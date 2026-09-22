// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `DissolvedReaction`, mirroring `micm::ProcessSet`.
  /// @details Non-templated data-holding class populated once at Finalize time (which is
  ///          when the sparse Jacobian pattern is available). Every device-facing solve-time
  ///          method captures only trivially-copyable scalars into `MICM_LAMBDA` and iterates
  ///          the phase-instance / reactant / product indices on the host \u2014 matching MICM's
  ///          `stub_aerosol_1` external model. Rate-constant evaluation is handled separately
  ///          by the Layer-D `ConstantsBucket`; this class only reads the resolved `k` column
  ///          from the state-parameters matrix.
  class DissolvedReactionSet
  {
   public:
    DissolvedReactionSet() = default;

    template<typename SparseMatrixPolicy>
    DissolvedReactionSet(
        const DissolvedReaction& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,
        const auto& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      k_state_parameter_index_ = LookupParameterIndex(config, state_parameter_indices);
      num_reactants_ = config.reactants_.size();
      num_products_ = config.products_.size();
      solvent_floor_ = config.solvent_floor_;
      min_halflife_ = config.min_halflife_;

      auto phase_it = phase_prefixes.find(config.phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedReactionSet: phase " + config.phase_.name_ + " not found in phase_prefixes for process " +
                config.uuid_);
      const auto& prefixes = phase_it->second;
      num_phases_ = prefixes.size();

      // Layout: reactant_indices_[phase * num_reactants_ + r], product_indices_[phase * num_products_ + p],
      //         solvent_indices_[phase], jacobian_flat_ids_[phase * pairs_per_phase + pair].
      const std::size_t pairs_per_phase = (num_reactants_ + 1) * (num_reactants_ + num_products_);
      reactant_indices_.assign(num_phases_ * num_reactants_, 0);
      product_indices_.assign(num_phases_ * num_products_, 0);
      solvent_indices_.assign(num_phases_, 0);
      jacobian_flat_ids_.assign(num_phases_ * pairs_per_phase, 0);

      std::size_t i_phase = 0;
      for (const auto& prefix : prefixes)
      {
        for (std::size_t r = 0; r < num_reactants_; ++r)
          reactant_indices_[i_phase * num_reactants_ + r] =
              LookupSpecies(state_variable_indices, prefix, config.phase_.name_, config.reactants_[r].name_);
        for (std::size_t p = 0; p < num_products_; ++p)
          product_indices_[i_phase * num_products_ + p] =
              LookupSpecies(state_variable_indices, prefix, config.phase_.name_, config.products_[p].name_);
        solvent_indices_[i_phase] =
            LookupSpecies(state_variable_indices, prefix, config.phase_.name_, config.solvent_.name_);

        std::size_t pair = 0;
        for (std::size_t i_ind = 0; i_ind < num_reactants_; ++i_ind)
        {
          const std::size_t ind_idx = reactant_indices_[i_phase * num_reactants_ + i_ind];
          for (std::size_t i_dep = 0; i_dep < num_reactants_; ++i_dep)
            jacobian_flat_ids_[i_phase * pairs_per_phase + pair++] =
                jacobian.VectorIndex(0, reactant_indices_[i_phase * num_reactants_ + i_dep], ind_idx);
          for (std::size_t i_dep = 0; i_dep < num_products_; ++i_dep)
            jacobian_flat_ids_[i_phase * pairs_per_phase + pair++] =
                jacobian.VectorIndex(0, product_indices_[i_phase * num_products_ + i_dep], ind_idx);
        }
        const std::size_t solv_idx = solvent_indices_[i_phase];
        for (std::size_t i_dep = 0; i_dep < num_reactants_; ++i_dep)
          jacobian_flat_ids_[i_phase * pairs_per_phase + pair++] =
              jacobian.VectorIndex(0, reactant_indices_[i_phase * num_reactants_ + i_dep], solv_idx);
        for (std::size_t i_dep = 0; i_dep < num_products_; ++i_dep)
          jacobian_flat_ids_[i_phase * pairs_per_phase + pair++] =
              jacobian.VectorIndex(0, product_indices_[i_phase * num_products_ + i_dep], solv_idx);
        ++i_phase;
      }
    }

    /// @brief Adds forcing contributions from every phase instance of this reaction into `forcing`.
    template<typename DenseMatrixPolicy>
    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const
    {
      const std::size_t num_reactants = num_reactants_;
      const std::size_t num_products = num_products_;
      const std::size_t k_index = k_state_parameter_index_;
      const double eps = solvent_floor_;
      const double t_half = min_halflife_;
      const bool capped = t_half > 0.0;

      for (std::size_t phase = 0; phase < num_phases_; ++phase)
      {
        const std::size_t solvent_idx = solvent_indices_[phase];
        DenseMatrixPolicy::Function(
            [this, phase, k_index, solvent_idx, num_reactants, num_products, eps, t_half, capped](
                auto&& params, auto&& vars, auto&& forcing_view)
            {
              auto rate = forcing_view.GetRowVariable();
              params.ForEachRow(
                  [num_reactants, eps](const double& k, const double& solvent, double& out)
                  { out = k * solvent / std::pow(solvent + eps, num_reactants); },
                  params.GetConstColumnView(k_index),
                  vars.GetConstColumnView(solvent_idx),
                  rate);
              for (std::size_t r = 0; r < num_reactants; ++r)
              {
                const std::size_t reactant_idx = reactant_indices_[phase * num_reactants + r];
                params.ForEachRow(
                    [](const double& reactant, double& out) { out *= reactant; },
                    vars.GetConstColumnView(reactant_idx),
                    rate);
              }

              if (capped)
              {
                auto accum = forcing_view.GetRowVariable();
                {
                  const std::size_t r0_idx = reactant_indices_[phase * num_reactants + 0];
                  params.ForEachRow(
                      [](const double& R, double& acc) { acc = std::pow(std::max(R, kSoftMinFloor), -kSoftMinP); },
                      vars.GetConstColumnView(r0_idx),
                      accum);
                }
                for (std::size_t r = 1; r < num_reactants; ++r)
                {
                  const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                  params.ForEachRow(
                      [](const double& R, double& acc) { acc += std::pow(std::max(R, kSoftMinFloor), -kSoftMinP); },
                      vars.GetConstColumnView(r_idx),
                      accum);
                }
                params.ForEachRow(
                    [t_half](double& out, double& acc)
                    {
                      const double c_min = std::pow(acc, -1.0 / kSoftMinP);
                      const double r_max = c_min / t_half;
                      if (r_max > kSoftMinFloor)
                        out = r_max * std::tanh(out / r_max);
                    },
                    rate,
                    accum);
              }

              for (std::size_t r = 0; r < num_reactants; ++r)
              {
                const std::size_t reactant_idx = reactant_indices_[phase * num_reactants + r];
                params.ForEachRow(
                    [](const double& rate, double& forcing) { forcing -= rate; },
                    rate,
                    forcing_view.GetColumnView(reactant_idx));
              }
              for (std::size_t p = 0; p < num_products; ++p)
              {
                const std::size_t product_idx = product_indices_[phase * num_products + p];
                params.ForEachRow(
                    [](const double& rate, double& forcing) { forcing += rate; },
                    rate,
                    forcing_view.GetColumnView(product_idx));
              }
            },
            state_parameters,
            state_variables,
            forcing)(state_parameters, state_variables, forcing);
      }
    }

    /// @brief Subtracts Jacobian contributions from every phase instance of this reaction into `jacobian`.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      const std::size_t num_reactants = num_reactants_;
      const std::size_t num_products = num_products_;
      const std::size_t k_index = k_state_parameter_index_;
      const std::size_t pairs_per_phase = (num_reactants_ + 1) * (num_reactants_ + num_products_);
      const double eps = solvent_floor_;
      const double t_half = min_halflife_;
      const bool capped = t_half > 0.0;

      for (std::size_t phase = 0; phase < num_phases_; ++phase)
      {
        const std::size_t solvent_idx = solvent_indices_[phase];
        SparseMatrixPolicy::Function(
            [this, phase, k_index, solvent_idx, num_reactants, num_products, pairs_per_phase, eps, t_half, capped](
                auto&& params, auto&& vars, auto&& jacobian_values)
            {
              auto d_rate_d_ind = jacobian_values.GetBlockVariable();
              auto raw_rate = jacobian_values.GetBlockVariable();
              auto sech2_var = jacobian_values.GetBlockVariable();
              auto corr_var = jacobian_values.GetBlockVariable();
              auto c_min_var = jacobian_values.GetBlockVariable();
              std::size_t pair = phase * pairs_per_phase;

              if (capped)
              {
                jacobian_values.ForEachBlock(
                    [num_reactants, eps](const double& k, const double& solvent, double& rr)
                    { rr = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params.GetConstColumnView(k_index),
                    vars.GetConstColumnView(solvent_idx),
                    raw_rate);
                for (std::size_t r = 0; r < num_reactants; ++r)
                {
                  const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                  jacobian_values.ForEachBlock(
                      [](const double& reactant, double& rr) { rr *= reactant; },
                      vars.GetConstColumnView(r_idx),
                      raw_rate);
                }
                {
                  const std::size_t r0_idx = reactant_indices_[phase * num_reactants + 0];
                  jacobian_values.ForEachBlock(
                      [](const double& R, double& cm) { cm = std::pow(std::max(R, kSoftMinFloor), -kSoftMinP); },
                      vars.GetConstColumnView(r0_idx),
                      c_min_var);
                }
                for (std::size_t r = 1; r < num_reactants; ++r)
                {
                  const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                  jacobian_values.ForEachBlock(
                      [](const double& R, double& cm) { cm += std::pow(std::max(R, kSoftMinFloor), -kSoftMinP); },
                      vars.GetConstColumnView(r_idx),
                      c_min_var);
                }
                jacobian_values.ForEachBlock(
                    [t_half](double& rr, double& cm, double& s2, double& cr)
                    {
                      cm = std::pow(cm, -1.0 / kSoftMinP);
                      const double r_max = cm / t_half;
                      if (r_max > kSoftMinFloor)
                      {
                        const double u = rr / r_max;
                        const double th = std::tanh(u);
                        s2 = 1.0 - th * th;
                        cr = (th - u * s2) / t_half;
                      }
                      else
                      {
                        s2 = 1.0;
                        cr = 0.0;
                      }
                    },
                    raw_rate,
                    c_min_var,
                    sech2_var,
                    corr_var);
              }

              for (std::size_t i_ind = 0; i_ind < num_reactants; ++i_ind)
              {
                jacobian_values.ForEachBlock(
                    [num_reactants, eps](const double& k, const double& solvent, double& partial)
                    { partial = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params.GetConstColumnView(k_index),
                    vars.GetConstColumnView(solvent_idx),
                    d_rate_d_ind);
                for (std::size_t r = 0; r < num_reactants; ++r)
                {
                  if (r == i_ind)
                    continue;
                  const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                  jacobian_values.ForEachBlock(
                      [](const double& reactant, double& partial) { partial *= reactant; },
                      vars.GetConstColumnView(r_idx),
                      d_rate_d_ind);
                }
                if (capped)
                {
                  const std::size_t i_ind_idx = reactant_indices_[phase * num_reactants + i_ind];
                  jacobian_values.ForEachBlock(
                      [](const double& s2, const double& cr, const double& cm, const double& R, double& partial)
                      {
                        const double ratio = cm / std::max(R, kSoftMinFloor);
                        partial = s2 * partial + cr * std::pow(ratio, kSoftMinP + 1.0);
                      },
                      sech2_var,
                      corr_var,
                      c_min_var,
                      vars.GetConstColumnView(i_ind_idx),
                      d_rate_d_ind);
                }
                for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                  const std::size_t flat = jacobian_flat_ids_[pair++];
                  jacobian_values.ForEachBlock(
                      [](const double& partial, double& jac) { jac += partial; },
                      d_rate_d_ind,
                      jacobian_values.GetBlockView(flat));
                }
                for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                  const std::size_t flat = jacobian_flat_ids_[pair++];
                  jacobian_values.ForEachBlock(
                      [](const double& partial, double& jac) { jac -= partial; },
                      d_rate_d_ind,
                      jacobian_values.GetBlockView(flat));
                }
              }

              jacobian_values.ForEachBlock(
                  [num_reactants, eps](const double& k, const double& solvent, double& partial) {
                    partial = k * (eps + (1.0 - static_cast<int>(num_reactants)) * solvent) /
                              std::pow(solvent + eps, num_reactants + 1);
                  },
                  params.GetConstColumnView(k_index),
                  vars.GetConstColumnView(solvent_idx),
                  d_rate_d_ind);
              for (std::size_t r = 0; r < num_reactants; ++r)
              {
                const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                jacobian_values.ForEachBlock(
                    [](const double& reactant, double& partial) { partial *= reactant; },
                    vars.GetConstColumnView(r_idx),
                    d_rate_d_ind);
              }
              if (capped)
                jacobian_values.ForEachBlock(
                    [](const double& s2, double& partial) { partial *= s2; }, sech2_var, d_rate_d_ind);
              for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
              {
                const std::size_t flat = jacobian_flat_ids_[pair++];
                jacobian_values.ForEachBlock(
                    [](const double& partial, double& jac) { jac += partial; },
                    d_rate_d_ind,
                    jacobian_values.GetBlockView(flat));
              }
              for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
              {
                const std::size_t flat = jacobian_flat_ids_[pair++];
                jacobian_values.ForEachBlock(
                    [](const double& partial, double& jac) { jac -= partial; },
                    d_rate_d_ind,
                    jacobian_values.GetBlockView(flat));
              }
            },
            state_parameters,
            state_variables,
            jacobian)(state_parameters, state_variables, jacobian);
      }
    }

   private:
    /// Soft-min exponent for rate capping; matches `DissolvedReaction::kSoftMinP`.
    static constexpr double kSoftMinP = 10.0;
    /// Tiny floor to prevent `pow(0, -p)` overflow; matches `DissolvedReaction::kSoftMinFloor`.
    static constexpr double kSoftMinFloor = 1.0e-300;

    static std::size_t LookupParameterIndex(const DissolvedReaction& config, const auto& state_parameter_indices)
    {
      const std::string key = config.phase_.name_ + "." + config.uuid_ + ".k";
      auto it = state_parameter_indices.find(key);
      if (it == state_parameter_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "DissolvedReactionSet: rate constant parameter " + key + " not found");
      return it->second;
    }

    static std::size_t LookupSpecies(
        const auto& state_variable_indices,
        const std::string& prefix,
        const std::string& phase_name,
        const std::string& species_name)
    {
      const std::string key = prefix + "." + phase_name + "." + species_name;
      auto it = state_variable_indices.find(key);
      if (it == state_variable_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_VARIABLE,
            "DissolvedReactionSet: state variable " + key + " not found");
      return it->second;
    }

    std::vector<std::size_t> reactant_indices_;
    std::vector<std::size_t> product_indices_;
    std::vector<std::size_t> solvent_indices_;
    std::vector<std::size_t> jacobian_flat_ids_;
    std::size_t k_state_parameter_index_ = 0;
    std::size_t num_phases_ = 0;
    std::size_t num_reactants_ = 0;
    std::size_t num_products_ = 0;
    double solvent_floor_ = 1.0e-20;
    double min_halflife_ = 0.0;
  };
}  // namespace miam
