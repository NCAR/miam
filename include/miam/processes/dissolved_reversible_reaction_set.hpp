// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `DissolvedReversibleReaction`, mirroring `micm::ProcessSet`.
  /// @details Non-templated data-holding class populated once at Finalize time (which is
  ///          when the sparse Jacobian pattern is available). Every device-facing solve-time
  ///          method captures only trivially-copyable scalars into `MICM_LAMBDA` and iterates
  ///          the phase-instance / reactant / product indices on the host — matching MICM's
  ///          `stub_aerosol_1` external model.
  class DissolvedReversibleReactionSet
  {
   public:
    DissolvedReversibleReactionSet() = default;

    template<typename SparseMatrixPolicy>
    DissolvedReversibleReactionSet(
        const DissolvedReversibleReaction& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,
        const auto& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      const std::string fwd_key = config.phase_.name_ + "." + config.uuid_ + ".k_forward";
      const std::string rev_key = config.phase_.name_ + "." + config.uuid_ + ".k_reverse";
      k_forward_parameter_index_ = LookupParameter(state_parameter_indices, fwd_key);
      k_reverse_parameter_index_ = LookupParameter(state_parameter_indices, rev_key);
      num_reactants_ = config.reactants_.size();
      num_products_ = config.products_.size();
      solvent_floor_ = config.solvent_floor_;

      auto phase_it = phase_prefixes.find(config.phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedReversibleReactionSet: phase " + config.phase_.name_ +
                " not found in phase_prefixes for process " + config.uuid_);
      const auto& prefixes = phase_it->second;
      num_phases_ = prefixes.size();

      const std::size_t pairs_per_phase = (num_reactants_ + num_products_ + 1) * (num_reactants_ + num_products_);
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
        for (std::size_t i_ind = 0; i_ind < num_products_; ++i_ind)
        {
          const std::size_t ind_idx = product_indices_[i_phase * num_products_ + i_ind];
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
      const std::size_t k_fwd = k_forward_parameter_index_;
      const std::size_t k_rev = k_reverse_parameter_index_;
      const double eps = solvent_floor_;

      for (std::size_t phase = 0; phase < num_phases_; ++phase)
      {
        const std::size_t solvent_idx = solvent_indices_[phase];
        DenseMatrixPolicy::Function(
            [this, phase, k_fwd, k_rev, solvent_idx, num_reactants, num_products, eps](
                auto&& params, auto&& vars, auto&& forcing_view)
            {
              auto forward_rate = forcing_view.GetRowVariable();
              auto reverse_rate = forcing_view.GetRowVariable();
              params.ForEachRow(
                  [num_reactants, num_products, eps](
                      const double& k_f, const double& k_r, const double& solvent, double& fwd, double& rev)
                  {
                    fwd = k_f * solvent / std::pow(solvent + eps, num_reactants);
                    rev = k_r * solvent / std::pow(solvent + eps, num_products);
                  },
                  params.GetConstColumnView(k_fwd),
                  params.GetConstColumnView(k_rev),
                  vars.GetConstColumnView(solvent_idx),
                  forward_rate,
                  reverse_rate);
              for (std::size_t r = 0; r < num_reactants; ++r)
              {
                const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                params.ForEachRow(
                    [](const double& reactant, double& fwd) { fwd *= reactant; },
                    vars.GetConstColumnView(r_idx),
                    forward_rate);
              }
              for (std::size_t p = 0; p < num_products; ++p)
              {
                const std::size_t p_idx = product_indices_[phase * num_products + p];
                params.ForEachRow(
                    [](const double& product, double& rev) { rev *= product; },
                    vars.GetConstColumnView(p_idx),
                    reverse_rate);
              }
              for (std::size_t r = 0; r < num_reactants; ++r)
              {
                const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                params.ForEachRow(
                    [](const double& fwd, const double& rev, double& forcing)
                    {
                      forcing -= fwd;
                      forcing += rev;
                    },
                    forward_rate,
                    reverse_rate,
                    forcing_view.GetColumnView(r_idx));
              }
              for (std::size_t p = 0; p < num_products; ++p)
              {
                const std::size_t p_idx = product_indices_[phase * num_products + p];
                params.ForEachRow(
                    [](const double& fwd, const double& rev, double& forcing)
                    {
                      forcing += fwd;
                      forcing -= rev;
                    },
                    forward_rate,
                    reverse_rate,
                    forcing_view.GetColumnView(p_idx));
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
      const std::size_t k_fwd = k_forward_parameter_index_;
      const std::size_t k_rev = k_reverse_parameter_index_;
      const std::size_t pairs_per_phase = (num_reactants_ + num_products_ + 1) * (num_reactants_ + num_products_);
      const double eps = solvent_floor_;

      for (std::size_t phase = 0; phase < num_phases_; ++phase)
      {
        const std::size_t solvent_idx = solvent_indices_[phase];
        SparseMatrixPolicy::Function(
            [this, phase, k_fwd, k_rev, solvent_idx, num_reactants, num_products, pairs_per_phase, eps](
                auto&& params, auto&& vars, auto&& jacobian_values)
            {
              auto d_fwd = jacobian_values.GetBlockVariable();
              auto d_rev = jacobian_values.GetBlockVariable();
              std::size_t pair = phase * pairs_per_phase;

              for (std::size_t i_ind = 0; i_ind < num_reactants; ++i_ind)
              {
                jacobian_values.ForEachBlock(
                    [num_reactants, eps](const double& k_f, const double& solvent, double& partial)
                    { partial = k_f * solvent / std::pow(solvent + eps, num_reactants); },
                    params.GetConstColumnView(k_fwd),
                    vars.GetConstColumnView(solvent_idx),
                    d_fwd);
                for (std::size_t r = 0; r < num_reactants; ++r)
                {
                  if (r == i_ind)
                    continue;
                  const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                  jacobian_values.ForEachBlock(
                      [](const double& reactant, double& partial) { partial *= reactant; },
                      vars.GetConstColumnView(r_idx),
                      d_fwd);
                }
                for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                  const std::size_t flat = jacobian_flat_ids_[pair++];
                  jacobian_values.ForEachBlock(
                      [](const double& partial, double& jac) { jac += partial; },
                      d_fwd,
                      jacobian_values.GetBlockView(flat));
                }
                for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                  const std::size_t flat = jacobian_flat_ids_[pair++];
                  jacobian_values.ForEachBlock(
                      [](const double& partial, double& jac) { jac -= partial; },
                      d_fwd,
                      jacobian_values.GetBlockView(flat));
                }
              }

              for (std::size_t i_ind = 0; i_ind < num_products; ++i_ind)
              {
                jacobian_values.ForEachBlock(
                    [num_products, eps](const double& k_r, const double& solvent, double& partial)
                    { partial = k_r * solvent / std::pow(solvent + eps, num_products); },
                    params.GetConstColumnView(k_rev),
                    vars.GetConstColumnView(solvent_idx),
                    d_rev);
                for (std::size_t p = 0; p < num_products; ++p)
                {
                  if (p == i_ind)
                    continue;
                  const std::size_t p_idx = product_indices_[phase * num_products + p];
                  jacobian_values.ForEachBlock(
                      [](const double& product, double& partial) { partial *= product; },
                      vars.GetConstColumnView(p_idx),
                      d_rev);
                }
                for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                  const std::size_t flat = jacobian_flat_ids_[pair++];
                  jacobian_values.ForEachBlock(
                      [](const double& partial, double& jac) { jac -= partial; },
                      d_rev,
                      jacobian_values.GetBlockView(flat));
                }
                for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                  const std::size_t flat = jacobian_flat_ids_[pair++];
                  jacobian_values.ForEachBlock(
                      [](const double& partial, double& jac) { jac += partial; },
                      d_rev,
                      jacobian_values.GetBlockView(flat));
                }
              }

              jacobian_values.ForEachBlock(
                  [num_reactants, num_products, eps](
                      const double& k_f, const double& k_r, const double& solvent, double& fwd_p, double& rev_p) {
                    fwd_p = k_f * (eps + (1.0 - static_cast<int>(num_reactants)) * solvent) /
                            std::pow(solvent + eps, num_reactants + 1);
                    rev_p = k_r * (eps + (1.0 - static_cast<int>(num_products)) * solvent) /
                            std::pow(solvent + eps, num_products + 1);
                  },
                  params.GetConstColumnView(k_fwd),
                  params.GetConstColumnView(k_rev),
                  vars.GetConstColumnView(solvent_idx),
                  d_fwd,
                  d_rev);
              for (std::size_t r = 0; r < num_reactants; ++r)
              {
                const std::size_t r_idx = reactant_indices_[phase * num_reactants + r];
                jacobian_values.ForEachBlock(
                    [](const double& reactant, double& fwd_p) { fwd_p *= reactant; },
                    vars.GetConstColumnView(r_idx),
                    d_fwd);
              }
              for (std::size_t p = 0; p < num_products; ++p)
              {
                const std::size_t p_idx = product_indices_[phase * num_products + p];
                jacobian_values.ForEachBlock(
                    [](const double& product, double& rev_p) { rev_p *= product; },
                    vars.GetConstColumnView(p_idx),
                    d_rev);
              }
              for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
              {
                const std::size_t flat = jacobian_flat_ids_[pair++];
                jacobian_values.ForEachBlock(
                    [](const double& fwd_p, const double& rev_p, double& jac)
                    {
                      jac += fwd_p;
                      jac -= rev_p;
                    },
                    d_fwd,
                    d_rev,
                    jacobian_values.GetBlockView(flat));
              }
              for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
              {
                const std::size_t flat = jacobian_flat_ids_[pair++];
                jacobian_values.ForEachBlock(
                    [](const double& fwd_p, const double& rev_p, double& jac)
                    {
                      jac -= fwd_p;
                      jac += rev_p;
                    },
                    d_fwd,
                    d_rev,
                    jacobian_values.GetBlockView(flat));
              }
            },
            state_parameters,
            state_variables,
            jacobian)(state_parameters, state_variables, jacobian);
      }
    }

   private:
    static std::size_t LookupParameter(const auto& state_parameter_indices, const std::string& key)
    {
      auto it = state_parameter_indices.find(key);
      if (it == state_parameter_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "DissolvedReversibleReactionSet: state parameter " + key + " not found");
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
            "DissolvedReversibleReactionSet: state variable " + key + " not found");
      return it->second;
    }

    std::vector<std::size_t> reactant_indices_;
    std::vector<std::size_t> product_indices_;
    std::vector<std::size_t> solvent_indices_;
    std::vector<std::size_t> jacobian_flat_ids_;
    std::size_t k_forward_parameter_index_ = 0;
    std::size_t k_reverse_parameter_index_ = 0;
    std::size_t num_phases_ = 0;
    std::size_t num_reactants_ = 0;
    std::size_t num_products_ = 0;
    double solvent_floor_ = 1.0e-20;
  };
}  // namespace miam
