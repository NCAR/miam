// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/util/types.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `DissolvedReaction`, mirroring `micm::ProcessSet`.
  /// @details Indices live in `SparseMatrixPolicy::VectorType<Index>` so they are
  ///          device-accessible via `.GetView()`. A POD `Views` bundle is captured by
  ///          value in every `MICM_LAMBDA`.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class DissolvedReactionSet
  {
   public:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;
    template<class U>
    using VectorView = typename Vector<U>::ConstViewType;

    struct Views
    {
      VectorView<micm::Index> reactant_indices_;
      VectorView<micm::Index> product_indices_;
      VectorView<micm::Index> solvent_indices_;
      VectorView<micm::Index> jacobian_flat_ids_;
      micm::Index k_state_parameter_index_;
      micm::Index num_phases_;
      micm::Index num_reactants_;
      micm::Index num_products_;
      micm::Real solvent_floor_;
      micm::Real min_halflife_;
      bool capped_;

      Views() = default;

      Views(
          const Vector<micm::Index>& reactant_indices,
          const Vector<micm::Index>& product_indices,
          const Vector<micm::Index>& solvent_indices,
          const Vector<micm::Index>& jacobian_flat_ids,
          micm::Index k_state_parameter_index,
          micm::Index num_phases,
          micm::Index num_reactants,
          micm::Index num_products,
          micm::Real solvent_floor,
          micm::Real min_halflife)
          : reactant_indices_(reactant_indices.GetView()),
            product_indices_(product_indices.GetView()),
            solvent_indices_(solvent_indices.GetView()),
            jacobian_flat_ids_(jacobian_flat_ids.GetView()),
            k_state_parameter_index_(k_state_parameter_index),
            num_phases_(num_phases),
            num_reactants_(num_reactants),
            num_products_(num_products),
            solvent_floor_(solvent_floor),
            min_halflife_(min_halflife),
            capped_(min_halflife > 0.0)
      {
      }
    };

    DissolvedReactionSet() = default;

    DissolvedReactionSet(
        const DissolvedReaction& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      const micm::Index k_state_parameter_index =
          static_cast<micm::Index>(LookupParameterIndex(config, state_parameter_indices));
      const micm::Index num_reactants = static_cast<micm::Index>(config.reactants_.size());
      const micm::Index num_products = static_cast<micm::Index>(config.products_.size());
      const micm::Real solvent_floor = config.solvent_floor_;
      const micm::Real min_halflife = config.min_halflife_;

      auto phase_it = phase_prefixes.find(config.phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedReactionSet: phase " + config.phase_.name_ + " not found in phase_prefixes for process " +
                config.uuid_);
      const auto& prefixes = phase_it->second;
      const micm::Index num_phases = static_cast<micm::Index>(prefixes.size());

      const micm::Index pairs_per_phase = (num_reactants + 1) * (num_reactants + num_products);

      std::vector<micm::Index> reactant_indices_host(num_phases * num_reactants, 0);
      std::vector<micm::Index> product_indices_host(num_phases * num_products, 0);
      std::vector<micm::Index> solvent_indices_host(num_phases, 0);
      std::vector<micm::Index> jacobian_flat_ids_host(num_phases * pairs_per_phase, 0);

      micm::Index i_phase = 0;
      for (const auto& prefix : prefixes)
      {
        for (micm::Index r = 0; r < num_reactants; ++r)
          reactant_indices_host[i_phase * num_reactants + r] = static_cast<micm::Index>(
              LookupSpecies(state_variable_indices, prefix, config.phase_.name_, config.reactants_[r].name_));
        for (micm::Index p = 0; p < num_products; ++p)
          product_indices_host[i_phase * num_products + p] = static_cast<micm::Index>(
              LookupSpecies(state_variable_indices, prefix, config.phase_.name_, config.products_[p].name_));
        solvent_indices_host[i_phase] = static_cast<micm::Index>(
            LookupSpecies(state_variable_indices, prefix, config.phase_.name_, config.solvent_.name_));

        micm::Index pair = 0;
        for (micm::Index i_ind = 0; i_ind < num_reactants; ++i_ind)
        {
          const micm::Index ind_idx = reactant_indices_host[i_phase * num_reactants + i_ind];
          for (micm::Index i_dep = 0; i_dep < num_reactants; ++i_dep)
            jacobian_flat_ids_host[i_phase * pairs_per_phase + pair++] = static_cast<micm::Index>(
                jacobian.VectorIndex(0, reactant_indices_host[i_phase * num_reactants + i_dep], ind_idx));
          for (micm::Index i_dep = 0; i_dep < num_products; ++i_dep)
            jacobian_flat_ids_host[i_phase * pairs_per_phase + pair++] = static_cast<micm::Index>(
                jacobian.VectorIndex(0, product_indices_host[i_phase * num_products + i_dep], ind_idx));
        }
        const micm::Index solv_idx = solvent_indices_host[i_phase];
        for (micm::Index i_dep = 0; i_dep < num_reactants; ++i_dep)
          jacobian_flat_ids_host[i_phase * pairs_per_phase + pair++] = static_cast<micm::Index>(
              jacobian.VectorIndex(0, reactant_indices_host[i_phase * num_reactants + i_dep], solv_idx));
        for (micm::Index i_dep = 0; i_dep < num_products; ++i_dep)
          jacobian_flat_ids_host[i_phase * pairs_per_phase + pair++] = static_cast<micm::Index>(
              jacobian.VectorIndex(0, product_indices_host[i_phase * num_products + i_dep], solv_idx));
        ++i_phase;
      }

      reactant_indices_ = Vector<micm::Index>(std::move(reactant_indices_host));
      product_indices_ = Vector<micm::Index>(std::move(product_indices_host));
      solvent_indices_ = Vector<micm::Index>(std::move(solvent_indices_host));
      jacobian_flat_ids_ = Vector<micm::Index>(std::move(jacobian_flat_ids_host));

      reactant_indices_.CopyToDevice();
      product_indices_.CopyToDevice();
      solvent_indices_.CopyToDevice();
      jacobian_flat_ids_.CopyToDevice();

      views_ = Views(
          reactant_indices_,
          product_indices_,
          solvent_indices_,
          jacobian_flat_ids_,
          k_state_parameter_index,
          num_phases,
          num_reactants,
          num_products,
          solvent_floor,
          min_halflife);
    }

    /// @brief Adds forcing contributions from every phase instance of this reaction into `forcing`.
    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const
    {
      const auto& views = views_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& params_view,
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ViewType& forcing_view)
          {
            const micm::Index num_reactants = views.num_reactants_;
            const micm::Index num_products = views.num_products_;
            const micm::Index k_index = views.k_state_parameter_index_;
            const micm::Real eps = views.solvent_floor_;
            const micm::Real t_half = views.min_halflife_;
            const bool capped = views.capped_;

            for (micm::Index phase = 0; phase < views.num_phases_; ++phase)
            {
              const micm::Index solvent_idx = views.solvent_indices_[phase];
              auto rate = forcing_view.GetRowVariable();
              forcing_view.ForEachRowStrict(
                  [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& out)
                  { out = k * solvent / std::pow(solvent + eps, num_reactants); },
                  params_view.GetConstColumnView(k_index),
                  state_view.GetConstColumnView(solvent_idx),
                  rate);
              for (micm::Index r = 0; r < num_reactants; ++r)
              {
                const micm::Index reactant_idx = views.reactant_indices_[phase * num_reactants + r];
                forcing_view.ForEachRowStrict(
                    [](const micm::Real& reactant, micm::Real& out) { out *= reactant; },
                    state_view.GetConstColumnView(reactant_idx),
                    rate);
              }

              if (capped)
              {
                auto accum = forcing_view.GetRowVariable();
                {
                  const micm::Index r0_idx = views.reactant_indices_[phase * num_reactants + 0];
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& R, micm::Real& acc)
                      { acc = std::pow(std::max(R, micm::Real(1.0e-300)), -micm::Real(10.0)); },
                      state_view.GetConstColumnView(r0_idx),
                      accum);
                }
                for (micm::Index r = 1; r < num_reactants; ++r)
                {
                  const micm::Index r_idx = views.reactant_indices_[phase * num_reactants + r];
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& R, micm::Real& acc)
                      { acc += std::pow(std::max(R, micm::Real(1.0e-300)), -micm::Real(10.0)); },
                      state_view.GetConstColumnView(r_idx),
                      accum);
                }
                forcing_view.ForEachRowStrict(
                    [t_half](micm::Real& out, micm::Real& acc)
                    {
                      const micm::Real c_min = std::pow(acc, -1.0 / micm::Real(10.0));
                      const micm::Real r_max = c_min / t_half;
                      if (r_max > micm::Real(1.0e-300))
                        out = r_max * std::tanh(out / r_max);
                    },
                    rate,
                    accum);
              }

              for (micm::Index r = 0; r < num_reactants; ++r)
              {
                const micm::Index reactant_idx = views.reactant_indices_[phase * num_reactants + r];
                forcing_view.ForEachRowStrict(
                    [](const micm::Real& rate, micm::Real& forcing) { forcing -= rate; },
                    rate,
                    forcing_view.GetColumnView(reactant_idx));
              }
              for (micm::Index p = 0; p < num_products; ++p)
              {
                const micm::Index product_idx = views.product_indices_[phase * num_products + p];
                forcing_view.ForEachRowStrict(
                    [](const micm::Real& rate, micm::Real& forcing) { forcing += rate; },
                    rate,
                    forcing_view.GetColumnView(product_idx));
              }
            }
          },
          state_parameters,
          state_variables,
          forcing)(state_parameters, state_variables, forcing);
    }

    /// @brief Subtracts Jacobian contributions from every phase instance of this reaction into `jacobian`.
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      const auto& views = views_;
      SparseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& params_view,
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename SparseMatrixPolicy::ViewType& jac_view)
          {
            const micm::Index num_reactants = views.num_reactants_;
            const micm::Index num_products = views.num_products_;
            const micm::Index k_index = views.k_state_parameter_index_;
            const micm::Index pairs_per_phase = (num_reactants + 1) * (num_reactants + num_products);
            const micm::Real eps = views.solvent_floor_;
            const micm::Real t_half = views.min_halflife_;
            const bool capped = views.capped_;

            for (micm::Index phase = 0; phase < views.num_phases_; ++phase)
            {
              const micm::Index solvent_idx = views.solvent_indices_[phase];
              auto d_rate_d_ind = jac_view.GetBlockVariable();
              auto raw_rate = jac_view.GetBlockVariable();
              auto sech2_var = jac_view.GetBlockVariable();
              auto corr_var = jac_view.GetBlockVariable();
              auto c_min_var = jac_view.GetBlockVariable();
              micm::Index pair = phase * pairs_per_phase;

              if (capped)
              {
                jac_view.ForEachBlockStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& rr)
                    { rr = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                    raw_rate);
                for (micm::Index r = 0; r < num_reactants; ++r)
                {
                  const micm::Index r_idx = views.reactant_indices_[phase * num_reactants + r];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& reactant, micm::Real& rr) { rr *= reactant; },
                      state_view.GetConstColumnView(r_idx),
                      raw_rate);
                }
                {
                  const micm::Index r0_idx = views.reactant_indices_[phase * num_reactants + 0];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& R, micm::Real& cm)
                      { cm = std::pow(std::max(R, micm::Real(1.0e-300)), -micm::Real(10.0)); },
                      state_view.GetConstColumnView(r0_idx),
                      c_min_var);
                }
                for (micm::Index r = 1; r < num_reactants; ++r)
                {
                  const micm::Index r_idx = views.reactant_indices_[phase * num_reactants + r];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& R, micm::Real& cm)
                      { cm += std::pow(std::max(R, micm::Real(1.0e-300)), -micm::Real(10.0)); },
                      state_view.GetConstColumnView(r_idx),
                      c_min_var);
                }
                jac_view.ForEachBlockStrict(
                    [t_half](micm::Real& rr, micm::Real& cm, micm::Real& s2, micm::Real& cr)
                    {
                      cm = std::pow(cm, -1.0 / micm::Real(10.0));
                      const micm::Real r_max = cm / t_half;
                      if (r_max > micm::Real(1.0e-300))
                      {
                        const micm::Real u = rr / r_max;
                        const micm::Real th = std::tanh(u);
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

              for (micm::Index i_ind = 0; i_ind < num_reactants; ++i_ind)
              {
                jac_view.ForEachBlockStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& partial)
                    { partial = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                    d_rate_d_ind);
                for (micm::Index r = 0; r < num_reactants; ++r)
                {
                  if (r == i_ind)
                    continue;
                  const micm::Index r_idx = views.reactant_indices_[phase * num_reactants + r];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                      state_view.GetConstColumnView(r_idx),
                      d_rate_d_ind);
                }
                if (capped)
                {
                  const micm::Index i_ind_idx = views.reactant_indices_[phase * num_reactants + i_ind];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& s2,
                         const micm::Real& cr,
                         const micm::Real& cm,
                         const micm::Real& R,
                         micm::Real& partial)
                      {
                        const micm::Real ratio = cm / std::max(R, micm::Real(1.0e-300));
                        partial = s2 * partial + cr * std::pow(ratio, micm::Real(10.0) + 1.0);
                      },
                      sech2_var,
                      corr_var,
                      c_min_var,
                      state_view.GetConstColumnView(i_ind_idx),
                      d_rate_d_ind);
                }
                for (micm::Index i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                  const micm::Index flat = views.jacobian_flat_ids_[pair++];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                      d_rate_d_ind,
                      jac_view.GetBlockView(flat));
                }
                for (micm::Index i_dep = 0; i_dep < num_products; ++i_dep)
                {
                  const micm::Index flat = views.jacobian_flat_ids_[pair++];
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                      d_rate_d_ind,
                      jac_view.GetBlockView(flat));
                }
              }

              jac_view.ForEachBlockStrict(
                  [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& partial)
                  {
                    partial = k * (eps + (1.0 - static_cast<micm::Real>(num_reactants)) * solvent) /
                              std::pow(solvent + eps, num_reactants + 1);
                  },
                  params_view.GetConstColumnView(k_index),
                  state_view.GetConstColumnView(solvent_idx),
                  d_rate_d_ind);
              for (micm::Index r = 0; r < num_reactants; ++r)
              {
                const micm::Index r_idx = views.reactant_indices_[phase * num_reactants + r];
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                    state_view.GetConstColumnView(r_idx),
                    d_rate_d_ind);
              }
              if (capped)
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& s2, micm::Real& partial) { partial *= s2; }, sech2_var, d_rate_d_ind);
              for (micm::Index i_dep = 0; i_dep < num_reactants; ++i_dep)
              {
                const micm::Index flat = views.jacobian_flat_ids_[pair++];
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                    d_rate_d_ind,
                    jac_view.GetBlockView(flat));
              }
              for (micm::Index i_dep = 0; i_dep < num_products; ++i_dep)
              {
                const micm::Index flat = views.jacobian_flat_ids_[pair++];
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                    d_rate_d_ind,
                    jac_view.GetBlockView(flat));
              }
            }
          },
          state_parameters,
          state_variables,
          jacobian)(state_parameters, state_variables, jacobian);
    }

   private:
    static std::size_t LookupParameterIndex(
        const DissolvedReaction& config,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices)
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
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
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

    Vector<micm::Index> reactant_indices_{};
    Vector<micm::Index> product_indices_{};
    Vector<micm::Index> solvent_indices_{};
    Vector<micm::Index> jacobian_flat_ids_{};
    Views views_{};
  };
}  // namespace miam
