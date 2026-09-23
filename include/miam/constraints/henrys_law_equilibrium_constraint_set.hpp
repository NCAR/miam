// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/henrys_law_equilibrium_constraint.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/util/types.hpp>

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `HenrysLawEquilibriumConstraint`, mirroring `micm::ProcessSet`.
  /// @details G = HLC * R * T * f_v * [A_g] - [A_aq]; f_v = [S] * solvent_molecular_weight / solvent_density.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class HenrysLawEquilibriumConstraintSet
  {
   public:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;
    template<class U>
    using VectorView = typename Vector<U>::ConstViewType;

    struct Views
    {
      VectorView<micm::Index> aq_indices_;
      VectorView<micm::Index> solvent_indices_;
      VectorView<micm::Index> hlc_rt_indices_;
      VectorView<micm::Index> gas_jac_ids_;
      VectorView<micm::Index> aq_jac_ids_;
      VectorView<micm::Index> solvent_jac_ids_;
      micm::Index num_phases_;
      micm::Index gas_idx_;
      micm::Real molar_volume_;

      Views() = default;

      Views(
          const Vector<micm::Index>& aq_indices,
          const Vector<micm::Index>& solvent_indices,
          const Vector<micm::Index>& hlc_rt_indices,
          const Vector<micm::Index>& gas_jac_ids,
          const Vector<micm::Index>& aq_jac_ids,
          const Vector<micm::Index>& solvent_jac_ids,
          micm::Index num_phases,
          micm::Index gas_idx,
          micm::Real molar_volume)
          : aq_indices_(aq_indices.GetView()),
            solvent_indices_(solvent_indices.GetView()),
            hlc_rt_indices_(hlc_rt_indices.GetView()),
            gas_jac_ids_(gas_jac_ids.GetView()),
            aq_jac_ids_(aq_jac_ids.GetView()),
            solvent_jac_ids_(solvent_jac_ids.GetView()),
            num_phases_(num_phases),
            gas_idx_(gas_idx),
            molar_volume_(molar_volume)
      {
      }
    };

    HenrysLawEquilibriumConstraintSet() = default;

    HenrysLawEquilibriumConstraintSet(
        const HenrysLawEquilibriumConstraint& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      auto gas_it = state_variable_indices.find(config.gas_species_.name_);
      if (gas_it == state_variable_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_VARIABLE,
            "HenrysLawEquilibriumConstraintSet: gas species '" + config.gas_species_.name_ + "' not found");
      const auto gas_idx = static_cast<micm::Index>(gas_it->second);
      const micm::Real molar_volume = config.solvent_molecular_weight_ / config.solvent_density_;

      auto phase_it = phase_prefixes.find(config.condensed_phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "HenrysLawEquilibriumConstraintSet: phase " + config.condensed_phase_.name_ + " not found");

      const auto num_phases = static_cast<micm::Index>(phase_it->second.size());

      std::vector<micm::Index> aq_indices_host;
      std::vector<micm::Index> solvent_indices_host;
      std::vector<micm::Index> hlc_rt_indices_host;
      std::vector<micm::Index> gas_jac_ids_host;
      std::vector<micm::Index> aq_jac_ids_host;
      std::vector<micm::Index> solvent_jac_ids_host;

      for (const auto& prefix : phase_it->second)
      {
        const auto aq_idx = static_cast<micm::Index>(state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.condensed_species_.name_));
        const auto sol_idx = static_cast<micm::Index>(state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.solvent_.name_));
        aq_indices_host.push_back(aq_idx);
        solvent_indices_host.push_back(sol_idx);
        hlc_rt_indices_host.push_back(static_cast<micm::Index>(state_parameter_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.uuid_ + ".hlc_rt")));
        gas_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, gas_idx)));
        aq_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, aq_idx)));
        solvent_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, sol_idx)));
      }

      aq_indices_ = Vector<micm::Index>(std::move(aq_indices_host));
      solvent_indices_ = Vector<micm::Index>(std::move(solvent_indices_host));
      hlc_rt_indices_ = Vector<micm::Index>(std::move(hlc_rt_indices_host));
      gas_jac_ids_ = Vector<micm::Index>(std::move(gas_jac_ids_host));
      aq_jac_ids_ = Vector<micm::Index>(std::move(aq_jac_ids_host));
      solvent_jac_ids_ = Vector<micm::Index>(std::move(solvent_jac_ids_host));

      aq_indices_.CopyToDevice();
      solvent_indices_.CopyToDevice();
      hlc_rt_indices_.CopyToDevice();
      gas_jac_ids_.CopyToDevice();
      aq_jac_ids_.CopyToDevice();
      solvent_jac_ids_.CopyToDevice();

      views_ = Views(
          aq_indices_,
          solvent_indices_,
          hlc_rt_indices_,
          gas_jac_ids_,
          aq_jac_ids_,
          solvent_jac_ids_,
          num_phases,
          gas_idx,
          molar_volume);
    }

    /// @brief Add G(y) into the algebraic rows of `residual`.
    void AddResidual(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& residual) const
    {
      const auto& views = views_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ViewType& residual_view,
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ConstViewType& params_view)
          {
            const micm::Real molar_volume = views.molar_volume_;
            const micm::Index gas_idx = views.gas_idx_;
            for (micm::Index i_phase = 0; i_phase < views.num_phases_; ++i_phase)
            {
              const micm::Index hlc_rt_idx = views.hlc_rt_indices_[i_phase];
              const micm::Index aq_idx = views.aq_indices_[i_phase];
              const micm::Index sol_idx = views.solvent_indices_[i_phase];
              residual_view.ForEachRowStrict(
                  [molar_volume](
                      const micm::Real& hlc_rt,
                      const micm::Real& gas,
                      const micm::Real& aq,
                      const micm::Real& sol,
                      micm::Real& r) { r = hlc_rt * (sol * molar_volume) * gas - aq; },
                  params_view.GetConstColumnView(hlc_rt_idx),
                  state_view.GetConstColumnView(gas_idx),
                  state_view.GetConstColumnView(aq_idx),
                  state_view.GetConstColumnView(sol_idx),
                  residual_view.GetColumnView(aq_idx));
            }
          },
          residual,
          state_variables,
          state_parameters)(residual, state_variables, state_parameters);
    }

    /// @brief Subtract dG/dy from `jacobian`.
    void SubtractJacobian(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        SparseMatrixPolicy& jacobian) const
    {
      const auto& views = views_;
      SparseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ConstViewType& params_view,
              const typename SparseMatrixPolicy::ViewType& jac_view)
          {
            const micm::Real molar_volume = views.molar_volume_;
            const micm::Index gas_idx = views.gas_idx_;
            for (micm::Index i_phase = 0; i_phase < views.num_phases_; ++i_phase)
            {
              const micm::Index hlc_rt_idx = views.hlc_rt_indices_[i_phase];
              const micm::Index sol_idx = views.solvent_indices_[i_phase];
              const micm::Index gas_jac_id = views.gas_jac_ids_[i_phase];
              const micm::Index aq_jac_id = views.aq_jac_ids_[i_phase];
              const micm::Index solvent_jac_id = views.solvent_jac_ids_[i_phase];
              auto bv_gas = jac_view.GetBlockView(gas_jac_id);
              auto bv_aq = jac_view.GetBlockView(aq_jac_id);
              auto bv_sol = jac_view.GetBlockView(solvent_jac_id);
              jac_view.ForEachBlockStrict(
                  [molar_volume](
                      const micm::Real& hlc_rt,
                      const micm::Real& gas,
                      const micm::Real& sol,
                      micm::Real& j_gas,
                      micm::Real& j_aq,
                      micm::Real& j_sol)
                  {
                    const micm::Real f_v = sol * molar_volume;
                    j_gas -= hlc_rt * f_v;
                    j_aq -= (-1.0);
                    j_sol -= hlc_rt * molar_volume * gas;
                  },
                  params_view.GetConstColumnView(hlc_rt_idx),
                  state_view.GetConstColumnView(gas_idx),
                  state_view.GetConstColumnView(sol_idx),
                  bv_gas,
                  bv_aq,
                  bv_sol);
            }
          },
          state_variables,
          state_parameters,
          jacobian)(state_variables, state_parameters, jacobian);
    }

   private:
    Vector<micm::Index> aq_indices_{};
    Vector<micm::Index> solvent_indices_{};
    Vector<micm::Index> hlc_rt_indices_{};
    Vector<micm::Index> gas_jac_ids_{};
    Vector<micm::Index> aq_jac_ids_{};
    Vector<micm::Index> solvent_jac_ids_{};
    Views views_{};
  };
}  // namespace miam
