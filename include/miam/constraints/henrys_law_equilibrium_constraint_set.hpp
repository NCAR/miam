// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/henrys_law_equilibrium_constraint.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `HenrysLawEquilibriumConstraint`.
  /// @details Residual: G = HLC * R * T * f_v * [A_g] - [A_aq]
  ///          where f_v = [S] * solvent_molecular_weight / solvent_density.
  class HenrysLawEquilibriumConstraintSet
  {
   public:
    HenrysLawEquilibriumConstraintSet() = default;

    template<typename SparseMatrixPolicy>
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
      gas_idx_ = gas_it->second;
      molar_volume_ = config.solvent_molecular_weight_ / config.solvent_density_;

      auto phase_it = phase_prefixes.find(config.condensed_phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "HenrysLawEquilibriumConstraintSet: phase " + config.condensed_phase_.name_ + " not found");

      num_phases_ = phase_it->second.size();
      for (const auto& prefix : phase_it->second)
      {
        std::size_t aq_idx = state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.condensed_species_.name_);
        std::size_t solvent_idx = state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.solvent_.name_);
        aq_indices_.push_back(aq_idx);
        solvent_indices_.push_back(solvent_idx);
        hlc_rt_indices_.push_back(
            state_parameter_indices.at(prefix + "." + config.condensed_phase_.name_ + "." + config.uuid_ + ".hlc_rt"));
        gas_jac_ids_.push_back(jacobian.VectorIndex(0, aq_idx, gas_idx_));
        aq_jac_ids_.push_back(jacobian.VectorIndex(0, aq_idx, aq_idx));
        solvent_jac_ids_.push_back(jacobian.VectorIndex(0, aq_idx, solvent_idx));
      }
    }

    /// @brief Add G(y) into the algebraic rows of `residual`.
    template<typename DenseMatrixPolicy>
    void AddResidual(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& residual) const
    {
      const double molar_volume = molar_volume_;
      const std::size_t gas_idx = gas_idx_;
      for (std::size_t i_phase = 0; i_phase < num_phases_; ++i_phase)
      {
        const std::size_t hlc_rt_idx = hlc_rt_indices_[i_phase];
        const std::size_t aq_idx = aq_indices_[i_phase];
        const std::size_t sol_idx = solvent_indices_[i_phase];

        DenseMatrixPolicy::Function(
            [molar_volume, gas_idx, hlc_rt_idx, aq_idx, sol_idx](auto&& sv, auto&& sp, auto&& res)
            {
              res.ForEachRow(
                  [molar_volume](const double& hlc_rt, const double& gas, const double& aq, const double& sol, double& r)
                  { r = hlc_rt * (sol * molar_volume) * gas - aq; },
                  sp.GetConstColumnView(hlc_rt_idx),
                  sv.GetConstColumnView(gas_idx),
                  sv.GetConstColumnView(aq_idx),
                  sv.GetConstColumnView(sol_idx),
                  res.GetColumnView(aq_idx));
            },
            state_variables,
            state_parameters,
            residual)(state_variables, state_parameters, residual);
      }
    }

    /// @brief Subtract dG/dy from `jacobian`.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractJacobian(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        SparseMatrixPolicy& jacobian) const
    {
      const double molar_volume = molar_volume_;
      const std::size_t gas_idx = gas_idx_;
      for (std::size_t i_phase = 0; i_phase < num_phases_; ++i_phase)
      {
        const std::size_t hlc_rt_idx = hlc_rt_indices_[i_phase];
        const std::size_t sol_idx = solvent_indices_[i_phase];
        const std::size_t gas_jac_id = gas_jac_ids_[i_phase];
        const std::size_t aq_jac_id = aq_jac_ids_[i_phase];
        const std::size_t solvent_jac_id = solvent_jac_ids_[i_phase];

        SparseMatrixPolicy::Function(
            [molar_volume, gas_idx, hlc_rt_idx, sol_idx, gas_jac_id, aq_jac_id, solvent_jac_id](
                auto&& sv, auto&& sp, auto&& jac)
            {
              auto bv_gas = jac.GetBlockView(gas_jac_id);
              auto bv_aq = jac.GetBlockView(aq_jac_id);
              auto bv_sol = jac.GetBlockView(solvent_jac_id);
              jac.ForEachBlock(
                  [molar_volume](
                      const double& hlc_rt, const double& gas, const double& sol,
                      double& j_gas, double& j_aq, double& j_sol)
                  {
                    double f_v = sol * molar_volume;
                    j_gas -= hlc_rt * f_v;
                    j_aq -= (-1.0);
                    j_sol -= hlc_rt * molar_volume * gas;
                  },
                  sp.GetConstColumnView(hlc_rt_idx),
                  sv.GetConstColumnView(gas_idx),
                  sv.GetConstColumnView(sol_idx),
                  bv_gas,
                  bv_aq,
                  bv_sol);
            },
            state_variables,
            state_parameters,
            jacobian)(state_variables, state_parameters, jacobian);
      }
    }

   private:
    std::size_t num_phases_{ 0 };
    std::size_t gas_idx_{ 0 };
    double molar_volume_{ 0.0 };
    std::vector<std::size_t> aq_indices_{};
    std::vector<std::size_t> solvent_indices_{};
    std::vector<std::size_t> hlc_rt_indices_{};
    std::vector<std::size_t> gas_jac_ids_{};
    std::vector<std::size_t> aq_jac_ids_{};
    std::vector<std::size_t> solvent_jac_ids_{};
  };
}  // namespace miam
