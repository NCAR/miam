// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/condensation_rate.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace miam
{
  namespace constraint
  {
    /// @brief A Henry's Law equilibrium constraint
    /// @details Replaces the ODE row for the condensed-phase species with the
    ///          steady-state Henry's Law equilibrium condition:
    ///
    ///          \f$ G = \text{HLC} \cdot R \cdot T \cdot f_v \cdot [A_g] - [A_{aq}] = 0 \f$
    ///
    ///          where \f$ f_v = [S] \cdot M_{w,S} / \rho_S \f$ is the volume fraction of liquid
    ///          water, HLC is the Henry's Law constant [mol m⁻³ Pa⁻¹], and the algebraic variable
    ///          is always the condensed-phase species (gas disallowed to avoid overconstrained
    ///          systems with multiple phase instances sharing the same gas species).
    class HenryLawEquilibriumConstraint
    {
     public:
      std::function<double(const micm::Conditions& conditions)> henry_law_constant_;  ///< HLC(T) function [mol m⁻³ Pa⁻¹]
      micm::Species gas_species_;                                                     ///< Gas-phase species
      micm::Species condensed_species_;                                               ///< Condensed-phase solute species
      micm::Species solvent_;                                                         ///< Condensed-phase solvent species
      micm::Phase condensed_phase_;                                                   ///< The condensed phase
      double Mw_solvent_;   ///< Solvent molecular weight [kg mol⁻¹]
      double rho_solvent_;  ///< Solvent density [kg m⁻³]
      std::string uuid_;    ///< Unique identifier

      /// @brief Shared mutable HLC*R*T values, bridging UpdateConstraintParametersFunction → constraint functions.
      ///        One value per grid cell.
      std::shared_ptr<std::vector<double>> hlc_rt_values_;

      HenryLawEquilibriumConstraint() = delete;

      /// @brief Constructor
      HenryLawEquilibriumConstraint(
          std::function<double(const micm::Conditions& conditions)> henry_law_constant,
          const micm::Species& gas_species,
          const micm::Species& condensed_species,
          const micm::Species& solvent,
          const micm::Phase& condensed_phase,
          double Mw_solvent,
          double rho_solvent)
          : henry_law_constant_(henry_law_constant),
            gas_species_(gas_species),
            condensed_species_(condensed_species),
            solvent_(solvent),
            condensed_phase_(condensed_phase),
            Mw_solvent_(Mw_solvent),
            rho_solvent_(rho_solvent),
            uuid_(miam::util::generate_uuid_v4()),
            hlc_rt_values_(std::make_shared<std::vector<double>>())
      {
      }

      /// @brief Create a copy with a new UUID
      HenryLawEquilibriumConstraint CopyWithNewUuid() const
      {
        return HenryLawEquilibriumConstraint(
            henry_law_constant_, gas_species_, condensed_species_, solvent_, condensed_phase_, Mw_solvent_, rho_solvent_);
      }

      /// @brief Returns the names of algebraic variables (one per phase instance)
      /// @details The algebraic variable is always the condensed-phase species.
      std::set<std::string> ConstraintAlgebraicVariableNames(
          const std::map<std::string, std::set<std::string>>& phase_prefixes) const
      {
        std::set<std::string> names;
        auto phase_it = phase_prefixes.find(condensed_phase_.name_);
        if (phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : phase_it->second)
          {
            names.insert(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
          }
        }
        return names;
      }

      /// @brief Returns all species the constraint depends on
      std::set<std::string> ConstraintSpeciesDependencies(
          const std::map<std::string, std::set<std::string>>& phase_prefixes) const
      {
        std::set<std::string> species_names;
        // Gas species is standalone (no phase prefix)
        species_names.insert(gas_species_.name_);
        auto phase_it = phase_prefixes.find(condensed_phase_.name_);
        if (phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : phase_it->second)
          {
            species_names.insert(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
            species_names.insert(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
          }
        }
        return species_names;
      }

      /// @brief Returns non-zero constraint Jacobian element positions
      /// @details For each phase instance, the algebraic row (condensed species) depends on:
      ///          the gas species, the condensed species, and the solvent.
      std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        std::set<std::pair<std::size_t, std::size_t>> elements;
        auto gas_it = state_variable_indices.find(gas_species_.name_);
        if (gas_it == state_variable_indices.end())
          throw std::runtime_error(
              "Internal Error: Gas species " + gas_species_.name_ + " not found in state_variable_indices");
        std::size_t gas_idx = gas_it->second;

        auto phase_it = phase_prefixes.find(condensed_phase_.name_);
        if (phase_it == phase_prefixes.end())
          return elements;

        for (const auto& prefix : phase_it->second)
        {
          std::size_t aq_idx =
              state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
          std::size_t solvent_idx =
              state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);

          // dG/d[A_g], dG/d[A_aq], dG/d[S]
          elements.insert({ aq_idx, gas_idx });
          elements.insert({ aq_idx, aq_idx });
          elements.insert({ aq_idx, solvent_idx });
        }
        return elements;
      }

      /// @brief Returns a function that updates HLC*R*T for each grid cell (called during UpdateStateParameters)
      template<typename DenseMatrixPolicy>
      std::function<void(const std::vector<micm::Conditions>&)> UpdateConstraintParametersFunction() const
      {
        auto hlc_rt_values = hlc_rt_values_;
        auto hlc_fn = henry_law_constant_;
        return [hlc_rt_values, hlc_fn](const std::vector<micm::Conditions>& conditions)
        {
          hlc_rt_values->resize(conditions.size());
          for (std::size_t i = 0; i < conditions.size(); ++i)
          {
            (*hlc_rt_values)[i] = hlc_fn(conditions[i]) * util::R_gas * conditions[i].temperature_;
          }
        };
      }

      /// @brief Returns a function that computes constraint residuals G(y) = 0
      /// @details G = HLC * R * T * f_v * [A_g] - [A_aq]
      ///          where f_v = [S] * Mw_solvent / rho_solvent
      template<typename DenseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
        auto hlc_rt_values = hlc_rt_values_;
        double Mw_rho = Mw_solvent_ / rho_solvent_;

        DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };
        std::vector<double> dummy_hlc;

        auto inner = DenseMatrixPolicy::Function(
            [indices, Mw_rho](auto&& hlc_rt_vec, auto&& state_variables, auto&& residual)
            {
              for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
              {
                // G = HLC*R*T * f_v * [A_g] - [A_aq],  where f_v = [S] * Mw/rho
                residual.ForEachRow(
                    [Mw_rho](const double& hlc_rt, const double& gas, const double& aq, const double& sol, double& res)
                    { res = hlc_rt * (sol * Mw_rho) * gas - aq; },
                    hlc_rt_vec,
                    state_variables.GetConstColumnView(indices.gas_idx_),
                    state_variables.GetConstColumnView(indices.aq_indices_[i_phase]),
                    state_variables.GetConstColumnView(indices.solvent_indices_[i_phase]),
                    residual.GetColumnView(indices.aq_indices_[i_phase]));
              }
            },
            dummy_hlc, dummy_state, dummy_state);

        return [inner = std::move(inner), hlc_rt_values](
                   const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& residual) mutable
        { inner(*hlc_rt_values, state_variables, residual); };
      }

      /// @brief Returns a function that computes constraint Jacobian entries (subtracts dG/dy)
      /// @details Follows MICM convention: jac[row][col] -= dG/dy
      ///          dG/d[A_g] = HLC * R * T * f_v
      ///          dG/d[A_aq] = -1
      ///          dG/d[S] = HLC * R * T * Mw/rho * [A_g]
      template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices,
          const SparseMatrixPolicy& jacobian) const
      {
        auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
        auto hlc_rt_values = hlc_rt_values_;
        double Mw_rho = Mw_solvent_ / rho_solvent_;

        // Pre-compute block-0 VectorIndex offsets per instance
        struct PerInstanceJacData
        {
          std::size_t gas_vec;
          std::size_t aq_vec;
          std::size_t solvent_vec;
        };
        std::vector<PerInstanceJacData> jac_data(indices.number_of_phase_instances_);
        for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
        {
          std::size_t aq_row = indices.aq_indices_[i_phase];
          jac_data[i_phase].gas_vec = jacobian.VectorIndex(0, aq_row, indices.gas_idx_);
          jac_data[i_phase].aq_vec = jacobian.VectorIndex(0, aq_row, aq_row);
          jac_data[i_phase].solvent_vec = jacobian.VectorIndex(0, aq_row, indices.solvent_indices_[i_phase]);
        }

        DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };
        std::vector<double> dummy_hlc;

        auto inner = SparseMatrixPolicy::Function(
            [indices, jac_data, Mw_rho](auto&& hlc_rt_vec, auto&& state_variables, auto&& jacobian_values)
            {
              for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
              {
                const auto& jd = jac_data[i_phase];

                auto bv_gas = jacobian_values.GetBlockView(jd.gas_vec);
                auto bv_aq = jacobian_values.GetBlockView(jd.aq_vec);
                auto bv_sol = jacobian_values.GetBlockView(jd.solvent_vec);

                // jac -= dG/d[A_g] = HLC*R*T * f_v
                // jac -= dG/d[A_aq] = -1
                // jac -= dG/d[S] = HLC*R*T * Mw/rho * [A_g]
                jacobian_values.ForEachBlock(
                    [Mw_rho](const double& hlc_rt, const double& gas, const double& sol,
                             double& j_gas, double& j_aq, double& j_sol)
                    {
                      double f_v = sol * Mw_rho;
                      j_gas -= hlc_rt * f_v;
                      j_aq -= (-1.0);
                      j_sol -= hlc_rt * Mw_rho * gas;
                    },
                    hlc_rt_vec,
                    state_variables.GetConstColumnView(indices.gas_idx_),
                    state_variables.GetConstColumnView(indices.solvent_indices_[i_phase]),
                    bv_gas, bv_aq, bv_sol);
              }
            },
            dummy_hlc, dummy_state, jacobian);

        return [inner = std::move(inner), hlc_rt_values](
                   const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian_values) mutable
        { inner(*hlc_rt_values, state_variables, jacobian_values); };
      }

     private:
      /// @brief Helper struct for state variable indices across phase instances
      struct StateVariableIndices
      {
        std::size_t number_of_phase_instances_;
        std::size_t gas_idx_;                      ///< Gas species index (shared across instances)
        std::vector<std::size_t> aq_indices_;      ///< Condensed species index per instance
        std::vector<std::size_t> solvent_indices_;  ///< Solvent index per instance
      };

      /// @brief Helper struct for Jacobian sparse matrix indices
      struct JacobianIndices
      {
        std::vector<std::size_t> gas_jac_indices_;      ///< [aq_row, gas_col] per instance (block 0)
        std::vector<std::size_t> aq_jac_indices_;       ///< [aq_row, aq_col] per instance (block 0)
        std::vector<std::size_t> solvent_jac_indices_;   ///< [aq_row, solvent_col] per instance (block 0)
        std::size_t block_stride_;                       ///< Flat vector stride between blocks
      };

      /// @brief Build state variable indices for all phase instances
      StateVariableIndices GetStateVariableIndices(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        StateVariableIndices indices;
        auto gas_it = state_variable_indices.find(gas_species_.name_);
        if (gas_it == state_variable_indices.end())
          throw std::runtime_error(
              "Internal Error: HenryLawEquilibriumConstraint: Gas species " + gas_species_.name_ +
              " not found in state_variable_indices");
        indices.gas_idx_ = gas_it->second;

        auto phase_it = phase_prefixes.find(condensed_phase_.name_);
        if (phase_it == phase_prefixes.end())
        {
          throw std::runtime_error(
              "Internal Error: HenryLawEquilibriumConstraint: Phase " + condensed_phase_.name_ +
              " not found in phase_prefixes");
        }
        const auto& prefixes = phase_it->second;
        indices.number_of_phase_instances_ = prefixes.size();
        indices.aq_indices_.resize(prefixes.size());
        indices.solvent_indices_.resize(prefixes.size());

        std::size_t i_phase = 0;
        for (const auto& prefix : prefixes)
        {
          indices.aq_indices_[i_phase] =
              state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
          indices.solvent_indices_[i_phase] =
              state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
          ++i_phase;
        }
        return indices;
      }

      /// @brief Build Jacobian sparse matrix indices for all phase instances
      JacobianIndices GetJacobianIndices(
          const StateVariableIndices& var_indices,
          const auto& jacobian) const
      {
        JacobianIndices jac_indices;
        std::size_t num_blocks = jacobian.NumberOfBlocks();
        std::size_t block_stride = jacobian.FlatBlockSize();
        jac_indices.gas_jac_indices_.resize(var_indices.number_of_phase_instances_);
        jac_indices.aq_jac_indices_.resize(var_indices.number_of_phase_instances_);
        jac_indices.solvent_jac_indices_.resize(var_indices.number_of_phase_instances_);
        jac_indices.block_stride_ = block_stride;

        for (std::size_t i_phase = 0; i_phase < var_indices.number_of_phase_instances_; ++i_phase)
        {
          std::size_t aq_row = var_indices.aq_indices_[i_phase];
          jac_indices.gas_jac_indices_[i_phase] = jacobian.VectorIndex(0, aq_row, var_indices.gas_idx_);
          jac_indices.aq_jac_indices_[i_phase] = jacobian.VectorIndex(0, aq_row, aq_row);
          jac_indices.solvent_jac_indices_[i_phase] =
              jacobian.VectorIndex(0, aq_row, var_indices.solvent_indices_[i_phase]);
        }
        return jac_indices;
      }
    };
  }  // namespace constraint
}  // namespace miam
