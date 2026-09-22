// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/math/condensation_rate.hpp>
#include <miam/processes/constants/rate_expression.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>

#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
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
  class HenrysLawEquilibriumConstraint
  {
   public:
    HenrysLawConstantExpression henrys_law_constant_;   ///< HLC(T) expression [mol m⁻³ Pa⁻¹]
    micm::Species gas_species_;                         ///< Gas-phase species
    micm::Species condensed_species_;                   ///< Condensed-phase solute species
    micm::Species solvent_;                             ///< Condensed-phase solvent species
    micm::Phase condensed_phase_;                       ///< The condensed phase
    double solvent_molecular_weight_;  ///< Solvent molecular weight [kg mol⁻¹]
    double solvent_density_;           ///< Solvent density [kg m⁻³]
    std::string uuid_;                 ///< Unique identifier

    HenrysLawEquilibriumConstraint() = delete;

    /// @brief Constructor
    HenrysLawEquilibriumConstraint(
        HenrysLawConstantExpression henrys_law_constant,
        const micm::Species& gas_species,
        const micm::Species& condensed_species,
        const micm::Species& solvent,
        const micm::Phase& condensed_phase,
        double solvent_molecular_weight,
        double solvent_density)
        : henrys_law_constant_(std::move(henrys_law_constant)),
          gas_species_(gas_species),
          condensed_species_(condensed_species),
          solvent_(solvent),
          condensed_phase_(condensed_phase),
          solvent_molecular_weight_(solvent_molecular_weight),
          solvent_density_(solvent_density),
          uuid_(GenerateUuid())
    {
    }

    /// @brief Create a copy with a new UUID
    HenrysLawEquilibriumConstraint CopyWithNewUuid() const
    {
      return HenrysLawEquilibriumConstraint(
          henrys_law_constant_,
          gas_species_,
          condensed_species_,
          solvent_,
          condensed_phase_,
          solvent_molecular_weight_,
          solvent_density_);
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
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_VARIABLE,
            "Internal Error: Gas species " + gas_species_.name_ + " not found in state_variable_indices");
      std::size_t gas_idx = gas_it->second;

      auto phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (phase_it == phase_prefixes.end())
        return elements;

      for (const auto& prefix : phase_it->second)
      {
        std::size_t aq_idx =
            state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
        std::size_t solvent_idx = state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);

        // dG/d[A_g], dG/d[A_aq], dG/d[S]
        elements.insert({ aq_idx, gas_idx });
        elements.insert({ aq_idx, aq_idx });
        elements.insert({ aq_idx, solvent_idx });
      }
      return elements;
    }

    /// @brief Returns the names of state parameters owned by this constraint (one per phase instance).
    /// @details Each phase instance writes \f$ \text{HLC}(T) \cdot R \cdot T \f$ to a dedicated
    ///          column of the state parameter matrix every time conditions change.
    std::set<std::string> ConstraintStateParameterNames(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      auto phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          names.insert(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc_rt");
      }
      return names;
    }

    /// @brief Returns a function that writes \f$ \text{HLC}(T) \cdot R \cdot T \f$ per grid cell
    ///        into the state parameter matrix.
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateConstraintParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
    {
      std::vector<std::size_t> hlc_rt_indices;
      auto phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          hlc_rt_indices.push_back(
              state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc_rt"));
      }
      auto hlc_expr = henrys_law_constant_;

      DenseMatrixPolicy state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> conditions_vector;

      return DenseMatrixPolicy::Function(
          [hlc_rt_indices, hlc_expr](auto&& conditions, auto&& params)
          {
            for (const auto& hlc_rt_idx : hlc_rt_indices)
              params.ForEachRowStrict(
                  [hlc_expr](const micm::Conditions& cond, double& hlc_rt)
                  { hlc_rt = EvaluateExpression(hlc_expr, cond) * micm::constants::GAS_CONSTANT * cond.temperature_; },
                  conditions,
                  params.GetColumnView(hlc_rt_idx));
          },
          conditions_vector,
          state_parameters);
    }

   private:
    /// @brief Helper struct for state variable indices across phase instances
    struct StateVariableIndices
    {
      std::size_t number_of_phase_instances_;
      std::size_t gas_idx_;                       ///< Gas species index (shared across instances)
      std::vector<std::size_t> aq_indices_;       ///< Condensed species index per instance
      std::vector<std::size_t> solvent_indices_;  ///< Solvent index per instance
    };

    /// @brief Helper struct for Jacobian sparse matrix indices
    struct JacobianIndices
    {
      std::vector<std::size_t> gas_jac_indices_;      ///< [aq_row, gas_col] per instance (block 0)
      std::vector<std::size_t> aq_jac_indices_;       ///< [aq_row, aq_col] per instance (block 0)
      std::vector<std::size_t> solvent_jac_indices_;  ///< [aq_row, solvent_col] per instance (block 0)
      std::size_t block_stride_;                      ///< Flat vector stride between blocks
    };

    /// @brief Build state variable indices for all phase instances
    StateVariableIndices GetStateVariableIndices(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      StateVariableIndices indices;
      auto gas_it = state_variable_indices.find(gas_species_.name_);
      if (gas_it == state_variable_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_VARIABLE,
            "HenrysLawEquilibriumConstraint: Gas species " + gas_species_.name_ + " not found in state_variable_indices");
      indices.gas_idx_ = gas_it->second;

      auto phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (phase_it == phase_prefixes.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "HenrysLawEquilibriumConstraint: Phase " + condensed_phase_.name_ + " not found in phase_prefixes");
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
    JacobianIndices GetJacobianIndices(const StateVariableIndices& var_indices, const auto& jacobian) const
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
        jac_indices.solvent_jac_indices_[i_phase] = jacobian.VectorIndex(0, aq_row, var_indices.solvent_indices_[i_phase]);
      }
      return jac_indices;
    }
  };
}  // namespace miam
