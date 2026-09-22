// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/math/condensation_rate.hpp>
#include <miam/processes/constants/rate_expression.hpp>
#include <miam/representations/aerosol_property.hpp>
#include <miam/representations/aerosol_property_descriptor.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  /// @brief Henry's Law phase transfer process
  /// @details Represents the transfer of a gas-phase species into a condensed phase
  ///          (and its re-evaporation) governed by Henry's Law equilibrium. The net rate is:
  ///
  ///          d[A]_gas/dt  = -φ_p · k_cond · [A]_gas + φ_p · k_evap · [A]_aq / f_v
  ///          d[A]_aq/dt   = +φ_p · k_cond · [A]_gas - φ_p · k_evap · [A]_aq / f_v
  ///
  ///          where k_evap = k_cond / (HLC · R · T), f_v = [solvent] · solvent_molecular_weight / solvent_density  [m³
  ///          mol⁻¹], and φ_p is the phase volume fraction.
  class HenrysLawPhaseTransfer
  {
   public:
    HenrysLawConstantExpression henrys_law_constant_;   ///< HLC(T) expression [mol m⁻³ Pa⁻¹]
    micm::Species gas_species_;                         ///< Gas-phase species
    micm::Species condensed_species_;                   ///< Condensed-phase solute species
    micm::Species solvent_;                             ///< Condensed-phase solvent species
    micm::Phase condensed_phase_;                       ///< The condensed phase
    double diffusion_coefficient_;      ///< Gas-phase diffusion coefficient [m² s⁻¹]
    double accommodation_coefficient_;  ///< Mass accommodation coefficient [dimensionless]
    double gas_molecular_weight_;       ///< Gas-phase molecular weight [kg mol⁻¹]
    double solvent_molecular_weight_;   ///< Solvent molecular weight [kg mol⁻¹]
    double solvent_density_;            ///< Solvent density [kg m⁻³]
    std::string uuid_;                  ///< Unique identifier

    HenrysLawPhaseTransfer() = delete;

    /// @brief Constructor
    HenrysLawPhaseTransfer(
        HenrysLawConstantExpression henrys_law_constant,
        const micm::Species& gas_species,
        const micm::Species& condensed_species,
        const micm::Species& solvent,
        const micm::Phase& condensed_phase,
        double diffusion_coefficient,
        double accommodation_coefficient,
        double gas_molecular_weight,
        double solvent_molecular_weight,
        double solvent_density)
        : henrys_law_constant_(std::move(henrys_law_constant)),
          gas_species_(gas_species),
          condensed_species_(condensed_species),
          solvent_(solvent),
          condensed_phase_(condensed_phase),
          diffusion_coefficient_(diffusion_coefficient),
          accommodation_coefficient_(accommodation_coefficient),
          gas_molecular_weight_(gas_molecular_weight),
          solvent_molecular_weight_(solvent_molecular_weight),
          solvent_density_(solvent_density),
          uuid_(GenerateUuid())
    {
    }

    /// @brief Create a copy with a new UUID
    HenrysLawPhaseTransfer CopyWithNewUuid() const
    {
      return HenrysLawPhaseTransfer(
          henrys_law_constant_,
          gas_species_,
          condensed_species_,
          solvent_,
          condensed_phase_,
          diffusion_coefficient_,
          accommodation_coefficient_,
          gas_molecular_weight_,
          solvent_molecular_weight_,
          solvent_density_);
    }

    /// @brief Returns unique parameter names for this process
    std::set<std::string> ProcessParameterNames(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      auto it = phase_prefixes.find(condensed_phase_.name_);
      if (it != phase_prefixes.end())
      {
        for (const auto& prefix : it->second)
        {
          names.insert(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
          names.insert(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
        }
      }
      return names;
    }

    /// @brief Returns participating species' unique state names
    std::set<std::string> SpeciesUsed(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> species_names;
      // Gas species is a standalone state variable
      species_names.insert(gas_species_.name_);
      // Condensed-phase species are per instance
      auto it = phase_prefixes.find(condensed_phase_.name_);
      if (it != phase_prefixes.end())
      {
        for (const auto& prefix : it->second)
        {
          species_names.insert(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
          species_names.insert(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
        }
      }
      else
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: Phase " + condensed_phase_.name_ + " not found in phase_prefixes for process " + uuid_);
      }
      return species_names;
    }

    /// @brief Returns the aerosol properties required by this process
    std::map<std::string, std::vector<AerosolProperty>> RequiredAerosolProperties() const
    {
      return {
        { condensed_phase_.name_,
          { AerosolProperty::EffectiveRadius, AerosolProperty::NumberConcentration, AerosolProperty::PhaseVolumeFraction } }
      };
    }

    /// @brief Returns non-zero Jacobian element positions
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
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
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: Phase " + condensed_phase_.name_ + " not found in phase_prefixes for process " + uuid_);
      }

      // We need provider-dependent indices, but at this stage we don't have providers yet.
      // Conservatively include all variables in each representation prefix as potential
      // indirect dependencies (through EffectiveRadius, NumberConcentration, PhaseVolumeFraction).
      for (const auto& prefix : phase_it->second)
      {
        std::size_t aq_idx =
            state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
        std::size_t solvent_idx = state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);

        // Direct dependencies
        elements.insert({ gas_idx, gas_idx });
        elements.insert({ gas_idx, aq_idx });
        elements.insert({ gas_idx, solvent_idx });
        elements.insert({ aq_idx, gas_idx });
        elements.insert({ aq_idx, aq_idx });
        elements.insert({ aq_idx, solvent_idx });

        // Indirect dependencies: any variable under this prefix may affect aerosol properties
        std::string prefix_dot = prefix + ".";
        for (const auto& [var_name, var_idx] : state_variable_indices)
        {
          if (var_name.substr(0, prefix_dot.size()) == prefix_dot)
          {
            elements.insert({ gas_idx, var_idx });
            elements.insert({ aq_idx, var_idx });
          }
        }
      }
      return elements;
    }

    /// @brief Returns non-zero Jacobian elements including indirect dependencies through aerosol descriptors.
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const auto& descriptors) const
    {
      auto elements = NonZeroJacobianElements(phase_prefixes, state_variable_indices);
      auto gas_idx = state_variable_indices.at(gas_species_.name_);

      for (const auto& [prefix, desc_map] : descriptors)
      {
        std::size_t aq_idx =
            state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);

        for (const auto& [prop, descriptor] : desc_map)
        {
          for (std::size_t var_j : DependentVariableIndices(descriptor))
          {
            elements.insert({ gas_idx, var_j });
            elements.insert({ aq_idx, var_j });
          }
        }
      }
      return elements;
    }

    /// @brief Returns a function that updates state parameters (HLC and temperature)
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateStateParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices) const
    {
      std::vector<std::size_t> hlc_indices;
      std::vector<std::size_t> temp_indices;
      auto it = phase_prefixes.find(condensed_phase_.name_);
      if (it != phase_prefixes.end())
      {
        for (const auto& prefix : it->second)
        {
          std::string hlc_param = prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc";
          std::string temp_param = prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature";
          if (state_parameter_indices.find(hlc_param) == state_parameter_indices.end())
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_MISSING_STATE_PARAMETER,
                "Internal Error: HLC parameter " + hlc_param + " not found");
          if (state_parameter_indices.find(temp_param) == state_parameter_indices.end())
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_MISSING_STATE_PARAMETER,
                "Internal Error: Temperature parameter " + temp_param + " not found");
          hlc_indices.push_back(state_parameter_indices.at(hlc_param));
          temp_indices.push_back(state_parameter_indices.at(temp_param));
        }
      }

      DenseMatrixPolicy dummy{ 1, state_parameter_indices.size(), 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> dummy_conditions;

      return DenseMatrixPolicy::Function(
          [this, hlc_indices, temp_indices](auto&& conditions, auto&& params)
          {
            for (std::size_t i = 0; i < hlc_indices.size(); ++i)
            {
              params.ForEachRowStrict(
                  [&](const micm::Conditions& cond, double& hlc, double& T)
                  {
                    hlc = EvaluateExpression(henrys_law_constant_, cond);
                    T = cond.temperature_;
                  },
                  conditions,
                  params.GetColumnView(hlc_indices[i]),
                  params.GetColumnView(temp_indices[i]));
            }
          },
          dummy_conditions,
          dummy);
    }

  };
}  // namespace miam
