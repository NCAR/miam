// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/math/condensation_rate.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <miam/representations/aerosol_property.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/types.hpp>

#include <algorithm>
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
    HenrysLawConstant henrys_law_constant_;  ///< HLC(T) function [mol m⁻³ Pa⁻¹]
    micm::Species gas_species_;                                                      ///< Gas-phase species
    micm::Species condensed_species_;                                                ///< Condensed-phase solute species
    micm::Species solvent_;                                                          ///< Condensed-phase solvent species
    micm::Phase condensed_phase_;                                                    ///< The condensed phase
    double diffusion_coefficient_;      ///< Gas-phase diffusion coefficient [m² s⁻¹]
    double accommodation_coefficient_;  ///< Mass accommodation coefficient [dimensionless]
    double gas_molecular_weight_;       ///< Gas-phase molecular weight [kg mol⁻¹]
    double solvent_molecular_weight_;   ///< Solvent molecular weight [kg mol⁻¹]
    double solvent_density_;            ///< Solvent density [kg m⁻³]
    std::string uuid_;                  ///< Unique identifier

    HenrysLawPhaseTransfer() = delete;

    /// @brief Constructor
    HenrysLawPhaseTransfer(
        HenrysLawConstant henrys_law_constant,
        const micm::Species& gas_species,
        const micm::Species& condensed_species,
        const micm::Species& solvent,
        const micm::Phase& condensed_phase,
        double diffusion_coefficient,
        double accommodation_coefficient,
        double gas_molecular_weight,
        double solvent_molecular_weight,
        double solvent_density)
        : henrys_law_constant_(henrys_law_constant),
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

    /// @brief Returns non-zero Jacobian elements (common interface overload with providers)
    template<typename DenseMatrixPolicy>
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>>& providers) const
    {
      auto elements = NonZeroJacobianElements(phase_prefixes, state_variable_indices);
      auto gas_idx = state_variable_indices.at(gas_species_.name_);

      // Add indirect dependencies through aerosol property providers
      for (const auto& [prefix, prov_map] : providers)
      {
        std::size_t aq_idx =
            state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);

        for (const auto& [prop, provider] : prov_map)
        {
          for (std::size_t var_j : provider.dependent_variable_indices)
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
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
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

      const HenrysLawConstant henrys_law_constant = henrys_law_constant_;
      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      struct Storage
      {
        Vector hlc_indices, temp_indices;
      };
      auto storage = std::make_shared<Storage>(Storage{ Vector(hlc_indices), Vector(temp_indices) });
      storage->hlc_indices.CopyToDevice();
      storage->temp_indices.CopyToDevice();
      auto hlc_view = storage->hlc_indices.GetView();
      auto temp_view = storage->temp_indices.GetView();
      const std::size_t num_instances = hlc_indices.size();
      DenseMatrixPolicy dummy{ 1, state_parameter_indices.size(), 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> dummy_conditions;

      auto function = DenseMatrixPolicy::Function(
              MICM_LAMBDA(
                  const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
                  const typename DenseMatrixPolicy::ViewType& params_view)
          {
                for (std::size_t i = 0; i < num_instances; ++i)
                params_view.ForEachRowStrict(
                    [henrys_law_constant](const micm::Conditions& cond, micm::Real& hlc, micm::Real& T)
            {
                      hlc = Calculate(henrys_law_constant, cond);
                    T = cond.temperature_;
                  },
                    conditions_view,
                      params_view.GetColumnView(hlc_view[i]),
                      params_view.GetColumnView(temp_view[i]));
          },
          dummy_conditions,
          dummy);

      return [storage, function](
                 const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
                 DenseMatrixPolicy& params) mutable { function(conditions, params); };
    }

    /// @brief Returns a function that calculates the forcing terms (common interface with providers)
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> providers) const
    {
      auto gas_idx = state_variable_indices.at(gas_species_.name_);

      struct InstanceData
      {
        std::size_t aq_species_idx;
        std::size_t solvent_species_idx;
        std::size_t hlc_param_idx;
        std::size_t temperature_param_idx;
        double molar_volume;  ///< Solvent molar volume [m³ mol⁻¹] = solvent_molecular_weight / solvent_density
        AerosolPropertyProvider<DenseMatrixPolicy> r_eff_provider;
        AerosolPropertyProvider<DenseMatrixPolicy> N_provider;
        AerosolPropertyProvider<DenseMatrixPolicy> phi_provider;
        CondensationRateProvider cond_rate_provider;
      };

      std::vector<InstanceData> instances;
      auto my_phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (my_phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : my_phase_it->second)
        {
          auto prov_it = providers.find(prefix);
          if (prov_it == providers.end())
            continue;
          const auto& prov_map = prov_it->second;
          InstanceData inst;
          inst.aq_species_idx =
              state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
          inst.solvent_species_idx = state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
          inst.hlc_param_idx = state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
          inst.temperature_param_idx =
              state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
          inst.molar_volume = solvent_molecular_weight_ / solvent_density_;
          inst.r_eff_provider = prov_map.at(AerosolProperty::EffectiveRadius);
          inst.N_provider = prov_map.at(AerosolProperty::NumberConcentration);
          inst.phi_provider = prov_map.at(AerosolProperty::PhaseVolumeFraction);
          inst.cond_rate_provider =
              MakeCondensationRateProvider(diffusion_coefficient_, accommodation_coefficient_, gas_molecular_weight_);
          instances.push_back(std::move(inst));
        }
      }

      // Build one Function per instance at setup time
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_buf{ 1, 1, 0.0 };

      using InnerFuncType = std::function<void(
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&)>;
      std::vector<InnerFuncType> inner_functions;

      for (const auto& inst : instances)
      {
          // Copy the instance data that the kernel uses into plain values
          const std::size_t aq_idx = inst.aq_species_idx;
          const std::size_t solvent_idx = inst.solvent_species_idx;
          const std::size_t hlc_idx = inst.hlc_param_idx;
          const std::size_t temp_idx = inst.temperature_param_idx;
          const double molar_volume = inst.molar_volume;
          const CondensationRateProvider cond = inst.cond_rate_provider;

          inner_functions.push_back(DenseMatrixPolicy::Function(
              MICM_LAMBDA(
                  const typename DenseMatrixPolicy::ConstViewType& params_view,
                  const typename DenseMatrixPolicy::ConstViewType& state_view,
                  const typename DenseMatrixPolicy::ViewType& forcing_view,
                  const typename DenseMatrixPolicy::ConstViewType& r_eff_view,
                  const typename DenseMatrixPolicy::ConstViewType& N_view,
                  const typename DenseMatrixPolicy::ConstViewType& phi_view)
            {
                auto net = forcing_view.GetRowVariable();

              // Compute net transfer rate
                forcing_view.ForEachRowStrict(
                    [molar_volume, cond](
                        const micm::Real& r_eff,
                        const micm::Real& N,
                        const micm::Real& phi,
                        const micm::Real& hlc,
                        const micm::Real& T,
                        const micm::Real& gas,
                        const micm::Real& aq,
                        const micm::Real& solvent,
                        micm::Real& net_val)
                  {
                      const micm::Real kc = cond.ComputeValue(r_eff, N, T);
                      const micm::Real kc_eff = phi * kc;
                      const micm::Real ke_eff = kc_eff / (hlc * micm::constants::GAS_CONSTANT * T);
                      const micm::Real fv = solvent * molar_volume;
                    net_val = kc_eff * gas - ke_eff * aq / fv;
                  },
                  r_eff_view.GetConstColumnView(0),
                  N_view.GetConstColumnView(0),
                  phi_view.GetConstColumnView(0),
                    params_view.GetConstColumnView(hlc_idx),
                    params_view.GetConstColumnView(temp_idx),
                    state_view.GetConstColumnView(gas_idx),
                    state_view.GetConstColumnView(aq_idx),
                    state_view.GetConstColumnView(solvent_idx),
                  net);

              // Apply to gas forcing (subtract)
                forcing_view.ForEachRowStrict(
                    [](const micm::Real& net_val, micm::Real& f_gas) { f_gas -= net_val; },
                    net,
                    forcing_view.GetColumnView(gas_idx));

              // Apply to aq forcing (add)
                forcing_view.ForEachRowStrict(
                    [](const micm::Real& net_val, micm::Real& f_aq) { f_aq += net_val; },
                  net,
                    forcing_view.GetColumnView(aq_idx));
            },
            dummy_state_parameters,
            dummy_state_variables,
            dummy_state_variables,
            dummy_buf,
            dummy_buf,
              dummy_buf));
      }

      return [instances = std::move(instances), inner_functions = std::move(inner_functions)](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms)
      {
        std::size_t num_rows = state_parameters.NumRows();

        for (std::size_t i = 0; i < instances.size(); ++i)
        {
          const auto& inst = instances[i];

          // Pre-compute aerosol properties into temp buffers (full-matrix calls)
          DenseMatrixPolicy r_eff_buf{ num_rows, 1, 0.0 };
          DenseMatrixPolicy N_buf{ num_rows, 1, 0.0 };
          DenseMatrixPolicy phi_buf{ num_rows, 1, 0.0 };
          inst.r_eff_provider.ComputeValue(state_parameters, state_variables, r_eff_buf);
          inst.N_provider.ComputeValue(state_parameters, state_variables, N_buf);
          inst.phi_provider.ComputeValue(state_parameters, state_variables, phi_buf);

          // Call the Function that was built for this instance
          inner_functions[i](state_parameters, state_variables, forcing_terms, r_eff_buf, N_buf, phi_buf);
        }
      };
    }

    /// @brief Returns a function that calculates Jacobian contributions (common interface with providers)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian,
        std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> providers) const
    {
      auto gas_idx = state_variable_indices.at(gas_species_.name_);
      using Vector = typename SparseMatrixPolicy::template VectorType<std::size_t>;

      struct InstanceData
      {
        std::size_t aq_species_idx;
        std::size_t solvent_species_idx;
        std::size_t hlc_param_idx;
        std::size_t temperature_param_idx;
        double molar_volume;  ///< Solvent molar volume [m³ mol⁻¹] = solvent_molecular_weight / solvent_density
        AerosolPropertyProvider<DenseMatrixPolicy> r_eff_provider;
        AerosolPropertyProvider<DenseMatrixPolicy> N_provider;
        AerosolPropertyProvider<DenseMatrixPolicy> phi_provider;
        CondensationRateProvider cond_rate_provider;
        std::size_t n_r_eff_deps;
        std::size_t n_N_deps;
        std::size_t n_phi_deps;
        // Jacobian flat ids stored in a device-ready vector for use with GetBlockView
        // Layout: [6 direct] [2*n_r_eff_deps indirect_r_eff] [2*n_N_deps indirect_N] [2*n_phi_deps indirect_phi]
        Vector jac_indices;
      };

      std::vector<InstanceData> jac_instances;
      auto my_jac_phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (my_jac_phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : my_jac_phase_it->second)
        {
          auto prov_it = providers.find(prefix);
          if (prov_it == providers.end())
            continue;
          const auto& prov_map = prov_it->second;
          InstanceData inst;
          inst.aq_species_idx =
              state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
          inst.solvent_species_idx = state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
          inst.hlc_param_idx = state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
          inst.temperature_param_idx =
              state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
          inst.molar_volume = solvent_molecular_weight_ / solvent_density_;
          inst.r_eff_provider = prov_map.at(AerosolProperty::EffectiveRadius);
          inst.N_provider = prov_map.at(AerosolProperty::NumberConcentration);
          inst.phi_provider = prov_map.at(AerosolProperty::PhaseVolumeFraction);
          inst.cond_rate_provider =
              MakeCondensationRateProvider(diffusion_coefficient_, accommodation_coefficient_, gas_molecular_weight_);
          inst.n_r_eff_deps = prov_map.at(AerosolProperty::EffectiveRadius).dependent_variable_indices.size();
          inst.n_N_deps = prov_map.at(AerosolProperty::NumberConcentration).dependent_variable_indices.size();
          inst.n_phi_deps = prov_map.at(AerosolProperty::PhaseVolumeFraction).dependent_variable_indices.size();

          std::size_t aq_idx = inst.aq_species_idx;
          std::size_t solvent_idx = inst.solvent_species_idx;

          // Build flat Jacobian index vector
          std::size_t total_indices = 6 + 2 * inst.n_r_eff_deps + 2 * inst.n_N_deps + 2 * inst.n_phi_deps;
          std::vector<std::size_t> jac_ids(total_indices);
          std::size_t idx = 0;

          // Direct entries (6 total)
          jac_ids[idx++] = jacobian.VectorIndex(0, gas_idx, gas_idx);
          jac_ids[idx++] = jacobian.VectorIndex(0, gas_idx, aq_idx);
          jac_ids[idx++] = jacobian.VectorIndex(0, gas_idx, solvent_idx);
          jac_ids[idx++] = jacobian.VectorIndex(0, aq_idx, gas_idx);
          jac_ids[idx++] = jacobian.VectorIndex(0, aq_idx, aq_idx);
          jac_ids[idx++] = jacobian.VectorIndex(0, aq_idx, solvent_idx);

          // Indirect through r_eff
          for (std::size_t var_j : prov_map.at(AerosolProperty::EffectiveRadius).dependent_variable_indices)
          {
            jac_ids[idx++] = jacobian.VectorIndex(0, gas_idx, var_j);
            jac_ids[idx++] = jacobian.VectorIndex(0, aq_idx, var_j);
          }
          // Indirect through N
          for (std::size_t var_j : prov_map.at(AerosolProperty::NumberConcentration).dependent_variable_indices)
          {
            jac_ids[idx++] = jacobian.VectorIndex(0, gas_idx, var_j);
            jac_ids[idx++] = jacobian.VectorIndex(0, aq_idx, var_j);
          }
          // Indirect through phi
          for (std::size_t var_j : prov_map.at(AerosolProperty::PhaseVolumeFraction).dependent_variable_indices)
          {
            jac_ids[idx++] = jacobian.VectorIndex(0, gas_idx, var_j);
            jac_ids[idx++] = jacobian.VectorIndex(0, aq_idx, var_j);
          }
          inst.jac_indices = Vector(jac_ids);
          inst.jac_indices.CopyToDevice();

          jac_instances.push_back(std::move(inst));
        }
      }

      // The kernels hold views into the Jacobian index vectors, so keep the instances on the heap
      auto storage = std::make_shared<std::vector<InstanceData>>(std::move(jac_instances));

      // Build one Function per instance at setup time
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_buf{ 1, 1, 0.0 };
      DenseMatrixPolicy dummy_partials{ 1, 1, 0.0 };

      using InnerJacFuncType = std::function<void(
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          SparseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&,
          const DenseMatrixPolicy&)>;
      std::vector<InnerJacFuncType> inner_jac_functions;

      for (const auto& inst : *storage)
      {
          // Copy the instance data that the kernel uses into plain values and views
          const std::size_t aq_idx = inst.aq_species_idx;
          const std::size_t solvent_idx = inst.solvent_species_idx;
          const std::size_t hlc_idx = inst.hlc_param_idx;
          const std::size_t temp_idx = inst.temperature_param_idx;
          const double molar_volume = inst.molar_volume;
          const CondensationRateProvider cond = inst.cond_rate_provider;
          const std::size_t n_r_eff_deps = inst.n_r_eff_deps;
          const std::size_t n_N_deps = inst.n_N_deps;
          const std::size_t n_phi_deps = inst.n_phi_deps;
          auto jac_id_view = inst.jac_indices.GetView();

        inner_jac_functions.push_back(SparseMatrixPolicy::Function(
              MICM_LAMBDA(
                  const typename DenseMatrixPolicy::ConstViewType& params_view,
                  const typename DenseMatrixPolicy::ConstViewType& state_view,
                  const typename SparseMatrixPolicy::ViewType& jac_view,
                  const typename DenseMatrixPolicy::ConstViewType& r_eff_view,
                  const typename DenseMatrixPolicy::ConstViewType& N_view,
                  const typename DenseMatrixPolicy::ConstViewType& phi_view,
                  const typename DenseMatrixPolicy::ConstViewType& r_eff_partials_view,
                  const typename DenseMatrixPolicy::ConstViewType& N_partials_view,
                  const typename DenseMatrixPolicy::ConstViewType& phi_partials_view)
            {
                std::size_t idx = 0;

              // Pre-extract BlockViews sequentially to avoid unspecified argument evaluation order
                auto bv_gg = jac_view.GetBlockView(jac_id_view[idx++]);
                auto bv_ga = jac_view.GetBlockView(jac_id_view[idx++]);
                auto bv_gs = jac_view.GetBlockView(jac_id_view[idx++]);
                auto bv_ag = jac_view.GetBlockView(jac_id_view[idx++]);
                auto bv_aa = jac_view.GetBlockView(jac_id_view[idx++]);
                auto bv_as = jac_view.GetBlockView(jac_id_view[idx++]);

              // Read inputs and compute direct Jacobian entries
                jac_view.ForEachBlockStrict(
                    [molar_volume, cond](
                        const micm::Real& r_eff,
                        const micm::Real& N,
                        const micm::Real& phi,
                        const micm::Real& hlc,
                        const micm::Real& T,
                        const micm::Real& gas,
                        const micm::Real& aq,
                        const micm::Real& solvent,
                        micm::Real& j_gg,
                        micm::Real& j_ga,
                        micm::Real& j_gs,
                        micm::Real& j_ag,
                        micm::Real& j_aa,
                        micm::Real& j_as)
                  {
                      const micm::Real kc = cond.ComputeValue(r_eff, N, T);
                      const micm::Real ke = kc / (hlc * micm::constants::GAS_CONSTANT * T);
                      const micm::Real fv = solvent * molar_volume;
                    // -J[gas, gas] = +φ · k_cond
                    j_gg += phi * kc;
                    // -J[gas, aq] = -φ · k_evap / f_v
                    j_ga -= phi * ke / fv;
                    // -J[gas, solvent] = +φ · k_evap · [aq] / (f_v · [solvent])
                    j_gs += phi * ke * aq / (fv * solvent);
                    // -J[aq, gas] = -φ · k_cond
                    j_ag -= phi * kc;
                    // -J[aq, aq] = +φ · k_evap / f_v
                    j_aa += phi * ke / fv;
                    // -J[aq, solvent] = -φ · k_evap · [aq] / (f_v · [solvent])
                    j_as -= phi * ke * aq / (fv * solvent);
                  },
                  r_eff_view.GetConstColumnView(0),
                  N_view.GetConstColumnView(0),
                  phi_view.GetConstColumnView(0),
                    params_view.GetConstColumnView(hlc_idx),
                    params_view.GetConstColumnView(temp_idx),
                    state_view.GetConstColumnView(gas_idx),
                    state_view.GetConstColumnView(aq_idx),
                    state_view.GetConstColumnView(solvent_idx),
                  bv_gg,
                  bv_ga,
                  bv_gs,
                  bv_ag,
                  bv_aa,
                  bv_as);

              // Indirect entries through r_eff
                for (std::size_t k = 0; k < n_r_eff_deps; ++k)
              {
                  auto bv_r_gas = jac_view.GetBlockView(jac_id_view[idx++]);
                  auto bv_r_aq = jac_view.GetBlockView(jac_id_view[idx++]);
                  jac_view.ForEachBlockStrict(
                      [molar_volume, cond](
                          const micm::Real& r_eff,
                          const micm::Real& N,
                          const micm::Real& phi,
                          const micm::Real& hlc,
                          const micm::Real& T,
                          const micm::Real& gas,
                          const micm::Real& aq,
                          const micm::Real& solvent,
                          const micm::Real& dr_dvar,
                          micm::Real& j_gas,
                          micm::Real& j_aq)
                    {
                        micm::Real kc_dummy, dk_dr, dk_dN_unused;
                        cond.ComputeValueAndDerivatives(r_eff, N, T, kc_dummy, dk_dr, dk_dN_unused);
                        const micm::Real dke_dr = dk_dr / (hlc * micm::constants::GAS_CONSTANT * T);
                        const micm::Real fv = solvent * molar_volume;
                        const micm::Real eff = phi * (dk_dr * dr_dvar * gas - dke_dr * dr_dvar * aq / fv);
                      j_gas += eff;
                      j_aq -= eff;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                      params_view.GetConstColumnView(hlc_idx),
                      params_view.GetConstColumnView(temp_idx),
                      state_view.GetConstColumnView(gas_idx),
                      state_view.GetConstColumnView(aq_idx),
                      state_view.GetConstColumnView(solvent_idx),
                    r_eff_partials_view.GetConstColumnView(k),
                    bv_r_gas,
                    bv_r_aq);
              }

              // Indirect entries through N
                for (std::size_t k = 0; k < n_N_deps; ++k)
              {
                  auto bv_N_gas = jac_view.GetBlockView(jac_id_view[idx++]);
                  auto bv_N_aq = jac_view.GetBlockView(jac_id_view[idx++]);
                  jac_view.ForEachBlockStrict(
                      [molar_volume, cond](
                          const micm::Real& r_eff,
                          const micm::Real& N,
                          const micm::Real& phi,
                          const micm::Real& hlc,
                          const micm::Real& T,
                          const micm::Real& gas,
                          const micm::Real& aq,
                          const micm::Real& solvent,
                          const micm::Real& dN_dvar,
                          micm::Real& j_gas,
                          micm::Real& j_aq)
                    {
                        micm::Real kc_dummy, dk_dr_unused, dk_dN;
                        cond.ComputeValueAndDerivatives(r_eff, N, T, kc_dummy, dk_dr_unused, dk_dN);
                        const micm::Real dke_dN = dk_dN / (hlc * micm::constants::GAS_CONSTANT * T);
                        const micm::Real fv = solvent * molar_volume;
                        const micm::Real eff = phi * (dk_dN * dN_dvar * gas - dke_dN * dN_dvar * aq / fv);
                      j_gas += eff;
                      j_aq -= eff;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                      params_view.GetConstColumnView(hlc_idx),
                      params_view.GetConstColumnView(temp_idx),
                      state_view.GetConstColumnView(gas_idx),
                      state_view.GetConstColumnView(aq_idx),
                      state_view.GetConstColumnView(solvent_idx),
                    N_partials_view.GetConstColumnView(k),
                    bv_N_gas,
                    bv_N_aq);
              }

              // Indirect entries through φ_p (negated: MICM solver expects -J)
                for (std::size_t k = 0; k < n_phi_deps; ++k)
              {
                  auto bv_phi_gas = jac_view.GetBlockView(jac_id_view[idx++]);
                  auto bv_phi_aq = jac_view.GetBlockView(jac_id_view[idx++]);
                  jac_view.ForEachBlockStrict(
                      [molar_volume, cond](
                          const micm::Real& r_eff,
                          const micm::Real& N,
                          const micm::Real& phi,
                          const micm::Real& hlc,
                          const micm::Real& T,
                          const micm::Real& gas,
                          const micm::Real& aq,
                          const micm::Real& solvent,
                          const micm::Real& dphi_dvar,
                          micm::Real& j_gas,
                          micm::Real& j_aq)
                    {
                        const micm::Real kc = cond.ComputeValue(r_eff, N, T);
                        const micm::Real ke = kc / (hlc * micm::constants::GAS_CONSTANT * T);
                        const micm::Real fv = solvent * molar_volume;
                        const micm::Real R = kc * gas - ke * aq / fv;
                      j_gas += R * dphi_dvar;
                      j_aq -= R * dphi_dvar;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                      params_view.GetConstColumnView(hlc_idx),
                      params_view.GetConstColumnView(temp_idx),
                      state_view.GetConstColumnView(gas_idx),
                      state_view.GetConstColumnView(aq_idx),
                      state_view.GetConstColumnView(solvent_idx),
                    phi_partials_view.GetConstColumnView(k),
                    bv_phi_gas,
                    bv_phi_aq);
              }
            },
            dummy_state_parameters,
            dummy_state_variables,
            jacobian,
            dummy_buf,
            dummy_buf,
            dummy_buf,
            dummy_partials,
            dummy_partials,
            dummy_partials));
      }

      return [storage, inner_jac_functions = std::move(inner_jac_functions)](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian_matrix)
      {
        std::size_t num_blocks = jacobian_matrix.NumberOfBlocks();

        for (std::size_t i = 0; i < storage->size(); ++i)
        {
          const auto& inst = (*storage)[i];

          // Pre-compute aerosol properties (full-matrix calls)
          DenseMatrixPolicy r_eff_buf{ num_blocks, 1, 0.0 };
          DenseMatrixPolicy N_buf{ num_blocks, 1, 0.0 };
          DenseMatrixPolicy phi_buf{ num_blocks, 1, 0.0 };
          inst.r_eff_provider.ComputeValue(state_parameters, state_variables, r_eff_buf);
          inst.N_provider.ComputeValue(state_parameters, state_variables, N_buf);
          inst.phi_provider.ComputeValue(state_parameters, state_variables, phi_buf);

          // Pre-compute partials
          DenseMatrixPolicy r_eff_partials{ num_blocks, std::max(inst.n_r_eff_deps, std::size_t(1)), 0.0 };
          if (inst.n_r_eff_deps > 0)
            inst.r_eff_provider.ComputeValueAndDerivatives(state_parameters, state_variables, r_eff_buf, r_eff_partials);

          DenseMatrixPolicy N_partials{ num_blocks, std::max(inst.n_N_deps, std::size_t(1)), 0.0 };
          if (inst.n_N_deps > 0)
            inst.N_provider.ComputeValueAndDerivatives(state_parameters, state_variables, N_buf, N_partials);

          DenseMatrixPolicy phi_partials{ num_blocks, std::max(inst.n_phi_deps, std::size_t(1)), 0.0 };
          if (inst.n_phi_deps > 0)
            inst.phi_provider.ComputeValueAndDerivatives(state_parameters, state_variables, phi_buf, phi_partials);

          // Call the Function that was built for this instance
          inner_jac_functions[i](
              state_parameters,
              state_variables,
              jacobian_matrix,
              r_eff_buf,
              N_buf,
              phi_buf,
              r_eff_partials,
              N_partials,
              phi_partials);
        }
      };
    }
  };
}  // namespace miam
