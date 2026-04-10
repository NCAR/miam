// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/aerosol_property.hpp>
#include <miam/util/condensation_rate.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>

#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  namespace process
  {
    /// @brief Henry's Law phase transfer process
    /// @details Represents the transfer of a gas-phase species into a condensed phase
    ///          (and its re-evaporation) governed by Henry's Law equilibrium. The net rate is:
    ///
    ///          d[A]_gas/dt  = -φ_p · k_cond · [A]_gas + φ_p · k_evap · [A]_aq / f_v
    ///          d[A]_aq/dt   = +φ_p · k_cond · [A]_gas - φ_p · k_evap · [A]_aq / f_v
    ///
    ///          where k_evap = k_cond / (HLC · R · T), f_v = [solvent] · Mw/ρ,
    ///          and φ_p is the phase volume fraction.
    class HenryLawPhaseTransfer
    {
     public:
      std::function<double(const micm::Conditions& conditions)> henry_law_constant_;  ///< HLC(T) function [mol m⁻³ Pa⁻¹]
      micm::Species gas_species_;                                                     ///< Gas-phase species
      micm::Species condensed_species_;                                               ///< Condensed-phase solute species
      micm::Species solvent_;                                                         ///< Condensed-phase solvent species
      micm::Phase condensed_phase_;                                                   ///< The condensed phase
      double D_g_;          ///< Gas-phase diffusion coefficient [m² s⁻¹]
      double alpha_;        ///< Mass accommodation coefficient [dimensionless]
      double Mw_gas_;       ///< Gas-phase molecular weight [kg mol⁻¹]
      double Mw_solvent_;   ///< Solvent molecular weight [kg mol⁻¹]
      double rho_solvent_;  ///< Solvent density [kg m⁻³]
      std::string uuid_;    ///< Unique identifier

      HenryLawPhaseTransfer() = delete;

      /// @brief Constructor
      HenryLawPhaseTransfer(
          std::function<double(const micm::Conditions& conditions)> henry_law_constant,
          const micm::Species& gas_species,
          const micm::Species& condensed_species,
          const micm::Species& solvent,
          const micm::Phase& condensed_phase,
          double D_g,
          double alpha,
          double Mw_gas,
          double Mw_solvent,
          double rho_solvent)
          : henry_law_constant_(henry_law_constant),
            gas_species_(gas_species),
            condensed_species_(condensed_species),
            solvent_(solvent),
            condensed_phase_(condensed_phase),
            D_g_(D_g),
            alpha_(alpha),
            Mw_gas_(Mw_gas),
            Mw_solvent_(Mw_solvent),
            rho_solvent_(rho_solvent),
            uuid_(miam::util::generate_uuid_v4())
      {
      }

      /// @brief Create a copy with a new UUID
      HenryLawPhaseTransfer CopyWithNewUuid() const
      {
        return HenryLawPhaseTransfer(
            henry_law_constant_,
            gas_species_,
            condensed_species_,
            solvent_,
            condensed_phase_,
            D_g_,
            alpha_,
            Mw_gas_,
            Mw_solvent_,
            rho_solvent_);
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
          throw std::runtime_error(
              "Internal Error: Phase " + condensed_phase_.name_ + " not found in phase_prefixes for process " + uuid_);
        }
        return species_names;
      }

      /// @brief Returns the aerosol properties required by this process
      std::map<std::string, std::vector<AerosolProperty>> RequiredAerosolProperties() const
      {
        return { { condensed_phase_.name_,
                   { AerosolProperty::EffectiveRadius,
                     AerosolProperty::NumberConcentration,
                     AerosolProperty::PhaseVolumeFraction } } };
      }

      /// @brief Returns non-zero Jacobian element positions
      std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
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
          throw std::runtime_error("Internal Error: Phase " + condensed_phase_.name_ + " not found in phase_prefixes");

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
          const std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>>& providers)
          const
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
      std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
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
              throw std::runtime_error("Internal Error: HLC parameter " + hlc_param + " not found");
            if (state_parameter_indices.find(temp_param) == state_parameter_indices.end())
              throw std::runtime_error("Internal Error: Temperature parameter " + temp_param + " not found");
            hlc_indices.push_back(state_parameter_indices.at(hlc_param));
            temp_indices.push_back(state_parameter_indices.at(temp_param));
          }
        }

        DenseMatrixPolicy dummy{ 1, state_parameter_indices.size(), 0.0 };
        std::vector<micm::Conditions> dummy_conditions;

        return DenseMatrixPolicy::Function(
            [this, hlc_indices, temp_indices](auto&& conditions, auto&& params)
            {
              for (std::size_t i = 0; i < hlc_indices.size(); ++i)
              {
                params.ForEachRow(
                    [&](const micm::Conditions& cond, double& hlc, double& T)
                    {
                      hlc = henry_law_constant_(cond);
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

      /// @brief Returns a function that calculates the forcing terms (common interface with providers)
      template<typename DenseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const auto& state_parameter_indices,
          const auto& state_variable_indices,
          std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> providers) const
      {
        auto gas_idx = state_variable_indices.at(gas_species_.name_);

        struct InstanceData
        {
          std::size_t aq_species_idx;
          std::size_t solvent_species_idx;
          std::size_t hlc_param_idx;
          std::size_t temperature_param_idx;
          double Mw_rho;
          AerosolPropertyProvider<DenseMatrixPolicy> r_eff_provider;
          AerosolPropertyProvider<DenseMatrixPolicy> N_provider;
          AerosolPropertyProvider<DenseMatrixPolicy> phi_provider;
          util::CondensationRateProvider cond_rate_provider;
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
            inst.solvent_species_idx =
                state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
            inst.hlc_param_idx = state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
            inst.temperature_param_idx =
                state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
            inst.Mw_rho = Mw_solvent_ / rho_solvent_;
            inst.r_eff_provider = prov_map.at(AerosolProperty::EffectiveRadius);
            inst.N_provider = prov_map.at(AerosolProperty::NumberConcentration);
            inst.phi_provider = prov_map.at(AerosolProperty::PhaseVolumeFraction);
            inst.cond_rate_provider = util::MakeCondensationRateProvider(D_g_, alpha_, Mw_gas_);
            instances.push_back(std::move(inst));
          }
        }

        // Create Function-wrapped inner loops — one per instance, built at setup time.
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
          auto inner = DenseMatrixPolicy::Function(
              [inst, gas_idx](auto&& state_parameters, auto&& state_variables, auto&& forcing_terms,
                              auto&& r_eff_view, auto&& N_view, auto&& phi_view)
              {
                auto net = forcing_terms.GetRowVariable();

                // Compute net transfer rate
                state_parameters.ForEachRow(
                    [&inst](const double& r_eff, const double& N, const double& phi,
                            const double& hlc, const double& T,
                            const double& gas, const double& aq, const double& solvent,
                            double& net_val)
                    {
                      double kc = inst.cond_rate_provider.ComputeValue(r_eff, N, T);
                      double kc_eff = phi * kc;
                      double ke_eff = kc_eff / (hlc * util::R_gas * T);
                      double fv = solvent * inst.Mw_rho;
                      net_val = kc_eff * gas - ke_eff * aq / fv;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                    state_parameters.GetConstColumnView(inst.hlc_param_idx),
                    state_parameters.GetConstColumnView(inst.temperature_param_idx),
                    state_variables.GetConstColumnView(gas_idx),
                    state_variables.GetConstColumnView(inst.aq_species_idx),
                    state_variables.GetConstColumnView(inst.solvent_species_idx),
                    net);

                // Apply to gas forcing (subtract)
                state_parameters.ForEachRow(
                    [](const double& net_val, double& f_gas) { f_gas -= net_val; },
                    net,
                    forcing_terms.GetColumnView(gas_idx));

                // Apply to aq forcing (add)
                state_parameters.ForEachRow(
                    [](const double& net_val, double& f_aq) { f_aq += net_val; },
                    net,
                    forcing_terms.GetColumnView(inst.aq_species_idx));
              },
              dummy_state_parameters,
              dummy_state_variables,
              dummy_state_variables,
              dummy_buf,
              dummy_buf,
              dummy_buf);
          inner_functions.push_back(std::move(inner));
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

            // Call the Function-wrapped inner loop with actual data
            inner_functions[i](state_parameters, state_variables, forcing_terms, r_eff_buf, N_buf, phi_buf);
          }
        };
      }

      /// @brief Returns a function that calculates Jacobian contributions (common interface with providers)
      template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const auto& state_parameter_indices,
          const auto& state_variable_indices,
          const SparseMatrixPolicy& jacobian,
          std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> providers) const
      {
        auto gas_idx = state_variable_indices.at(gas_species_.name_);

        struct InstanceData
        {
          std::size_t aq_species_idx;
          std::size_t solvent_species_idx;
          std::size_t hlc_param_idx;
          std::size_t temperature_param_idx;
          double Mw_rho;
          AerosolPropertyProvider<DenseMatrixPolicy> r_eff_provider;
          AerosolPropertyProvider<DenseMatrixPolicy> N_provider;
          AerosolPropertyProvider<DenseMatrixPolicy> phi_provider;
          util::CondensationRateProvider cond_rate_provider;
          std::size_t n_r_eff_deps;
          std::size_t n_N_deps;
          std::size_t n_phi_deps;
          // Jacobian indices stored in a flat Matrix<std::size_t> (1 x N) for use with GetBlockView via *jac_id++
          // Layout: [6 direct] [2*n_r_eff_deps indirect_r_eff] [2*n_N_deps indirect_N] [2*n_phi_deps indirect_phi]
          micm::Matrix<std::size_t> jac_indices;
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
            inst.solvent_species_idx =
                state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
            inst.hlc_param_idx =
                state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
            inst.temperature_param_idx =
                state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
            inst.Mw_rho = Mw_solvent_ / rho_solvent_;
            inst.r_eff_provider = prov_map.at(AerosolProperty::EffectiveRadius);
            inst.N_provider = prov_map.at(AerosolProperty::NumberConcentration);
            inst.phi_provider = prov_map.at(AerosolProperty::PhaseVolumeFraction);
            inst.cond_rate_provider = util::MakeCondensationRateProvider(D_g_, alpha_, Mw_gas_);
            inst.n_r_eff_deps = prov_map.at(AerosolProperty::EffectiveRadius).dependent_variable_indices.size();
            inst.n_N_deps = prov_map.at(AerosolProperty::NumberConcentration).dependent_variable_indices.size();
            inst.n_phi_deps = prov_map.at(AerosolProperty::PhaseVolumeFraction).dependent_variable_indices.size();

            std::size_t aq_idx = inst.aq_species_idx;
            std::size_t solvent_idx = inst.solvent_species_idx;

            // Build flat Jacobian index vector
            std::size_t total_indices = 6 + 2 * inst.n_r_eff_deps + 2 * inst.n_N_deps + 2 * inst.n_phi_deps;
            inst.jac_indices = micm::Matrix<std::size_t>(1, total_indices);
            std::size_t idx = 0;

            // Direct entries (6 total)
            inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, gas_idx, gas_idx);
            inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, gas_idx, aq_idx);
            inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, gas_idx, solvent_idx);
            inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, aq_idx, gas_idx);
            inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, aq_idx, aq_idx);
            inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, aq_idx, solvent_idx);

            // Indirect through r_eff
            for (std::size_t var_j : prov_map.at(AerosolProperty::EffectiveRadius).dependent_variable_indices)
            {
              inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, gas_idx, var_j);
              inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, aq_idx, var_j);
            }
            // Indirect through N
            for (std::size_t var_j : prov_map.at(AerosolProperty::NumberConcentration).dependent_variable_indices)
            {
              inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, gas_idx, var_j);
              inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, aq_idx, var_j);
            }
            // Indirect through phi
            for (std::size_t var_j : prov_map.at(AerosolProperty::PhaseVolumeFraction).dependent_variable_indices)
            {
              inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, gas_idx, var_j);
              inst.jac_indices[0][idx++] = jacobian.VectorIndex(0, aq_idx, var_j);
            }

            jac_instances.push_back(std::move(inst));
          }
        }

        // Build Function-wrapped inner loops per instance
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

        for (const auto& inst : jac_instances)
        {
          inner_jac_functions.push_back(SparseMatrixPolicy::Function(
              [inst, gas_idx](auto&& state_parameters, auto&& state_variables, auto&& jacobian_values,
                              auto&& r_eff_view, auto&& N_view, auto&& phi_view,
                              auto&& r_eff_partials_view, auto&& N_partials_view, auto&& phi_partials_view)
              {
                auto jac_id = inst.jac_indices.AsVector().begin();

                // Pre-extract BlockViews sequentially to avoid unspecified argument evaluation order
                auto bv_gg = jacobian_values.GetBlockView(*jac_id++);
                auto bv_ga = jacobian_values.GetBlockView(*jac_id++);
                auto bv_gs = jacobian_values.GetBlockView(*jac_id++);
                auto bv_ag = jacobian_values.GetBlockView(*jac_id++);
                auto bv_aa = jacobian_values.GetBlockView(*jac_id++);
                auto bv_as = jacobian_values.GetBlockView(*jac_id++);

                // Read inputs and compute direct Jacobian entries
                jacobian_values.ForEachBlock(
                    [&inst](const double& r_eff, const double& N, const double& phi,
                            const double& hlc, const double& T,
                            const double& gas, const double& aq, const double& solvent,
                            double& j_gg, double& j_ga, double& j_gs,
                            double& j_ag, double& j_aa, double& j_as)
                    {
                      double kc = inst.cond_rate_provider.ComputeValue(r_eff, N, T);
                      double ke = kc / (hlc * util::R_gas * T);
                      double fv = solvent * inst.Mw_rho;
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
                    state_parameters.GetConstColumnView(inst.hlc_param_idx),
                    state_parameters.GetConstColumnView(inst.temperature_param_idx),
                    state_variables.GetConstColumnView(gas_idx),
                    state_variables.GetConstColumnView(inst.aq_species_idx),
                    state_variables.GetConstColumnView(inst.solvent_species_idx),
                    bv_gg, bv_ga, bv_gs, bv_ag, bv_aa, bv_as);

                // Indirect entries through r_eff
                for (std::size_t k = 0; k < inst.n_r_eff_deps; ++k)
                {
                  auto bv_r_gas = jacobian_values.GetBlockView(*jac_id++);
                  auto bv_r_aq = jacobian_values.GetBlockView(*jac_id++);
                  jacobian_values.ForEachBlock(
                      [&inst](const double& r_eff, const double& N, const double& phi,
                              const double& hlc, const double& T,
                              const double& gas, const double& aq, const double& solvent,
                              const double& dr_dvar,
                              double& j_gas, double& j_aq)
                      {
                        double kc_dummy, dk_dr, dk_dN_unused;
                        inst.cond_rate_provider.ComputeValueAndDerivatives(r_eff, N, T, kc_dummy, dk_dr, dk_dN_unused);
                        double dke_dr = dk_dr / (hlc * util::R_gas * T);
                        double fv = solvent * inst.Mw_rho;
                        double eff = phi * (dk_dr * dr_dvar * gas - dke_dr * dr_dvar * aq / fv);
                        j_gas += eff;
                        j_aq -= eff;
                      },
                      r_eff_view.GetConstColumnView(0),
                      N_view.GetConstColumnView(0),
                      phi_view.GetConstColumnView(0),
                      state_parameters.GetConstColumnView(inst.hlc_param_idx),
                      state_parameters.GetConstColumnView(inst.temperature_param_idx),
                      state_variables.GetConstColumnView(gas_idx),
                      state_variables.GetConstColumnView(inst.aq_species_idx),
                      state_variables.GetConstColumnView(inst.solvent_species_idx),
                      r_eff_partials_view.GetConstColumnView(k),
                      bv_r_gas, bv_r_aq);
                }

                // Indirect entries through N
                for (std::size_t k = 0; k < inst.n_N_deps; ++k)
                {
                  auto bv_N_gas = jacobian_values.GetBlockView(*jac_id++);
                  auto bv_N_aq = jacobian_values.GetBlockView(*jac_id++);
                  jacobian_values.ForEachBlock(
                      [&inst](const double& r_eff, const double& N, const double& phi,
                              const double& hlc, const double& T,
                              const double& gas, const double& aq, const double& solvent,
                              const double& dN_dvar,
                              double& j_gas, double& j_aq)
                      {
                        double kc_dummy, dk_dr_unused, dk_dN;
                        inst.cond_rate_provider.ComputeValueAndDerivatives(r_eff, N, T, kc_dummy, dk_dr_unused, dk_dN);
                        double dke_dN = dk_dN / (hlc * util::R_gas * T);
                        double fv = solvent * inst.Mw_rho;
                        double eff = phi * (dk_dN * dN_dvar * gas - dke_dN * dN_dvar * aq / fv);
                        j_gas += eff;
                        j_aq -= eff;
                      },
                      r_eff_view.GetConstColumnView(0),
                      N_view.GetConstColumnView(0),
                      phi_view.GetConstColumnView(0),
                      state_parameters.GetConstColumnView(inst.hlc_param_idx),
                      state_parameters.GetConstColumnView(inst.temperature_param_idx),
                      state_variables.GetConstColumnView(gas_idx),
                      state_variables.GetConstColumnView(inst.aq_species_idx),
                      state_variables.GetConstColumnView(inst.solvent_species_idx),
                      N_partials_view.GetConstColumnView(k),
                      bv_N_gas, bv_N_aq);
                }

                // Indirect entries through φ_p (negated: MICM solver expects -J)
                for (std::size_t k = 0; k < inst.n_phi_deps; ++k)
                {
                  auto bv_phi_gas = jacobian_values.GetBlockView(*jac_id++);
                  auto bv_phi_aq = jacobian_values.GetBlockView(*jac_id++);
                  jacobian_values.ForEachBlock(
                      [&inst](const double& r_eff, const double& N, const double& phi,
                              const double& hlc, const double& T,
                              const double& gas, const double& aq, const double& solvent,
                              const double& dphi_dvar,
                              double& j_gas, double& j_aq)
                      {
                        double kc = inst.cond_rate_provider.ComputeValue(r_eff, N, T);
                        double ke = kc / (hlc * util::R_gas * T);
                        double fv = solvent * inst.Mw_rho;
                        double R = kc * gas - ke * aq / fv;
                        j_gas += R * dphi_dvar;
                        j_aq -= R * dphi_dvar;
                      },
                      r_eff_view.GetConstColumnView(0),
                      N_view.GetConstColumnView(0),
                      phi_view.GetConstColumnView(0),
                      state_parameters.GetConstColumnView(inst.hlc_param_idx),
                      state_parameters.GetConstColumnView(inst.temperature_param_idx),
                      state_variables.GetConstColumnView(gas_idx),
                      state_variables.GetConstColumnView(inst.aq_species_idx),
                      state_variables.GetConstColumnView(inst.solvent_species_idx),
                      phi_partials_view.GetConstColumnView(k),
                      bv_phi_gas, bv_phi_aq);
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

        return [jac_instances = std::move(jac_instances), inner_jac_functions = std::move(inner_jac_functions), gas_idx](
                   const DenseMatrixPolicy& state_parameters,
                   const DenseMatrixPolicy& state_variables,
                   SparseMatrixPolicy& jacobian_matrix)
        {
          std::size_t num_blocks = jacobian_matrix.NumberOfBlocks();

          for (std::size_t i = 0; i < jac_instances.size(); ++i)
          {
            const auto& inst = jac_instances[i];

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

            // Call the Function-wrapped inner loop
            inner_jac_functions[i](
                state_parameters, state_variables, jacobian_matrix,
                r_eff_buf, N_buf, phi_buf,
                r_eff_partials, N_partials, phi_partials);
          }
        };
      }
    };
  }  // namespace process
}  // namespace miam
