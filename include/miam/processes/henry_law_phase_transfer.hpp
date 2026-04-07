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
        // Return the direct dependencies; the provider-based overload adds indirect ones.
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

        return [instances = std::move(instances), gas_idx](
                   const DenseMatrixPolicy& state_parameters,
                   const DenseMatrixPolicy& state_variables,
                   DenseMatrixPolicy& forcing_terms)
        {
          std::size_t num_rows = state_parameters.NumRows();

          for (const auto& inst : instances)
          {
            // Pre-compute aerosol properties into temp buffers (full-matrix calls)
            DenseMatrixPolicy r_eff_buf{ num_rows, 1, 0.0 };
            DenseMatrixPolicy N_buf{ num_rows, 1, 0.0 };
            DenseMatrixPolicy phi_buf{ num_rows, 1, 0.0 };
            inst.r_eff_provider.ComputeValue(state_parameters, state_variables, r_eff_buf);
            inst.N_provider.ComputeValue(state_parameters, state_variables, N_buf);
            inst.phi_provider.ComputeValue(state_parameters, state_variables, phi_buf);

            // Compute net rate and apply forcing for each grid cell
            for (std::size_t row = 0; row < num_rows; ++row)
            {
              double r_eff_val = r_eff_buf[row][0];
              double N_val = N_buf[row][0];
              double phi_val = phi_buf[row][0];
              double hlc = state_parameters[row][inst.hlc_param_idx];
              double T = state_parameters[row][inst.temperature_param_idx];
              double gas_conc = state_variables[row][gas_idx];
              double aq_conc = state_variables[row][inst.aq_species_idx];
              double solvent_conc = state_variables[row][inst.solvent_species_idx];

              double kc = inst.cond_rate_provider.ComputeValue(r_eff_val, N_val, T);
              double kc_eff = phi_val * kc;
              double ke_eff = kc_eff / (hlc * util::R_gas * T);
              double f_v = solvent_conc * inst.Mw_rho;
              double net = kc_eff * gas_conc - ke_eff * aq_conc / f_v;

              forcing_terms[row][gas_idx] -= net;
              forcing_terms[row][inst.aq_species_idx] += net;
            }
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
        };

        // Jacobian element VectorIndex offsets (block 0) — add block offset at runtime
        struct JacobianInstanceData
        {
          InstanceData instance;
          // Direct entry indices: [gas,gas], [gas,aq], [gas,solvent], [aq,gas], [aq,aq], [aq,solvent]
          std::vector<std::size_t> direct_jac_indices;
          // Indirect r_eff entries: pairs of [gas,var_j], [aq,var_j] per dependent variable
          std::vector<std::size_t> r_eff_jac_indices;
          // Indirect N entries: pairs of [gas,var_j], [aq,var_j] per dependent variable
          std::vector<std::size_t> N_jac_indices;
          // Indirect phi entries: pairs of [gas,var_j], [aq,var_j] per dependent variable
          std::vector<std::size_t> phi_jac_indices;
        };

        std::vector<JacobianInstanceData> jac_instances;
        auto my_jac_phase_it = phase_prefixes.find(condensed_phase_.name_);
        if (my_jac_phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : my_jac_phase_it->second)
          {
            auto prov_it = providers.find(prefix);
            if (prov_it == providers.end())
              continue;
            const auto& prov_map = prov_it->second;
            JacobianInstanceData jac_inst;
            jac_inst.instance.aq_species_idx =
                state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
            jac_inst.instance.solvent_species_idx =
                state_variable_indices.at(prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
            jac_inst.instance.hlc_param_idx =
                state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
            jac_inst.instance.temperature_param_idx =
                state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
            jac_inst.instance.Mw_rho = Mw_solvent_ / rho_solvent_;
            jac_inst.instance.r_eff_provider = prov_map.at(AerosolProperty::EffectiveRadius);
            jac_inst.instance.N_provider = prov_map.at(AerosolProperty::NumberConcentration);
            jac_inst.instance.phi_provider = prov_map.at(AerosolProperty::PhaseVolumeFraction);
            jac_inst.instance.cond_rate_provider = util::MakeCondensationRateProvider(D_g_, alpha_, Mw_gas_);

            std::size_t aq_idx = jac_inst.instance.aq_species_idx;
            std::size_t solvent_idx = jac_inst.instance.solvent_species_idx;

            // Direct entries (6 total)
            jac_inst.direct_jac_indices.push_back(jacobian.VectorIndex(0, gas_idx, gas_idx));
            jac_inst.direct_jac_indices.push_back(jacobian.VectorIndex(0, gas_idx, aq_idx));
            jac_inst.direct_jac_indices.push_back(jacobian.VectorIndex(0, gas_idx, solvent_idx));
            jac_inst.direct_jac_indices.push_back(jacobian.VectorIndex(0, aq_idx, gas_idx));
            jac_inst.direct_jac_indices.push_back(jacobian.VectorIndex(0, aq_idx, aq_idx));
            jac_inst.direct_jac_indices.push_back(jacobian.VectorIndex(0, aq_idx, solvent_idx));

            // Indirect through r_eff
            for (std::size_t var_j : prov_map.at(AerosolProperty::EffectiveRadius).dependent_variable_indices)
            {
              jac_inst.r_eff_jac_indices.push_back(jacobian.VectorIndex(0, gas_idx, var_j));
              jac_inst.r_eff_jac_indices.push_back(jacobian.VectorIndex(0, aq_idx, var_j));
            }
            // Indirect through N
            for (std::size_t var_j : prov_map.at(AerosolProperty::NumberConcentration).dependent_variable_indices)
            {
              jac_inst.N_jac_indices.push_back(jacobian.VectorIndex(0, gas_idx, var_j));
              jac_inst.N_jac_indices.push_back(jacobian.VectorIndex(0, aq_idx, var_j));
            }
            // Indirect through phi
            for (std::size_t var_j : prov_map.at(AerosolProperty::PhaseVolumeFraction).dependent_variable_indices)
            {
              jac_inst.phi_jac_indices.push_back(jacobian.VectorIndex(0, gas_idx, var_j));
              jac_inst.phi_jac_indices.push_back(jacobian.VectorIndex(0, aq_idx, var_j));
            }

            jac_instances.push_back(std::move(jac_inst));
          }
        }

        // Compute the block stride: difference between VectorIndex(0,...) and VectorIndex(1,...)
        // For standard ordering (L=1), this is the number of non-zero elements per block
        std::size_t block_stride = 0;
        if (jacobian.NumberOfBlocks() > 1)
          block_stride = jacobian.VectorIndex(1, 0, 0) - jacobian.VectorIndex(0, 0, 0);
        else if (jacobian.NumberOfBlocks() == 1)
          block_stride = jacobian.AsVector().size();  // Only one block

        return [jac_instances = std::move(jac_instances), gas_idx, block_stride](
                   const DenseMatrixPolicy& state_parameters,
                   const DenseMatrixPolicy& state_variables,
                   SparseMatrixPolicy& jacobian_matrix)
        {
          std::size_t num_blocks = jacobian_matrix.NumberOfBlocks();
          auto& jac_vec = jacobian_matrix.AsVector();

          for (const auto& jac_inst : jac_instances)
          {
            const auto& inst = jac_inst.instance;

            // Pre-compute aerosol properties (full-matrix calls)
            DenseMatrixPolicy r_eff_buf{ num_blocks, 1, 0.0 };
            DenseMatrixPolicy N_buf{ num_blocks, 1, 0.0 };
            DenseMatrixPolicy phi_buf{ num_blocks, 1, 0.0 };
            inst.r_eff_provider.ComputeValue(state_parameters, state_variables, r_eff_buf);
            inst.N_provider.ComputeValue(state_parameters, state_variables, N_buf);
            inst.phi_provider.ComputeValue(state_parameters, state_variables, phi_buf);

            // Pre-compute r_eff partials if needed
            std::size_t n_r_eff_deps = inst.r_eff_provider.dependent_variable_indices.size();
            DenseMatrixPolicy r_eff_partials{ num_blocks, std::max(n_r_eff_deps, std::size_t(1)), 0.0 };
            if (n_r_eff_deps > 0)
              inst.r_eff_provider.ComputeValueAndDerivatives(state_parameters, state_variables, r_eff_buf, r_eff_partials);

            // Pre-compute N partials if needed
            std::size_t n_N_deps = inst.N_provider.dependent_variable_indices.size();
            DenseMatrixPolicy N_partials{ num_blocks, std::max(n_N_deps, std::size_t(1)), 0.0 };
            if (n_N_deps > 0)
              inst.N_provider.ComputeValueAndDerivatives(state_parameters, state_variables, N_buf, N_partials);

            // Pre-compute phi partials if needed
            std::size_t n_phi_deps = inst.phi_provider.dependent_variable_indices.size();
            DenseMatrixPolicy phi_partials{ num_blocks, std::max(n_phi_deps, std::size_t(1)), 0.0 };
            if (n_phi_deps > 0)
              inst.phi_provider.ComputeValueAndDerivatives(state_parameters, state_variables, phi_buf, phi_partials);

            // Iterate over blocks (grid cells)
            for (std::size_t block = 0; block < num_blocks; ++block)
            {
              std::size_t jac_offset = block * block_stride;

              double r_eff_val = r_eff_buf[block][0];
              double N_val = N_buf[block][0];
              double phi_val = phi_buf[block][0];
              double hlc = state_parameters[block][inst.hlc_param_idx];
              double T = state_parameters[block][inst.temperature_param_idx];
              double gas_conc = state_variables[block][gas_idx];
              double aq_conc = state_variables[block][inst.aq_species_idx];
              double solvent_conc = state_variables[block][inst.solvent_species_idx];

              double kc = inst.cond_rate_provider.ComputeValue(r_eff_val, N_val, T);
              double ke = kc / (hlc * util::R_gas * T);
              double fv = solvent_conc * inst.Mw_rho;
              double R = kc * gas_conc - ke * aq_conc / fv;

              // Direct Jacobian entries (negated: MICM solver expects -J)
              // -J[gas, gas] = +φ · k_cond
              jac_vec[jac_inst.direct_jac_indices[0] + jac_offset] += phi_val * kc;
              // -J[gas, aq] = -φ · k_evap / f_v
              jac_vec[jac_inst.direct_jac_indices[1] + jac_offset] -= phi_val * ke / fv;
              // -J[gas, solvent] = +φ · k_evap · [aq] / (f_v · [solvent])
              jac_vec[jac_inst.direct_jac_indices[2] + jac_offset] += phi_val * ke * aq_conc / (fv * solvent_conc);
              // -J[aq, gas] = -φ · k_cond
              jac_vec[jac_inst.direct_jac_indices[3] + jac_offset] -= phi_val * kc;
              // -J[aq, aq] = +φ · k_evap / f_v
              jac_vec[jac_inst.direct_jac_indices[4] + jac_offset] += phi_val * ke / fv;
              // -J[aq, solvent] = -φ · k_evap · [aq] / (f_v · [solvent])
              jac_vec[jac_inst.direct_jac_indices[5] + jac_offset] -= phi_val * ke * aq_conc / (fv * solvent_conc);

              // Indirect entries through r_eff
              if (n_r_eff_deps > 0)
              {
                double dk_dr, dk_dN_unused;
                double kc_unused;
                inst.cond_rate_provider.ComputeValueAndDerivatives(r_eff_val, N_val, T, kc_unused, dk_dr, dk_dN_unused);
                double dke_dr = dk_dr / (hlc * util::R_gas * T);

                for (std::size_t k = 0; k < n_r_eff_deps; ++k)
                {
                  double dr_dvar = r_eff_partials[block][k];
                  double effect = phi_val * (dk_dr * dr_dvar * gas_conc - dke_dr * dr_dvar * aq_conc / fv);
                  jac_vec[jac_inst.r_eff_jac_indices[2 * k] + jac_offset] += effect;
                  jac_vec[jac_inst.r_eff_jac_indices[2 * k + 1] + jac_offset] -= effect;
                }
              }

              // Indirect entries through N
              if (n_N_deps > 0)
              {
                double dk_dr_unused, dk_dN;
                double kc_unused;
                inst.cond_rate_provider.ComputeValueAndDerivatives(r_eff_val, N_val, T, kc_unused, dk_dr_unused, dk_dN);
                double dke_dN = dk_dN / (hlc * util::R_gas * T);

                for (std::size_t k = 0; k < n_N_deps; ++k)
                {
                  double dN_dvar = N_partials[block][k];
                  double effect = phi_val * (dk_dN * dN_dvar * gas_conc - dke_dN * dN_dvar * aq_conc / fv);
                  jac_vec[jac_inst.N_jac_indices[2 * k] + jac_offset] += effect;
                  jac_vec[jac_inst.N_jac_indices[2 * k + 1] + jac_offset] -= effect;
                }
              }

              // Indirect entries through φ_p (negated: MICM solver expects -J)
              for (std::size_t k = 0; k < n_phi_deps; ++k)
              {
                double dphi_dvar = phi_partials[block][k];
                jac_vec[jac_inst.phi_jac_indices[2 * k] + jac_offset] += R * dphi_dvar;
                jac_vec[jac_inst.phi_jac_indices[2 * k + 1] + jac_offset] -= R * dphi_dvar;
              }
            }
          }
        };
      }
    };
  }  // namespace process
}  // namespace miam
