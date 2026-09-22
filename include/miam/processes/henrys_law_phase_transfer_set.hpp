// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/math/condensation_rate.hpp>
#include <miam/processes/henrys_law_phase_transfer.hpp>
#include <miam/representations/aerosol_property.hpp>
#include <miam/representations/aerosol_property_descriptor.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/util/constants.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `HenrysLawPhaseTransfer`, mirroring `micm::ProcessSet`.
  /// @details Non-templated data-holding class populated once at Finalize time with a
  ///          reference (host-side) descriptor map so the Jacobian flat-ID layout, which
  ///          depends on `n_deps` for `r_eff` / `N` / `phi`, can be pre-computed. At solve
  ///          time the caller passes a policy-typed descriptor map back in; the Set's
  ///          methods evaluate aerosol properties into per-instance buffers on the host,
  ///          then dispatch a per-instance kernel that reads from those buffers.
  class HenrysLawPhaseTransferSet
  {
   public:
    HenrysLawPhaseTransferSet() = default;

    template<typename SparseMatrixPolicy, typename DescriptorMap>
    HenrysLawPhaseTransferSet(
        const HenrysLawPhaseTransfer& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,
        const auto& state_variable_indices,
        const SparseMatrixPolicy& jacobian,
        const DescriptorMap& reference_descriptors)
    {
      gas_species_index_ = LookupSpecies(state_variable_indices, "", "", config.gas_species_.name_);
      molar_volume_ = config.solvent_molecular_weight_ / config.solvent_density_;
      cond_rate_provider_ = MakeCondensationRateProvider(
          config.diffusion_coefficient_, config.accommodation_coefficient_, config.gas_molecular_weight_);

      auto phase_it = phase_prefixes.find(config.condensed_phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "HenrysLawPhaseTransferSet: phase " + config.condensed_phase_.name_ +
                " not found in phase_prefixes for process " + config.uuid_);

      instances_.clear();
      for (const auto& prefix : phase_it->second)
      {
        auto desc_it = reference_descriptors.find(prefix);
        if (desc_it == reference_descriptors.end())
          continue;
        const auto& desc_map = desc_it->second;

        InstanceData inst;
        inst.prefix = prefix;
        inst.aq_species_idx = state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.condensed_species_.name_);
        inst.solvent_species_idx =
            state_variable_indices.at(prefix + "." + config.condensed_phase_.name_ + "." + config.solvent_.name_);
        inst.hlc_param_idx =
            state_parameter_indices.at(prefix + "." + config.condensed_phase_.name_ + "." + config.uuid_ + ".hlc");
        inst.temperature_param_idx = state_parameter_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.uuid_ + ".temperature");

        const auto& r_eff_deps = DependentVariableIndices(desc_map.at(AerosolProperty::EffectiveRadius));
        const auto& N_deps = DependentVariableIndices(desc_map.at(AerosolProperty::NumberConcentration));
        const auto& phi_deps = DependentVariableIndices(desc_map.at(AerosolProperty::PhaseVolumeFraction));
        inst.n_r_eff_deps = r_eff_deps.size();
        inst.n_N_deps = N_deps.size();
        inst.n_phi_deps = phi_deps.size();

        const std::size_t gas_idx = gas_species_index_;
        const std::size_t aq_idx = inst.aq_species_idx;
        const std::size_t solvent_idx = inst.solvent_species_idx;

        inst.jac_flat_ids.reserve(6 + 2 * inst.n_r_eff_deps + 2 * inst.n_N_deps + 2 * inst.n_phi_deps);
        inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, gas_idx, gas_idx));
        inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, gas_idx, aq_idx));
        inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, gas_idx, solvent_idx));
        inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, aq_idx, gas_idx));
        inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, aq_idx, aq_idx));
        inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, aq_idx, solvent_idx));
        for (std::size_t var_j : r_eff_deps)
        {
          inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, gas_idx, var_j));
          inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, aq_idx, var_j));
        }
        for (std::size_t var_j : N_deps)
        {
          inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, gas_idx, var_j));
          inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, aq_idx, var_j));
        }
        for (std::size_t var_j : phi_deps)
        {
          inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, gas_idx, var_j));
          inst.jac_flat_ids.push_back(jacobian.VectorIndex(0, aq_idx, var_j));
        }
        instances_.push_back(std::move(inst));
      }
    }

    template<typename DenseMatrixPolicy, typename DescriptorMap>
    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing,
        const DescriptorMap& descriptors) const
    {
      const std::size_t num_rows = state_parameters.NumRows();
      const std::size_t gas_idx = gas_species_index_;
      const double molar_volume = molar_volume_;
      const CondensationRateProvider cond = cond_rate_provider_;

      for (const auto& inst : instances_)
      {
        const auto desc_it = descriptors.find(inst.prefix);
        if (desc_it == descriptors.end())
          continue;
        const auto& desc_map = desc_it->second;

        DenseMatrixPolicy r_eff_buf{ num_rows, 1, 0.0 };
        DenseMatrixPolicy N_buf{ num_rows, 1, 0.0 };
        DenseMatrixPolicy phi_buf{ num_rows, 1, 0.0 };
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::EffectiveRadius), state_parameters, state_variables, r_eff_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::NumberConcentration), state_parameters, state_variables, N_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::PhaseVolumeFraction), state_parameters, state_variables, phi_buf);

        const std::size_t aq_idx = inst.aq_species_idx;
        const std::size_t solvent_idx = inst.solvent_species_idx;
        const std::size_t hlc_idx = inst.hlc_param_idx;
        const std::size_t temp_idx = inst.temperature_param_idx;

        DenseMatrixPolicy::Function(
            [gas_idx, aq_idx, solvent_idx, hlc_idx, temp_idx, molar_volume, cond](
                auto&& params, auto&& vars, auto&& forcing_view, auto&& r_eff_view, auto&& N_view, auto&& phi_view)
            {
              auto net = forcing_view.GetRowVariable();

              params.ForEachRow(
                  [molar_volume, cond](
                      const double& r_eff,
                      const double& N,
                      const double& phi,
                      const double& hlc,
                      const double& T,
                      const double& gas,
                      const double& aq,
                      const double& solvent,
                      double& net_val)
                  {
                    double kc = cond.ComputeValue(r_eff, N, T);
                    double kc_eff = phi * kc;
                    double ke_eff = kc_eff / (hlc * micm::constants::GAS_CONSTANT * T);
                    double fv = solvent * molar_volume;
                    net_val = kc_eff * gas - ke_eff * aq / fv;
                  },
                  r_eff_view.GetConstColumnView(0),
                  N_view.GetConstColumnView(0),
                  phi_view.GetConstColumnView(0),
                  params.GetConstColumnView(hlc_idx),
                  params.GetConstColumnView(temp_idx),
                  vars.GetConstColumnView(gas_idx),
                  vars.GetConstColumnView(aq_idx),
                  vars.GetConstColumnView(solvent_idx),
                  net);

              params.ForEachRow(
                  [](const double& n, double& f) { f -= n; }, net, forcing_view.GetColumnView(gas_idx));
              params.ForEachRow(
                  [](const double& n, double& f) { f += n; }, net, forcing_view.GetColumnView(aq_idx));
            },
            state_parameters,
            state_variables,
            forcing,
            r_eff_buf,
            N_buf,
            phi_buf)(state_parameters, state_variables, forcing, r_eff_buf, N_buf, phi_buf);
      }
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy, typename DescriptorMap>
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian,
        const DescriptorMap& descriptors) const
    {
      const std::size_t num_blocks = jacobian.NumberOfBlocks();
      const std::size_t gas_idx = gas_species_index_;
      const double molar_volume = molar_volume_;
      const CondensationRateProvider cond = cond_rate_provider_;

      for (const auto& inst : instances_)
      {
        const auto desc_it = descriptors.find(inst.prefix);
        if (desc_it == descriptors.end())
          continue;
        const auto& desc_map = desc_it->second;

        DenseMatrixPolicy r_eff_buf{ num_blocks, 1, 0.0 };
        DenseMatrixPolicy N_buf{ num_blocks, 1, 0.0 };
        DenseMatrixPolicy phi_buf{ num_blocks, 1, 0.0 };
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::EffectiveRadius), state_parameters, state_variables, r_eff_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::NumberConcentration), state_parameters, state_variables, N_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::PhaseVolumeFraction), state_parameters, state_variables, phi_buf);

        DenseMatrixPolicy r_eff_partials{ num_blocks, std::max(inst.n_r_eff_deps, std::size_t(1)), 0.0 };
        if (inst.n_r_eff_deps > 0)
          EvaluateAerosolPropertyAndDerivatives(
              desc_map.at(AerosolProperty::EffectiveRadius), state_parameters, state_variables, r_eff_buf, r_eff_partials);
        DenseMatrixPolicy N_partials{ num_blocks, std::max(inst.n_N_deps, std::size_t(1)), 0.0 };
        if (inst.n_N_deps > 0)
          EvaluateAerosolPropertyAndDerivatives(
              desc_map.at(AerosolProperty::NumberConcentration), state_parameters, state_variables, N_buf, N_partials);
        DenseMatrixPolicy phi_partials{ num_blocks, std::max(inst.n_phi_deps, std::size_t(1)), 0.0 };
        if (inst.n_phi_deps > 0)
          EvaluateAerosolPropertyAndDerivatives(
              desc_map.at(AerosolProperty::PhaseVolumeFraction), state_parameters, state_variables, phi_buf, phi_partials);

        const std::size_t aq_idx = inst.aq_species_idx;
        const std::size_t solvent_idx = inst.solvent_species_idx;
        const std::size_t hlc_idx = inst.hlc_param_idx;
        const std::size_t temp_idx = inst.temperature_param_idx;
        const std::size_t n_r_eff = inst.n_r_eff_deps;
        const std::size_t n_N = inst.n_N_deps;
        const std::size_t n_phi = inst.n_phi_deps;
        const std::size_t* flat_ids = inst.jac_flat_ids.data();

        SparseMatrixPolicy::Function(
            [gas_idx, aq_idx, solvent_idx, hlc_idx, temp_idx, molar_volume, cond, n_r_eff, n_N, n_phi, flat_ids](
                auto&& params,
                auto&& vars,
                auto&& jacobian_values,
                auto&& r_eff_view,
                auto&& N_view,
                auto&& phi_view,
                auto&& r_eff_partials_view,
                auto&& N_partials_view,
                auto&& phi_partials_view)
            {
              std::size_t idx = 0;

              auto bv_gg = jacobian_values.GetBlockView(flat_ids[idx++]);
              auto bv_ga = jacobian_values.GetBlockView(flat_ids[idx++]);
              auto bv_gs = jacobian_values.GetBlockView(flat_ids[idx++]);
              auto bv_ag = jacobian_values.GetBlockView(flat_ids[idx++]);
              auto bv_aa = jacobian_values.GetBlockView(flat_ids[idx++]);
              auto bv_as = jacobian_values.GetBlockView(flat_ids[idx++]);

              jacobian_values.ForEachBlock(
                  [molar_volume, cond](
                      const double& r_eff,
                      const double& N,
                      const double& phi,
                      const double& hlc,
                      const double& T,
                      const double& gas,
                      const double& aq,
                      const double& solvent,
                      double& j_gg,
                      double& j_ga,
                      double& j_gs,
                      double& j_ag,
                      double& j_aa,
                      double& j_as)
                  {
                    double kc = cond.ComputeValue(r_eff, N, T);
                    double ke = kc / (hlc * micm::constants::GAS_CONSTANT * T);
                    double fv = solvent * molar_volume;
                    j_gg += phi * kc;
                    j_ga -= phi * ke / fv;
                    j_gs += phi * ke * aq / (fv * solvent);
                    j_ag -= phi * kc;
                    j_aa += phi * ke / fv;
                    j_as -= phi * ke * aq / (fv * solvent);
                  },
                  r_eff_view.GetConstColumnView(0),
                  N_view.GetConstColumnView(0),
                  phi_view.GetConstColumnView(0),
                  params.GetConstColumnView(hlc_idx),
                  params.GetConstColumnView(temp_idx),
                  vars.GetConstColumnView(gas_idx),
                  vars.GetConstColumnView(aq_idx),
                  vars.GetConstColumnView(solvent_idx),
                  bv_gg,
                  bv_ga,
                  bv_gs,
                  bv_ag,
                  bv_aa,
                  bv_as);

              for (std::size_t k = 0; k < n_r_eff; ++k)
              {
                auto bv_r_gas = jacobian_values.GetBlockView(flat_ids[idx++]);
                auto bv_r_aq = jacobian_values.GetBlockView(flat_ids[idx++]);
                jacobian_values.ForEachBlock(
                    [molar_volume, cond](
                        const double& r_eff,
                        const double& N,
                        const double& phi,
                        const double& hlc,
                        const double& T,
                        const double& gas,
                        const double& aq,
                        const double& solvent,
                        const double& dr_dvar,
                        double& j_gas,
                        double& j_aq)
                    {
                      double kc_dummy, dk_dr, dk_dN_unused;
                      cond.ComputeValueAndDerivatives(r_eff, N, T, kc_dummy, dk_dr, dk_dN_unused);
                      double dke_dr = dk_dr / (hlc * micm::constants::GAS_CONSTANT * T);
                      double fv = solvent * molar_volume;
                      double eff = phi * (dk_dr * dr_dvar * gas - dke_dr * dr_dvar * aq / fv);
                      j_gas += eff;
                      j_aq -= eff;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                    params.GetConstColumnView(hlc_idx),
                    params.GetConstColumnView(temp_idx),
                    vars.GetConstColumnView(gas_idx),
                    vars.GetConstColumnView(aq_idx),
                    vars.GetConstColumnView(solvent_idx),
                    r_eff_partials_view.GetConstColumnView(k),
                    bv_r_gas,
                    bv_r_aq);
              }

              for (std::size_t k = 0; k < n_N; ++k)
              {
                auto bv_N_gas = jacobian_values.GetBlockView(flat_ids[idx++]);
                auto bv_N_aq = jacobian_values.GetBlockView(flat_ids[idx++]);
                jacobian_values.ForEachBlock(
                    [molar_volume, cond](
                        const double& r_eff,
                        const double& N,
                        const double& phi,
                        const double& hlc,
                        const double& T,
                        const double& gas,
                        const double& aq,
                        const double& solvent,
                        const double& dN_dvar,
                        double& j_gas,
                        double& j_aq)
                    {
                      double kc_dummy, dk_dr_unused, dk_dN;
                      cond.ComputeValueAndDerivatives(r_eff, N, T, kc_dummy, dk_dr_unused, dk_dN);
                      double dke_dN = dk_dN / (hlc * micm::constants::GAS_CONSTANT * T);
                      double fv = solvent * molar_volume;
                      double eff = phi * (dk_dN * dN_dvar * gas - dke_dN * dN_dvar * aq / fv);
                      j_gas += eff;
                      j_aq -= eff;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                    params.GetConstColumnView(hlc_idx),
                    params.GetConstColumnView(temp_idx),
                    vars.GetConstColumnView(gas_idx),
                    vars.GetConstColumnView(aq_idx),
                    vars.GetConstColumnView(solvent_idx),
                    N_partials_view.GetConstColumnView(k),
                    bv_N_gas,
                    bv_N_aq);
              }

              for (std::size_t k = 0; k < n_phi; ++k)
              {
                auto bv_phi_gas = jacobian_values.GetBlockView(flat_ids[idx++]);
                auto bv_phi_aq = jacobian_values.GetBlockView(flat_ids[idx++]);
                jacobian_values.ForEachBlock(
                    [molar_volume, cond](
                        const double& r_eff,
                        const double& N,
                        const double& phi,
                        const double& hlc,
                        const double& T,
                        const double& gas,
                        const double& aq,
                        const double& solvent,
                        const double& dphi_dvar,
                        double& j_gas,
                        double& j_aq)
                    {
                      double kc = cond.ComputeValue(r_eff, N, T);
                      double ke = kc / (hlc * micm::constants::GAS_CONSTANT * T);
                      double fv = solvent * molar_volume;
                      double R = kc * gas - ke * aq / fv;
                      j_gas += R * dphi_dvar;
                      j_aq -= R * dphi_dvar;
                    },
                    r_eff_view.GetConstColumnView(0),
                    N_view.GetConstColumnView(0),
                    phi_view.GetConstColumnView(0),
                    params.GetConstColumnView(hlc_idx),
                    params.GetConstColumnView(temp_idx),
                    vars.GetConstColumnView(gas_idx),
                    vars.GetConstColumnView(aq_idx),
                    vars.GetConstColumnView(solvent_idx),
                    phi_partials_view.GetConstColumnView(k),
                    bv_phi_gas,
                    bv_phi_aq);
              }
            },
            state_parameters,
            state_variables,
            jacobian,
            r_eff_buf,
            N_buf,
            phi_buf,
            r_eff_partials,
            N_partials,
            phi_partials)(
            state_parameters,
            state_variables,
            jacobian,
            r_eff_buf,
            N_buf,
            phi_buf,
            r_eff_partials,
            N_partials,
            phi_partials);
      }
    }

   private:
    struct InstanceData
    {
      std::string prefix;
      std::size_t aq_species_idx = 0;
      std::size_t solvent_species_idx = 0;
      std::size_t hlc_param_idx = 0;
      std::size_t temperature_param_idx = 0;
      std::size_t n_r_eff_deps = 0;
      std::size_t n_N_deps = 0;
      std::size_t n_phi_deps = 0;
      /// Layout: [6 direct] [2 * n_r_eff_deps] [2 * n_N_deps] [2 * n_phi_deps].
      std::vector<std::size_t> jac_flat_ids;
    };

    static std::size_t LookupSpecies(
        const auto& state_variable_indices,
        const std::string& /*prefix*/,
        const std::string& /*phase_name*/,
        const std::string& species_name)
    {
      auto it = state_variable_indices.find(species_name);
      if (it == state_variable_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_VARIABLE,
            "HenrysLawPhaseTransferSet: state variable " + species_name + " not found");
      return it->second;
    }

    std::vector<InstanceData> instances_;
    std::size_t gas_species_index_ = 0;
    double molar_volume_ = 0.0;
    CondensationRateProvider cond_rate_provider_{};
  };
}  // namespace miam
