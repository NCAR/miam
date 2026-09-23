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
  /// @brief Solve-time companion to `HenrysLawPhaseTransfer`, mirroring `micm::ProcessSet`.
  /// @details Per-phase-instance scalars and jacobian flat-IDs are stored in
  ///          `SparseMatrixPolicy::VectorType<Index>` containers. Aerosol-property evaluation
  ///          runs on host for each instance; the resulting buffers plus per-instance scalars
  ///          are handed to a `DP::Function` / `SP::Function` kernel via `MICM_LAMBDA`.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class HenrysLawPhaseTransferSet
  {
   public:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;
    template<class U>
    using VectorView = typename Vector<U>::ConstViewType;

    struct Views
    {
      VectorView<micm::Index> aq_species_indices_;
      VectorView<micm::Index> solvent_species_indices_;
      VectorView<micm::Index> hlc_param_indices_;
      VectorView<micm::Index> temperature_param_indices_;
      VectorView<micm::Index> jac_flat_id_offsets_;
      VectorView<micm::Index> jac_flat_ids_;
      micm::Index gas_species_index_;
      micm::Real molar_volume_;
      CondensationRateProvider cond_;

      Views() = default;

      Views(
          const Vector<micm::Index>& aq_species_indices,
          const Vector<micm::Index>& solvent_species_indices,
          const Vector<micm::Index>& hlc_param_indices,
          const Vector<micm::Index>& temperature_param_indices,
          const Vector<micm::Index>& jac_flat_id_offsets,
          const Vector<micm::Index>& jac_flat_ids,
          micm::Index gas_species_index,
          micm::Real molar_volume,
          CondensationRateProvider cond)
          : aq_species_indices_(aq_species_indices.GetView()),
            solvent_species_indices_(solvent_species_indices.GetView()),
            hlc_param_indices_(hlc_param_indices.GetView()),
            temperature_param_indices_(temperature_param_indices.GetView()),
            jac_flat_id_offsets_(jac_flat_id_offsets.GetView()),
            jac_flat_ids_(jac_flat_ids.GetView()),
            gas_species_index_(gas_species_index),
            molar_volume_(molar_volume),
            cond_(cond)
      {
      }
    };

    HenrysLawPhaseTransferSet() = default;

    template<class DescriptorMap>
    HenrysLawPhaseTransferSet(
        const HenrysLawPhaseTransfer& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian,
        const DescriptorMap& reference_descriptors)
    {
      const auto gas_species_index =
          static_cast<micm::Index>(LookupSpecies(state_variable_indices, config.gas_species_.name_));
      const micm::Real molar_volume = config.solvent_molecular_weight_ / config.solvent_density_;
      cond_rate_provider_ = MakeCondensationRateProvider(
          config.diffusion_coefficient_, config.accommodation_coefficient_, config.gas_molecular_weight_);

      auto phase_it = phase_prefixes.find(config.condensed_phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "HenrysLawPhaseTransferSet: phase " + config.condensed_phase_.name_ +
                " not found in phase_prefixes for process " + config.uuid_);

      prefixes_.clear();
      n_r_eff_deps_.clear();
      n_N_deps_.clear();
      n_phi_deps_.clear();

      std::vector<micm::Index> aq_species_indices_host;
      std::vector<micm::Index> solvent_species_indices_host;
      std::vector<micm::Index> hlc_param_indices_host;
      std::vector<micm::Index> temperature_param_indices_host;
      std::vector<micm::Index> jac_flat_id_offsets_host;
      std::vector<micm::Index> jac_flat_ids_host;

      for (const auto& prefix : phase_it->second)
      {
        auto desc_it = reference_descriptors.find(prefix);
        if (desc_it == reference_descriptors.end())
          continue;
        const auto& desc_map = desc_it->second;

        const auto aq_idx = static_cast<micm::Index>(state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.condensed_species_.name_));
        const auto solvent_idx = static_cast<micm::Index>(state_variable_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.solvent_.name_));
        const auto hlc_idx = static_cast<micm::Index>(state_parameter_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.uuid_ + ".hlc"));
        const auto temp_idx = static_cast<micm::Index>(state_parameter_indices.at(
            prefix + "." + config.condensed_phase_.name_ + "." + config.uuid_ + ".temperature"));

        aq_species_indices_host.push_back(aq_idx);
        solvent_species_indices_host.push_back(solvent_idx);
        hlc_param_indices_host.push_back(hlc_idx);
        temperature_param_indices_host.push_back(temp_idx);

        const auto& r_eff_deps = DependentVariableIndices(desc_map.at(AerosolProperty::EffectiveRadius));
        const auto& N_deps = DependentVariableIndices(desc_map.at(AerosolProperty::NumberConcentration));
        const auto& phi_deps = DependentVariableIndices(desc_map.at(AerosolProperty::PhaseVolumeFraction));

        prefixes_.push_back(prefix);
        n_r_eff_deps_.push_back(r_eff_deps.size());
        n_N_deps_.push_back(N_deps.size());
        n_phi_deps_.push_back(phi_deps.size());

        jac_flat_id_offsets_host.push_back(static_cast<micm::Index>(jac_flat_ids_host.size()));
        jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, gas_species_index, gas_species_index)));
        jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, gas_species_index, aq_idx)));
        jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, gas_species_index, solvent_idx)));
        jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, gas_species_index)));
        jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, aq_idx)));
        jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, solvent_idx)));
        for (std::size_t var_j : r_eff_deps)
        {
          jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, gas_species_index, var_j)));
          jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, var_j)));
        }
        for (std::size_t var_j : N_deps)
        {
          jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, gas_species_index, var_j)));
          jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, var_j)));
        }
        for (std::size_t var_j : phi_deps)
        {
          jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, gas_species_index, var_j)));
          jac_flat_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, aq_idx, var_j)));
        }
      }

      aq_species_indices_ = Vector<micm::Index>(std::move(aq_species_indices_host));
      solvent_species_indices_ = Vector<micm::Index>(std::move(solvent_species_indices_host));
      hlc_param_indices_ = Vector<micm::Index>(std::move(hlc_param_indices_host));
      temperature_param_indices_ = Vector<micm::Index>(std::move(temperature_param_indices_host));
      jac_flat_id_offsets_ = Vector<micm::Index>(std::move(jac_flat_id_offsets_host));
      jac_flat_ids_ = Vector<micm::Index>(std::move(jac_flat_ids_host));

      aq_species_indices_.CopyToDevice();
      solvent_species_indices_.CopyToDevice();
      hlc_param_indices_.CopyToDevice();
      temperature_param_indices_.CopyToDevice();
      jac_flat_id_offsets_.CopyToDevice();
      jac_flat_ids_.CopyToDevice();

      views_ = Views(
          aq_species_indices_,
          solvent_species_indices_,
          hlc_param_indices_,
          temperature_param_indices_,
          jac_flat_id_offsets_,
          jac_flat_ids_,
          gas_species_index,
          molar_volume,
          cond_rate_provider_);
    }

    template<class DescriptorMap>
    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing,
        const DescriptorMap& descriptors) const
    {
      const auto& views = views_;
      const std::size_t num_rows = state_parameters.NumRows();
      for (std::size_t i_inst = 0; i_inst < prefixes_.size(); ++i_inst)
      {
        const auto& prefix = prefixes_[i_inst];
        auto desc_it = descriptors.find(prefix);
        if (desc_it == descriptors.end())
          continue;
        const auto& desc_map = desc_it->second;

        DenseMatrixPolicy r_eff_buf{ num_rows, 1, 0.0 };
        DenseMatrixPolicy N_buf{ num_rows, 1, 0.0 };
        DenseMatrixPolicy phi_buf{ num_rows, 1, 0.0 };
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::EffectiveRadius), state_parameters, state_variables, r_eff_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::NumberConcentration), state_parameters, state_variables, N_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::PhaseVolumeFraction), state_parameters, state_variables, phi_buf);

        const micm::Index inst = static_cast<micm::Index>(i_inst);
        DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ViewType& forcing_view,
                const typename DenseMatrixPolicy::ConstViewType& r_eff_view,
                const typename DenseMatrixPolicy::ConstViewType& N_view,
                const typename DenseMatrixPolicy::ConstViewType& phi_view)
            {
              const micm::Index gas_idx = views.gas_species_index_;
              const micm::Index aq_idx = views.aq_species_indices_[inst];
              const micm::Index solvent_idx = views.solvent_species_indices_[inst];
              const micm::Index hlc_idx = views.hlc_param_indices_[inst];
              const micm::Index temp_idx = views.temperature_param_indices_[inst];
              const micm::Real molar_volume = views.molar_volume_;
              const CondensationRateProvider cond = views.cond_;

              auto net = forcing_view.GetRowVariable();
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

              forcing_view.ForEachRowStrict(
                  [](const micm::Real& n, micm::Real& f) { f -= n; }, net, forcing_view.GetColumnView(gas_idx));
              forcing_view.ForEachRowStrict(
                  [](const micm::Real& n, micm::Real& f) { f += n; }, net, forcing_view.GetColumnView(aq_idx));
            },
            state_parameters,
            state_variables,
            forcing,
            r_eff_buf,
            N_buf,
            phi_buf)(state_parameters, state_variables, forcing, r_eff_buf, N_buf, phi_buf);
      }
    }

    template<class DescriptorMap>
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian,
        const DescriptorMap& descriptors) const
    {
      const auto& views = views_;
      const std::size_t num_blocks = jacobian.NumberOfBlocks();
      for (std::size_t i_inst = 0; i_inst < prefixes_.size(); ++i_inst)
      {
        const auto& prefix = prefixes_[i_inst];
        auto desc_it = descriptors.find(prefix);
        if (desc_it == descriptors.end())
          continue;
        const auto& desc_map = desc_it->second;

        const std::size_t n_r_eff = n_r_eff_deps_[i_inst];
        const std::size_t n_N = n_N_deps_[i_inst];
        const std::size_t n_phi = n_phi_deps_[i_inst];

        DenseMatrixPolicy r_eff_buf{ num_blocks, 1, 0.0 };
        DenseMatrixPolicy N_buf{ num_blocks, 1, 0.0 };
        DenseMatrixPolicy phi_buf{ num_blocks, 1, 0.0 };
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::EffectiveRadius), state_parameters, state_variables, r_eff_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::NumberConcentration), state_parameters, state_variables, N_buf);
        EvaluateAerosolProperty(desc_map.at(AerosolProperty::PhaseVolumeFraction), state_parameters, state_variables, phi_buf);

        DenseMatrixPolicy r_eff_partials{ num_blocks, std::max(n_r_eff, std::size_t(1)), 0.0 };
        if (n_r_eff > 0)
          EvaluateAerosolPropertyAndDerivatives(
              desc_map.at(AerosolProperty::EffectiveRadius), state_parameters, state_variables, r_eff_buf, r_eff_partials);
        DenseMatrixPolicy N_partials{ num_blocks, std::max(n_N, std::size_t(1)), 0.0 };
        if (n_N > 0)
          EvaluateAerosolPropertyAndDerivatives(
              desc_map.at(AerosolProperty::NumberConcentration), state_parameters, state_variables, N_buf, N_partials);
        DenseMatrixPolicy phi_partials{ num_blocks, std::max(n_phi, std::size_t(1)), 0.0 };
        if (n_phi > 0)
          EvaluateAerosolPropertyAndDerivatives(
              desc_map.at(AerosolProperty::PhaseVolumeFraction), state_parameters, state_variables, phi_buf, phi_partials);

        const micm::Index inst = static_cast<micm::Index>(i_inst);
        const micm::Index n_r_eff_i = static_cast<micm::Index>(n_r_eff);
        const micm::Index n_N_i = static_cast<micm::Index>(n_N);
        const micm::Index n_phi_i = static_cast<micm::Index>(n_phi);

        SparseMatrixPolicy::Function(
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
              const micm::Index gas_idx = views.gas_species_index_;
              const micm::Index aq_idx = views.aq_species_indices_[inst];
              const micm::Index solvent_idx = views.solvent_species_indices_[inst];
              const micm::Index hlc_idx = views.hlc_param_indices_[inst];
              const micm::Index temp_idx = views.temperature_param_indices_[inst];
              const micm::Real molar_volume = views.molar_volume_;
              const CondensationRateProvider cond = views.cond_;
              micm::Index idx = views.jac_flat_id_offsets_[inst];

              auto bv_gg = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
              auto bv_ga = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
              auto bv_gs = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
              auto bv_ag = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
              auto bv_aa = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
              auto bv_as = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);

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

              for (micm::Index k = 0; k < n_r_eff_i; ++k)
              {
                auto bv_r_gas = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
                auto bv_r_aq = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
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

              for (micm::Index k = 0; k < n_N_i; ++k)
              {
                auto bv_N_gas = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
                auto bv_N_aq = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
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

              for (micm::Index k = 0; k < n_phi_i; ++k)
              {
                auto bv_phi_gas = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
                auto bv_phi_aq = jac_view.GetBlockView(views.jac_flat_ids_[idx++]);
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
    static std::size_t LookupSpecies(
        const std::unordered_map<std::string, std::size_t>& state_variable_indices, const std::string& species_name)
    {
      auto it = state_variable_indices.find(species_name);
      if (it == state_variable_indices.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_VARIABLE,
            "HenrysLawPhaseTransferSet: state variable " + species_name + " not found");
      return it->second;
    }

    Vector<micm::Index> aq_species_indices_{};
    Vector<micm::Index> solvent_species_indices_{};
    Vector<micm::Index> hlc_param_indices_{};
    Vector<micm::Index> temperature_param_indices_{};
    Vector<micm::Index> jac_flat_id_offsets_{};
    Vector<micm::Index> jac_flat_ids_{};
    Views views_{};
    CondensationRateProvider cond_rate_provider_{};

    // Host-only bookkeeping for per-instance descriptor lookup and buffer allocation.
    std::vector<std::string> prefixes_{};
    std::vector<std::size_t> n_r_eff_deps_{};
    std::vector<std::size_t> n_N_deps_{};
    std::vector<std::size_t> n_phi_deps_{};
  };
}  // namespace miam
