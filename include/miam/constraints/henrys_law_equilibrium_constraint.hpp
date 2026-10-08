// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/math/condensation_rate.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
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
    HenrysLawConstant henrys_law_constant_;  ///< HLC(T) function [mol m⁻³ Pa⁻¹]
    micm::Species gas_species_;                                                      ///< Gas-phase species
    micm::Species condensed_species_;                                                ///< Condensed-phase solute species
    micm::Species solvent_;                                                          ///< Condensed-phase solvent species
    micm::Phase condensed_phase_;                                                    ///< The condensed phase
    double solvent_molecular_weight_;  ///< Solvent molecular weight [kg mol⁻¹]
    double solvent_density_;           ///< Solvent density [kg m⁻³]
    std::string uuid_;                 ///< Unique identifier

    HenrysLawEquilibriumConstraint() = delete;

    /// @brief Constructor
    HenrysLawEquilibriumConstraint(
        HenrysLawConstant henrys_law_constant,
        const micm::Species& gas_species,
        const micm::Species& condensed_species,
        const micm::Species& solvent,
        const micm::Phase& condensed_phase,
        double solvent_molecular_weight,
        double solvent_density)
        : henrys_law_constant_(henrys_law_constant),
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
      const HenrysLawConstant henrys_law_constant = henrys_law_constant_;

      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      auto storage = std::make_shared<Vector>(hlc_rt_indices);
      storage->CopyToDevice();
      auto hlc_rt_view = storage->GetView();
      const std::size_t num_hlc_rt = hlc_rt_indices.size();
      DenseMatrixPolicy dummy_params{ 1, state_parameter_indices.size(), 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> dummy_conditions;

      auto function = DenseMatrixPolicy::Function(
              MICM_LAMBDA(
                  const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
                  const typename DenseMatrixPolicy::ViewType& params_view)
          {
                for (std::size_t i = 0; i < num_hlc_rt; ++i)
                params_view.ForEachRowStrict(
                    [henrys_law_constant](const micm::Conditions& cond, micm::Real& hlc_rt)
                  { hlc_rt = Calculate(henrys_law_constant, cond) * micm::constants::GAS_CONSTANT * cond.temperature_; },
                    conditions_view,
                      params_view.GetColumnView(hlc_rt_view[i]));
          },
          dummy_conditions,
          dummy_params);

      return [storage, function](
                 const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
                 DenseMatrixPolicy& params) mutable { function(conditions, params); };
    }

    /// @brief Returns a function that computes constraint residuals G(y) = 0
    /// @details G = HLC * R * T * f_v * [A_g] - [A_aq]
    ///          where f_v = [S] * solvent_molecular_weight / solvent_density
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      double molar_volume = solvent_molecular_weight_ / solvent_density_;  // [m³ mol⁻¹]

      std::vector<std::size_t> hlc_rt_indices;
      auto phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          hlc_rt_indices.push_back(
              state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc_rt"));
      }

      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      Vector aq_indices(indices.aq_indices_);
      Vector solvent_indices(indices.solvent_indices_);
      Vector hlc_rt_indices_vec(hlc_rt_indices);
      aq_indices.CopyToDevice();
      solvent_indices.CopyToDevice();
      hlc_rt_indices_vec.CopyToDevice();
      std::size_t num_phases = indices.number_of_phase_instances_;
      std::size_t gas_idx = indices.gas_idx_;

      struct Storage
      {
        Vector aq_indices, solvent_indices, hlc_rt_indices_vec;
      };
      auto storage = std::make_shared<Storage>(
          Storage{ std::move(aq_indices), std::move(solvent_indices), std::move(hlc_rt_indices_vec) });
      auto aq_view = storage->aq_indices.GetView();
      auto solvent_view = storage->solvent_indices.GetView();
      auto hlc_rt_view = storage->hlc_rt_indices_vec.GetView();
      DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_params{ 1, std::max(state_parameter_indices.size(), std::size_t{ 1 }), 0.0 };

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ViewType& residual_view)
          {
              for (std::size_t i_phase = 0; i_phase < num_phases; ++i_phase)
            {
                const std::size_t aq_idx = aq_view[i_phase];
              // G = HLC*R*T * f_v * [A_g] - [A_aq],  where f_v = [S] * solvent_molecular_weight / solvent_density [m³
              // mol⁻¹]
                residual_view.ForEachRowStrict(
                    [molar_volume](
                        const micm::Real& hlc_rt, const micm::Real& gas, const micm::Real& aq, const micm::Real& sol, micm::Real& res)
                  { res = hlc_rt * (sol * molar_volume) * gas - aq; },
                    params_view.GetConstColumnView(hlc_rt_view[i_phase]),
                    state_view.GetConstColumnView(gas_idx),
                    state_view.GetConstColumnView(aq_idx),
                    state_view.GetConstColumnView(solvent_view[i_phase]),
                    residual_view.GetColumnView(aq_idx));
            }
          },
          dummy_state,
          dummy_params,
          dummy_state);

      return [storage, function](
                 const DenseMatrixPolicy& state_variables,
                 const DenseMatrixPolicy& state_parameters,
                 DenseMatrixPolicy& residual) mutable { function(state_variables, state_parameters, residual); };
    }

    /// @brief Returns a function that computes constraint Jacobian entries (subtracts dG/dy)
    /// @details Follows MICM convention: jac[row][col] -= dG/dy
    ///          dG/d[A_g] = HLC * R * T * f_v
    ///          dG/d[A_aq] = -1
    ///          dG/d[S] = HLC * R * T * (solvent_molecular_weight / solvent_density) [m³ mol⁻¹] * [A_g]
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      double molar_volume = solvent_molecular_weight_ / solvent_density_;  // [m³ mol⁻¹]

      std::vector<std::size_t> hlc_rt_indices;
      auto phase_it = phase_prefixes.find(condensed_phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          hlc_rt_indices.push_back(
              state_parameter_indices.at(prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc_rt"));
      }

      // Pre-compute block-0 VectorIndex offsets per instance
      std::vector<std::size_t> gas_jac_ids(indices.number_of_phase_instances_);
      std::vector<std::size_t> aq_jac_ids(indices.number_of_phase_instances_);
      std::vector<std::size_t> solvent_jac_ids(indices.number_of_phase_instances_);
      for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
      {
        std::size_t aq_row = indices.aq_indices_[i_phase];
        gas_jac_ids[i_phase] = jacobian.VectorIndex(0, aq_row, indices.gas_idx_);
        aq_jac_ids[i_phase] = jacobian.VectorIndex(0, aq_row, aq_row);
        solvent_jac_ids[i_phase] = jacobian.VectorIndex(0, aq_row, indices.solvent_indices_[i_phase]);
      }

      using Vector = typename SparseMatrixPolicy::template VectorType<std::size_t>;
      Vector solvent_indices(indices.solvent_indices_);
      Vector hlc_rt_indices_vec(hlc_rt_indices);
      Vector gas_jac_ids_vec(gas_jac_ids);
      Vector aq_jac_ids_vec(aq_jac_ids);
      Vector solvent_jac_ids_vec(solvent_jac_ids);
      solvent_indices.CopyToDevice();
      hlc_rt_indices_vec.CopyToDevice();
      gas_jac_ids_vec.CopyToDevice();
      aq_jac_ids_vec.CopyToDevice();
      solvent_jac_ids_vec.CopyToDevice();
      std::size_t num_phases = indices.number_of_phase_instances_;
      std::size_t gas_idx = indices.gas_idx_;

      struct Storage
      {
        Vector solvent_indices, hlc_rt_indices_vec, gas_jac_ids_vec, aq_jac_ids_vec, solvent_jac_ids_vec;
      };
      auto storage = std::make_shared<Storage>(Storage{ std::move(solvent_indices),
                                                        std::move(hlc_rt_indices_vec),
                                                        std::move(gas_jac_ids_vec),
                                                        std::move(aq_jac_ids_vec),
                                                        std::move(solvent_jac_ids_vec) });
      auto solvent_view = storage->solvent_indices.GetView();
      auto hlc_rt_view = storage->hlc_rt_indices_vec.GetView();
      auto gas_jac_view = storage->gas_jac_ids_vec.GetView();
      auto aq_jac_view = storage->aq_jac_ids_vec.GetView();
      auto solvent_jac_view = storage->solvent_jac_ids_vec.GetView();
      DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_params{ 1, std::max(state_parameter_indices.size(), std::size_t{ 1 }), 0.0 };

      auto function = SparseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename SparseMatrixPolicy::ViewType& jac_view)
          {
              for (std::size_t i_phase = 0; i_phase < num_phases; ++i_phase)
            {
                auto bv_gas = jac_view.GetBlockView(gas_jac_view[i_phase]);
                auto bv_aq = jac_view.GetBlockView(aq_jac_view[i_phase]);
                auto bv_sol = jac_view.GetBlockView(solvent_jac_view[i_phase]);

              // jac -= dG/d[A_g] = HLC*R*T * f_v
              // jac -= dG/d[A_aq] = -1
              // jac -= dG/d[S] = HLC*R*T * molar_volume [m³ mol⁻¹] * [A_g]
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
                    params_view.GetConstColumnView(hlc_rt_view[i_phase]),
                    state_view.GetConstColumnView(gas_idx),
                    state_view.GetConstColumnView(solvent_view[i_phase]),
                  bv_gas,
                  bv_aq,
                  bv_sol);
            }
          },
          dummy_state,
          dummy_params,
          jacobian);

      return [storage, function](
                 const DenseMatrixPolicy& state_variables,
                 const DenseMatrixPolicy& state_parameters,
                 SparseMatrixPolicy& jacobian_values) mutable { function(state_variables, state_parameters, jacobian_values); };
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
