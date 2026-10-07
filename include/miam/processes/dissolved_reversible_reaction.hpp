// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/rate_constant.hpp>
#include <miam/representations/aerosol_property.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/types.hpp>

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace miam
{
  /// @brief A dissolved reversible reaction
  /// @details Dissolved reversible reactions involve reactants and products in solution, and
  ///          are characterized by both a forward and reverse rate constant. If an equilibrium
  ///          constant is provided, only one of the forward or reverse rate constants needs to be
  ///          specified, and the other is computed from the equilibrium constant.
  ///
  ///          The reaction takes the form:
  ///          \f$ \mathrm{Reactants} \leftrightarrow \mathrm{Products} \f$
  ///
  ///          with forward rate constant \f$ k_f \f$ and reverse rate constant \f$ k_r \f$. The
  ///          relationship between the rate constants and the equilibrium constant \f$ K_{eq} \f$ is:
  ///          \f$ K_{eq} = \frac{k_f}{k_r} = \frac{\prod[\mathrm{Products}]}{\prod[\mathrm{Reactants}]} \f$.
  ///
  ///          Both \f$ k_f \f$ and \f$ k_r \f$ are in units of s⁻¹ after the
  ///          solvent-normalization conversion from literature molar units:
  ///          \f$ k_f = k_{f,lit} \times c_{H_2O}^{n_r - 1} \f$,
  ///          \f$ k_r = k_{r,lit} \times c_{H_2O}^{n_p - 1} \f$
  ///          where \f$ c_{H_2O} = 55.51 \f$ mol/L (1000 g/L ÷ 18.015 g/mol).
  ///          \f$ K_{eq} \f$ is dimensionless:
  ///          \f$ K_{eq} = K_{lit} / c_{H_2O}^{n_p - n_r} \f$.
  ///
  ///          The forward and reverse rate constants are stored per representation prefix, so the
  ///          same reaction can carry different kinetics in different aerosol representations (for
  ///          example, droplet-water-molarity-dependent rates). Every representation holding the
  ///          reaction's phase must have a forward and reverse rate constant configured.
  class DissolvedReversibleReaction
  {
   public:
    // Exactly two of the three constants are set. The third follows from k_f = K_eq * k_r.
    std::optional<RateConstant> forward_rate_constant_;         ///< Forward rate constant
    std::optional<RateConstant> reverse_rate_constant_;         ///< Reverse rate constant
    std::optional<EquilibriumConstant> equilibrium_constant_;  ///< Equilibrium constant
    std::vector<micm::Species> reactants_;                                             ///< Reactant species
    std::vector<micm::Species> products_;                                              ///< Product species
    micm::Species solvent_;                                                            ///< Solvent species
    micm::Phase phase_;  ///< Phase in which the reaction occurs
    std::string uuid_;   ///< Unique identifier for the reaction
    double solvent_floor_{
      1.0e-20
    };  ///< Floor [mol m⁻³] added to [S] in ([S]+δ)^n denominator to prevent singularity as [S] → 0

    DissolvedReversibleReaction() = delete;

    /// @brief Constructor
    DissolvedReversibleReaction(
        std::optional<RateConstant> forward_rate_constant,
        std::optional<RateConstant> reverse_rate_constant,
        const std::vector<micm::Species>& reactants,
        const std::vector<micm::Species>& products,
        micm::Species solvent,
        micm::Phase phase,
        double solvent_floor = 1.0e-20,
        std::optional<EquilibriumConstant> equilibrium_constant = std::nullopt)
        : forward_rate_constant_(forward_rate_constant),
          reverse_rate_constant_(reverse_rate_constant),
          equilibrium_constant_(equilibrium_constant),
          reactants_(reactants),
          products_(products),
          solvent_(solvent),
          phase_(phase),
          uuid_(GenerateUuid()),
          solvent_floor_(solvent_floor)
    {
    }

    /// @brief Create a copy of this reaction with a new UUID
    /// @return A new DissolvedReversibleReaction with the same properties but a unique UUID
    DissolvedReversibleReaction CopyWithNewUuid() const
    {
      return DissolvedReversibleReaction(
          forward_rate_constant_,
          reverse_rate_constant_,
          reactants_,
          products_,
          solvent_,
          phase_,
          solvent_floor_,
          equilibrium_constant_);
    }

    /// @brief Returns a set of unique parameter names for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @return Set of unique parameter names for this process
    std::set<std::string> ProcessParameterNames(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      // The conditions are shared by the whole system, so we just need one value each for the forward and reverse rate
      // constants. We can use the phase name and uuid to create unique parameter names.
      std::set<std::string> parameter_names;
      parameter_names.insert(phase_.name_ + "." + uuid_ + ".k_forward");
      parameter_names.insert(phase_.name_ + "." + uuid_ + ".k_reverse");
      return parameter_names;
    }

    /// @brief Returns participating species' unique state names
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @return Set of unique state variable names for all species involved in the reaction
    std::set<std::string> SpeciesUsed(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> species_names;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        const auto& prefixes = phase_it->second;
        for (const auto& prefix : prefixes)
        {
          for (const auto& reactant : reactants_)
          {
            species_names.insert(prefix + "." + phase_.name_ + "." + reactant.name_);
          }
          for (const auto& product : products_)
          {
            species_names.insert(prefix + "." + phase_.name_ + "." + product.name_);
          }
          species_names.insert(prefix + "." + phase_.name_ + "." + solvent_.name_);
        }
      }
      else
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: Phase " + phase_.name_ + " not found in phase_prefixes map for process " + uuid_);
      }
      return species_names;
    }

    /// @brief Returns the aerosol properties required by this process
    /// @details DissolvedReversibleReaction does not depend on aerosol properties.
    /// @return Empty map
    std::map<std::string, std::vector<AerosolProperty>> RequiredAerosolProperties() const
    {
      return {};
    }

    /// @brief Returns a set of Jacobian index pairs for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_variable_indices Map of state variable names to their corresponding indices in the Jacobian
    /// @return Set of index pairs representing the Jacobian entries affected by this process (dependent, independent)
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_variable_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      std::set<std::pair<std::size_t, std::size_t>> jacobian_indices;

      StateVariableIndices variable_indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      // Get pairs for each phase instance
      for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
      {
        // Reactant contributions
        for (std::size_t r = 0; r < variable_indices.reactant_indices_.NumColumns(); ++r)
        {
          std::size_t independent_index = variable_indices.reactant_indices_[i_phase][r];
          // Each reactant affects all reactants
          for (std::size_t r2 = 0; r2 < variable_indices.reactant_indices_.NumColumns(); ++r2)
          {
            std::size_t dependent_index = variable_indices.reactant_indices_[i_phase][r2];
            jacobian_indices.insert({ dependent_index, independent_index });
          }
          // Each reactant affects all products
          for (std::size_t p = 0; p < variable_indices.product_indices_.NumColumns(); ++p)
          {
            std::size_t dependent_index = variable_indices.product_indices_[i_phase][p];
            jacobian_indices.insert({ dependent_index, independent_index });
          }
        }
        // Product contributions
        for (std::size_t p = 0; p < variable_indices.product_indices_.NumColumns(); ++p)
        {
          std::size_t independent_index = variable_indices.product_indices_[i_phase][p];
          // Each product affects all reactants
          for (std::size_t r = 0; r < variable_indices.reactant_indices_.NumColumns(); ++r)
          {
            std::size_t dependent_index = variable_indices.reactant_indices_[i_phase][r];
            jacobian_indices.insert({ dependent_index, independent_index });
          }
          // Each product affects all products
          for (std::size_t p2 = 0; p2 < variable_indices.product_indices_.NumColumns(); ++p2)
          {
            std::size_t dependent_index = variable_indices.product_indices_[i_phase][p2];
            jacobian_indices.insert({ dependent_index, independent_index });
          }
        }
        // Solvent contributions (affects all reactants and products)
        std::size_t independent_index = variable_indices.solvent_indices_[i_phase];
        for (std::size_t r = 0; r < variable_indices.reactant_indices_.NumColumns(); ++r)
        {
          std::size_t dependent_index = variable_indices.reactant_indices_[i_phase][r];
          jacobian_indices.insert({ dependent_index, independent_index });
        }
        for (std::size_t p = 0; p < variable_indices.product_indices_.NumColumns(); ++p)
        {
          std::size_t dependent_index = variable_indices.product_indices_[i_phase][p];
          jacobian_indices.insert({ dependent_index, independent_index });
        }
      }

      return jacobian_indices;
    }

    /// @brief Returns non-zero Jacobian elements (common interface overload accepting providers)
    /// @details Delegates to the existing two-argument version; providers are unused.
    template<typename DenseMatrixPolicy>
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>>& /* providers */)
        const
    {
      return NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    }

    /// @brief Returns a function that updates state parameters for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_parameter_indices Map of state parameter names to their corresponding indices in the state parameter
    /// vector
    /// @return Function that updates state parameters for this process
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateStateParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      // throw an error if the expected parameters don't exist
      std::string forward_param = phase_.name_ + "." + uuid_ + ".k_forward";
      std::string reverse_param = phase_.name_ + "." + uuid_ + ".k_reverse";
      if (state_parameter_indices.find(forward_param) == state_parameter_indices.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "Internal Error: UpdateStateParametersFunction: Forward rate constant parameter " + forward_param +
                " not found in state_parameter_indices");
      }
      if (state_parameter_indices.find(reverse_param) == state_parameter_indices.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "Internal Error: UpdateStateParametersFunction: Reverse rate constant parameter " + reverse_param +
                " not found in state_parameter_indices");
      }
      std::size_t forward_index = state_parameter_indices.at(forward_param);
      std::size_t reverse_index = state_parameter_indices.at(reverse_param);

      std::size_t num_params = state_parameter_indices.size();

      // A rate constant that is not given is a placeholder. The kernel computes it from K_eq.
      const RateConstant forward = forward_rate_constant_.value_or(0.0);
      const RateConstant reverse = reverse_rate_constant_.value_or(0.0);
      const EquilibriumConstant equilibrium = equilibrium_constant_.value_or(EquilibriumConstant{});
      const bool has_forward = forward_rate_constant_.has_value();
      const bool has_reverse = reverse_rate_constant_.has_value();

      // NVCC forbids an extended __host__ __device__ lambda inside the generic std::visit lambda.
      return std::visit(
          [&](const auto& forward_form, const auto& reverse_form)
          {
            return MakeReversibleUpdateFn<DenseMatrixPolicy>(
                forward_form, reverse_form, equilibrium, has_forward, has_reverse, forward_index, reverse_index, num_params);
          },
          forward,
          reverse);
    }

    /// @brief Builds the forward and reverse rate constant update function for one pair of rate constant types
    /// @details NVCC requires this function to be public.
    template<typename DenseMatrixPolicy, typename ForwardT, typename ReverseT>
    static std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    MakeReversibleUpdateFn(
        const ForwardT& forward,
        const ReverseT& reverse,
        const EquilibriumConstant& equilibrium,
        bool has_forward,
        bool has_reverse,
        std::size_t forward_index,
        std::size_t reverse_index,
        std::size_t num_params)
    {
      // Set up dummy arguments to build the function
      DenseMatrixPolicy state_parameters{ 1, num_params, 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> conditions_vector;

      // return a function that updates the forward and reverse rate constant parameters based on the current conditions
      return DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
              const typename DenseMatrixPolicy::ViewType& params_view)
          {
            params_view.ForEachRowStrict(
                [forward, reverse, equilibrium, has_forward, has_reverse](
                    const micm::Conditions& condition, micm::Real& k_forward, micm::Real& k_reverse)
                {
                  k_forward = has_forward ? Calculate(forward, condition)
                                          : Calculate(equilibrium, condition) * Calculate(reverse, condition);
                  k_reverse = has_reverse ? Calculate(reverse, condition)
                                          : Calculate(forward, condition) / Calculate(equilibrium, condition);
                },
                conditions_view,
                params_view.GetColumnView(forward_index),
                params_view.GetColumnView(reverse_index));
          },
          conditions_vector,
          state_parameters);
    }

    /// @brief Returns a function that calculates the forcing terms for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_parameter_indices Map of state parameter names to their corresponding indices in the state parameter
    /// vector
    /// @param state_variable_indices Map of state variable names to their corresponding indices in the state variable
    /// vector
    /// @return Function that calculates the forcing terms for this process
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices,   // acts like std::unordered_map<std::string, std::size_t>
        std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> /* providers */
    ) const
    {
      return ForcingFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
    }

    /// @brief Returns a function that calculates the forcing terms for this process (original)
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices    // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      StateVariableIndices variable_indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      auto [forward_index, reverse_index] = GetParameterIndices(state_parameter_indices);
      auto storage = std::make_shared<IndexStorage<DenseMatrixPolicy>>(variable_indices);
      auto reactant_view = storage->reactant_indices_.GetView();
      auto product_view = storage->product_indices_.GetView();
      auto solvent_view = storage->solvent_indices_.GetView();
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      const std::size_t num_phases = variable_indices.number_of_phase_instances_;
      const std::size_t num_reactants = reactants_.size();
      const std::size_t num_products = products_.size();
      const double eps = solvent_floor_;
      const std::size_t k_fwd = forward_index;
      const std::size_t k_rev = reverse_index;

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ViewType& forcing_view)
          {
              // For each phase instance, calculate the reaction rate and update the forcing terms for reactants and
              // products
              for (std::size_t phase = 0; phase < num_phases; ++phase)
            {
                const std::size_t solvent_idx = solvent_view[phase];
                auto forward_rate = forcing_view.GetRowVariable();
                auto reverse_rate = forcing_view.GetRowVariable();
              // Calculate the damped forward and reverse rates
                forcing_view.ForEachRowStrict(
                    [num_reactants, num_products, eps](
                        const micm::Real& k_f,
                        const micm::Real& k_r,
                        const micm::Real& solvent,
                        micm::Real& fwd,
                        micm::Real& rev)
                  {
                      fwd = k_f * solvent / std::pow(solvent + eps, num_reactants);
                      rev = k_r * solvent / std::pow(solvent + eps, num_products);
                  },
                    params_view.GetConstColumnView(k_fwd),
                    params_view.GetConstColumnView(k_rev),
                    state_view.GetConstColumnView(solvent_idx),
                  forward_rate,
                  reverse_rate);
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& reactant, micm::Real& fwd) { fwd *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    forward_rate);
              }
                for (std::size_t p = 0; p < num_products; ++p)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& product, micm::Real& rev) { rev *= product; },
                      state_view.GetConstColumnView(product_view[phase * num_products + p]),
                    reverse_rate);
              }

              // Apply the reaction rates to the forcing terms for reactants and products
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& fwd, const micm::Real& rev, micm::Real& forcing)
                    {
                        forcing -= fwd;
                        forcing += rev;
                    },
                    forward_rate,
                    reverse_rate,
                      forcing_view.GetColumnView(reactant_view[phase * num_reactants + r]));
              }
                for (std::size_t p = 0; p < num_products; ++p)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& fwd, const micm::Real& rev, micm::Real& forcing)
                    {
                        forcing += fwd;
                        forcing -= rev;
                    },
                    forward_rate,
                    reverse_rate,
                      forcing_view.GetColumnView(product_view[phase * num_products + p]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          dummy_state_variables);

      return [storage, function](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms) mutable { function(state_parameters, state_variables, forcing_terms); };
    }

    /// @brief Returns a function that calculates the Jacobian contributions for this process (common interface overload)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices,   // acts like std::unordered_map<std::string, std::size_t>
        const SparseMatrixPolicy& jacobian,
        std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> /* providers */) const
    {
      return JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
          phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    }

    /// @brief Returns a function that calculates the Jacobian contributions for this process (original)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices,   // acts like std::unordered_map<std::string, std::size_t>
        const SparseMatrixPolicy& jacobian) const
    {
      StateVariableIndices variable_indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      JacobianIndices jacobian_indices = GetJacobianIndices(variable_indices, jacobian);
      auto [forward_index, reverse_index] = GetParameterIndices(state_parameter_indices);
      auto storage = std::make_shared<IndexStorage<SparseMatrixPolicy>>(variable_indices, jacobian_indices);
      auto reactant_view = storage->reactant_indices_.GetView();
      auto product_view = storage->product_indices_.GetView();
      auto solvent_view = storage->solvent_indices_.GetView();
      auto jac_id_view = storage->jacobian_flat_ids_.GetView();
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      const std::size_t num_phases = variable_indices.number_of_phase_instances_;
      const std::size_t num_reactants = reactants_.size();
      const std::size_t num_products = products_.size();
      const double eps = solvent_floor_;
      const std::size_t k_fwd = forward_index;
      const std::size_t k_rev = reverse_index;

      auto function = SparseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename SparseMatrixPolicy::ViewType& jac_view)
          {
              const std::size_t pairs_per_phase = (num_reactants + num_products + 1) * (num_reactants + num_products);

            // For each phase instance, calculate the partial derivatives for the Jacobian entries
              for (std::size_t phase = 0; phase < num_phases; ++phase)
            {
                const std::size_t solvent_idx = solvent_view[phase];
                auto d_forward_rate_d_ind = jac_view.GetBlockVariable();
                auto d_reverse_rate_d_ind = jac_view.GetBlockVariable();
                std::size_t pair = phase * pairs_per_phase;

              // Calculate partials for independent reactants
                for (std::size_t i_ind = 0; i_ind < num_reactants; ++i_ind)
              {
                // dr_fwd/d[R_i] = k_f * [S] / ([S]+eps)^n_r * prod(R_j, j!=i)
                  jac_view.ForEachBlockStrict(
                      [num_reactants, eps](const micm::Real& k_f, const micm::Real& solvent, micm::Real& partial)
                      { partial = k_f * solvent / std::pow(solvent + eps, num_reactants); },
                      params_view.GetConstColumnView(k_fwd),
                      state_view.GetConstColumnView(solvent_idx),
                    d_forward_rate_d_ind);
                // add contributions to the partial from the other reactants
                  for (std::size_t r = 0; r < num_reactants; ++r)
                {
                  if (r == i_ind)
                    continue;  // Skip the variable we're taking the derivative with respect to
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                        state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                      d_forward_rate_d_ind);
                }
                // apply partial to dependent reactants (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                      d_forward_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
                // apply partial to dependent products (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                      d_forward_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
              }
              // Calculate partials for independent products
                for (std::size_t i_ind = 0; i_ind < num_products; ++i_ind)
              {
                // dr_rev/d[P_i] = k_r * [S] / ([S]+eps)^n_p * prod(P_j, j!=i)
                  jac_view.ForEachBlockStrict(
                      [num_products, eps](const micm::Real& k_r, const micm::Real& solvent, micm::Real& partial)
                      { partial = k_r * solvent / std::pow(solvent + eps, num_products); },
                      params_view.GetConstColumnView(k_rev),
                      state_view.GetConstColumnView(solvent_idx),
                    d_reverse_rate_d_ind);
                // add contributions to the partial from the other products
                  for (std::size_t p = 0; p < num_products; ++p)
                {
                  if (p == i_ind)
                    continue;  // Skip the variable we're taking the derivative with respect to
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& product, micm::Real& partial) { partial *= product; },
                        state_view.GetConstColumnView(product_view[phase * num_products + p]),
                      d_reverse_rate_d_ind);
                }
                // apply partial to dependent reactants (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                      d_reverse_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
                // apply partial to dependent products (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                      d_reverse_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
              }
              // Calculate partials for independent solvent
              // dr/d[S] = k * (eps + (1-n)*[S]) / ([S]+eps)^(n+1) * prod([species])
                jac_view.ForEachBlockStrict(
                    [num_reactants, num_products, eps](
                        const micm::Real& k_f,
                        const micm::Real& k_r,
                        const micm::Real& solvent,
                        micm::Real& forward_partial,
                        micm::Real& reverse_partial)
                  {
                      forward_partial = k_f * (eps + (1.0 - static_cast<micm::Real>(num_reactants)) * solvent) /
                                        std::pow(solvent + eps, num_reactants + 1);
                      reverse_partial = k_r * (eps + (1.0 - static_cast<micm::Real>(num_products)) * solvent) /
                                        std::pow(solvent + eps, num_products + 1);
                  },
                    params_view.GetConstColumnView(k_fwd),
                    params_view.GetConstColumnView(k_rev),
                    state_view.GetConstColumnView(solvent_idx),
                  d_forward_rate_d_ind,
                  d_reverse_rate_d_ind);
              // add contributions to the partial from the reactants/products
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& reactant, micm::Real& forward_partial) { forward_partial *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    d_forward_rate_d_ind);
              }
                for (std::size_t p = 0; p < num_products; ++p)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& product, micm::Real& reverse_partial) { reverse_partial *= product; },
                      state_view.GetConstColumnView(product_view[phase * num_products + p]),
                    d_reverse_rate_d_ind);
              }
              // apply partials to dependent reactants (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& forward_partial, const micm::Real& reverse_partial, micm::Real& jac)
                    {
                        jac += forward_partial;
                        jac -= reverse_partial;
                    },
                    d_forward_rate_d_ind,
                    d_reverse_rate_d_ind,
                      jac_view.GetBlockView(jac_id_view[pair++]));
              }
              // apply partials to dependent products (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& forward_partial, const micm::Real& reverse_partial, micm::Real& jac)
                    {
                        jac -= forward_partial;
                        jac += reverse_partial;
                    },
                    d_forward_rate_d_ind,
                    d_reverse_rate_d_ind,
                      jac_view.GetBlockView(jac_id_view[pair++]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          jacobian);

      return [storage, function](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian_values) mutable { function(state_parameters, state_variables, jacobian_values); };
    }

   private:
    /// @brief Helper struct for keeping track of state varible indices for reactants, products, and solvent across
    /// multiple phase instances (e.g. grid cells)
    struct StateVariableIndices
    {
      std::size_t number_of_phase_instances_;  ///< Number of instances of the phase in the system (e.g. number of grid
                                               ///< cells containing this phase)
      micm::Matrix<std::size_t>
          reactant_indices_;  ///< Matrix of state variable indices for reactants (num_reactants x num_prefixes)
      micm::Matrix<std::size_t>
          product_indices_;  ///< Matrix of state variable indices for products (num_products x num_prefixes)
      std::vector<std::size_t> solvent_indices_;  ///< Vector of state variable indices for solvent (num_prefixes)
    };

    /// @brief Helper struct for keeping track of Jacobian sparse matrix elements
    struct JacobianIndices
    {
      micm::Matrix<std::size_t>
          indices_;  // Index in sparse matrix for each dependent/independent pair (num_pairs x num_prefixes)
    };

    /// @brief Device-ready copies of the indices that the solve-time functions read
    /// @details The solve-time functions hold a shared pointer to this storage. Their kernels capture views into it.
    template<typename MatrixPolicy>
    struct IndexStorage
    {
      using Vector = typename MatrixPolicy::template VectorType<std::size_t>;
      Vector reactant_indices_;
      Vector product_indices_;
      Vector solvent_indices_;
      Vector jacobian_flat_ids_;

      IndexStorage(const StateVariableIndices& variable_indices, const JacobianIndices& jacobian_indices = {})
          : reactant_indices_(variable_indices.reactant_indices_.AsVector()),
            product_indices_(variable_indices.product_indices_.AsVector()),
            solvent_indices_(variable_indices.solvent_indices_),
            jacobian_flat_ids_(jacobian_indices.indices_.AsVector())
      {
        reactant_indices_.CopyToDevice();
        product_indices_.CopyToDevice();
        solvent_indices_.CopyToDevice();
        jacobian_flat_ids_.CopyToDevice();
      }
    };

    /// @brief Helper function to return parameter indices for the forward and reverse rate constants
    std::pair<std::size_t, std::size_t> GetParameterIndices(
        const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      std::string forward_param = phase_.name_ + "." + uuid_ + ".k_forward";
      std::string reverse_param = phase_.name_ + "." + uuid_ + ".k_reverse";
      if (state_parameter_indices.find(forward_param) == state_parameter_indices.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "Internal Error: GetParameterIndices: Forward rate constant parameter " + forward_param +
                " not found in state_parameter_indices");
      }
      if (state_parameter_indices.find(reverse_param) == state_parameter_indices.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "Internal Error: GetParameterIndices: Reverse rate constant parameter " + reverse_param +
                " not found in state_parameter_indices");
      }
      return { state_parameter_indices.at(forward_param), state_parameter_indices.at(reverse_param) };
    }

    /// @brief Helper function to return variable indices for all species involved in the reaction
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_variable_indices Map of state variable names to their corresponding indices in the state variable
    /// vector
    /// @return StateVariableIndices struct containing matrices of indices for reactants, products, and solvent
    StateVariableIndices GetStateVariableIndices(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_variable_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      StateVariableIndices indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it == phase_prefixes.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: GetStateVariableIndices: Phase " + phase_.name_ + " not found in phase_prefixes");
      }
      const auto& prefixes = phase_it->second;
      indices.number_of_phase_instances_ = prefixes.size();
      indices.reactant_indices_ = micm::Matrix<std::size_t>(prefixes.size(), reactants_.size());
      indices.product_indices_ = micm::Matrix<std::size_t>(prefixes.size(), products_.size());
      indices.solvent_indices_ = std::vector<std::size_t>(prefixes.size());
      std::size_t i_phase = 0;
      for (const auto& prefix : prefixes)
      {
        for (std::size_t i_reactant = 0; i_reactant < reactants_.size(); ++i_reactant)
        {
          std::string reactant_var = prefix + "." + phase_.name_ + "." + reactants_[i_reactant].name_;
          if (state_variable_indices.find(reactant_var) == state_variable_indices.end())
          {
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_MISSING_STATE_VARIABLE,
                "Internal Error: GetStateVariableIndices: Reactant variable " + reactant_var +
                    " not found in state_variable_indices");
          }
          indices.reactant_indices_[i_phase][i_reactant] = state_variable_indices.at(reactant_var);
        }
        for (std::size_t i_product = 0; i_product < products_.size(); ++i_product)
        {
          std::string product_var = prefix + "." + phase_.name_ + "." + products_[i_product].name_;
          if (state_variable_indices.find(product_var) == state_variable_indices.end())
          {
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_MISSING_STATE_VARIABLE,
                "Internal Error: GetStateVariableIndices: Product variable " + product_var +
                    " not found in state_variable_indices");
          }
          indices.product_indices_[i_phase][i_product] = state_variable_indices.at(product_var);
        }
        std::string solvent_var = prefix + "." + phase_.name_ + "." + solvent_.name_;
        if (state_variable_indices.find(solvent_var) == state_variable_indices.end())
        {
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_STATE_VARIABLE,
              "Internal Error: GetStateVariableIndices: Solvent variable " + solvent_var +
                  " not found in state_variable_indices");
        }
        indices.solvent_indices_[i_phase] = state_variable_indices.at(solvent_var);
        ++i_phase;
      }
      return indices;
    }

    /// @brief Helper function to return Jacobian sparse matrix indices for all pairs of species involved in the reaction
    /// @param variable_indices StateVariableIndices struct containing matrices of indices for reactants, products, and
    /// solvent
    /// @param jacobian Sparse matrix policy object for the Jacobian structure
    /// @return JacobianIndices struct containing a matrix of sparse matrix indices for each dependent/independent pair of
    /// species
    JacobianIndices GetJacobianIndices(
        const StateVariableIndices& variable_indices,
        const auto& jacobian  // sparse matrix policy object for the Jacobian structure
    ) const
    {
      // Each reactant and each product depends on all reactants, all products, and the solvent
      std::size_t num_pairs =
          (reactants_.size() + products_.size()) * (reactants_.size() + products_.size() + 1);  // +1 for solvent
      JacobianIndices jacobian_indices;
      jacobian_indices.indices_ = micm::Matrix<std::size_t>(variable_indices.number_of_phase_instances_, num_pairs);
      for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
      {
        std::size_t pair_index = 0;
        // add terms for independent reactants
        for (std::size_t i_ind = 0; i_ind < reactants_.size(); ++i_ind)
        {
          // ... and dependent reactants
          for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
          {
            jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
                0, variable_indices.reactant_indices_[i_phase][i_dep], variable_indices.reactant_indices_[i_phase][i_ind]);
          }
          // ... and dependent products
          for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
          {
            jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
                0, variable_indices.product_indices_[i_phase][i_dep], variable_indices.reactant_indices_[i_phase][i_ind]);
          }
        }
        // add terms for independent products
        for (std::size_t i_ind = 0; i_ind < products_.size(); ++i_ind)
        {
          // ... and dependent reactants
          for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
          {
            jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
                0, variable_indices.reactant_indices_[i_phase][i_dep], variable_indices.product_indices_[i_phase][i_ind]);
          }
          // ... and dependent products
          for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
          {
            jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
                0, variable_indices.product_indices_[i_phase][i_dep], variable_indices.product_indices_[i_phase][i_ind]);
          }
        }
        // add terms for independent solvent
        // ... and dependent reactants
        for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
        {
          jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
              0, variable_indices.reactant_indices_[i_phase][i_dep], variable_indices.solvent_indices_[i_phase]);
        }
        // ... and dependent products
        for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
        {
          jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
              0, variable_indices.product_indices_[i_phase][i_dep], variable_indices.solvent_indices_[i_phase]);
        }
      }
      return jacobian_indices;
    }
  };
}  // namespace miam
