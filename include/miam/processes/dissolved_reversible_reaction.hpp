// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/representations/aerosol_property.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
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
  ///          where \f$ c_{H_2O} = 55.556 \f$ mol/L.
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
    std::map<std::string, std::function<double(const micm::Conditions& conditions)>>
        forward_rate_constants_;  ///< Forward rate constant functions keyed by representation prefix
    std::map<std::string, std::function<double(const micm::Conditions& conditions)>>
        reverse_rate_constants_;             ///< Reverse rate constant functions keyed by representation prefix
    std::vector<micm::Species> reactants_;   ///< Reactant species
    std::vector<micm::Species> products_;    ///< Product species
    micm::Species solvent_;                  ///< Solvent species
    micm::Phase phase_;                      ///< Phase in which the reaction occurs
    std::string uuid_;                       ///< Unique identifier for the reaction
    double solvent_floor_{
      1.0e-20
    };  ///< Floor [mol m⁻³] added to [S] in ([S]+δ)^n denominator to prevent singularity as [S] → 0

    DissolvedReversibleReaction() = delete;

    /// @brief Constructor
    DissolvedReversibleReaction(
        std::map<std::string, std::function<double(const micm::Conditions& conditions)>> forward_rate_constants,
        std::map<std::string, std::function<double(const micm::Conditions& conditions)>> reverse_rate_constants,
        const std::vector<micm::Species>& reactants,
        const std::vector<micm::Species>& products,
        micm::Species solvent,
        micm::Phase phase,
        double solvent_floor = 1.0e-20)
        : forward_rate_constants_(std::move(forward_rate_constants)),
          reverse_rate_constants_(std::move(reverse_rate_constants)),
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
          forward_rate_constants_, reverse_rate_constants_, reactants_, products_, solvent_, phase_, solvent_floor_);
    }

    /// @brief Returns a set of unique parameter names for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @return Set of unique parameter names for this process
    std::set<std::string> ProcessParameterNames(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      // One forward and one reverse rate constant parameter per representation instance of this phase.
      std::set<std::string> parameter_names;
      auto it = phase_prefixes.find(phase_.name_);
      if (it != phase_prefixes.end())
      {
        for (const auto& prefix : it->second)
        {
          parameter_names.insert(prefix + "." + phase_.name_ + "." + uuid_ + ".k_forward");
          parameter_names.insert(prefix + "." + phase_.name_ + "." + uuid_ + ".k_reverse");
        }
      }
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
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      // Build per-prefix parameter slots paired with the matching rate constant functions
      std::vector<std::pair<std::size_t, std::function<double(const micm::Conditions&)>>> forward_slots;
      std::vector<std::pair<std::size_t, std::function<double(const micm::Conditions&)>>> reverse_slots;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: UpdateStateParametersFunction: Phase " + phase_.name_ + " not found in phase_prefixes");
      for (const auto& prefix : phase_it->second)
      {
        std::string forward_param = prefix + "." + phase_.name_ + "." + uuid_ + ".k_forward";
        std::string reverse_param = prefix + "." + phase_.name_ + "." + uuid_ + ".k_reverse";
        if (state_parameter_indices.find(forward_param) == state_parameter_indices.end())
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_STATE_PARAMETER,
              "Internal Error: UpdateStateParametersFunction: Forward rate constant parameter " + forward_param +
                  " not found in state_parameter_indices");
        if (state_parameter_indices.find(reverse_param) == state_parameter_indices.end())
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_STATE_PARAMETER,
              "Internal Error: UpdateStateParametersFunction: Reverse rate constant parameter " + reverse_param +
                  " not found in state_parameter_indices");
        auto forward_it = forward_rate_constants_.find(prefix);
        if (forward_it == forward_rate_constants_.end())
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
              "DissolvedReversibleReaction: No forward rate constant configured for representation prefix '" + prefix +
                  "'");
        auto reverse_it = reverse_rate_constants_.find(prefix);
        if (reverse_it == reverse_rate_constants_.end())
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
              "DissolvedReversibleReaction: No reverse rate constant configured for representation prefix '" + prefix +
                  "'");
        forward_slots.push_back({ state_parameter_indices.at(forward_param), forward_it->second });
        reverse_slots.push_back({ state_parameter_indices.at(reverse_param), reverse_it->second });
      }

      // Set up dummy arguments to build the function
      DenseMatrixPolicy state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      std::vector<micm::Conditions> conditions_vector;

      // return a function that updates the forward and reverse rate constant parameters based on the current conditions
      return DenseMatrixPolicy::Function(
          [forward_slots, reverse_slots](auto&& conditions, auto&& params)
          {
            for (const auto& [param_index, rate_fn] : forward_slots)
              params.ForEachRow(
                  [&](const micm::Conditions& condition, double& parameter) { parameter = rate_fn(condition); },
                  conditions,
                  params.GetColumnView(param_index));
            for (const auto& [param_index, rate_fn] : reverse_slots)
              params.ForEachRow(
                  [&](const micm::Conditions& condition, double& parameter) { parameter = rate_fn(condition); },
                  conditions,
                  params.GetColumnView(param_index));
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
      auto [forward_indices, reverse_indices] = GetParameterIndices(phase_prefixes, state_parameter_indices);
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      return DenseMatrixPolicy::Function(
          [this, variable_indices, forward_indices, reverse_indices](
              auto&& state_parameters, auto&& state_variables, auto&& forcing_terms)
          {
            auto forward_rate = forcing_terms.GetRowVariable();
            auto reverse_rate = forcing_terms.GetRowVariable();
            const double eps = solvent_floor_;
            const std::size_t n_r = reactants_.size();
            const std::size_t n_p = products_.size();

            // For each phase instance, calculate the reaction rate and update the forcing terms for reactants and products
            for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
            {
              // Calculate the damped forward and reverse rates
              state_parameters.ForEachRow(
                  [&](const double& forward_rate_constant,
                      const double& reverse_rate_constant,
                      const double& solvent,
                      double& forward_rate,
                      double& reverse_rate)
                  {
                    forward_rate = forward_rate_constant * solvent / std::pow(solvent + eps, n_r);
                    reverse_rate = reverse_rate_constant * solvent / std::pow(solvent + eps, n_p);
                  },
                  state_parameters.GetConstColumnView(forward_indices[i_phase]),
                  state_parameters.GetConstColumnView(reverse_indices[i_phase]),
                  state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                  forward_rate,
                  reverse_rate);
              for (std::size_t r = 0; r < reactants_.size(); ++r)
              {
                state_variables.ForEachRow(
                    [&](const double& reactant, double& forward_rate) { forward_rate *= reactant; },
                    state_variables.GetConstColumnView(variable_indices.reactant_indices_[i_phase][r]),
                    forward_rate);
              }
              for (std::size_t p = 0; p < products_.size(); ++p)
              {
                state_variables.ForEachRow(
                    [&](const double& product, double& reverse_rate) { reverse_rate *= product; },
                    state_variables.GetConstColumnView(variable_indices.product_indices_[i_phase][p]),
                    reverse_rate);
              }

              // Apply the reaction rates to the forcing terms for reactants and products
              for (std::size_t r = 0; r < reactants_.size(); ++r)
              {
                state_variables.ForEachRow(
                    [&](const double& forward_rate, const double& reverse_rate, double& forcing)
                    {
                      forcing -= forward_rate;
                      forcing += reverse_rate;
                    },
                    forward_rate,
                    reverse_rate,
                    forcing_terms.GetColumnView(variable_indices.reactant_indices_[i_phase][r]));
              }
              for (std::size_t p = 0; p < products_.size(); ++p)
              {
                state_variables.ForEachRow(
                    [&](const double& forward_rate, const double& reverse_rate, double& forcing)
                    {
                      forcing += forward_rate;
                      forcing -= reverse_rate;
                    },
                    forward_rate,
                    reverse_rate,
                    forcing_terms.GetColumnView(variable_indices.product_indices_[i_phase][p]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          dummy_state_variables);
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
      auto [forward_indices, reverse_indices] = GetParameterIndices(phase_prefixes, state_parameter_indices);
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      return SparseMatrixPolicy::Function(
          [this, variable_indices, jacobian_indices, forward_indices, reverse_indices](
              auto&& state_parameters, auto&& state_variables, auto&& jacobian_values)
          {
            auto d_forward_rate_d_ind = jacobian_values.GetBlockVariable();
            auto d_reverse_rate_d_ind = jacobian_values.GetBlockVariable();
            auto jac_id = jacobian_indices.indices_.AsVector().begin();
            const double eps = solvent_floor_;
            const std::size_t n_r = reactants_.size();
            const std::size_t n_p = products_.size();

            // For each phase instance, calculate the partial derivatives for the Jacobian entries
            for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
            {
              // Calculate partials for independent reactants
              for (std::size_t i_ind = 0; i_ind < reactants_.size(); ++i_ind)
              {
                // dr_fwd/d[R_i] = k_f * [S] / ([S]+eps)^n_r * prod(R_j, j!=i)
                jacobian_values.ForEachBlock(
                    [&](const double& forward_rate_constant, const double& solvent, double& partial)
                    { partial = forward_rate_constant * solvent / std::pow(solvent + eps, n_r); },
                    state_parameters.GetConstColumnView(forward_indices[i_phase]),
                    state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                    d_forward_rate_d_ind);
                // add contributions to the partial from the other reactants
                for (std::size_t r = 0; r < reactants_.size(); ++r)
                {
                  if (r == i_ind)
                    continue;  // Skip the variable we're taking the derivative with respect to
                  jacobian_values.ForEachBlock(
                      [&](const double& reactant, double& partial) { partial *= reactant; },
                      state_variables.GetConstColumnView(variable_indices.reactant_indices_[i_phase][r]),
                      d_forward_rate_d_ind);
                }
                // apply partial to dependent reactants (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& partial, double& jacobian) { jacobian += partial; },
                      d_forward_rate_d_ind,
                      jacobian_values.GetBlockView(*jac_id++));
                }
                // apply partial to dependent products (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& partial, double& jacobian) { jacobian -= partial; },
                      d_forward_rate_d_ind,
                      jacobian_values.GetBlockView(*jac_id++));
                }
              }
              // Calculate partials for independent products
              for (std::size_t i_ind = 0; i_ind < products_.size(); ++i_ind)
              {
                // dr_rev/d[P_i] = k_r * [S] / ([S]+eps)^n_p * prod(P_j, j!=i)
                jacobian_values.ForEachBlock(
                    [&](const double& reverse_rate_constant, const double& solvent, double& partial)
                    { partial = reverse_rate_constant * solvent / std::pow(solvent + eps, n_p); },
                    state_parameters.GetConstColumnView(reverse_indices[i_phase]),
                    state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                    d_reverse_rate_d_ind);
                // add contributions to the partial from the other products
                for (std::size_t p = 0; p < products_.size(); ++p)
                {
                  if (p == i_ind)
                    continue;  // Skip the variable we're taking the derivative with respect to
                  jacobian_values.ForEachBlock(
                      [&](const double& product, double& partial) { partial *= product; },
                      state_variables.GetConstColumnView(variable_indices.product_indices_[i_phase][p]),
                      d_reverse_rate_d_ind);
                }
                // apply partial to dependent reactants (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& partial, double& jacobian) { jacobian -= partial; },
                      d_reverse_rate_d_ind,
                      jacobian_values.GetBlockView(*jac_id++));
                }
                // apply partial to dependent products (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& partial, double& jacobian) { jacobian += partial; },
                      d_reverse_rate_d_ind,
                      jacobian_values.GetBlockView(*jac_id++));
                }
              }
              // Calculate partials for independent solvent
              // dr/d[S] = k * (eps + (1-n)*[S]) / ([S]+eps)^(n+1) * prod([species])
              jacobian_values.ForEachBlock(
                  [&](const double& forward_rate_constant,
                      const double& reverse_rate_constant,
                      const double& solvent,
                      double& forward_partial,
                      double& reverse_partial)
                  {
                    forward_partial = forward_rate_constant * (eps + (1.0 - static_cast<int>(n_r)) * solvent) /
                                      std::pow(solvent + eps, n_r + 1);
                    reverse_partial = reverse_rate_constant * (eps + (1.0 - static_cast<int>(n_p)) * solvent) /
                                      std::pow(solvent + eps, n_p + 1);
                  },
                  state_parameters.GetConstColumnView(forward_indices[i_phase]),
                  state_parameters.GetConstColumnView(reverse_indices[i_phase]),
                  state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                  d_forward_rate_d_ind,
                  d_reverse_rate_d_ind);
              // add contributions to the partial from the reactants/products
              for (std::size_t r = 0; r < reactants_.size(); ++r)
              {
                jacobian_values.ForEachBlock(
                    [&](const double& reactant, double& forward_partial) { forward_partial *= reactant; },
                    state_variables.GetConstColumnView(variable_indices.reactant_indices_[i_phase][r]),
                    d_forward_rate_d_ind);
              }
              for (std::size_t p = 0; p < products_.size(); ++p)
              {
                jacobian_values.ForEachBlock(
                    [&](const double& product, double& reverse_partial) { reverse_partial *= product; },
                    state_variables.GetConstColumnView(variable_indices.product_indices_[i_phase][p]),
                    d_reverse_rate_d_ind);
              }
              // apply partials to dependent reactants (subtract: -J convention)
              for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
              {
                jacobian_values.ForEachBlock(
                    [&](const double& forward_partial, const double& reverse_partial, double& jacobian)
                    {
                      jacobian += forward_partial;
                      jacobian -= reverse_partial;
                    },
                    d_forward_rate_d_ind,
                    d_reverse_rate_d_ind,
                    jacobian_values.GetBlockView(*jac_id++));
              }
              // apply partials to dependent products (subtract: -J convention)
              for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
              {
                jacobian_values.ForEachBlock(
                    [&](const double& forward_partial, const double& reverse_partial, double& jacobian)
                    {
                      jacobian -= forward_partial;
                      jacobian += reverse_partial;
                    },
                    d_forward_rate_d_ind,
                    d_reverse_rate_d_ind,
                    jacobian_values.GetBlockView(*jac_id++));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          jacobian);
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

    /// @brief Helper function to return parameter indices for the forward and reverse rate constants
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_parameter_indices Map of state parameter names to their corresponding indices in the state parameter
    /// vector
    /// @return Pair of (forward, reverse) parameter index vectors, one entry per phase instance in prefix-sorted order
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>> GetParameterIndices(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      std::vector<std::size_t> forward_indices;
      std::vector<std::size_t> reverse_indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: GetParameterIndices: Phase " + phase_.name_ + " not found in phase_prefixes");
      for (const auto& prefix : phase_it->second)
      {
        std::string forward_param = prefix + "." + phase_.name_ + "." + uuid_ + ".k_forward";
        std::string reverse_param = prefix + "." + phase_.name_ + "." + uuid_ + ".k_reverse";
        if (state_parameter_indices.find(forward_param) == state_parameter_indices.end())
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_STATE_PARAMETER,
              "Internal Error: GetParameterIndices: Forward rate constant parameter " + forward_param +
                  " not found in state_parameter_indices");
        if (state_parameter_indices.find(reverse_param) == state_parameter_indices.end())
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_STATE_PARAMETER,
              "Internal Error: GetParameterIndices: Reverse rate constant parameter " + reverse_param +
                  " not found in state_parameter_indices");
        forward_indices.push_back(state_parameter_indices.at(forward_param));
        reverse_indices.push_back(state_parameter_indices.at(reverse_param));
      }
      return { forward_indices, reverse_indices };
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
