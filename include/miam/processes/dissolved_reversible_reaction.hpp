// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/aerosol_property.hpp>
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
  namespace process
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
    struct DissolvedReversibleReaction
    {
      std::function<double(const micm::Conditions& conditions)> forward_rate_constant_;  ///< Forward rate constant function
      std::function<double(const micm::Conditions& conditions)> reverse_rate_constant_;  ///< Reverse rate constant function
      std::vector<micm::Species> reactants_;                                             ///< Reactant species
      std::vector<micm::Species> products_;                                              ///< Product species
      micm::Species solvent_;                                                            ///< Solvent species
      micm::Phase phase_;  ///< Phase in which the reaction occurs
      std::string uuid_;   ///< Unique identifier for the reaction

      DissolvedReversibleReaction() = delete;

      /// @brief Constructor
      DissolvedReversibleReaction(
          std::function<double(const micm::Conditions& conditions)> forward_rate_constant,
          std::function<double(const micm::Conditions& conditions)> reverse_rate_constant,
          const std::vector<micm::Species>& reactants,
          const std::vector<micm::Species>& products,
          micm::Species solvent,
          micm::Phase phase)
          : forward_rate_constant_(forward_rate_constant),
            reverse_rate_constant_(reverse_rate_constant),
            reactants_(reactants),
            products_(products),
            solvent_(solvent),
            phase_(phase),
            uuid_(generate_uuid_v4())
      {
      }

      /// @brief Create a copy of this reaction with a new UUID
      /// @return A new DissolvedReversibleReaction with the same properties but a unique UUID
      DissolvedReversibleReaction CopyWithNewUuid() const
      {
        return DissolvedReversibleReaction(
            forward_rate_constant_, reverse_rate_constant_, reactants_, products_, solvent_, phase_);
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
          throw std::runtime_error(
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
          const std::map<std::string, std::vector<AerosolPropertyProvider<DenseMatrixPolicy>>>& /* providers */) const
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
        // throw an error if the expected parameters don't exist
        std::string forward_param = phase_.name_ + "." + uuid_ + ".k_forward";
        std::string reverse_param = phase_.name_ + "." + uuid_ + ".k_reverse";
        if (state_parameter_indices.find(forward_param) == state_parameter_indices.end())
        {
          throw std::runtime_error(
              "Internal Error: UpdateStateParametersFunction: Forward rate constant parameter " + forward_param +
              " not found in state_parameter_indices");
        }
        if (state_parameter_indices.find(reverse_param) == state_parameter_indices.end())
        {
          throw std::runtime_error(
              "Internal Error: UpdateStateParametersFunction: Reverse rate constant parameter " + reverse_param +
              " not found in state_parameter_indices");
        }
        std::size_t forward_index = state_parameter_indices.at(forward_param);
        std::size_t reverse_index = state_parameter_indices.at(reverse_param);

        // Set up dummy arguments to build the function
        DenseMatrixPolicy state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        std::vector<micm::Conditions> conditions_vector;

        // return a function that updates the forward and reverse rate constant parameters based on the current conditions
        return DenseMatrixPolicy::Function(
            [this, forward_index, reverse_index](auto&& conditions, auto&& params)
            {
              params.ForEachRow(
                  [&](const micm::Conditions& condition, double& parameter)
                  { parameter = forward_rate_constant_(condition); },
                  conditions,
                  params.GetColumnView(forward_index));
              params.ForEachRow(
                  [&](const micm::Conditions& condition, double& parameter)
                  { parameter = reverse_rate_constant_(condition); },
                  conditions,
                  params.GetColumnView(reverse_index));
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
          std::map<std::string, std::vector<AerosolPropertyProvider<DenseMatrixPolicy>>> /* providers */
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
        DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
        return DenseMatrixPolicy::Function(
            [this, variable_indices, forward_index, reverse_index](
                auto&& state_parameters, auto&& state_variables, auto&& forcing_terms)
            {
              auto forward_rate = forcing_terms.GetRowVariable();
              auto reverse_rate = forcing_terms.GetRowVariable();

              // For each phase instance, calculate the reaction rate and update the forcing terms for reactants and products
              for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
              {
                // Calculate the forward and reverse rates
                state_parameters.ForEachRow(
                    [&](const double& forward_rate_constant,
                        const double& reverse_rate_constant,
                        const double& solvent,
                        double& forward_rate,
                        double& reverse_rate)
                    {
                      // We want to end up with rates in units of mol m-3 s-1 (same as state variable units per second)
                      forward_rate = forward_rate_constant / std::pow(solvent, reactants_.size() - 1);
                      reverse_rate = reverse_rate_constant / std::pow(solvent, products_.size() - 1);
                    },
                    state_parameters.GetConstColumnView(forward_index),
                    state_parameters.GetConstColumnView(reverse_index),
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
          std::map<std::string, std::vector<AerosolPropertyProvider<DenseMatrixPolicy>>> /* providers */) const
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
        DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
        return SparseMatrixPolicy::Function(
            [this, variable_indices, jacobian_indices, forward_index, reverse_index](
                auto&& state_parameters, auto&& state_variables, auto&& jacobian_values)
            {
              auto d_forward_rate_d_ind = jacobian_values.GetBlockVariable();
              auto d_reverse_rate_d_ind = jacobian_values.GetBlockVariable();
              auto jac_id = jacobian_indices.indices_.AsVector().begin();

              // For each phase instance, calculate the partial derivatives for the Jacobian entries
              for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
              {
                // Calculate the partial derivatives of the forward and reverse rates with respect to the independent
                // variable Calculate partials for independent reactants
                for (std::size_t i_ind = 0; i_ind < reactants_.size(); ++i_ind)
                {
                  // Start the rate calculation with the rate constant and solvent
                  jacobian_values.ForEachBlock(
                      [&](const double& forward_rate_constant, const double& solvent, double& partial)
                      { partial = forward_rate_constant / std::pow(solvent, reactants_.size() - 1); },
                      state_parameters.GetConstColumnView(forward_index),
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
                  // apply partial to dependent reactants
                  for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
                  {
                    jacobian_values.ForEachBlock(
                        [&](const double& partial, double& jacobian) { jacobian -= partial; },
                        d_forward_rate_d_ind,
                        jacobian_values.GetBlockView(*jac_id++));
                  }
                  // apply partial to dependent products
                  for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
                  {
                    jacobian_values.ForEachBlock(
                        [&](const double& partial, double& jacobian) { jacobian += partial; },
                        d_forward_rate_d_ind,
                        jacobian_values.GetBlockView(*jac_id++));
                  }
                }
                // Calculate partials for independent products
                for (std::size_t i_ind = 0; i_ind < products_.size(); ++i_ind)
                {
                  // Start the rate calculation with the rate constant and solvent
                  jacobian_values.ForEachBlock(
                      [&](const double& reverse_rate_constant, const double& solvent, double& partial)
                      { partial = reverse_rate_constant / std::pow(solvent, products_.size() - 1); },
                      state_parameters.GetConstColumnView(reverse_index),
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
                  // apply partial to dependent reactants
                  for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
                  {
                    jacobian_values.ForEachBlock(
                        [&](const double& partial, double& jacobian) { jacobian += partial; },
                        d_reverse_rate_d_ind,
                        jacobian_values.GetBlockView(*jac_id++));
                  }
                  // apply partial to dependent products
                  for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
                  {
                    jacobian_values.ForEachBlock(
                        [&](const double& partial, double& jacobian) { jacobian -= partial; },
                        d_reverse_rate_d_ind,
                        jacobian_values.GetBlockView(*jac_id++));
                  }
                }
                // Calculate partials for independent solvent
                jacobian_values.ForEachBlock(
                    [&](const double& forward_rate_constant,
                        const double& reverse_rate_constant,
                        const double& solvent,
                        double& forward_partial,
                        double& reverse_partial)
                    {
                      forward_partial = forward_rate_constant * (1 - static_cast<int>(reactants_.size())) /
                                        std::pow(solvent, reactants_.size());
                      reverse_partial = reverse_rate_constant * (1 - static_cast<int>(products_.size())) /
                                        std::pow(solvent, products_.size());
                    },
                    state_parameters.GetConstColumnView(forward_index),
                    state_parameters.GetConstColumnView(reverse_index),
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
                // apply partials to dependent reactants
                for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
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
                // apply partials to dependent products
                for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
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
      std::pair<std::size_t, std::size_t> GetParameterIndices(
          const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
      ) const
      {
        std::string forward_param = phase_.name_ + "." + uuid_ + ".k_forward";
        std::string reverse_param = phase_.name_ + "." + uuid_ + ".k_reverse";
        if (state_parameter_indices.find(forward_param) == state_parameter_indices.end())
        {
          throw std::runtime_error(
              "Internal Error: GetParameterIndices: Forward rate constant parameter " + forward_param +
              " not found in state_parameter_indices");
        }
        if (state_parameter_indices.find(reverse_param) == state_parameter_indices.end())
        {
          throw std::runtime_error(
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
          throw std::runtime_error(
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
              throw std::runtime_error(
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
              throw std::runtime_error(
                  "Internal Error: GetStateVariableIndices: Product variable " + product_var +
                  " not found in state_variable_indices");
            }
            indices.product_indices_[i_phase][i_product] = state_variable_indices.at(product_var);
          }
          std::string solvent_var = prefix + "." + phase_.name_ + "." + solvent_.name_;
          if (state_variable_indices.find(solvent_var) == state_variable_indices.end())
          {
            throw std::runtime_error(
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
  }  // namespace process
}  // namespace miam