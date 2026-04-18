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
    /// @brief A dissolved (irreversible) reaction
    /// @details Dissolved reactions involve reactants and products in solution, and
    ///          are characterized by a single rate constant.
    ///
    ///          The reaction takes the form:
    ///          \f$ \mathrm{Reactants} \rightarrow \mathrm{Products} \f$
    ///
    ///          with rate constant \f$ k \f$.
    ///
    ///          The rate expression is:
    ///          \f$ r = \frac{k}{[S]^{n_r - 1}} \prod [R_i] \f$
    ///
    ///          where [S] is the solvent concentration (mol m⁻³ air) and
    ///          \f$ n_r \f$ is the number of reactants. The solvent-normalization
    ///          factor absorbs concentration dimensions, so \f$ k \f$ always has
    ///          units of s⁻¹. To convert from a literature rate constant
    ///          \f$ k_{lit} \f$ in molar units:
    ///          \f$ k = k_{lit} \times c_{H_2O}^{n_r - 1} \f$
    ///          where \f$ c_{H_2O} = 55.556 \f$ mol/L.
    class DissolvedReaction
    {
     public:
      std::function<double(const micm::Conditions& conditions)> rate_constant_;  ///< Rate constant function
      std::vector<micm::Species> reactants_;                                     ///< Reactant species
      std::vector<micm::Species> products_;                                      ///< Product species
      micm::Species solvent_;                                                    ///< Solvent species
      micm::Phase phase_;  ///< Phase in which the reaction occurs
      std::string uuid_;   ///< Unique identifier for the reaction
      double solvent_damping_epsilon_{ 1.0e-10 };  ///< Regularization parameter to prevent singularity as solvent → 0

      DissolvedReaction() = delete;

      /// @brief Constructor
      DissolvedReaction(
          std::function<double(const micm::Conditions& conditions)> rate_constant,
          const std::vector<micm::Species>& reactants,
          const std::vector<micm::Species>& products,
          micm::Species solvent,
          micm::Phase phase,
          double solvent_damping_epsilon = 1.0e-10)
          : rate_constant_(rate_constant),
            reactants_(reactants),
            products_(products),
            solvent_(solvent),
            phase_(phase),
            uuid_(miam::util::generate_uuid_v4()),
            solvent_damping_epsilon_(solvent_damping_epsilon)
      {
      }

      /// @brief Create a copy of this reaction with a new UUID
      /// @return A new DissolvedReaction with the same properties but a unique UUID
      DissolvedReaction CopyWithNewUuid() const
      {
        return DissolvedReaction(rate_constant_, reactants_, products_, solvent_, phase_, solvent_damping_epsilon_);
      }

      /// @brief Returns a set of unique parameter names for this process
      /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
      /// species names)
      /// @return Set of unique parameter names for this process
      std::set<std::string> ProcessParameterNames(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
      {
        std::set<std::string> parameter_names;
        parameter_names.insert(phase_.name_ + "." + uuid_ + ".k");
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
      /// @details DissolvedReaction does not depend on aerosol properties.
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
          // For a forward-only reaction, the rate depends only on reactants and solvent.
          // Products do NOT appear as independent variables.

          // Reactant contributions (each reactant is an independent variable)
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
          const std::
              map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>>& /* providers */) const
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
        std::string k_param = phase_.name_ + "." + uuid_ + ".k";
        if (state_parameter_indices.find(k_param) == state_parameter_indices.end())
        {
          throw std::runtime_error(
              "Internal Error: UpdateStateParametersFunction: Rate constant parameter " + k_param +
              " not found in state_parameter_indices");
        }
        std::size_t k_index = state_parameter_indices.at(k_param);

        // Set up dummy arguments to build the function
        DenseMatrixPolicy state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        std::vector<micm::Conditions> conditions_vector;

        // return a function that updates the rate constant parameter based on the current conditions
        return DenseMatrixPolicy::Function(
            [this, k_index](auto&& conditions, auto&& params)
            {
              params.ForEachRow(
                  [&](const micm::Conditions& condition, double& parameter)
                  { parameter = rate_constant_(condition); },
                  conditions,
                  params.GetColumnView(k_index));
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
        std::size_t k_index = GetParameterIndex(state_parameter_indices);
        DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
        return DenseMatrixPolicy::Function(
            [this, variable_indices, k_index](auto&& state_parameters, auto&& state_variables, auto&& forcing_terms)
            {
              auto rate = forcing_terms.GetRowVariable();
              const double eps = solvent_damping_epsilon_;
              const std::size_t n_r = reactants_.size();

              // For each phase instance, calculate the reaction rate and update the forcing terms
              for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
              {
                // Calculate the damped rate: k * [S] / ([S] + eps)^n_r * prod([reactants])
                state_parameters.ForEachRow(
                    [&](const double& rate_constant, const double& solvent, double& rate)
                    {
                      rate = rate_constant * solvent / std::pow(solvent + eps, n_r);
                    },
                    state_parameters.GetConstColumnView(k_index),
                    state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                    rate);
                for (std::size_t r = 0; r < reactants_.size(); ++r)
                {
                  state_variables.ForEachRow(
                      [&](const double& reactant, double& rate) { rate *= reactant; },
                      state_variables.GetConstColumnView(variable_indices.reactant_indices_[i_phase][r]),
                      rate);
                }

                // Apply the reaction rate to the forcing terms for reactants and products
                for (std::size_t r = 0; r < reactants_.size(); ++r)
                {
                  state_variables.ForEachRow(
                      [&](const double& rate, double& forcing) { forcing -= rate; },
                      rate,
                      forcing_terms.GetColumnView(variable_indices.reactant_indices_[i_phase][r]));
                }
                for (std::size_t p = 0; p < products_.size(); ++p)
                {
                  state_variables.ForEachRow(
                      [&](const double& rate, double& forcing) { forcing += rate; },
                      rate,
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
        std::size_t k_index = GetParameterIndex(state_parameter_indices);
        DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
        return SparseMatrixPolicy::Function(
            [this, variable_indices, jacobian_indices, k_index](
                auto&& state_parameters, auto&& state_variables, auto&& jacobian_values)
            {
              auto d_rate_d_ind = jacobian_values.GetBlockVariable();
              auto jac_id = jacobian_indices.indices_.AsVector().begin();
              const double eps = solvent_damping_epsilon_;
              const std::size_t n_r = reactants_.size();

              // For each phase instance, calculate the partial derivatives for the Jacobian entries
              for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
              {
                // Calculate partials for independent reactants
                for (std::size_t i_ind = 0; i_ind < reactants_.size(); ++i_ind)
                {
                  // dr/d[R_i] = k * [S] / ([S]+eps)^n_r * prod(R_j, j!=i)
                  jacobian_values.ForEachBlock(
                      [&](const double& rate_constant, const double& solvent, double& partial)
                      { partial = rate_constant * solvent / std::pow(solvent + eps, n_r); },
                      state_parameters.GetConstColumnView(k_index),
                      state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                      d_rate_d_ind);
                  // add contributions to the partial from the other reactants
                  for (std::size_t r = 0; r < reactants_.size(); ++r)
                  {
                    if (r == i_ind)
                      continue;  // Skip the variable we're taking the derivative with respect to
                    jacobian_values.ForEachBlock(
                        [&](const double& reactant, double& partial) { partial *= reactant; },
                        state_variables.GetConstColumnView(variable_indices.reactant_indices_[i_phase][r]),
                        d_rate_d_ind);
                  }
                  // apply partial to dependent reactants (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
                  {
                    jacobian_values.ForEachBlock(
                        [&](const double& partial, double& jacobian) { jacobian += partial; },
                        d_rate_d_ind,
                        jacobian_values.GetBlockView(*jac_id++));
                  }
                  // apply partial to dependent products (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
                  {
                    jacobian_values.ForEachBlock(
                        [&](const double& partial, double& jacobian) { jacobian -= partial; },
                        d_rate_d_ind,
                        jacobian_values.GetBlockView(*jac_id++));
                  }
                }
                // Calculate partials for independent solvent
                // dr/d[S] = k * (eps + (1-n_r)*[S]) / ([S]+eps)^(n_r+1) * prod([R_i])
                jacobian_values.ForEachBlock(
                    [&](const double& rate_constant, const double& solvent, double& partial)
                    {
                      partial = rate_constant * (eps + (1.0 - static_cast<int>(n_r)) * solvent) /
                                std::pow(solvent + eps, n_r + 1);
                    },
                    state_parameters.GetConstColumnView(k_index),
                    state_variables.GetConstColumnView(variable_indices.solvent_indices_[i_phase]),
                    d_rate_d_ind);
                // add contributions to the partial from the reactants
                for (std::size_t r = 0; r < reactants_.size(); ++r)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& reactant, double& partial) { partial *= reactant; },
                      state_variables.GetConstColumnView(variable_indices.reactant_indices_[i_phase][r]),
                      d_rate_d_ind);
                }
                // apply partials to dependent reactants (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& partial, double& jacobian) { jacobian += partial; },
                      d_rate_d_ind,
                      jacobian_values.GetBlockView(*jac_id++));
                }
                // apply partials to dependent products (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
                {
                  jacobian_values.ForEachBlock(
                      [&](const double& partial, double& jacobian) { jacobian -= partial; },
                      d_rate_d_ind,
                      jacobian_values.GetBlockView(*jac_id++));
                }
              }
            },
            dummy_state_parameters,
            dummy_state_variables,
            jacobian);
      }

     private:
      /// @brief Helper struct for keeping track of state variable indices for reactants, products, and solvent across
      /// multiple phase instances (e.g. grid cells)
      struct StateVariableIndices
      {
        std::size_t number_of_phase_instances_;  ///< Number of instances of the phase in the system
        micm::Matrix<std::size_t>
            reactant_indices_;  ///< Matrix of state variable indices for reactants (num_prefixes x num_reactants)
        micm::Matrix<std::size_t>
            product_indices_;  ///< Matrix of state variable indices for products (num_prefixes x num_products)
        std::vector<std::size_t> solvent_indices_;  ///< Vector of state variable indices for solvent (num_prefixes)
      };

      /// @brief Helper struct for keeping track of Jacobian sparse matrix elements
      struct JacobianIndices
      {
        micm::Matrix<std::size_t>
            indices_;  // Index in sparse matrix for each dependent/independent pair (num_pairs x num_prefixes)
      };

      /// @brief Helper function to return parameter index for the rate constant
      std::size_t GetParameterIndex(
          const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
      ) const
      {
        std::string k_param = phase_.name_ + "." + uuid_ + ".k";
        if (state_parameter_indices.find(k_param) == state_parameter_indices.end())
        {
          throw std::runtime_error(
              "Internal Error: GetParameterIndex: Rate constant parameter " + k_param +
              " not found in state_parameter_indices");
        }
        return state_parameter_indices.at(k_param);
      }

      /// @brief Helper function to return variable indices for all species involved in the reaction
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
      JacobianIndices GetJacobianIndices(
          const StateVariableIndices& variable_indices,
          const auto& jacobian  // sparse matrix policy object for the Jacobian structure
      ) const
      {
        // For a forward-only reaction, independent variables are only reactants and solvent (not products).
        // Each reactant and each product depends on all reactants and the solvent.
        std::size_t num_pairs =
            (reactants_.size() + products_.size()) * (reactants_.size() + 1);  // +1 for solvent
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
