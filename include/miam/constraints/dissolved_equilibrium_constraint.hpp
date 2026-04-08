// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace miam
{
  namespace constraint
  {
    /// @brief A dissolved equilibrium constraint
    /// @details Replaces the ODE row for a designated algebraic species with the
    ///          steady-state equilibrium condition:
    ///
    ///          \f$ G = K_{eq} \frac{\prod[R_i]}{[S]^{n_r - 1}}
    ///                       - \frac{\prod[P_j]}{[S]^{n_p - 1}} = 0 \f$
    ///
    ///          where \f$ K_{eq} = k_f / k_r \f$, \f$ R_i \f$ are reactants, \f$ P_j \f$ are
    ///          products, and \f$ S \f$ is the solvent. One of the product species is designated
    ///          as the algebraic variable whose ODE row is replaced by this constraint.
    class DissolvedEquilibriumConstraint
    {
     public:
      std::function<double(const micm::Conditions& conditions)> equilibrium_constant_;  ///< K_eq function
      std::vector<micm::Species> reactants_;                                             ///< Reactant species
      std::vector<micm::Species> products_;                                              ///< Product species
      micm::Species algebraic_species_;  ///< Product species whose ODE row is replaced
      micm::Species solvent_;            ///< Solvent species
      micm::Phase phase_;                ///< Phase in which the reaction occurs
      std::string uuid_;                 ///< Unique identifier

      /// @brief Shared mutable K_eq values, bridging UpdateStateParametersFunction → constraint functions.
      ///        One value per grid cell.
      std::shared_ptr<std::vector<double>> k_eq_values_;

      DissolvedEquilibriumConstraint() = delete;

      /// @brief Constructor
      DissolvedEquilibriumConstraint(
          std::function<double(const micm::Conditions& conditions)> equilibrium_constant,
          const std::vector<micm::Species>& reactants,
          const std::vector<micm::Species>& products,
          const micm::Species& algebraic_species,
          micm::Species solvent,
          micm::Phase phase)
          : equilibrium_constant_(equilibrium_constant),
            reactants_(reactants),
            products_(products),
            algebraic_species_(algebraic_species),
            solvent_(solvent),
            phase_(phase),
            uuid_(miam::util::generate_uuid_v4()),
            k_eq_values_(std::make_shared<std::vector<double>>())
      {
        // Validate that the algebraic species is one of the products
        bool found = false;
        for (const auto& product : products_)
        {
          if (product.name_ == algebraic_species_.name_)
          {
            found = true;
            break;
          }
        }
        if (!found)
        {
          throw std::invalid_argument(
              "DissolvedEquilibriumConstraint: algebraic species '" + algebraic_species_.name_ +
              "' must be one of the products.");
        }
      }

      /// @brief Create a copy with a new UUID
      DissolvedEquilibriumConstraint CopyWithNewUuid() const
      {
        return DissolvedEquilibriumConstraint(
            equilibrium_constant_, reactants_, products_, algebraic_species_, solvent_, phase_);
      }

      /// @brief Returns the names of algebraic variables (one per phase instance)
      std::set<std::string> ConstraintAlgebraicVariableNames(
          const std::map<std::string, std::set<std::string>>& phase_prefixes) const
      {
        std::set<std::string> names;
        auto phase_it = phase_prefixes.find(phase_.name_);
        if (phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : phase_it->second)
          {
            names.insert(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);
          }
        }
        return names;
      }

      /// @brief Returns all species the constraint depends on
      std::set<std::string> ConstraintSpeciesDependencies(
          const std::map<std::string, std::set<std::string>>& phase_prefixes) const
      {
        std::set<std::string> species_names;
        auto phase_it = phase_prefixes.find(phase_.name_);
        if (phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : phase_it->second)
          {
            for (const auto& reactant : reactants_)
              species_names.insert(prefix + "." + phase_.name_ + "." + reactant.name_);
            for (const auto& product : products_)
              species_names.insert(prefix + "." + phase_.name_ + "." + product.name_);
            species_names.insert(prefix + "." + phase_.name_ + "." + solvent_.name_);
          }
        }
        return species_names;
      }

      /// @brief Returns non-zero constraint Jacobian element positions
      /// @details For each phase instance, the algebraic row depends on all reactants, all products, and the solvent.
      std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        std::set<std::pair<std::size_t, std::size_t>> elements;
        auto phase_it = phase_prefixes.find(phase_.name_);
        if (phase_it == phase_prefixes.end())
          return elements;

        for (const auto& prefix : phase_it->second)
        {
          std::size_t alg_row = state_variable_indices.at(
              prefix + "." + phase_.name_ + "." + algebraic_species_.name_);

          for (const auto& reactant : reactants_)
          {
            std::size_t col = state_variable_indices.at(prefix + "." + phase_.name_ + "." + reactant.name_);
            elements.insert({ alg_row, col });
          }
          for (const auto& product : products_)
          {
            std::size_t col = state_variable_indices.at(prefix + "." + phase_.name_ + "." + product.name_);
            elements.insert({ alg_row, col });
          }
          std::size_t solvent_col = state_variable_indices.at(prefix + "." + phase_.name_ + "." + solvent_.name_);
          elements.insert({ alg_row, solvent_col });
        }
        return elements;
      }

      /// @brief Returns a function that updates K_eq for each grid cell (called during UpdateStateParameters)
      template<typename DenseMatrixPolicy>
      std::function<void(const std::vector<micm::Conditions>&)> UpdateConstraintParametersFunction() const
      {
        auto k_eq_values = k_eq_values_;
        auto eq_const_fn = equilibrium_constant_;
        return [k_eq_values, eq_const_fn](const std::vector<micm::Conditions>& conditions)
        {
          k_eq_values->resize(conditions.size());
          for (std::size_t i = 0; i < conditions.size(); ++i)
          {
            (*k_eq_values)[i] = eq_const_fn(conditions[i]);
          }
        };
      }

      /// @brief Returns a function that computes constraint residuals G(y) = 0
      /// @details G = K_eq * prod([R_i]) / [S]^(n_r-1) - prod([P_j]) / [S]^(n_p-1)
      template<typename DenseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
        auto k_eq_values = k_eq_values_;
        std::size_t n_reactants = reactants_.size();
        std::size_t n_products = products_.size();

        return [indices, k_eq_values, n_reactants, n_products](
                   const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& residual)
        {
          for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
          {
            std::size_t alg_idx = indices.algebraic_indices_[i_phase];
            for (std::size_t row = 0; row < state_variables.NumRows(); ++row)
            {
              double k_eq = (*k_eq_values)[row];
              double solvent = state_variables[row][indices.solvent_indices_[i_phase]];

              // Forward part: K_eq * prod([R_i]) / [S]^(n_r - 1)
              double forward = k_eq / std::pow(solvent, static_cast<double>(n_reactants) - 1.0);
              for (std::size_t r = 0; r < n_reactants; ++r)
                forward *= state_variables[row][indices.reactant_indices_[i_phase][r]];

              // Reverse part: prod([P_j]) / [S]^(n_p - 1)
              double reverse = 1.0 / std::pow(solvent, static_cast<double>(n_products) - 1.0);
              for (std::size_t p = 0; p < n_products; ++p)
                reverse *= state_variables[row][indices.product_indices_[i_phase][p]];

              residual[row][alg_idx] = forward - reverse;
            }
          }
        };
      }

      /// @brief Returns a function that computes constraint Jacobian entries (subtracts dG/dy)
      /// @details Follows MICM convention: jac[row][col] -= dG/dy
      template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices,
          const SparseMatrixPolicy& jacobian) const
      {
        auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
        auto jac_indices = GetJacobianIndices(indices, jacobian);
        auto k_eq_values = k_eq_values_;
        std::size_t n_reactants = reactants_.size();
        std::size_t n_products = products_.size();

        return [indices, jac_indices, k_eq_values, n_reactants, n_products](
                   const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian_values)
        {
          for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
          {
            auto jac_it = jac_indices.indices_[i_phase].begin();

            for (std::size_t row = 0; row < state_variables.NumRows(); ++row)
            {
              double k_eq = (*k_eq_values)[row];
              double solvent = state_variables[row][indices.solvent_indices_[i_phase]];

              // dG/d[R_i]: partial of forward part w.r.t. each reactant
              for (std::size_t i_ind = 0; i_ind < n_reactants; ++i_ind)
              {
                double dG_dRi = k_eq / std::pow(solvent, static_cast<double>(n_reactants) - 1.0);
                for (std::size_t r = 0; r < n_reactants; ++r)
                {
                  if (r != i_ind)
                    dG_dRi *= state_variables[row][indices.reactant_indices_[i_phase][r]];
                }
                jacobian_values.AsVector()[*(jac_it + i_ind)] -= dG_dRi;
              }
              jac_it += n_reactants;

              // dG/d[P_j]: partial of reverse part w.r.t. each product
              for (std::size_t i_ind = 0; i_ind < n_products; ++i_ind)
              {
                double dG_dPj = -1.0 / std::pow(solvent, static_cast<double>(n_products) - 1.0);
                for (std::size_t p = 0; p < n_products; ++p)
                {
                  if (p != i_ind)
                    dG_dPj *= state_variables[row][indices.product_indices_[i_phase][p]];
                }
                jacobian_values.AsVector()[*(jac_it + i_ind)] -= dG_dPj;
              }
              jac_it += n_products;

              // dG/d[S]: partial w.r.t. solvent
              double forward_term = k_eq;
              for (std::size_t r = 0; r < n_reactants; ++r)
                forward_term *= state_variables[row][indices.reactant_indices_[i_phase][r]];
              double dG_dS_forward = forward_term * (1.0 - static_cast<double>(n_reactants)) /
                                     std::pow(solvent, static_cast<double>(n_reactants));

              double reverse_term = 1.0;
              for (std::size_t p = 0; p < n_products; ++p)
                reverse_term *= state_variables[row][indices.product_indices_[i_phase][p]];
              double dG_dS_reverse = reverse_term * (1.0 - static_cast<double>(n_products)) /
                                     std::pow(solvent, static_cast<double>(n_products));

              jacobian_values.AsVector()[*jac_it] -= (dG_dS_forward - dG_dS_reverse);
              ++jac_it;
            }
          }
        };
      }

     private:
      /// @brief Helper struct for state variable indices across phase instances
      struct StateVariableIndices
      {
        std::size_t number_of_phase_instances_;
        micm::Matrix<std::size_t> reactant_indices_;   // [n_instances x n_reactants]
        micm::Matrix<std::size_t> product_indices_;    // [n_instances x n_products]
        std::vector<std::size_t> solvent_indices_;      // [n_instances]
        std::vector<std::size_t> algebraic_indices_;    // [n_instances]
      };

      /// @brief Helper struct for Jacobian sparse matrix indices
      struct JacobianIndices
      {
        // [n_instances][n_reactants + n_products + 1 (solvent)] per grid row
        std::vector<std::vector<std::size_t>> indices_;
      };

      /// @brief Build state variable indices for all phase instances
      StateVariableIndices GetStateVariableIndices(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        StateVariableIndices indices;
        auto phase_it = phase_prefixes.find(phase_.name_);
        if (phase_it == phase_prefixes.end())
        {
          throw std::runtime_error(
              "Internal Error: DissolvedEquilibriumConstraint: Phase " + phase_.name_ + " not found in phase_prefixes");
        }
        const auto& prefixes = phase_it->second;
        indices.number_of_phase_instances_ = prefixes.size();
        indices.reactant_indices_ = micm::Matrix<std::size_t>(prefixes.size(), reactants_.size());
        indices.product_indices_ = micm::Matrix<std::size_t>(prefixes.size(), products_.size());
        indices.solvent_indices_.resize(prefixes.size());
        indices.algebraic_indices_.resize(prefixes.size());

        std::size_t i_phase = 0;
        for (const auto& prefix : prefixes)
        {
          for (std::size_t r = 0; r < reactants_.size(); ++r)
          {
            std::string var = prefix + "." + phase_.name_ + "." + reactants_[r].name_;
            indices.reactant_indices_[i_phase][r] = state_variable_indices.at(var);
          }
          for (std::size_t p = 0; p < products_.size(); ++p)
          {
            std::string var = prefix + "." + phase_.name_ + "." + products_[p].name_;
            indices.product_indices_[i_phase][p] = state_variable_indices.at(var);
          }
          indices.solvent_indices_[i_phase] =
              state_variable_indices.at(prefix + "." + phase_.name_ + "." + solvent_.name_);
          indices.algebraic_indices_[i_phase] =
              state_variable_indices.at(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);
          ++i_phase;
        }
        return indices;
      }

      /// @brief Build Jacobian sparse matrix indices for all phase instances
      JacobianIndices GetJacobianIndices(
          const StateVariableIndices& var_indices,
          const auto& jacobian) const
      {
        JacobianIndices jac_indices;
        jac_indices.indices_.resize(var_indices.number_of_phase_instances_);

        for (std::size_t i_phase = 0; i_phase < var_indices.number_of_phase_instances_; ++i_phase)
        {
          std::size_t alg_row = var_indices.algebraic_indices_[i_phase];
          auto& inst = jac_indices.indices_[i_phase];

          // Reactant columns
          for (std::size_t r = 0; r < reactants_.size(); ++r)
            inst.push_back(jacobian.VectorIndex(0, alg_row, var_indices.reactant_indices_[i_phase][r]));

          // Product columns
          for (std::size_t p = 0; p < products_.size(); ++p)
            inst.push_back(jacobian.VectorIndex(0, alg_row, var_indices.product_indices_[i_phase][p]));

          // Solvent column
          inst.push_back(jacobian.VectorIndex(0, alg_row, var_indices.solvent_indices_[i_phase]));
        }
        return jac_indices;
      }
    };
  }  // namespace constraint
}  // namespace miam
