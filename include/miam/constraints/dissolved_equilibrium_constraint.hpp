// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/rate_expression.hpp>
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
#include <unordered_map>
#include <vector>

namespace miam
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
  ///
  ///          To prevent singularity as \f$[S] \to 0\f$, the solvent denominator is
  ///          regularized by a small floor \f$\delta\f$ (\c solvent_floor_):
  ///
  ///          \f$ G = K_{eq} \frac{[S]\prod[R_i]}{([S]+\delta)^{n_r}}
  ///                       - \frac{[S]\prod[P_j]}{([S]+\delta)^{n_p}} = 0 \f$
  class DissolvedEquilibriumConstraint
  {
   public:
    EquilibriumConstantExpression equilibrium_constant_;  ///< K_eq expression
    std::vector<micm::Species> reactants_;                ///< Reactant species
    std::vector<micm::Species> products_;                 ///< Product species
    micm::Species algebraic_species_;  ///< Product species whose ODE row is replaced
    micm::Species solvent_;            ///< Solvent species
    micm::Phase phase_;                ///< Phase in which the reaction occurs
    std::string uuid_;                 ///< Unique identifier
    double solvent_floor_{ 1.0e-20 };  ///< Floor \f$\delta\f$ [mol m⁻³] added to \f$[S]\f$ in \f$([S]+\delta)^n\f$
                                       ///< denominator to prevent singularity as \f$[S] \to 0\f$

    DissolvedEquilibriumConstraint() = delete;

    /// @brief Constructor
    DissolvedEquilibriumConstraint(
        EquilibriumConstantExpression equilibrium_constant,
        const std::vector<micm::Species>& reactants,
        const std::vector<micm::Species>& products,
        const micm::Species& algebraic_species,
        micm::Species solvent,
        micm::Phase phase,
        double solvent_floor = 1.0e-20)
        : equilibrium_constant_(std::move(equilibrium_constant)),
          reactants_(reactants),
          products_(products),
          algebraic_species_(algebraic_species),
          solvent_(solvent),
          phase_(phase),
          uuid_(GenerateUuid()),
          solvent_floor_(solvent_floor)
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
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_ALGEBRAIC_SPECIES_NOT_FOUND_IN_PRODUCTS,
            "DissolvedEquilibriumConstraint: algebraic species '" + algebraic_species_.name_ +
                "' must be one of the products.");
      }
    }

    /// @brief Create a copy with a new UUID
    DissolvedEquilibriumConstraint CopyWithNewUuid() const
    {
      return DissolvedEquilibriumConstraint(
          equilibrium_constant_, reactants_, products_, algebraic_species_, solvent_, phase_, solvent_floor_);
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
        std::size_t alg_row = state_variable_indices.at(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);

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

    /// @brief Returns the names of state parameters owned by this constraint (one per phase instance).
    /// @details Each phase instance writes \f$ K_{eq}(T) \f$ to a dedicated column of the
    ///          state parameter matrix every time conditions change.
    std::set<std::string> ConstraintStateParameterNames(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          names.insert(prefix + "." + phase_.name_ + "." + uuid_ + ".k_eq");
      }
      return names;
    }

    /// @brief Returns a function that writes \f$ K_{eq}(T) \f$ per grid cell into the state parameter matrix.
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateConstraintParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
    {
      std::vector<std::size_t> k_eq_indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          k_eq_indices.push_back(state_parameter_indices.at(prefix + "." + phase_.name_ + "." + uuid_ + ".k_eq"));
      }
      auto eq_const_expr = equilibrium_constant_;

      DenseMatrixPolicy state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> conditions_vector;

      return DenseMatrixPolicy::Function(
          [k_eq_indices, eq_const_expr](auto&& conditions, auto&& params)
          {
            for (const auto& k_eq_idx : k_eq_indices)
              params.ForEachRowStrict(
                  [eq_const_expr](const micm::Conditions& cond, double& k_eq)
                  { k_eq = EvaluateExpression(eq_const_expr, cond); },
                  conditions,
                  params.GetColumnView(k_eq_idx));
          },
          conditions_vector,
          state_parameters);
    }

   private:
    /// @brief Helper struct for state variable indices across phase instances
    struct StateVariableIndices
    {
      std::size_t number_of_phase_instances_;
      micm::Matrix<std::size_t> reactant_indices_;  // [n_instances x n_reactants]
      micm::Matrix<std::size_t> product_indices_;   // [n_instances x n_products]
      std::vector<std::size_t> solvent_indices_;    // [n_instances]
      std::vector<std::size_t> algebraic_indices_;  // [n_instances]
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
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedEquilibriumConstraint: Phase " + phase_.name_ + " not found in phase_prefixes");
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
        indices.solvent_indices_[i_phase] = state_variable_indices.at(prefix + "." + phase_.name_ + "." + solvent_.name_);
        indices.algebraic_indices_[i_phase] =
            state_variable_indices.at(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);
        ++i_phase;
      }
      return indices;
    }

    /// @brief Build Jacobian sparse matrix indices for all phase instances
    JacobianIndices GetJacobianIndices(const StateVariableIndices& var_indices, const auto& jacobian) const
    {
      JacobianIndices jac_indices;
      jac_indices.indices_.resize(var_indices.number_of_phase_instances_);

      std::size_t num_blocks = jacobian.NumberOfBlocks();
      for (std::size_t i_phase = 0; i_phase < var_indices.number_of_phase_instances_; ++i_phase)
      {
        std::size_t alg_row = var_indices.algebraic_indices_[i_phase];
        auto& inst = jac_indices.indices_[i_phase];

        for (std::size_t block = 0; block < num_blocks; ++block)
        {
          // Reactant columns
          for (std::size_t r = 0; r < reactants_.size(); ++r)
            inst.push_back(jacobian.VectorIndex(block, alg_row, var_indices.reactant_indices_[i_phase][r]));

          // Product columns
          for (std::size_t p = 0; p < products_.size(); ++p)
            inst.push_back(jacobian.VectorIndex(block, alg_row, var_indices.product_indices_[i_phase][p]));

          // Solvent column
          inst.push_back(jacobian.VectorIndex(block, alg_row, var_indices.solvent_indices_[i_phase]));
        }
      }
      return jac_indices;
    }
  };
}  // namespace miam
