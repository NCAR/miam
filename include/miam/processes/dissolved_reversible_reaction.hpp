// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/rate_expression.hpp>
#include <miam/representations/aerosol_property.hpp>
#include <miam/representations/aerosol_property_descriptor.hpp>
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
    RateConstantExpression forward_rate_constant_;  ///< Forward rate constant expression
    RateConstantExpression reverse_rate_constant_;  ///< Reverse rate constant expression
    std::vector<micm::Species> reactants_;          ///< Reactant species
    std::vector<micm::Species> products_;           ///< Product species
    micm::Species solvent_;                         ///< Solvent species
    micm::Phase phase_;                             ///< Phase in which the reaction occurs
    std::string uuid_;                              ///< Unique identifier for the reaction
    double solvent_floor_{
      1.0e-20
    };  ///< Floor [mol m⁻³] added to [S] in ([S]+δ)^n denominator to prevent singularity as [S] → 0

    DissolvedReversibleReaction() = delete;

    /// @brief Constructor
    DissolvedReversibleReaction(
        RateConstantExpression forward_rate_constant,
        RateConstantExpression reverse_rate_constant,
        const std::vector<micm::Species>& reactants,
        const std::vector<micm::Species>& products,
        micm::Species solvent,
        micm::Phase phase,
        double solvent_floor = 1.0e-20)
        : forward_rate_constant_(std::move(forward_rate_constant)),
          reverse_rate_constant_(std::move(reverse_rate_constant)),
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
          forward_rate_constant_, reverse_rate_constant_, reactants_, products_, solvent_, phase_, solvent_floor_);
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
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const auto& /* providers */)
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

      // Hoist both variant visits to the host so the concrete expression PODs are captured by
      // value into MICM_LAMBDA (built inside `MakeReversibleUpdateFn`, a function template, so the
      // extended lambda is not nested in the generic std::visit lambda, which NVCC forbids).
      return std::visit(
          [forward_index, reverse_index, num_params](const auto& fwd_expr, const auto& rev_expr)
          { return MakeReversibleUpdateFn<DenseMatrixPolicy>(fwd_expr, rev_expr, forward_index, reverse_index, num_params); },
          forward_rate_constant_,
          reverse_rate_constant_);
    }

    /// @brief Builds the forward/reverse rate-constant update callable for concrete expression alternatives.
    /// @details Public because NVCC forbids private functions that define an extended lambda.
    template<typename DenseMatrixPolicy, typename FwdT, typename RevT>
    static std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    MakeReversibleUpdateFn(
        const FwdT& fwd_expr, const RevT& rev_expr, std::size_t forward_index, std::size_t reverse_index, std::size_t num_params)
    {
      const FwdT fwd_copy = fwd_expr;
      const RevT rev_copy = rev_expr;
      DenseMatrixPolicy state_parameters{ 1, num_params, 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> conditions_vector;
      return DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
              const typename DenseMatrixPolicy::ViewType& params_view)
          {
            params_view.ForEachRowStrict(
                [fwd_copy](const micm::Conditions& condition, micm::Real& parameter)
                { parameter = fwd_copy.Calculate(condition); },
                conditions_view,
                params_view.GetColumnView(forward_index));
            params_view.ForEachRowStrict(
                [rev_copy](const micm::Conditions& condition, micm::Real& parameter)
                { parameter = rev_copy.Calculate(condition); },
                conditions_view,
                params_view.GetColumnView(reverse_index));
          },
          conditions_vector,
          state_parameters);
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

  };
}  // namespace miam
