// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/rate_expression.hpp>
#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <concepts>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <variant>

namespace miam
{
  /// @brief A dissolved reversible reaction builder
  /// @details Builder class for constructing DissolvedReversibleReaction objects.
  ///
  ///          Exactly two of {forward rate constant, reverse rate constant,
  ///          equilibrium constant} must be provided; the third is derived at
  ///          `Build()` time by wrapping the two into a `CombinedExpression`.
  class DissolvedReversibleReactionBuilder
  {
   public:
    DissolvedReversibleReactionBuilder() = default;

    /// @brief Sets the phase in which the reaction occurs
    DissolvedReversibleReactionBuilder& SetPhase(const micm::Phase& phase)
    {
      phase_ = phase;
      phase_is_set_ = true;
      return *this;
    }

    /// @brief Sets the reactant species
    DissolvedReversibleReactionBuilder& SetReactants(const std::vector<micm::Species>& reactants)
    {
      reactants_ = reactants;
      return *this;
    }

    /// @brief Sets the product species
    DissolvedReversibleReactionBuilder& SetProducts(const std::vector<micm::Species>& products)
    {
      products_ = products;
      return *this;
    }

    /// @brief Sets the solvent species
    DissolvedReversibleReactionBuilder& SetSolvent(const micm::Species& solvent)
    {
      solvent_ = solvent;
      solvent_is_set_ = true;
      return *this;
    }

    /// @brief Sets the floor \f$\delta\f$ [mol m⁻³] added to the solvent in the denominator
    ///        to prevent singularity as \f$[S] \to 0\f$. Default: 1e-20.
    DissolvedReversibleReactionBuilder& SetSolventFloor(double solvent_floor)
    {
      solvent_floor_ = solvent_floor;
      return *this;
    }

    /// @brief Sets the forward rate constant from any operand-expression alternative
    template<class Expression>
      requires std::constructible_from<detail::CombinedExpressionOperand, Expression>
    DissolvedReversibleReactionBuilder& SetForwardRateConstant(Expression expression)
    {
      forward_rate_constant_ = detail::CombinedExpressionOperand{ std::move(expression) };
      return *this;
    }

    /// @brief Sets the forward rate constant from MICM Arrhenius parameters
    DissolvedReversibleReactionBuilder& SetForwardRateConstant(const micm::ArrheniusRateConstantParameters& params)
    {
      forward_rate_constant_ = detail::CombinedExpressionOperand{ ArrheniusExpression{ params } };
      return *this;
    }

    /// @brief Sets the reverse rate constant from any operand-expression alternative
    template<class Expression>
      requires std::constructible_from<detail::CombinedExpressionOperand, Expression>
    DissolvedReversibleReactionBuilder& SetReverseRateConstant(Expression expression)
    {
      reverse_rate_constant_ = detail::CombinedExpressionOperand{ std::move(expression) };
      return *this;
    }

    /// @brief Sets the reverse rate constant from MICM Arrhenius parameters
    DissolvedReversibleReactionBuilder& SetReverseRateConstant(const micm::ArrheniusRateConstantParameters& params)
    {
      reverse_rate_constant_ = detail::CombinedExpressionOperand{ ArrheniusExpression{ params } };
      return *this;
    }

    /// @brief Sets the equilibrium constant from any operand-expression alternative
    template<class Expression>
      requires std::constructible_from<detail::CombinedExpressionOperand, Expression>
    DissolvedReversibleReactionBuilder& SetEquilibriumConstant(Expression expression)
    {
      equilibrium_constant_ = detail::CombinedExpressionOperand{ std::move(expression) };
      return *this;
    }

    /// @brief Builds and returns the DissolvedReversibleReaction object
    DissolvedReversibleReaction Build() const
    {
      if (reactants_.empty())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires at least one reactant species.");
      }
      if (products_.empty())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires at least one product species.");
      }
      if (!phase_is_set_)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires the phase to be set.");
      }
      if (!solvent_is_set_)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires the solvent to be set.");
      }

      const bool has_fwd = forward_rate_constant_.has_value();
      const bool has_rev = reverse_rate_constant_.has_value();
      const bool has_eq = equilibrium_constant_.has_value();
      const int num_set = (has_fwd ? 1 : 0) + (has_rev ? 1 : 0) + (has_eq ? 1 : 0);
      if (num_set != 2)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder: exactly two of forward rate constant, reverse rate constant, or "
            "equilibrium constant must be set.");
      }

      RateConstantExpression forward = ResolveForward();
      RateConstantExpression reverse = ResolveReverse();
      return DissolvedReversibleReaction(
          std::move(forward), std::move(reverse), reactants_, products_, solvent_, phase_, solvent_floor_);
    }

   private:
    /// @brief Lift an operand variant into the outer `RateConstantExpression` variant.
    static RateConstantExpression OperandToRateConstant(const detail::CombinedExpressionOperand& operand)
    {
      return std::visit([](const auto& expr) -> RateConstantExpression { return expr; }, operand);
    }

    /// @brief Return the forward rate constant, deriving it from `K_eq * k_r` if unset.
    RateConstantExpression ResolveForward() const
    {
      if (forward_rate_constant_.has_value())
        return OperandToRateConstant(*forward_rate_constant_);
      return CombinedExpression{ *equilibrium_constant_, *reverse_rate_constant_, CombinedExpressionOp::Multiply };
    }

    /// @brief Return the reverse rate constant, deriving it from `k_f / K_eq` if unset.
    RateConstantExpression ResolveReverse() const
    {
      if (reverse_rate_constant_.has_value())
        return OperandToRateConstant(*reverse_rate_constant_);
      return CombinedExpression{ *forward_rate_constant_, *equilibrium_constant_, CombinedExpressionOp::Divide };
    }

    micm::Phase phase_;
    bool phase_is_set_ = false;
    std::vector<micm::Species> reactants_;
    std::vector<micm::Species> products_;
    micm::Species solvent_;
    bool solvent_is_set_ = false;
    std::optional<detail::CombinedExpressionOperand> forward_rate_constant_;
    std::optional<detail::CombinedExpressionOperand> reverse_rate_constant_;
    std::optional<detail::CombinedExpressionOperand> equilibrium_constant_;
    double solvent_floor_{ 1.0e-20 };
  };
}  // namespace miam
