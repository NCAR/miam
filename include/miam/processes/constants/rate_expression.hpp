// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <cmath>
#include <cstdint>
#include <variant>

namespace miam
{
  /// @brief Constant-valued expression: returns `value_` regardless of conditions.
  /// @details Used for temperature- and pressure-independent rate constants,
  ///          equilibrium constants, and Henry's law constants.
  class UserDefinedConstantExpression
  {
   public:
    micm::Real value_ = 0.0;

    UserDefinedConstantExpression() = default;
    explicit UserDefinedConstantExpression(micm::Real value)
        : value_(value)
    {
    }

    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& /*conditions*/) const
    {
      return value_;
    }
  };

  /// @brief Arrhenius-form expression: k = A * T^B * exp(C/T + D + E*P)
  /// @details Wraps a MICM Arrhenius parameter set. Used for temperature-
  ///          and pressure-dependent rate constants.
  class ArrheniusExpression
  {
   public:
    micm::ArrheniusRateConstantParameters params_;

    ArrheniusExpression() = default;
    explicit ArrheniusExpression(const micm::ArrheniusRateConstantParameters& params)
        : params_(params)
    {
    }

    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const
    {
      return micm::CalculateArrhenius(params_, conditions.temperature_, conditions.pressure_);
    }
  };

  /// @brief van 't Hoff-form expression: f(T) = A * exp(C * (1/T0 - 1/T))
  /// @details Used for MIAM equilibrium constants (K_eq) and for rate constants
  ///          whose temperature dependence is expressed in this form directly.
  class VantHoffExpression
  {
   public:
    micm::Real A_ = 1.0;    ///< Value at the reference temperature `T0_`
    micm::Real C_ = 0.0;    ///< Temperature-dependence parameter [K] (e.g. Ea/R)
    micm::Real T0_ = 298.15;  ///< Reference temperature [K]

    VantHoffExpression() = default;
    VantHoffExpression(micm::Real A, micm::Real C, micm::Real T0 = 298.15)
        : A_(A),
          C_(C),
          T0_(T0)
    {
    }

    /// @brief Implicit conversion from the legacy `EquilibriumConstant` class.
    VantHoffExpression(const EquilibriumConstant& eq)
        : A_(eq.parameters_.A_),
          C_(eq.parameters_.C_),
          T0_(eq.parameters_.T0_)
    {
    }

    /// @brief Implicit conversion from the legacy `EquilibriumConstantParameters` struct.
    VantHoffExpression(const EquilibriumConstantParameters& params)
        : A_(params.A_),
          C_(params.C_),
          T0_(params.T0_)
    {
    }

    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const
    {
      return A_ * std::exp(C_ * (1.0 / T0_ - 1.0 / conditions.temperature_));
    }
  };

  /// @brief Henry's Law-form expression: HLC(T) = HLC_ref * exp(C * (1/T - 1/T0))
  /// @details Same functional form as `VantHoffExpression` with the sign of C
  ///          flipped, matching the standard Henry's law temperature convention.
  ///          Used for Henry's law constants (mol m⁻³ Pa⁻¹).
  class HenrysLawExpression
  {
   public:
    micm::Real HLC_ref_ = 1.0;  ///< HLC at the reference temperature `T0_` [mol m⁻³ Pa⁻¹]
    micm::Real C_ = 0.0;        ///< Temperature-dependence parameter [K]
    micm::Real T0_ = 298.15;    ///< Reference temperature [K]

    HenrysLawExpression() = default;
    HenrysLawExpression(micm::Real HLC_ref, micm::Real C, micm::Real T0 = 298.15)
        : HLC_ref_(HLC_ref),
          C_(C),
          T0_(T0)
    {
    }

    /// @brief Implicit conversion from the legacy `HenrysLawConstant` class.
    HenrysLawExpression(const HenrysLawConstant& hlc)
        : HLC_ref_(hlc.parameters_.HLC_ref_),
          C_(hlc.parameters_.C_),
          T0_(hlc.parameters_.T0_)
    {
    }

    /// @brief Implicit conversion from the legacy `HenrysLawConstantParameters` struct.
    HenrysLawExpression(const HenrysLawConstantParameters& params)
        : HLC_ref_(params.HLC_ref_),
          C_(params.C_),
          T0_(params.T0_)
    {
    }

    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const
    {
      return HLC_ref_ * std::exp(C_ * (1.0 / conditions.temperature_ - 1.0 / T0_));
    }
  };

  /// @brief Binary operator selector for `CombinedExpression`.
  enum class CombinedExpressionOp : std::uint8_t
  {
    Multiply,
    Divide
  };

  namespace detail
  {
    /// @brief Non-recursive operand set for `CombinedExpression`.
    /// @details Excludes `CombinedExpression` itself to prevent recursive variant
    ///          types and to bound the compile-time expression tree at one level.
    using CombinedExpressionOperand =
        std::variant<ArrheniusExpression, VantHoffExpression, UserDefinedConstantExpression>;

    /// @brief Device-callable form of `CombinedExpressionOperand`.
    /// @details `std::visit` cannot run in device code, so the operand stores every alternative
    ///          and a tag, and `Calculate` selects the alternative with a `switch`.
    class DeviceExpressionOperand
    {
     public:
      enum class Kind : std::uint8_t
      {
        Arrhenius,
        VantHoff,
        UserDefinedConstant
      };

      Kind kind_ = Kind::UserDefinedConstant;
      ArrheniusExpression arrhenius_{};
      VantHoffExpression vant_hoff_{};
      UserDefinedConstantExpression user_defined_constant_{};

      DeviceExpressionOperand() = default;

      /// @brief Host-side conversion from the variant operand.
      DeviceExpressionOperand(const CombinedExpressionOperand& operand)
      {
        if (const auto* arrhenius = std::get_if<ArrheniusExpression>(&operand))
        {
          kind_ = Kind::Arrhenius;
          arrhenius_ = *arrhenius;
        }
        else if (const auto* vant_hoff = std::get_if<VantHoffExpression>(&operand))
        {
          kind_ = Kind::VantHoff;
          vant_hoff_ = *vant_hoff;
        }
        else
        {
          kind_ = Kind::UserDefinedConstant;
          user_defined_constant_ = std::get<UserDefinedConstantExpression>(operand);
        }
      }

      MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const
      {
        switch (kind_)
        {
          case Kind::Arrhenius: return arrhenius_.Calculate(conditions);
          case Kind::VantHoff: return vant_hoff_.Calculate(conditions);
          case Kind::UserDefinedConstant: break;
        }
        return user_defined_constant_.Calculate(conditions);
      }
    };
  }  // namespace detail

  /// @brief Product or quotient of two non-composite expressions.
  /// @details Emitted by `DissolvedReversibleReactionBuilder::Build()` when the
  ///          missing forward or reverse rate constant must be derived from the
  ///          other rate constant and the equilibrium constant.
  ///
  ///          For \f$k_f = K_{eq} \cdot k_r\f$: `CombinedExpression(K_eq, k_r, Multiply)`.
  ///          For \f$k_r = k_f / K_{eq}\f$: `CombinedExpression(k_f, K_eq, Divide)`.
  class CombinedExpression
  {
   public:
    detail::DeviceExpressionOperand left_;
    detail::DeviceExpressionOperand right_;
    CombinedExpressionOp op_ = CombinedExpressionOp::Multiply;

    CombinedExpression() = default;
    CombinedExpression(
        detail::CombinedExpressionOperand left,
        detail::CombinedExpressionOperand right,
        CombinedExpressionOp op)
        : left_(left),
          right_(right),
          op_(op)
    {
    }

    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const
    {
      const micm::Real l = left_.Calculate(conditions);
      const micm::Real r = right_.Calculate(conditions);
      return op_ == CombinedExpressionOp::Multiply ? l * r : l / r;
    }
  };

  /// @brief Compile-time list of rate-constant expressions.
  /// @details Configuration-side type stored on `DissolvedReaction`,
  ///          `DissolvedReversibleReaction`, and any other class whose
  ///          per-reaction rate constant may take on multiple forms.
  ///          Host-only — visited exactly once at `Model::FinalizeProcessSetup`
  ///          time to populate `ConstantsBucket`; never crosses to device.
  using RateConstantExpression = std::
      variant<ArrheniusExpression, VantHoffExpression, UserDefinedConstantExpression, CombinedExpression>;

  /// @brief Compile-time list of equilibrium-constant expressions.
  using EquilibriumConstantExpression =
      std::variant<VantHoffExpression, UserDefinedConstantExpression>;

  /// @brief Compile-time list of Henry's law-constant expressions.
  using HenrysLawConstantExpression =
      std::variant<HenrysLawExpression, UserDefinedConstantExpression>;

  /// @brief Evaluate an expression variant on the host.
  /// @details Convenience wrapper for callsites that still need a `double` from
  ///          a variant (e.g. the OLD `std::function`-returning factories on
  ///          `Model` and the per-process/per-constraint classes). Removed
  ///          entirely in Layer C/D once the buckets take over on-device
  ///          evaluation.
  template<class ExpressionVariant>
  inline micm::Real EvaluateExpression(const ExpressionVariant& expression, const micm::Conditions& conditions)
  {
    return std::visit([&](const auto& expr) { return expr.Calculate(conditions); }, expression);
  }
}  // namespace miam
