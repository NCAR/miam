// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <functional>
#include <string>

namespace miam
{
  /// @brief A dissolved reversible reaction builder
  /// @details Builder class for constructing DissolvedReversibleReaction objects.
  ///          Exactly two of SetForwardRateConstant, SetReverseRateConstant, and
  ///          SetEquilibriumConstant must be called; the third is derived automatically.
  ///          Rate constants are shared across all representations of the phase.
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

    /// @brief Sets the forward rate constant function (shared across all representations)
    DissolvedReversibleReactionBuilder& SetForwardRateConstant(
        std::function<double(const micm::Conditions&)> forward_rate_constant)
    {
      forward_rate_constant_ = std::move(forward_rate_constant);
      return *this;
    }

    /// @brief Sets the forward rate constant from an object with a Calculate() method
    template<typename RateConstantT>
      requires requires(const RateConstantT& kc, const micm::Conditions& c) {
        { kc.Calculate(c) } -> std::convertible_to<double>;
      }
    DissolvedReversibleReactionBuilder& SetForwardRateConstant(const RateConstantT& forward_rate_constant)
    {
      forward_rate_constant_ = [forward_rate_constant](const micm::Conditions& c)
      { return forward_rate_constant.Calculate(c); };
      return *this;
    }

    /// @brief Sets the forward rate constant from Arrhenius parameters
    DissolvedReversibleReactionBuilder& SetForwardRateConstant(
        const micm::ArrheniusRateConstantParameters& forward_rate_constant)
    {
      forward_rate_constant_ = [forward_rate_constant](const micm::Conditions& c)
      { return micm::CalculateArrhenius(forward_rate_constant, c.temperature_, c.pressure_); };
      return *this;
    }

    /// @brief Sets the reverse rate constant function (shared across all representations)
    DissolvedReversibleReactionBuilder& SetReverseRateConstant(
        std::function<double(const micm::Conditions&)> reverse_rate_constant)
    {
      reverse_rate_constant_ = std::move(reverse_rate_constant);
      return *this;
    }

    /// @brief Sets the reverse rate constant from an object with a Calculate() method
    template<typename RateConstantT>
      requires requires(const RateConstantT& kc, const micm::Conditions& c) {
        { kc.Calculate(c) } -> std::convertible_to<double>;
      }
    DissolvedReversibleReactionBuilder& SetReverseRateConstant(const RateConstantT& reverse_rate_constant)
    {
      reverse_rate_constant_ = [reverse_rate_constant](const micm::Conditions& c)
      { return reverse_rate_constant.Calculate(c); };
      return *this;
    }

    /// @brief Sets the reverse rate constant from Arrhenius parameters
    DissolvedReversibleReactionBuilder& SetReverseRateConstant(
        const micm::ArrheniusRateConstantParameters& reverse_rate_constant)
    {
      reverse_rate_constant_ = [reverse_rate_constant](const micm::Conditions& c)
      { return micm::CalculateArrhenius(reverse_rate_constant, c.temperature_, c.pressure_); };
      return *this;
    }

    /// @brief Sets the equilibrium constant function (shared across all representations)
    /// @details When provided, the missing rate constant is derived: k_f = K_eq * k_r or k_r = k_f / K_eq
    DissolvedReversibleReactionBuilder& SetEquilibriumConstant(
        std::function<double(const micm::Conditions&)> equilibrium_constant)
    {
      equilibrium_constant_ = std::move(equilibrium_constant);
      return *this;
    }

    /// @brief Sets the equilibrium constant from an object with a Calculate() method
    template<typename RateConstantT>
      requires requires(const RateConstantT& kc, const micm::Conditions& c) {
        { kc.Calculate(c) } -> std::convertible_to<double>;
      }
    DissolvedReversibleReactionBuilder& SetEquilibriumConstant(const RateConstantT& equilibrium_constant)
    {
      equilibrium_constant_ = [equilibrium_constant](const micm::Conditions& c)
      { return equilibrium_constant.Calculate(c); };
      return *this;
    }

    /// @brief Sets the equilibrium constant from Arrhenius parameters
    DissolvedReversibleReactionBuilder& SetEquilibriumConstant(
        const micm::ArrheniusRateConstantParameters& equilibrium_constant)
    {
      equilibrium_constant_ = [equilibrium_constant](const micm::Conditions& c)
      { return micm::CalculateArrhenius(equilibrium_constant, c.temperature_, c.pressure_); };
      return *this;
    }

    /// @brief Builds and returns the DissolvedReversibleReaction object
    DissolvedReversibleReaction Build() const
    {
      int num_set = 0;
      if (forward_rate_constant_) ++num_set;
      if (reverse_rate_constant_) ++num_set;
      if (equilibrium_constant_) ++num_set;
      if (num_set != 2)
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires exactly two of SetForwardRateConstant, "
            "SetReverseRateConstant, or SetEquilibriumConstant to be called.");
      if (reactants_.empty())
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires at least one reactant species.");
      if (products_.empty())
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires at least one product species.");
      if (!phase_is_set_)
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires the phase to be set.");
      if (!solvent_is_set_)
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires the solvent to be set.");

      auto forward = forward_rate_constant_;
      auto reverse = reverse_rate_constant_;
      if (equilibrium_constant_)
      {
        auto keq = equilibrium_constant_;
        if (!forward)
        {
          auto rev = reverse;
          forward = [keq, rev](const micm::Conditions& c) { return keq(c) * rev(c); };
        }
        else
        {
          auto fwd = forward;
          reverse = [keq, fwd](const micm::Conditions& c) { return fwd(c) / keq(c); };
        }
      }
      return DissolvedReversibleReaction(forward, reverse, reactants_, products_, solvent_, phase_, solvent_floor_);
    }

   private:
    micm::Phase phase_;                     ///< Phase in which the reaction occurs
    bool phase_is_set_ = false;             ///< Flag to track if the phase has been set
    std::vector<micm::Species> reactants_;  ///< Reactant species
    std::vector<micm::Species> products_;   ///< Product species
    micm::Species solvent_;                 ///< Solvent species
    bool solvent_is_set_ = false;           ///< Flag to track if the solvent has been set
    double solvent_floor_{ 1.0e-20 };       ///< Floor δ [mol m⁻³]
    std::function<double(const micm::Conditions&)> forward_rate_constant_;  ///< Shared forward rate constant
    std::function<double(const micm::Conditions&)> reverse_rate_constant_;  ///< Shared reverse rate constant
    std::function<double(const micm::Conditions&)> equilibrium_constant_;   ///< Shared equilibrium constant
  };
}  // namespace miam
