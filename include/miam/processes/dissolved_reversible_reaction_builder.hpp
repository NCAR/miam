// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/rate_constant.hpp>
#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <map>
#include <optional>
#include <set>
#include <string>

namespace miam
{
  /// @brief A dissolved reversible reaction builder
  /// @details Builder class for constructing DissolvedReversibleReaction objects.
  ///
  ///          Exactly two of {forward rate constant, reverse rate
  ///          constant, equilibrium constant} must be determinable; the third is derived.
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

    /// @brief Sets the forward rate constant
    DissolvedReversibleReactionBuilder& SetForwardRateConstant(const RateConstant& forward_rate_constant)
    {
      forward_rate_constant_ = forward_rate_constant;
      return *this;
    }

    /// @brief Sets the reverse rate constant
    DissolvedReversibleReactionBuilder& SetReverseRateConstant(const RateConstant& reverse_rate_constant)
    {
      reverse_rate_constant_ = reverse_rate_constant;
      return *this;
    }

    /// @brief Sets the equilibrium constant
    DissolvedReversibleReactionBuilder& SetEquilibriumConstant(const EquilibriumConstant& equilibrium_constant)
    {
      equilibrium_constant_ = equilibrium_constant;
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

      int num_set = 0;
      if (forward_rate_constant_)
        ++num_set;
      if (reverse_rate_constant_)
        ++num_set;
      if (equilibrium_constant_)
        ++num_set;
      if (num_set != 2)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder: exactly two of forward rate constant, reverse rate constant, or "
            "equilibrium constant must be set.");
      }

      return DissolvedReversibleReaction(
          forward_rate_constant_, reverse_rate_constant_, reactants_, products_, solvent_, phase_, solvent_floor_, equilibrium_constant_);
    }

   private:
    micm::Phase phase_;                     ///< Phase in which the reaction occurs
    bool phase_is_set_ = false;             ///< Flag to track if the phase has been set
    std::vector<micm::Species> reactants_;  ///< Reactant species
    std::vector<micm::Species> products_;   ///< Product species
    micm::Species solvent_;                 ///< Solvent species
    bool solvent_is_set_ = false;           ///< Flag to track if the solvent has been set
    std::optional<RateConstant> forward_rate_constant_;         ///< Forward rate constant
    std::optional<RateConstant> reverse_rate_constant_;         ///< Reverse rate constant
    std::optional<EquilibriumConstant> equilibrium_constant_;  ///< Equilibrium constant
    double solvent_floor_{ 1.0e-20 };  ///< Floor δ [mol m⁻³] added to [S] in ([S]+δ)^n denominator; see SetSolventFloor()
  };
}  // namespace miam
