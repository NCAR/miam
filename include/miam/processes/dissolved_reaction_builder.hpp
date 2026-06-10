// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reaction.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/system/conditions.hpp>

#include <functional>
#include <map>
#include <string>

namespace miam
{
  /// @brief A dissolved reaction builder
  /// @details Builder class for constructing DissolvedReaction objects.
  class DissolvedReactionBuilder
  {
   public:
    DissolvedReactionBuilder() = default;

    /// @brief Sets the phase in which the reaction occurs
    DissolvedReactionBuilder& SetPhase(const micm::Phase& phase)
    {
      phase_ = phase;
      phase_is_set_ = true;
      return *this;
    }

    /// @brief Sets the reactant species
    DissolvedReactionBuilder& SetReactants(const std::vector<micm::Species>& reactants)
    {
      reactants_ = reactants;
      return *this;
    }

    /// @brief Sets the product species
    DissolvedReactionBuilder& SetProducts(const std::vector<micm::Species>& products)
    {
      products_ = products;
      return *this;
    }

    /// @brief Sets the solvent species
    DissolvedReactionBuilder& SetSolvent(const micm::Species& solvent)
    {
      solvent_ = solvent;
      solvent_is_set_ = true;
      return *this;
    }

    /// @brief Sets the floor \f$\delta\f$ [mol m⁻³] added to the solvent in the denominator
    ///        to prevent singularity as \f$[S] \to 0\f$. Default: 1e-20.
    DissolvedReactionBuilder& SetSolventFloor(double solvent_floor)
    {
      solvent_floor_ = solvent_floor;
      return *this;
    }

    /// @brief Sets the minimum half-life for rate capping [s]
    /// @details When positive, the reaction rate is smoothly capped so that no
    ///          reactant is consumed faster than the specified half-life. When zero
    ///          (default), rate capping is disabled with zero runtime overhead.
    DissolvedReactionBuilder& SetMinHalflife(double min_halflife)
    {
      min_halflife_ = min_halflife;
      return *this;
    }

    /// @brief Adds a rate constant for a specific representation prefix
    DissolvedReactionBuilder& AddRateConstant(
        const std::string& prefix,
        std::function<double(const micm::Conditions&)> rate_constant)
    {
      rate_constants_[prefix] = std::move(rate_constant);
      return *this;
    }

    /// @brief Builds and returns the DissolvedReaction object
    DissolvedReaction Build() const
    {
      if (rate_constants_.empty())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReactionBuilder requires at least one rate constant via AddRateConstant.");
      }
      if (reactants_.empty())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReactionBuilder requires at least one reactant species.");
      }
      if (products_.empty())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReactionBuilder requires at least one product species.");
      }
      if (!phase_is_set_)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReactionBuilder requires the phase to be set.");
      }
      if (!solvent_is_set_)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReactionBuilder requires the solvent to be set.");
      }
      return DissolvedReaction(rate_constants_, reactants_, products_, solvent_, phase_, solvent_floor_, min_halflife_);
    }

   private:
    micm::Phase phase_;                                                        ///< Phase in which the reaction occurs
    bool phase_is_set_ = false;                                                ///< Flag to track if the phase has been set
    std::vector<micm::Species> reactants_;                                     ///< Reactant species
    std::vector<micm::Species> products_;                                      ///< Product species
    micm::Species solvent_;                                                    ///< Solvent species
    bool solvent_is_set_ = false;                                              ///< Flag to track if the solvent has been set
    std::map<std::string, std::function<double(const micm::Conditions& conditions)>> rate_constants_;  ///< Per-prefix rate constants
    double solvent_floor_{ 1.0e-20 };  ///< Floor δ [mol m⁻³] added to [S] in ([S]+δ)^n denominator; see SetSolventFloor()
    double min_halflife_{ 0.0 };       ///< Minimum half-life for rate capping [s]
  };
}  // namespace miam
