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
#include <map>
#include <set>
#include <string>

namespace miam
{
  /// @brief A dissolved reversible reaction builder
  /// @details Builder class for constructing DissolvedReversibleReaction objects.
  ///
  ///          Forward and reverse rate constants are configured per representation prefix (via
  ///          AddForwardRateConstant / AddReverseRateConstant), mirroring DissolvedReactionBuilder,
  ///          because the kinetics may differ between aerosol representations. The equilibrium
  ///          constant, by contrast, is an intrinsic thermodynamic property and is set once for the
  ///          whole reaction (via SetEquilibriumConstant).
  ///
  ///          For each representation prefix, exactly two of {forward rate constant, reverse rate
  ///          constant, equilibrium constant} must be determinable; the third is derived. Because
  ///          the equilibrium constant is shared, supplying it plus either a forward or reverse rate
  ///          constant per prefix is the common case.
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

    /// @brief Adds a forward rate constant for a specific representation prefix
    DissolvedReversibleReactionBuilder& AddForwardRateConstant(const std::string& prefix, const auto& forward_rate_constant)
    {
      forward_rate_constants_[prefix] = [forward_rate_constant](const micm::Conditions& conditions)
      { return forward_rate_constant.Calculate(conditions); };
      return *this;
    }

    /// @brief Adds a reverse rate constant for a specific representation prefix
    DissolvedReversibleReactionBuilder& AddReverseRateConstant(const std::string& prefix, const auto& reverse_rate_constant)
    {
      reverse_rate_constants_[prefix] = [reverse_rate_constant](const micm::Conditions& conditions)
      { return reverse_rate_constant.Calculate(conditions); };
      return *this;
    }

    /// @brief Adds a reverse rate constant from Arrhenius parameters for a specific representation prefix
    DissolvedReversibleReactionBuilder& AddReverseRateConstant(
        const std::string& prefix,
        const micm::ArrheniusRateConstantParameters& params)
    {
      reverse_rate_constants_[prefix] = [params](const micm::Conditions& conditions)
      { return micm::CalculateArrhenius(params, conditions.temperature_, conditions.pressure_); };
      return *this;
    }

    /// @brief Sets the (shared, intrinsic) equilibrium constant function
    DissolvedReversibleReactionBuilder& SetEquilibriumConstant(const auto& equilibrium_constant)
    {
      equilibrium_constant_ = [equilibrium_constant](const micm::Conditions& conditions)
      { return equilibrium_constant.Calculate(conditions); };
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

      // Collect the set of representation prefixes that have at least one rate constant configured
      std::set<std::string> prefixes;
      for (const auto& [prefix, fn] : forward_rate_constants_)
        prefixes.insert(prefix);
      for (const auto& [prefix, fn] : reverse_rate_constants_)
        prefixes.insert(prefix);
      if (prefixes.empty())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "DissolvedReversibleReactionBuilder requires at least one forward or reverse rate constant via "
            "AddForwardRateConstant/AddReverseRateConstant.");
      }

      // For each prefix, require exactly two of {forward, reverse, equilibrium} and derive the third.
      // The equilibrium constant is shared across all prefixes.
      std::map<std::string, std::function<double(const micm::Conditions&)>> forward = forward_rate_constants_;
      std::map<std::string, std::function<double(const micm::Conditions&)>> reverse = reverse_rate_constants_;
      const bool has_eq = static_cast<bool>(equilibrium_constant_);
      for (const auto& prefix : prefixes)
      {
        const bool has_fwd = forward.count(prefix) > 0;
        const bool has_rev = reverse.count(prefix) > 0;
        const int num_set = (has_fwd ? 1 : 0) + (has_rev ? 1 : 0) + (has_eq ? 1 : 0);
        if (num_set != 2)
        {
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
              "DissolvedReversibleReactionBuilder: for representation '" + prefix +
                  "', exactly two of forward rate constant, reverse rate constant, or equilibrium constant must be set.");
        }
        if (has_eq)
        {
          if (!has_fwd)
          {
            auto eq_const = equilibrium_constant_;
            auto rev_const = reverse.at(prefix);
            forward[prefix] = [eq_const, rev_const](const micm::Conditions& conditions)
            { return eq_const(conditions) * rev_const(conditions); };
          }
          else if (!has_rev)
          {
            auto eq_const = equilibrium_constant_;
            auto fwd_const = forward.at(prefix);
            reverse[prefix] = [eq_const, fwd_const](const micm::Conditions& conditions)
            { return fwd_const(conditions) / eq_const(conditions); };
          }
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
    std::map<std::string, std::function<double(const micm::Conditions& conditions)>>
        forward_rate_constants_;  ///< Per-prefix forward rate constant functions
    std::map<std::string, std::function<double(const micm::Conditions& conditions)>>
        reverse_rate_constants_;  ///< Per-prefix reverse rate constant functions
    std::function<double(const micm::Conditions& conditions)>
        equilibrium_constant_;         ///< Shared equilibrium constant function (intrinsic, representation-independent)
    double solvent_floor_{ 1.0e-20 };  ///< Floor δ [mol m⁻³] added to [S] in ([S]+δ)^n denominator; see SetSolventFloor()
  };
}  // namespace miam
