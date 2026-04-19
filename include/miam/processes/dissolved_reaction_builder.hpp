// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reaction.hpp>

#include <micm/system/conditions.hpp>

#include <functional>

namespace miam
{
  namespace process
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

      /// @brief Sets the solvent damping epsilon for regularization near zero solvent
      DissolvedReactionBuilder& SetSolventDampingEpsilon(double epsilon)
      {
        solvent_damping_epsilon_ = epsilon;
        return *this;
      }

      /// @brief Sets the maximum half-life for rate capping [s]
      /// @details When positive, the reaction rate is smoothly capped so that no
      ///          reactant is consumed faster than the specified half-life. When zero
      ///          (default), rate capping is disabled with zero runtime overhead.
      DissolvedReactionBuilder& SetMaxHalflife(double max_halflife)
      {
        max_halflife_ = max_halflife;
        return *this;
      }

      /// @brief Sets the rate constant from an object with a Calculate method
      template<typename T>
          requires requires(const T& t, const micm::Conditions& c) { { t.Calculate(c) }; }
      DissolvedReactionBuilder& SetRateConstant(const T& rate_constant)
      {
        rate_constant_ = [rate_constant](const micm::Conditions& conditions)
        { return rate_constant.Calculate(conditions); };
        return *this;
      }

      /// @brief Sets the rate constant from a std::function or lambda
      DissolvedReactionBuilder& SetRateConstant(std::function<double(const micm::Conditions&)> rate_constant)
      {
        rate_constant_ = std::move(rate_constant);
        return *this;
      }

      /// @brief Builds and returns the DissolvedReaction object
      process::DissolvedReaction Build() const
      {
        if (!rate_constant_)
        {
          throw std::runtime_error("DissolvedReactionBuilder requires the rate constant to be set.");
        }
        if (reactants_.empty())
        {
          throw std::runtime_error("DissolvedReactionBuilder requires at least one reactant species.");
        }
        if (products_.empty())
        {
          throw std::runtime_error("DissolvedReactionBuilder requires at least one product species.");
        }
        if (!phase_is_set_)
        {
          throw std::runtime_error("DissolvedReactionBuilder requires the phase to be set.");
        }
        if (!solvent_is_set_)
        {
          throw std::runtime_error("DissolvedReactionBuilder requires the solvent to be set.");
        }
        return process::DissolvedReaction(
            rate_constant_, reactants_, products_, solvent_, phase_, solvent_damping_epsilon_, max_halflife_);
      }

     private:
      micm::Phase phase_;                     ///< Phase in which the reaction occurs
      bool phase_is_set_ = false;             ///< Flag to track if the phase has been set
      std::vector<micm::Species> reactants_;  ///< Reactant species
      std::vector<micm::Species> products_;   ///< Product species
      micm::Species solvent_;                 ///< Solvent species
      bool solvent_is_set_ = false;           ///< Flag to track if the solvent has been set
      std::function<double(const micm::Conditions& conditions)> rate_constant_;  ///< Rate constant function
      double solvent_damping_epsilon_{ 1.0e-20 };  ///< Regularization parameter for solvent damping
      double max_halflife_{ 0.0 };                    ///< Maximum half-life for rate capping [s]
    };
  }  // namespace process
}  // namespace miam
