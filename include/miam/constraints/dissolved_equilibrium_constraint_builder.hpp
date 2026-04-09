// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>

#include <micm/system/conditions.hpp>

#include <functional>
#include <stdexcept>

namespace miam
{
  namespace constraint
  {
    /// @brief Builder for DissolvedEquilibriumConstraint
    class DissolvedEquilibriumConstraintBuilder
    {
     public:
      DissolvedEquilibriumConstraintBuilder() = default;

      DissolvedEquilibriumConstraintBuilder& SetPhase(const micm::Phase& phase)
      {
        phase_ = phase;
        phase_is_set_ = true;
        return *this;
      }

      DissolvedEquilibriumConstraintBuilder& SetReactants(const std::vector<micm::Species>& reactants)
      {
        reactants_ = reactants;
        return *this;
      }

      DissolvedEquilibriumConstraintBuilder& SetProducts(const std::vector<micm::Species>& products)
      {
        products_ = products;
        return *this;
      }

      DissolvedEquilibriumConstraintBuilder& SetAlgebraicSpecies(const micm::Species& species)
      {
        algebraic_species_ = species;
        algebraic_species_is_set_ = true;
        return *this;
      }

      DissolvedEquilibriumConstraintBuilder& SetSolvent(const micm::Species& solvent)
      {
        solvent_ = solvent;
        solvent_is_set_ = true;
        return *this;
      }

      template<typename T>
          requires requires(const T& t, const micm::Conditions& c) { { t.Calculate(c) }; }
      DissolvedEquilibriumConstraintBuilder& SetEquilibriumConstant(const T& equilibrium_constant)
      {
        equilibrium_constant_ = [equilibrium_constant](const micm::Conditions& conditions)
        { return equilibrium_constant.Calculate(conditions); };
        return *this;
      }

      DissolvedEquilibriumConstraintBuilder& SetEquilibriumConstant(
          std::function<double(const micm::Conditions&)> equilibrium_constant)
      {
        equilibrium_constant_ = std::move(equilibrium_constant);
        return *this;
      }

      DissolvedEquilibriumConstraint Build() const
      {
        if (!phase_is_set_)
          throw std::runtime_error("DissolvedEquilibriumConstraintBuilder requires the phase to be set.");
        if (reactants_.empty())
          throw std::runtime_error("DissolvedEquilibriumConstraintBuilder requires at least one reactant.");
        if (products_.empty())
          throw std::runtime_error("DissolvedEquilibriumConstraintBuilder requires at least one product.");
        if (!algebraic_species_is_set_)
          throw std::runtime_error("DissolvedEquilibriumConstraintBuilder requires the algebraic species to be set.");
        if (!solvent_is_set_)
          throw std::runtime_error("DissolvedEquilibriumConstraintBuilder requires the solvent to be set.");
        if (!equilibrium_constant_)
          throw std::runtime_error("DissolvedEquilibriumConstraintBuilder requires the equilibrium constant to be set.");

        return DissolvedEquilibriumConstraint(
            equilibrium_constant_, reactants_, products_, algebraic_species_, solvent_, phase_);
      }

     private:
      micm::Phase phase_;
      bool phase_is_set_ = false;
      std::vector<micm::Species> reactants_;
      std::vector<micm::Species> products_;
      micm::Species algebraic_species_;
      bool algebraic_species_is_set_ = false;
      micm::Species solvent_;
      bool solvent_is_set_ = false;
      std::function<double(const micm::Conditions& conditions)> equilibrium_constant_;
    };
  }  // namespace constraint
}  // namespace miam
