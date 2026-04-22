// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/henry_law_equilibrium_constraint.hpp>

#include <micm/system/conditions.hpp>

#include <functional>
#include <stdexcept>

namespace miam
{
  namespace constraint
  {
    /// @brief Builder for HenryLawEquilibriumConstraint
    class HenryLawEquilibriumConstraintBuilder
    {
     public:
      HenryLawEquilibriumConstraintBuilder() = default;

      HenryLawEquilibriumConstraintBuilder& SetGasSpecies(const micm::Species& species)
      {
        gas_species_ = species;
        gas_species_is_set_ = true;
        return *this;
      }

      HenryLawEquilibriumConstraintBuilder& SetCondensedSpecies(const micm::Species& species)
      {
        condensed_species_ = species;
        condensed_species_is_set_ = true;
        return *this;
      }

      HenryLawEquilibriumConstraintBuilder& SetSolvent(const micm::Species& solvent)
      {
        solvent_ = solvent;
        solvent_is_set_ = true;
        return *this;
      }

      HenryLawEquilibriumConstraintBuilder& SetCondensedPhase(const micm::Phase& phase)
      {
        condensed_phase_ = phase;
        condensed_phase_is_set_ = true;
        return *this;
      }

      template<typename T>
          requires requires(const T& t, const micm::Conditions& c) { { t.Calculate(c) }; }
      HenryLawEquilibriumConstraintBuilder& SetHenryLawConstant(const T& henry_law_constant)
      {
        henry_law_constant_ = [henry_law_constant](const micm::Conditions& conditions)
        { return henry_law_constant.Calculate(conditions); };
        return *this;
      }

      HenryLawEquilibriumConstraintBuilder& SetHenryLawConstant(
          std::function<double(const micm::Conditions&)> henry_law_constant)
      {
        henry_law_constant_ = std::move(henry_law_constant);
        return *this;
      }

      HenryLawEquilibriumConstraintBuilder& SetMwSolvent(double Mw_solvent)
      {
        Mw_solvent_ = Mw_solvent;
        Mw_solvent_is_set_ = true;
        return *this;
      }

      HenryLawEquilibriumConstraintBuilder& SetRhoSolvent(double rho_solvent)
      {
        rho_solvent_ = rho_solvent;
        rho_solvent_is_set_ = true;
        return *this;
      }

      HenryLawEquilibriumConstraint Build() const
      {
        if (!gas_species_is_set_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the gas species to be set.");
        if (!condensed_species_is_set_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the condensed species to be set.");
        if (!solvent_is_set_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the solvent to be set.");
        if (!condensed_phase_is_set_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the condensed phase to be set.");
        if (!henry_law_constant_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the Henry's Law constant to be set.");
        if (!Mw_solvent_is_set_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the solvent molecular weight to be set.");
        if (!rho_solvent_is_set_)
          throw std::runtime_error("HenryLawEquilibriumConstraintBuilder requires the solvent density to be set.");

        return HenryLawEquilibriumConstraint(
            henry_law_constant_, gas_species_, condensed_species_, solvent_, condensed_phase_, Mw_solvent_, rho_solvent_);
      }

     private:
      micm::Species gas_species_;
      bool gas_species_is_set_ = false;
      micm::Species condensed_species_;
      bool condensed_species_is_set_ = false;
      micm::Species solvent_;
      bool solvent_is_set_ = false;
      micm::Phase condensed_phase_;
      bool condensed_phase_is_set_ = false;
      std::function<double(const micm::Conditions& conditions)> henry_law_constant_;
      double Mw_solvent_ = 0.0;
      bool Mw_solvent_is_set_ = false;
      double rho_solvent_ = 0.0;
      bool rho_solvent_is_set_ = false;
    };
  }  // namespace constraint
}  // namespace miam
