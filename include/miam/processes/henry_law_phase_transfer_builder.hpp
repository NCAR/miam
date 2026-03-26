// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/henry_law_phase_transfer.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <functional>
#include <stdexcept>
#include <string>

namespace miam
{
  namespace process
  {
    /// @brief Builder for HenryLawPhaseTransfer processes
    class HenryLawPhaseTransferBuilder
    {
     public:
      HenryLawPhaseTransferBuilder() = default;

      /// @brief Sets the condensed phase in which the solute dissolves
      HenryLawPhaseTransferBuilder& SetCondensedPhase(const micm::Phase& phase)
      {
        condensed_phase_ = phase;
        condensed_phase_is_set_ = true;
        return *this;
      }

      /// @brief Sets the gas-phase species
      HenryLawPhaseTransferBuilder& SetGasSpecies(const micm::Species& species)
      {
        gas_species_ = species;
        gas_species_is_set_ = true;
        return *this;
      }

      /// @brief Sets the condensed-phase solute species
      HenryLawPhaseTransferBuilder& SetCondensedSpecies(const micm::Species& species)
      {
        condensed_species_ = species;
        condensed_species_is_set_ = true;
        return *this;
      }

      /// @brief Sets the condensed-phase solvent species
      HenryLawPhaseTransferBuilder& SetSolvent(const micm::Species& solvent)
      {
        solvent_ = solvent;
        solvent_is_set_ = true;
        return *this;
      }

      /// @brief Sets the Henry's Law constant
      HenryLawPhaseTransferBuilder& SetHenrysLawConstant(const auto& henry_law_constant)
      {
        henry_law_constant_ = [henry_law_constant](const micm::Conditions& conditions)
        { return henry_law_constant.Calculate(conditions); };
        return *this;
      }

      /// @brief Sets the gas-phase diffusion coefficient [m² s⁻¹]
      HenryLawPhaseTransferBuilder& SetDiffusionCoefficient(double D_g)
      {
        D_g_ = D_g;
        D_g_is_set_ = true;
        return *this;
      }

      /// @brief Sets the mass accommodation coefficient [dimensionless, 0-1]
      HenryLawPhaseTransferBuilder& SetAccommodationCoefficient(double alpha)
      {
        alpha_ = alpha;
        alpha_is_set_ = true;
        return *this;
      }

      /// @brief Builds and returns the HenryLawPhaseTransfer object
      HenryLawPhaseTransfer Build() const
      {
        if (!condensed_phase_is_set_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the condensed phase to be set.");
        if (!gas_species_is_set_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the gas species to be set.");
        if (!condensed_species_is_set_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the condensed species to be set.");
        if (!solvent_is_set_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the solvent to be set.");
        if (!henry_law_constant_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the Henry's Law constant to be set.");
        if (!D_g_is_set_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the diffusion coefficient to be set.");
        if (!alpha_is_set_)
          throw std::runtime_error("HenryLawPhaseTransferBuilder requires the accommodation coefficient to be set.");

        double Mw_gas = gas_species_.GetProperty<double>("molecular weight [kg mol-1]");
        double Mw_solvent = solvent_.GetProperty<double>("molecular weight [kg mol-1]");
        double rho_solvent = solvent_.GetProperty<double>("density [kg m-3]");

        return HenryLawPhaseTransfer(
            henry_law_constant_,
            gas_species_,
            condensed_species_,
            solvent_,
            condensed_phase_,
            D_g_,
            alpha_,
            Mw_gas,
            Mw_solvent,
            rho_solvent);
      }

     private:
      micm::Phase condensed_phase_;
      bool condensed_phase_is_set_ = false;
      micm::Species gas_species_;
      bool gas_species_is_set_ = false;
      micm::Species condensed_species_;
      bool condensed_species_is_set_ = false;
      micm::Species solvent_;
      bool solvent_is_set_ = false;
      std::function<double(const micm::Conditions& conditions)> henry_law_constant_;
      double D_g_ = 0.0;
      bool D_g_is_set_ = false;
      double alpha_ = 0.0;
      bool alpha_is_set_ = false;
    };
  }  // namespace process
}  // namespace miam
