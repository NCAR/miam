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
    HenryLawPhaseTransferBuilder& SetDiffusionCoefficient(double diffusion_coefficient)
    {
      if (diffusion_coefficient <= 0)
        throw std::invalid_argument("Diffusion coefficient must be positive.");
      diffusion_coefficient_ = diffusion_coefficient;
      diffusion_coefficient_is_set_ = true;
      return *this;
    }

    /// @brief Sets the mass accommodation coefficient [dimensionless, 0–1]
    HenryLawPhaseTransferBuilder& SetAccommodationCoefficient(double accommodation_coefficient)
    {
      if (accommodation_coefficient <= 0)
        throw std::invalid_argument("Accommodation coefficient must be positive.");
      accommodation_coefficient_ = accommodation_coefficient;
      accommodation_coefficient_is_set_ = true;
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
      if (!diffusion_coefficient_is_set_)
        throw std::runtime_error("HenryLawPhaseTransferBuilder requires the diffusion coefficient to be set.");
      if (!accommodation_coefficient_is_set_)
        throw std::runtime_error("HenryLawPhaseTransferBuilder requires the accommodation coefficient to be set.");

      double gas_molecular_weight = gas_species_.GetProperty<double>("molecular weight [kg mol-1]");
      double solvent_molecular_weight = solvent_.GetProperty<double>("molecular weight [kg mol-1]");
      double solvent_density = solvent_.GetProperty<double>("density [kg m-3]");

      return HenryLawPhaseTransfer(
          henry_law_constant_,
          gas_species_,
          condensed_species_,
          solvent_,
          condensed_phase_,
          diffusion_coefficient_,
          accommodation_coefficient_,
          gas_molecular_weight,
          solvent_molecular_weight,
          solvent_density);
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
    double diffusion_coefficient_ = 0.0;  ///< Gas-phase diffusion coefficient [m² s⁻¹]
    bool diffusion_coefficient_is_set_ = false;
    double accommodation_coefficient_ = 0.0;  ///< Mass accommodation coefficient [dimensionless, 0–1]
    bool accommodation_coefficient_is_set_ = false;
  };
}  // namespace miam
