// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/linear_constraint.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <stdexcept>

namespace miam
{
  /// @brief Builder for LinearConstraint
  class LinearConstraintBuilder
  {
   public:
    LinearConstraintBuilder() = default;

    LinearConstraintBuilder& SetAlgebraicSpecies(const micm::Phase& phase, const micm::Species& species)
    {
      algebraic_phase_ = phase;
      algebraic_species_ = species;
      algebraic_is_set_ = true;
      return *this;
    }

    LinearConstraintBuilder& AddTerm(const micm::Phase& phase, const micm::Species& species, double coefficient)
    {
      terms_.push_back({ phase, species, coefficient, "" });
      return *this;
    }

    /// @brief Adds a term restricted to a single aerosol representation.
    /// @details Only meaningful for instanced (condensed) phases. The term contributes only the named
    ///          representation's variable instead of summing over all representations holding the phase.
    LinearConstraintBuilder& AddTerm(
        const micm::Phase& phase,
        const micm::Species& species,
        double coefficient,
        const std::string& representation)
    {
      terms_.push_back({ phase, species, coefficient, representation });
      return *this;
    }

    LinearConstraintBuilder& SetConstant(double constant)
    {
      if (diagnose_from_state_)
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MUTUALLY_EXCLUSIVE_PARAMETERS,
            "LinearConstraintBuilder: SetConstant() and DiagnoseConstantFromState() are mutually exclusive.");
      constant_ = constant;
      constant_is_set_ = true;
      return *this;
    }

    /// @brief Diagnose the constant from the current state at the beginning of each solve step.
    ///        The constant C will be computed as C = sum(c_i * [species_i]) from the state variables.
    ///        This is useful for mass conservation constraints where the total mass varies by grid cell
    ///        and may change between solve steps due to emissions, transport, or deposition.
    LinearConstraintBuilder& DiagnoseConstantFromState()
    {
      if (constant_is_set_)
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MUTUALLY_EXCLUSIVE_PARAMETERS,
            "LinearConstraintBuilder: SetConstant() and DiagnoseConstantFromState() are mutually exclusive.");
      diagnose_from_state_ = true;
      return *this;
    }

    LinearConstraint Build() const
    {
      if (!algebraic_is_set_)
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "LinearConstraintBuilder requires the algebraic species to be set.");
      if (terms_.empty())
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_REQUIRED_PARAMETER,
            "LinearConstraintBuilder requires at least one term.");

      return LinearConstraint(algebraic_phase_, algebraic_species_, terms_, constant_, diagnose_from_state_);
    }

   private:
    micm::Phase algebraic_phase_;
    micm::Species algebraic_species_;
    bool algebraic_is_set_ = false;
    std::vector<LinearConstraint::Term> terms_;
    double constant_{ 0.0 };
    bool constant_is_set_ = false;
    bool diagnose_from_state_ = false;
  };
}  // namespace miam
