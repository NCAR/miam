// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/linear_constraint.hpp>

#include <stdexcept>

namespace miam
{
  namespace constraint
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
        terms_.push_back({ phase, species, coefficient });
        return *this;
      }

      LinearConstraintBuilder& SetConstant(double constant)
      {
        constant_ = constant;
        return *this;
      }

      LinearConstraint Build() const
      {
        if (!algebraic_is_set_)
          throw std::runtime_error("LinearConstraintBuilder requires the algebraic species to be set.");
        if (terms_.empty())
          throw std::runtime_error("LinearConstraintBuilder requires at least one term.");

        return LinearConstraint(algebraic_phase_, algebraic_species_, terms_, constant_);
      }

     private:
      micm::Phase algebraic_phase_;
      micm::Species algebraic_species_;
      bool algebraic_is_set_ = false;
      std::vector<LinearConstraint::Term> terms_;
      double constant_{ 0.0 };
    };
  }  // namespace constraint
}  // namespace miam
