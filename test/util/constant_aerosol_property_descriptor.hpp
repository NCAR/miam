// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <vector>

namespace miam::test
{
  /// @brief Test-only descriptor that returns a fixed value and zero partial derivatives.
  /// @details Used to isolate the consumer-side (`HenrysLawPhaseTransfer`) Jacobian machinery
  ///          from the specifics of any real representation's aerosol-property formulas.
  template<typename DenseMatrixPolicy>
  class ConstantAerosolPropertyDescriptor
  {
   public:
    ConstantAerosolPropertyDescriptor() = default;

    ConstantAerosolPropertyDescriptor(double value, std::vector<std::size_t> dependent_variable_indices)
        : value_(value),
          dependent_variable_indices_(std::move(dependent_variable_indices))
    {
    }

    const std::vector<std::size_t>& DependentVariableIndices() const
    {
      return dependent_variable_indices_;
    }

    void Evaluate(
        const DenseMatrixPolicy& /*state_parameters*/,
        const DenseMatrixPolicy& /*state_variables*/,
        DenseMatrixPolicy& result) const
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
        result[row][0] = value_;
    }

    void EvaluateAndDerivatives(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result,
        DenseMatrixPolicy& /*partials*/) const
    {
      Evaluate(state_parameters, state_variables, result);
    }

   private:
    double value_ = 0.0;
    std::vector<std::size_t> dependent_variable_indices_{};
  };
}  // namespace miam::test
