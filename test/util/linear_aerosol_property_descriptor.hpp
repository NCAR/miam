// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace miam::test
{
  /// @brief Test-only descriptor: `value = base_value + Σ_k coeffs[k] * variables[dep_indices[k]]`.
  /// @details Exists to exercise the Jacobian-through-aerosol-properties code paths in
  ///          `HenrysLawPhaseTransfer` with an analytic linear dependency the tests can verify.
  template<typename DenseMatrixPolicy>
  class LinearAerosolPropertyDescriptor
  {
   public:
    LinearAerosolPropertyDescriptor() = default;

    LinearAerosolPropertyDescriptor(
        double base_value,
        std::vector<std::size_t> dependent_variable_indices,
        std::vector<double> coefficients)
        : base_value_(base_value),
          dependent_variable_indices_(std::move(dependent_variable_indices)),
          coefficients_(std::move(coefficients))
    {
    }

    const std::vector<std::size_t>& DependentVariableIndices() const
    {
      return dependent_variable_indices_;
    }

    void Evaluate(
        const DenseMatrixPolicy& /*state_parameters*/,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result) const
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
      {
        double val = base_value_;
        for (std::size_t k = 0; k < dependent_variable_indices_.size(); ++k)
          val += coefficients_[k] * state_variables[row][dependent_variable_indices_[k]];
        result[row][0] = val;
      }
    }

    void EvaluateAndDerivatives(
        const DenseMatrixPolicy& /*state_parameters*/,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result,
        DenseMatrixPolicy& partials) const
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
      {
        double val = base_value_;
        for (std::size_t k = 0; k < dependent_variable_indices_.size(); ++k)
        {
          val += coefficients_[k] * state_variables[row][dependent_variable_indices_[k]];
          partials[row][k] = coefficients_[k];
        }
        result[row][0] = val;
      }
    }

   private:
    double base_value_ = 0.0;
    std::vector<std::size_t> dependent_variable_indices_{};
    std::vector<double> coefficients_{};
  };
}  // namespace miam::test
