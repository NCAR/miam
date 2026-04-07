// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace miam
{
  /// @brief Aerosol/cloud particle properties that representations can provide
  enum class AerosolProperty
  {
    EffectiveRadius,      // [m]
    NumberConcentration,  // [# m^-3]
    PhaseVolumeFraction   // [dimensionless, 0-1]
  };

  /// @brief A provider for a single aerosol property, created at setup time by a representation instance
  /// @details Captures all needed parameter/variable column indices internally. Operates on
  ///          ForEachRow-compatible column views — no per-cell indexing. Partial derivatives are
  ///          written into columns of a pre-allocated DenseMatrixPolicy, one column per dependent variable.
  /// @tparam DenseMatrixPolicy The dense matrix type used for state data
  template<typename DenseMatrixPolicy>
  struct AerosolPropertyProvider
  {
    /// @brief State variable indices that this property has non-zero partial derivatives with respect to
    /// @details Fixed at creation time. Used by processes to:
    ///   1. Determine Jacobian sparsity (NonZeroJacobianElements)
    ///   2. Know how many columns the partials matrix needs
    ///   3. Map partials columns back to state variable indices
    std::vector<std::size_t> dependent_variable_indices;

    /// @brief Compute the property value for all grid cells in the current group
    /// @details Called inside a ForEachRow loop — receives column views and writes into a RowVariable.
    ///   The provider internally calls params_view.GetConstColumnView(...) and
    ///   vars_view.GetConstColumnView(...) using its captured indices.
    ///
    ///   Parameters:
    ///     params_view: const GroupView of state parameters
    ///     vars_view:   const GroupView of state variables
    ///     result:      mutable RowVariable to write the property value into
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ComputeValue;

    /// @brief Compute the property value AND partial derivatives for all grid cells in the current group
    /// @details Called inside a ForEachRow loop.
    ///
    ///   Parameters:
    ///     params_view:      const GroupView of state parameters
    ///     vars_view:        const GroupView of state variables
    ///     result:           mutable RowVariable for the property value
    ///     partials_matrix:  mutable GroupView of the partials DenseMatrix
    ///                       (num_cells x num_dependent_variables)
    ///                       Column k corresponds to d(property)/d(var[dependent_variable_indices[k]])
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&, DenseMatrixPolicy&)>
        ComputeValueAndDerivatives;
  };
}  // namespace miam
