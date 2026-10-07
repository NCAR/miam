// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/util/types.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
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
  ///          ForEachRow-compatible column views - no per-cell indexing. Partial derivatives are
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
    /// @details Called inside a ForEachRow loop - receives column views and writes into a RowVariable.
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

  /// @brief Creates a provider for phi = V_phase / V_total, where V = sum([species_k] * molar_volume_k)
  /// @param species_indices State variable indices. The first phase_species_count belong to the target phase.
  /// @param molar_volumes Molar volume [m3 mol-1] of each species in species_indices
  /// @param phase_species_count Number of target phase species at the start of species_indices
  /// @param num_state_parameters Number of columns in the state parameter matrix
  /// @param num_state_variables Number of columns in the state variable matrix
  /// @details An empty species_indices gives phi = 1 and no dependent variables.
  template<typename DenseMatrixPolicy>
  AerosolPropertyProvider<DenseMatrixPolicy> MakePhaseVolumeFractionProvider(
      const std::vector<std::size_t>& species_indices,
      const std::vector<double>& molar_volumes,
      std::size_t phase_species_count,
      std::size_t num_state_parameters,
      std::size_t num_state_variables)
  {
    using IndexVector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
    using DoubleVector = typename DenseMatrixPolicy::template VectorType<double>;
    auto storage = std::make_shared<std::pair<IndexVector, DoubleVector>>(species_indices, molar_volumes);
    storage->first.CopyToDevice();
    storage->second.CopyToDevice();
    const auto species_view = storage->first.GetView();
    const auto volumes_view = storage->second.GetView();
    DenseMatrixPolicy example_params{ 1, num_state_parameters, 0.0 };
    DenseMatrixPolicy example_vars{ 1, num_state_variables, 0.0 };
    DenseMatrixPolicy example_result{ 1, 1, 0.0 };
    DenseMatrixPolicy example_partials{ 1, species_indices.size(), 0.0 };

    AerosolPropertyProvider<DenseMatrixPolicy> provider;
    provider.dependent_variable_indices = species_indices;
    auto value_function =
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& params_view,
              const typename DenseMatrixPolicy::ConstViewType& vars_view,
              const typename DenseMatrixPolicy::ViewType& result_view) {
            auto phi = result_view.GetColumnView(0);
            if (species_view.size() == 0)
            {
              params_view.ForEachRowStrict([](double& v) { v = 1.0; }, phi);
              return;
            }
            auto V_phase = result_view.GetRowVariable();
            params_view.ForEachRowStrict(
                [](double& vt, double& vp)
                {
                  vt = 0.0;
                  vp = 0.0;
                },
                phi,
                V_phase);
            const std::size_t n = species_view.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = volumes_view[k];
              if (k < phase_species_count)
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt, double& vp)
                    {
                      const double vol = c * mv;
                      vt += vol;
                      vp += vol;
                    },
                    vars_view.GetConstColumnView(species_view[k]),
                    phi,
                    V_phase);
              else
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt) { vt += c * mv; }, vars_view.GetConstColumnView(species_view[k]), phi);
            }
            params_view.ForEachRowStrict(
                [](double& vt, const double& vp) { vt = (vt > 0.0) ? vp / vt : 1.0; }, phi, V_phase);
          },
          example_params,
          example_vars,
          example_result);
    provider.ComputeValue = [storage, value_function](
                                const DenseMatrixPolicy& params,
                                const DenseMatrixPolicy& vars,
                                DenseMatrixPolicy& result) mutable { value_function(params, vars, result); };
    if (species_indices.empty())
    {
      provider.ComputeValueAndDerivatives = [compute_value = provider.ComputeValue](
                                                const DenseMatrixPolicy& params,
                                                const DenseMatrixPolicy& vars,
                                                DenseMatrixPolicy& result,
                                                DenseMatrixPolicy& /*partials*/) { compute_value(params, vars, result); };
      return provider;
    }
    auto value_and_derivatives_function =
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& params_view,
              const typename DenseMatrixPolicy::ConstViewType& vars_view,
              const typename DenseMatrixPolicy::ViewType& result_view,
              const typename DenseMatrixPolicy::ViewType& partials_view) {
            auto phi = result_view.GetColumnView(0);
            auto V_phase = result_view.GetRowVariable();
            auto V_total = result_view.GetRowVariable();
            params_view.ForEachRowStrict(
                [](double& vt, double& vp)
                {
                  vt = 0.0;
                  vp = 0.0;
                },
                V_total,
                V_phase);
            const std::size_t n = species_view.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = volumes_view[k];
              if (k < phase_species_count)
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt, double& vp)
                    {
                      const double vol = c * mv;
                      vt += vol;
                      vp += vol;
                    },
                    vars_view.GetConstColumnView(species_view[k]),
                    V_total,
                    V_phase);
              else
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt) { vt += c * mv; },
                    vars_view.GetConstColumnView(species_view[k]),
                    V_total);
            }
            params_view.ForEachRowStrict(
                [](const double& vp, const double& vt, double& phi_out) { phi_out = (vt > 0.0) ? vp / vt : 1.0; },
                V_phase,
                V_total,
                phi);
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = volumes_view[k];
              if (k < phase_species_count)
                params_view.ForEachRowStrict(
                    [mv](const double& phi_row, const double& vt, double& dphi)
                    { dphi = (vt > 0.0) ? mv * (1.0 - phi_row) / vt : 0.0; },
                    phi,
                    V_total,
                    partials_view.GetColumnView(k));
              else
                params_view.ForEachRowStrict(
                    [mv](const double& phi_row, const double& vt, double& dphi)
                    { dphi = (vt > 0.0) ? -mv * phi_row / vt : 0.0; },
                    phi,
                    V_total,
                    partials_view.GetColumnView(k));
            }
          },
          example_params,
          example_vars,
          example_result,
          example_partials);
    provider.ComputeValueAndDerivatives = [storage, value_and_derivatives_function](
                                              const DenseMatrixPolicy& params,
                                              const DenseMatrixPolicy& vars,
                                              DenseMatrixPolicy& result,
                                              DenseMatrixPolicy& partials) mutable
    { value_and_derivatives_function(params, vars, result, partials); };
    return provider;
  }
}  // namespace miam
