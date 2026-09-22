// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/representations/aerosol_property.hpp>

#include <micm/util/types.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Phase volume fraction: `phi = V_target_phase / V_total_across_all_phases`.
  /// @details The formula depends only on a per-mode list of species indices and molar
  ///          volumes plus a split point separating target-phase species from other-phase
  ///          species, so every representation (SingleMomentMode, TwoMomentMode, UniformSection)
  ///          instantiates this same descriptor type.
  template<typename DenseMatrixPolicy>
  class PhaseVolumeFractionDescriptor
  {
   public:
    template<typename U>
    using Vector = typename DenseMatrixPolicy::template VectorType<U>;

    /// @brief Trivially-copyable snapshot suitable for capture by value into `MICM_LAMBDA`.
    struct View
    {
      typename Vector<std::size_t>::ConstViewType species_variable_indices_{};
      typename Vector<double>::ConstViewType species_molar_volumes_{};
      /// Number of species at the start of the two views that belong to the target phase.
      std::size_t phase_species_count_ = 0;
    };

    PhaseVolumeFractionDescriptor() = default;

    PhaseVolumeFractionDescriptor(
        std::vector<std::size_t> species_variable_indices,
        std::vector<double> species_molar_volumes,
        std::size_t phase_species_count,
        std::vector<std::size_t> dependent_variable_indices)
        : species_variable_indices_(std::move(species_variable_indices)),
          species_molar_volumes_(std::move(species_molar_volumes)),
          phase_species_count_(phase_species_count),
          dependent_variable_indices_(std::move(dependent_variable_indices))
    {
      species_variable_indices_.CopyToDevice();
      species_molar_volumes_.CopyToDevice();
    }

    View GetView() const
    {
      return { species_variable_indices_.GetView(), species_molar_volumes_.GetView(), phase_species_count_ };
    }

    const std::vector<std::size_t>& DependentVariableIndices() const
    {
      return dependent_variable_indices_;
    }

    void Evaluate(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result) const
    {
      const auto view = GetView();
      DenseMatrixPolicy::Function(
          [view](auto&& params_view, auto&& vars_view, auto&& result_view)
          {
            auto phi = result_view.GetColumnView(0);
            if (view.species_variable_indices_.size() == 0)
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
            const std::size_t n = view.species_variable_indices_.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              if (k < view.phase_species_count_)
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt, double& vp)
                    {
                      const double vol = c * mv;
                      vt += vol;
                      vp += vol;
                    },
                    vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                    phi,
                    V_phase);
              else
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt) { vt += c * mv; },
                    vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                    phi);
            }
            params_view.ForEachRowStrict(
                [](double& vt, const double& vp) { vt = (vt > 0.0) ? vp / vt : 1.0; }, phi, V_phase);
          },
          state_parameters,
          state_variables,
          result)(state_parameters, state_variables, result);
    }

    void EvaluateAndDerivatives(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result,
        DenseMatrixPolicy& partials) const
    {
      if (dependent_variable_indices_.empty())
      {
        Evaluate(state_parameters, state_variables, result);
        return;
      }
      const auto view = GetView();
      DenseMatrixPolicy::Function(
          [view](auto&& params_view, auto&& vars_view, auto&& result_view, auto&& partials_view)
          {
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
            const std::size_t n = view.species_variable_indices_.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              if (k < view.phase_species_count_)
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt, double& vp)
                    {
                      const double vol = c * mv;
                      vt += vol;
                      vp += vol;
                    },
                    vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                    V_total,
                    V_phase);
              else
                params_view.ForEachRowStrict(
                    [mv](const double& c, double& vt) { vt += c * mv; },
                    vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                    V_total);
            }
            params_view.ForEachRowStrict(
                [](const double& vp, const double& vt, double& phi_out) { phi_out = (vt > 0.0) ? vp / vt : 1.0; },
                V_phase,
                V_total,
                phi);
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              if (k < view.phase_species_count_)
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
          state_parameters,
          state_variables,
          result,
          partials)(state_parameters, state_variables, result, partials);
    }

   private:
    Vector<std::size_t> species_variable_indices_{};
    Vector<double> species_molar_volumes_{};
    std::size_t phase_species_count_ = 0;
    std::vector<std::size_t> dependent_variable_indices_{};
  };
}  // namespace miam
