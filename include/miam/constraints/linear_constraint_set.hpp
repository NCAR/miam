// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/linear_constraint.hpp>

#include <micm/util/types.hpp>

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `LinearConstraint`, mirroring `micm::ProcessSet`.
  /// @details Indices/coefficients are stored in `SparseMatrixPolicy::VectorType<T>` so they
  ///          are device-accessible via `.GetView()`. A POD `Views` bundle is captured by
  ///          value in every `MICM_LAMBDA`.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class LinearConstraintSet
  {
   public:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;
    template<class U>
    using VectorView = typename Vector<U>::ConstViewType;

    struct Views
    {
      VectorView<micm::Index> alg_indices_;
      VectorView<micm::Index> param_indices_;
      VectorView<micm::Index> counts_per_instance_;
      VectorView<micm::Index> flat_term_indices_;
      VectorView<micm::Real> flat_term_coeffs_;
      VectorView<micm::Index> flat_jac_ids_;
      micm::Index num_instances_;
      micm::Real constant_;

      Views() = default;

      Views(
          const Vector<micm::Index>& alg_indices,
          const Vector<micm::Index>& param_indices,
          const Vector<micm::Index>& counts_per_instance,
          const Vector<micm::Index>& flat_term_indices,
          const Vector<micm::Real>& flat_term_coeffs,
          const Vector<micm::Index>& flat_jac_ids,
          micm::Index num_instances,
          micm::Real constant)
          : alg_indices_(alg_indices.GetView()),
            param_indices_(param_indices.GetView()),
            counts_per_instance_(counts_per_instance.GetView()),
            flat_term_indices_(flat_term_indices.GetView()),
            flat_term_coeffs_(flat_term_coeffs.GetView()),
            flat_jac_ids_(flat_jac_ids.GetView()),
            num_instances_(num_instances),
            constant_(constant)
      {
      }
    };

    LinearConstraintSet() = default;

    LinearConstraintSet(
        const LinearConstraint& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      diagnose_from_state_ = config.diagnose_from_state_;
      const bool is_global = (phase_prefixes.find(config.algebraic_phase_.name_) == phase_prefixes.end());

      std::vector<micm::Index> alg_indices_host;
      std::vector<micm::Index> param_indices_host;
      std::vector<micm::Index> counts_per_instance_host;
      std::vector<micm::Index> flat_term_indices_host;
      std::vector<micm::Real> flat_term_coeffs_host;
      std::vector<micm::Index> flat_jac_ids_host;

      if (is_global)
      {
        auto resolved = ResolveGlobalTerms(config, phase_prefixes, state_variable_indices);
        const auto alg_idx = static_cast<micm::Index>(state_variable_indices.at(config.algebraic_species_.name_));
        alg_indices_host.push_back(alg_idx);
        counts_per_instance_host.push_back(static_cast<micm::Index>(resolved.size()));
        for (const auto& [idx, coeff] : resolved)
        {
          flat_term_indices_host.push_back(static_cast<micm::Index>(idx));
          flat_term_coeffs_host.push_back(coeff);
          flat_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, alg_idx, idx)));
        }
        if (config.diagnose_from_state_)
          param_indices_host.push_back(
              static_cast<micm::Index>(state_parameter_indices.at("LC_" + config.uuid_ + "_constant")));
      }
      else
      {
        const auto& alg_prefixes = phase_prefixes.at(config.algebraic_phase_.name_);
        for (const auto& alg_prefix : alg_prefixes)
        {
          const auto alg_idx = static_cast<micm::Index>(state_variable_indices.at(
              alg_prefix + "." + config.algebraic_phase_.name_ + "." + config.algebraic_species_.name_));
          alg_indices_host.push_back(alg_idx);
          if (config.diagnose_from_state_)
            param_indices_host.push_back(static_cast<micm::Index>(
                state_parameter_indices.at("LC_" + config.uuid_ + "_" + alg_prefix + "_constant")));

          micm::Index count = 0;
          for (const auto& term : config.terms_)
          {
            auto phase_it = phase_prefixes.find(term.phase.name_);
            const std::size_t idx =
                (phase_it != phase_prefixes.end())
                    ? state_variable_indices.at(alg_prefix + "." + term.phase.name_ + "." + term.species.name_)
                    : state_variable_indices.at(term.species.name_);
            flat_term_indices_host.push_back(static_cast<micm::Index>(idx));
            flat_term_coeffs_host.push_back(term.coefficient);
            flat_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, alg_idx, idx)));
            ++count;
          }
          counts_per_instance_host.push_back(count);
        }
      }

      alg_indices_ = Vector<micm::Index>(std::move(alg_indices_host));
      param_indices_ = Vector<micm::Index>(std::move(param_indices_host));
      counts_per_instance_ = Vector<micm::Index>(std::move(counts_per_instance_host));
      flat_term_indices_ = Vector<micm::Index>(std::move(flat_term_indices_host));
      flat_term_coeffs_ = Vector<micm::Real>(std::move(flat_term_coeffs_host));
      flat_jac_ids_ = Vector<micm::Index>(std::move(flat_jac_ids_host));

      alg_indices_.CopyToDevice();
      param_indices_.CopyToDevice();
      counts_per_instance_.CopyToDevice();
      flat_term_indices_.CopyToDevice();
      flat_term_coeffs_.CopyToDevice();
      flat_jac_ids_.CopyToDevice();

      views_ = Views(
          alg_indices_,
          param_indices_,
          counts_per_instance_,
          flat_term_indices_,
          flat_term_coeffs_,
          flat_jac_ids_,
          static_cast<micm::Index>(alg_indices_.size()),
          config.constant_);
    }

    /// @brief Add G(y) into the algebraic rows of `residual`.
    /// @details G = sum(coeff_i * [species_i]) - C (or - state_parameters[C_idx] when diagnosed).
    void AddResidual(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& residual) const
    {
      const auto& views = views_;
      if (diagnose_from_state_)
      {
        DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ViewType& residual_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ConstViewType& params_view)
            {
              micm::Index term_offset = 0;
              for (micm::Index i = 0; i < views.num_instances_; ++i)
              {
                const micm::Index count = views.counts_per_instance_[i];
                const micm::Index alg_idx = views.alg_indices_[i];
                const micm::Index param_idx = views.param_indices_[i];
                auto sum = residual_view.GetRowVariable();
                residual_view.ForEachRowStrict(
                    [](const micm::Real& p, micm::Real& s) { s = -p; },
                    params_view.GetConstColumnView(param_idx),
                    sum);
                for (micm::Index k = 0; k < count; ++k)
                {
                  const micm::Index term_idx = views.flat_term_indices_[term_offset + k];
                  const micm::Real coeff = views.flat_term_coeffs_[term_offset + k];
                  residual_view.ForEachRowStrict(
                      [coeff](const micm::Real& v, micm::Real& s) { s += coeff * v; },
                      state_view.GetConstColumnView(term_idx),
                      sum);
                }
                residual_view.ForEachRowStrict(
                    [](const micm::Real& s, micm::Real& r) { r = s; },
                    sum,
                    residual_view.GetColumnView(alg_idx));
                term_offset += count;
              }
            },
            residual,
            state_variables,
            state_parameters)(residual, state_variables, state_parameters);
      }
      else
      {
        DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ViewType& residual_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view)
            {
              micm::Index term_offset = 0;
              for (micm::Index i = 0; i < views.num_instances_; ++i)
              {
                const micm::Index count = views.counts_per_instance_[i];
                const micm::Index alg_idx = views.alg_indices_[i];
                const micm::Real c = views.constant_;
                auto sum = residual_view.GetRowVariable();
                residual_view.ForEachRowStrict([c](micm::Real& s) { s = -c; }, sum);
                for (micm::Index k = 0; k < count; ++k)
                {
                  const micm::Index term_idx = views.flat_term_indices_[term_offset + k];
                  const micm::Real coeff = views.flat_term_coeffs_[term_offset + k];
                  residual_view.ForEachRowStrict(
                      [coeff](const micm::Real& v, micm::Real& s) { s += coeff * v; },
                      state_view.GetConstColumnView(term_idx),
                      sum);
                }
                residual_view.ForEachRowStrict(
                    [](const micm::Real& s, micm::Real& r) { r = s; },
                    sum,
                    residual_view.GetColumnView(alg_idx));
                term_offset += count;
              }
            },
            residual,
            state_variables)(residual, state_variables);
        (void)state_parameters;
      }
    }

    /// @brief Subtract dG/dy from `jacobian`. dG/d[species_i] = coeff_i.
    void SubtractJacobian(
        const DenseMatrixPolicy& /*state_variables*/,
        const DenseMatrixPolicy& /*state_parameters*/,
        SparseMatrixPolicy& jacobian) const
    {
      const auto& views = views_;
      SparseMatrixPolicy::Function(
          MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jac_view)
          {
            const micm::Index total = static_cast<micm::Index>(views.flat_jac_ids_.size());
            for (micm::Index k = 0; k < total; ++k)
            {
              const micm::Real coeff = views.flat_term_coeffs_[k];
              const micm::Index vec_idx = views.flat_jac_ids_[k];
              jac_view.ForEachBlockStrict(
                  [coeff](micm::Real& j) { j -= coeff; }, jac_view.GetBlockView(vec_idx));
            }
          },
          jacobian)(jacobian);
    }

   private:
    static std::vector<std::pair<std::size_t, double>> ResolveGlobalTerms(
        const LinearConstraint& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices)
    {
      std::vector<std::pair<std::size_t, double>> resolved;
      for (const auto& term : config.terms_)
      {
        auto phase_it = phase_prefixes.find(term.phase.name_);
        if (phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : phase_it->second)
          {
            std::size_t idx = state_variable_indices.at(prefix + "." + term.phase.name_ + "." + term.species.name_);
            resolved.push_back({ idx, term.coefficient });
          }
        }
        else
        {
          std::size_t idx = state_variable_indices.at(term.species.name_);
          resolved.push_back({ idx, term.coefficient });
        }
      }
      return resolved;
    }

    Vector<micm::Index> alg_indices_{};
    Vector<micm::Index> param_indices_{};
    Vector<micm::Index> counts_per_instance_{};
    Vector<micm::Index> flat_term_indices_{};
    Vector<micm::Real> flat_term_coeffs_{};
    Vector<micm::Index> flat_jac_ids_{};
    Views views_{};
    bool diagnose_from_state_{ false };  //!< host-only; selects AddResidual kernel variant
  };
}  // namespace miam
