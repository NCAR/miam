// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/linear_constraint.hpp>

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `LinearConstraint`, mirroring the process Set pattern.
  /// @details Non-templated data-holding class populated once at Finalize time. Solve-time
  ///          methods build the DP/SP kernel per call from trivially-copyable captures.
  ///          Handles all four (is_global, diagnose_from_state) combinations at construction
  ///          time; solve-time methods just iterate the resolved indices.
  class LinearConstraintSet
  {
   public:
    LinearConstraintSet() = default;

    template<typename SparseMatrixPolicy>
    LinearConstraintSet(
        const LinearConstraint& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      is_global_ = (phase_prefixes.find(config.algebraic_phase_.name_) == phase_prefixes.end());
      diagnose_from_state_ = config.diagnose_from_state_;
      constant_ = config.constant_;

      if (is_global_)
      {
        auto resolved = ResolveGlobalTerms(config, phase_prefixes, state_variable_indices);
        alg_indices_.push_back(state_variable_indices.at(config.algebraic_species_.name_));
        counts_per_instance_.push_back(resolved.size());
        for (const auto& [idx, coeff] : resolved)
        {
          flat_term_indices_.push_back(idx);
          flat_term_coeffs_.push_back(coeff);
          flat_jac_ids_.push_back(jacobian.VectorIndex(0, alg_indices_.back(), idx));
        }
        if (diagnose_from_state_)
          param_indices_.push_back(state_parameter_indices.at("LC_" + config.uuid_ + "_constant"));
      }
      else
      {
        const auto& alg_prefixes = phase_prefixes.at(config.algebraic_phase_.name_);
        for (const auto& alg_prefix : alg_prefixes)
        {
          std::size_t alg_idx = state_variable_indices.at(
              alg_prefix + "." + config.algebraic_phase_.name_ + "." + config.algebraic_species_.name_);
          alg_indices_.push_back(alg_idx);
          if (diagnose_from_state_)
            param_indices_.push_back(state_parameter_indices.at("LC_" + config.uuid_ + "_" + alg_prefix + "_constant"));

          std::size_t count = 0;
          for (const auto& term : config.terms_)
          {
            auto phase_it = phase_prefixes.find(term.phase.name_);
            std::size_t idx = (phase_it != phase_prefixes.end())
                                  ? state_variable_indices.at(alg_prefix + "." + term.phase.name_ + "." + term.species.name_)
                                  : state_variable_indices.at(term.species.name_);
            flat_term_indices_.push_back(idx);
            flat_term_coeffs_.push_back(term.coefficient);
            flat_jac_ids_.push_back(jacobian.VectorIndex(0, alg_idx, idx));
            ++count;
          }
          counts_per_instance_.push_back(count);
        }
      }
    }

    /// @brief Add G(y) into the algebraic rows of `residual`.
    /// @details G = sum(coeff_i * [species_i]) - C, written into residual[alg_row].
    template<typename DenseMatrixPolicy>
    void AddResidual(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& residual) const
    {
      const std::size_t n_inst = alg_indices_.size();
      std::size_t offset = 0;
      for (std::size_t i = 0; i < n_inst; ++i)
      {
        const std::size_t count = counts_per_instance_[i];
        const std::size_t alg_idx = alg_indices_[i];
        if (diagnose_from_state_)
        {
          const std::size_t param_idx = param_indices_[i];
          DenseMatrixPolicy::Function(
              [this, offset, count, alg_idx, param_idx](auto&& sv, auto&& sp, auto&& res)
              {
                auto sum = res.GetRowVariable();
                res.ForEachRow(
                    [](const double& p, double& s) { s = -p; }, sp.GetConstColumnView(param_idx), sum);
                for (std::size_t k = 0; k < count; ++k)
                {
                  const std::size_t term_idx = flat_term_indices_[offset + k];
                  const double coeff = flat_term_coeffs_[offset + k];
                  res.ForEachRow(
                      [coeff](const double& val, double& s) { s += coeff * val; },
                      sv.GetConstColumnView(term_idx),
                      sum);
                }
                res.ForEachRow([](const double& s, double& r) { r = s; }, sum, res.GetColumnView(alg_idx));
              },
              state_variables,
              state_parameters,
              residual)(state_variables, state_parameters, residual);
        }
        else
        {
          const double constant = constant_;
          DenseMatrixPolicy::Function(
              [this, offset, count, alg_idx, constant](auto&& sv, auto&& res)
              {
                auto sum = res.GetRowVariable();
                res.ForEachRow([constant](double& s) { s = -constant; }, sum);
                for (std::size_t k = 0; k < count; ++k)
                {
                  const std::size_t term_idx = flat_term_indices_[offset + k];
                  const double coeff = flat_term_coeffs_[offset + k];
                  res.ForEachRow(
                      [coeff](const double& val, double& s) { s += coeff * val; },
                      sv.GetConstColumnView(term_idx),
                      sum);
                }
                res.ForEachRow([](const double& s, double& r) { r = s; }, sum, res.GetColumnView(alg_idx));
              },
              state_variables,
              residual)(state_variables, residual);
        }
        offset += count;
      }
    }

    /// @brief Subtract dG/dy from `jacobian` (jac -= dG/dy). dG/d[species_i] = coeff_i.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractJacobian(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& /*state_parameters*/,
        SparseMatrixPolicy& jacobian) const
    {
      const std::size_t total = flat_jac_ids_.size();
      SparseMatrixPolicy::Function(
          [this, total](auto&& /*sv*/, auto&& jac)
          {
            for (std::size_t k = 0; k < total; ++k)
            {
              const double coeff = flat_term_coeffs_[k];
              const std::size_t vec_idx = flat_jac_ids_[k];
              auto bv = jac.GetBlockView(vec_idx);
              jac.ForEachBlock([coeff](double& j) { j -= coeff; }, bv);
            }
          },
          state_variables,
          jacobian)(state_variables, jacobian);
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

    bool is_global_{ true };
    bool diagnose_from_state_{ false };
    double constant_{ 0.0 };
    std::vector<std::size_t> alg_indices_{};        ///< One per constraint instance (size 1 if global).
    std::vector<std::size_t> param_indices_{};      ///< Diagnosed-constant param index per instance; empty if not diagnosed.
    std::vector<std::size_t> counts_per_instance_{};///< Number of terms per instance.
    std::vector<std::size_t> flat_term_indices_{};  ///< Flattened state-variable indices across all instances.
    std::vector<double> flat_term_coeffs_{};        ///< Coefficient for each flat term.
    std::vector<std::size_t> flat_jac_ids_{};       ///< VectorIndex per flat term for Jacobian (block 0).
  };
}  // namespace miam
