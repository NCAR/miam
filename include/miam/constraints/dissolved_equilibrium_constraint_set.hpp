// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `DissolvedEquilibriumConstraint`.
  /// @details Non-templated data-holding class populated once at Finalize time.
  ///          Residual: G = k_eq * [S] * prod([R_i]) / ([S]+eps)^n_r - [S] * prod([P_j]) / ([S]+eps)^n_p.
  class DissolvedEquilibriumConstraintSet
  {
   public:
    DissolvedEquilibriumConstraintSet() = default;

    template<typename SparseMatrixPolicy>
    DissolvedEquilibriumConstraintSet(
        const DissolvedEquilibriumConstraint& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      n_reactants_ = config.reactants_.size();
      n_products_ = config.products_.size();
      solvent_floor_ = config.solvent_floor_;

      auto phase_it = phase_prefixes.find(config.phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedEquilibriumConstraintSet: phase " + config.phase_.name_ + " not found");

      num_phases_ = phase_it->second.size();
      for (const auto& prefix : phase_it->second)
      {
        k_eq_indices_.push_back(
            state_parameter_indices.at(prefix + "." + config.phase_.name_ + "." + config.uuid_ + ".k_eq"));
        std::size_t alg_idx = state_variable_indices.at(
            prefix + "." + config.phase_.name_ + "." + config.algebraic_species_.name_);
        algebraic_indices_.push_back(alg_idx);
        solvent_indices_.push_back(
            state_variable_indices.at(prefix + "." + config.phase_.name_ + "." + config.solvent_.name_));

        for (const auto& r : config.reactants_)
        {
          std::size_t idx = state_variable_indices.at(prefix + "." + config.phase_.name_ + "." + r.name_);
          reactant_indices_.push_back(idx);
          reactant_jac_ids_.push_back(jacobian.VectorIndex(0, alg_idx, idx));
        }
        for (const auto& p : config.products_)
        {
          std::size_t idx = state_variable_indices.at(prefix + "." + config.phase_.name_ + "." + p.name_);
          product_indices_.push_back(idx);
          product_jac_ids_.push_back(jacobian.VectorIndex(0, alg_idx, idx));
        }
        solvent_jac_ids_.push_back(jacobian.VectorIndex(0, alg_idx, solvent_indices_.back()));
      }
    }

    /// @brief Add G(y) into the algebraic rows of `residual`.
    template<typename DenseMatrixPolicy>
    void AddResidual(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& residual) const
    {
      const std::size_t nr = n_reactants_;
      const std::size_t np = n_products_;
      const double eps = solvent_floor_;
      for (std::size_t i_phase = 0; i_phase < num_phases_; ++i_phase)
      {
        const std::size_t alg_idx = algebraic_indices_[i_phase];
        const std::size_t solvent_idx = solvent_indices_[i_phase];
        const std::size_t k_eq_idx = k_eq_indices_[i_phase];
        const std::size_t r_base = i_phase * nr;
        const std::size_t p_base = i_phase * np;

        DenseMatrixPolicy::Function(
            [this, r_base, p_base, nr, np, eps, alg_idx, solvent_idx, k_eq_idx](
                auto&& sv, auto&& sp, auto&& res)
            {
              auto forward = res.GetRowVariable();
              res.ForEachRow(
                  [](const double& keq, double& fwd) { fwd = keq; },
                  sp.GetConstColumnView(k_eq_idx),
                  forward);
              for (std::size_t r = 0; r < nr; ++r)
                res.ForEachRow(
                    [](const double& conc, double& fwd) { fwd *= conc; },
                    sv.GetConstColumnView(reactant_indices_[r_base + r]),
                    forward);
              res.ForEachRow(
                  [nr, eps](const double& sol, double& fwd)
                  { fwd *= sol / std::pow(sol + eps, static_cast<double>(nr)); },
                  sv.GetConstColumnView(solvent_idx),
                  forward);

              auto reverse = res.GetRowVariable();
              res.ForEachRow([](double& rev) { rev = 1.0; }, reverse);
              for (std::size_t p = 0; p < np; ++p)
                res.ForEachRow(
                    [](const double& conc, double& rev) { rev *= conc; },
                    sv.GetConstColumnView(product_indices_[p_base + p]),
                    reverse);
              res.ForEachRow(
                  [np, eps](const double& sol, double& rev)
                  { rev *= sol / std::pow(sol + eps, static_cast<double>(np)); },
                  sv.GetConstColumnView(solvent_idx),
                  reverse);

              res.ForEachRow(
                  [](const double& fwd, const double& rev, double& r) { r = fwd - rev; },
                  forward,
                  reverse,
                  res.GetColumnView(alg_idx));
            },
            state_variables,
            state_parameters,
            residual)(state_variables, state_parameters, residual);
      }
    }

    /// @brief Subtract dG/dy from `jacobian`.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractJacobian(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        SparseMatrixPolicy& jacobian) const
    {
      const std::size_t nr = n_reactants_;
      const std::size_t np = n_products_;
      const double eps = solvent_floor_;
      for (std::size_t i_phase = 0; i_phase < num_phases_; ++i_phase)
      {
        const std::size_t solvent_idx = solvent_indices_[i_phase];
        const std::size_t k_eq_idx = k_eq_indices_[i_phase];
        const std::size_t solvent_jac_id = solvent_jac_ids_[i_phase];
        const std::size_t r_base = i_phase * nr;
        const std::size_t p_base = i_phase * np;

        SparseMatrixPolicy::Function(
            [this, r_base, p_base, nr, np, eps, solvent_idx, k_eq_idx, solvent_jac_id](
                auto&& sv, auto&& sp, auto&& jac)
            {
              // dG/d[R_i] = k_eq * prod_{j!=i}[R_j] * [S] / ([S]+eps)^n_r  →  jac -= dG
              for (std::size_t r = 0; r < nr; ++r)
              {
                auto d = jac.GetBlockVariable();
                jac.ForEachBlock(
                    [](const double& keq, double& v) { v = keq; },
                    sp.GetConstColumnView(k_eq_idx),
                    d);
                for (std::size_t j = 0; j < nr; ++j)
                  if (j != r)
                    jac.ForEachBlock(
                        [](const double& conc, double& v) { v *= conc; },
                        sv.GetConstColumnView(reactant_indices_[r_base + j]),
                        d);
                jac.ForEachBlock(
                    [nr, eps](const double& sol, double& v)
                    { v *= sol / std::pow(sol + eps, static_cast<double>(nr)); },
                    sv.GetConstColumnView(solvent_idx),
                    d);
                auto bv = jac.GetBlockView(reactant_jac_ids_[r_base + r]);
                jac.ForEachBlock([](const double& v, double& j_val) { j_val -= v; }, d, bv);
              }

              // dG/d[P_j] = -prod_{k!=j}[P_k] * [S] / ([S]+eps)^n_p  →  jac -= -dG = +dG
              for (std::size_t p = 0; p < np; ++p)
              {
                auto d = jac.GetBlockVariable();
                jac.ForEachBlock([](double& v) { v = 1.0; }, d);
                for (std::size_t k = 0; k < np; ++k)
                  if (k != p)
                    jac.ForEachBlock(
                        [](const double& conc, double& v) { v *= conc; },
                        sv.GetConstColumnView(product_indices_[p_base + k]),
                        d);
                jac.ForEachBlock(
                    [np, eps](const double& sol, double& v)
                    { v *= sol / std::pow(sol + eps, static_cast<double>(np)); },
                    sv.GetConstColumnView(solvent_idx),
                    d);
                auto bv = jac.GetBlockView(product_jac_ids_[p_base + p]);
                jac.ForEachBlock([](const double& v, double& j_val) { j_val += v; }, d, bv);
              }

              // dG/d[S]: damped solvent derivative
              auto fwd_d = jac.GetBlockVariable();
              jac.ForEachBlock(
                  [](const double& keq, double& v) { v = keq; }, sp.GetConstColumnView(k_eq_idx), fwd_d);
              for (std::size_t r = 0; r < nr; ++r)
                jac.ForEachBlock(
                    [](const double& conc, double& v) { v *= conc; },
                    sv.GetConstColumnView(reactant_indices_[r_base + r]),
                    fwd_d);
              jac.ForEachBlock(
                  [nr, eps](const double& sol, double& v)
                  {
                    v *= (eps + (1.0 - static_cast<double>(nr)) * sol) /
                         std::pow(sol + eps, static_cast<double>(nr) + 1.0);
                  },
                  sv.GetConstColumnView(solvent_idx),
                  fwd_d);

              auto rev_d = jac.GetBlockVariable();
              jac.ForEachBlock([](double& v) { v = 1.0; }, rev_d);
              for (std::size_t p = 0; p < np; ++p)
                jac.ForEachBlock(
                    [](const double& conc, double& v) { v *= conc; },
                    sv.GetConstColumnView(product_indices_[p_base + p]),
                    rev_d);
              jac.ForEachBlock(
                  [np, eps](const double& sol, double& v)
                  {
                    v *= (eps + (1.0 - static_cast<double>(np)) * sol) /
                         std::pow(sol + eps, static_cast<double>(np) + 1.0);
                  },
                  sv.GetConstColumnView(solvent_idx),
                  rev_d);

              auto bv = jac.GetBlockView(solvent_jac_id);
              jac.ForEachBlock(
                  [](const double& fwd, const double& rev, double& j_val) { j_val -= (fwd - rev); },
                  fwd_d,
                  rev_d,
                  bv);
            },
            state_variables,
            state_parameters,
            jacobian)(state_variables, state_parameters, jacobian);
      }
    }

   private:
    std::size_t num_phases_{ 0 };
    std::size_t n_reactants_{ 0 };
    std::size_t n_products_{ 0 };
    double solvent_floor_{ 1.0e-20 };
    std::vector<std::size_t> algebraic_indices_{};  ///< [num_phases]
    std::vector<std::size_t> solvent_indices_{};    ///< [num_phases]
    std::vector<std::size_t> k_eq_indices_{};       ///< [num_phases]
    std::vector<std::size_t> reactant_indices_{};   ///< [num_phases * n_reactants], row-major
    std::vector<std::size_t> product_indices_{};    ///< [num_phases * n_products],  row-major
    std::vector<std::size_t> reactant_jac_ids_{};   ///< VectorIndex per reactant, [num_phases * n_reactants]
    std::vector<std::size_t> product_jac_ids_{};    ///< VectorIndex per product,  [num_phases * n_products]
    std::vector<std::size_t> solvent_jac_ids_{};    ///< VectorIndex per solvent, [num_phases]
  };
}  // namespace miam
