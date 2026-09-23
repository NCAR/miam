// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/util/types.hpp>

#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  /// @brief Solve-time companion to `DissolvedEquilibriumConstraint`, mirroring `micm::ProcessSet`.
  /// @details G = k_eq * [S] * prod([R_i]) / ([S]+eps)^n_r - [S] * prod([P_j]) / ([S]+eps)^n_p.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class DissolvedEquilibriumConstraintSet
  {
   public:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;
    template<class U>
    using VectorView = typename Vector<U>::ConstViewType;

    struct Views
    {
      VectorView<micm::Index> algebraic_indices_;
      VectorView<micm::Index> solvent_indices_;
      VectorView<micm::Index> k_eq_indices_;
      VectorView<micm::Index> reactant_indices_;
      VectorView<micm::Index> product_indices_;
      VectorView<micm::Index> reactant_jac_ids_;
      VectorView<micm::Index> product_jac_ids_;
      VectorView<micm::Index> solvent_jac_ids_;
      micm::Index num_phases_;
      micm::Index n_reactants_;
      micm::Index n_products_;
      micm::Real solvent_floor_;

      Views() = default;

      Views(
          const Vector<micm::Index>& algebraic_indices,
          const Vector<micm::Index>& solvent_indices,
          const Vector<micm::Index>& k_eq_indices,
          const Vector<micm::Index>& reactant_indices,
          const Vector<micm::Index>& product_indices,
          const Vector<micm::Index>& reactant_jac_ids,
          const Vector<micm::Index>& product_jac_ids,
          const Vector<micm::Index>& solvent_jac_ids,
          micm::Index num_phases,
          micm::Index n_reactants,
          micm::Index n_products,
          micm::Real solvent_floor)
          : algebraic_indices_(algebraic_indices.GetView()),
            solvent_indices_(solvent_indices.GetView()),
            k_eq_indices_(k_eq_indices.GetView()),
            reactant_indices_(reactant_indices.GetView()),
            product_indices_(product_indices.GetView()),
            reactant_jac_ids_(reactant_jac_ids.GetView()),
            product_jac_ids_(product_jac_ids.GetView()),
            solvent_jac_ids_(solvent_jac_ids.GetView()),
            num_phases_(num_phases),
            n_reactants_(n_reactants),
            n_products_(n_products),
            solvent_floor_(solvent_floor)
      {
      }
    };

    DissolvedEquilibriumConstraintSet() = default;

    DissolvedEquilibriumConstraintSet(
        const DissolvedEquilibriumConstraint& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      const micm::Index nr = static_cast<micm::Index>(config.reactants_.size());
      const micm::Index np = static_cast<micm::Index>(config.products_.size());
      const micm::Real eps = config.solvent_floor_;

      auto phase_it = phase_prefixes.find(config.phase_.name_);
      if (phase_it == phase_prefixes.end())
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedEquilibriumConstraintSet: phase " + config.phase_.name_ + " not found");

      const auto num_phases = static_cast<micm::Index>(phase_it->second.size());

      std::vector<micm::Index> algebraic_indices_host;
      std::vector<micm::Index> solvent_indices_host;
      std::vector<micm::Index> k_eq_indices_host;
      std::vector<micm::Index> reactant_indices_host;
      std::vector<micm::Index> product_indices_host;
      std::vector<micm::Index> reactant_jac_ids_host;
      std::vector<micm::Index> product_jac_ids_host;
      std::vector<micm::Index> solvent_jac_ids_host;

      for (const auto& prefix : phase_it->second)
      {
        k_eq_indices_host.push_back(static_cast<micm::Index>(
            state_parameter_indices.at(prefix + "." + config.phase_.name_ + "." + config.uuid_ + ".k_eq")));
        const auto alg_idx = static_cast<micm::Index>(state_variable_indices.at(
            prefix + "." + config.phase_.name_ + "." + config.algebraic_species_.name_));
        algebraic_indices_host.push_back(alg_idx);
        const auto sol_idx = static_cast<micm::Index>(
            state_variable_indices.at(prefix + "." + config.phase_.name_ + "." + config.solvent_.name_));
        solvent_indices_host.push_back(sol_idx);

        for (const auto& r : config.reactants_)
        {
          const auto idx =
              static_cast<micm::Index>(state_variable_indices.at(prefix + "." + config.phase_.name_ + "." + r.name_));
          reactant_indices_host.push_back(idx);
          reactant_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, alg_idx, idx)));
        }
        for (const auto& p : config.products_)
        {
          const auto idx =
              static_cast<micm::Index>(state_variable_indices.at(prefix + "." + config.phase_.name_ + "." + p.name_));
          product_indices_host.push_back(idx);
          product_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, alg_idx, idx)));
        }
        solvent_jac_ids_host.push_back(static_cast<micm::Index>(jacobian.VectorIndex(0, alg_idx, sol_idx)));
      }

      algebraic_indices_ = Vector<micm::Index>(std::move(algebraic_indices_host));
      solvent_indices_ = Vector<micm::Index>(std::move(solvent_indices_host));
      k_eq_indices_ = Vector<micm::Index>(std::move(k_eq_indices_host));
      reactant_indices_ = Vector<micm::Index>(std::move(reactant_indices_host));
      product_indices_ = Vector<micm::Index>(std::move(product_indices_host));
      reactant_jac_ids_ = Vector<micm::Index>(std::move(reactant_jac_ids_host));
      product_jac_ids_ = Vector<micm::Index>(std::move(product_jac_ids_host));
      solvent_jac_ids_ = Vector<micm::Index>(std::move(solvent_jac_ids_host));

      algebraic_indices_.CopyToDevice();
      solvent_indices_.CopyToDevice();
      k_eq_indices_.CopyToDevice();
      reactant_indices_.CopyToDevice();
      product_indices_.CopyToDevice();
      reactant_jac_ids_.CopyToDevice();
      product_jac_ids_.CopyToDevice();
      solvent_jac_ids_.CopyToDevice();

      views_ = Views(
          algebraic_indices_,
          solvent_indices_,
          k_eq_indices_,
          reactant_indices_,
          product_indices_,
          reactant_jac_ids_,
          product_jac_ids_,
          solvent_jac_ids_,
          num_phases,
          nr,
          np,
          eps);
    }

    /// @brief Add G(y) into the algebraic rows of `residual`.
    void AddResidual(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& residual) const
    {
      const auto& views = views_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ViewType& residual_view,
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ConstViewType& params_view)
          {
            const micm::Index nr = views.n_reactants_;
            const micm::Index np = views.n_products_;
            const micm::Real eps = views.solvent_floor_;
            for (micm::Index i_phase = 0; i_phase < views.num_phases_; ++i_phase)
            {
              const micm::Index alg_idx = views.algebraic_indices_[i_phase];
              const micm::Index solvent_idx = views.solvent_indices_[i_phase];
              const micm::Index k_eq_idx = views.k_eq_indices_[i_phase];
              const micm::Index r_base = i_phase * nr;
              const micm::Index p_base = i_phase * np;

              auto forward = residual_view.GetRowVariable();
              residual_view.ForEachRowStrict(
                  [](const micm::Real& keq, micm::Real& fwd) { fwd = keq; },
                  params_view.GetConstColumnView(k_eq_idx),
                  forward);
              for (micm::Index r = 0; r < nr; ++r)
                residual_view.ForEachRowStrict(
                    [](const micm::Real& conc, micm::Real& fwd) { fwd *= conc; },
                    state_view.GetConstColumnView(views.reactant_indices_[r_base + r]),
                    forward);
              residual_view.ForEachRowStrict(
                  [nr, eps](const micm::Real& sol, micm::Real& fwd)
                  { fwd *= sol / std::pow(sol + eps, static_cast<micm::Real>(nr)); },
                  state_view.GetConstColumnView(solvent_idx),
                  forward);

              auto reverse = residual_view.GetRowVariable();
              residual_view.ForEachRowStrict([](micm::Real& rev) { rev = 1.0; }, reverse);
              for (micm::Index p = 0; p < np; ++p)
                residual_view.ForEachRowStrict(
                    [](const micm::Real& conc, micm::Real& rev) { rev *= conc; },
                    state_view.GetConstColumnView(views.product_indices_[p_base + p]),
                    reverse);
              residual_view.ForEachRowStrict(
                  [np, eps](const micm::Real& sol, micm::Real& rev)
                  { rev *= sol / std::pow(sol + eps, static_cast<micm::Real>(np)); },
                  state_view.GetConstColumnView(solvent_idx),
                  reverse);

              residual_view.ForEachRowStrict(
                  [](const micm::Real& fwd, const micm::Real& rev, micm::Real& r) { r = fwd - rev; },
                  forward,
                  reverse,
                  residual_view.GetColumnView(alg_idx));
            }
          },
          residual,
          state_variables,
          state_parameters)(residual, state_variables, state_parameters);
    }

    /// @brief Subtract dG/dy from `jacobian`.
    void SubtractJacobian(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        SparseMatrixPolicy& jacobian) const
    {
      const auto& views = views_;
      SparseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ConstViewType& params_view,
              const typename SparseMatrixPolicy::ViewType& jac_view)
          {
            const micm::Index nr = views.n_reactants_;
            const micm::Index np = views.n_products_;
            const micm::Real eps = views.solvent_floor_;
            for (micm::Index i_phase = 0; i_phase < views.num_phases_; ++i_phase)
            {
              const micm::Index solvent_idx = views.solvent_indices_[i_phase];
              const micm::Index k_eq_idx = views.k_eq_indices_[i_phase];
              const micm::Index solvent_jac_id = views.solvent_jac_ids_[i_phase];
              const micm::Index r_base = i_phase * nr;
              const micm::Index p_base = i_phase * np;

              // dG/d[R_i] = k_eq * prod_{j!=i}[R_j] * [S] / ([S]+eps)^n_r  →  jac -= dG
              for (micm::Index r = 0; r < nr; ++r)
              {
                auto d = jac_view.GetBlockVariable();
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& keq, micm::Real& v) { v = keq; },
                    params_view.GetConstColumnView(k_eq_idx),
                    d);
                for (micm::Index j = 0; j < nr; ++j)
                  if (j != r)
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& conc, micm::Real& v) { v *= conc; },
                        state_view.GetConstColumnView(views.reactant_indices_[r_base + j]),
                        d);
                jac_view.ForEachBlockStrict(
                    [nr, eps](const micm::Real& sol, micm::Real& v)
                    { v *= sol / std::pow(sol + eps, static_cast<micm::Real>(nr)); },
                    state_view.GetConstColumnView(solvent_idx),
                    d);
                auto bv = jac_view.GetBlockView(views.reactant_jac_ids_[r_base + r]);
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& v, micm::Real& j_val) { j_val -= v; }, d, bv);
              }

              // dG/d[P_j] = -prod_{k!=j}[P_k] * [S] / ([S]+eps)^n_p  →  jac -= -dG = +dG
              for (micm::Index p = 0; p < np; ++p)
              {
                auto d = jac_view.GetBlockVariable();
                jac_view.ForEachBlockStrict([](micm::Real& v) { v = 1.0; }, d);
                for (micm::Index k = 0; k < np; ++k)
                  if (k != p)
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& conc, micm::Real& v) { v *= conc; },
                        state_view.GetConstColumnView(views.product_indices_[p_base + k]),
                        d);
                jac_view.ForEachBlockStrict(
                    [np, eps](const micm::Real& sol, micm::Real& v)
                    { v *= sol / std::pow(sol + eps, static_cast<micm::Real>(np)); },
                    state_view.GetConstColumnView(solvent_idx),
                    d);
                auto bv = jac_view.GetBlockView(views.product_jac_ids_[p_base + p]);
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& v, micm::Real& j_val) { j_val += v; }, d, bv);
              }

              // dG/d[S]: damped solvent derivative
              auto fwd_d = jac_view.GetBlockVariable();
              jac_view.ForEachBlockStrict(
                  [](const micm::Real& keq, micm::Real& v) { v = keq; },
                  params_view.GetConstColumnView(k_eq_idx),
                  fwd_d);
              for (micm::Index r = 0; r < nr; ++r)
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& conc, micm::Real& v) { v *= conc; },
                    state_view.GetConstColumnView(views.reactant_indices_[r_base + r]),
                    fwd_d);
              jac_view.ForEachBlockStrict(
                  [nr, eps](const micm::Real& sol, micm::Real& v)
                  {
                    v *= (eps + (1.0 - static_cast<micm::Real>(nr)) * sol) /
                         std::pow(sol + eps, static_cast<micm::Real>(nr) + 1.0);
                  },
                  state_view.GetConstColumnView(solvent_idx),
                  fwd_d);

              auto rev_d = jac_view.GetBlockVariable();
              jac_view.ForEachBlockStrict([](micm::Real& v) { v = 1.0; }, rev_d);
              for (micm::Index p = 0; p < np; ++p)
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& conc, micm::Real& v) { v *= conc; },
                    state_view.GetConstColumnView(views.product_indices_[p_base + p]),
                    rev_d);
              jac_view.ForEachBlockStrict(
                  [np, eps](const micm::Real& sol, micm::Real& v)
                  {
                    v *= (eps + (1.0 - static_cast<micm::Real>(np)) * sol) /
                         std::pow(sol + eps, static_cast<micm::Real>(np) + 1.0);
                  },
                  state_view.GetConstColumnView(solvent_idx),
                  rev_d);

              auto bv = jac_view.GetBlockView(solvent_jac_id);
              jac_view.ForEachBlockStrict(
                  [](const micm::Real& fwd, const micm::Real& rev, micm::Real& j_val) { j_val -= (fwd - rev); },
                  fwd_d,
                  rev_d,
                  bv);
            }
          },
          state_variables,
          state_parameters,
          jacobian)(state_variables, state_parameters, jacobian);
    }

   private:
    Vector<micm::Index> algebraic_indices_{};
    Vector<micm::Index> solvent_indices_{};
    Vector<micm::Index> k_eq_indices_{};
    Vector<micm::Index> reactant_indices_{};
    Vector<micm::Index> product_indices_{};
    Vector<micm::Index> reactant_jac_ids_{};
    Vector<micm::Index> product_jac_ids_{};
    Vector<micm::Index> solvent_jac_ids_{};
    Views views_{};
  };
}  // namespace miam
