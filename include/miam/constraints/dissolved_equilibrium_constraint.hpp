// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/types.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  /// @brief A dissolved equilibrium constraint
  /// @details Replaces the ODE row for a designated algebraic species with the
  ///          steady-state equilibrium condition:
  ///
  ///          \f$ G = K_{eq} \frac{\prod[R_i]}{[S]^{n_r - 1}}
  ///                       - \frac{\prod[P_j]}{[S]^{n_p - 1}} = 0 \f$
  ///
  ///          where \f$ K_{eq} = k_f / k_r \f$, \f$ R_i \f$ are reactants, \f$ P_j \f$ are
  ///          products, and \f$ S \f$ is the solvent. One of the product species is designated
  ///          as the algebraic variable whose ODE row is replaced by this constraint.
  ///
  ///          To prevent singularity as \f$[S] \to 0\f$, the solvent denominator is
  ///          regularized by a small floor \f$\delta\f$ (\c solvent_floor_):
  ///
  ///          \f$ G = K_{eq} \frac{[S]\prod[R_i]}{([S]+\delta)^{n_r}}
  ///                       - \frac{[S]\prod[P_j]}{([S]+\delta)^{n_p}} = 0 \f$
  class DissolvedEquilibriumConstraint
  {
   public:
    EquilibriumConstant equilibrium_constant_;  ///< K_eq
    std::vector<micm::Species> reactants_;                                            ///< Reactant species
    std::vector<micm::Species> products_;                                             ///< Product species
    micm::Species algebraic_species_;  ///< Product species whose ODE row is replaced
    micm::Species solvent_;            ///< Solvent species
    micm::Phase phase_;                ///< Phase in which the reaction occurs
    std::string uuid_;                 ///< Unique identifier
    double solvent_floor_{ 1.0e-20 };  ///< Floor \f$\delta\f$ [mol m⁻³] added to \f$[S]\f$ in \f$([S]+\delta)^n\f$
                                       ///< denominator to prevent singularity as \f$[S] \to 0\f$

    DissolvedEquilibriumConstraint() = delete;

    /// @brief Constructor
    DissolvedEquilibriumConstraint(
        EquilibriumConstant equilibrium_constant,
        const std::vector<micm::Species>& reactants,
        const std::vector<micm::Species>& products,
        const micm::Species& algebraic_species,
        micm::Species solvent,
        micm::Phase phase,
        double solvent_floor = 1.0e-20)
        : equilibrium_constant_(equilibrium_constant),
          reactants_(reactants),
          products_(products),
          algebraic_species_(algebraic_species),
          solvent_(solvent),
          phase_(phase),
          uuid_(GenerateUuid()),
          solvent_floor_(solvent_floor)
    {
      // Validate that the algebraic species is one of the products
      bool found = false;
      for (const auto& product : products_)
      {
        if (product.name_ == algebraic_species_.name_)
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_ALGEBRAIC_SPECIES_NOT_FOUND_IN_PRODUCTS,
            "DissolvedEquilibriumConstraint: algebraic species '" + algebraic_species_.name_ +
                "' must be one of the products.");
      }
    }

    /// @brief Create a copy with a new UUID
    DissolvedEquilibriumConstraint CopyWithNewUuid() const
    {
      return DissolvedEquilibriumConstraint(
          equilibrium_constant_, reactants_, products_, algebraic_species_, solvent_, phase_, solvent_floor_);
    }

    /// @brief Returns the names of algebraic variables (one per phase instance)
    std::set<std::string> ConstraintAlgebraicVariableNames(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
        {
          names.insert(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);
        }
      }
      return names;
    }

    /// @brief Returns all species the constraint depends on
    std::set<std::string> ConstraintSpeciesDependencies(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> species_names;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
        {
          for (const auto& reactant : reactants_)
            species_names.insert(prefix + "." + phase_.name_ + "." + reactant.name_);
          for (const auto& product : products_)
            species_names.insert(prefix + "." + phase_.name_ + "." + product.name_);
          species_names.insert(prefix + "." + phase_.name_ + "." + solvent_.name_);
        }
      }
      return species_names;
    }

    /// @brief Returns non-zero constraint Jacobian element positions
    /// @details For each phase instance, the algebraic row depends on all reactants, all products, and the solvent.
    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      std::set<std::pair<std::size_t, std::size_t>> elements;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it == phase_prefixes.end())
        return elements;

      for (const auto& prefix : phase_it->second)
      {
        std::size_t alg_row = state_variable_indices.at(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);

        for (const auto& reactant : reactants_)
        {
          std::size_t col = state_variable_indices.at(prefix + "." + phase_.name_ + "." + reactant.name_);
          elements.insert({ alg_row, col });
        }
        for (const auto& product : products_)
        {
          std::size_t col = state_variable_indices.at(prefix + "." + phase_.name_ + "." + product.name_);
          elements.insert({ alg_row, col });
        }
        std::size_t solvent_col = state_variable_indices.at(prefix + "." + phase_.name_ + "." + solvent_.name_);
        elements.insert({ alg_row, solvent_col });
      }
      return elements;
    }

    /// @brief Returns the names of state parameters owned by this constraint (one per phase instance).
    /// @details Each phase instance writes \f$ K_{eq}(T) \f$ to a dedicated column of the
    ///          state parameter matrix every time conditions change.
    std::set<std::string> ConstraintStateParameterNames(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          names.insert(prefix + "." + phase_.name_ + "." + uuid_ + ".k_eq");
      }
      return names;
    }

    /// @brief Returns a function that writes \f$ K_{eq}(T) \f$ per grid cell into the state parameter matrix.
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateConstraintParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
    {
      std::vector<std::size_t> k_eq_indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          k_eq_indices.push_back(state_parameter_indices.at(prefix + "." + phase_.name_ + "." + uuid_ + ".k_eq"));
      }
      const EquilibriumConstant equilibrium_constant = equilibrium_constant_;

      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      auto storage = std::make_shared<Vector>(k_eq_indices);
      storage->CopyToDevice();
      auto k_eq_view = storage->GetView();
      const std::size_t num_k_eq = k_eq_indices.size();
      DenseMatrixPolicy dummy_params{ 1, state_parameter_indices.size(), 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> dummy_conditions;

      auto function = DenseMatrixPolicy::Function(
              MICM_LAMBDA(
                  const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
                  const typename DenseMatrixPolicy::ViewType& params_view)
          {
                for (std::size_t i = 0; i < num_k_eq; ++i)
                params_view.ForEachRowStrict(
                    [equilibrium_constant](const micm::Conditions& cond, micm::Real& k_eq)
                  { k_eq = Calculate(equilibrium_constant, cond); },
                    conditions_view,
                      params_view.GetColumnView(k_eq_view[i]));
          },
          dummy_conditions,
          dummy_params);

      return [storage, function](
                 const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
                 DenseMatrixPolicy& params) mutable { function(conditions, params); };
    }

    /// @brief Returns a function that computes constraint residuals G(y) = 0
    /// @details G = K_eq * [S] * prod([R_i]) / ([S]+δ)^n_r - [S] * prod([P_j]) / ([S]+δ)^n_p
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      std::size_t n_reactants = reactants_.size();
      std::size_t n_products = products_.size();
      double eps = solvent_floor_;

      std::vector<std::size_t> k_eq_indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          k_eq_indices.push_back(state_parameter_indices.at(prefix + "." + phase_.name_ + "." + uuid_ + ".k_eq"));
      }

      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      Vector reactant_indices(indices.reactant_indices_.AsVector());
      Vector product_indices(indices.product_indices_.AsVector());
      Vector solvent_indices(indices.solvent_indices_);
      Vector algebraic_indices(indices.algebraic_indices_);
      Vector k_eq_indices_vec(k_eq_indices);
      reactant_indices.CopyToDevice();
      product_indices.CopyToDevice();
      solvent_indices.CopyToDevice();
      algebraic_indices.CopyToDevice();
      k_eq_indices_vec.CopyToDevice();
      std::size_t num_phases = indices.number_of_phase_instances_;

      struct Storage
      {
        Vector reactant_indices, product_indices, solvent_indices, algebraic_indices, k_eq_indices_vec;
      };
      auto storage = std::make_shared<Storage>(Storage{ std::move(reactant_indices),
                                                        std::move(product_indices),
                                                        std::move(solvent_indices),
                                                        std::move(algebraic_indices),
                                                        std::move(k_eq_indices_vec) });
      auto reactant_view = storage->reactant_indices.GetView();
      auto product_view = storage->product_indices.GetView();
      auto solvent_view = storage->solvent_indices.GetView();
      auto algebraic_view = storage->algebraic_indices.GetView();
      auto k_eq_view = storage->k_eq_indices_vec.GetView();
      DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_params{ 1, std::max(state_parameter_indices.size(), std::size_t{ 1 }), 0.0 };

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ViewType& residual_view)
          {
              for (std::size_t i_phase = 0; i_phase < num_phases; ++i_phase)
            {
                const std::size_t alg_idx = algebraic_view[i_phase];
                const std::size_t solvent_idx = solvent_view[i_phase];

              // Forward part: K_eq * prod([R_i]) * [S] / ([S]+eps)^n_r
                auto forward = residual_view.GetRowVariable();
                residual_view.ForEachRowStrict(
                    [](const micm::Real& keq, micm::Real& fwd) { fwd = keq; },
                    params_view.GetConstColumnView(k_eq_view[i_phase]),
                  forward);
              for (std::size_t r = 0; r < n_reactants; ++r)
                  residual_view.ForEachRowStrict(
                      [](const micm::Real& conc, micm::Real& fwd) { fwd *= conc; },
                      state_view.GetConstColumnView(reactant_view[i_phase * n_reactants + r]),
                    forward);
                residual_view.ForEachRowStrict(
                    [n_reactants, eps](const micm::Real& sol, micm::Real& fwd)
                    { fwd *= sol / std::pow(sol + eps, static_cast<micm::Real>(n_reactants)); },
                    state_view.GetConstColumnView(solvent_idx),
                  forward);

              // Reverse part: prod([P_j]) * [S] / ([S]+eps)^n_p
                auto reverse = residual_view.GetRowVariable();
                residual_view.ForEachRowStrict([](micm::Real& rev) { rev = 1.0; }, reverse);
              for (std::size_t p = 0; p < n_products; ++p)
                  residual_view.ForEachRowStrict(
                      [](const micm::Real& conc, micm::Real& rev) { rev *= conc; },
                      state_view.GetConstColumnView(product_view[i_phase * n_products + p]),
                    reverse);
                residual_view.ForEachRowStrict(
                    [n_products, eps](const micm::Real& sol, micm::Real& rev)
                    { rev *= sol / std::pow(sol + eps, static_cast<micm::Real>(n_products)); },
                    state_view.GetConstColumnView(solvent_idx),
                  reverse);

              // G = forward - reverse
                residual_view.ForEachRowStrict(
                    [](const micm::Real& fwd, const micm::Real& rev, micm::Real& res) { res = fwd - rev; },
                  forward,
                  reverse,
                    residual_view.GetColumnView(alg_idx));
            }
          },
          dummy_state,
          dummy_params,
          dummy_state);

      return [storage, function](
                 const DenseMatrixPolicy& state_variables,
                 const DenseMatrixPolicy& state_parameters,
                 DenseMatrixPolicy& residual) mutable { function(state_variables, state_parameters, residual); };
    }

    /// @brief Returns a function that computes constraint Jacobian entries (subtracts dG/dy)
    /// @details Follows MICM convention: jac[row][col] -= dG/dy
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      auto indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      std::size_t n_reactants = reactants_.size();
      std::size_t n_products = products_.size();
      double eps = solvent_floor_;

      std::vector<std::size_t> k_eq_indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        for (const auto& prefix : phase_it->second)
          k_eq_indices.push_back(state_parameter_indices.at(prefix + "." + phase_.name_ + "." + uuid_ + ".k_eq"));
      }

      // Pre-compute block-0 VectorIndex values per instance
      std::vector<std::size_t> reactant_jac_ids;
      std::vector<std::size_t> product_jac_ids;
      std::vector<std::size_t> solvent_jac_ids;
      for (std::size_t i_phase = 0; i_phase < indices.number_of_phase_instances_; ++i_phase)
      {
        std::size_t alg_row = indices.algebraic_indices_[i_phase];
        for (std::size_t r = 0; r < n_reactants; ++r)
          reactant_jac_ids.push_back(jacobian.VectorIndex(0, alg_row, indices.reactant_indices_[i_phase][r]));
        for (std::size_t p = 0; p < n_products; ++p)
          product_jac_ids.push_back(jacobian.VectorIndex(0, alg_row, indices.product_indices_[i_phase][p]));
        solvent_jac_ids.push_back(jacobian.VectorIndex(0, alg_row, indices.solvent_indices_[i_phase]));
      }

      using Vector = typename SparseMatrixPolicy::template VectorType<std::size_t>;
      Vector reactant_indices(indices.reactant_indices_.AsVector());
      Vector product_indices(indices.product_indices_.AsVector());
      Vector solvent_indices(indices.solvent_indices_);
      Vector k_eq_indices_vec(k_eq_indices);
      Vector reactant_jac_ids_vec(reactant_jac_ids);
      Vector product_jac_ids_vec(product_jac_ids);
      Vector solvent_jac_ids_vec(solvent_jac_ids);
      reactant_indices.CopyToDevice();
      product_indices.CopyToDevice();
      solvent_indices.CopyToDevice();
      k_eq_indices_vec.CopyToDevice();
      reactant_jac_ids_vec.CopyToDevice();
      product_jac_ids_vec.CopyToDevice();
      solvent_jac_ids_vec.CopyToDevice();
      std::size_t num_phases = indices.number_of_phase_instances_;

      struct Storage
      {
        Vector reactant_indices, product_indices, solvent_indices, k_eq_indices_vec;
        Vector reactant_jac_ids_vec, product_jac_ids_vec, solvent_jac_ids_vec;
      };
      auto storage = std::make_shared<Storage>(Storage{ std::move(reactant_indices),
                                                        std::move(product_indices),
                                                        std::move(solvent_indices),
                                                        std::move(k_eq_indices_vec),
                                                        std::move(reactant_jac_ids_vec),
                                                        std::move(product_jac_ids_vec),
                                                        std::move(solvent_jac_ids_vec) });
      auto reactant_view = storage->reactant_indices.GetView();
      auto product_view = storage->product_indices.GetView();
      auto solvent_view = storage->solvent_indices.GetView();
      auto k_eq_view = storage->k_eq_indices_vec.GetView();
      auto reactant_jac_view = storage->reactant_jac_ids_vec.GetView();
      auto product_jac_view = storage->product_jac_ids_vec.GetView();
      auto solvent_jac_view = storage->solvent_jac_ids_vec.GetView();
      DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_params{ 1, std::max(state_parameter_indices.size(), std::size_t{ 1 }), 0.0 };

      auto function = SparseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename SparseMatrixPolicy::ViewType& jac_view)
          {
              for (std::size_t i_phase = 0; i_phase < num_phases; ++i_phase)
            {
                const std::size_t solvent_idx = solvent_view[i_phase];
                const std::size_t k_eq_idx = k_eq_view[i_phase];
                const std::size_t r_base = i_phase * n_reactants;
                const std::size_t p_base = i_phase * n_products;

              // dG/d[R_i] = K_eq * prod([R_j], j!=i) * [S] / ([S]+eps)^n_r
              for (std::size_t r = 0; r < n_reactants; ++r)
              {
                  auto deriv = jac_view.GetBlockVariable();
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& keq, micm::Real& d) { d = keq; },
                      params_view.GetConstColumnView(k_eq_idx),
                    deriv);
                for (std::size_t j = 0; j < n_reactants; ++j)
                  if (j != r)
                      jac_view.ForEachBlockStrict(
                          [](const micm::Real& conc, micm::Real& d) { d *= conc; },
                          state_view.GetConstColumnView(reactant_view[r_base + j]),
                        deriv);
                  jac_view.ForEachBlockStrict(
                      [n_reactants, eps](const micm::Real& sol, micm::Real& d)
                      { d *= sol / std::pow(sol + eps, static_cast<micm::Real>(n_reactants)); },
                      state_view.GetConstColumnView(solvent_idx),
                    deriv);
                  auto bv = jac_view.GetBlockView(reactant_jac_view[r_base + r]);
                  jac_view.ForEachBlockStrict([](const micm::Real& d, micm::Real& j) { j -= d; }, deriv, bv);
              }

              // dG/d[P_j] = -prod([P_k], k!=j) * [S] / ([S]+eps)^n_p
              for (std::size_t p = 0; p < n_products; ++p)
              {
                  auto deriv = jac_view.GetBlockVariable();
                  jac_view.ForEachBlockStrict([](micm::Real& d) { d = 1.0; }, deriv);
                for (std::size_t k = 0; k < n_products; ++k)
                  if (k != p)
                      jac_view.ForEachBlockStrict(
                          [](const micm::Real& conc, micm::Real& d) { d *= conc; },
                          state_view.GetConstColumnView(product_view[p_base + k]),
                        deriv);
                  jac_view.ForEachBlockStrict(
                      [n_products, eps](const micm::Real& sol, micm::Real& d)
                      { d *= sol / std::pow(sol + eps, static_cast<micm::Real>(n_products)); },
                      state_view.GetConstColumnView(solvent_idx),
                    deriv);
                  auto bv = jac_view.GetBlockView(product_jac_view[p_base + p]);
                  jac_view.ForEachBlockStrict([](const micm::Real& d, micm::Real& j) { j += d; }, deriv, bv);
              }

              // dG/d[S]: damped solvent derivative
              {
                  auto forward_deriv = jac_view.GetBlockVariable();
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& keq, micm::Real& fwd) { fwd = keq; },
                      params_view.GetConstColumnView(k_eq_idx),
                    forward_deriv);
                for (std::size_t r = 0; r < n_reactants; ++r)
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& conc, micm::Real& fwd) { fwd *= conc; },
                        state_view.GetConstColumnView(reactant_view[r_base + r]),
                      forward_deriv);
                  jac_view.ForEachBlockStrict(
                      [n_reactants, eps](const micm::Real& sol, micm::Real& fwd)
                    {
                        fwd *= (eps + (1.0 - static_cast<micm::Real>(n_reactants)) * sol) /
                               std::pow(sol + eps, static_cast<micm::Real>(n_reactants) + 1.0);
                    },
                      state_view.GetConstColumnView(solvent_idx),
                    forward_deriv);

                  auto reverse_deriv = jac_view.GetBlockVariable();
                  jac_view.ForEachBlockStrict([](micm::Real& rev) { rev = 1.0; }, reverse_deriv);
                for (std::size_t p = 0; p < n_products; ++p)
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& conc, micm::Real& rev) { rev *= conc; },
                        state_view.GetConstColumnView(product_view[p_base + p]),
                      reverse_deriv);
                  jac_view.ForEachBlockStrict(
                      [n_products, eps](const micm::Real& sol, micm::Real& rev)
                    {
                        rev *= (eps + (1.0 - static_cast<micm::Real>(n_products)) * sol) /
                               std::pow(sol + eps, static_cast<micm::Real>(n_products) + 1.0);
                    },
                      state_view.GetConstColumnView(solvent_idx),
                    reverse_deriv);

                  auto bv = jac_view.GetBlockView(solvent_jac_view[i_phase]);
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& fwd_d, const micm::Real& rev_d, micm::Real& j) { j -= (fwd_d - rev_d); },
                    forward_deriv,
                    reverse_deriv,
                    bv);
              }
            }
          },
          dummy_state,
          dummy_params,
          jacobian);

      return [storage, function](
                 const DenseMatrixPolicy& state_variables,
                 const DenseMatrixPolicy& state_parameters,
                 SparseMatrixPolicy& jacobian_values) mutable { function(state_variables, state_parameters, jacobian_values); };
    }

   private:
    /// @brief Helper struct for state variable indices across phase instances
    struct StateVariableIndices
    {
      std::size_t number_of_phase_instances_;
      micm::Matrix<std::size_t> reactant_indices_;  // [n_instances x n_reactants]
      micm::Matrix<std::size_t> product_indices_;   // [n_instances x n_products]
      std::vector<std::size_t> solvent_indices_;    // [n_instances]
      std::vector<std::size_t> algebraic_indices_;  // [n_instances]
    };

    /// @brief Helper struct for Jacobian sparse matrix indices
    struct JacobianIndices
    {
      // [n_instances][n_reactants + n_products + 1 (solvent)] per grid row
      std::vector<std::vector<std::size_t>> indices_;
    };

    /// @brief Build state variable indices for all phase instances
    StateVariableIndices GetStateVariableIndices(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      StateVariableIndices indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it == phase_prefixes.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "DissolvedEquilibriumConstraint: Phase " + phase_.name_ + " not found in phase_prefixes");
      }
      const auto& prefixes = phase_it->second;
      indices.number_of_phase_instances_ = prefixes.size();
      indices.reactant_indices_ = micm::Matrix<std::size_t>(prefixes.size(), reactants_.size());
      indices.product_indices_ = micm::Matrix<std::size_t>(prefixes.size(), products_.size());
      indices.solvent_indices_.resize(prefixes.size());
      indices.algebraic_indices_.resize(prefixes.size());

      std::size_t i_phase = 0;
      for (const auto& prefix : prefixes)
      {
        for (std::size_t r = 0; r < reactants_.size(); ++r)
        {
          std::string var = prefix + "." + phase_.name_ + "." + reactants_[r].name_;
          indices.reactant_indices_[i_phase][r] = state_variable_indices.at(var);
        }
        for (std::size_t p = 0; p < products_.size(); ++p)
        {
          std::string var = prefix + "." + phase_.name_ + "." + products_[p].name_;
          indices.product_indices_[i_phase][p] = state_variable_indices.at(var);
        }
        indices.solvent_indices_[i_phase] = state_variable_indices.at(prefix + "." + phase_.name_ + "." + solvent_.name_);
        indices.algebraic_indices_[i_phase] =
            state_variable_indices.at(prefix + "." + phase_.name_ + "." + algebraic_species_.name_);
        ++i_phase;
      }
      return indices;
    }

    /// @brief Build Jacobian sparse matrix indices for all phase instances
    JacobianIndices GetJacobianIndices(const StateVariableIndices& var_indices, const auto& jacobian) const
    {
      JacobianIndices jac_indices;
      jac_indices.indices_.resize(var_indices.number_of_phase_instances_);

      std::size_t num_blocks = jacobian.NumberOfBlocks();
      for (std::size_t i_phase = 0; i_phase < var_indices.number_of_phase_instances_; ++i_phase)
      {
        std::size_t alg_row = var_indices.algebraic_indices_[i_phase];
        auto& inst = jac_indices.indices_[i_phase];

        for (std::size_t block = 0; block < num_blocks; ++block)
        {
          // Reactant columns
          for (std::size_t r = 0; r < reactants_.size(); ++r)
            inst.push_back(jacobian.VectorIndex(block, alg_row, var_indices.reactant_indices_[i_phase][r]));

          // Product columns
          for (std::size_t p = 0; p < products_.size(); ++p)
            inst.push_back(jacobian.VectorIndex(block, alg_row, var_indices.product_indices_[i_phase][p]));

          // Solvent column
          inst.push_back(jacobian.VectorIndex(block, alg_row, var_indices.solvent_indices_[i_phase]));
        }
      }
      return jac_indices;
    }
  };
}  // namespace miam
