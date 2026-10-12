// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/types.hpp>

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
  /// @brief A linear algebraic constraint of the form  sum(coeff_i * [species_i]) = C
  /// @details Replaces the ODE row for a designated algebraic species with a linear
  ///          constraint equation:
  ///
  ///          \f$ G = \sum_i c_i \cdot [\text{species}_i] - C = 0 \f$
  ///
  ///          Typical uses include mass conservation and charge balance.
  ///
  ///          Resolution rules:
  ///          - If the algebraic species is in a non-instanced phase (not in phase_prefixes),
  ///            a single global constraint is generated. Instanced phase terms are summed
  ///            across all instances.
  ///          - If the algebraic species is in an instanced phase, one constraint is generated
  ///            per phase instance. Each constraint refers only to that instance's variables
  ///            (plus any non-instanced terms that are shared).
  ///
  ///          Limitation with multiple representations:
  ///          When diagnose_from_state_ is false, the same constant_ is broadcast to every
  ///          representation instance. A LinearConstraint with a fixed constant_ therefore
  ///          cannot encode different conserved totals across representations (e.g., two
  ///          droplets initialized to different amounts). In that case either set
  ///          diagnose_from_state_ = true so the solver diagnoses a per-instance constant
  ///          from the initial state, or omit the constraint if the system is already fully
  ///          determined by the remaining equilibrium and charge constraints.
  class LinearConstraint
  {
   public:
    /// @brief A term in the linear sum
    struct Term
    {
      micm::Phase phase;
      micm::Species species;
      double coefficient;
    };

    micm::Phase algebraic_phase_;        ///< Phase of the algebraic variable
    micm::Species algebraic_species_;    ///< Species whose ODE row is replaced
    std::vector<Term> terms_;            ///< Linear combination terms
    double constant_{ 0.0 };             ///< RHS constant C (used when diagnose_from_state_ is false)
    bool diagnose_from_state_{ false };  ///< If true, C is diagnosed from state at start of each Solve()
    std::string uuid_;                   ///< Unique identifier

    LinearConstraint() = delete;

    /// @brief Constructor
    LinearConstraint(
        const micm::Phase& algebraic_phase,
        const micm::Species& algebraic_species,
        const std::vector<Term>& terms,
        double constant,
        bool diagnose_from_state = false)
        : algebraic_phase_(algebraic_phase),
          algebraic_species_(algebraic_species),
          terms_(terms),
          constant_(constant),
          diagnose_from_state_(diagnose_from_state),
          uuid_(GenerateUuid())
    {
    }

    /// @brief Create a copy with a new UUID
    LinearConstraint CopyWithNewUuid() const
    {
      return LinearConstraint(algebraic_phase_, algebraic_species_, terms_, constant_, diagnose_from_state_);
    }

    /// @brief Returns the names of algebraic variables
    std::set<std::string> ConstraintAlgebraicVariableNames(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      auto phase_it = phase_prefixes.find(algebraic_phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        // Instanced phase: one algebraic variable per instance
        for (const auto& prefix : phase_it->second)
        {
          names.insert(prefix + "." + algebraic_phase_.name_ + "." + algebraic_species_.name_);
        }
      }
      else
      {
        // Non-instanced (gas) phase: single global algebraic variable
        names.insert(algebraic_species_.name_);
      }
      return names;
    }

    /// @brief Returns all species the constraint depends on
    std::set<std::string> ConstraintSpeciesDependencies(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> species_names;
      for (const auto& term : terms_)
      {
        auto phase_it = phase_prefixes.find(term.phase.name_);
        if (phase_it != phase_prefixes.end())
        {
          for (const auto& prefix : phase_it->second)
          {
            species_names.insert(prefix + "." + term.phase.name_ + "." + term.species.name_);
          }
        }
        else
        {
          species_names.insert(term.species.name_);
        }
      }
      return species_names;
    }

    /// @brief Returns non-zero constraint Jacobian element positions
    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      std::set<std::pair<std::size_t, std::size_t>> elements;
      bool is_global = (phase_prefixes.find(algebraic_phase_.name_) == phase_prefixes.end());

      if (is_global)
      {
        // Single algebraic row
        std::size_t alg_row = state_variable_indices.at(algebraic_species_.name_);
        for (const auto& term : terms_)
        {
          auto phase_it = phase_prefixes.find(term.phase.name_);
          if (phase_it != phase_prefixes.end())
          {
            for (const auto& prefix : phase_it->second)
            {
              std::size_t col = state_variable_indices.at(prefix + "." + term.phase.name_ + "." + term.species.name_);
              elements.insert({ alg_row, col });
            }
          }
          else
          {
            std::size_t col = state_variable_indices.at(term.species.name_);
            elements.insert({ alg_row, col });
          }
        }
      }
      else
      {
        // Per-instance algebraic rows
        const auto& alg_prefixes = phase_prefixes.at(algebraic_phase_.name_);
        for (const auto& prefix : alg_prefixes)
        {
          std::size_t alg_row =
              state_variable_indices.at(prefix + "." + algebraic_phase_.name_ + "." + algebraic_species_.name_);
          for (const auto& term : terms_)
          {
            auto phase_it = phase_prefixes.find(term.phase.name_);
            if (phase_it != phase_prefixes.end())
            {
              // Same instanced phase as algebraic: use only this instance
              std::size_t col = state_variable_indices.at(prefix + "." + term.phase.name_ + "." + term.species.name_);
              elements.insert({ alg_row, col });
            }
            else
            {
              // Non-instanced (gas) term: shared across all instances
              std::size_t col = state_variable_indices.at(term.species.name_);
              elements.insert({ alg_row, col });
            }
          }
        }
      }
      return elements;
    }

    /// @brief Returns a no-op constraint parameter update function (linear constraints have no parameters)
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateConstraintParametersFunction(
        const std::map<std::string, std::set<std::string>>& /*phase_prefixes*/,
        const std::unordered_map<std::string, std::size_t>& /*state_parameter_indices*/) const
    {
      return [](const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&) {};
    }

    /// @brief Returns state parameter names for constraints that update via conditions
    /// @details Diagnosed-from-state parameters are handled by InitializeConstraintParameterNames()
    ///          instead, so this always returns an empty set.
    std::set<std::string> ConstraintStateParameterNames(
        const std::map<std::string, std::set<std::string>>& /*phase_prefixes*/) const
    {
      return {};
    }

    /// @brief Returns parameter names that need initialization from state variables
    std::set<std::string> InitializeConstraintParameterNames(
        const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      if (!diagnose_from_state_)
        return {};
      return DiagnoseParamNames(phase_prefixes);
    }

    /// @brief Returns a function that diagnoses constraint constants from the current state
    /// @details Computes C = sum(c_i * [species_i]) for each grid cell at the start of each Solve()
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> InitializeConstraintParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      if (!diagnose_from_state_)
        return [](const DenseMatrixPolicy&, DenseMatrixPolicy&) {};

      FlatTerms flat = GetFlatTerms(phase_prefixes, state_parameter_indices, state_variable_indices);
      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      Vector param_indices(flat.param_indices);
      Vector counts(flat.counts);
      Vector term_indices(flat.term_indices);
      typename DenseMatrixPolicy::template VectorType<double> coefficients(flat.coefficients);
      param_indices.CopyToDevice();
      counts.CopyToDevice();
      term_indices.CopyToDevice();
      coefficients.CopyToDevice();
      std::size_t num_instances = flat.alg_indices.size();

      struct Storage
      {
        Vector param_indices, counts, term_indices;
        typename DenseMatrixPolicy::template VectorType<double> coefficients;
      };
      auto storage = std::make_shared<Storage>(
          Storage{ std::move(param_indices), std::move(counts), std::move(term_indices), std::move(coefficients) });
      auto param_view = storage->param_indices.GetView();
      auto count_view = storage->counts.GetView();
      auto term_view = storage->term_indices.GetView();
      auto coeff_view = storage->coefficients.GetView();
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ViewType& params_view)
      {
              std::size_t term_offset = 0;
              for (std::size_t i_inst = 0; i_inst < num_instances; ++i_inst)
            {
                auto total = params_view.GetRowVariable();
                params_view.ForEachRowStrict([](micm::Real& t) { t = 0.0; }, total);
                for (std::size_t k = 0; k < count_view[i_inst]; ++k)
                {
                  const micm::Real coeff = coeff_view[term_offset + k];
                  params_view.ForEachRowStrict(
                      [coeff](const micm::Real& val, micm::Real& t) { t += coeff * val; },
                      state_view.GetConstColumnView(term_view[term_offset + k]),
                    total);
      }
                params_view.ForEachRowStrict(
                    [](const micm::Real& t, micm::Real& param) { param = t; },
                    total,
                    params_view.GetColumnView(param_view[i_inst]));
                term_offset += count_view[i_inst];
              }
            },
            dummy_state_variables,
            dummy_state_parameters);

      return [storage, function](const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& state_parameters) mutable
      { function(state_variables, state_parameters); };
    }

    /// @brief Returns a function that computes constraint residuals G(y) = 0
    /// @details G = sum(coeff_i * [species_i]) - C
    ///          When diagnose_from_state_ is true, C is read from state_parameters per grid cell.
    ///          Otherwise, C is the compile-time constant_.
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      bool diagnose = diagnose_from_state_;
      double constant = constant_;

      FlatTerms flat = GetFlatTerms(phase_prefixes, state_parameter_indices, state_variable_indices);
      using Vector = typename DenseMatrixPolicy::template VectorType<std::size_t>;
      Vector alg_indices(flat.alg_indices);
      Vector param_indices(flat.param_indices);
      Vector counts(flat.counts);
      Vector term_indices(flat.term_indices);
      typename DenseMatrixPolicy::template VectorType<double> coefficients(flat.coefficients);
      alg_indices.CopyToDevice();
      param_indices.CopyToDevice();
      counts.CopyToDevice();
      term_indices.CopyToDevice();
      coefficients.CopyToDevice();
      std::size_t num_instances = flat.alg_indices.size();

      struct Storage
      {
        Vector alg_indices, param_indices, counts, term_indices;
        typename DenseMatrixPolicy::template VectorType<double> coefficients;
      };
      auto storage = std::make_shared<Storage>(Storage{ std::move(alg_indices),
                                                        std::move(param_indices),
                                                        std::move(counts),
                                                        std::move(term_indices),
                                                        std::move(coefficients) });
      auto alg_view = storage->alg_indices.GetView();
      auto param_view = storage->param_indices.GetView();
      auto count_view = storage->counts.GetView();
      auto term_view = storage->term_indices.GetView();
      auto coeff_view = storage->coefficients.GetView();
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };

        if (diagnose)
        {
          DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
        auto function = DenseMatrixPolicy::Function(
              MICM_LAMBDA(
                  const typename DenseMatrixPolicy::ConstViewType& state_view,
                  const typename DenseMatrixPolicy::ConstViewType& params_view,
                  const typename DenseMatrixPolicy::ViewType& residual_view)
              {
                std::size_t term_offset = 0;
                for (std::size_t i_inst = 0; i_inst < num_instances; ++i_inst)
                {
                  auto sum = residual_view.GetRowVariable();
                  residual_view.ForEachRowStrict(
                      [](const micm::Real& param, micm::Real& s) { s = -param; },
                      params_view.GetConstColumnView(param_view[i_inst]),
                      sum);
                  for (std::size_t k = 0; k < count_view[i_inst]; ++k)
                  {
                    const micm::Real coeff = coeff_view[term_offset + k];
                    residual_view.ForEachRowStrict(
                        [coeff](const micm::Real& val, micm::Real& s) { s += coeff * val; },
                        state_view.GetConstColumnView(term_view[term_offset + k]),
                        sum);
        }
                  residual_view.ForEachRowStrict(
                      [](const micm::Real& s, micm::Real& res) { res = s; },
                      sum,
                      residual_view.GetColumnView(alg_view[i_inst]));
                  term_offset += count_view[i_inst];
                }
              },
              dummy_state_variables,
              dummy_state_parameters,
              dummy_state_variables);

        return [storage, function](
                     const DenseMatrixPolicy& state_variables,
                     const DenseMatrixPolicy& state_parameters,
                   DenseMatrixPolicy& residual) mutable { function(state_variables, state_parameters, residual); };
        }

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ViewType& residual_view)
        {
              std::size_t term_offset = 0;
              for (std::size_t i_inst = 0; i_inst < num_instances; ++i_inst)
              {
                auto sum = residual_view.GetRowVariable();
                residual_view.ForEachRowStrict([constant](micm::Real& s) { s = -constant; }, sum);
                for (std::size_t k = 0; k < count_view[i_inst]; ++k)
                {
                  const micm::Real coeff = coeff_view[term_offset + k];
                  residual_view.ForEachRowStrict(
                      [coeff](const micm::Real& val, micm::Real& s) { s += coeff * val; },
                      state_view.GetConstColumnView(term_view[term_offset + k]),
                        sum);
        }
                residual_view.ForEachRowStrict(
                    [](const micm::Real& s, micm::Real& res) { res = s; }, sum, residual_view.GetColumnView(alg_view[i_inst]));
                term_offset += count_view[i_inst];
                }
              },
              dummy_state_variables,
              dummy_state_variables);

      return [storage, function](
                     const DenseMatrixPolicy& state_variables,
                     const DenseMatrixPolicy& /*state_parameters*/,
                 DenseMatrixPolicy& residual) mutable { function(state_variables, residual); };
    }

    /// @brief Returns a function that computes constraint Jacobian entries (subtracts dG/dy)
    /// @details dG/d[species_i] = coeff_i. Follows MICM convention: jac -= dG/dy
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      FlatTerms flat = GetFlatTerms(phase_prefixes, state_parameter_indices, state_variable_indices);

      // Pre-compute Jacobian VectorIndex offsets per instance (block 0)
      std::vector<std::size_t> jac_ids;
      std::size_t term_offset = 0;
      for (std::size_t i_inst = 0; i_inst < flat.alg_indices.size(); ++i_inst)
      {
        for (std::size_t k = 0; k < flat.counts[i_inst]; ++k)
          jac_ids.push_back(jacobian.VectorIndex(0, flat.alg_indices[i_inst], flat.term_indices[term_offset + k]));
        term_offset += flat.counts[i_inst];
      }

      typename SparseMatrixPolicy::template VectorType<std::size_t> jac_ids_vec(jac_ids);
      typename SparseMatrixPolicy::template VectorType<double> coefficients(flat.coefficients);
      jac_ids_vec.CopyToDevice();
      coefficients.CopyToDevice();
      std::size_t num_terms = jac_ids.size();

      struct Storage
      {
        typename SparseMatrixPolicy::template VectorType<std::size_t> jac_ids_vec;
        typename SparseMatrixPolicy::template VectorType<double> coefficients;
      };
      auto storage = std::make_shared<Storage>(Storage{ std::move(jac_ids_vec), std::move(coefficients) });
      auto jac_id_view = storage->jac_ids_vec.GetView();
      auto coeff_view = storage->coefficients.GetView();

      auto function = SparseMatrixPolicy::Function(
            MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jac_view)
            {
              for (std::size_t k = 0; k < num_terms; ++k)
              {
                const micm::Real coeff = coeff_view[k];
                jac_view.ForEachBlockStrict([coeff](micm::Real& j) { j -= coeff; }, jac_view.GetBlockView(jac_id_view[k]));
              }
            },
            jacobian);

      return [storage, function](
                 const DenseMatrixPolicy& /*state_variables*/,
                 const DenseMatrixPolicy& /*state_parameters*/,
                 SparseMatrixPolicy& jacobian_values) mutable { function(jacobian_values); };
    }

   private:
    /// @brief Flattened (index, coefficient) terms of all constraint instances
    /// @details A global constraint has one instance.
    struct FlatTerms
    {
      std::vector<std::size_t> alg_indices;    ///< Algebraic row per instance
      std::vector<std::size_t> param_indices;  ///< Diagnosed constant per instance (empty if not diagnosed)
      std::vector<std::size_t> counts;         ///< Number of terms per instance
      std::vector<std::size_t> term_indices;   ///< State variable index per term
      std::vector<double> coefficients;        ///< Coefficient per term
    };

    /// @brief Resolve the terms, algebraic rows, and diagnosed parameters of all instances
    FlatTerms GetFlatTerms(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      FlatTerms flat;
      std::vector<std::vector<std::pair<std::size_t, double>>> per_instance;
      bool is_global = (phase_prefixes.find(algebraic_phase_.name_) == phase_prefixes.end());
      if (is_global)
      {
        per_instance.push_back(ResolveGlobalTerms(phase_prefixes, state_variable_indices));
        flat.alg_indices.push_back(state_variable_indices.at(algebraic_species_.name_));
        if (diagnose_from_state_)
          flat.param_indices.push_back(state_parameter_indices.at("LC_" + uuid_ + "_constant"));
      }
      else
      {
        per_instance = ResolvePerInstanceTerms(phase_prefixes, state_variable_indices);
        for (const auto& prefix : phase_prefixes.at(algebraic_phase_.name_))
        {
          flat.alg_indices.push_back(
              state_variable_indices.at(prefix + "." + algebraic_phase_.name_ + "." + algebraic_species_.name_));
          if (diagnose_from_state_)
            flat.param_indices.push_back(state_parameter_indices.at("LC_" + uuid_ + "_" + prefix + "_constant"));
        }
      }
      for (const auto& terms : per_instance)
      {
        flat.counts.push_back(terms.size());
        for (const auto& [idx, coeff] : terms)
        {
          flat.term_indices.push_back(idx);
          flat.coefficients.push_back(coeff);
        }
      }
      return flat;
        }

    /// @brief Returns parameter names for diagnosed constants
    std::set<std::string> DiagnoseParamNames(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> names;
      bool is_global = (phase_prefixes.find(algebraic_phase_.name_) == phase_prefixes.end());
      if (is_global)
      {
        names.insert("LC_" + uuid_ + "_constant");
      }
      else
      {
        for (const auto& prefix : phase_prefixes.at(algebraic_phase_.name_))
          names.insert("LC_" + uuid_ + "_" + prefix + "_constant");
      }
      return names;
    }

    /// @brief Resolve terms for a global (non-instanced algebraic) constraint.
    ///        All instanced terms are expanded into (index, coefficient) pairs.
    std::vector<std::pair<std::size_t, double>> ResolveGlobalTerms(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      std::vector<std::pair<std::size_t, double>> resolved;
      for (const auto& term : terms_)
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

    /// @brief Resolve terms for per-instance constraints.
    ///        Returns a vector of (index, coefficient) pairs per algebraic instance.
    std::vector<std::vector<std::pair<std::size_t, double>>> ResolvePerInstanceTerms(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      const auto& alg_prefixes = phase_prefixes.at(algebraic_phase_.name_);
      std::vector<std::vector<std::pair<std::size_t, double>>> per_instance(alg_prefixes.size());

      std::size_t i_inst = 0;
      for (const auto& alg_prefix : alg_prefixes)
      {
        for (const auto& term : terms_)
        {
          auto phase_it = phase_prefixes.find(term.phase.name_);
          if (phase_it != phase_prefixes.end())
          {
            // Match the instance prefix for this term's phase
            std::size_t idx = state_variable_indices.at(alg_prefix + "." + term.phase.name_ + "." + term.species.name_);
            per_instance[i_inst].push_back({ idx, term.coefficient });
          }
          else
          {
            // Non-instanced term: shared gas-phase variable
            std::size_t idx = state_variable_indices.at(term.species.name_);
            per_instance[i_inst].push_back({ idx, term.coefficient });
          }
        }
        ++i_inst;
      }
      return per_instance;
    }
  };
}  // namespace miam