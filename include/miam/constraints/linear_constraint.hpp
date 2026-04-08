// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>

#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace miam
{
  namespace constraint
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

      micm::Phase algebraic_phase_;      ///< Phase of the algebraic variable
      micm::Species algebraic_species_;   ///< Species whose ODE row is replaced
      std::vector<Term> terms_;           ///< Linear combination terms
      double constant_{ 0.0 };            ///< RHS constant C
      std::string uuid_;                  ///< Unique identifier

      LinearConstraint() = delete;

      /// @brief Constructor
      LinearConstraint(
          const micm::Phase& algebraic_phase,
          const micm::Species& algebraic_species,
          const std::vector<Term>& terms,
          double constant)
          : algebraic_phase_(algebraic_phase),
            algebraic_species_(algebraic_species),
            terms_(terms),
            constant_(constant),
            uuid_(miam::util::generate_uuid_v4())
      {
      }

      /// @brief Create a copy with a new UUID
      LinearConstraint CopyWithNewUuid() const
      {
        return LinearConstraint(algebraic_phase_, algebraic_species_, terms_, constant_);
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
                std::size_t col =
                    state_variable_indices.at(prefix + "." + term.phase.name_ + "." + term.species.name_);
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
                std::size_t col =
                    state_variable_indices.at(prefix + "." + term.phase.name_ + "." + term.species.name_);
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
      std::function<void(const std::vector<micm::Conditions>&)> UpdateConstraintParametersFunction() const
      {
        return [](const std::vector<micm::Conditions>&) {};
      }

      /// @brief Returns a function that computes constraint residuals G(y) = 0
      /// @details G = sum(coeff_i * [species_i]) - C
      template<typename DenseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
      {
        bool is_global = (phase_prefixes.find(algebraic_phase_.name_) == phase_prefixes.end());
        double constant = constant_;

        DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };

        if (is_global)
        {
          auto resolved = ResolveGlobalTerms(phase_prefixes, state_variable_indices);
          std::size_t alg_idx = state_variable_indices.at(algebraic_species_.name_);

          auto inner = DenseMatrixPolicy::Function(
              [resolved, alg_idx, constant](auto&& state_variables, auto&& residual)
              {
                auto sum = residual.GetRowVariable();
                residual.ForEachRow([constant](double& s) { s = -constant; }, sum);
                for (const auto& [idx, coeff] : resolved)
                {
                  residual.ForEachRow(
                      [coeff](const double& val, double& s) { s += coeff * val; },
                      state_variables.GetConstColumnView(idx),
                      sum);
                }
                residual.ForEachRow(
                    [](const double& s, double& res) { res = s; },
                    sum,
                    residual.GetColumnView(alg_idx));
              },
              dummy_state, dummy_state);

          return [inner = std::move(inner)](const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& residual) mutable
          { inner(state_variables, residual); };
        }
        else
        {
          auto per_instance = ResolvePerInstanceTerms(phase_prefixes, state_variable_indices);
          std::vector<std::size_t> alg_indices;
          for (const auto& prefix : phase_prefixes.at(algebraic_phase_.name_))
          {
            alg_indices.push_back(
                state_variable_indices.at(prefix + "." + algebraic_phase_.name_ + "." + algebraic_species_.name_));
          }

          auto inner = DenseMatrixPolicy::Function(
              [per_instance, alg_indices, constant](auto&& state_variables, auto&& residual)
              {
                for (std::size_t i_inst = 0; i_inst < alg_indices.size(); ++i_inst)
                {
                  auto sum = residual.GetRowVariable();
                  residual.ForEachRow([constant](double& s) { s = -constant; }, sum);
                  for (const auto& [idx, coeff] : per_instance[i_inst])
                  {
                    residual.ForEachRow(
                        [coeff](const double& val, double& s) { s += coeff * val; },
                        state_variables.GetConstColumnView(idx),
                        sum);
                  }
                  residual.ForEachRow(
                      [](const double& s, double& res) { res = s; },
                      sum,
                      residual.GetColumnView(alg_indices[i_inst]));
                }
              },
              dummy_state, dummy_state);

          return [inner = std::move(inner)](const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& residual) mutable
          { inner(state_variables, residual); };
        }
      }

      /// @brief Returns a function that computes constraint Jacobian entries (subtracts dG/dy)
      /// @details dG/d[species_i] = coeff_i. Follows MICM convention: jac -= dG/dy
      template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_variable_indices,
          const SparseMatrixPolicy& jacobian) const
      {
        bool is_global = (phase_prefixes.find(algebraic_phase_.name_) == phase_prefixes.end());

        DenseMatrixPolicy dummy_state{ 1, state_variable_indices.size(), 0.0 };

        if (is_global)
        {
          auto resolved = ResolveGlobalTerms(phase_prefixes, state_variable_indices);
          std::size_t alg_row = state_variable_indices.at(algebraic_species_.name_);

          // Pre-compute Jacobian VectorIndex offsets (block 0)
          std::vector<std::pair<std::size_t, double>> jac_entries;
          for (const auto& [col_idx, coeff] : resolved)
          {
            jac_entries.push_back({ jacobian.VectorIndex(0, alg_row, col_idx), coeff });
          }

          auto inner = SparseMatrixPolicy::Function(
              [jac_entries](auto&& state_variables, auto&& jacobian_values)
              {
                for (const auto& [vec_idx, coeff] : jac_entries)
                {
                  auto bv = jacobian_values.GetBlockView(vec_idx);
                  jacobian_values.ForEachBlock(
                      [coeff](double& j) { j -= coeff; },
                      bv);
                }
              },
              dummy_state, jacobian);

          return [inner = std::move(inner)](const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian_values) mutable
          { inner(state_variables, jacobian_values); };
        }
        else
        {
          auto per_instance = ResolvePerInstanceTerms(phase_prefixes, state_variable_indices);
          std::vector<std::size_t> alg_rows;
          for (const auto& prefix : phase_prefixes.at(algebraic_phase_.name_))
          {
            alg_rows.push_back(
                state_variable_indices.at(prefix + "." + algebraic_phase_.name_ + "." + algebraic_species_.name_));
          }

          // Pre-compute Jacobian VectorIndex offsets per instance (block 0)
          std::vector<std::vector<std::pair<std::size_t, double>>> jac_entries_per_instance;
          for (std::size_t i_inst = 0; i_inst < alg_rows.size(); ++i_inst)
          {
            std::vector<std::pair<std::size_t, double>> inst_entries;
            for (const auto& [col_idx, coeff] : per_instance[i_inst])
            {
              inst_entries.push_back({ jacobian.VectorIndex(0, alg_rows[i_inst], col_idx), coeff });
            }
            jac_entries_per_instance.push_back(std::move(inst_entries));
          }

          auto inner = SparseMatrixPolicy::Function(
              [jac_entries_per_instance](auto&& state_variables, auto&& jacobian_values)
              {
                for (const auto& inst_entries : jac_entries_per_instance)
                {
                  for (const auto& [vec_idx, coeff] : inst_entries)
                  {
                    auto bv = jacobian_values.GetBlockView(vec_idx);
                    jacobian_values.ForEachBlock(
                        [coeff](double& j) { j -= coeff; },
                        bv);
                  }
                }
              },
              dummy_state, jacobian);

          return [inner = std::move(inner)](const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian_values) mutable
          { inner(state_variables, jacobian_values); };
        }
      }

     private:
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
              std::size_t idx =
                  state_variable_indices.at(prefix + "." + term.phase.name_ + "." + term.species.name_);
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
              std::size_t idx =
                  state_variable_indices.at(alg_prefix + "." + term.phase.name_ + "." + term.species.name_);
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

    /// @brief Builder for LinearConstraint
    class LinearConstraintBuilder
    {
     public:
      LinearConstraintBuilder() = default;

      LinearConstraintBuilder& SetAlgebraicSpecies(const micm::Phase& phase, const micm::Species& species)
      {
        algebraic_phase_ = phase;
        algebraic_species_ = species;
        algebraic_is_set_ = true;
        return *this;
      }

      LinearConstraintBuilder& AddTerm(const micm::Phase& phase, const micm::Species& species, double coefficient)
      {
        terms_.push_back({ phase, species, coefficient });
        return *this;
      }

      LinearConstraintBuilder& SetConstant(double constant)
      {
        constant_ = constant;
        return *this;
      }

      LinearConstraint Build() const
      {
        if (!algebraic_is_set_)
          throw std::runtime_error("LinearConstraintBuilder requires the algebraic species to be set.");
        if (terms_.empty())
          throw std::runtime_error("LinearConstraintBuilder requires at least one term.");

        return LinearConstraint(algebraic_phase_, algebraic_species_, terms_, constant_);
      }

     private:
      micm::Phase algebraic_phase_;
      micm::Species algebraic_species_;
      bool algebraic_is_set_ = false;
      std::vector<LinearConstraint::Term> terms_;
      double constant_{ 0.0 };
    };
  }  // namespace constraint
}  // namespace miam
