// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/constraints/dissolved_equilibrium_constraint_set.hpp>
#include <miam/constraints/henrys_law_equilibrium_constraint.hpp>
#include <miam/constraints/henrys_law_equilibrium_constraint_set.hpp>
#include <miam/constraints/linear_constraint.hpp>
#include <miam/constraints/linear_constraint_set.hpp>
#include <miam/processes.hpp>
#include <miam/processes/dissolved_reaction_set.hpp>
#include <miam/processes/dissolved_reversible_reaction_set.hpp>
#include <miam/processes/henrys_law_phase_transfer_set.hpp>
#include <miam/representations.hpp>
#include <miam/util/error.hpp>
#include <miam/util/matching_sparse.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/system/conditions.hpp>
#include <micm/util/matrix.hpp>

#include <algorithm>
#include <any>
#include <concepts>
#include <cstddef>
#include <functional>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <variant>
#include <vector>

namespace miam
{
  /// @brief Aerosol/Cloud Model
  /// @details Model is a collection of representations that collectively define an aerosol
  ///          and/or cloud system. Model is compatible with the micm::ExternalModelSystem
  ///          and micm::ExternalModelProcessSet interfaces.
  class Model
  {
   public:
    using RepresentationVariant = std::variant<SingleMomentMode, TwoMomentMode, UniformSection>;
    using ProcessVariant = std::variant<DissolvedReaction, DissolvedReversibleReaction, HenrysLawPhaseTransfer>;
    using ConstraintVariant = std::variant<DissolvedEquilibriumConstraint, HenrysLawEquilibriumConstraint, LinearConstraint>;

    std::string name_;
    std::vector<RepresentationVariant> representations_;
    std::vector<ProcessVariant> processes_{};
    std::vector<ConstraintVariant> constraints_{};

    std::unordered_map<std::string, std::size_t> state_parameter_indices_{};
    std::unordered_map<std::string, std::size_t> state_variable_indices_{};

    // Lazy caches for per-process / per-constraint std::function objects.
    mutable std::any cached_process_update_fns_{};
    mutable std::any cached_constraint_update_fns_{};
    mutable std::any cached_constraint_init_fns_{};

    // Cached DP-typed aerosol descriptor map used by cached AddForcingTerms /
    // SubtractJacobianTerms. Rebuilt lazily when the stored DP does not match.
    mutable std::any cached_descriptors_{};

    // Solve-time companion Sets. Populated in FinalizeProcessSetup (which is when
    // the sparse Jacobian pattern is available); consumed on-device by AddForcingTerms /
    // SubtractJacobianTerms.
    mutable std::any dissolved_reaction_sets_any_{};
    mutable std::any dissolved_reversible_reaction_sets_any_{};
    mutable std::any henrys_law_phase_transfer_sets_any_{};

    // Solve-time companion Sets for constraints. Populated in FinalizeConstraintSetup.
    // Each holds `std::vector<ConstraintSet<DP, SP>>` for the SP MICM finalized us with
    // (DP recovered via `detail::MatchingDenseT<SP>`).
    mutable std::any linear_constraint_sets_any_{};
    mutable std::any dissolved_equilibrium_constraint_sets_any_{};
    mutable std::any henrys_law_equilibrium_constraint_sets_any_{};

    /// @brief Returns the total state size (number of variables, number of parameters)
    std::tuple<std::size_t, std::size_t> StateSize() const
    {
      std::size_t num_variables = 0;
      std::size_t num_parameters = 0;
      for (const auto& repr : representations_)
      {
        std::visit(
            [&](const auto& r)
            {
              auto [vars, params] = r.StateSize();
              num_variables += vars;
              num_parameters += params;
            },
            repr);
      }
      // Add parameters for each process
      auto phase_prefixes = CollectPhaseStatePrefixes();
      ForEachProcess(
          [&](const auto& process)
          {
            auto process_params = process.ProcessParameterNames(phase_prefixes);
            num_parameters += process_params.size();
          });
      return { num_variables, num_parameters };
    }

    /// @brief Returns unique names for all state variables
    std::set<std::string> StateVariableNames() const
    {
      std::set<std::string> names;
      for (const auto& repr : representations_)
      {
        std::visit(
            [&](const auto& r)
            {
              auto repr_names = r.StateVariableNames();
              names.insert(repr_names.begin(), repr_names.end());
            },
            repr);
      }
      return names;
    }

    /// @brief Returns unique names for all state parameters
    std::set<std::string> StateParameterNames() const
    {
      std::set<std::string> names;
      // Collect parameter names from all representations
      for (const auto& repr : representations_)
      {
        std::visit(
            [&](const auto& r)
            {
              auto repr_names = r.StateParameterNames();
              names.insert(repr_names.begin(), repr_names.end());
            },
            repr);
      }
      // Add parameters for each process
      auto phase_prefixes = CollectPhaseStatePrefixes();
      ForEachProcess(
          [&](const auto& process)
          {
            auto process_params = process.ProcessParameterNames(phase_prefixes);
            names.insert(process_params.begin(), process_params.end());
          });
      return names;
    }

    /// @brief Returns names of all species used in the model's processes and constraints
    std::set<std::string> SpeciesUsed() const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      // Collect participating species' unique state names for all processes
      std::set<std::string> species_names;
      ForEachProcess(
          [&](const auto& process)
          {
            auto process_species = process.SpeciesUsed(phase_prefixes);
            species_names.insert(process_species.begin(), process_species.end());
          });
      // Collect species from constraints
      ForEachConstraint(
          [&](const auto& c)
          {
            auto deps = c.ConstraintSpeciesDependencies(phase_prefixes);
            species_names.insert(deps.begin(), deps.end());
          });
      return species_names;
    }

    /// @brief Add processes to the model
    /// @details Accepts a vector of any process type stored in ProcessVariant.
    ///          Each process is copied with a new UUID to ensure uniqueness across models.
    template<typename ProcessType>
    void AddProcesses(const std::vector<ProcessType>& new_processes)
    {
      for (const auto& process : new_processes)
      {
        processes_.push_back(ProcessVariant{ process.CopyWithNewUuid() });
      }
    }

    /// @brief Add processes to the model from an initializer list
    template<typename ProcessType>
    void AddProcesses(std::initializer_list<ProcessType> new_processes)
    {
      for (const auto& process : new_processes)
      {
        processes_.push_back(ProcessVariant{ process.CopyWithNewUuid() });
      }
    }

    /// @brief Add processes to the model (variadic form for mixed types)
    template<typename... ProcessTypes>
      requires(sizeof...(ProcessTypes) >= 1 && (std::constructible_from<ProcessVariant, std::decay_t<ProcessTypes>> && ...))
    void AddProcesses(ProcessTypes&&... processes)
    {
      (processes_.push_back(ProcessVariant{ processes.CopyWithNewUuid() }), ...);
    }

    /// @brief Add constraints to the model
    template<typename ConstraintType>
    void AddConstraints(const std::vector<ConstraintType>& new_constraints)
    {
      for (const auto& c : new_constraints)
      {
        constraints_.push_back(ConstraintVariant{ c.CopyWithNewUuid() });
      }
    }

    /// @brief Add constraints to the model from an initializer list
    template<typename ConstraintType>
    void AddConstraints(std::initializer_list<ConstraintType> new_constraints)
    {
      for (const auto& c : new_constraints)
      {
        constraints_.push_back(ConstraintVariant{ c.CopyWithNewUuid() });
      }
    }

    /// @brief Add constraints to the model (variadic form for mixed types)
    template<typename... ConstraintTypes>
      requires(
          sizeof...(ConstraintTypes) >= 1 &&
          (std::constructible_from<ConstraintVariant, std::decay_t<ConstraintTypes>> && ...))
    void AddConstraints(ConstraintTypes&&... constraints)
    {
      (constraints_.push_back(ConstraintVariant{ constraints.CopyWithNewUuid() }), ...);
    }

    /// @brief Returns non-zero Jacobian element positions
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& state_indices) const
    {
      // Collect needed Jacobian element indices from all processes
      std::set<std::pair<std::size_t, std::size_t>> elements;
      auto phase_prefixes = CollectPhaseStatePrefixes();
      ForEachProcess(
          [&](const auto& process)
          {
            auto process_elements = process.NonZeroJacobianElements(phase_prefixes, state_indices);
            elements.insert(process_elements.begin(), process_elements.end());
          });
      return elements;
    }

    /// @brief Returns a function that updates state parameters
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateStateParametersFunction(const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
    {
      // Collect parameter update functions from all processes and return a combined function
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::vector<
          std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>>
          update_functions;
      ForEachProcess(
          [&](const auto& process)
          {
            auto update_fn =
                process.template UpdateStateParametersFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices);
            update_functions.push_back(update_fn);
          });
      return [update_functions](
                 const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
                 DenseMatrixPolicy& state_parameters)
      {
        for (const auto& fn : update_functions)
        {
          fn(conditions, state_parameters);
        }
      };
    }

    /// @brief Returns a function that calculates forcing terms.
    /// @details `SparseMatrixPolicy` is required so that process Sets can pre-compute Jacobian flat IDs at build time.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      auto descriptors =
          BuildDescriptors<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
      auto nz_elements = NonZeroJacobianElements(state_variable_indices);
      auto jacobian_builder = SparseMatrixPolicy::Create(state_variable_indices.size()).InitialValue(0.0);
      for (const auto& elem : nz_elements)
        jacobian_builder = jacobian_builder.WithElement(elem.first, elem.second);
      SparseMatrixPolicy jacobian_pattern(jacobian_builder);

      std::vector<DissolvedReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>> dr_sets;
      std::vector<DissolvedReversibleReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>> drr_sets;
      std::vector<HenrysLawPhaseTransferSet<DenseMatrixPolicy, SparseMatrixPolicy>> hlpt_sets;
      ForEachProcess(
          [&](const auto& process)
          {
            using P = std::decay_t<decltype(process)>;
            if constexpr (std::is_same_v<P, DissolvedReaction>)
              dr_sets.emplace_back(process, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian_pattern);
            else if constexpr (std::is_same_v<P, DissolvedReversibleReaction>)
              drr_sets.emplace_back(process, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian_pattern);
            else if constexpr (std::is_same_v<P, HenrysLawPhaseTransfer>)
              hlpt_sets.emplace_back(
                  process, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian_pattern, descriptors);
          });
      return [dr_sets = std::move(dr_sets),
              drr_sets = std::move(drr_sets),
              hlpt_sets = std::move(hlpt_sets),
              descriptors = std::move(descriptors)](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms)
      {
        for (const auto& set : dr_sets)
          set.AddForcingTerms(state_parameters, state_variables, forcing_terms);
        for (const auto& set : drr_sets)
          set.AddForcingTerms(state_parameters, state_variables, forcing_terms);
        for (const auto& set : hlpt_sets)
          set.AddForcingTerms(state_parameters, state_variables, forcing_terms, descriptors);
      };
    }

    /// @brief Returns a function that calculates Jacobian contributions
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      auto descriptors =
          BuildDescriptors<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
      std::vector<DissolvedReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>> dr_sets;
      std::vector<DissolvedReversibleReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>> drr_sets;
      std::vector<HenrysLawPhaseTransferSet<DenseMatrixPolicy, SparseMatrixPolicy>> hlpt_sets;
      ForEachProcess(
          [&](const auto& process)
          {
            using P = std::decay_t<decltype(process)>;
            if constexpr (std::is_same_v<P, DissolvedReaction>)
              dr_sets.emplace_back(process, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
            else if constexpr (std::is_same_v<P, DissolvedReversibleReaction>)
              drr_sets.emplace_back(process, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
            else if constexpr (std::is_same_v<P, HenrysLawPhaseTransfer>)
              hlpt_sets.emplace_back(
                  process, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, descriptors);
          });
      return [dr_sets = std::move(dr_sets),
              drr_sets = std::move(drr_sets),
              hlpt_sets = std::move(hlpt_sets),
              descriptors = std::move(descriptors)](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian)
      {
        for (const auto& set : dr_sets)
          set.SubtractJacobianTerms(state_parameters, state_variables, jacobian);
        for (const auto& set : drr_sets)
          set.SubtractJacobianTerms(state_parameters, state_variables, jacobian);
        for (const auto& set : hlpt_sets)
          set.SubtractJacobianTerms(state_parameters, state_variables, jacobian, descriptors);
      };
    }

    // ── HasConstraints concept methods ──

    /// @brief Returns unique names for constraint-specific state parameters
    ///        Includes parameters for diagnosed constants (mass conservation totals)
    std::set<std::string> ConstraintStateParameterNames() const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::set<std::string> names;
      ForEachConstraint(
          [&](const auto& c)
          {
            if constexpr (requires { c.ConstraintStateParameterNames(phase_prefixes); })
            {
              auto c_names = c.ConstraintStateParameterNames(phase_prefixes);
              names.insert(c_names.begin(), c_names.end());
            }
          });
      return names;
    }

    // ── HasInitializeConstraintParameters concept methods ──

    /// @brief Returns parameter names that need initialization from state variables
    std::set<std::string> InitializeConstraintParameterNames() const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::set<std::string> names;
      ForEachConstraint(
          [&](const auto& c)
          {
            if constexpr (requires { c.InitializeConstraintParameterNames(phase_prefixes); })
            {
              auto c_names = c.InitializeConstraintParameterNames(phase_prefixes);
              names.insert(c_names.begin(), c_names.end());
            }
          });
      return names;
    }

    /// @brief Returns a function that diagnoses constraint parameters from current state
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> InitializeConstraintParametersFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::vector<std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>> init_fns;
      ForEachConstraint(
          [&](const auto& c)
          {
            if constexpr (requires {
                            c.template InitializeConstraintParametersFunction<DenseMatrixPolicy>(
                                phase_prefixes, state_parameter_indices, state_variable_indices);
                          })
            {
              init_fns.push_back(c.template InitializeConstraintParametersFunction<DenseMatrixPolicy>(
                  phase_prefixes, state_parameter_indices, state_variable_indices));
            }
          });
      return [init_fns](const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& state_parameters)
      {
        for (const auto& fn : init_fns)
          fn(state_variables, state_parameters);
      };
    }

    /// @brief Returns a function that updates constraint parameters based on conditions
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::vector<
          std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>>
          update_fns;
      ForEachConstraint(
          [&](const auto& c)
          {
            update_fns.push_back(
                c.template UpdateConstraintParametersFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices));
          });
      return [update_fns](
                 const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
                 DenseMatrixPolicy& state_parameters) mutable
      {
        for (auto& fn : update_fns)
          fn(conditions, state_parameters);
      };
    }

    /// @brief Returns names of all algebraic variables across all constraints
    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::set<std::string> names;
      ForEachConstraint(
          [&](const auto& c)
          {
            auto c_names = c.ConstraintAlgebraicVariableNames(phase_prefixes);
            names.insert(c_names.begin(), c_names.end());
          });
      return names;
    }

    /// @brief Returns all species that constraints depend on
    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::set<std::string> species_names;
      ForEachConstraint(
          [&](const auto& c)
          {
            auto deps = c.ConstraintSpeciesDependencies(phase_prefixes);
            species_names.insert(deps.begin(), deps.end());
          });
      return species_names;
    }

    /// @brief Returns non-zero constraint Jacobian element positions
    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& state_indices) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::set<std::pair<std::size_t, std::size_t>> elements;
      ForEachConstraint(
          [&](const auto& c)
          {
            auto c_elements = c.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
            elements.insert(c_elements.begin(), c_elements.end());
          });
      return elements;
    }

    /// @brief Cache build-time indices for solve-time direct-dispatch methods
    /// @details Called once by micm::SolverBuilder after the parameter map,
    ///          species map, and Jacobian sparsity pattern are finalized.
    template<typename SparseMatrixPolicy>
    void FinalizeProcessSetup(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      state_parameter_indices_ = state_parameter_indices;
      state_variable_indices_ = state_variable_indices;
      cached_process_update_fns_.reset();
      cached_descriptors_.reset();

      auto phase_prefixes = CollectPhaseStatePrefixes();
      // Reference (CPU) descriptor map used only to enumerate `n_deps` counts for
      // Jacobian flat-ID layout in HenrysLawPhaseTransferSet.
      using ReferenceDense = micm::Matrix<double>;
      auto reference_descriptors =
          BuildDescriptors<ReferenceDense>(phase_prefixes, state_parameter_indices, state_variable_indices);

      using MatchingDP = detail::MatchingDenseT<SparseMatrixPolicy>;
      using DissolvedReactionSetT = DissolvedReactionSet<MatchingDP, SparseMatrixPolicy>;
      using DissolvedReversibleReactionSetT = DissolvedReversibleReactionSet<MatchingDP, SparseMatrixPolicy>;
      using HenrysLawPhaseTransferSetT = HenrysLawPhaseTransferSet<MatchingDP, SparseMatrixPolicy>;
      std::vector<DissolvedReactionSetT> dissolved_reaction_sets;
      std::vector<DissolvedReversibleReactionSetT> dissolved_reversible_reaction_sets;
      std::vector<HenrysLawPhaseTransferSetT> henrys_law_phase_transfer_sets;
      for (const auto& process : processes_)
      {
        std::visit(
            [&](const auto& p)
            {
              using P = std::decay_t<decltype(p)>;
              if constexpr (std::is_same_v<P, DissolvedReaction>)
              {
                dissolved_reaction_sets.emplace_back(
                    p, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              }
              else if constexpr (std::is_same_v<P, DissolvedReversibleReaction>)
              {
                dissolved_reversible_reaction_sets.emplace_back(
                    p, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              }
              else if constexpr (std::is_same_v<P, HenrysLawPhaseTransfer>)
              {
                henrys_law_phase_transfer_sets.emplace_back(
                    p, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, reference_descriptors);
              }
            },
            process);
      }
      dissolved_reaction_sets_any_ = std::move(dissolved_reaction_sets);
      dissolved_reversible_reaction_sets_any_ = std::move(dissolved_reversible_reaction_sets);
      henrys_law_phase_transfer_sets_any_ = std::move(henrys_law_phase_transfer_sets);
    }

    /// @brief Cache build-time indices for constraint solve-time methods
    template<typename SparseMatrixPolicy>
    void FinalizeConstraintSetup(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      state_parameter_indices_ = state_parameter_indices;
      state_variable_indices_ = state_variable_indices;
      cached_constraint_update_fns_.reset();
      cached_constraint_init_fns_.reset();

      auto phase_prefixes = CollectPhaseStatePrefixes();

      using MatchingDP = detail::MatchingDenseT<SparseMatrixPolicy>;
      using LinearSetT = LinearConstraintSet<MatchingDP, SparseMatrixPolicy>;
      using DissEqSetT = DissolvedEquilibriumConstraintSet<MatchingDP, SparseMatrixPolicy>;
      using HenrysEqSetT = HenrysLawEquilibriumConstraintSet<MatchingDP, SparseMatrixPolicy>;
      std::vector<LinearSetT> linear_sets;
      std::vector<DissEqSetT> dissolved_equilibrium_sets;
      std::vector<HenrysEqSetT> henrys_law_equilibrium_sets;
      for (const auto& constraint : constraints_)
      {
        std::visit(
            [&](const auto& c)
            {
              using C = std::decay_t<decltype(c)>;
              if constexpr (std::is_same_v<C, LinearConstraint>)
                linear_sets.emplace_back(
                    c, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              else if constexpr (std::is_same_v<C, DissolvedEquilibriumConstraint>)
                dissolved_equilibrium_sets.emplace_back(
                    c, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              else if constexpr (std::is_same_v<C, HenrysLawEquilibriumConstraint>)
                henrys_law_equilibrium_sets.emplace_back(
                    c, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
            },
            constraint);
      }
      linear_constraint_sets_any_ = std::move(linear_sets);
      dissolved_equilibrium_constraint_sets_any_ = std::move(dissolved_equilibrium_sets);
      henrys_law_equilibrium_constraint_sets_any_ = std::move(henrys_law_equilibrium_sets);
    }

    /// @brief Solve-time: refresh temperature-/pressure-dependent process parameters
    template<typename DenseMatrixPolicy>
    void UpdateStateParameters(
        const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
        DenseMatrixPolicy& state_parameters) const
    {
      using FnType = std::function<void(
          const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>;
      using CacheType = std::vector<FnType>;
      if (!cached_process_update_fns_.has_value())
      {
        CacheType cache;
        auto phase_prefixes = CollectPhaseStatePrefixes();
        ForEachProcess(
            [&](const auto& process)
            {
              cache.push_back(process.template UpdateStateParametersFunction<DenseMatrixPolicy>(
                  phase_prefixes, state_parameter_indices_));
            });
        cached_process_update_fns_ = std::move(cache);
      }
      for (auto& fn : std::any_cast<CacheType&>(cached_process_update_fns_))
        fn(conditions, state_parameters);
    }

    /// @brief Solve-time: add process forcing (tendency) contributions
    template<typename DenseMatrixPolicy>
    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const
    {
      using SP = detail::MatchingSparseT<DenseMatrixPolicy>;
      using DissolvedReactionSetT = DissolvedReactionSet<DenseMatrixPolicy, SP>;
      using DissolvedReversibleReactionSetT = DissolvedReversibleReactionSet<DenseMatrixPolicy, SP>;
      using HenrysLawPhaseTransferSetT = HenrysLawPhaseTransferSet<DenseMatrixPolicy, SP>;
      if (dissolved_reaction_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<DissolvedReactionSetT>&>(dissolved_reaction_sets_any_);
        for (const auto& set : sets)
          set.AddForcingTerms(state_parameters, state_variables, forcing);
      }
      if (dissolved_reversible_reaction_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<DissolvedReversibleReactionSetT>&>(dissolved_reversible_reaction_sets_any_);
        for (const auto& set : sets)
          set.AddForcingTerms(state_parameters, state_variables, forcing);
      }
      if (henrys_law_phase_transfer_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<HenrysLawPhaseTransferSetT>&>(henrys_law_phase_transfer_sets_any_);
        if (!sets.empty())
        {
          const auto& descriptors = GetCachedDescriptors<DenseMatrixPolicy>();
          for (const auto& set : sets)
            set.AddForcingTerms(state_parameters, state_variables, forcing, descriptors);
        }
      }
    }

    /// @brief Solve-time: subtract process Jacobian contributions (-J convention)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      using DissolvedReactionSetT = DissolvedReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>;
      using DissolvedReversibleReactionSetT = DissolvedReversibleReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>;
      using HenrysLawPhaseTransferSetT = HenrysLawPhaseTransferSet<DenseMatrixPolicy, SparseMatrixPolicy>;
      if (dissolved_reaction_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<DissolvedReactionSetT>&>(dissolved_reaction_sets_any_);
        for (const auto& set : sets)
          set.SubtractJacobianTerms(state_parameters, state_variables, jacobian);
      }
      if (dissolved_reversible_reaction_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<DissolvedReversibleReactionSetT>&>(dissolved_reversible_reaction_sets_any_);
        for (const auto& set : sets)
          set.SubtractJacobianTerms(state_parameters, state_variables, jacobian);
      }
      if (henrys_law_phase_transfer_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<HenrysLawPhaseTransferSetT>&>(henrys_law_phase_transfer_sets_any_);
        if (!sets.empty())
        {
          const auto& descriptors = GetCachedDescriptors<DenseMatrixPolicy>();
          for (const auto& set : sets)
            set.SubtractJacobianTerms(state_parameters, state_variables, jacobian, descriptors);
        }
      }
    }

    /// @brief Solve-time: refresh temperature-/pressure-dependent constraint parameters
    template<typename DenseMatrixPolicy>
    void UpdateConstraintStateParameters(
        const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
        DenseMatrixPolicy& state_parameters) const
    {
      using FnType = std::function<void(
          const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>;
      using CacheType = std::vector<FnType>;
      if (!cached_constraint_update_fns_.has_value())
      {
        CacheType cache;
        auto phase_prefixes = CollectPhaseStatePrefixes();
        ForEachConstraint(
            [&](const auto& c)
            {
              cache.push_back(c.template UpdateConstraintParametersFunction<DenseMatrixPolicy>(
                  phase_prefixes, state_parameter_indices_));
            });
        cached_constraint_update_fns_ = std::move(cache);
      }
      for (auto& fn : std::any_cast<CacheType&>(cached_constraint_update_fns_))
        fn(conditions, state_parameters);
    }

    /// @brief Solve-time: diagnose constraint parameters from current state
    template<typename DenseMatrixPolicy>
    void InitializeConstraintParameters(
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& state_parameters) const
    {
      using FnType = std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>;
      using CacheType = std::vector<FnType>;
      if (!cached_constraint_init_fns_.has_value())
      {
        CacheType cache;
        auto phase_prefixes = CollectPhaseStatePrefixes();
        ForEachConstraint(
            [&](const auto& c)
            {
              if constexpr (requires {
                              c.template InitializeConstraintParametersFunction<DenseMatrixPolicy>(
                                  phase_prefixes, state_parameter_indices_, state_variable_indices_);
                            })
              {
                cache.push_back(c.template InitializeConstraintParametersFunction<DenseMatrixPolicy>(
                    phase_prefixes, state_parameter_indices_, state_variable_indices_));
              }
            });
        cached_constraint_init_fns_ = std::move(cache);
      }
      for (auto& fn : std::any_cast<CacheType&>(cached_constraint_init_fns_))
        fn(state_variables, state_parameters);
    }

    /// @brief Solve-time: add constraint residual G(y) to algebraic forcing rows
    template<typename DenseMatrixPolicy>
    void AddConstraintResidual(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const
    {
      using SP = detail::MatchingSparseT<DenseMatrixPolicy>;
      using LinearSetT = LinearConstraintSet<DenseMatrixPolicy, SP>;
      using DissEqSetT = DissolvedEquilibriumConstraintSet<DenseMatrixPolicy, SP>;
      using HenrysEqSetT = HenrysLawEquilibriumConstraintSet<DenseMatrixPolicy, SP>;
      if (linear_constraint_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<LinearSetT>&>(linear_constraint_sets_any_);
        for (const auto& set : sets)
          set.AddResidual(state_variables, state_parameters, forcing);
      }
      if (dissolved_equilibrium_constraint_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<DissEqSetT>&>(dissolved_equilibrium_constraint_sets_any_);
        for (const auto& set : sets)
          set.AddResidual(state_variables, state_parameters, forcing);
      }
      if (henrys_law_equilibrium_constraint_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<HenrysEqSetT>&>(henrys_law_equilibrium_constraint_sets_any_);
        for (const auto& set : sets)
          set.AddResidual(state_variables, state_parameters, forcing);
      }
    }

    /// @brief Solve-time: subtract dG/dy from algebraic Jacobian rows (-J convention)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractConstraintJacobian(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      using LinearSetT = LinearConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>;
      using DissEqSetT = DissolvedEquilibriumConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>;
      using HenrysEqSetT = HenrysLawEquilibriumConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>;
      if (linear_constraint_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<LinearSetT>&>(linear_constraint_sets_any_);
        for (const auto& set : sets)
          set.SubtractJacobian(state_variables, state_parameters, jacobian);
      }
      if (dissolved_equilibrium_constraint_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<DissEqSetT>&>(dissolved_equilibrium_constraint_sets_any_);
        for (const auto& set : sets)
          set.SubtractJacobian(state_variables, state_parameters, jacobian);
      }
      if (henrys_law_equilibrium_constraint_sets_any_.has_value())
      {
        auto& sets = std::any_cast<std::vector<HenrysEqSetT>&>(henrys_law_equilibrium_constraint_sets_any_);
        for (const auto& set : sets)
          set.SubtractJacobian(state_variables, state_parameters, jacobian);
      }
    }

   private:
    /// @brief Iterate over all registered processes with a generic callable
    template<typename Func>
    void ForEachProcess(Func&& fn) const
    {
      for (const auto& process : processes_)
      {
        std::visit([&](const auto& p) { fn(p); }, process);
      }
    }

    /// @brief Iterate over all registered constraints with a generic callable
    template<typename Func>
    void ForEachConstraint(Func&& fn) const
    {
      for (const auto& c : constraints_)
      {
        std::visit([&](const auto& cv) { fn(cv); }, c);
      }
    }

    /// @brief Return a reference to the lazily-built, DP-typed aerosol descriptor map used by
    ///        the cached AddForcingTerms / SubtractJacobianTerms paths. Rebuilds if the currently
    ///        cached DP does not match `DenseMatrixPolicy`.
    template<typename DenseMatrixPolicy>
    const std::map<std::string, std::map<AerosolProperty, AerosolPropertyDescriptor<DenseMatrixPolicy>>>&
    GetCachedDescriptors() const
    {
      using DescriptorMap = std::map<std::string, std::map<AerosolProperty, AerosolPropertyDescriptor<DenseMatrixPolicy>>>;
      if (!cached_descriptors_.has_value() || cached_descriptors_.type() != typeid(DescriptorMap))
      {
        auto phase_prefixes = CollectPhaseStatePrefixes();
        cached_descriptors_ =
            BuildDescriptors<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices_, state_variable_indices_);
      }
      return std::any_cast<const DescriptorMap&>(cached_descriptors_);
    }

    /// @brief Build aerosol property descriptors for all processes.
    /// @details Queries RequiredAerosolProperties() on each process, finds the representation
    ///          that owns each phase prefix, and calls GetPropertyDescriptor() to create descriptors.
    template<typename DenseMatrixPolicy>
    std::map<std::string, std::map<AerosolProperty, AerosolPropertyDescriptor<DenseMatrixPolicy>>> BuildDescriptors(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      // Collect all required properties across all processes
      std::map<std::string, std::vector<AerosolProperty>> required;
      ForEachProcess(
          [&](const auto& process)
          {
            auto process_required = process.RequiredAerosolProperties();
            for (const auto& [phase_name, properties] : process_required)
              for (const auto& prop : properties)
              {
                auto& existing = required[phase_name];
                if (std::find(existing.begin(), existing.end(), prop) == existing.end())
                  existing.push_back(prop);
              }
          });

      std::map<std::string, std::map<AerosolProperty, AerosolPropertyDescriptor<DenseMatrixPolicy>>> result;
      if (required.empty())
        return result;

      for (const auto& [phase_name, properties] : required)
      {
        auto pp_it = phase_prefixes.find(phase_name);
        if (pp_it == phase_prefixes.end())
        {
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_PHASE_PREFIX,
              "BuildDescriptors: phase not found: " + phase_name);
        }

        for (const auto& prefix : pp_it->second)
        {
          for (const auto& repr : representations_)
          {
            std::visit(
                [&](const auto& r)
                {
                  auto repr_prefixes = r.PhaseStatePrefixes();
                  auto phase_it = repr_prefixes.find(phase_name);
                  if (phase_it != repr_prefixes.end() && phase_it->second.count(prefix))
                  {
                    for (const auto& prop : properties)
                    {
                      result[prefix][prop] = WidenDescriptorVariant<DenseMatrixPolicy>(
                          r.template GetPropertyDescriptor<DenseMatrixPolicy>(
                              prop, state_parameter_indices, state_variable_indices, phase_name));
                    }
                  }
                },
                repr);
          }
        }
      }
      return result;
    }

    std::map<std::string, std::size_t> CountPhaseInstances() const
    {
      std::map<std::string, std::size_t> phase_instance_counts;
      for (const auto& repr : representations_)
      {
        std::visit(
            [&](const auto& r)
            {
              auto counts = r.NumPhaseInstances();
              for (const auto& [phase_name, count] : counts)
              {
                phase_instance_counts[phase_name] += count;  // Sum instances across representations
              }
            },
            repr);
      }
      return phase_instance_counts;
    }

    std::map<std::string, std::set<std::string>> CollectPhaseStatePrefixes() const
    {
      std::map<std::string, std::set<std::string>> phase_prefixes;
      for (const auto& repr : representations_)
      {
        std::visit(
            [&](const auto& r)
            {
              auto prefixes = r.PhaseStatePrefixes();
              for (const auto& [phase_name, prefix_set] : prefixes)
              {
                phase_prefixes[phase_name].insert(prefix_set.begin(), prefix_set.end());
              }
            },
            repr);
      }
      /// Validate against expected count of phase instances to ensure uniqueness
      auto expected_counts = CountPhaseInstances();
      for (const auto& [phase_name, prefixes] : phase_prefixes)
      {
        auto expected_count_it = expected_counts.find(phase_name);
        if (expected_count_it != expected_counts.end())
        {
          if (prefixes.size() != expected_count_it->second)
          {
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_DUPLICATE_STATE_PREFIX,
                "PhaseStatePrefixes: Non-unique state variable prefixes detected for phase " + phase_name +
                    ". Check your aerosol representation and phase names for duplicates.");
          }
        }
      }
      return phase_prefixes;
    }
  };
}  // namespace miam