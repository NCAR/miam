// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/constraints/dissolved_equilibrium_constraint_evaluator.hpp>
#include <miam/constraints/henrys_law_equilibrium_constraint.hpp>
#include <miam/constraints/henrys_law_equilibrium_constraint_evaluator.hpp>
#include <miam/constraints/linear_constraint.hpp>
#include <miam/constraints/linear_constraint_evaluator.hpp>
#include <miam/processes.hpp>
#include <miam/processes/dissolved_reaction_evaluator.hpp>
#include <miam/processes/dissolved_reversible_reaction_evaluator.hpp>
#include <miam/processes/henrys_law_phase_transfer_evaluator.hpp>
#include <miam/representations.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/system/conditions.hpp>
#include <micm/util/matrix.hpp>

#include <algorithm>
#include <any>
#include <concepts>
#include <cstddef>
#include <functional>
#include <map>
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

    // Per-process / per-constraint parameter update closures, built at Finalize for the solver's
    // DenseMatrixPolicy. Each holds a `std::vector<UpdateFn<DP>>` or `std::vector<InitFn<DP>>`.
    std::any process_update_fns_{};
    std::any constraint_update_fns_{};
    std::any constraint_init_fns_{};

    // Solve-time evaluators, built by FinalizeProcessSetup / FinalizeConstraintSetup for the
    // solver's <DenseMatrixPolicy, SparseMatrixPolicy> pair. `process_evaluators_` and
    // `constraint_evaluators_` hold a `std::shared_ptr<const ProcessEvaluators<DP, SP>>` /
    // `ConstraintEvaluators<DP, SP>`, so copies of the Model share them. The `*_fn_` members hold
    // DP-only closures over the same evaluators for the solve-time methods that MICM calls with
    // only DP.
    std::any process_evaluators_{};
    std::any process_forcing_fn_{};
    std::any constraint_evaluators_{};
    std::any constraint_residual_fn_{};

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
    /// @details `SparseMatrixPolicy` is required so that process evaluators can pre-compute Jacobian flat IDs at build time.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      auto nz_elements = NonZeroJacobianElements(state_variable_indices);
      auto jacobian_builder = SparseMatrixPolicy::Create(state_variable_indices.size()).InitialValue(0.0);
      for (const auto& elem : nz_elements)
        jacobian_builder = jacobian_builder.WithElement(elem.first, elem.second);
      SparseMatrixPolicy jacobian_pattern(jacobian_builder);

      auto evaluators = BuildProcessEvaluators<DenseMatrixPolicy>(
          phase_prefixes, state_parameter_indices, state_variable_indices, jacobian_pattern);
      return [evaluators](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms)
      { evaluators->AddForcingTerms(state_parameters, state_variables, forcing_terms); };
    }

    /// @brief Returns a function that calculates Jacobian contributions
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      auto evaluators = BuildProcessEvaluators<DenseMatrixPolicy>(
          CollectPhaseStatePrefixes(), state_parameter_indices, state_variable_indices, jacobian);
      return [evaluators](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian)
      { evaluators->SubtractJacobianTerms(state_parameters, state_variables, jacobian); };
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
    ///          species map, and Jacobian sparsity pattern are finalized. The builder
    ///          supplies DenseMatrixPolicy; SparseMatrixPolicy is deduced from `jacobian`.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void FinalizeProcessSetup(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      state_parameter_indices_ = state_parameter_indices;
      state_variable_indices_ = state_variable_indices;

      auto phase_prefixes = CollectPhaseStatePrefixes();

      std::vector<UpdateFn<DenseMatrixPolicy>> update_fns;
      ForEachProcess(
          [&](const auto& process)
          {
            update_fns.push_back(process.template UpdateStateParametersFunction<DenseMatrixPolicy>(
                phase_prefixes, state_parameter_indices));
          });
      process_update_fns_ = std::move(update_fns);

      auto shared_evaluators = BuildProcessEvaluators<DenseMatrixPolicy>(
          phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
      process_forcing_fn_ = ForcingFn<DenseMatrixPolicy>(
          [shared_evaluators](
              const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing)
          { shared_evaluators->AddForcingTerms(state_parameters, state_variables, forcing); });
      process_evaluators_ = std::move(shared_evaluators);
    }

    /// @brief Cache build-time indices for constraint solve-time methods
    /// @details The builder supplies DenseMatrixPolicy; SparseMatrixPolicy is deduced from `jacobian`.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void FinalizeConstraintSetup(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      state_parameter_indices_ = state_parameter_indices;
      state_variable_indices_ = state_variable_indices;

      auto phase_prefixes = CollectPhaseStatePrefixes();

      std::vector<UpdateFn<DenseMatrixPolicy>> update_fns;
      std::vector<InitFn<DenseMatrixPolicy>> init_fns;
      ForEachConstraint(
          [&](const auto& c)
          {
            update_fns.push_back(
                c.template UpdateConstraintParametersFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices));
            if constexpr (requires {
                            c.template InitializeConstraintParametersFunction<DenseMatrixPolicy>(
                                phase_prefixes, state_parameter_indices, state_variable_indices);
                          })
            {
              init_fns.push_back(c.template InitializeConstraintParametersFunction<DenseMatrixPolicy>(
                  phase_prefixes, state_parameter_indices, state_variable_indices));
            }
          });
      constraint_update_fns_ = std::move(update_fns);
      constraint_init_fns_ = std::move(init_fns);

      auto evaluators = std::make_shared<ConstraintEvaluators<DenseMatrixPolicy, SparseMatrixPolicy>>();
      for (const auto& constraint : constraints_)
      {
        std::visit(
            [&](const auto& c)
            {
              using C = std::decay_t<decltype(c)>;
              if constexpr (std::is_same_v<C, LinearConstraint>)
                evaluators->linear_.emplace_back(c, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              else if constexpr (std::is_same_v<C, DissolvedEquilibriumConstraint>)
                evaluators->dissolved_equilibrium_.emplace_back(
                    c, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              else if constexpr (std::is_same_v<C, HenrysLawEquilibriumConstraint>)
                evaluators->henrys_law_equilibrium_.emplace_back(
                    c, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
            },
            constraint);
      }

      std::shared_ptr<const ConstraintEvaluators<DenseMatrixPolicy, SparseMatrixPolicy>> shared_evaluators = std::move(evaluators);
      constraint_residual_fn_ = ForcingFn<DenseMatrixPolicy>(
          [shared_evaluators](
              const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing)
          { shared_evaluators->AddResidual(state_parameters, state_variables, forcing); });
      constraint_evaluators_ = std::move(shared_evaluators);
    }

    /// @brief Solve-time: refresh temperature-/pressure-dependent process parameters
    template<typename DenseMatrixPolicy>
    void UpdateStateParameters(
        const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
        DenseMatrixPolicy& state_parameters) const
    {
      if (process_update_fns_.has_value())
        for (const auto& fn : std::any_cast<const std::vector<UpdateFn<DenseMatrixPolicy>>&>(process_update_fns_))
          fn(conditions, state_parameters);
    }

    /// @brief Solve-time: add process forcing (tendency) contributions
    template<typename DenseMatrixPolicy>
    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const
    {
      if (process_forcing_fn_.has_value())
        std::any_cast<const ForcingFn<DenseMatrixPolicy>&>(process_forcing_fn_)(state_parameters, state_variables, forcing);
    }

    /// @brief Solve-time: subtract process Jacobian contributions (-J convention)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      using EvaluatorsPtr = std::shared_ptr<const ProcessEvaluators<DenseMatrixPolicy, SparseMatrixPolicy>>;
      if (process_evaluators_.has_value())
        std::any_cast<const EvaluatorsPtr&>(process_evaluators_)->SubtractJacobianTerms(state_parameters, state_variables, jacobian);
    }

    /// @brief Solve-time: refresh temperature-/pressure-dependent constraint parameters
    template<typename DenseMatrixPolicy>
    void UpdateConstraintStateParameters(
        const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
        DenseMatrixPolicy& state_parameters) const
    {
      if (constraint_update_fns_.has_value())
        for (const auto& fn : std::any_cast<const std::vector<UpdateFn<DenseMatrixPolicy>>&>(constraint_update_fns_))
          fn(conditions, state_parameters);
    }

    /// @brief Solve-time: diagnose constraint parameters from current state
    template<typename DenseMatrixPolicy>
    void InitializeConstraintParameters(
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& state_parameters) const
    {
      if (constraint_init_fns_.has_value())
        for (const auto& fn : std::any_cast<const std::vector<InitFn<DenseMatrixPolicy>>&>(constraint_init_fns_))
          fn(state_variables, state_parameters);
    }

    /// @brief Solve-time: add constraint residual G(y) to algebraic forcing rows
    template<typename DenseMatrixPolicy>
    void AddConstraintResidual(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const
    {
      if (constraint_residual_fn_.has_value())
        std::any_cast<const ForcingFn<DenseMatrixPolicy>&>(constraint_residual_fn_)(
            state_parameters, state_variables, forcing);
    }

    /// @brief Solve-time: subtract dG/dy from algebraic Jacobian rows (-J convention)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    void SubtractConstraintJacobian(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      using EvaluatorsPtr = std::shared_ptr<const ConstraintEvaluators<DenseMatrixPolicy, SparseMatrixPolicy>>;
      if (constraint_evaluators_.has_value())
        std::any_cast<const EvaluatorsPtr&>(constraint_evaluators_)->SubtractJacobian(state_parameters, state_variables, jacobian);
    }

   private:
    template<typename DenseMatrixPolicy>
    using DescriptorMap = std::map<std::string, std::map<AerosolProperty, AerosolPropertyDescriptor<DenseMatrixPolicy>>>;

    template<typename DenseMatrixPolicy>
    using UpdateFn =
        std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>;

    template<typename DenseMatrixPolicy>
    using InitFn = std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>;

    template<typename DenseMatrixPolicy>
    using ForcingFn = std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>;

    /// @brief Process evaluators built by FinalizeProcessSetup for one <DP, SP> pair
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    struct ProcessEvaluators
    {
      std::vector<DissolvedReactionEvaluator<DenseMatrixPolicy, SparseMatrixPolicy>> dissolved_reactions_;
      std::vector<DissolvedReversibleReactionEvaluator<DenseMatrixPolicy, SparseMatrixPolicy>> dissolved_reversible_reactions_;
      std::vector<HenrysLawPhaseTransferEvaluator<DenseMatrixPolicy, SparseMatrixPolicy>> henrys_law_phase_transfers_;
      DescriptorMap<DenseMatrixPolicy> descriptors_;

      void AddForcingTerms(
          const DenseMatrixPolicy& state_parameters,
          const DenseMatrixPolicy& state_variables,
          DenseMatrixPolicy& forcing) const
      {
        for (const auto& evaluator : dissolved_reactions_)
          evaluator.AddForcingTerms(state_parameters, state_variables, forcing);
        for (const auto& evaluator : dissolved_reversible_reactions_)
          evaluator.AddForcingTerms(state_parameters, state_variables, forcing);
        for (const auto& evaluator : henrys_law_phase_transfers_)
          evaluator.AddForcingTerms(state_parameters, state_variables, forcing, descriptors_);
      }

      void SubtractJacobianTerms(
          const DenseMatrixPolicy& state_parameters,
          const DenseMatrixPolicy& state_variables,
          SparseMatrixPolicy& jacobian) const
      {
        for (const auto& evaluator : dissolved_reactions_)
          evaluator.SubtractJacobianTerms(state_parameters, state_variables, jacobian);
        for (const auto& evaluator : dissolved_reversible_reactions_)
          evaluator.SubtractJacobianTerms(state_parameters, state_variables, jacobian);
        for (const auto& evaluator : henrys_law_phase_transfers_)
          evaluator.SubtractJacobianTerms(state_parameters, state_variables, jacobian, descriptors_);
      }
    };

    /// @brief Constraint evaluators built by FinalizeConstraintSetup for one <DP, SP> pair
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    struct ConstraintEvaluators
    {
      std::vector<LinearConstraintEvaluator<DenseMatrixPolicy, SparseMatrixPolicy>> linear_;
      std::vector<DissolvedEquilibriumConstraintEvaluator<DenseMatrixPolicy, SparseMatrixPolicy>> dissolved_equilibrium_;
      std::vector<HenrysLawEquilibriumConstraintEvaluator<DenseMatrixPolicy, SparseMatrixPolicy>> henrys_law_equilibrium_;

      void AddResidual(
          const DenseMatrixPolicy& state_parameters,
          const DenseMatrixPolicy& state_variables,
          DenseMatrixPolicy& forcing) const
      {
        for (const auto& evaluator : linear_)
          evaluator.AddResidual(state_variables, state_parameters, forcing);
        for (const auto& evaluator : dissolved_equilibrium_)
          evaluator.AddResidual(state_variables, state_parameters, forcing);
        for (const auto& evaluator : henrys_law_equilibrium_)
          evaluator.AddResidual(state_variables, state_parameters, forcing);
      }

      void SubtractJacobian(
          const DenseMatrixPolicy& state_parameters,
          const DenseMatrixPolicy& state_variables,
          SparseMatrixPolicy& jacobian) const
      {
        for (const auto& evaluator : linear_)
          evaluator.SubtractJacobian(state_variables, state_parameters, jacobian);
        for (const auto& evaluator : dissolved_equilibrium_)
          evaluator.SubtractJacobian(state_variables, state_parameters, jacobian);
        for (const auto& evaluator : henrys_law_equilibrium_)
          evaluator.SubtractJacobian(state_variables, state_parameters, jacobian);
      }
    };

    /// @brief Builds the process evaluators for one <DP, SP> pair. The result is shared, because
    ///        the evaluators hold views into their own index storage and must not be copied.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::shared_ptr<const ProcessEvaluators<DenseMatrixPolicy, SparseMatrixPolicy>> BuildProcessEvaluators(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      // Reference (CPU) descriptor map used only to enumerate `n_deps` counts for
      // Jacobian flat-ID layout in HenrysLawPhaseTransferEvaluator.
      using ReferenceDense = micm::Matrix<double>;
      auto reference_descriptors =
          BuildDescriptors<ReferenceDense>(phase_prefixes, state_parameter_indices, state_variable_indices);

      auto evaluators = std::make_shared<ProcessEvaluators<DenseMatrixPolicy, SparseMatrixPolicy>>();
      for (const auto& process : processes_)
      {
        std::visit(
            [&](const auto& p)
            {
              using P = std::decay_t<decltype(p)>;
              if constexpr (std::is_same_v<P, DissolvedReaction>)
              {
                evaluators->dissolved_reactions_.emplace_back(
                    p, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              }
              else if constexpr (std::is_same_v<P, DissolvedReversibleReaction>)
              {
                evaluators->dissolved_reversible_reactions_.emplace_back(
                    p, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
              }
              else if constexpr (std::is_same_v<P, HenrysLawPhaseTransfer>)
              {
                evaluators->henrys_law_phase_transfers_.emplace_back(
                    p, phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, reference_descriptors);
              }
            },
            process);
      }
      if (!evaluators->henrys_law_phase_transfers_.empty())
        evaluators->descriptors_ =
            BuildDescriptors<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);

      return evaluators;
    }

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