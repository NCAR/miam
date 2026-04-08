// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/aerosol_property.hpp>
#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/constraints/henry_law_equilibrium_constraint.hpp>
#include <miam/constraints/linear_constraint.hpp>
#include <miam/process.hpp>
#include <miam/representation.hpp>

#include <micm/system/conditions.hpp>

#include <algorithm>
#include <any>
#include <concepts>
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
    using RepresentationVariant =
        std::variant<representation::SingleMomentMode, representation::TwoMomentMode, representation::UniformSection>;
    using ProcessVariant = std::variant<process::DissolvedReaction, process::DissolvedReversibleReaction, process::HenryLawPhaseTransfer>;
    using ConstraintVariant = std::variant<
        constraint::DissolvedEquilibriumConstraint,
        constraint::HenryLawEquilibriumConstraint,
        constraint::LinearConstraint>;

    std::string name_;
    std::vector<RepresentationVariant> representations_;
    std::vector<ProcessVariant> processes_{};
    std::vector<ConstraintVariant> constraints_{};

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
      requires(sizeof...(ProcessTypes) >= 1 &&
               (std::constructible_from<ProcessVariant, std::decay_t<ProcessTypes>> && ...))
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
      requires(sizeof...(ConstraintTypes) >= 1 &&
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
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
    {
      // Collect parameter update functions from all processes and return a combined function
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>> update_functions;
      ForEachProcess(
          [&](const auto& process)
          {
            auto update_fn =
                process.template UpdateStateParametersFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices);
            update_functions.push_back(update_fn);
          });
      // Collect constraint parameter update functions (e.g. K_eq, HLC*R*T)
      std::vector<std::function<void(const std::vector<micm::Conditions>&)>> constraint_update_functions;
      ForEachConstraint(
          [&](const auto& c)
          {
            constraint_update_functions.push_back(c.template UpdateConstraintParametersFunction<DenseMatrixPolicy>());
          });
      return [update_functions, constraint_update_functions](
                 const std::vector<micm::Conditions>& conditions, DenseMatrixPolicy& state_parameters)
      {
        for (const auto& fn : update_functions)
        {
          fn(conditions, state_parameters);
        }
        for (const auto& fn : constraint_update_functions)
        {
          fn(conditions);
        }
      };
    }

    /// @brief Returns a function that calculates forcing terms
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      // Collect forcing functions from all processes and return a combined function
      auto phase_prefixes = CollectPhaseStatePrefixes();
      auto providers = BuildProviders<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
      std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>>
          forcing_functions;
      ForEachProcess(
          [&](const auto& process)
          {
            auto forcing_fn = process.template ForcingFunction<DenseMatrixPolicy>(
                phase_prefixes, state_parameter_indices, state_variable_indices, providers);
            forcing_functions.push_back(forcing_fn);
          });
      return [forcing_functions](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms)
      {
        for (const auto& fn : forcing_functions)
        {
          fn(state_parameters, state_variables, forcing_terms);
        }
      };
    }

    /// @brief Returns a function that calculates Jacobian contributions
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      // Collect Jacobian functions from all processes and return a combined function
      auto phase_prefixes = CollectPhaseStatePrefixes();
      auto providers = BuildProviders<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
      std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>>
          jacobian_functions;
      ForEachProcess(
          [&](const auto& process)
          {
            auto jacobian_fn = process.template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
                phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);
            jacobian_functions.push_back(jacobian_fn);
          });
      return [jacobian_functions](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian)
      {
        for (const auto& fn : jacobian_functions)
        {
          fn(state_parameters, state_variables, jacobian);
        }
      };
    }

    // ── HasConstraints concept methods ──

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
    /// @details Also includes process Jacobian elements for algebraic rows, since the solver builder
    ///          removes process elements from algebraic rows but the process JacobianFunction still
    ///          needs those elements to exist in the sparse matrix for VectorIndex lookups.
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

      // Include process Jacobian elements whose rows are algebraic.
      // The solver builder filters process elements from algebraic rows, then merges
      // constraint elements. By including process elements here, they survive the filter.
      auto algebraic_names = ConstraintAlgebraicVariableNames();
      std::set<std::size_t> algebraic_rows;
      for (const auto& name : algebraic_names)
      {
        auto it = state_indices.find(name);
        if (it != state_indices.end())
          algebraic_rows.insert(it->second);
      }
      if (!algebraic_rows.empty())
      {
        ForEachProcess(
            [&](const auto& process)
            {
              auto p_elements = process.NonZeroJacobianElements(phase_prefixes, state_indices);
              for (const auto& elem : p_elements)
              {
                if (algebraic_rows.count(elem.first) > 0)
                  elements.insert(elem);
              }
            });
      }

      return elements;
    }

    /// @brief Returns combined constraint residual function G(y) = 0
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::vector<std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>> residual_fns;
      ForEachConstraint(
          [&](const auto& c)
          {
            residual_fns.push_back(
                c.template ConstraintResidualFunction<DenseMatrixPolicy>(phase_prefixes, state_variable_indices));
          });
      return [residual_fns](const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& residual)
      {
        for (const auto& fn : residual_fns)
          fn(state_variables, residual);
      };
    }

    /// @brief Returns combined constraint Jacobian function (subtracts dG/dy)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian) const
    {
      auto phase_prefixes = CollectPhaseStatePrefixes();
      std::vector<std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)>> jac_fns;
      ForEachConstraint(
          [&](const auto& c)
          {
            jac_fns.push_back(
                c.template ConstraintJacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
                    phase_prefixes, state_variable_indices, jacobian));
          });
      return [jac_fns](const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian_values)
      {
        for (const auto& fn : jac_fns)
          fn(state_variables, jacobian_values);
      };
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

    /// @brief Build aerosol property providers for all processes
    /// @details Queries RequiredAerosolProperties() on each process, finds the representation
    ///          that owns each phase prefix, and calls GetPropertyProvider() to create providers.
    template<typename DenseMatrixPolicy>
    std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> BuildProviders(
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

      std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> result;
      if (required.empty())
        return result;

      for (const auto& [phase_name, properties] : required)
      {
        auto pp_it = phase_prefixes.find(phase_name);
        if (pp_it == phase_prefixes.end())
          throw std::runtime_error("BuildProviders: phase not found: " + phase_name);

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
                      result[prefix][prop] = r.template GetPropertyProvider<DenseMatrixPolicy>(
                          prop, state_parameter_indices, state_variable_indices, phase_name);
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
            throw std::runtime_error(
                "Internal Error: PhaseStatePrefixes: Non-unique state variable prefixes detected for phase " + phase_name);
          }
        }
      }
      return phase_prefixes;
    }
  };
}  // namespace miam