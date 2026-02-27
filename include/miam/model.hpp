// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/representation.hpp>
#include <miam/process.hpp>
#include <micm/system/conditions.hpp>

#include <functional>
#include <memory>
#include <set>
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
        using RepresentationVariant = std::variant<
            representation::SingleMomentMode,
            representation::TwoMomentMode,
            representation::UniformSection
        >;

        std::string name_;
        std::vector<RepresentationVariant> representations_;
        std::vector<process::DissolvedReversibleReaction> dissolved_reactions_{};

        /// @brief Returns the total state size (number of variables, number of parameters)
        std::tuple<std::size_t, std::size_t> StateSize() const
        {
            std::size_t num_variables = 0;
            std::size_t num_parameters = 0;
            for (const auto& repr : representations_)
            {
                std::visit([&](const auto& r) {
                    auto [vars, params] = r.StateSize();
                    num_variables += vars;
                    num_parameters += params;
                }, repr);
            }
            // Add parameters for each process
            auto phase_prefixes = CollectPhaseStatePrefixes();
            for (const auto& reaction : dissolved_reactions_)
            {
                auto process_params = reaction.ProcessParameterNames(phase_prefixes);
                num_parameters += process_params.size();
            }
            return { num_variables, num_parameters };
        }

        /// @brief Returns unique names for all state variables
        std::set<std::string> StateVariableNames() const
        {
            std::set<std::string> names;
            for (const auto& repr : representations_)
            {
                std::visit([&](const auto& r) {
                    auto repr_names = r.StateVariableNames();
                    names.insert(repr_names.begin(), repr_names.end());
                }, repr);
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
                std::visit([&](const auto& r) {
                    auto repr_names = r.StateParameterNames();
                    names.insert(repr_names.begin(), repr_names.end());
                }, repr);
            }
            // Add parameters for each process
            auto phase_prefixes = CollectPhaseStatePrefixes();
            for (const auto& reaction : dissolved_reactions_)            {
                auto process_params = reaction.ProcessParameterNames(phase_prefixes);
                names.insert(process_params.begin(), process_params.end());
            }
            return names;
        }

        /// @brief Returns names of all species used in the model's processes
        std::set<std::string> SpeciesUsed() const
        {
            auto phase_prefixes = CollectPhaseStatePrefixes();
            // Collect participating species' unique state names for all processes
            std::set<std::string> species_names;
            for (const auto& reaction : dissolved_reactions_)
            {
                auto process_species = reaction.SpeciesUsed(phase_prefixes);
                species_names.insert(process_species.begin(), process_species.end());
            }
            return species_names;
        }

        /// @brief Add dissolved reversible reactions to the model
        void AddProcesses(const std::vector<process::DissolvedReversibleReaction>& new_reactions)
        {
            // Create copies with new UUIDs for each reaction to ensure uniqueness across models
            for (const auto& reaction : new_reactions)
            {
                dissolved_reactions_.push_back(reaction.CopyWithNewUuid());
            }
        }

        /// @brief Returns non-zero Jacobian element positions
        std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
            const std::unordered_map<std::string, std::size_t>& state_indices) const
        {
            // Collect needed Jacobian element indices from all processes
            std::set<std::pair<std::size_t, std::size_t>> elements;
            auto phase_prefixes = CollectPhaseStatePrefixes();
            for (const auto& reaction : dissolved_reactions_)            {
                auto process_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);
                elements.insert(process_elements.begin(), process_elements.end());
            }
            return elements;
        }

        /// @brief Returns a function that updates state parameters
        template<typename DenseMatrixPolicy>
        std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> 
        UpdateStateParametersFunction(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
        {
            // Collect parameter update functions from all processes and return a combined function
            auto phase_prefixes = CollectPhaseStatePrefixes();
            std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>> update_functions;
            for (const auto& reaction : dissolved_reactions_)            {
                auto update_fn = reaction.template UpdateStateParametersFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices);
                update_functions.push_back(update_fn);
            }
            return [update_functions](const std::vector<micm::Conditions>& conditions, DenseMatrixPolicy& state_parameters) {
                for (const auto& fn : update_functions)
                {
                    fn(conditions, state_parameters);
                }
            };
        }

        /// @brief Returns a function that calculates forcing terms
        template<typename DenseMatrixPolicy>
        std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> 
        ForcingFunction(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
            const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
        {
            // Collect forcing functions from all processes and return a combined function
            auto phase_prefixes = CollectPhaseStatePrefixes();
            std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>> forcing_functions;
            for (const auto& reaction : dissolved_reactions_)            {
                auto forcing_fn = reaction.template ForcingFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
                forcing_functions.push_back(forcing_fn);
            }
            return [forcing_functions](const DenseMatrixPolicy& state_parameters, 
                     const DenseMatrixPolicy& state_variables, 
                     DenseMatrixPolicy& forcing_terms) {
                for (const auto& fn : forcing_functions)
                {
                    fn(state_parameters, state_variables, forcing_terms);
                }
            };
        }

        /// @brief Returns a function that calculates Jacobian contributions
        template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
        std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> 
        JacobianFunction(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
            const std::unordered_map<std::string, std::size_t>& state_variable_indices,
            const SparseMatrixPolicy& jacobian) const
        {
            // Collect Jacobian functions from all processes and return a combined function
            auto phase_prefixes = CollectPhaseStatePrefixes();
            std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>> jacobian_functions;
            for (const auto& reaction : dissolved_reactions_)            {
                auto jacobian_fn = reaction.template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
                jacobian_functions.push_back(jacobian_fn);
            }
            return [jacobian_functions](const DenseMatrixPolicy& state_parameters, 
                     const DenseMatrixPolicy& state_variables, 
                     SparseMatrixPolicy& jacobian) {
                for (const auto& fn : jacobian_functions)
                {
                    fn(state_parameters, state_variables, jacobian);
                }
            };
        }
    private:

        std::map<std::string, std::size_t> CountPhaseInstances() const
        {
            std::map<std::string, std::size_t> phase_instance_counts;
            for (const auto& repr : representations_)
            {
                std::visit([&](const auto& r) {
                    auto counts = r.NumPhaseInstances();
                    for (const auto& [phase_name, count] : counts)
                    {
                        phase_instance_counts[phase_name] += count; // Sum instances across representations
                    }
                }, repr);
            }
            return phase_instance_counts;
        }

        std::map<std::string, std::set<std::string>> CollectPhaseStatePrefixes() const
        {
            std::map<std::string, std::set<std::string>> phase_prefixes;
            for (const auto& repr : representations_)
            {
                std::visit([&](const auto& r) {
                    auto prefixes = r.PhaseStatePrefixes();
                    for (const auto& [phase_name, prefix_set] : prefixes)
                    {
                        phase_prefixes[phase_name].insert(prefix_set.begin(), prefix_set.end());
                    }
                }, repr);
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
                        throw std::runtime_error("Internal Error: PhaseStatePrefixes: Non-unique state variable prefixes detected for phase " + phase_name);
                    }
                }
            }
            return phase_prefixes;
        }
    };
}