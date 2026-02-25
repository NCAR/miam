// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace miam
{
    namespace process
    {
        /// @brief A dissolved reversible reaction
        /// @details Dissolved reversible reactions involve reactants and products in solution, and
        ///          are characterized by both a forward and reverse rate constant. If an equilibrium
        ///          constant is provided, only one of the forward or reverse rate constants needs to be
        ///          specified, and the other is computed from the equilibrium constant.
        ///
        ///          The reaction takes the form:
        ///          \f$ \mathrm{Reactants} \leftrightarrow \mathrm{Products} \f$
        ///
        ///          with forward rate constant \f$ k_f \f$ and reverse rate constant \f$ k_r \f$. The
        ///          relationship between the rate constants and the equilibrium constant \f$ K_{eq} \f$ is:
        ///          \f$ K_{eq} = \frac{k_f}{k_r} = \frac{\prod[\mathrm{Products}]}{\prod[\mathrm{Reactants}]} \f$.
        struct DissolvedReversibleReaction
        {
            std::function<double(const micm::Conditions& conditions)> forward_rate_constant_; ///< Forward rate constant function
            std::function<double(const micm::Conditions& conditions)> reverse_rate_constant_; ///< Reverse rate constant function
            std::vector<micm::Species> reactants_;  ///< Reactant species
            std::vector<micm::Species> products_;   ///< Product species
            micm::Species solvent_;                 ///< Solvent species
            micm::Phase phase_;                     ///< Phase in which the reaction occurs

            DissolvedReversibleReaction() = delete;

            /// @brief Constructor
            DissolvedReversibleReaction(
                std::function<double(const micm::Conditions& conditions)> forward_rate_constant,
                std::function<double(const micm::Conditions& conditions)> reverse_rate_constant,
                const std::vector<micm::Species>& reactants,
                const std::vector<micm::Species>& products,
                micm::Species solvent,
                micm::Phase phase)
                : forward_rate_constant_(forward_rate_constant),
                  reverse_rate_constant_(reverse_rate_constant),
                  reactants_(reactants),
                  products_(products),
                  solvent_(solvent),
                  phase_(phase)
            {}

            /// @brief Returns participating species' unique state names
            /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or species names)
            /// @return Set of unique state variable names for all species involved in the reaction
            std::set<std::string> SpeciesUsed(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
            {
                std::set<std::string> species_names;
                auto phase_it = phase_prefixes.find(phase_.name_);
                if (phase_it != phase_prefixes.end())
                {
                    const auto& prefixes = phase_it->second;
                    for (const auto& prefix : prefixes)
                    {
                        for (const auto& reactant : reactants_)
                        {
                            species_names.insert(prefix + "." + phase_.name_ + "." + reactant.name_);
                        }
                        for (const auto& product : products_)
                        {
                            species_names.insert(prefix + "." + phase_.name_ + "." + product.name_);
                        }
                        species_names.insert(prefix + "." + phase_.name_ + "." + solvent_.name_);
                    }
                }
                return species_names;
            }

            /// @brief Returns a set of Jacobian index pairs for this process
            /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or species names)
            /// @param state_variable_indices Map of state variable names to their corresponding indices in the Jacobian
            /// @return Set of index pairs representing the Jacobian entries affected by this process (dependent, independent)
            std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
                const std::map<std::string, std::set<std::string>>& phase_prefixes,
                const auto& state_variable_indices // acts like std::unordered_map<std::string, std::size_t>
            ) const
            {
                std::set<std::pair<std::size_t, std::size_t>> jacobian_indices;
                
                // Check if this reaction's phase exists in the phase_prefixes
                auto phase_it = phase_prefixes.find(phase_.name_);
                if (phase_it == phase_prefixes.end())
                {
                    // Phase not present in this model, return empty set
                    return jacobian_indices;
                }
                
                for(const auto& prefix : phase_it->second)
                {
                    // Reversible reactions will have contributions from each reactant on each reactant and product variable,
                    // and from each product on each reactant and product variable, and from the solvent on all reactant and product variables.
                    for(const auto& reactant : reactants_)
                    {
                        std::string reactant_var = prefix + "." + phase_.name_ + "." + reactant.name_;
                        for(const auto& other_reactant : reactants_)
                        {
                            std::string other_reactant_var = prefix + "." + phase_.name_ + "." + other_reactant.name_;
                            jacobian_indices.insert({state_variable_indices.at(reactant_var), state_variable_indices.at(other_reactant_var)});
                        }
                        for(const auto& product : products_)
                        {
                            std::string product_var = prefix + "." + phase_.name_ + "." + product.name_;
                            jacobian_indices.insert({state_variable_indices.at(reactant_var), state_variable_indices.at(product_var)});
                        }
                        std::string solvent_var = prefix + "." + phase_.name_ + "." + solvent_.name_;
                        jacobian_indices.insert({state_variable_indices.at(reactant_var), state_variable_indices.at(solvent_var)});
                    }
                    for(const auto& product : products_)
                    {                        
                        std::string product_var = prefix + "." + phase_.name_ + "." + product.name_;
                        for(const auto& reactant : reactants_)                        {
                            std::string reactant_var = prefix + "." + phase_.name_ + "." + reactant.name_;
                            jacobian_indices.insert({state_variable_indices.at(product_var), state_variable_indices.at(reactant_var)});
                        }
                        for(const auto& other_product : products_)                        {
                            std::string other_product_var = prefix + "." + phase_.name_ + "." + other_product.name_;
                            jacobian_indices.insert({state_variable_indices.at(product_var), state_variable_indices.at(other_product_var)});
                        }
                        std::string solvent_var = prefix + "." + phase_.name_ + "." + solvent_.name_;
                        jacobian_indices.insert({state_variable_indices.at(product_var), state_variable_indices.at(solvent_var)});
                    }
                }
                return jacobian_indices;
            }                        

        };
    }
}