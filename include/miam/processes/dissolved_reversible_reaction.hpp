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

        };
    }
}