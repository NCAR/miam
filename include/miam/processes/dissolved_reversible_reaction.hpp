// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/conditions.hpp>

#include <functional>

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

        };
    }
}