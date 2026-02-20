// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <micm/system/conditions.hpp>

#include <functional>

namespace miam
{
    namespace process
    {
        /// @brief A dissolved reversible reaction builder
        /// @details Builder class for constructing DissolvedReversibleReaction objects.
        class DissolvedReversibleReactionBuilder
        {
        public:
            DissolvedReversibleReactionBuilder() = default;

            /// @brief Sets the phase in which the reaction occurs
            DissolvedReversibleReactionBuilder& SetPhase(const micm::Phase& phase)
            {
                phase_ = phase;
                return *this;
            }

            /// @brief Sets the reactant species
            DissolvedReversibleReactionBuilder& SetReactants(const std::vector<micm::Species>& reactants)
            {
                reactants_ = reactants;
                return *this;
            }

            /// @brief Sets the product species
            DissolvedReversibleReactionBuilder& SetProducts(const std::vector<micm::Species>& products)
            {
                products_ = products;
                return *this;
            }

            /// @brief Sets the solvent species
            DissolvedReversibleReactionBuilder& SetSolvent(const micm::Species& solvent)
            {
                solvent_ = solvent;
                return *this;
            }

            /// @brief Sets the forward rate constant function
            DissolvedReversibleReactionBuilder& SetForwardRateConstant(const auto& forward_rate_constant)
            {
                forward_rate_constant_ = [forward_rate_constant](const micm::Conditions& conditions)
                {
                    return forward_rate_constant.Calculate(conditions);
                };
                return *this;
            }

            /// @brief Sets the reverse rate constant function
            DissolvedReversibleReactionBuilder& SetReverseRateConstant(const auto& reverse_rate_constant)
            {
                reverse_rate_constant_ = [reverse_rate_constant](const micm::Conditions& conditions)
                {
                    return reverse_rate_constant.Calculate(conditions);
                };
                return *this;
            }

            /// @brief Sets the equilibrium constant function
            DissolvedReversibleReactionBuilder& SetEquilibriumConstant(const auto& equilibrium_constant)
            {
                equilibrium_constant_ = [equilibrium_constant](const micm::Conditions& conditions)
                {
                    return equilibrium_constant.Calculate(conditions);
                };
                return *this;
            }

            /// @brief Builds and returns the DissolvedReversibleReaction object
            process::DissolvedReversibleReaction Build() const
            {
                // Check that exactly two of the three rate constant/equilibrium constant functions are set
                int num_set = 0;
                if (forward_rate_constant_) ++num_set;
                if (reverse_rate_constant_) ++num_set;
                if (equilibrium_constant_) ++num_set;
                if (num_set != 2)
                {
                    throw std::runtime_error("DissolvedReversibleReactionBuilder requires exactly two of forward rate constant, reverse rate constant, or equilibrium constant to be set.");
                }
                // If equilibrium constant is set, compute the missing rate constant
                if (equilibrium_constant_)
                {
                    if (!forward_rate_constant_)
                    {
                        forward_rate_constant_ = [this](const micm::Conditions& conditions)
                        {
                            double K_eq = equilibrium_constant_(conditions);
                            double k_r = reverse_rate_constant_(conditions);
                            return K_eq * k_r;
                        };
                    }
                    else if (!reverse_rate_constant_)
                    {
                        reverse_rate_constant_ = [this](const micm::Conditions& conditions)
                        {
                            double K_eq = equilibrium_constant_(conditions);
                            double k_f = forward_rate_constant_(conditions);
                            return k_f / K_eq;
                        };
                    }
                }
                return process::DissolvedReversibleReaction(
                    forward_rate_constant_,
                    reverse_rate_constant_,
                    reactants_,
                    products_,
                    solvent_,
                    phase_);
            }

        private:
            micm::Phase phase_;                     ///< Phase in which the reaction occurs
            std::vector<micm::Species> reactants_;  ///< Reactant species
            std::vector<micm::Species> products_;   ///< Product species
            micm::Species solvent_;                 ///< Solvent species
            mutable std::function<double(const micm::Conditions& conditions)> forward_rate_constant_; ///< Forward rate constant function
            mutable std::function<double(const micm::Conditions& conditions)> reverse_rate_constant_; ///< Reverse rate constant function
            mutable std::function<double(const micm::Conditions& conditions)> equilibrium_constant_;  ///< Equilibrium constant function
        };
    }
}