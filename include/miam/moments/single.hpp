// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/phase.hpp>

#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace miam
{
    namespace moment
    {

        /// @brief A single-moment scheme for tracking aerosol distribution state
        /// @details In a single-moment scheme, only the total mass concentration of the distribution is tracked in the state.
        ///          All other properties are derived from the mass concentration and fixed distribution parameters.
        class Single
        {
        public:
            Single() = default;

            /// @brief Returns the number of state variables and parameters needed to describe the distribution
            /// @param phases Phases associated with the distribution
            /// @param moments Moments associated with the distribution (the first moment must be "VOLUME")
            /// @return Number of state variables
            static constexpr std::tuple<std::size_t, std::size_t> StateSize(const std::vector<micm::Phase>& phases, std::vector<std::string> moments)
            {
                std::size_t size = 0;
                if (moments.size() < 2)
                {
                    throw std::invalid_argument("At least two possible moments must be specified for a single moment scheme.");
                }
                if (moments[0] != "VOLUME")
                {
                    throw std::invalid_argument("The first moment must be 'VOLUME' for Single moment scheme.");
                }
                for (const auto& phase : phases)
                {
                    size += phase.StateSize();
                }
                return { size, moments.size() - 2 }; // Exclude first two moments
            }

            /// @brief Returns a set of strings with unique names for each state variable
            /// @param prefix Prefix for variable names (e.g., distribution name)
            /// @param phases Phases associated with the distribution
            /// @param moments Moments associated with the distribution (the first moment must be "VOLUME")
            /// @return Set of variable names
            static std::set<std::string> StateVariableNames(const std::string& prefix, const std::vector<micm::Phase>& phases, const std::vector<std::string>& moments)
            {
                std::set<std::string> names;
                if (moments.size() < 2)
                {
                    throw std::invalid_argument("At least two possible moments must be specified for a single moment scheme.");
                }
                if (moments[0] != "VOLUME")
                {
                    throw std::invalid_argument("The first moment must be 'VOLUME' for Single moment scheme.");
                }
                for (const auto& phase : phases)
                {
                    for (const auto& species : phase.UniqueNames())
                    {
                        names.insert(prefix + "." + species);
                    }
                }
                return names;
            }

            /// @brief Returns a set of strings with unique names for each state parameter
            /// @param prefix Prefix for parameter names (e.g., distribution name)
            /// @param phases Phases associated with the distribution
            /// @param moments Moments associated with the distribution (the first moment must be "VOLUME")
            /// @return Set of parameter names
            static std::set<std::string> StateParameterNames(const std::string& prefix, const std::vector<micm::Phase>& phases, const std::vector<std::string>& moments)
            {
                std::set<std::string> names;
                if (moments.size() < 2)
                {
                    throw std::invalid_argument("At least two possible moments must be specified for a single moment scheme.");
                }
                if (moments[0] != "VOLUME")
                {
                    throw std::invalid_argument("The first moment must be 'VOLUME' for Single moment scheme.");
                }
                for (std::size_t i = 2; i < moments.size(); ++i)
                {
                    names.insert(prefix + "." + moments[i]);
                }
                return names;
            }

            /// @brief Returns the state variable name for a specific species in a phase
            /// @param prefix Prefix for variable names (e.g., distribution name)
            /// @param phase Phase associated with the species
            /// @param species Species name
            /// @return Variable name
            static std::string Species(const std::string& prefix, const micm::Phase& phase, const micm::Species& species)
            {
                return prefix + "." + phase.name_ + "." + species.name_;
            }
        };
    }
}