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

        /// @brief A two-moment scheme for tracking aerosol distribution state
        /// @details In a two-moment scheme, both the total mass concentration
        ///          and number concentration of the distribution are tracked in
        ///          the state. All other properties are derived from these two
        ///          moments and fixed distribution parameters.
        class Two
        {
        public:
            Two() = default;

            /// @brief Returns the number of state variables and parameters needed to describe the distribution
            /// @param phases Phases associated with the distribution
            /// @param moments Moments associated with the distribution (the first moment must be "VOLUME")
            /// @return Number of state variables
            static constexpr std::tuple<std::size_t, std::size_t> StateSize(const std::vector<micm::Phase>& phases, const std::vector<std::string>& moments)
            {
                std::size_t size = 0;
                if (moments.size() < 3)
                {
                    throw std::invalid_argument("At least three possible moments must be specified for a two moment scheme.");
                }
                if (moments[0] != "VOLUME")
                {
                    throw std::invalid_argument("The first moment must be 'VOLUME' for Two moment scheme.");
                }
                for (const auto& phase : phases)
                {
                    size += phase.StateSize();
                }
                // Add one state variable for second moment
                size += 1;
                return { size, moments.size() - 3 }; // Exclude first three moments
            }

            /// @brief Returns a set of strings with unique names for each state variable
            /// @param prefix Prefix for variable names (e.g., distribution name)
            /// @param phases Phases associated with the distribution
            /// @param moments Moments associated with the distribution (the first moment must be "VOLUME")
            /// @return Set of variable names
            static std::set<std::string> StateVariableNames(const std::string& prefix, const std::vector<micm::Phase>& phases, const std::vector<std::string>& moments)
            {
                std::set<std::string> names;
                if (moments.size() < 3)
                {
                    throw std::invalid_argument("At least three possible moments must be specified for Two moment scheme.");
                }
                if (moments[0] != "VOLUME")
                {
                    throw std::invalid_argument("The first moment must be 'VOLUME' for Two moment scheme.");
                }
                for (const auto& phase : phases)
                {
                    for (const auto& species : phase.UniqueNames())
                    {
                        names.insert(prefix + "." + species);
                    }
                }
                names.insert(prefix + "." + moments[1]);
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
                if (moments.size() < 3)
                {
                    throw std::invalid_argument("At least three possible moments must be specified for Two moment scheme.");
                }
                if (moments[0] != "VOLUME")
                {
                    throw std::invalid_argument("The first moment must be 'VOLUME' for Two moment scheme.");
                }
                for (std::size_t i = 3; i < moments.size(); ++i)
                {
                    names.insert(prefix + "." + moments[i]);
                }
                return names;
            }

            /// @brief Returns the state variable name for a specific species in a phase
            /// @param prefix Prefix for variable names (e.g., distribution name)
            /// @param phase Phase associated with the species
            /// @param species Species name
            /// @return State variable name
            static std::string Species(const std::string& prefix, const micm::Phase& phase, const micm::Species& species)
            {
                return prefix + "." + phase.name_ + "." + species.name_;
            }
        };
    }
}