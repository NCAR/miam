// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/phase.hpp>

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
            ~Single() = default;

            /// @brief Returns the number of state variables needed to describe the distribution
            /// @param phases Phases associated with the distribution
            /// @return Number of state variables
            static constexpr std::size_t StateSize(const std::vector<micm::Phase>& phases)
            {
                std::size_t size = 0;
                for (const auto& phase : phases)
                {
                    size += phase.StateSize();
                }
                return size;
            }

            /// @brief Returns a vector of strings with unique names for each state variable
            /// @param prefix Prefix for variable names (e.g., distribution name)
            /// @param phases Phases associated with the distribution
            /// @return Vector of variable names
            static std::vector<std::string> StateVariableNames(const std::string& prefix, const std::vector<micm::Phase>& phases)
            {
                std::vector<std::string> names;
                for (const auto& phase : phases)
                {
                    for (const auto& species : phase.UniqueNames())
                    {
                        names.push_back(prefix + "." + species);
                    }
                }
                return names;
            }
        };
    }
}