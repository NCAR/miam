// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace miam
{
    namespace shape
    {

        /// @brief Delta-function size distribution shape
        /// @details Represents a delta-function distribution shape for aerosol or cloud particle size distributions.
        ///          All particles are assumed to have the same size.
        class DeltaFunction
        {
        public:
            DeltaFunction() = delete;
            
            /// @brief Constructor with specified prefix
            /// @param prefix Prefix to apply to shape properties in state variable names
            explicit DeltaFunction(const std::string& prefix) : prefix_(prefix) {}

            /// @brief Returns a vector of strings with names for each potential moment for the delta-function shape
            /// @return Vector of moment names
            static std::vector<std::string> PossibleMoments()
            {
                return { VOLUME_MOMENT, NUMBER_CONCENTRATION_MOMENT, RADIUS_MOMENT };
            }

            /// @brief Returns the state name for Number Concentration
            /// @return State variable name for Number Concentration
            std::string NumberConcentration() const
            {
                return prefix_ + "." + NUMBER_CONCENTRATION_MOMENT;
            }

            /// @brief Returns the state name for Radius
            /// @return State variable name for Radius
            std::string Radius() const
            {
                return prefix_ + "." + RADIUS_MOMENT;
            }

        private:
            std::string prefix_; // State name prefix to apply to shape properties
            static constexpr char VOLUME_MOMENT[] = "VOLUME";
            static constexpr char NUMBER_CONCENTRATION_MOMENT[] = "NUMBER_CONCENTRATION";
            static constexpr char RADIUS_MOMENT[] = "RADIUS";
        };
    }
}