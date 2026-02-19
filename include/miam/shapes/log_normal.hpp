// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace miam
{
    namespace shape
    {

        /// @brief Log-normal size distribution shape
        /// @details Represents a log-normal distribution shape for aerosol or cloud particle size distributions.
        ///          Characterized by a geometric mean diameter and geometric standard deviation.
        class LogNormal
        {
        public:
            LogNormal() = delete;
            
            /// @brief Constructor with specified prefix
            /// @param prefix Prefix to apply to shape properties in state variable names
            explicit LogNormal(const std::string& prefix) : prefix_(prefix) {}

            /// @brief Returns a vector of strings with names for each potential moment for the log-normal shape
            /// @return Vector of moment names
            static std::vector<std::string> PossibleMoments()
            {
                return { VOLUME_MOMENT, NUMBER_CONCENTRATION_MOMENT, GEOMETRIC_MEAN_RADIUS_MOMENT, GEOMETRIC_STANDARD_DEVIATION_MOMENT };
            }

            /// @brief Returns the state name for Number Concentration
            /// @return State variable name for Number Concentration
            std::string NumberConcentration() const
            {
                return prefix_ + "." + NUMBER_CONCENTRATION_MOMENT;
            }

            /// @brief Returns the state name for Gemoetric Mean Radius
            /// @return State variable name for Geometric Mean Radius
            std::string GeometricMeanRadius() const
            {
                return prefix_ + "." + GEOMETRIC_MEAN_RADIUS_MOMENT;
            }

            /// @brief Returns the state name for Geometric Standard Deviation
            /// @return State variable name for Geometric Standard Deviation
            std::string GeometricStandardDeviation() const
            {
                return prefix_ + "." + GEOMETRIC_STANDARD_DEVIATION_MOMENT;
            }

        private:
            std::string prefix_; // State name prefix to apply to shape properties
            static constexpr char VOLUME_MOMENT[] = "VOLUME";
            static constexpr char NUMBER_CONCENTRATION_MOMENT[] = "NUMBER_CONCENTRATION";
            static constexpr char GEOMETRIC_MEAN_RADIUS_MOMENT[] = "GEOMETRIC_MEAN_RADIUS";
            static constexpr char GEOMETRIC_STANDARD_DEVIATION_MOMENT[] = "GEOMETRIC_STANDARD_DEVIATION";
        };
    }
}