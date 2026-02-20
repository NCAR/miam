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
    namespace representation
    {
        /// @brief Two moment log-normal particle size distribution representation
        /// @details Represents a two moment log-normal distribution for aerosol or cloud particle size distributions.
        ///          Characterized by number concentration, geometric mean radius, and geometric standard deviation.
        class TwoMomentMode
        {
        public:
            TwoMomentMode() = delete;

            TwoMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases)
                : prefix_(prefix), phases_(phases), default_geometric_standard_deviation_(1.0) {}

            TwoMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases, const double geometric_standard_deviation)
                : prefix_(prefix), phases_(phases), default_geometric_standard_deviation_(geometric_standard_deviation) {}

            std::tuple<std::size_t, std::size_t> StateSize() const
            {
                std::size_t size = 0;
                for (const auto& phase : phases_)
                {
                    size += phase.StateSize();
                }
                size++; // Number concentration
                return { size, 1 }; // One parameter: geometric standard deviation
            }

            std::set<std::string> StateVariableNames() const
            {
                std::set<std::string> names;
                for (const auto& phase : phases_)
                {
                    for (const auto& species : phase.UniqueNames())
                    {
                        names.insert(prefix_ + "." + species);
                    }
                }
                names.insert(prefix_ + ".NUMBER_CONCENTRATION");
                return names;
            }

            std::set<std::string> StateParameterNames() const
            {
                std::set<std::string> names;
                names.insert(prefix_ + ".GEOMETRIC_STANDARD_DEVIATION");
                return names;
            }

            std::string Species(const micm::Phase& phase, const micm::Species& species) const
            {
                return prefix_ + "." + phase.name_ + "." + species.name_;
            }

            std::map<std::string, double> DefaultParameters() const
            {
                return {
                    { prefix_ + ".GEOMETRIC_STANDARD_DEVIATION", default_geometric_standard_deviation_ }
                };
            }

            std::string NumberConcentration() const
            {
                return prefix_ + ".NUMBER_CONCENTRATION";
            }

            std::string GeometricStandardDeviation() const
            {
                return prefix_ + ".GEOMETRIC_STANDARD_DEVIATION";
            }

            void SetDefaultParameters(auto& state) const
            {
                auto gsd_it = state.custom_rate_parameter_map_.find(GeometricStandardDeviation());
                if (gsd_it == state.custom_rate_parameter_map_.end())
                {
                    throw std::runtime_error("TwoMomentMode::SetDefaultParameters: Geometric standard deviation parameter not found in state.");
                }
                for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
                {
                    state.custom_rate_parameters_[i_cell][gsd_it->second] = default_geometric_standard_deviation_;
                }
            }

        private:
            std::string prefix_; // State name prefix to apply to shape properties
            std::vector<micm::Phase> phases_; // Phases associated with the mode
            double default_geometric_standard_deviation_; // Default geometric standard deviation
        };
    }
}