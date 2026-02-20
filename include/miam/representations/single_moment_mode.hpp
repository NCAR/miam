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
        /// @brief Single moment log-normal particle size distribution representation
        /// @details Represents a single moment log-normal distribution for aerosol or cloud particle size distributions.
        ///          Characterized by a geometric mean radius and geometric standard deviation.
        class SingleMomentMode
        {
        public:
            SingleMomentMode() = delete;

            SingleMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases)
                : prefix_(prefix), phases_(phases), default_geometric_mean_radius_(0.0), default_geometric_standard_deviation_(1.0) {}

            SingleMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases, const double geometric_mean_radius, const double geometric_standard_deviation)
                : prefix_(prefix), phases_(phases), default_geometric_mean_radius_(geometric_mean_radius), default_geometric_standard_deviation_(geometric_standard_deviation) {}

            std::tuple<std::size_t, std::size_t> StateSize() const
            {
                std::size_t size = 0;
                for (const auto& phase : phases_)
                {
                    size += phase.StateSize();
                }
                return { size, 2 }; // Two parameters: geometric mean radius and geometric standard deviation
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
                return names;
            }

            std::set<std::string> StateParameterNames() const
            {
                std::set<std::string> names;
                names.insert(prefix_ + ".GEOMETRIC_MEAN_RADIUS");
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
                    { prefix_ + ".GEOMETRIC_MEAN_RADIUS", default_geometric_mean_radius_ },
                    { prefix_ + ".GEOMETRIC_STANDARD_DEVIATION", default_geometric_standard_deviation_ }
                };
            }

            std::string GeometricMeanRadius() const
            {
                return prefix_ + ".GEOMETRIC_MEAN_RADIUS";
            }

            std::string GeometricStandardDeviation() const
            {
                return prefix_ + ".GEOMETRIC_STANDARD_DEVIATION";
            }

            void SetDefaultParameters(auto& state) const
            {
                auto gmd_it = state.custom_rate_parameter_map_.find(GeometricMeanRadius());
                if (gmd_it == state.custom_rate_parameter_map_.end())
                {
                    throw std::runtime_error("Geometric mean radius parameter not found in state for " + prefix_);
                }
                auto gsd_it = state.custom_rate_parameter_map_.find(GeometricStandardDeviation());
                if (gsd_it == state.custom_rate_parameter_map_.end())
                {
                    throw std::runtime_error("Geometric standard deviation parameter not found in state for " + prefix_);
                }
                for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
                {
                    state.custom_rate_parameters_[i_cell][gmd_it->second] = default_geometric_mean_radius_;
                    state.custom_rate_parameters_[i_cell][gsd_it->second] = default_geometric_standard_deviation_;
                }
            }

        private:
            std::string prefix_; // State name prefix to apply to mode properties
            std::vector<micm::Phase> phases_; // Phases associated with the mode
            double default_geometric_mean_radius_; // Default geometric mean radius
            double default_geometric_standard_deviation_; // Default geometric standard deviation
        };
    }
}