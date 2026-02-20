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
        /// @brief Sectional particle size distribution representation with uniform sections
        /// @details Represents a sectional distribution with uniform sections for aerosol or cloud particle size distributions.
        ///          Each section is characterized by a fixed size range and variable total volume. Number concentrations are
        ///          derived from the total volume and section size.
        class UniformSection
        {
        public:
            UniformSection() = delete;

            UniformSection(const std::string& prefix, const std::vector<micm::Phase>& phases)
                : prefix_(prefix), phases_(phases), default_min_radius_(0.0), default_max_radius_(0.0) {}

            UniformSection(const std::string& prefix, const std::vector<micm::Phase>& phases, const double minimum_radius, const double maximum_radius)
                : prefix_(prefix), phases_(phases), default_min_radius_(minimum_radius), default_max_radius_(maximum_radius) {}

            std::tuple<std::size_t, std::size_t> StateSize() const
            {
                std::size_t size = 0;
                for (const auto& phase : phases_)
                {
                    size += phase.StateSize();
                }
                return { size, 2 }; // Two parameters: min and max radius
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
                names.insert(prefix_ + ".MIN_RADIUS");
                names.insert(prefix_ + ".MAX_RADIUS");
                return names;
            }

            std::string Species(const micm::Phase& phase, const micm::Species& species) const
            {
                return prefix_ + "." + phase.name_ + "." + species.name_;
            }

            std::map<std::string, double> DefaultParameters() const
            {
                return {
                    { prefix_ + ".MIN_RADIUS", default_min_radius_ },
                    { prefix_ + ".MAX_RADIUS", default_max_radius_ }
                };
            }

            std::string MinRadius() const
            {
                return prefix_ + ".MIN_RADIUS";
            }

            std::string MaxRadius() const
            {
                return prefix_ + ".MAX_RADIUS";
            }

            void SetDefaultParameters(auto& state) const
            {
                auto min_radius_it = state.custom_rate_parameter_map_.find(MinRadius());
                if (min_radius_it == state.custom_rate_parameter_map_.end())
                {
                    throw std::runtime_error("Custom parameter map missing MIN_RADIUS for " + prefix_);
                }
                auto max_radius_it = state.custom_rate_parameter_map_.find(MaxRadius());
                if (max_radius_it == state.custom_rate_parameter_map_.end())
                {
                    throw std::runtime_error("Custom parameter map missing MAX_RADIUS for " + prefix_);
                }
                for (std::size_t cell = 0; cell < state.variables_.NumRows(); ++cell)
                {
                    state.custom_rate_parameters_[cell][min_radius_it->second] = default_min_radius_;
                    state.custom_rate_parameters_[cell][max_radius_it->second] = default_max_radius_;
                }
            }

        private:
            std::string prefix_; // State name prefix to apply to section properties
            std::vector<micm::Phase> phases_; // Phases associated with the section
            double default_min_radius_; // Minimum radius of the section
            double default_max_radius_; // Maximum radius of the section
        };
    }
}