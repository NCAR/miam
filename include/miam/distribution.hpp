// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/phase.hpp>

#include <vector>

namespace miam
{

    /// @brief A homogeneous population of particles characterized by a size distribution
    /// @details Particles within a distribution are considered to be of the same composition and physical properties,
    ///          differing only in volume (size). The shape of the distribution is used to determine average properties
    ///          such as effective radius and surface area. The moment determines which properties (e.g., mass, number
    ///          concentration, effective radius) are tracked in the state distribution, and which are fixed in time.
    ///
    ///          Aerosol distributions are often combined into sets to represent the full aerosol particle population
    ///          in a system, and can be used to represent any suspension of fine solid or liquid particles in a gas,
    ///          including cloud, rain, or ice droplets, as well as smaller dry or aqueous aerosol particles.
    template<typename ShapeType, typename MomentType>
    class Distribution
    {
        // ShapeType is reserved for future use to represent the distribution's shape
        static_assert(
            sizeof(ShapeType) >= 0,
            "ShapeType template parameter is reserved for future use to represent the distribution shape.");
        /// @brief Name for this distribution
        std::string name_;
        /// @brief Phases associated with this distribution (e.g., aqueous, organic)
        std::vector<micm::Phase> phases_;

    public:
        Distribution() = delete;

        /// @brief Construct a Distribution with specified phases, shape, and moment
        /// @param name Name of the distribution
        /// @param phases Phases associated with this distribution
        Distribution(const std::string& name, const std::vector<micm::Phase>& phases)
            : name_(name), phases_(phases)
        {}

        ~Distribution() = default;

        /// @brief Returns the number of state variables needed to describe the distribution
        /// @return Number of state variables
        std::size_t StateSize() const
        {
            return MomentType::StateSize(phases_);
        }

        /// @brief Returns a vector of strings with unique names for each state variable
        /// @return Vector of variable names
        std::vector<std::string> UniqueNames() const
        {
            return MomentType::StateVariableNames(name_, phases_);
        }
    };
}