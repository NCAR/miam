// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/phase.hpp>

#include <set>
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
        // @brief Shape
        ShapeType shape_;
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
            : shape_(name), name_(name), phases_(phases)
        {}

        ~Distribution() = default;

        /// @brief Returns the number of state variables and parameters needed to describe the distribution
        /// @return Number of state variables
        std::tuple<std::size_t, std::size_t> StateSize() const
        {
            return MomentType::StateSize(phases_, ShapeType::PossibleMoments());
        }

        /// @brief Returns a set of strings with unique names for each state variable
        /// @return Set of variable names
        std::set<std::string> StateVariableNames() const
        {
            return MomentType::StateVariableNames(name_, phases_, ShapeType::PossibleMoments());
        }

        /// @brief Returns a set of strings with unique names for each state parameter
        std::set<std::string> StateParameterNames() const
        {
            return MomentType::StateParameterNames(name_, phases_, ShapeType::PossibleMoments());
        }

        /// @brief Returns a const reference to the shape associated with this distribution
        const ShapeType& Shape() const
        {
            return shape_;
        }

        /// @brief Returns the state variable name for a specific species in a specific phase
        /// @param phase Phase associated with the species
        /// @param species Species for which to get the state variable name
        /// @return State variable name
        std::string Species(const micm::Phase& phase, const micm::Species& species) const
        {
            return MomentType::Species(name_, phase, species);
        }

        /// @brief Sets a state parameter value
        /// @param state Reference to the state array
        /// @param parameter_name Name of the parameter to set
        /// @param value Value to set
        void SetParameter(auto& state, const std::string& parameter_name, double value) const
        {
            auto param_it = state.custom_rate_parameter_map_.find(parameter_name);
            if (param_it == state.custom_rate_parameter_map_.end())
            {
                throw std::invalid_argument("Parameter name " + parameter_name + " not found in state.");
            }
            if (state.custom_rate_parameters_.NumRows() != 1)
            {
                throw std::invalid_argument("Cannot apply scalar value to multiple grid cell state parameters.");
            }
            state.custom_rate_parameters_[0][param_it->second] = value;
        }

        /// @brief Sets a state parameter value for multiple grid cells
        /// @param state Reference to the state array
        /// @param parameter_name Name of the parameter to set
        /// @param values Vector of values to set for each grid cell
        void SetParameter(auto& state, const std::string& parameter_name, const std::vector<double>& values) const
        {
            auto param_it = state.custom_rate_parameter_map_.find(parameter_name);
            if (param_it == state.custom_rate_parameter_map_.end())
            {
                throw std::invalid_argument("Parameter name " + parameter_name + " not found in state.");
            }
            if (state.custom_rate_parameters_.NumRows() != values.size())
            {
                throw std::invalid_argument("Size of values vector does not match number of grid cells in state.");
            }
            for (std::size_t cell = 0; cell < state.custom_rate_parameters_.NumRows(); ++cell)
            {
                state.custom_rate_parameters_[cell][param_it->second] = values[cell];
            }
        }
    };
}