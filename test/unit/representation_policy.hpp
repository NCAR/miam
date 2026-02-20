// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>

#include <micm/system/phase.hpp>

#include <gtest/gtest.h>

constexpr std::size_t minStateSize = 3;
const std::string test_model_name = "TEST_MODEL";

std::vector<micm::Phase> getTestPhases()
{
    micm::Phase phase1{ "PHASE1", { { micm::Species{ "SPECIES_A" } }, { micm::Species{ "SPECIES_B" } } } };
    micm::Phase phase2{ "PHASE2", { { micm::Species{ "SPECIES_C" } } } };
    return { phase1, phase2 };
}

template<typename ModelPolicy>
std::tuple<std::size_t, std::size_t> testStateSize()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    EXPECT_GE(std::get<0>(model.StateSize()), minStateSize);
    EXPECT_GE(std::get<1>(model.StateSize()), 0);
    return model.StateSize();
}

template<typename ModelPolicy>
std::set<std::string> testStateVariableNames()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    const auto names = model.StateVariableNames();
    EXPECT_EQ(names.size(), std::get<0>(model.StateSize()));

    return names;
}

template<typename ModelPolicy>
std::set<std::string> testStateParameterNames()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    const auto names = model.StateParameterNames();
    EXPECT_EQ(names.size(), std::get<1>(model.StateSize()));

    return names;
}

template<typename ModelPolicy>
std::map<std::string, double> testDefaultParameters()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    auto params = model.DefaultParameters();
    auto param_names = model.StateParameterNames();
    
    // Make sure the number of default parameters matches the number of state parameters
    EXPECT_EQ(params.size(), param_names.size());
    
    // Make sure all state parameter names have default values
    for (const auto& name : param_names)
    {
        EXPECT_TRUE(params.find(name) != params.end()); // name should exist in params
    }
    return params;
}