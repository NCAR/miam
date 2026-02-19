// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>

#include <micm/system/phase.hpp>

#include <gtest/gtest.h>

constexpr std::size_t minStateSize = 3;

std::vector<micm::Phase> getTestPhases()
{
    micm::Phase phase1{ "PHASE1", { { micm::Species{ "SPECIES_A" } }, { micm::Species{ "SPECIES_B" } } } };
    micm::Phase phase2{ "PHASE2", { { micm::Species{ "SPECIES_C" } } } };
    return { phase1, phase2 };
}

std::vector<std::string> getTestMoments()
{
    return { "VOLUME", "NUMBER_CONCENTRATION", "GEOMETRIC_MEAN_RADIUS", "GEOMETRIC_STANDARD_DEVIATION" };
}

template<typename MomentPolicy>
std::tuple<std::size_t, std::size_t> testStateSize()
{
    const auto phases = getTestPhases();
    const auto moments = getTestMoments();
    EXPECT_GE(std::get<0>(MomentPolicy::StateSize(phases, moments)), minStateSize);
    EXPECT_GE(std::get<1>(MomentPolicy::StateSize(phases, moments)), 0);
    return MomentPolicy::StateSize(phases, moments);
}

template<typename MomentPolicy>
std::set<std::string> testStateVariableNames()
{
    const auto phases = getTestPhases();
    const auto moments = getTestMoments();
    const auto names = MomentPolicy::StateVariableNames("TEST_DISTRIBUTION", phases, moments);
    EXPECT_EQ(names.size(), std::get<0>(MomentPolicy::StateSize(phases, moments)));
    
    return names;
}

template<typename MomentPolicy>
std::set<std::string> testStateParameterNames()
{
    const auto phases = getTestPhases();
    const auto moments = getTestMoments();
    const auto names = MomentPolicy::StateParameterNames("TEST_DISTRIBUTION", phases, moments);
    EXPECT_EQ(names.size(), std::get<1>(MomentPolicy::StateSize(phases, moments)));

    return names;
}
