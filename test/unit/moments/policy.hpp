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

template<typename MomentPolicy>
std::size_t testStateSize()
{
    const auto phases = getTestPhases();
    EXPECT_GE(MomentPolicy::StateSize(phases), minStateSize);
    return MomentPolicy::StateSize(phases);
}

template<typename MomentPolicy>
std::vector<std::string> testUniqueNames()
{
    const auto phases = getTestPhases();
    const auto names = MomentPolicy::StateVariableNames("TEST_DISTRIBUTION", phases);
    EXPECT_EQ(names.size(), MomentPolicy::StateSize(phases));
    
    // Check that names are unique
    std::set<std::string> unique_names(names.begin(), names.end());
    EXPECT_EQ(unique_names.size(), names.size());
    return names;
}
