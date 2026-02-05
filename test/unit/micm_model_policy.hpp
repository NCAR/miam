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
std::size_t testStateSize()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    EXPECT_GE(model.StateSize(), minStateSize);
    return model.StateSize();
}

template<typename ModelPolicy>
std::vector<std::string> testUniqueNames()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    const auto names = model.UniqueNames();
    EXPECT_EQ(names.size(), model.StateSize());

    // Check that names are unique
    std::set<std::string> unique_names(names.begin(), names.end());
    EXPECT_EQ(unique_names.size(), names.size());
    return names;
}
   