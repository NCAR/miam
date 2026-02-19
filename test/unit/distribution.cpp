// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "micm_model_policy.hpp"

#include <miam/distribution.hpp>
#include <miam/shape.hpp>
#include <miam/moment.hpp>

#include <gtest/gtest.h>

TEST(Distribution, StateSize)
{
    auto size = testStateSize<miam::Distribution<miam::shape::LogNormal, miam::moment::Two>>();
    EXPECT_EQ(std::get<0>(size), 4);
    EXPECT_EQ(std::get<1>(size), 1);

    size = testStateSize<miam::Distribution<miam::shape::DeltaFunction, miam::moment::Single>>();
    EXPECT_EQ(std::get<0>(size), 3);
    EXPECT_EQ(std::get<1>(size), 1);
}

TEST(Distribution, StateVariableNames)
{
    auto names = testStateVariableNames<miam::Distribution<miam::shape::LogNormal, miam::moment::Two>>();
    EXPECT_EQ(names.size(), 4);
    EXPECT_TRUE(names.contains("TEST_MODEL.PHASE1.SPECIES_A"));
    EXPECT_TRUE(names.contains("TEST_MODEL.PHASE1.SPECIES_B"));
    EXPECT_TRUE(names.contains("TEST_MODEL.PHASE2.SPECIES_C"));
    EXPECT_TRUE(names.contains("TEST_MODEL.NUMBER_CONCENTRATION"));

    names = testStateVariableNames<miam::Distribution<miam::shape::DeltaFunction, miam::moment::Single>>();
    EXPECT_EQ(names.size(), 3);
    EXPECT_TRUE(names.contains("TEST_MODEL.PHASE1.SPECIES_A"));
    EXPECT_TRUE(names.contains("TEST_MODEL.PHASE1.SPECIES_B"));
    EXPECT_TRUE(names.contains("TEST_MODEL.PHASE2.SPECIES_C"));
}

TEST(Distribution, StateParameterNames)
{
    auto names = testStateParameterNames<miam::Distribution<miam::shape::LogNormal, miam::moment::Two>>();
    EXPECT_EQ(names.size(), 1);
    EXPECT_TRUE(names.contains("TEST_MODEL.GEOMETRIC_STANDARD_DEVIATION"));

    names = testStateParameterNames<miam::Distribution<miam::shape::DeltaFunction, miam::moment::Single>>();
    EXPECT_EQ(names.size(), 1);
    EXPECT_TRUE(names.contains("TEST_MODEL.RADIUS"));
}