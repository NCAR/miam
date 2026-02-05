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
    EXPECT_EQ(size, 4);

    size = testStateSize<miam::Distribution<miam::shape::DeltaFunction, miam::moment::Single>>();
    EXPECT_EQ(size, 3);
}

TEST(Distribution, UniqueNames)
{
    auto names = testUniqueNames<miam::Distribution<miam::shape::LogNormal, miam::moment::Two>>();
    EXPECT_EQ(names.size(), 4);
    EXPECT_EQ(names[0], "TEST_MODEL.PHASE1.SPECIES_A");
    EXPECT_EQ(names[1], "TEST_MODEL.PHASE1.SPECIES_B");
    EXPECT_EQ(names[2], "TEST_MODEL.PHASE2.SPECIES_C");
    EXPECT_EQ(names[3], "TEST_MODEL.number_concentration");

    names = testUniqueNames<miam::Distribution<miam::shape::DeltaFunction, miam::moment::Single>>();
    EXPECT_EQ(names.size(), 3);
    EXPECT_EQ(names[0], "TEST_MODEL.PHASE1.SPECIES_A");
    EXPECT_EQ(names[1], "TEST_MODEL.PHASE1.SPECIES_B");
    EXPECT_EQ(names[2], "TEST_MODEL.PHASE2.SPECIES_C");
}