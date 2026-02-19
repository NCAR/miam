#include "policy.hpp"

#include <miam/moments/two.hpp>

#include <gtest/gtest.h>

TEST(MomentsTwo, StateSize)
{
    auto size = testStateSize<miam::moment::Two>();
    EXPECT_EQ(std::get<0>(size), 4);
    EXPECT_EQ(std::get<1>(size), 1);
}

TEST(MomentsTwo, StateVariableNames)
{
    auto names = testStateVariableNames<miam::moment::Two>();
    EXPECT_EQ(names.size(), 4);
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.PHASE1.SPECIES_A"));
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.PHASE1.SPECIES_B"));
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.PHASE2.SPECIES_C"));
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.NUMBER_CONCENTRATION"));
}

TEST(MomentsTwo, StateParameterNames)
{
    auto names = testStateParameterNames<miam::moment::Two>();
    EXPECT_EQ(names.size(), 1);
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.GEOMETRIC_STANDARD_DEVIATION"));
}