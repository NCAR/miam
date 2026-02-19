#include "policy.hpp"

#include <miam/moments/single.hpp>

#include <gtest/gtest.h>

TEST(MomentsSingle, StateSize)
{
    auto size = testStateSize<miam::moment::Single>();
    EXPECT_EQ(std::get<0>(size), 3);
    EXPECT_EQ(std::get<1>(size), 2);
}

TEST(MomentsSingle, StateVariableNames)
{
    auto names = testStateVariableNames<miam::moment::Single>();
    EXPECT_EQ(names.size(), 3);
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.PHASE1.SPECIES_A"));
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.PHASE1.SPECIES_B"));
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.PHASE2.SPECIES_C"));
}

TEST(MomentsSingle, StateParameterNames)
{
    auto names = testStateParameterNames<miam::moment::Single>();
    EXPECT_EQ(names.size(), 2);
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.GEOMETRIC_MEAN_RADIUS"));
    EXPECT_TRUE(names.contains("TEST_DISTRIBUTION.GEOMETRIC_STANDARD_DEVIATION"));
}