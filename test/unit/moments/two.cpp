#include "policy.hpp"

#include <miam/moments/two.hpp>

#include <gtest/gtest.h>

TEST(MomentsTwo, StateSize)
{
    auto size = testStateSize<miam::moment::Two>();
    EXPECT_EQ(size, 4);
}

TEST(MomentsTwo, UniqueNames)
{
    auto names = testUniqueNames<miam::moment::Two>();
    EXPECT_EQ(names.size(), 4);
    EXPECT_EQ(names[0], "TEST_DISTRIBUTION.PHASE1.SPECIES_A");
    EXPECT_EQ(names[1], "TEST_DISTRIBUTION.PHASE1.SPECIES_B");
    EXPECT_EQ(names[2], "TEST_DISTRIBUTION.PHASE2.SPECIES_C");
    EXPECT_EQ(names[3], "TEST_DISTRIBUTION.number_concentration");
}