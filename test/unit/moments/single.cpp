#include "policy.hpp"

#include <miam/moments/single.hpp>

#include <gtest/gtest.h>

TEST(MomentsSingle, StateSize)
{
    auto size = testStateSize<miam::moment::Single>();
    EXPECT_EQ(size, 3);
}

TEST(MomentsSingle, UniqueNames)
{
    auto names = testUniqueNames<miam::moment::Single>();
    EXPECT_EQ(names.size(), 3);
    EXPECT_EQ(names[0], "TEST_DISTRIBUTION.PHASE1.SPECIES_A");
    EXPECT_EQ(names[1], "TEST_DISTRIBUTION.PHASE1.SPECIES_B");
    EXPECT_EQ(names[2], "TEST_DISTRIBUTION.PHASE2.SPECIES_C");
}