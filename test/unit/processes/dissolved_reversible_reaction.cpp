// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

using namespace miam::process;

TEST(DissolvedReversibleReaction, SpeciesUsedWithSinglePrefix)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },      // reactants
        { hp, ohm },  // products
        h2o,          // solvent
        aqueous_phase
    };
    
    // Create phase prefixes map - single representation with one prefix
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
    
    auto species_used = reaction.SpeciesUsed(phase_prefixes);
    
    // Should have 3 species: reactant (H2O), 2 products (H+, OH-), and solvent (H2O - already in reactants)
    EXPECT_EQ(species_used.size(), 3);
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H+") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.OH-") != species_used.end());
}

TEST(DissolvedReversibleReaction, SpeciesUsedWithMultiplePrefixes)
{
    auto co2 = micm::Species{ "CO2" };
    auto h2o = micm::Species{ "H2O" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { co2 }, { h2o }, { h2co3 } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { co2, h2o },  // reactants
        { h2co3 },     // products
        h2o,           // solvent
        aqueous_phase
    };
    
    // Create phase prefixes map - multiple representations with different prefixes
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
    phase_prefixes["AQUEOUS"].insert("LARGE_DROP");
    
    auto species_used = reaction.SpeciesUsed(phase_prefixes);
    
    // Should have 3 species in each of 2 representations = 6 total
    EXPECT_EQ(species_used.size(), 6);
    
    // Check SMALL_DROP species
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.CO2") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H2CO3") != species_used.end());
    
    // Check LARGE_DROP species
    EXPECT_TRUE(species_used.find("LARGE_DROP.AQUEOUS.CO2") != species_used.end());
    EXPECT_TRUE(species_used.find("LARGE_DROP.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("LARGE_DROP.AQUEOUS.H2CO3") != species_used.end());
}

TEST(DissolvedReversibleReaction, SpeciesUsedWithNoMatchingPhase)
{
    auto co2 = micm::Species{ "CO2" };
    auto h2o = micm::Species{ "H2O" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { co2 }, { h2o }, { h2co3 } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { co2, h2o },  // reactants
        { h2co3 },     // products
        h2o,           // solvent
        aqueous_phase
    };
    
    // Create phase prefixes map without the AQUEOUS phase
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["ORGANIC"].insert("AEROSOL_MODE");
    
    auto species_used = reaction.SpeciesUsed(phase_prefixes);
    
    // Should be empty since no matching phase
    EXPECT_EQ(species_used.size(), 0);
}

TEST(DissolvedReversibleReaction, SpeciesUsedDuplicateHandling)
{
    // Test case where reactant and product overlap
    auto hp = micm::Species{ "H+" };
    auto hco3m = micm::Species{ "HCO3-" };
    auto co32m = micm::Species{ "CO32-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { hp }, { hco3m }, { co32m } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-10; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { hco3m },         // reactants
        { hp, co32m },     // products
        hco3m,             // solvent (same as reactant)
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    auto species_used = reaction.SpeciesUsed(phase_prefixes);
    
    // Should have 3 unique species (solvent is same as reactant, so no duplicate)
    EXPECT_EQ(species_used.size(), 3);
    EXPECT_TRUE(species_used.find("MODE1.AQUEOUS.HCO3-") != species_used.end());
    EXPECT_TRUE(species_used.find("MODE1.AQUEOUS.H+") != species_used.end());
    EXPECT_TRUE(species_used.find("MODE1.AQUEOUS.CO32-") != species_used.end());
}

TEST(DissolvedReversibleReaction, SpeciesUsedComplexReaction)
{
    // More complex reaction with multiple reactants and products
    auto a = micm::Species{ "A" };
    auto b = micm::Species{ "B" };
    auto c = micm::Species{ "C" };
    auto d = micm::Species{ "D" };
    auto solvent = micm::Species{ "SOLVENT" };
    
    auto phase = micm::Phase{ "LIQUID", { { a }, { b }, { c }, { d }, { solvent } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 2.0; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { a, b },      // 2 reactants
        { c, d },      // 2 products
        solvent,
        phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["LIQUID"].insert("REP1");
    phase_prefixes["LIQUID"].insert("REP2");
    phase_prefixes["LIQUID"].insert("REP3");
    
    auto species_used = reaction.SpeciesUsed(phase_prefixes);
    
    // Should have 5 species (A, B, C, D, SOLVENT) × 3 representations = 15 total
    EXPECT_EQ(species_used.size(), 15);
    
    // Verify all combinations exist
    for (const auto& prefix : {"REP1", "REP2", "REP3"})
    {
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".LIQUID.A") != species_used.end());
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".LIQUID.B") != species_used.end());
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".LIQUID.C") != species_used.end());
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".LIQUID.D") != species_used.end());
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".LIQUID.SOLVENT") != species_used.end());
    }
}

TEST(DissolvedReversibleReaction, NonZeroJacobianElementsBasic)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    // H2O <-> H+ + OH-
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },      // reactants
        { hp, ohm },  // products
        h2o,          // solvent
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_indices["MODE1.AQUEOUS.H+"] = 1;
    state_indices["MODE1.AQUEOUS.OH-"] = 2;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);
    
    // For reversible reaction H2O <-> H+ + OH-:
    // d[H2O]/d[H2O], d[H2O]/d[H+], d[H2O]/d[OH-]
    // d[H+]/d[H2O], d[H+]/d[H+], d[H+]/d[OH-]
    // d[OH-]/d[H2O], d[OH-]/d[H+], d[OH-]/d[OH-]
    // Total: 9 elements
    EXPECT_EQ(jacobian_elements.size(), 9);
    
    // Check specific elements (row, col)
    EXPECT_TRUE(jacobian_elements.find({0, 0}) != jacobian_elements.end()); // d[H2O]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({0, 1}) != jacobian_elements.end()); // d[H2O]/d[H+]
    EXPECT_TRUE(jacobian_elements.find({0, 2}) != jacobian_elements.end()); // d[H2O]/d[OH-]
    EXPECT_TRUE(jacobian_elements.find({1, 0}) != jacobian_elements.end()); // d[H+]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({1, 1}) != jacobian_elements.end()); // d[H+]/d[H+]
    EXPECT_TRUE(jacobian_elements.find({1, 2}) != jacobian_elements.end()); // d[H+]/d[OH-]
    EXPECT_TRUE(jacobian_elements.find({2, 0}) != jacobian_elements.end()); // d[OH-]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({2, 1}) != jacobian_elements.end()); // d[OH-]/d[H+]
    EXPECT_TRUE(jacobian_elements.find({2, 2}) != jacobian_elements.end()); // d[OH-]/d[OH-]
}

TEST(DissolvedReversibleReaction, NonZeroJacobianElementsMultipleReactants)
{
    auto co2 = micm::Species{ "CO2" };
    auto h2o = micm::Species{ "H2O" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { co2 }, { h2o }, { h2co3 } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    // CO2 + H2O <-> H2CO3
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { co2, h2o },  // 2 reactants
        { h2co3 },     // 1 product
        h2o,           // solvent
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("DROP");
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["DROP.AQUEOUS.CO2"] = 0;
    state_indices["DROP.AQUEOUS.H2O"] = 1;
    state_indices["DROP.AQUEOUS.H2CO3"] = 2;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);
    
    // For CO2 + H2O <-> H2CO3:
    // Reactant CO2 affects: CO2, H2O (other reactant), H2CO3 (product), H2O (solvent) = 3 unique
    // Reactant H2O affects: H2O, CO2 (other reactant), H2CO3 (product), H2O (solvent) = 3 unique
    // Product H2CO3 affects: H2CO3, CO2 (reactant), H2O (reactant), H2O (solvent) = 3 unique
    // Total: 9 elements
    EXPECT_EQ(jacobian_elements.size(), 9);
}

TEST(DissolvedReversibleReaction, NonZeroJacobianElementsMultiplePrefixes)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    // Two representations of the same phase
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
    phase_prefixes["AQUEOUS"].insert("LARGE_DROP");
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["SMALL_DROP.AQUEOUS.H2O"] = 0;
    state_indices["SMALL_DROP.AQUEOUS.H+"] = 1;
    state_indices["SMALL_DROP.AQUEOUS.OH-"] = 2;
    state_indices["LARGE_DROP.AQUEOUS.H2O"] = 3;
    state_indices["LARGE_DROP.AQUEOUS.H+"] = 4;
    state_indices["LARGE_DROP.AQUEOUS.OH-"] = 5;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);
    
    // 9 elements per representation × 2 representations = 18 elements
    EXPECT_EQ(jacobian_elements.size(), 18);
    
    // Check a few elements from each representation
    EXPECT_TRUE(jacobian_elements.find({0, 0}) != jacobian_elements.end()); // SMALL_DROP
    EXPECT_TRUE(jacobian_elements.find({1, 2}) != jacobian_elements.end()); // SMALL_DROP
    EXPECT_TRUE(jacobian_elements.find({3, 3}) != jacobian_elements.end()); // LARGE_DROP
    EXPECT_TRUE(jacobian_elements.find({4, 5}) != jacobian_elements.end()); // LARGE_DROP
}

TEST(DissolvedReversibleReaction, NonZeroJacobianElementsComplexReaction)
{
    // HCO3- <-> H+ + CO32-
    auto hco3m = micm::Species{ "HCO3-" };
    auto hp = micm::Species{ "H+" };
    auto co32m = micm::Species{ "CO32-" };
    auto h2o = micm::Species{ "H2O" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { hco3m }, { hp }, { co32m }, { h2o } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-10; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { hco3m },
        { hp, co32m },
        h2o,
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE");
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["MODE.AQUEOUS.HCO3-"] = 0;
    state_indices["MODE.AQUEOUS.H+"] = 1;
    state_indices["MODE.AQUEOUS.CO32-"] = 2;
    state_indices["MODE.AQUEOUS.H2O"] = 3;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);
    
    // Reactant HCO3- affects itself, H+, CO32-, and H2O (solvent) = 4
    // Product H+ affects all reactants and products and solvent = 4
    // Product CO32- affects all reactants and products and solvent = 4
    // Total unique: 12 elements
    EXPECT_EQ(jacobian_elements.size(), 12);
    
    // Verify some specific elements
    EXPECT_TRUE(jacobian_elements.find({0, 0}) != jacobian_elements.end()); // d[HCO3-]/d[HCO3-]
    EXPECT_TRUE(jacobian_elements.find({0, 1}) != jacobian_elements.end()); // d[HCO3-]/d[H+]
    EXPECT_TRUE(jacobian_elements.find({0, 2}) != jacobian_elements.end()); // d[HCO3-]/d[CO32-]
    EXPECT_TRUE(jacobian_elements.find({0, 3}) != jacobian_elements.end()); // d[HCO3-]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({1, 0}) != jacobian_elements.end()); // d[H+]/d[HCO3-]
    EXPECT_TRUE(jacobian_elements.find({2, 0}) != jacobian_elements.end()); // d[CO32-]/d[HCO3-]
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
