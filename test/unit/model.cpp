// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/model.hpp>
#include <miam/representations/single_moment_mode.hpp>
#include <miam/representations/two_moment_mode.hpp>
#include <miam/representations/uniform_section.hpp>
#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

using namespace miam;

TEST(Model, SpeciesUsedWithNoProcesses)
{
    auto h2o = micm::Species{ "H2O" };
    auto co2 = micm::Species{ "CO2" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { co2 } } };
    
    auto mode = representation::SingleMomentMode{ "MODE1", { aqueous_phase } };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    
    auto species_used = model.SpeciesUsed();
    
    // Should be empty since no processes
    EXPECT_EQ(species_used.size(), 0);
}

TEST(Model, SpeciesUsedWithSingleRepresentationAndProcess)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto mode = representation::SingleMomentMode{ "SMALL_DROP", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    process::DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },      // reactants
        { hp, ohm },  // products
        h2o,          // solvent
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    model.processes_.push_back(reaction);
    
    auto species_used = model.SpeciesUsed();
    
    // Should have 3 species from the reaction
    EXPECT_EQ(species_used.size(), 3);
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H+") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.OH-") != species_used.end());
}

TEST(Model, SpeciesUsedWithMultipleRepresentations)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto small_drop = representation::SingleMomentMode{ "SMALL_DROP", { aqueous_phase } };
    auto large_drop = representation::SingleMomentMode{ "LARGE_DROP", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    process::DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },      // reactants
        { hp, ohm },  // products
        h2o,          // solvent
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(small_drop);
    model.representations_.push_back(large_drop);
    model.processes_.push_back(reaction);
    
    auto species_used = model.SpeciesUsed();
    
    // Should have 3 species × 2 representations = 6 species
    EXPECT_EQ(species_used.size(), 6);
    
    // Small drop species
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.H+") != species_used.end());
    EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.OH-") != species_used.end());
    
    // Large drop species
    EXPECT_TRUE(species_used.find("LARGE_DROP.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("LARGE_DROP.AQUEOUS.H+") != species_used.end());
    EXPECT_TRUE(species_used.find("LARGE_DROP.AQUEOUS.OH-") != species_used.end());
}

TEST(Model, SpeciesUsedWithMultipleProcesses)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    auto co2 = micm::Species{ "CO2" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm }, { co2 }, { h2co3 } } };
    
    auto mode = representation::TwoMomentMode{ "AITKEN", { aqueous_phase } };
    
    auto h2o_forward = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto h2o_reverse = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    process::DissolvedReversibleReaction h2o_dissociation{
        h2o_forward,
        h2o_reverse,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    auto co2_forward = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto co2_reverse = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    process::DissolvedReversibleReaction co2_hydration{
        co2_forward,
        co2_reverse,
        { co2, h2o },
        { h2co3 },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    model.processes_.push_back(h2o_dissociation);
    model.processes_.push_back(co2_hydration);
    
    auto species_used = model.SpeciesUsed();
    
    // Should have unique species from both reactions: H2O, H+, OH-, CO2, H2CO3 = 5 species
    EXPECT_EQ(species_used.size(), 5);
    EXPECT_TRUE(species_used.find("AITKEN.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("AITKEN.AQUEOUS.H+") != species_used.end());
    EXPECT_TRUE(species_used.find("AITKEN.AQUEOUS.OH-") != species_used.end());
    EXPECT_TRUE(species_used.find("AITKEN.AQUEOUS.CO2") != species_used.end());
    EXPECT_TRUE(species_used.find("AITKEN.AQUEOUS.H2CO3") != species_used.end());
}

TEST(Model, SpeciesUsedMixedRepresentationTypes)
{
    auto h2o = micm::Species{ "H2O" };
    auto co2 = micm::Species{ "CO2" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { co2 }, { h2co3 } } };
    
    auto mode1 = representation::SingleMomentMode{ "MODE1", { aqueous_phase } };
    auto mode2 = representation::TwoMomentMode{ "MODE2", { aqueous_phase } };
    auto section = representation::UniformSection{ "BIN_01", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    process::DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { co2, h2o },
        { h2co3 },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode1);
    model.representations_.push_back(mode2);
    model.representations_.push_back(section);
    model.processes_.push_back(reaction);
    
    auto species_used = model.SpeciesUsed();
    
    // Should have 3 species × 3 representations = 9 species
    EXPECT_EQ(species_used.size(), 9);
    
    // Check each representation
    for (const auto& prefix : {"MODE1", "MODE2", "BIN_01"})
    {
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.H2O") != species_used.end());
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.CO2") != species_used.end());
        EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.H2CO3") != species_used.end());
    }
}

TEST(Model, AddProcesses)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto mode = representation::SingleMomentMode{ "MODE1", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    process::DissolvedReversibleReaction reaction1{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    process::DissolvedReversibleReaction reaction2{
        reverse_rate,
        forward_rate,
        { hp, ohm },
        { h2o },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    
    EXPECT_EQ(model.processes_.size(), 0);
    
    model.AddProcesses({ reaction1, reaction2 });
    
    EXPECT_EQ(model.processes_.size(), 2);
}

TEST(Model, CollectPhaseStatePrefixesValidation)
{
    // Test that non-unique prefixes throw an error
    auto h2o = micm::Species{ "H2O" };
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o } } };
    
    // Create two representations with the same prefix (this should cause an error)
    auto mode1 = representation::SingleMomentMode{ "DUPLICATE", { aqueous_phase } };
    auto mode2 = representation::SingleMomentMode{ "DUPLICATE", { aqueous_phase } };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode1);
    model.representations_.push_back(mode2);
    
    // Calling SpeciesUsed should throw because of duplicate prefixes
    EXPECT_THROW(model.SpeciesUsed(), std::runtime_error);
}

TEST(Model, SpeciesUsedWithDifferentPhasesInReactions)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto c6h14 = micm::Species{ "C6H14" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp } } };
    auto organic_phase = micm::Phase{ "ORGANIC", { { c6h14 } } };
    
    auto aqueous_mode = representation::SingleMomentMode{ "DROPLET", { aqueous_phase } };
    auto organic_mode = representation::UniformSection{ "PARTICLE", { organic_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 2.0; };
    
    // Reaction only in aqueous phase
    process::DissolvedReversibleReaction aqueous_reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(aqueous_mode);
    model.representations_.push_back(organic_mode);
    model.processes_.push_back(aqueous_reaction);
    
    auto species_used = model.SpeciesUsed();
    
    // Should only have species from aqueous phase (2 species)
    EXPECT_EQ(species_used.size(), 2);
    EXPECT_TRUE(species_used.find("DROPLET.AQUEOUS.H2O") != species_used.end());
    EXPECT_TRUE(species_used.find("DROPLET.AQUEOUS.H+") != species_used.end());
    
    // Organic species should not be present
    EXPECT_TRUE(species_used.find("PARTICLE.ORGANIC.C6H14") == species_used.end());
}

TEST(Model, NonZeroJacobianElementsWithNoProcesses)
{
    auto h2o = micm::Species{ "H2O" };
    auto co2 = micm::Species{ "CO2" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { co2 } } };
    
    auto mode = representation::SingleMomentMode{ "MODE1", { aqueous_phase } };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_indices["MODE1.AQUEOUS.CO2"] = 1;
    
    auto jacobian_elements = model.NonZeroJacobianElements(state_indices);
    
    // Should be empty since no processes
    EXPECT_EQ(jacobian_elements.size(), 0);
}

TEST(Model, NonZeroJacobianElementsWithSingleProcess)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto mode = representation::SingleMomentMode{ "SMALL_DROP", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    // H2O <-> H+ + OH-
    process::DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    model.processes_.push_back(reaction);
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["SMALL_DROP.AQUEOUS.H2O"] = 0;
    state_indices["SMALL_DROP.AQUEOUS.H+"] = 1;
    state_indices["SMALL_DROP.AQUEOUS.OH-"] = 2;
    
    auto jacobian_elements = model.NonZeroJacobianElements(state_indices);
    
    // 9 elements for H2O <-> H+ + OH- reaction
    EXPECT_EQ(jacobian_elements.size(), 9);
    
    // Verify all elements are present
    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            EXPECT_TRUE(jacobian_elements.find({i, j}) != jacobian_elements.end());
        }
    }
}

TEST(Model, NonZeroJacobianElementsWithMultipleProcesses)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    auto co2 = micm::Species{ "CO2" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm }, { co2 }, { h2co3 } } };
    
    auto mode = representation::TwoMomentMode{ "AITKEN", { aqueous_phase } };
    
    auto h2o_forward = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto h2o_reverse = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    process::DissolvedReversibleReaction h2o_dissociation{
        h2o_forward,
        h2o_reverse,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    auto co2_forward = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto co2_reverse = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    process::DissolvedReversibleReaction co2_hydration{
        co2_forward,
        co2_reverse,
        { co2, h2o },
        { h2co3 },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode);
    model.processes_.push_back(h2o_dissociation);
    model.processes_.push_back(co2_hydration);
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["AITKEN.AQUEOUS.H2O"] = 0;
    state_indices["AITKEN.AQUEOUS.H+"] = 1;
    state_indices["AITKEN.AQUEOUS.OH-"] = 2;
    state_indices["AITKEN.AQUEOUS.CO2"] = 3;
    state_indices["AITKEN.AQUEOUS.H2CO3"] = 4;
    
    auto jacobian_elements = model.NonZeroJacobianElements(state_indices);
    
    // Should have contributions from both reactions
    // H2O dissociation affects H2O, H+, OH- (9 elements)
    // CO2 hydration affects CO2, H2O, H2CO3 (9 elements)
    // But there's overlap in H2O row/column, so total < 18
    EXPECT_GT(jacobian_elements.size(), 9);  // More than just one reaction
    EXPECT_LE(jacobian_elements.size(), 18); // Less than or equal to sum of both
    
    // Check some specific elements
    EXPECT_TRUE(jacobian_elements.find({0, 0}) != jacobian_elements.end()); // d[H2O]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({1, 0}) != jacobian_elements.end()); // d[H+]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({3, 0}) != jacobian_elements.end()); // d[CO2]/d[H2O]
    EXPECT_TRUE(jacobian_elements.find({4, 3}) != jacobian_elements.end()); // d[H2CO3]/d[CO2]
}

TEST(Model, NonZeroJacobianElementsWithMultipleRepresentations)
{
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    auto small_drop = representation::SingleMomentMode{ "SMALL_DROP", { aqueous_phase } };
    auto large_drop = representation::SingleMomentMode{ "LARGE_DROP", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-14; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e11; };
    
    process::DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(small_drop);
    model.representations_.push_back(large_drop);
    model.processes_.push_back(reaction);
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["SMALL_DROP.AQUEOUS.H2O"] = 0;
    state_indices["SMALL_DROP.AQUEOUS.H+"] = 1;
    state_indices["SMALL_DROP.AQUEOUS.OH-"] = 2;
    state_indices["LARGE_DROP.AQUEOUS.H2O"] = 3;
    state_indices["LARGE_DROP.AQUEOUS.H+"] = 4;
    state_indices["LARGE_DROP.AQUEOUS.OH-"] = 5;
    
    auto jacobian_elements = model.NonZeroJacobianElements(state_indices);
    
    // 9 elements per representation × 2 representations = 18 elements
    EXPECT_EQ(jacobian_elements.size(), 18);
    
    // Check elements from both representations
    EXPECT_TRUE(jacobian_elements.find({0, 0}) != jacobian_elements.end()); // SMALL_DROP
    EXPECT_TRUE(jacobian_elements.find({1, 2}) != jacobian_elements.end()); // SMALL_DROP
    EXPECT_TRUE(jacobian_elements.find({3, 3}) != jacobian_elements.end()); // LARGE_DROP
    EXPECT_TRUE(jacobian_elements.find({4, 5}) != jacobian_elements.end()); // LARGE_DROP
}

TEST(Model, NonZeroJacobianElementsMixedTypes)
{
    auto h2o = micm::Species{ "H2O" };
    auto co2 = micm::Species{ "CO2" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { co2 }, { h2co3 } } };
    
    auto mode1 = representation::SingleMomentMode{ "MODE1", { aqueous_phase } };
    auto mode2 = representation::TwoMomentMode{ "MODE2", { aqueous_phase } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { return 1.0e-3; };
    auto reverse_rate = [](const micm::Conditions& conditions) { return 1.0e2; };
    
    process::DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { co2, h2o },
        { h2co3 },
        h2o,
        aqueous_phase
    };
    
    Model model;
    model.name_ = "TEST_MODEL";
    model.representations_.push_back(mode1);
    model.representations_.push_back(mode2);
    model.processes_.push_back(reaction);
    
    std::unordered_map<std::string, std::size_t> state_indices;
    state_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_indices["MODE1.AQUEOUS.CO2"] = 1;
    state_indices["MODE1.AQUEOUS.H2CO3"] = 2;
    state_indices["MODE2.AQUEOUS.H2O"] = 3;
    state_indices["MODE2.AQUEOUS.CO2"] = 4;
    state_indices["MODE2.AQUEOUS.H2CO3"] = 5;
    
    auto jacobian_elements = model.NonZeroJacobianElements(state_indices);
    
    // 9 elements per representation × 2 representations = 18 elements
    EXPECT_EQ(jacobian_elements.size(), 18);
    
    // Verify elements exist for both representations
    EXPECT_TRUE(jacobian_elements.find({0, 1}) != jacobian_elements.end()); // MODE1
    EXPECT_TRUE(jacobian_elements.find({2, 0}) != jacobian_elements.end()); // MODE1
    EXPECT_TRUE(jacobian_elements.find({3, 4}) != jacobian_elements.end()); // MODE2
    EXPECT_TRUE(jacobian_elements.find({5, 3}) != jacobian_elements.end()); // MODE2
}
