// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

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

TEST(DissolvedReversibleReaction, UpdateStateParametersFunctionBasic)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    // Create constant rate constants
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    auto update_func = reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);
    
    // Create state parameters matrix (3 grid cells, 2 parameters)
    MatrixPolicy state_parameters(3, 2, 0.0);
    
    // Create conditions for 3 grid cells
    std::vector<micm::Conditions> conditions(3);
    conditions[0].temperature_ = 298.15;
    conditions[1].temperature_ = 310.0;
    conditions[2].temperature_ = 285.0;
    
    // Apply the update function
    update_func(conditions, state_parameters);
    
    // Verify the rate constants were set correctly
    for (std::size_t i_cell = 0; i_cell < 3; ++i_cell)
    {
        EXPECT_EQ(state_parameters[i_cell][0], k_forward);
        EXPECT_EQ(state_parameters[i_cell][1], k_reverse);
    }
}

TEST(DissolvedReversibleReaction, UpdateStateParametersFunctionTemperatureDependent)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    // Create temperature-dependent rate constants (Arrhenius-like)
    auto forward_rate = [](const micm::Conditions& conditions) { 
        return 1.0e-14 * std::exp(-3000.0 / conditions.temperature_); 
    };
    auto reverse_rate = [](const micm::Conditions& conditions) { 
        return 1.0e11 * std::exp(-2000.0 / conditions.temperature_); 
    };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("DROPLET");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    auto update_func = reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);
    
    // Create state parameters matrix (3 grid cells, 2 parameters)
    MatrixPolicy state_parameters(3, 2, 0.0);
    
    // Create conditions with different temperatures
    std::vector<micm::Conditions> conditions(3);
    conditions[0].temperature_ = 298.15;
    conditions[1].temperature_ = 310.0;
    conditions[2].temperature_ = 273.15;
    
    // Apply the update function
    update_func(conditions, state_parameters);
    
    // Verify the rate constants vary with temperature
    double k_f_298 = 1.0e-14 * std::exp(-3000.0 / 298.15);
    double k_f_310 = 1.0e-14 * std::exp(-3000.0 / 310.0);
    double k_f_273 = 1.0e-14 * std::exp(-3000.0 / 273.15);
    
    double k_r_298 = 1.0e11 * std::exp(-2000.0 / 298.15);
    double k_r_310 = 1.0e11 * std::exp(-2000.0 / 310.0);
    double k_r_273 = 1.0e11 * std::exp(-2000.0 / 273.15);
    
    EXPECT_NEAR(state_parameters[0][0], k_f_298, 1e-20);
    EXPECT_NEAR(state_parameters[1][0], k_f_310, 1e-20);
    EXPECT_NEAR(state_parameters[2][0], k_f_273, 1e-20);
    
    EXPECT_NEAR(state_parameters[0][1], k_r_298, 1e5);
    EXPECT_NEAR(state_parameters[1][1], k_r_310, 1e5);
    EXPECT_NEAR(state_parameters[2][1], k_r_273, 1e5);
}

TEST(DissolvedReversibleReaction, UpdateStateParametersFunctionMissingParameter)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
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
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    // Only include forward parameter, not reverse
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    
    // Should throw because reverse parameter is missing
    EXPECT_THROW(
        reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices),
        std::runtime_error
    );
}

TEST(DissolvedReversibleReaction, UpdateStateParametersFunctionMultipleCells)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto co2 = micm::Species{ "CO2" };
    auto h2o = micm::Species{ "H2O" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { co2 }, { h2o }, { h2co3 } } };
    
    auto forward_rate = [](const micm::Conditions& conditions) { 
        return 1.0e-3 * conditions.pressure_ / 101325.0; // Pressure-dependent
    };
    auto reverse_rate = [](const micm::Conditions& conditions) { 
        return 1.0e2;
    };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { co2, h2o },
        { h2co3 },
        h2o,
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("CLOUD");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    auto update_func = reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);
    
    // Create state parameters matrix (5 grid cells, 2 parameters)
    MatrixPolicy state_parameters(5, 2, 0.0);
    
    // Create conditions with different pressures
    std::vector<micm::Conditions> conditions(5);
    for (std::size_t i = 0; i < 5; ++i)
    {
        conditions[i].temperature_ = 298.15;
        conditions[i].pressure_ = 101325.0 * (1.0 + 0.1 * i); // Varying pressure
    }
    
    // Apply the update function
    update_func(conditions, state_parameters);
    
    // Verify the rate constants
    for (std::size_t i_cell = 0; i_cell < 5; ++i_cell)
    {
        double expected_k_f = 1.0e-3 * (1.0 + 0.1 * i_cell);
        EXPECT_NEAR(state_parameters[i_cell][0], expected_k_f, 1e-7);
        EXPECT_EQ(state_parameters[i_cell][1], 1.0e2);
    }
}

// ============================================================================
// ForcingFunction Tests
// ============================================================================

TEST(DissolvedReversibleReaction, ForcingFunctionBasicRates)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
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
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_variable_indices["MODE1.AQUEOUS.H+"] = 1;
    state_variable_indices["MODE1.AQUEOUS.OH-"] = 2;
    
    auto forcing_func = reaction.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices);
    
    // Create state parameters (1 cell, 2 parameters)
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;  // k_forward
    state_parameters[0][1] = k_reverse;  // k_reverse
    
    // Create state variables (1 cell, 3 species)
    MatrixPolicy state_variables(1, 3);
    state_variables[0][0] = 55.0;  // [H2O] = 55 M (solvent)
    state_variables[0][1] = 1.0e-7;  // [H+] = 1e-7 M
    state_variables[0][2] = 1.0e-7;  // [OH-] = 1e-7 M
    
    // Create forcing terms
    MatrixPolicy forcing_terms(1, 3, 0.0);
    
    // Apply the forcing function
    forcing_func(state_parameters, state_variables, forcing_terms);
    
    // Calculate expected rates
    // Forward rate = k_f / [H2O]^(n_reactants-1) * [H2O] = k_f * [H2O] / [H2O]^0 = k_f * [H2O]
    double forward_rate_val = k_forward * 55.0 / std::pow(55.0, 0);  // n_reactants = 1, so pow(solvent, 0) = 1
    // Reverse rate = k_r / [H2O]^(n_products-1) * [H+] * [OH-] = k_r / [H2O] * [H+] * [OH-]
    double reverse_rate_val = k_reverse / 55.0 * 1.0e-7 * 1.0e-7;
    
    // For H2O (reactant): -forward + reverse
    EXPECT_NEAR(forcing_terms[0][0], -forward_rate_val + reverse_rate_val, 1e-20);
    // For H+ (product): +forward - reverse
    EXPECT_NEAR(forcing_terms[0][1], forward_rate_val - reverse_rate_val, 1e-20);
    // For OH- (product): +forward - reverse
    EXPECT_NEAR(forcing_terms[0][2], forward_rate_val - reverse_rate_val, 1e-20);
}

TEST(DissolvedReversibleReaction, ForcingFunctionSolventNormalization)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto co2 = micm::Species{ "CO2" };
    auto h2o = micm::Species{ "H2O" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { co2 }, { h2o }, { h2co3 } } };
    
    double k_forward = 1.0e-3;
    double k_reverse = 1.0e2;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
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
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["DROP.AQUEOUS.CO2"] = 0;
    state_variable_indices["DROP.AQUEOUS.H2O"] = 1;
    state_variable_indices["DROP.AQUEOUS.H2CO3"] = 2;
    
    auto forcing_func = reaction.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 3);
    state_variables[0][0] = 0.001;  // [CO2]
    state_variables[0][1] = 50.0;   // [H2O] (solvent)
    state_variables[0][2] = 0.0005; // [H2CO3]
    
    MatrixPolicy forcing_terms(1, 3, 0.0);
    
    forcing_func(state_parameters, state_variables, forcing_terms);
    
    // Forward rate = k_f / [H2O]^(2-1) * [CO2] * [H2O] = k_f / [H2O] * [CO2] * [H2O]
    double forward_rate_val = k_forward / std::pow(50.0, 1) * 0.001 * 50.0;
    // Reverse rate = k_r / [H2O]^(1-1) * [H2CO3] = k_r * [H2CO3]
    double reverse_rate_val = k_reverse / std::pow(50.0, 0) * 0.0005;
    
    // For CO2 (reactant): -forward + reverse
    EXPECT_NEAR(forcing_terms[0][0], -forward_rate_val + reverse_rate_val, 1e-10);
    // For H2O (reactant): -forward + reverse
    EXPECT_NEAR(forcing_terms[0][1], -forward_rate_val + reverse_rate_val, 1e-10);
    // For H2CO3 (product): +forward - reverse
    EXPECT_NEAR(forcing_terms[0][2], forward_rate_val - reverse_rate_val, 1e-10);
}

TEST(DissolvedReversibleReaction, ForcingFunctionMultipleReactantsProducts)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto a = micm::Species{ "A" };
    auto b = micm::Species{ "B" };
    auto c = micm::Species{ "C" };
    auto d = micm::Species{ "D" };
    auto solvent = micm::Species{ "SOLVENT" };
    
    auto phase = micm::Phase{ "LIQUID", { { a }, { b }, { c }, { d }, { solvent } } };
    
    double k_forward = 2.0;
    double k_reverse = 3.0;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    // A + B <-> C + D
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
    
    std::string forward_param = phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["REP1.LIQUID.A"] = 0;
    state_variable_indices["REP1.LIQUID.B"] = 1;
    state_variable_indices["REP1.LIQUID.C"] = 2;
    state_variable_indices["REP1.LIQUID.D"] = 3;
    state_variable_indices["REP1.LIQUID.SOLVENT"] = 4;
    
    auto forcing_func = reaction.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 5);
    state_variables[0][0] = 2.0;   // [A]
    state_variables[0][1] = 3.0;   // [B]
    state_variables[0][2] = 1.0;   // [C]
    state_variables[0][3] = 0.5;   // [D]
    state_variables[0][4] = 40.0;  // [SOLVENT]
    
    MatrixPolicy forcing_terms(1, 5, 0.0);
    
    forcing_func(state_parameters, state_variables, forcing_terms);
    
    // Forward rate = k_f / [solvent]^(2-1) * [A] * [B] = k_f / [solvent] * [A] * [B]
    double forward_rate_val = k_forward / 40.0 * 2.0 * 3.0;
    // Reverse rate = k_r / [solvent]^(2-1) * [C] * [D] = k_r / [solvent] * [C] * [D]
    double reverse_rate_val = k_reverse / 40.0 * 1.0 * 0.5;
    
    // Reactants: -forward + reverse
    EXPECT_NEAR(forcing_terms[0][0], -forward_rate_val + reverse_rate_val, 1e-10);  // A
    EXPECT_NEAR(forcing_terms[0][1], -forward_rate_val + reverse_rate_val, 1e-10);  // B
    // Products: +forward - reverse
    EXPECT_NEAR(forcing_terms[0][2], forward_rate_val - reverse_rate_val, 1e-10);  // C
    EXPECT_NEAR(forcing_terms[0][3], forward_rate_val - reverse_rate_val, 1e-10);  // D
    // Solvent is not updated
    EXPECT_EQ(forcing_terms[0][4], 0.0);
}

TEST(DissolvedReversibleReaction, ForcingFunctionMultipleCells)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_variable_indices["MODE1.AQUEOUS.H+"] = 1;
    state_variable_indices["MODE1.AQUEOUS.OH-"] = 2;
    
    auto forcing_func = reaction.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices);
    
    // Test with 3 cells
    MatrixPolicy state_parameters(3, 2);
    MatrixPolicy state_variables(3, 3);
    MatrixPolicy forcing_terms(3, 3, 0.0);
    
    // Different conditions for each cell
    for (std::size_t i = 0; i < 3; ++i)
    {
        state_parameters[i][0] = k_forward;
        state_parameters[i][1] = k_reverse;
        state_variables[i][0] = 55.0 + i * 5.0;  // Different [H2O]
        state_variables[i][1] = 1.0e-7 * (i + 1);  // Different [H+]
        state_variables[i][2] = 1.0e-7 * (i + 1);  // Different [OH-]
    }
    
    forcing_func(state_parameters, state_variables, forcing_terms);
    
    // Verify each cell independently
    for (std::size_t i = 0; i < 3; ++i)
    {
        double h2o_conc = 55.0 + i * 5.0;
        double hp_conc = 1.0e-7 * (i + 1);
        double ohm_conc = 1.0e-7 * (i + 1);
        
        double forward_rate_val = k_forward * h2o_conc;
        double reverse_rate_val = k_reverse / h2o_conc * hp_conc * ohm_conc;
        
        EXPECT_NEAR(forcing_terms[i][0], -forward_rate_val + reverse_rate_val, 1e-20);
        EXPECT_NEAR(forcing_terms[i][1], forward_rate_val - reverse_rate_val, 1e-20);
        EXPECT_NEAR(forcing_terms[i][2], forward_rate_val - reverse_rate_val, 1e-20);
    }
}

TEST(DissolvedReversibleReaction, ForcingFunctionMultiplePhaseInstances)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    // Two phase instances (e.g., small and large drops)
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
    phase_prefixes["AQUEOUS"].insert("LARGE_DROP");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["SMALL_DROP.AQUEOUS.H2O"] = 0;
    state_variable_indices["SMALL_DROP.AQUEOUS.H+"] = 1;
    state_variable_indices["SMALL_DROP.AQUEOUS.OH-"] = 2;
    state_variable_indices["LARGE_DROP.AQUEOUS.H2O"] = 3;
    state_variable_indices["LARGE_DROP.AQUEOUS.H+"] = 4;
    state_variable_indices["LARGE_DROP.AQUEOUS.OH-"] = 5;
    
    auto forcing_func = reaction.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 6);
    // Small drop
    state_variables[0][0] = 50.0;  // [H2O]
    state_variables[0][1] = 1.0e-7;  // [H+]
    state_variables[0][2] = 1.0e-7;  // [OH-]
    // Large drop
    state_variables[0][3] = 60.0;  // [H2O]
    state_variables[0][4] = 2.0e-7;  // [H+]
    state_variables[0][5] = 2.0e-7;  // [OH-]
    
    MatrixPolicy forcing_terms(1, 6, 0.0);
    
    forcing_func(state_parameters, state_variables, forcing_terms);
    
    // Verify small drop
    double forward_small = k_forward * 50.0;
    double reverse_small = k_reverse / 50.0 * 1.0e-7 * 1.0e-7;
    EXPECT_NEAR(forcing_terms[0][0], -forward_small + reverse_small, 1e-20);
    EXPECT_NEAR(forcing_terms[0][1], forward_small - reverse_small, 1e-20);
    EXPECT_NEAR(forcing_terms[0][2], forward_small - reverse_small, 1e-20);
    
    // Verify large drop
    double forward_large = k_forward * 60.0;
    double reverse_large = k_reverse / 60.0 * 2.0e-7 * 2.0e-7;
    EXPECT_NEAR(forcing_terms[0][3], -forward_large + reverse_large, 1e-20);
    EXPECT_NEAR(forcing_terms[0][4], forward_large - reverse_large, 1e-20);
    EXPECT_NEAR(forcing_terms[0][5], forward_large - reverse_large, 1e-20);
}

// ============================================================================
// JacobianFunction Tests
// ============================================================================

TEST(DissolvedReversibleReaction, JacobianFunctionBasicPartials)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    // H2O <-> H+ + OH-
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },      // 1 reactant
        { hp, ohm },  // 2 products
        h2o,          // solvent
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_variable_indices["MODE1.AQUEOUS.H+"] = 1;
    state_variable_indices["MODE1.AQUEOUS.OH-"] = 2;
    
    // Build sparse jacobian structure
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 3);
    state_variables[0][0] = 55.0;  // [H2O]
    state_variables[0][1] = 1.0e-7;  // [H+]
    state_variables[0][2] = 1.0e-7;  // [OH-]
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // Expected partial derivatives:
    // d[H2O]/d[H2O] = -k_f + k_r * [H+] * [OH-] * (1-1) / [H2O]^1 = -k_f - k_r * [H+] * [OH-] / [H2O]^2
    // d[H2O]/d[H+] = k_r / [H2O] * [OH-]
    // d[H2O]/d[OH-] = k_r / [H2O] * [H+]
    // d[H+]/d[H2O] = k_f - k_r * [H+] * [OH-] * (1-2) / [H2O]^2 = k_f + k_r * [H+] * [OH-] / [H2O]^2
    // d[H+]/d[H+] = -k_r / [H2O] * [OH-]
    // d[H+]/d[OH-] = -k_r / [H2O] * [H+]
    // d[OH-]/d[H2O] = k_f - k_r * [H+] * [OH-] * (1-2) / [H2O]^2 = k_f + k_r * [H+] * [OH-] / [H2O]^2
    // d[OH-]/d[H+] = -k_r / [H2O] * [OH-]
    // d[OH-]/d[OH-] = -k_r / [H2O] * [H+]
    
    double d_forward_d_h2o = k_forward;
    double d_reverse_d_hp = k_reverse / 55.0 * 1.0e-7;
    double d_reverse_d_ohm = k_reverse / 55.0 * 1.0e-7;
    double d_forward_d_solvent = k_forward * 55.0 * (1 - 1) / std::pow(55.0, 1);  // = 0
    // d_reverse_d_solvent is the partial derivative of the reverse rate w.r.t. solvent
    // reverse_rate = k_r / [H2O]^1 * [H+] * [OH-], so d(reverse_rate)/d[H2O] = -k_r / [H2O]^2 * [H+] * [OH-]
    double d_reverse_d_solvent = k_reverse * (1 - 2) / std::pow(55.0, 2) * 1.0e-7 * 1.0e-7;
    
    // For H2O (reactant): d[H2O]/dt = -forward_rate + reverse_rate
    // d[H2O]/d[H2O] = -d(forward_rate)/d[H2O] + d(reverse_rate)/d[H2O]
    EXPECT_NEAR(jacobian[0][0][0], -d_forward_d_h2o + d_reverse_d_solvent, 1e-15);  // d[H2O]/d[H2O]
    EXPECT_NEAR(jacobian[0][0][1], d_reverse_d_ohm, 1e-20);  // d[H2O]/d[H+]
    EXPECT_NEAR(jacobian[0][0][2], d_reverse_d_hp, 1e-20);   // d[H2O]/d[OH-]
    // For H+ and OH- (products): d[H+]/dt = +forward_rate - reverse_rate
    // d[H+]/d[H2O] = +d(forward_rate)/d[H2O] - d(reverse_rate)/d[H2O]
    EXPECT_NEAR(jacobian[0][1][0], d_forward_d_h2o - d_reverse_d_solvent, 1e-15);  // d[H+]/d[H2O]
    EXPECT_NEAR(jacobian[0][1][1], -d_reverse_d_ohm, 1e-20);  // d[H+]/d[H+]
    EXPECT_NEAR(jacobian[0][1][2], -d_reverse_d_hp, 1e-20);   // d[H+]/d[OH-]
    EXPECT_NEAR(jacobian[0][2][0], d_forward_d_h2o - d_reverse_d_solvent, 1e-15);  // d[OH-]/d[H2O]
    EXPECT_NEAR(jacobian[0][2][1], -d_reverse_d_ohm, 1e-20);  // d[OH-]/d[H+]
    EXPECT_NEAR(jacobian[0][2][2], -d_reverse_d_hp, 1e-20);   // d[OH-]/d[OH-]
}

TEST(DissolvedReversibleReaction, JacobianFunctionMultipleReactantsProducts)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    auto co2 = micm::Species{ "CO2" };
    auto h2o = micm::Species{ "H2O" };
    auto h2co3 = micm::Species{ "H2CO3" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { co2 }, { h2o }, { h2co3 } } };
    
    double k_forward = 1.0e-3;
    double k_reverse = 1.0e2;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
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
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["DROP.AQUEOUS.CO2"] = 0;
    state_variable_indices["DROP.AQUEOUS.H2O"] = 1;
    state_variable_indices["DROP.AQUEOUS.H2CO3"] = 2;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 3);
    state_variables[0][0] = 0.001;  // [CO2]
    state_variables[0][1] = 50.0;   // [H2O]
    state_variables[0][2] = 0.0005; // [H2CO3]
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // Expected partials for CO2 + H2O <-> H2CO3:
    // Forward rate = k_f / [H2O] * [CO2] * [H2O]
    // Reverse rate = k_r * [H2CO3]
    
    // d[CO2]/d[CO2] = -k_f / [H2O] * [H2O] = -k_f
    double d_co2_d_co2 = -k_forward / 50.0 * 50.0;
    // d[CO2]/d[H2O] = -k_f / [H2O] * [CO2] + k_f * (1-2) / [H2O]^2 * [CO2] * [H2O] = -k_f * [CO2] / [H2O] - k_f * [CO2] / [H2O]
    double d_forward_d_h2o = k_forward / 50.0 * 0.001;
    double d_forward_d_solvent = k_forward * (1 - 2) / std::pow(50.0, 2) * 0.001 * 50.0;
    double d_co2_d_h2o = -d_forward_d_h2o - d_forward_d_solvent;
    // d[CO2]/d[H2CO3] = k_r
    double d_co2_d_h2co3 = k_reverse;
    
    // d[H2O]/d[CO2] = -k_f / [H2O] * [H2O] = -k_f
    double d_h2o_d_co2 = d_co2_d_co2;
    // d[H2O]/d[H2O] = same as d[CO2]/d[H2O]
    double d_h2o_d_h2o = d_co2_d_h2o;
    // d[H2O]/d[H2CO3] = k_r
    double d_h2o_d_h2co3 = d_co2_d_h2co3;
    
    // d[H2CO3]/d[CO2] = k_f / [H2O] * [H2O] = k_f
    double d_h2co3_d_co2 = -d_co2_d_co2;
    // d[H2CO3]/d[H2O] = -d[CO2]/d[H2O]
    double d_h2co3_d_h2o = -d_co2_d_h2o;
    // d[H2CO3]/d[H2CO3] = -k_r
    double d_h2co3_d_h2co3 = -k_reverse;
    
    EXPECT_NEAR(jacobian[0][0][0], d_co2_d_co2, 1e-10);
    EXPECT_NEAR(jacobian[0][0][1], d_co2_d_h2o, 1e-10);
    EXPECT_NEAR(jacobian[0][0][2], d_co2_d_h2co3, 1e-10);
    EXPECT_NEAR(jacobian[0][1][0], d_h2o_d_co2, 1e-10);
    EXPECT_NEAR(jacobian[0][1][1], d_h2o_d_h2o, 1e-10);
    EXPECT_NEAR(jacobian[0][1][2], d_h2o_d_h2co3, 1e-10);
    EXPECT_NEAR(jacobian[0][2][0], d_h2co3_d_co2, 1e-10);
    EXPECT_NEAR(jacobian[0][2][1], d_h2co3_d_h2o, 1e-10);
    EXPECT_NEAR(jacobian[0][2][2], d_h2co3_d_h2co3, 1e-10);
}

TEST(DissolvedReversibleReaction, JacobianFunctionSigns)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    auto hco3m = micm::Species{ "HCO3-" };
    auto hp = micm::Species{ "H+" };
    auto co32m = micm::Species{ "CO32-" };
    auto h2o = micm::Species{ "H2O" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { hco3m }, { hp }, { co32m }, { h2o } } };
    
    double k_forward = 1.0e-10;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    // HCO3- <-> H+ + CO32-
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { hco3m },      // 1 reactant
        { hp, co32m },  // 2 products
        h2o,            // solvent
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE.AQUEOUS.HCO3-"] = 0;
    state_variable_indices["MODE.AQUEOUS.H+"] = 1;
    state_variable_indices["MODE.AQUEOUS.CO32-"] = 2;
    state_variable_indices["MODE.AQUEOUS.H2O"] = 3;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    auto jacobian_builder = SparseMatrixPolicy::Create(4).SetNumberOfBlocks(1);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 4);
    state_variables[0][0] = 0.001;  // [HCO3-]
    state_variables[0][1] = 1.0e-7;  // [H+]
    state_variables[0][2] = 5.0e-8;  // [CO32-]
    state_variables[0][3] = 55.0;    // [H2O]
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // Check signs: 
    // - Reactant w.r.t. reactant: negative (forward dominates)
    // - Reactant w.r.t. product: positive (reverse adds back)
    // - Product w.r.t. reactant: positive (forward produces)
    // - Product w.r.t. product: negative (reverse consumes)
    
    // d[HCO3-]/d[HCO3-]: negative
    EXPECT_LT(jacobian[0][0][0], 0.0);
    // d[HCO3-]/d[H+]: positive
    EXPECT_GT(jacobian[0][0][1], 0.0);
    // d[HCO3-]/d[CO32-]: positive
    EXPECT_GT(jacobian[0][0][2], 0.0);
    
    // d[H+]/d[HCO3-]: positive
    EXPECT_GT(jacobian[0][1][0], 0.0);
    // d[H+]/d[H+]: negative
    EXPECT_LT(jacobian[0][1][1], 0.0);
    // d[H+]/d[CO32-]: negative
    EXPECT_LT(jacobian[0][1][2], 0.0);
    
    // d[CO32-]/d[HCO3-]: positive
    EXPECT_GT(jacobian[0][2][0], 0.0);
    // d[CO32-]/d[H+]: negative
    EXPECT_LT(jacobian[0][2][1], 0.0);
    // d[CO32-]/d[CO32-]: negative
    EXPECT_LT(jacobian[0][2][2], 0.0);
}

TEST(DissolvedReversibleReaction, JacobianFunctionMultipleCells)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE1.AQUEOUS.H2O"] = 0;
    state_variable_indices["MODE1.AQUEOUS.H+"] = 1;
    state_variable_indices["MODE1.AQUEOUS.OH-"] = 2;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(3);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    // Test with 3 cells
    MatrixPolicy state_parameters(3, 2);
    MatrixPolicy state_variables(3, 3);
    
    for (std::size_t i = 0; i < 3; ++i)
    {
        state_parameters[i][0] = k_forward;
        state_parameters[i][1] = k_reverse;
        state_variables[i][0] = 55.0 + i;  // Different [H2O]
        state_variables[i][1] = 1.0e-7 * (i + 1);  // Different [H+]
        state_variables[i][2] = 1.0e-7 * (i + 1);  // Different [OH-]
    }
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // Verify each cell independently
    for (std::size_t i = 0; i < 3; ++i)
    {
        double h2o_conc = 55.0 + i;
        double hp_conc = 1.0e-7 * (i + 1);
        double ohm_conc = 1.0e-7 * (i + 1);
        
        // d[H2O]/d[H2O] should be negative (forward dominates over reverse for typical water concentrations)
        EXPECT_LT(jacobian[i][0][0], 0.0);
        
        // d[H+]/d[H2O] should be positive (forward produces H+)
        EXPECT_GT(jacobian[i][1][0], 0.0);
        
        // d[OH-]/d[H2O] should be positive (forward produces OH-)
        EXPECT_GT(jacobian[i][2][0], 0.0);
    }
}

TEST(DissolvedReversibleReaction, JacobianFunctionMultiplePhaseInstances)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    auto h2o = micm::Species{ "H2O" };
    auto hp = micm::Species{ "H+" };
    auto ohm = micm::Species{ "OH-" };
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { hp }, { ohm } } };
    
    double k_forward = 1.0e-14;
    double k_reverse = 1.0e11;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { h2o },
        { hp, ohm },
        h2o,
        aqueous_phase
    };
    
    // Two phase instances
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
    phase_prefixes["AQUEOUS"].insert("LARGE_DROP");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["SMALL_DROP.AQUEOUS.H2O"] = 0;
    state_variable_indices["SMALL_DROP.AQUEOUS.H+"] = 1;
    state_variable_indices["SMALL_DROP.AQUEOUS.OH-"] = 2;
    state_variable_indices["LARGE_DROP.AQUEOUS.H2O"] = 3;
    state_variable_indices["LARGE_DROP.AQUEOUS.H+"] = 4;
    state_variable_indices["LARGE_DROP.AQUEOUS.OH-"] = 5;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    auto jacobian_builder = SparseMatrixPolicy::Create(6).SetNumberOfBlocks(1);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 6);
    // Small drop (indices 0-2)
    state_variables[0][0] = 50.0;  // [H2O]
    state_variables[0][1] = 1.0e-7;  // [H+]
    state_variables[0][2] = 1.0e-7;  // [OH-]
    // Large drop (indices 3-5)
    state_variables[0][3] = 60.0;  // [H2O]
    state_variables[0][4] = 2.0e-7;  // [H+]
    state_variables[0][5] = 2.0e-7;  // [OH-]
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // Verify small drop elements are populated
    EXPECT_NE(jacobian[0][0][0], 0.0);  // d[H2O_small]/d[H2O_small]
    EXPECT_NE(jacobian[0][1][0], 0.0);  // d[H+_small]/d[H2O_small]
    EXPECT_NE(jacobian[0][2][0], 0.0);  // d[OH-_small]/d[H2O_small]
    
    // Verify large drop elements are populated
    EXPECT_NE(jacobian[0][3][3], 0.0);  // d[H2O_large]/d[H2O_large]
    EXPECT_NE(jacobian[0][4][3], 0.0);  // d[H+_large]/d[H2O_large]
    EXPECT_NE(jacobian[0][5][3], 0.0);  // d[OH-_large]/d[H2O_large]
}

TEST(DissolvedReversibleReaction, JacobianFunctionSimpleDistinctSpecies)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    // Use completely distinct species to avoid any overlap
    auto foo = micm::Species{ "foo" };  // reactant
    auto bar = micm::Species{ "bar" };  // product
    auto baz = micm::Species{ "baz" };  // solvent
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { foo }, { bar }, { baz } } };
    
    double k_forward = 2.0;
    double k_reverse = 3.0;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    // foo <-> bar (with solvent baz)
    // With 1 reactant and 1 product, there's no solvent concentration dependence
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { foo },      // 1 reactant
        { bar },      // 1 product
        baz,          // solvent
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE1.AQUEOUS.foo"] = 0;
    state_variable_indices["MODE1.AQUEOUS.bar"] = 1;
    state_variable_indices["MODE1.AQUEOUS.baz"] = 2;
    
    // Build sparse jacobian structure
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    EXPECT_EQ(jacobian_elements.size(), 6);  // 2 species x 3 dependencies (self, other species, solvent)
    
    auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 3);
    state_variables[0][0] = 10.0;  // [foo]
    state_variables[0][1] = 5.0;   // [bar]
    state_variables[0][2] = 50.0;  // [baz]
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // For foo <-> bar:
    // Forward rate = k_f / [baz]^(N_reactants - 1) * [foo] = k_f / [baz]^0 * [foo] = k_f * [foo]
    // Reverse rate = k_r / [baz]^(N_products - 1) * [bar] = k_r / [baz]^0 * [bar] = k_r * [bar]
    // 
    // d[foo]/dt = -k_f * [foo] + k_r * [bar]
    // d[bar]/dt = +k_f * [foo] - k_r * [bar]
    //
    // Jacobian elements:
    // d[foo]/d[foo] = -k_f = -2.0
    // d[foo]/d[bar] = +k_r = +3.0
    // d[foo]/d[baz] = 0 (no solvent dependence when N_reactants = N_products = 1)
    //
    // d[bar]/d[foo] = +k_f = +2.0
    // d[bar]/d[bar] = -k_r = -3.0
    // d[bar]/d[baz] = 0 (no solvent dependence when N_reactants = N_products = 1)
    
    double d_forward_d_foo = k_forward;
    double d_reverse_d_bar = k_reverse;
    double d_foo_d_baz = 0.0;  // No solvent dependence for 1 reactant / 1 product
    double d_bar_d_baz = 0.0;  // No solvent dependence for 1 reactant / 1 product
    
    EXPECT_NEAR(jacobian[0][0][0], -d_forward_d_foo, 1e-10);  // d[foo]/d[foo]
    EXPECT_NEAR(jacobian[0][0][1], d_reverse_d_bar, 1e-10);    // d[foo]/d[bar]
    EXPECT_NEAR(jacobian[0][0][2], d_foo_d_baz, 1e-10);        // d[foo]/d[baz]
    EXPECT_NEAR(jacobian[0][1][0], d_forward_d_foo, 1e-10);    // d[bar]/d[foo]
    EXPECT_NEAR(jacobian[0][1][1], -d_reverse_d_bar, 1e-10);   // d[bar]/d[bar]
    EXPECT_NEAR(jacobian[0][1][2], d_bar_d_baz, 1e-10);        // d[bar]/d[baz]
}

TEST(DissolvedReversibleReaction, JacobianFunctionTwoReactantsWithSolventDependence)
{
    using MatrixPolicy = micm::Matrix<double>;
    using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
    
    auto foo = micm::Species{ "foo" };   // reactant 1
    auto qux = micm::Species{ "qux" };   // reactant 2
    auto bar = micm::Species{ "bar" };   // product
    auto baz = micm::Species{ "baz" };   // solvent
    
    auto aqueous_phase = micm::Phase{ "AQUEOUS", { { foo }, { qux }, { bar }, { baz } } };
    
    double k_forward = 2.0;
    double k_reverse = 3.0;
    
    auto forward_rate = [k_forward](const micm::Conditions& conditions) { return k_forward; };
    auto reverse_rate = [k_reverse](const micm::Conditions& conditions) { return k_reverse; };
    
    // foo + qux <-> bar (with solvent baz)
    // With 2 reactants, forward rate has solvent dependence: k_f / [baz]^1
    DissolvedReversibleReaction reaction{
        forward_rate,
        reverse_rate,
        { foo, qux },  // 2 reactants
        { bar },       // 1 product
        baz,           // solvent
        aqueous_phase
    };
    
    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("MODE1");
    
    std::string forward_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_forward";
    std::string reverse_param = aqueous_phase.name_ + "." + reaction.uuid_ + ".k_reverse";
    
    std::unordered_map<std::string, std::size_t> state_parameter_indices;
    state_parameter_indices[forward_param] = 0;
    state_parameter_indices[reverse_param] = 1;
    
    std::unordered_map<std::string, std::size_t> state_variable_indices;
    state_variable_indices["MODE1.AQUEOUS.foo"] = 0;
    state_variable_indices["MODE1.AQUEOUS.qux"] = 1;
    state_variable_indices["MODE1.AQUEOUS.bar"] = 2;
    state_variable_indices["MODE1.AQUEOUS.baz"] = 3;
    
    auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    EXPECT_EQ(jacobian_elements.size(), 12);  // 3 species x 4 dependencies (2 reactants, 1 product, 1 solvent)
    
    auto jacobian_builder = SparseMatrixPolicy::Create(4).SetNumberOfBlocks(1);
    for (const auto& elem : jacobian_elements)
    {
        jacobian_builder.WithElement(elem.first, elem.second);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);
    
    auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    
    MatrixPolicy state_parameters(1, 2);
    state_parameters[0][0] = k_forward;
    state_parameters[0][1] = k_reverse;
    
    MatrixPolicy state_variables(1, 4);
    state_variables[0][0] = 10.0;  // [foo]
    state_variables[0][1] = 8.0;   // [qux]
    state_variables[0][2] = 5.0;   // [bar]
    state_variables[0][3] = 50.0;  // [baz]
    
    jacobian_func(state_parameters, state_variables, jacobian);
    
    // Forward rate = k_f / [baz]^(2-1) * [foo] * [qux] = k_f / [baz] * [foo] * [qux]
    // Reverse rate = k_r / [baz]^(1-1) * [bar] = k_r * [bar]
    //
    // d[foo]/dt = -(k_f / [baz]) * [foo] * [qux] + k_r * [bar]
    // d[qux]/dt = -(k_f / [baz]) * [foo] * [qux] + k_r * [bar]
    // d[bar]/dt = +(k_f / [baz]) * [foo] * [qux] - k_r * [bar]
    //
    // Jacobian elements involving [baz]:
    // d(forward_rate)/d[baz] = k_f * [foo] * [qux] * d/d[baz]([baz]^-1) 
    //                        = k_f * [foo] * [qux] * (-1) / [baz]^2
    //                        = -k_f * [foo] * [qux] / [baz]^2
    //                        = -2.0 * 10.0 * 8.0 / 2500 = -0.064
    //
    // d[foo]/d[baz] = -d(forward_rate)/d[baz] = -(-0.064) = +0.064    // (because of the minus sign in the rate equation)
    // d[qux]/d[baz] = -d(forward_rate)/d[baz] = +0.064
    // d[bar]/d[baz] = +d(forward_rate)/d[baz] =  -0.064
    
    double forward_rate_val = k_forward / state_variables[0][3] * state_variables[0][0] * state_variables[0][1];
    double d_forward_rate_d_baz = -k_forward * state_variables[0][0] * state_variables[0][1] / std::pow(state_variables[0][3], 2);
    double d_foo_d_baz = -d_forward_rate_d_baz;   // Opposite sign because foo is consumed
    double d_qux_d_baz = -d_forward_rate_d_baz;   // Opposite sign because qux is consumed
    double d_bar_d_baz = d_forward_rate_d_baz;     // Same sign because bar is produced
    
    EXPECT_NEAR(jacobian[0][0][3], d_foo_d_baz, 1e-10);  // d[foo]/d[baz]
    EXPECT_NEAR(jacobian[0][1][3], d_qux_d_baz, 1e-10);  // d[qux]/d[baz]
    EXPECT_NEAR(jacobian[0][2][3], d_bar_d_baz, 1e-10);  // d[bar]/d[baz]
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
