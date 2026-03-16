// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../representation_policy.hpp"

#include <miam/representations/uniform_section.hpp>

using namespace miam::representation;

TEST(UniformSection, StateSize)
{
    auto size = testStateSize<UniformSection>();
    // Should have 3 species and 2 parameters (MIN_RADIUS, MAX_RADIUS)
    EXPECT_EQ(std::get<0>(size), 3);
    EXPECT_EQ(std::get<1>(size), 2);
}

TEST(UniformSection, StateVariableNames)
{
    auto names = testStateVariableNames<UniformSection>();
    // Should have 3 species with phase-qualified names
    EXPECT_EQ(names.size(), 3);
    EXPECT_TRUE(names.find(test_model_name + ".PHASE1.SPECIES_A") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".PHASE1.SPECIES_B") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".PHASE2.SPECIES_C") != names.end());
}

TEST(UniformSection, StateParameterNames)
{
    auto names = testStateParameterNames<UniformSection>();
    // Should have 2 parameters: MIN_RADIUS and MAX_RADIUS
    EXPECT_EQ(names.size(), 2);
    EXPECT_TRUE(names.find(test_model_name + ".MIN_RADIUS") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".MAX_RADIUS") != names.end());
}

TEST(UniformSection, DefaultParametersPolicy)
{
    auto params = testDefaultParameters<UniformSection>();
    
    // Check default values (0.0 for both when using default constructor)
    EXPECT_EQ(params[test_model_name + ".MIN_RADIUS"], 0.0);
    EXPECT_EQ(params[test_model_name + ".MAX_RADIUS"], 0.0);
}

TEST(UniformSection, DefaultParameters)
{
    const auto phases = getTestPhases();
    auto model = UniformSection{ test_model_name, phases };
    auto params = model.DefaultParameters();
    
    // Should have 2 default parameters
    EXPECT_EQ(params.size(), 2);
    
    // Check parameter names exist
    EXPECT_TRUE(params.find(test_model_name + ".MIN_RADIUS") != params.end());
    EXPECT_TRUE(params.find(test_model_name + ".MAX_RADIUS") != params.end());
    
    // Check default values (0.0 for both when using default constructor)
    EXPECT_EQ(params[test_model_name + ".MIN_RADIUS"], 0.0);
    EXPECT_EQ(params[test_model_name + ".MAX_RADIUS"], 0.0);
}

TEST(UniformSection, DefaultParametersWithCustomValues)
{
    const auto phases = getTestPhases();
    const double custom_min_radius = 1.0e-7;
    const double custom_max_radius = 1.0e-6;
    auto model = UniformSection{ test_model_name, phases, custom_min_radius, custom_max_radius };
    auto params = model.DefaultParameters();
    
    // Check custom default values
    EXPECT_EQ(params[test_model_name + ".MIN_RADIUS"], custom_min_radius);
    EXPECT_EQ(params[test_model_name + ".MAX_RADIUS"], custom_max_radius);
}

TEST(UniformSection, MinRadiusParameterName)
{
    const auto phases = getTestPhases();
    auto model = UniformSection{ test_model_name, phases };
    
    EXPECT_EQ(model.MinRadius(), test_model_name + ".MIN_RADIUS");
}

TEST(UniformSection, MaxRadiusParameterName)
{
    const auto phases = getTestPhases();
    auto model = UniformSection{ test_model_name, phases };
    
    EXPECT_EQ(model.MaxRadius(), test_model_name + ".MAX_RADIUS");
}

TEST(UniformSection, SpeciesNaming)
{
    micm::Species test_species{ "CO2" };
    micm::Phase test_phase{ "AQUEOUS", { { test_species } } };
    std::vector<micm::Phase> phases = { test_phase };
    
    auto model = UniformSection{ test_model_name, phases };
    
    std::string species_name = model.Species(test_phase, test_species);
    EXPECT_EQ(species_name, test_model_name + ".AQUEOUS.CO2");
}

TEST(UniformSection, MultiplePhases)
{
    micm::Phase aqueous{ "AQUEOUS", { { micm::Species{ "H2O" } }, { micm::Species{ "CO2" } } } };
    micm::Phase organic{ "ORGANIC", { { micm::Species{ "C6H14" } } } };
    std::vector<micm::Phase> phases = { aqueous, organic };
    
    auto model = UniformSection{ "SECTION1", phases };
    
    auto size = model.StateSize();
    // 3 species (2 from aqueous + 1 from organic) and 2 parameters
    EXPECT_EQ(std::get<0>(size), 3);
    EXPECT_EQ(std::get<1>(size), 2);
    
    auto var_names = model.StateVariableNames();
    EXPECT_EQ(var_names.size(), 3);
    // Names should include phase: SECTION1.AQUEOUS.H2O, SECTION1.AQUEOUS.CO2, SECTION1.ORGANIC.C6H14
    EXPECT_TRUE(var_names.find("SECTION1.AQUEOUS.H2O") != var_names.end());
    EXPECT_TRUE(var_names.find("SECTION1.AQUEOUS.CO2") != var_names.end());
    EXPECT_TRUE(var_names.find("SECTION1.ORGANIC.C6H14") != var_names.end());
}

TEST(UniformSection, PrefixConsistency)
{
    const auto phases = getTestPhases();
    const std::string custom_prefix = "CUSTOM_SECTION";
    auto model = UniformSection{ custom_prefix, phases };
    
    // Check that all method outputs use the custom prefix
    EXPECT_EQ(model.MinRadius(), custom_prefix + ".MIN_RADIUS");
    EXPECT_EQ(model.MaxRadius(), custom_prefix + ".MAX_RADIUS");
    
    auto var_names = model.StateVariableNames();
    for (const auto& name : var_names)
    {
        EXPECT_EQ(name.substr(0, custom_prefix.length()), custom_prefix);
    }
    
    auto param_names = model.StateParameterNames();
    for (const auto& name : param_names)
    {
        EXPECT_EQ(name.substr(0, custom_prefix.length()), custom_prefix);
    }
}

TEST(UniformSection, RadiusRangeValidation)
{
    const auto phases = getTestPhases();
    const double min_radius = 1.0e-7;
    const double max_radius = 1.0e-5;
    auto model = UniformSection{ test_model_name, phases, min_radius, max_radius };
    
    auto params = model.DefaultParameters();
    
    // Ensure both radius parameters are present and have correct values
    EXPECT_EQ(params[test_model_name + ".MIN_RADIUS"], min_radius);
    EXPECT_EQ(params[test_model_name + ".MAX_RADIUS"], max_radius);
    
    // In a physical system, max should be greater than min
    EXPECT_GT(params[test_model_name + ".MAX_RADIUS"], params[test_model_name + ".MIN_RADIUS"]);
}

TEST(UniformSection, NumPhaseInstances)
{
    const auto phases = getTestPhases();
    auto model = UniformSection{ test_model_name, phases };
    
    auto num_instances = model.NumPhaseInstances();
    
    // Should have one instance per phase (2 phases in test setup)
    EXPECT_EQ(num_instances.size(), 2);
    EXPECT_EQ(num_instances["PHASE1"], 1);
    EXPECT_EQ(num_instances["PHASE2"], 1);
}

TEST(UniformSection, PhaseStatePrefixes)
{
    const auto phases = getTestPhases();
    const std::string prefix = "SECTION1";
    auto model = UniformSection{ prefix, phases };
    
    auto phase_prefixes = model.PhaseStatePrefixes();
    
    // Should have one prefix per phase
    EXPECT_EQ(phase_prefixes.size(), 2);
    
    // Each phase should have exactly one prefix (the model's prefix)
    EXPECT_EQ(phase_prefixes["PHASE1"].size(), 1);
    EXPECT_EQ(phase_prefixes["PHASE2"].size(), 1);
    
    // The prefix should match what we provided
    EXPECT_TRUE(phase_prefixes["PHASE1"].find(prefix) != phase_prefixes["PHASE1"].end());
    EXPECT_TRUE(phase_prefixes["PHASE2"].find(prefix) != phase_prefixes["PHASE2"].end());
}

TEST(UniformSection, PhaseStatePrefixesWithMultiplePhases)
{
    micm::Phase aqueous{ "AQUEOUS", { { micm::Species{ "H2O" } } } };
    micm::Phase organic{ "ORGANIC", { { micm::Species{ "C6H14" } } } };
    micm::Phase gas{ "GAS", { { micm::Species{ "CO2" } } } };
    std::vector<micm::Phase> phases = { aqueous, organic, gas };
    
    const std::string prefix = "BIN_03";
    auto model = UniformSection{ prefix, phases };
    
    auto phase_prefixes = model.PhaseStatePrefixes();
    
    // Should have entries for all three phases
    EXPECT_EQ(phase_prefixes.size(), 3);
    EXPECT_EQ(phase_prefixes["AQUEOUS"].size(), 1);
    EXPECT_EQ(phase_prefixes["ORGANIC"].size(), 1);
    EXPECT_EQ(phase_prefixes["GAS"].size(), 1);
    
    // All should contain the same prefix
    EXPECT_TRUE(phase_prefixes["AQUEOUS"].find(prefix) != phase_prefixes["AQUEOUS"].end());
    EXPECT_TRUE(phase_prefixes["ORGANIC"].find(prefix) != phase_prefixes["ORGANIC"].end());
    EXPECT_TRUE(phase_prefixes["GAS"].find(prefix) != phase_prefixes["GAS"].end());
}
