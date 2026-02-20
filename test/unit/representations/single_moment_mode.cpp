// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../representation_policy.hpp"

#include <miam/representations/single_moment_mode.hpp>

using namespace miam::representation;

TEST(SingleMomentMode, StateSize)
{
    auto size = testStateSize<SingleMomentMode>();
    // Should have 3 species from test phases and 2 parameters (GMD and GSD)
    EXPECT_EQ(std::get<0>(size), 3);
    EXPECT_EQ(std::get<1>(size), 2);
}

TEST(SingleMomentMode, StateVariableNames)
{
    auto names = testStateVariableNames<SingleMomentMode>();
    // Should have 3 species with phase-qualified names
    EXPECT_EQ(names.size(), 3);
    EXPECT_TRUE(names.find(test_model_name + ".PHASE1.SPECIES_A") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".PHASE1.SPECIES_B") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".PHASE2.SPECIES_C") != names.end());
}

TEST(SingleMomentMode, StateParameterNames)
{
    auto names = testStateParameterNames<SingleMomentMode>();
    // Should have 2 parameters: GEOMETRIC_MEAN_RADIUS and GEOMETRIC_STANDARD_DEVIATION
    EXPECT_EQ(names.size(), 2);
    EXPECT_TRUE(names.find(test_model_name + ".GEOMETRIC_MEAN_RADIUS") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".GEOMETRIC_STANDARD_DEVIATION") != names.end());
}

TEST(SingleMomentMode, DefaultParametersPolicy)
{
    auto params = testDefaultParameters<SingleMomentMode>();
    
    // Check default values (0.0 for GMD, 1.0 for GSD when using default constructor)
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_MEAN_RADIUS"], 0.0);
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_STANDARD_DEVIATION"], 1.0);
}

TEST(SingleMomentMode, DefaultParameters)
{
    const auto phases = getTestPhases();
    auto model = SingleMomentMode{ test_model_name, phases };
    auto params = model.DefaultParameters();
    
    // Should have 2 default parameters
    EXPECT_EQ(params.size(), 2);
    
    // Check parameter names exist
    EXPECT_TRUE(params.find(test_model_name + ".GEOMETRIC_MEAN_RADIUS") != params.end());
    EXPECT_TRUE(params.find(test_model_name + ".GEOMETRIC_STANDARD_DEVIATION") != params.end());
    
    // Check default values (0.0 for GMD, 1.0 for GSD when using default constructor)
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_MEAN_RADIUS"], 0.0);
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_STANDARD_DEVIATION"], 1.0);
}

TEST(SingleMomentMode, DefaultParametersWithCustomValues)
{
    const auto phases = getTestPhases();
    const double custom_gmd = 1.5e-7;
    const double custom_gsd = 1.3;
    auto model = SingleMomentMode{ test_model_name, phases, custom_gmd, custom_gsd };
    auto params = model.DefaultParameters();
    
    // Check custom default values
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_MEAN_RADIUS"], custom_gmd);
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_STANDARD_DEVIATION"], custom_gsd);
}

TEST(SingleMomentMode, GeometricMeanRadiusParameterName)
{
    const auto phases = getTestPhases();
    auto model = SingleMomentMode{ test_model_name, phases };
    
    EXPECT_EQ(model.GeometricMeanRadius(), test_model_name + ".GEOMETRIC_MEAN_RADIUS");
}

TEST(SingleMomentMode, GeometricStandardDeviationParameterName)
{
    const auto phases = getTestPhases();
    auto model = SingleMomentMode{ test_model_name, phases };
    
    EXPECT_EQ(model.GeometricStandardDeviation(), test_model_name + ".GEOMETRIC_STANDARD_DEVIATION");
}

TEST(SingleMomentMode, SpeciesNaming)
{
    const auto phases = getTestPhases();
    auto model = SingleMomentMode{ test_model_name, phases };
    
    micm::Phase test_phase{ "AQUEOUS", { { micm::Species{ "H2O" } } } };
    micm::Species test_species{ "CO2" };
    
    std::string species_name = model.Species(test_phase, test_species);
    EXPECT_EQ(species_name, test_model_name + ".AQUEOUS.CO2");
}

TEST(SingleMomentMode, MultiplePhases)
{
    micm::Phase aqueous{ "AQUEOUS", { { micm::Species{ "H2O" } }, { micm::Species{ "CO2" } } } };
    micm::Phase organic{ "ORGANIC", { { micm::Species{ "C6H14" } } } };
    std::vector<micm::Phase> phases = { aqueous, organic };
    
    auto model = SingleMomentMode{ "MODE1", phases };
    
    auto size = model.StateSize();
    // 3 species (2 from aqueous + 1 from organic) and 2 parameters
    EXPECT_EQ(std::get<0>(size), 3);
    EXPECT_EQ(std::get<1>(size), 2);
    
    auto var_names = model.StateVariableNames();
    EXPECT_EQ(var_names.size(), 3);
    // Names should include phase: MODE1.AQUEOUS.H2O, MODE1.AQUEOUS.CO2, MODE1.ORGANIC.C6H14
    EXPECT_TRUE(var_names.find("MODE1.AQUEOUS.H2O") != var_names.end());
    EXPECT_TRUE(var_names.find("MODE1.AQUEOUS.CO2") != var_names.end());
    EXPECT_TRUE(var_names.find("MODE1.ORGANIC.C6H14") != var_names.end());
}

TEST(SingleMomentMode, PrefixConsistency)
{
    const auto phases = getTestPhases();
    const std::string custom_prefix = "CUSTOM_MODE";
    auto model = SingleMomentMode{ custom_prefix, phases };
    
    // Check that all method outputs use the custom prefix
    EXPECT_EQ(model.GeometricMeanRadius(), custom_prefix + ".GEOMETRIC_MEAN_RADIUS");
    EXPECT_EQ(model.GeometricStandardDeviation(), custom_prefix + ".GEOMETRIC_STANDARD_DEVIATION");
    
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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
