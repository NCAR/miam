// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../representation_policy.hpp"

#include <miam/representations/two_moment_mode.hpp>

using namespace miam::representation;

TEST(TwoMomentMode, StateSize)
{
    auto size = testStateSize<TwoMomentMode>();
    // Should have 3 species + 1 number concentration = 4 variables, and 1 parameter (GSD)
    EXPECT_EQ(std::get<0>(size), 4);
    EXPECT_EQ(std::get<1>(size), 1);
}

TEST(TwoMomentMode, StateVariableNames)
{
    auto names = testStateVariableNames<TwoMomentMode>();
    // Should have 3 species + number concentration = 4 variables
    EXPECT_EQ(names.size(), 4);
    EXPECT_TRUE(names.find(test_model_name + ".PHASE1.SPECIES_A") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".PHASE1.SPECIES_B") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".PHASE2.SPECIES_C") != names.end());
    EXPECT_TRUE(names.find(test_model_name + ".NUMBER_CONCENTRATION") != names.end());
}

TEST(TwoMomentMode, StateParameterNames)
{
    auto names = testStateParameterNames<TwoMomentMode>();
    // Should have 1 parameter: GEOMETRIC_STANDARD_DEVIATION
    EXPECT_EQ(names.size(), 1);
    EXPECT_TRUE(names.find(test_model_name + ".GEOMETRIC_STANDARD_DEVIATION") != names.end());
}

TEST(TwoMomentMode, DefaultParametersPolicy)
{
    auto params = testDefaultParameters<TwoMomentMode>();
    
    // Check default value (1.0 for GSD when using default constructor)
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_STANDARD_DEVIATION"], 1.0);
}

TEST(TwoMomentMode, DefaultParameters)
{
    const auto phases = getTestPhases();
    auto model = TwoMomentMode{ test_model_name, phases };
    auto params = model.DefaultParameters();
    
    // Should have 1 default parameter
    EXPECT_EQ(params.size(), 1);
    
    // Check parameter name exists
    EXPECT_TRUE(params.find(test_model_name + ".GEOMETRIC_STANDARD_DEVIATION") != params.end());
    
    // Check default value (1.0 for GSD when using default constructor)
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_STANDARD_DEVIATION"], 1.0);
}

TEST(TwoMomentMode, DefaultParametersWithCustomValues)
{
    const auto phases = getTestPhases();
    const double custom_gsd = 1.5;
    auto model = TwoMomentMode{ test_model_name, phases, custom_gsd };
    auto params = model.DefaultParameters();
    
    // Check custom default value
    EXPECT_EQ(params[test_model_name + ".GEOMETRIC_STANDARD_DEVIATION"], custom_gsd);
}

TEST(TwoMomentMode, GeometricStandardDeviationParameterName)
{
    const auto phases = getTestPhases();
    auto model = TwoMomentMode{ test_model_name, phases };
    
    EXPECT_EQ(model.GeometricStandardDeviation(), test_model_name + ".GEOMETRIC_STANDARD_DEVIATION");
}

TEST(TwoMomentMode, NumberConcentrationVariableName)
{
    const auto phases = getTestPhases();
    auto model = TwoMomentMode{ test_model_name, phases };
    
    EXPECT_EQ(model.NumberConcentration(), test_model_name + ".NUMBER_CONCENTRATION");
}

TEST(TwoMomentMode, SpeciesNaming)
{
    const auto phases = getTestPhases();
    auto model = TwoMomentMode{ test_model_name, phases };
    
    micm::Phase test_phase{ "AQUEOUS", { { micm::Species{ "H2O" } } } };
    micm::Species test_species{ "CO2" };
    
    std::string species_name = model.Species(test_phase, test_species);
    EXPECT_EQ(species_name, test_model_name + ".AQUEOUS.CO2");
}

TEST(TwoMomentMode, MultiplePhases)
{
    micm::Phase aqueous{ "AQUEOUS", { { micm::Species{ "H2O" } }, { micm::Species{ "CO2" } } } };
    micm::Phase organic{ "ORGANIC", { { micm::Species{ "C6H14" } } } };
    std::vector<micm::Phase> phases = { aqueous, organic };
    
    auto model = TwoMomentMode{ "MODE1", phases };
    
    auto size = model.StateSize();
    // 3 species + 1 number concentration = 4 variables, and 1 parameter
    EXPECT_EQ(std::get<0>(size), 4);
    EXPECT_EQ(std::get<1>(size), 1);
    
    auto var_names = model.StateVariableNames();
    EXPECT_EQ(var_names.size(), 4);
    // Names should include phase: MODE1.AQUEOUS.H2O, MODE1.AQUEOUS.CO2, MODE1.ORGANIC.C6H14, MODE1.NUMBER_CONCENTRATION
    EXPECT_TRUE(var_names.find("MODE1.AQUEOUS.H2O") != var_names.end());
    EXPECT_TRUE(var_names.find("MODE1.AQUEOUS.CO2") != var_names.end());
    EXPECT_TRUE(var_names.find("MODE1.ORGANIC.C6H14") != var_names.end());
    EXPECT_TRUE(var_names.find("MODE1.NUMBER_CONCENTRATION") != var_names.end());
}

TEST(TwoMomentMode, PrefixConsistency)
{
    const auto phases = getTestPhases();
    const std::string custom_prefix = "CUSTOM_MODE";
    auto model = TwoMomentMode{ custom_prefix, phases };
    
    // Check that all method outputs use the custom prefix
    EXPECT_EQ(model.GeometricStandardDeviation(), custom_prefix + ".GEOMETRIC_STANDARD_DEVIATION");
    EXPECT_EQ(model.NumberConcentration(), custom_prefix + ".NUMBER_CONCENTRATION");
    
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

TEST(TwoMomentMode, NumberConcentrationInStateVariables)
{
    const auto phases = getTestPhases();
    auto model = TwoMomentMode{ test_model_name, phases };
    
    auto var_names = model.StateVariableNames();
    
    // NUMBER_CONCENTRATION should be a state variable
    EXPECT_TRUE(var_names.find(test_model_name + ".NUMBER_CONCENTRATION") != var_names.end());
    
    auto param_names = model.StateParameterNames();
    
    // NUMBER_CONCENTRATION should NOT be a state parameter
    EXPECT_FALSE(param_names.find(test_model_name + ".NUMBER_CONCENTRATION") != param_names.end());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
