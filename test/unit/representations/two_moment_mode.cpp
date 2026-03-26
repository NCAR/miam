// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../representation_policy.hpp"

#include <miam/aerosol_property.hpp>
#include <miam/representations/two_moment_mode.hpp>

#include <micm/util/matrix.hpp>

#include <cmath>
#include <numbers>
#include <unordered_map>

using namespace miam;
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
    micm::Species test_species{ "CO2" };
    micm::Phase test_phase{ "AQUEOUS", { { test_species } } };
    std::vector<micm::Phase> phases = { test_phase };
    
    auto model = TwoMomentMode{ test_model_name, phases };
    
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

TEST(TwoMomentMode, NumPhaseInstances)
{
    const auto phases = getTestPhases();
    auto model = TwoMomentMode{ test_model_name, phases };
    
    auto num_instances = model.NumPhaseInstances();
    
    // Should have one instance per phase (2 phases in test setup)
    EXPECT_EQ(num_instances.size(), 2);
    EXPECT_EQ(num_instances["PHASE1"], 1);
    EXPECT_EQ(num_instances["PHASE2"], 1);
}

TEST(TwoMomentMode, PhaseStatePrefixes)
{
    const auto phases = getTestPhases();
    const std::string prefix = "MODE1";
    auto model = TwoMomentMode{ prefix, phases };
    
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

TEST(TwoMomentMode, PhaseStatePrefixesWithMultiplePhases)
{
    micm::Phase aqueous{ "AQUEOUS", { { micm::Species{ "H2O" } } } };
    micm::Phase organic{ "ORGANIC", { { micm::Species{ "C6H14" } } } };
    micm::Phase gas{ "GAS", { { micm::Species{ "CO2" } } } };
    std::vector<micm::Phase> phases = { aqueous, organic, gas };
    
    const std::string prefix = "AITKEN";
    auto model = TwoMomentMode{ prefix, phases };
    
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

// ============================================================================
// Property Provider Tests
// ============================================================================

namespace
{
  using Matrix = micm::Matrix<double>;

  micm::Species MakeSpecies(const std::string& name, double mw, double rho)
  {
    return micm::Species{ name, { { "molecular weight [kg mol-1]", mw }, { "density [kg m-3]", rho } } };
  }

  const double mw_A = 0.018;
  const double rho_A = 1000.0;
  const double mw_B = 0.044;
  const double rho_B = 2000.0;
  const double mwr_A = mw_A / rho_A;
  const double mwr_B = mw_B / rho_B;

  micm::Phase MakeProviderTestPhase()
  {
    return micm::Phase{ "aqueous", { { MakeSpecies("A", mw_A, rho_A) }, { MakeSpecies("B", mw_B, rho_B) } } };
  }
}  // namespace

TEST(TwoMomentMode, ProviderEffectiveRadius)
{
  auto phase = MakeProviderTestPhase();
  TwoMomentMode mode("MODE1", { phase }, 2.0);
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 0 } };
  std::unordered_map<std::string, std::size_t> var_idx = {
    { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 }, { "MODE1.NUMBER_CONCENTRATION", 2 }
  };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::EffectiveRadius, param_idx, var_idx);

  ASSERT_EQ(provider.dependent_variable_indices.size(), 3u);
  EXPECT_EQ(provider.dependent_variable_indices[0], 0u);
  EXPECT_EQ(provider.dependent_variable_indices[1], 1u);
  EXPECT_EQ(provider.dependent_variable_indices[2], 2u);

  double gsd = 2.0;
  double conc_A = 100.0;
  double conc_B = 200.0;
  double N = 1.0e12;

  Matrix params{ 1, 1, 0.0 };
  Matrix vars{ 1, 3, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  params[0][0] = gsd;
  vars[0][0] = conc_A;
  vars[0][1] = conc_B;
  vars[0][2] = N;

  provider.ComputeValue(params, vars, result);

  double V_total = conc_A * mwr_A + conc_B * mwr_B;
  double ln_gsd = std::log(gsd);
  double V_mean = V_total / N;
  double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
  double expected = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
  EXPECT_NEAR(result[0][0], expected, std::abs(expected) * 1e-10);

  Matrix result2{ 1, 1, 0.0 };
  Matrix partials{ 1, 3, 0.0 };
  provider.ComputeValueAndDerivatives(params, vars, result2, partials);
  EXPECT_NEAR(result2[0][0], expected, std::abs(expected) * 1e-10);

  double dr_dA = expected * mwr_A / (3.0 * V_total);
  EXPECT_NEAR(partials[0][0], dr_dA, std::abs(dr_dA) * 1e-10);

  double dr_dB = expected * mwr_B / (3.0 * V_total);
  EXPECT_NEAR(partials[0][1], dr_dB, std::abs(dr_dB) * 1e-10);

  double dr_dN = -expected / (3.0 * N);
  EXPECT_NEAR(partials[0][2], dr_dN, std::abs(dr_dN) * 1e-10);
}

TEST(TwoMomentMode, ProviderNumberConcentration)
{
  auto phase = MakeProviderTestPhase();
  TwoMomentMode mode("MODE1", { phase }, 2.0);
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 0 } };
  std::unordered_map<std::string, std::size_t> var_idx = {
    { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 }, { "MODE1.NUMBER_CONCENTRATION", 2 }
  };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::NumberConcentration, param_idx, var_idx);

  ASSERT_EQ(provider.dependent_variable_indices.size(), 1u);
  EXPECT_EQ(provider.dependent_variable_indices[0], 2u);

  double N_input = 1.5e12;
  Matrix params{ 1, 1, 0.0 };
  Matrix vars{ 1, 3, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  params[0][0] = 2.0;
  vars[0][2] = N_input;

  provider.ComputeValue(params, vars, result);
  EXPECT_DOUBLE_EQ(result[0][0], N_input);

  Matrix result2{ 1, 1, 0.0 };
  Matrix partials{ 1, 1, 0.0 };
  provider.ComputeValueAndDerivatives(params, vars, result2, partials);
  EXPECT_DOUBLE_EQ(result2[0][0], N_input);
  EXPECT_DOUBLE_EQ(partials[0][0], 1.0);
}

TEST(TwoMomentMode, ProviderPhaseVolumeFractionSinglePhase)
{
  auto phase = MakeProviderTestPhase();
  TwoMomentMode mode("MODE1", { phase }, 2.0);
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 0 } };
  std::unordered_map<std::string, std::size_t> var_idx = {
    { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 }, { "MODE1.NUMBER_CONCENTRATION", 2 }
  };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::PhaseVolumeFraction, param_idx, var_idx);
  EXPECT_TRUE(provider.dependent_variable_indices.empty());

  Matrix params{ 1, 1, 0.0 };
  Matrix vars{ 1, 3, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  vars[0][0] = 100.0;
  vars[0][1] = 200.0;
  vars[0][2] = 1e12;

  provider.ComputeValue(params, vars, result);
  EXPECT_DOUBLE_EQ(result[0][0], 1.0);
}

// ============================================================================
// Policy Tests — cross-representation property contracts
// ============================================================================

TEST(TwoMomentMode, PolicyPhaseStatePrefixes)
{
    testPhaseStatePrefixes<TwoMomentMode>();
}

TEST(TwoMomentMode, PolicyPhaseVolumeFractionSinglePhase)
{
    testPhaseVolumeFractionSinglePhase<TwoMomentMode>();
}

TEST(TwoMomentMode, PolicyPhaseVolumeFractionMultiPhase)
{
    testPhaseVolumeFractionMultiPhase<TwoMomentMode>();
}

TEST(TwoMomentMode, PolicyEffectiveRadius)
{
    testEffectiveRadiusProvider<TwoMomentMode>();
}

TEST(TwoMomentMode, PolicyNumberConcentration)
{
    testNumberConcentrationProvider<TwoMomentMode>();
}
