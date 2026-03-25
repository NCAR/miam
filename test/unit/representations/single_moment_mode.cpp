// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../representation_policy.hpp"

#include <miam/aerosol_property.hpp>
#include <miam/representations/single_moment_mode.hpp>

#include <micm/util/matrix.hpp>

#include <cmath>
#include <numbers>
#include <unordered_map>

using namespace miam;
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
    micm::Species test_species{ "CO2" };
    micm::Phase test_phase{ "AQUEOUS", { { test_species } } };
    std::vector<micm::Phase> phases = { test_phase };
    
    auto model = SingleMomentMode{ test_model_name, phases };
    
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

TEST(SingleMomentMode, NumPhaseInstances)
{
    const auto phases = getTestPhases();
    auto model = SingleMomentMode{ test_model_name, phases };
    
    auto num_instances = model.NumPhaseInstances();
    
    // Should have one instance per phase (2 phases in test setup)
    EXPECT_EQ(num_instances.size(), 2);
    EXPECT_EQ(num_instances["PHASE1"], 1);
    EXPECT_EQ(num_instances["PHASE2"], 1);
}

TEST(SingleMomentMode, PhaseStatePrefixes)
{
    const auto phases = getTestPhases();
    const std::string prefix = "MODE1";
    auto model = SingleMomentMode{ prefix, phases };
    
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

TEST(SingleMomentMode, PhaseStatePrefixesWithMultiplePhases)
{
    micm::Phase aqueous{ "AQUEOUS", { { micm::Species{ "H2O" } } } };
    micm::Phase organic{ "ORGANIC", { { micm::Species{ "C6H14" } } } };
    micm::Phase gas{ "GAS", { { micm::Species{ "CO2" } } } };
    std::vector<micm::Phase> phases = { aqueous, organic, gas };
    
    const std::string prefix = "SMALL_DROP";
    auto model = SingleMomentMode{ prefix, phases };
    
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

TEST(SingleMomentMode, ProviderEffectiveRadius)
{
  auto phase = MakeProviderTestPhase();
  SingleMomentMode mode("MODE1", { phase }, 1.0e-6, 2.0);
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_MEAN_RADIUS", 0 },
                                                             { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 1 } };
  std::unordered_map<std::string, std::size_t> var_idx = { { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 } };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::EffectiveRadius, param_idx, var_idx);
  EXPECT_TRUE(provider.dependent_variable_indices.empty());

  Matrix params{ 1, 2, 0.0 };
  Matrix vars{ 1, 2, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  params[0][0] = 1.0e-6;
  params[0][1] = 2.0;

  provider.ComputeValue(params, vars, result);

  double ln_gsd = std::log(2.0);
  double expected = 1.0e-6 * std::exp(2.5 * ln_gsd * ln_gsd);
  EXPECT_NEAR(result[0][0], expected, expected * 1e-10);

  Matrix result2{ 1, 1, 0.0 };
  Matrix partials{ 1, 0, 0.0 };
  provider.ComputeValueAndDerivatives(params, vars, result2, partials);
  EXPECT_NEAR(result2[0][0], expected, expected * 1e-10);
}

TEST(SingleMomentMode, ProviderNumberConcentration)
{
  auto phase = MakeProviderTestPhase();
  SingleMomentMode mode("MODE1", { phase }, 1.0e-6, 2.0);
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_MEAN_RADIUS", 0 },
                                                             { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 1 } };
  std::unordered_map<std::string, std::size_t> var_idx = { { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 } };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::NumberConcentration, param_idx, var_idx);
  ASSERT_EQ(provider.dependent_variable_indices.size(), 2u);

  double gmd = 1.0e-6;
  double gsd = 2.0;
  double conc_A = 100.0;
  double conc_B = 200.0;

  Matrix params{ 1, 2, 0.0 };
  Matrix vars{ 1, 2, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  params[0][0] = gmd;
  params[0][1] = gsd;
  vars[0][0] = conc_A;
  vars[0][1] = conc_B;

  provider.ComputeValue(params, vars, result);

  double V_total = conc_A * mwr_A + conc_B * mwr_B;
  double ln_gsd = std::log(gsd);
  double V_single = (4.0 / 3.0) * std::numbers::pi * gmd * gmd * gmd * std::exp(4.5 * ln_gsd * ln_gsd);
  double expected_N = V_total / V_single;
  EXPECT_NEAR(result[0][0], expected_N, std::abs(expected_N) * 1e-10);

  Matrix result2{ 1, 1, 0.0 };
  Matrix partials{ 1, 2, 0.0 };
  provider.ComputeValueAndDerivatives(params, vars, result2, partials);
  EXPECT_NEAR(result2[0][0], expected_N, std::abs(expected_N) * 1e-10);
  EXPECT_NEAR(partials[0][0], mwr_A / V_single, std::abs(mwr_A / V_single) * 1e-10);
  EXPECT_NEAR(partials[0][1], mwr_B / V_single, std::abs(mwr_B / V_single) * 1e-10);
}

TEST(SingleMomentMode, ProviderPhaseVolumeFractionSinglePhase)
{
  auto phase = MakeProviderTestPhase();
  SingleMomentMode mode("MODE1", { phase });
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_MEAN_RADIUS", 0 },
                                                             { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 1 } };
  std::unordered_map<std::string, std::size_t> var_idx = { { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 } };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::PhaseVolumeFraction, param_idx, var_idx);
  EXPECT_TRUE(provider.dependent_variable_indices.empty());

  Matrix params{ 1, 2, 0.0 };
  Matrix vars{ 1, 2, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  vars[0][0] = 100.0;
  vars[0][1] = 200.0;

  provider.ComputeValue(params, vars, result);
  EXPECT_DOUBLE_EQ(result[0][0], 1.0);

  Matrix result2{ 1, 1, 0.0 };
  Matrix partials{ 1, 0, 0.0 };
  provider.ComputeValueAndDerivatives(params, vars, result2, partials);
  EXPECT_DOUBLE_EQ(result2[0][0], 1.0);
}

TEST(SingleMomentMode, ProviderPhaseVolumeFractionMultiPhase)
{
  auto species_A = MakeSpecies("A", mw_A, rho_A);
  auto species_B = MakeSpecies("B", mw_B, rho_B);
  micm::Phase aqueous("aqueous", { { species_A } });
  micm::Phase organic("organic", { { species_B } });
  SingleMomentMode mode("MODE1", { aqueous, organic }, 1.0e-6, 2.0);

  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_MEAN_RADIUS", 0 },
                                                             { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 1 } };
  std::unordered_map<std::string, std::size_t> var_idx = { { "MODE1.aqueous.A", 0 }, { "MODE1.organic.B", 1 } };

  auto provider =
      mode.GetPropertyProvider<Matrix>(AerosolProperty::PhaseVolumeFraction, param_idx, var_idx, "aqueous");

  ASSERT_EQ(provider.dependent_variable_indices.size(), 2u);
  EXPECT_EQ(provider.dependent_variable_indices[0], 0u);
  EXPECT_EQ(provider.dependent_variable_indices[1], 1u);

  double conc_A = 100.0;
  double conc_B = 200.0;
  Matrix params{ 1, 2, 0.0 };
  Matrix vars{ 1, 2, 0.0 };
  Matrix result{ 1, 1, 0.0 };
  vars[0][0] = conc_A;
  vars[0][1] = conc_B;

  provider.ComputeValue(params, vars, result);

  double V_aqueous = conc_A * mwr_A;
  double V_organic = conc_B * mwr_B;
  double V_total = V_aqueous + V_organic;
  double expected_phi = V_aqueous / V_total;
  EXPECT_NEAR(result[0][0], expected_phi, 1e-10);

  Matrix result2{ 1, 1, 0.0 };
  Matrix partials{ 1, 2, 0.0 };
  provider.ComputeValueAndDerivatives(params, vars, result2, partials);
  EXPECT_NEAR(result2[0][0], expected_phi, 1e-10);

  double dphi_dA = mwr_A * (1.0 - expected_phi) / V_total;
  EXPECT_NEAR(partials[0][0], dphi_dA, std::abs(dphi_dA) * 1e-10);

  double dphi_dB = -mwr_B * expected_phi / V_total;
  EXPECT_NEAR(partials[0][1], dphi_dB, std::abs(dphi_dB) * 1e-10);
}

TEST(SingleMomentMode, ProviderMultiCell)
{
  auto phase = MakeProviderTestPhase();
  SingleMomentMode mode("MODE1", { phase }, 1.0e-6, 2.0);
  std::unordered_map<std::string, std::size_t> param_idx = { { "MODE1.GEOMETRIC_MEAN_RADIUS", 0 },
                                                             { "MODE1.GEOMETRIC_STANDARD_DEVIATION", 1 } };
  std::unordered_map<std::string, std::size_t> var_idx = { { "MODE1.aqueous.A", 0 }, { "MODE1.aqueous.B", 1 } };

  auto provider = mode.GetPropertyProvider<Matrix>(AerosolProperty::EffectiveRadius, param_idx, var_idx);

  const std::size_t num_cells = 3;
  Matrix params{ num_cells, 2, 0.0 };
  Matrix vars{ num_cells, 2, 0.0 };
  Matrix result{ num_cells, 1, 0.0 };

  double gmds[] = { 1.0e-6, 2.0e-6, 0.5e-6 };
  double gsds[] = { 1.5, 2.0, 1.2 };
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    params[i][0] = gmds[i];
    params[i][1] = gsds[i];
  }

  provider.ComputeValue(params, vars, result);

  for (std::size_t i = 0; i < num_cells; ++i)
  {
    double ln_gsd = std::log(gsds[i]);
    double expected = gmds[i] * std::exp(2.5 * ln_gsd * ln_gsd);
    EXPECT_NEAR(result[i][0], expected, expected * 1e-10) << "Cell " << i;
  }
}
