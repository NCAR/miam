// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/processes/henry_law_phase_transfer.hpp>
#include <miam/processes/henry_law_phase_transfer_builder.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <miam/util/condensation_rate.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace miam::process;
using MatrixPolicy = micm::Matrix<double>;
using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

namespace
{
  // Test constants
  constexpr double D_g = 1.5e-5;     // m^2 s^-1
  constexpr double alpha = 0.05;     // accommodation coefficient
  constexpr double Mw_gas = 0.044;   // kg mol^-1 (CO2)
  constexpr double Mw_solvent = 0.018;  // kg mol^-1 (H2O)
  constexpr double rho_solvent = 1000.0;  // kg m^-3 (H2O)
  constexpr double HLC_ref = 3.4e-2;  // mol m^-3 Pa^-1

  micm::Species MakeGasSpecies()
  {
    return micm::Species{ "CO2_g", { { "molecular weight [kg mol-1]", Mw_gas } } };
  }

  micm::Species MakeCondensedSpecies()
  {
    return micm::Species{ "CO2_aq", { { "molecular weight [kg mol-1]", Mw_gas },
                                       { "density [kg m-3]", 1800.0 } } };
  }

  micm::Species MakeSolvent()
  {
    return micm::Species{ "H2O", { { "molecular weight [kg mol-1]", Mw_solvent },
                                    { "density [kg m-3]", rho_solvent } } };
  }

  micm::Phase MakeAqueousPhase()
  {
    return micm::Phase{ "AQUEOUS", { { MakeCondensedSpecies() }, { MakeSolvent() } } };
  }

  /// Create a simple constant HLC function
  auto MakeConstantHLC(double hlc_value)
  {
    return [hlc_value](const micm::Conditions& conditions) { return hlc_value; };
  }

  /// Create a HenryLawPhaseTransfer process with default test settings
  HenryLawPhaseTransfer MakeTestProcess(double hlc_value = HLC_ref)
  {
    return HenryLawPhaseTransfer(
        MakeConstantHLC(hlc_value),
        MakeGasSpecies(),
        MakeCondensedSpecies(),
        MakeSolvent(),
        MakeAqueousPhase(),
        D_g,
        alpha,
        Mw_gas,
        Mw_solvent,
        rho_solvent);
  }

  /// Create a simple AerosolPropertyProvider that returns a constant value
  template<typename DenseMatrixPolicy>
  miam::AerosolPropertyProvider<DenseMatrixPolicy> MakeConstantProvider(
      double value,
      const std::vector<std::size_t>& dependent_variable_indices = {})
  {
    miam::AerosolPropertyProvider<DenseMatrixPolicy> provider;
    provider.dependent_variable_indices = dependent_variable_indices;
    provider.ComputeValue = [value](const DenseMatrixPolicy& params, const DenseMatrixPolicy& vars, DenseMatrixPolicy& result)
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
        result[row][0] = value;
    };
    provider.ComputeValueAndDerivatives =
        [value](const DenseMatrixPolicy& params, const DenseMatrixPolicy& vars,
                DenseMatrixPolicy& result, DenseMatrixPolicy& partials)
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
        result[row][0] = value;
      // partials remain 0 for constant providers
    };
    return provider;
  }

  /// Create providers with specified values for r_eff, N, phi
  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> MakeTestProviders(
      const std::string& prefix,
      double r_eff,
      double N,
      double phi,
      const std::vector<std::size_t>& r_eff_deps = {},
      const std::vector<std::size_t>& N_deps = {},
      const std::vector<std::size_t>& phi_deps = {})
  {
    std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
    providers[prefix][miam::AerosolProperty::EffectiveRadius] = MakeConstantProvider<MatrixPolicy>(r_eff, r_eff_deps);
    providers[prefix][miam::AerosolProperty::NumberConcentration] = MakeConstantProvider<MatrixPolicy>(N, N_deps);
    providers[prefix][miam::AerosolProperty::PhaseVolumeFraction] = MakeConstantProvider<MatrixPolicy>(phi, phi_deps);
    return providers;
  }
}  // namespace

// ======================== ProcessParameterNames ========================

TEST(HenryLawPhaseTransfer, ProcessParameterNamesSinglePrefix)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  auto names = process.ProcessParameterNames(phase_prefixes);
  EXPECT_EQ(names.size(), 2);
  EXPECT_TRUE(names.count("MODE1.AQUEOUS." + process.uuid_ + ".hlc"));
  EXPECT_TRUE(names.count("MODE1.AQUEOUS." + process.uuid_ + ".temperature"));
}

TEST(HenryLawPhaseTransfer, ProcessParameterNamesMultiplePrefixes)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  auto names = process.ProcessParameterNames(phase_prefixes);
  EXPECT_EQ(names.size(), 4);
  EXPECT_TRUE(names.count("MODE1.AQUEOUS." + process.uuid_ + ".hlc"));
  EXPECT_TRUE(names.count("MODE1.AQUEOUS." + process.uuid_ + ".temperature"));
  EXPECT_TRUE(names.count("MODE2.AQUEOUS." + process.uuid_ + ".hlc"));
  EXPECT_TRUE(names.count("MODE2.AQUEOUS." + process.uuid_ + ".temperature"));
}

TEST(HenryLawPhaseTransfer, ProcessParameterNamesNoMatchingPhase)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["ORGANIC"].insert("MODE1");

  auto names = process.ProcessParameterNames(phase_prefixes);
  EXPECT_EQ(names.size(), 0);
}

// ======================== SpeciesUsed ========================

TEST(HenryLawPhaseTransfer, SpeciesUsedSinglePrefix)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  auto species = process.SpeciesUsed(phase_prefixes);
  EXPECT_EQ(species.size(), 3);
  EXPECT_TRUE(species.count("CO2_g"));
  EXPECT_TRUE(species.count("MODE1.AQUEOUS.CO2_aq"));
  EXPECT_TRUE(species.count("MODE1.AQUEOUS.H2O"));
}

TEST(HenryLawPhaseTransfer, SpeciesUsedMultiplePrefixes)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  auto species = process.SpeciesUsed(phase_prefixes);
  // Gas species (1) + 2 condensed per prefix (4) = 5
  EXPECT_EQ(species.size(), 5);
  EXPECT_TRUE(species.count("CO2_g"));
  EXPECT_TRUE(species.count("MODE1.AQUEOUS.CO2_aq"));
  EXPECT_TRUE(species.count("MODE1.AQUEOUS.H2O"));
  EXPECT_TRUE(species.count("MODE2.AQUEOUS.CO2_aq"));
  EXPECT_TRUE(species.count("MODE2.AQUEOUS.H2O"));
}

TEST(HenryLawPhaseTransfer, SpeciesUsedNoMatchingPhaseThrows)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["ORGANIC"].insert("MODE1");

  EXPECT_THROW(process.SpeciesUsed(phase_prefixes), std::runtime_error);
}

// ======================== RequiredAerosolProperties ========================

TEST(HenryLawPhaseTransfer, RequiredAerosolProperties)
{
  auto process = MakeTestProcess();
  auto required = process.RequiredAerosolProperties();

  EXPECT_EQ(required.size(), 1);
  EXPECT_TRUE(required.count("AQUEOUS"));

  const auto& props = required.at("AQUEOUS");
  EXPECT_EQ(props.size(), 3);
  EXPECT_EQ(props[0], miam::AerosolProperty::EffectiveRadius);
  EXPECT_EQ(props[1], miam::AerosolProperty::NumberConcentration);
  EXPECT_EQ(props[2], miam::AerosolProperty::PhaseVolumeFraction);
}

// ======================== NonZeroJacobianElements ========================

TEST(HenryLawPhaseTransfer, NonZeroJacobianElementsSinglePrefix)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  auto elements = process.NonZeroJacobianElements(phase_prefixes, state_variable_indices);

  // 6 direct entries: (0,0), (0,1), (0,2), (1,0), (1,1), (1,2)
  EXPECT_EQ(elements.size(), 6);
  EXPECT_TRUE(elements.count({ 0, 0 }));  // gas,gas
  EXPECT_TRUE(elements.count({ 0, 1 }));  // gas,aq
  EXPECT_TRUE(elements.count({ 0, 2 }));  // gas,solvent
  EXPECT_TRUE(elements.count({ 1, 0 }));  // aq,gas
  EXPECT_TRUE(elements.count({ 1, 1 }));  // aq,aq
  EXPECT_TRUE(elements.count({ 1, 2 }));  // aq,solvent
}

TEST(HenryLawPhaseTransfer, NonZeroJacobianElementsWithProviders)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;
  state_variable_indices["MODE1.AQUEOUS.EXTRA"] = 3;

  // N provider depends on state variable index 1 and 2
  auto providers = MakeTestProviders("MODE1", 1e-6, 1e8, 1.0,
                                     {},          // r_eff: no deps
                                     { 1, 2 },    // N: depends on aq and solvent
                                     { 1 });      // phi: depends on aq

  auto elements = process.NonZeroJacobianElements<MatrixPolicy>(
      phase_prefixes, state_variable_indices, providers);

  // 6 direct + N indirect: 2*2 = 4 (gas,1),(aq,1),(gas,2),(aq,2) + phi indirect: 2*1 = 2 (gas,1),(aq,1)
  // But many overlap with direct entries. Let's just check key additions.
  EXPECT_GE(elements.size(), 6);  // At least 6 direct
  EXPECT_TRUE(elements.count({ 0, 1 }));  // gas,aq (N dep var + direct)
  EXPECT_TRUE(elements.count({ 0, 2 }));  // gas,solvent (N dep var + direct)
  EXPECT_TRUE(elements.count({ 1, 1 }));  // aq,aq (N dep var + direct)
  EXPECT_TRUE(elements.count({ 1, 2 }));  // aq,solvent (N dep var + direct)
}

// ======================== UpdateStateParametersFunction ========================

TEST(HenryLawPhaseTransfer, UpdateStateParametersFunction)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  auto update_func = process.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);

  MatrixPolicy state_parameters(1, 2, 0.0);

  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = 298.15;

  update_func(conditions, state_parameters);

  EXPECT_NEAR(state_parameters[0][0], HLC_ref, 1e-10);
  EXPECT_NEAR(state_parameters[0][1], 298.15, 1e-10);
}

TEST(HenryLawPhaseTransfer, UpdateStateParametersFunctionMultipleCells)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  auto update_func = process.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);

  MatrixPolicy state_parameters(3, 2, 0.0);

  std::vector<micm::Conditions> conditions(3);
  conditions[0].temperature_ = 280.0;
  conditions[1].temperature_ = 298.15;
  conditions[2].temperature_ = 310.0;

  update_func(conditions, state_parameters);

  EXPECT_NEAR(state_parameters[0][0], HLC_ref, 1e-10);
  EXPECT_NEAR(state_parameters[0][1], 280.0, 1e-10);
  EXPECT_NEAR(state_parameters[1][1], 298.15, 1e-10);
  EXPECT_NEAR(state_parameters[2][1], 310.0, 1e-10);
}

// ======================== CopyWithNewUuid ========================

TEST(HenryLawPhaseTransfer, CopyWithNewUuid)
{
  auto process = MakeTestProcess();
  auto copy = process.CopyWithNewUuid();

  EXPECT_NE(process.uuid_, copy.uuid_);
  EXPECT_EQ(process.gas_species_.name_, copy.gas_species_.name_);
  EXPECT_DOUBLE_EQ(process.D_g_, copy.D_g_);
  EXPECT_DOUBLE_EQ(process.alpha_, copy.alpha_);
  EXPECT_DOUBLE_EQ(process.Mw_gas_, copy.Mw_gas_);
}

// ======================== ForcingFunction ========================

TEST(HenryLawPhaseTransfer, ForcingFunctionBasicRates)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  double r_eff_val = 1.0e-6;   // 1 μm
  double N_val = 1.0e8;        // 10^8 m^-3
  double phi_val = 1.0e-6;     // volume fraction

  auto providers = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);

  auto forcing_func = process.ForcingFunction<MatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, providers);

  double T = 298.15;
  double hlc = HLC_ref;
  double gas_conc = 1.0e-3;       // mol m^-3
  double aq_conc = 1.0e-5;        // mol m^-3
  double solvent_conc = 55000.0;   // mol m^-3

  MatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = hlc;
  state_parameters[0][1] = T;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = aq_conc;
  state_variables[0][2] = solvent_conc;

  MatrixPolicy forcing_terms(1, 3, 0.0);

  forcing_func(state_parameters, state_variables, forcing_terms);

  // Compute expected net rate
  auto cond_rate_provider = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);
  double kc = cond_rate_provider.ComputeValue(r_eff_val, N_val, T);
  double kc_eff = phi_val * kc;
  double ke_eff = kc_eff / (hlc * miam::util::R_gas * T);
  double f_v = solvent_conc * Mw_solvent / rho_solvent;
  double expected_net = kc_eff * gas_conc - ke_eff * aq_conc / f_v;

  EXPECT_NEAR(forcing_terms[0][0], -expected_net, std::abs(expected_net) * 1e-10);
  EXPECT_NEAR(forcing_terms[0][1], expected_net, std::abs(expected_net) * 1e-10);
  EXPECT_NEAR(forcing_terms[0][2], 0.0, 1e-30);  // solvent unchanged
}

TEST(HenryLawPhaseTransfer, ForcingFunctionMultipleCells)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  double r_eff_val = 5.0e-6;
  double N_val = 1.0e9;
  double phi_val = 1.0e-5;

  auto providers = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);

  auto forcing_func = process.ForcingFunction<MatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, providers);

  std::size_t num_cells = 3;
  MatrixPolicy state_parameters(num_cells, 2);
  MatrixPolicy state_variables(num_cells, 3);
  MatrixPolicy forcing_terms(num_cells, 3, 0.0);

  double T = 298.15;
  double hlc = HLC_ref;
  double solvent = 55000.0;

  for (std::size_t i = 0; i < num_cells; ++i)
  {
    state_parameters[i][0] = hlc;
    state_parameters[i][1] = T;
    state_variables[i][0] = 1.0e-3 * (i + 1);  // varying gas concentration
    state_variables[i][1] = 1.0e-5;
    state_variables[i][2] = solvent;
  }

  forcing_func(state_parameters, state_variables, forcing_terms);

  auto cond_rate_provider = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);
  double kc = cond_rate_provider.ComputeValue(r_eff_val, N_val, T);
  double kc_eff = phi_val * kc;
  double ke_eff = kc_eff / (hlc * miam::util::R_gas * T);
  double f_v = solvent * Mw_solvent / rho_solvent;

  for (std::size_t i = 0; i < num_cells; ++i)
  {
    double gas = 1.0e-3 * (i + 1);
    double aq = 1.0e-5;
    double expected_net = kc_eff * gas - ke_eff * aq / f_v;
    EXPECT_NEAR(forcing_terms[i][0], -expected_net, std::abs(expected_net) * 1e-10);
    EXPECT_NEAR(forcing_terms[i][1], expected_net, std::abs(expected_net) * 1e-10);
  }
}

TEST(HenryLawPhaseTransfer, ForcingFunctionMassConservation)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  auto providers = MakeTestProviders("MODE1", 2.0e-6, 5.0e8, 1.0e-4);

  auto forcing_func = process.ForcingFunction<MatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, providers);

  MatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = HLC_ref;
  state_parameters[0][1] = 298.15;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = 0.5;      // substantial gas
  state_variables[0][1] = 0.001;    // some aq
  state_variables[0][2] = 55000.0;

  MatrixPolicy forcing_terms(1, 3, 0.0);
  forcing_func(state_parameters, state_variables, forcing_terms);

  // Mass conservation: gas forcing + aq forcing = 0
  EXPECT_NEAR(forcing_terms[0][0] + forcing_terms[0][1], 0.0, 1e-20);
}

// ======================== JacobianFunction ========================

TEST(HenryLawPhaseTransfer, JacobianFunctionDirectEntries)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  double r_eff_val = 1.0e-6;
  double N_val = 1.0e8;
  double phi_val = 1.0e-6;

  auto providers = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);

  // Build sparse jacobian
  auto elements = process.NonZeroJacobianElements<MatrixPolicy>(
      phase_prefixes, state_variable_indices, providers);
  auto builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
  for (const auto& elem : elements)
    builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian(builder);
  jacobian.Fill(0.0);

  auto jac_func = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);

  double T = 298.15;
  double hlc = HLC_ref;
  double gas_conc = 1.0e-3;
  double aq_conc = 1.0e-5;
  double solvent_conc = 55000.0;

  MatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = hlc;
  state_parameters[0][1] = T;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = aq_conc;
  state_variables[0][2] = solvent_conc;

  jac_func(state_parameters, state_variables, jacobian);

  auto cond_rate_provider = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);
  double kc = cond_rate_provider.ComputeValue(r_eff_val, N_val, T);
  double ke = kc / (hlc * miam::util::R_gas * T);
  double fv = solvent_conc * Mw_solvent / rho_solvent;

  // Stored as -J (MICM convention). -J[gas,gas] = +φ · k_cond
  EXPECT_NEAR(jacobian[0][0][0], phi_val * kc, std::abs(phi_val * kc) * 1e-10);
  // -J[gas,aq] = -φ · k_evap / f_v
  EXPECT_NEAR(jacobian[0][0][1], -phi_val * ke / fv, std::abs(phi_val * ke / fv) * 1e-10);
  // -J[gas,solvent] = +φ · k_evap · [aq] / (f_v · [solvent])
  double expected_gas_solvent = phi_val * ke * aq_conc / (fv * solvent_conc);
  EXPECT_NEAR(jacobian[0][0][2], expected_gas_solvent, std::abs(expected_gas_solvent) * 1e-10);
  // -J[aq,gas] = -φ · k_cond
  EXPECT_NEAR(jacobian[0][1][0], -phi_val * kc, std::abs(phi_val * kc) * 1e-10);
  // -J[aq,aq] = +φ · k_evap / f_v
  EXPECT_NEAR(jacobian[0][1][1], phi_val * ke / fv, std::abs(phi_val * ke / fv) * 1e-10);
  // -J[aq,solvent] = -φ · k_evap · [aq] / (f_v · [solvent])
  EXPECT_NEAR(jacobian[0][1][2], -expected_gas_solvent, std::abs(expected_gas_solvent) * 1e-10);
}

TEST(HenryLawPhaseTransfer, JacobianFunctionSymmetry)
{
  // J[gas,x] = -J[aq,x] for all x (mass conservation in Jacobian)
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  auto providers = MakeTestProviders("MODE1", 2.0e-6, 5.0e8, 1.0e-4);

  auto elements = process.NonZeroJacobianElements<MatrixPolicy>(
      phase_prefixes, state_variable_indices, providers);
  auto builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
  for (const auto& elem : elements)
    builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian(builder);
  jacobian.Fill(0.0);

  auto jac_func = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);

  MatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = HLC_ref;
  state_parameters[0][1] = 298.15;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = 0.1;
  state_variables[0][1] = 0.01;
  state_variables[0][2] = 55000.0;

  jac_func(state_parameters, state_variables, jacobian);

  // J[gas,x] = -J[aq,x] for all x
  EXPECT_NEAR(jacobian[0][0][0] + jacobian[0][1][0], 0.0, 1e-20);  // x = gas
  EXPECT_NEAR(jacobian[0][0][1] + jacobian[0][1][1], 0.0, 1e-20);  // x = aq
  EXPECT_NEAR(jacobian[0][0][2] + jacobian[0][1][2], 0.0, 1e-20);  // x = solvent
}

TEST(HenryLawPhaseTransfer, JacobianFunctionFiniteDifference)
{
  // Verify Jacobian against finite-difference approximation
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  double r_eff_val = 2.0e-6;
  double N_val = 5.0e8;
  double phi_val = 1.0e-4;

  auto providers = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);

  auto elements = process.NonZeroJacobianElements<MatrixPolicy>(
      phase_prefixes, state_variable_indices, providers);
  auto builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
  for (const auto& elem : elements)
    builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian(builder);
  jacobian.Fill(0.0);

  auto forcing_func = process.ForcingFunction<MatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, providers);
  auto jac_func = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);

  MatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = HLC_ref;
  state_parameters[0][1] = 298.15;

  double gas_conc = 1.0e-3;
  double aq_conc = 1.0e-5;
  double solvent_conc = 55000.0;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = aq_conc;
  state_variables[0][2] = solvent_conc;

  jac_func(state_parameters, state_variables, jacobian);

  // Finite difference for each state variable column j
  double eps = 1e-8;
  for (std::size_t j = 0; j < 3; ++j)
  {
    MatrixPolicy vars_plus(1, 3);
    MatrixPolicy vars_minus(1, 3);
    for (std::size_t k = 0; k < 3; ++k)
    {
      vars_plus[0][k] = state_variables[0][k];
      vars_minus[0][k] = state_variables[0][k];
    }
    double h = std::max(std::abs(state_variables[0][j]) * eps, eps);
    vars_plus[0][j] += h;
    vars_minus[0][j] -= h;

    MatrixPolicy forcing_plus(1, 3, 0.0);
    MatrixPolicy forcing_minus(1, 3, 0.0);

    // Need fresh providers for each call since they're consumed
    auto prov_plus = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);
    auto prov_minus = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);
    auto ff_plus = process.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, prov_plus);
    auto ff_minus = process.ForcingFunction<MatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, prov_minus);

    ff_plus(state_parameters, vars_plus, forcing_plus);
    ff_minus(state_parameters, vars_minus, forcing_minus);

    double fd_gas = (forcing_plus[0][0] - forcing_minus[0][0]) / (2.0 * h);
    double fd_aq = (forcing_plus[0][1] - forcing_minus[0][1]) / (2.0 * h);

    double tol = 1e-4;  // relative tolerance for finite difference
    // Jacobian stores -J (MICM convention), finite difference gives +J, so ratio is -1
    if (std::abs(jacobian[0][0][j]) > 1e-30)
      EXPECT_NEAR(fd_gas / jacobian[0][0][j], -1.0, tol) << "J[gas," << j << "]";
    if (std::abs(jacobian[0][1][j]) > 1e-30)
      EXPECT_NEAR(fd_aq / jacobian[0][1][j], -1.0, tol) << "J[aq," << j << "]";
  }
}

TEST(HenryLawPhaseTransfer, JacobianFunctionMultipleCells)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string hlc_param = "MODE1.AQUEOUS." + process.uuid_ + ".hlc";
  std::string temp_param = "MODE1.AQUEOUS." + process.uuid_ + ".temperature";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[hlc_param] = 0;
  state_parameter_indices[temp_param] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  double r_eff_val = 1.0e-6;
  double N_val = 1.0e8;
  double phi_val = 1.0e-6;

  auto providers = MakeTestProviders("MODE1", r_eff_val, N_val, phi_val);

  auto elements = process.NonZeroJacobianElements<MatrixPolicy>(
      phase_prefixes, state_variable_indices, providers);

  std::size_t num_cells = 2;
  auto builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(num_cells);
  for (const auto& elem : elements)
    builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian(builder);
  jacobian.Fill(0.0);

  auto jac_func = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);

  MatrixPolicy state_parameters(num_cells, 2);
  MatrixPolicy state_variables(num_cells, 3);

  double T = 298.15;
  double hlc = HLC_ref;
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    state_parameters[i][0] = hlc;
    state_parameters[i][1] = T;
    state_variables[i][0] = 1.0e-3 * (i + 1);
    state_variables[i][1] = 1.0e-5;
    state_variables[i][2] = 55000.0;
  }

  jac_func(state_parameters, state_variables, jacobian);

  auto cond_rate_provider = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);
  double kc = cond_rate_provider.ComputeValue(r_eff_val, N_val, T);

  // Both cells should have same -J[gas,gas] = +φ · k_cond (MICM convention)
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    EXPECT_NEAR(jacobian[i][0][0], phi_val * kc, std::abs(phi_val * kc) * 1e-10);
    // Mass conservation symmetry per cell
    EXPECT_NEAR(jacobian[i][0][0] + jacobian[i][1][0], 0.0, 1e-20);
  }
}

// ======================== Expanded Scenarios ========================

// ---------------------------------------------------------------------------
// Helpers: reusable FD Jacobian checker and linear-dependent provider factory
// ---------------------------------------------------------------------------
namespace
{
  /// @brief Build a sparse Jacobian matrix with declared sparsity for a single process
  SparseMatrixPolicy BuildJacobian(
      const HenryLawPhaseTransfer& process,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>>& providers,
      std::size_t num_blocks)
  {
    auto elements = process.NonZeroJacobianElements<MatrixPolicy>(
        phase_prefixes, state_variable_indices, providers);
    auto builder = SparseMatrixPolicy::Create(state_variable_indices.size())
                       .SetNumberOfBlocks(num_blocks)
                       .InitialValue(0.0);
    for (const auto& elem : elements)
      builder = builder.WithElement(elem.first, elem.second);
    return SparseMatrixPolicy(builder);
  }

  /// @brief Build a sparse Jacobian matrix with declared sparsity for multiple processes
  SparseMatrixPolicy BuildJacobian(
      const std::vector<std::reference_wrapper<const HenryLawPhaseTransfer>>& processes,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>>& providers,
      std::size_t num_blocks)
  {
    std::set<std::pair<std::size_t, std::size_t>> elements;
    for (const auto& proc : processes)
    {
      auto elts = proc.get().NonZeroJacobianElements<MatrixPolicy>(
          phase_prefixes, state_variable_indices, providers);
      elements.insert(elts.begin(), elts.end());
    }
    auto builder = SparseMatrixPolicy::Create(state_variable_indices.size())
                       .SetNumberOfBlocks(num_blocks)
                       .InitialValue(0.0);
    for (const auto& elem : elements)
      builder = builder.WithElement(elem.first, elem.second);
    return SparseMatrixPolicy(builder);
  }

  /// @brief Compare analytical Jacobian against central finite-difference approximation
  void CheckFiniteDifferenceJacobian(
      const HenryLawPhaseTransfer& process,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>>& providers,
      const MatrixPolicy& state_parameters,
      const MatrixPolicy& state_variables,
      double rel_tol = 1e-4)
  {
    std::size_t num_blocks = state_parameters.NumRows();
    std::size_t num_vars = state_variable_indices.size();

    // Build analytical Jacobian
    auto jacobian = BuildJacobian(process, phase_prefixes, state_variable_indices, providers, num_blocks);
    auto jac_func = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);
    jac_func(state_parameters, state_variables, jacobian);

    // Build forcing function for FD
    double eps = 1e-7;
    for (std::size_t j = 0; j < num_vars; ++j)
    {
      MatrixPolicy vars_plus(state_variables);
      MatrixPolicy vars_minus(state_variables);
      for (std::size_t b = 0; b < num_blocks; ++b)
      {
        double h = std::max(std::abs(state_variables[b][j]) * eps, eps);
        vars_plus[b][j] = state_variables[b][j] + h;
        vars_minus[b][j] = state_variables[b][j] - h;
      }

      MatrixPolicy forcing_plus(num_blocks, num_vars, 0.0);
      MatrixPolicy forcing_minus(num_blocks, num_vars, 0.0);

      auto ff_plus = process.ForcingFunction<MatrixPolicy>(
          phase_prefixes, state_parameter_indices, state_variable_indices, providers);
      auto ff_minus = process.ForcingFunction<MatrixPolicy>(
          phase_prefixes, state_parameter_indices, state_variable_indices, providers);

      ff_plus(state_parameters, vars_plus, forcing_plus);
      ff_minus(state_parameters, vars_minus, forcing_minus);

      for (std::size_t b = 0; b < num_blocks; ++b)
      {
        double h = std::max(std::abs(state_variables[b][j]) * eps, eps);
        for (std::size_t i = 0; i < num_vars; ++i)
        {
          double fd = (forcing_plus[b][i] - forcing_minus[b][i]) / (2.0 * h);
          // Jacobian stores -J, so analytical should be -fd
          double analytical;
          try
          {
            analytical = jacobian[b][i][j];
          }
          catch (...)
          {
            // Zero element in sparse matrix — fd should also be ~0
            if (std::abs(fd) > 1e-10)
              ADD_FAILURE() << "Missing Jacobian entry at block=" << b
                            << " row=" << i << " col=" << j << " fd=" << fd;
            continue;
          }

          double scale = std::max(std::abs(analytical), std::abs(fd));
          if (scale > 1e-20)
          {
            // analytical = -J[i,j], fd ≈ +J[i,j], so analytical + fd ≈ 0
            EXPECT_NEAR(analytical + fd, 0.0, scale * rel_tol)
                << "FD mismatch: block=" << b << " row=" << i << " col=" << j
                << " analytical(-J)=" << analytical << " fd(+J)=" << fd;
          }
        }
      }
    }
  }

  /// @brief Create a provider that varies linearly with given state variables:
  ///        value = base_value + sum(coeffs[k] * vars[dep_indices[k]])
  miam::AerosolPropertyProvider<MatrixPolicy> MakeLinearProvider(
      double base_value,
      const std::vector<std::size_t>& dep_indices,
      const std::vector<double>& coeffs)
  {
    miam::AerosolPropertyProvider<MatrixPolicy> provider;
    provider.dependent_variable_indices = dep_indices;
    provider.ComputeValue =
        [base_value, dep_indices, coeffs](
            const MatrixPolicy& params, const MatrixPolicy& vars, MatrixPolicy& result)
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
      {
        double val = base_value;
        for (std::size_t k = 0; k < dep_indices.size(); ++k)
          val += coeffs[k] * vars[row][dep_indices[k]];
        result[row][0] = val;
      }
    };
    provider.ComputeValueAndDerivatives =
        [base_value, dep_indices, coeffs](
            const MatrixPolicy& params, const MatrixPolicy& vars,
            MatrixPolicy& result, MatrixPolicy& partials)
    {
      for (std::size_t row = 0; row < result.NumRows(); ++row)
      {
        double val = base_value;
        for (std::size_t k = 0; k < dep_indices.size(); ++k)
        {
          val += coeffs[k] * vars[row][dep_indices[k]];
          partials[row][k] = coeffs[k];
        }
        result[row][0] = val;
      }
    };
    return provider;
  }
}  // namespace

// ---------------------------------------------------------------------------
// Test: Multiple phase instances (MODE1 + MODE2) — forcing
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, ForcingMultiplePhaseInstances)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".hlc"] = 2;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;
  state_variable_indices["MODE2.AQUEOUS.CO2_aq"] = 3;
  state_variable_indices["MODE2.AQUEOUS.H2O"] = 4;

  double r1 = 1.0e-6, N1 = 1.0e8, phi1 = 1.0e-6;
  double r2 = 5.0e-6, N2 = 1.0e7, phi2 = 1.0e-4;

  auto prov1 = MakeTestProviders("MODE1", r1, N1, phi1);
  auto prov2 = MakeTestProviders("MODE2", r2, N2, phi2);
  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers.insert(prov1.begin(), prov1.end());
  providers.insert(prov2.begin(), prov2.end());

  auto forcing_func = process.ForcingFunction<MatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, providers);

  double T = 298.15;
  double hlc = HLC_ref;
  double gas = 1.0e-3;
  double aq1 = 1.0e-5, solvent1 = 55000.0;
  double aq2 = 2.0e-5, solvent2 = 45000.0;

  MatrixPolicy params(1, 4);
  params[0][0] = hlc; params[0][1] = T;
  params[0][2] = hlc; params[0][3] = T;

  MatrixPolicy vars(1, 5);
  vars[0][0] = gas;
  vars[0][1] = aq1; vars[0][2] = solvent1;
  vars[0][3] = aq2; vars[0][4] = solvent2;

  MatrixPolicy forcing(1, 5, 0.0);
  forcing_func(params, vars, forcing);

  auto crp = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);

  double kc1 = crp.ComputeValue(r1, N1, T);
  double ke1 = kc1 / (hlc * miam::util::R_gas * T);
  double fv1 = solvent1 * Mw_solvent / rho_solvent;
  double net1 = phi1 * kc1 * gas - phi1 * ke1 * aq1 / fv1;

  double kc2 = crp.ComputeValue(r2, N2, T);
  double ke2 = kc2 / (hlc * miam::util::R_gas * T);
  double fv2 = solvent2 * Mw_solvent / rho_solvent;
  double net2 = phi2 * kc2 * gas - phi2 * ke2 * aq2 / fv2;

  // Gas forcing is sum of both instances (negative)
  EXPECT_NEAR(forcing[0][0], -(net1 + net2), std::abs(net1 + net2) * 1e-10);
  // Each aq mode gets its own forcing
  EXPECT_NEAR(forcing[0][1], net1, std::abs(net1) * 1e-10);
  EXPECT_NEAR(forcing[0][3], net2, std::abs(net2) * 1e-10);
  // Solvents unchanged
  EXPECT_NEAR(forcing[0][2], 0.0, 1e-30);
  EXPECT_NEAR(forcing[0][4], 0.0, 1e-30);
  // Mass conservation: gas + aq_mode1 + aq_mode2 = 0
  EXPECT_NEAR(forcing[0][0] + forcing[0][1] + forcing[0][3], 0.0, 1e-20);
}

// ---------------------------------------------------------------------------
// Test: Multiple phase instances — analytical Jacobian
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianMultiplePhaseInstances)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".hlc"] = 2;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;
  state_variable_indices["MODE2.AQUEOUS.CO2_aq"] = 3;
  state_variable_indices["MODE2.AQUEOUS.H2O"] = 4;

  double r1 = 1.0e-6, N1 = 1.0e8, phi1 = 1.0e-6;
  double r2 = 5.0e-6, N2 = 1.0e7, phi2 = 1.0e-4;

  auto prov1 = MakeTestProviders("MODE1", r1, N1, phi1);
  auto prov2 = MakeTestProviders("MODE2", r2, N2, phi2);
  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers.insert(prov1.begin(), prov1.end());
  providers.insert(prov2.begin(), prov2.end());

  auto jacobian = BuildJacobian(process, phase_prefixes, state_variable_indices, providers, 1);
  auto jac_func = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian, providers);

  double T = 298.15, hlc = HLC_ref;
  double gas = 1.0e-3;
  double aq1 = 1.0e-5, solvent1 = 55000.0;
  double aq2 = 2.0e-5, solvent2 = 45000.0;

  MatrixPolicy params(1, 4);
  params[0][0] = hlc; params[0][1] = T;
  params[0][2] = hlc; params[0][3] = T;

  MatrixPolicy vars(1, 5);
  vars[0][0] = gas;
  vars[0][1] = aq1; vars[0][2] = solvent1;
  vars[0][3] = aq2; vars[0][4] = solvent2;

  jac_func(params, vars, jacobian);

  auto crp = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);

  double kc1 = crp.ComputeValue(r1, N1, T);
  double ke1 = kc1 / (hlc * miam::util::R_gas * T);
  double fv1 = solvent1 * Mw_solvent / rho_solvent;

  double kc2 = crp.ComputeValue(r2, N2, T);
  double ke2 = kc2 / (hlc * miam::util::R_gas * T);
  double fv2 = solvent2 * Mw_solvent / rho_solvent;

  // -J[gas,gas] = phi1*kc1 + phi2*kc2  (both instances contribute)
  double j_gg = phi1 * kc1 + phi2 * kc2;
  EXPECT_NEAR(jacobian[0][0][0], j_gg, std::abs(j_gg) * 1e-10);

  // -J[gas, MODE1.aq] = -phi1*ke1/fv1
  EXPECT_NEAR(jacobian[0][0][1], -phi1 * ke1 / fv1, std::abs(phi1 * ke1 / fv1) * 1e-10);
  // -J[gas, MODE2.aq] = -phi2*ke2/fv2
  EXPECT_NEAR(jacobian[0][0][3], -phi2 * ke2 / fv2, std::abs(phi2 * ke2 / fv2) * 1e-10);

  // -J[MODE1.aq, gas] = -phi1*kc1
  EXPECT_NEAR(jacobian[0][1][0], -phi1 * kc1, std::abs(phi1 * kc1) * 1e-10);
  // -J[MODE2.aq, gas] = -phi2*kc2
  EXPECT_NEAR(jacobian[0][3][0], -phi2 * kc2, std::abs(phi2 * kc2) * 1e-10);

  // Cross-mode independence: MODE1 aq row should not depend on MODE2 aq variable
  EXPECT_THROW(jacobian[0][1][3], std::exception);
  EXPECT_THROW(jacobian[0][3][1], std::exception);

  // Mass conservation: sum of gas row and all aq rows for each column = 0
  for (std::size_t j = 0; j < 5; ++j)
  {
    double sum = 0.0;
    for (std::size_t i : { std::size_t(0), std::size_t(1), std::size_t(3) })
    {
      try { sum += jacobian[0][i][j]; } catch (...) {}
    }
    EXPECT_NEAR(sum, 0.0, 1e-15) << "Mass conservation violated for column " << j;
  }
}

// ---------------------------------------------------------------------------
// Test: Multiple phase instances — FD Jacobian check
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDMultiplePhaseInstances)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".hlc"] = 2;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;
  state_variable_indices["MODE2.AQUEOUS.CO2_aq"] = 3;
  state_variable_indices["MODE2.AQUEOUS.H2O"] = 4;

  double r1 = 1.0e-6, N1 = 1.0e8, phi1 = 1.0e-6;
  double r2 = 5.0e-6, N2 = 1.0e7, phi2 = 1.0e-4;

  auto prov1 = MakeTestProviders("MODE1", r1, N1, phi1);
  auto prov2 = MakeTestProviders("MODE2", r2, N2, phi2);
  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers.insert(prov1.begin(), prov1.end());
  providers.insert(prov2.begin(), prov2.end());

  MatrixPolicy params(1, 4);
  params[0][0] = HLC_ref; params[0][1] = 298.15;
  params[0][2] = HLC_ref; params[0][3] = 298.15;

  MatrixPolicy vars(1, 5);
  vars[0][0] = 1.0e-3;
  vars[0][1] = 1.0e-5; vars[0][2] = 55000.0;
  vars[0][3] = 2.0e-5; vars[0][4] = 45000.0;

  CheckFiniteDifferenceJacobian(process, phase_prefixes, state_parameter_indices,
                                state_variable_indices, providers, params, vars);
}

// ---------------------------------------------------------------------------
// Test: Multiple grid cells with varying conditions — FD Jacobian check
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDMultipleCellsVaryingConditions)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;

  auto providers = MakeTestProviders("MODE1", 2.0e-6, 5.0e8, 1.0e-4);

  std::size_t num_cells = 4;
  MatrixPolicy params(num_cells, 2);
  MatrixPolicy vars(num_cells, 3);

  double temperatures[] = { 270.0, 285.0, 298.15, 310.0 };
  double gas_concs[] = { 5.0e-4, 1.0e-3, 2.0e-3, 1.0e-2 };
  double aq_concs[] = { 1.0e-6, 1.0e-5, 5.0e-5, 1.0e-4 };
  double solvent_concs[] = { 55000.0, 50000.0, 45000.0, 40000.0 };

  for (std::size_t i = 0; i < num_cells; ++i)
  {
    params[i][0] = HLC_ref;
    params[i][1] = temperatures[i];
    vars[i][0] = gas_concs[i];
    vars[i][1] = aq_concs[i];
    vars[i][2] = solvent_concs[i];
  }

  CheckFiniteDifferenceJacobian(process, phase_prefixes, state_parameter_indices,
                                state_variable_indices, providers, params, vars);
}

// ---------------------------------------------------------------------------
// Test: Multiple grid cells + multiple instances — forcing analytical
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, ForcingMultiCellsMultiInstances)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".hlc"] = 2;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;
  state_variable_indices["MODE2.AQUEOUS.CO2_aq"] = 3;
  state_variable_indices["MODE2.AQUEOUS.H2O"] = 4;

  double r1 = 1.0e-6, N1 = 1.0e8, phi1 = 1.0e-6;
  double r2 = 5.0e-6, N2 = 1.0e7, phi2 = 1.0e-4;
  auto prov1 = MakeTestProviders("MODE1", r1, N1, phi1);
  auto prov2 = MakeTestProviders("MODE2", r2, N2, phi2);
  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers.insert(prov1.begin(), prov1.end());
  providers.insert(prov2.begin(), prov2.end());

  auto forcing_func = process.ForcingFunction<MatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, providers);

  std::size_t num_cells = 3;
  MatrixPolicy params(num_cells, 4);
  MatrixPolicy vars(num_cells, 5);

  double T = 298.15, hlc = HLC_ref;
  for (std::size_t c = 0; c < num_cells; ++c)
  {
    params[c][0] = hlc; params[c][1] = T;
    params[c][2] = hlc; params[c][3] = T;
    vars[c][0] = 1.0e-3 * (c + 1);
    vars[c][1] = 1.0e-5;   vars[c][2] = 55000.0;
    vars[c][3] = 2.0e-5;   vars[c][4] = 45000.0;
  }

  MatrixPolicy forcing(num_cells, 5, 0.0);
  forcing_func(params, vars, forcing);

  auto crp = miam::util::MakeCondensationRateProvider(D_g, alpha, Mw_gas);
  double kc1 = crp.ComputeValue(r1, N1, T);
  double ke1 = kc1 / (hlc * miam::util::R_gas * T);
  double fv1 = 55000.0 * Mw_solvent / rho_solvent;
  double kc2 = crp.ComputeValue(r2, N2, T);
  double ke2 = kc2 / (hlc * miam::util::R_gas * T);
  double fv2 = 45000.0 * Mw_solvent / rho_solvent;

  for (std::size_t c = 0; c < num_cells; ++c)
  {
    double gas_c = 1.0e-3 * (c + 1);
    double net1 = phi1 * kc1 * gas_c - phi1 * ke1 * 1.0e-5 / fv1;
    double net2 = phi2 * kc2 * gas_c - phi2 * ke2 * 2.0e-5 / fv2;
    EXPECT_NEAR(forcing[c][0], -(net1 + net2), std::abs(net1 + net2) * 1e-10)
        << "Gas forcing mismatch at cell " << c;
    EXPECT_NEAR(forcing[c][1], net1, std::abs(net1) * 1e-10)
        << "MODE1 aq forcing mismatch at cell " << c;
    EXPECT_NEAR(forcing[c][3], net2, std::abs(net2) * 1e-10)
        << "MODE2 aq forcing mismatch at cell " << c;
    EXPECT_NEAR(forcing[c][0] + forcing[c][1] + forcing[c][3], 0.0, 1e-20)
        << "Mass conservation at cell " << c;
  }
}

// ---------------------------------------------------------------------------
// Test: Multiple grid cells + multiple instances — FD Jacobian
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDMultiCellsMultiInstances)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  state_parameter_indices["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".hlc"] = 2;
  state_parameter_indices["MODE2.AQUEOUS." + process.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["CO2_g"] = 0;
  state_variable_indices["MODE1.AQUEOUS.CO2_aq"] = 1;
  state_variable_indices["MODE1.AQUEOUS.H2O"] = 2;
  state_variable_indices["MODE2.AQUEOUS.CO2_aq"] = 3;
  state_variable_indices["MODE2.AQUEOUS.H2O"] = 4;

  double r1 = 1.0e-6, N1 = 1.0e8, phi1 = 1.0e-6;
  double r2 = 5.0e-6, N2 = 1.0e7, phi2 = 1.0e-4;
  auto prov1 = MakeTestProviders("MODE1", r1, N1, phi1);
  auto prov2 = MakeTestProviders("MODE2", r2, N2, phi2);
  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers.insert(prov1.begin(), prov1.end());
  providers.insert(prov2.begin(), prov2.end());

  std::size_t num_cells = 3;
  MatrixPolicy params(num_cells, 4);
  MatrixPolicy vars(num_cells, 5);

  for (std::size_t c = 0; c < num_cells; ++c)
  {
    params[c][0] = HLC_ref; params[c][1] = 298.15;
    params[c][2] = HLC_ref; params[c][3] = 298.15;
    vars[c][0] = 1.0e-3 * (c + 1);
    vars[c][1] = 1.0e-5;   vars[c][2] = 55000.0;
    vars[c][3] = 2.0e-5;   vars[c][4] = 45000.0;
  }

  CheckFiniteDifferenceJacobian(process, phase_prefixes, state_parameter_indices,
                                state_variable_indices, providers, params, vars);
}

// ---------------------------------------------------------------------------
// Test: Two different gas species transferring into same condensed phase
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, ForcingMultipleTransferProcesses)
{
  auto gas_CO2 = micm::Species{ "CO2_g", { { "molecular weight [kg mol-1]", 0.044 } } };
  auto gas_SO2 = micm::Species{ "SO2_g", { { "molecular weight [kg mol-1]", 0.064 } } };
  auto aq_CO2 = micm::Species{ "CO2_aq", { { "molecular weight [kg mol-1]", 0.044 }, { "density [kg m-3]", 1800.0 } } };
  auto aq_SO2 = micm::Species{ "SO2_aq", { { "molecular weight [kg mol-1]", 0.064 }, { "density [kg m-3]", 1400.0 } } };
  auto solvent = micm::Species{ "H2O", { { "molecular weight [kg mol-1]", Mw_solvent }, { "density [kg m-3]", rho_solvent } } };
  auto aq_phase = micm::Phase{ "AQUEOUS", { { aq_CO2 }, { aq_SO2 }, { solvent } } };

  double HLC_CO2 = 3.4e-2, HLC_SO2 = 1.2;
  double D_CO2 = 1.5e-5, D_SO2 = 1.2e-5;

  auto proc_CO2 = HenryLawPhaseTransfer(
      [HLC_CO2](const micm::Conditions&) { return HLC_CO2; },
      gas_CO2, aq_CO2, solvent, aq_phase, D_CO2, alpha, 0.044, Mw_solvent, rho_solvent);
  auto proc_SO2 = HenryLawPhaseTransfer(
      [HLC_SO2](const micm::Conditions&) { return HLC_SO2; },
      gas_SO2, aq_SO2, solvent, aq_phase, D_SO2, alpha, 0.064, Mw_solvent, rho_solvent);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + proc_CO2.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + proc_CO2.uuid_ + ".temperature"] = 1;
  spi["MODE1.AQUEOUS." + proc_SO2.uuid_ + ".hlc"] = 2;
  spi["MODE1.AQUEOUS." + proc_SO2.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["SO2_g"] = 1;
  svi["MODE1.AQUEOUS.CO2_aq"] = 2;
  svi["MODE1.AQUEOUS.SO2_aq"] = 3;
  svi["MODE1.AQUEOUS.H2O"] = 4;

  double r = 2.0e-6, N = 5.0e8, phi = 1.0e-4;
  auto providers = MakeTestProviders("MODE1", r, N, phi);

  auto ff_CO2 = proc_CO2.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);
  auto ff_SO2 = proc_SO2.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);

  double T = 298.15;
  MatrixPolicy params(1, 4);
  params[0][0] = HLC_CO2; params[0][1] = T;
  params[0][2] = HLC_SO2; params[0][3] = T;

  double co2_g = 1.0e-3, so2_g = 5.0e-4;
  double co2_aq = 1.0e-5, so2_aq = 1.0e-4;
  double h2o = 55000.0;

  MatrixPolicy vars(1, 5);
  vars[0][0] = co2_g; vars[0][1] = so2_g;
  vars[0][2] = co2_aq; vars[0][3] = so2_aq; vars[0][4] = h2o;

  MatrixPolicy forcing(1, 5, 0.0);
  ff_CO2(params, vars, forcing);
  ff_SO2(params, vars, forcing);

  auto crp_CO2 = miam::util::MakeCondensationRateProvider(D_CO2, alpha, 0.044);
  auto crp_SO2 = miam::util::MakeCondensationRateProvider(D_SO2, alpha, 0.064);
  double fv = h2o * Mw_solvent / rho_solvent;

  double kc_co2 = crp_CO2.ComputeValue(r, N, T);
  double ke_co2 = kc_co2 / (HLC_CO2 * miam::util::R_gas * T);
  double net_co2 = phi * kc_co2 * co2_g - phi * ke_co2 * co2_aq / fv;

  double kc_so2 = crp_SO2.ComputeValue(r, N, T);
  double ke_so2 = kc_so2 / (HLC_SO2 * miam::util::R_gas * T);
  double net_so2 = phi * kc_so2 * so2_g - phi * ke_so2 * so2_aq / fv;

  EXPECT_NEAR(forcing[0][0], -net_co2, std::abs(net_co2) * 1e-10);
  EXPECT_NEAR(forcing[0][1], -net_so2, std::abs(net_so2) * 1e-10);
  EXPECT_NEAR(forcing[0][2], net_co2, std::abs(net_co2) * 1e-10);
  EXPECT_NEAR(forcing[0][3], net_so2, std::abs(net_so2) * 1e-10);

  // Per-species mass conservation
  EXPECT_NEAR(forcing[0][0] + forcing[0][2], 0.0, 1e-20);
  EXPECT_NEAR(forcing[0][1] + forcing[0][3], 0.0, 1e-20);
}

// ---------------------------------------------------------------------------
// Test: Two different transfer processes — FD Jacobian (combined)
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDMultipleTransferProcesses)
{
  auto gas_CO2 = micm::Species{ "CO2_g", { { "molecular weight [kg mol-1]", 0.044 } } };
  auto gas_SO2 = micm::Species{ "SO2_g", { { "molecular weight [kg mol-1]", 0.064 } } };
  auto aq_CO2 = micm::Species{ "CO2_aq", { { "molecular weight [kg mol-1]", 0.044 }, { "density [kg m-3]", 1800.0 } } };
  auto aq_SO2 = micm::Species{ "SO2_aq", { { "molecular weight [kg mol-1]", 0.064 }, { "density [kg m-3]", 1400.0 } } };
  auto solvent = micm::Species{ "H2O", { { "molecular weight [kg mol-1]", Mw_solvent }, { "density [kg m-3]", rho_solvent } } };
  auto aq_phase = micm::Phase{ "AQUEOUS", { { aq_CO2 }, { aq_SO2 }, { solvent } } };

  double HLC_CO2 = 3.4e-2, HLC_SO2 = 1.2;
  double D_CO2 = 1.5e-5, D_SO2 = 1.2e-5;

  auto proc_CO2 = HenryLawPhaseTransfer(
      [HLC_CO2](const micm::Conditions&) { return HLC_CO2; },
      gas_CO2, aq_CO2, solvent, aq_phase, D_CO2, alpha, 0.044, Mw_solvent, rho_solvent);
  auto proc_SO2 = HenryLawPhaseTransfer(
      [HLC_SO2](const micm::Conditions&) { return HLC_SO2; },
      gas_SO2, aq_SO2, solvent, aq_phase, D_SO2, alpha, 0.064, Mw_solvent, rho_solvent);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + proc_CO2.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + proc_CO2.uuid_ + ".temperature"] = 1;
  spi["MODE1.AQUEOUS." + proc_SO2.uuid_ + ".hlc"] = 2;
  spi["MODE1.AQUEOUS." + proc_SO2.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["SO2_g"] = 1;
  svi["MODE1.AQUEOUS.CO2_aq"] = 2;
  svi["MODE1.AQUEOUS.SO2_aq"] = 3;
  svi["MODE1.AQUEOUS.H2O"] = 4;

  double r = 2.0e-6, N = 5.0e8, phi = 1.0e-4;
  auto providers = MakeTestProviders("MODE1", r, N, phi);

  auto jacobian = BuildJacobian({ std::cref(proc_CO2), std::cref(proc_SO2) },
                                phase_prefixes, svi, providers, 1);

  auto jf_CO2 = proc_CO2.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, spi, svi, jacobian, providers);
  auto jf_SO2 = proc_SO2.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, spi, svi, jacobian, providers);

  double T = 298.15;
  MatrixPolicy params(1, 4);
  params[0][0] = HLC_CO2; params[0][1] = T;
  params[0][2] = HLC_SO2; params[0][3] = T;

  MatrixPolicy vars(1, 5);
  vars[0][0] = 1.0e-3;  vars[0][1] = 5.0e-4;
  vars[0][2] = 1.0e-5;  vars[0][3] = 1.0e-4;
  vars[0][4] = 55000.0;

  jf_CO2(params, vars, jacobian);
  jf_SO2(params, vars, jacobian);

  // FD check with central differences
  double eps = 1e-7;
  std::size_t nv = 5;
  for (std::size_t j = 0; j < nv; ++j)
  {
    MatrixPolicy vp(vars), vm(vars);
    double h = std::max(std::abs(vars[0][j]) * eps, eps);
    vp[0][j] += h;
    vm[0][j] -= h;

    MatrixPolicy fp(1, nv, 0.0), fm(1, nv, 0.0);
    auto fp_co2 = proc_CO2.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);
    auto fp_so2 = proc_SO2.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);
    auto fm_co2 = proc_CO2.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);
    auto fm_so2 = proc_SO2.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);
    fp_co2(params, vp, fp); fp_so2(params, vp, fp);
    fm_co2(params, vm, fm); fm_so2(params, vm, fm);

    for (std::size_t i = 0; i < nv; ++i)
    {
      double fd = (fp[0][i] - fm[0][i]) / (2.0 * h);
      double analytical;
      try { analytical = jacobian[0][i][j]; } catch (...) { continue; }
      double scale = std::max(std::abs(analytical), std::abs(fd));
      if (scale > 1e-20)
      {
        EXPECT_NEAR(analytical + fd, 0.0, scale * 1e-4)
            << "FD mismatch: row=" << i << " col=" << j
            << " analytical(-J)=" << analytical << " fd(+J)=" << fd;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Test: State-dependent providers (linear) — indirect Jacobian entries via FD
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDWithLinearProviders)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["MODE1.AQUEOUS.CO2_aq"] = 1;
  svi["MODE1.AQUEOUS.H2O"] = 2;
  svi["MODE1.AQUEOUS.EXTRA1"] = 3;
  svi["MODE1.AQUEOUS.EXTRA2"] = 4;

  // r_eff depends linearly on EXTRA1: r_eff = 1e-6 + 1e-8 * [EXTRA1]
  auto r_eff_prov = MakeLinearProvider(1.0e-6, { 3 }, { 1.0e-8 });
  // N depends linearly on EXTRA2: N = 1e8 + 1e6 * [EXTRA2]
  auto N_prov = MakeLinearProvider(1.0e8, { 4 }, { 1.0e6 });
  // phi depends linearly on H2O: phi = 1e-6 + 1e-10 * [H2O]
  auto phi_prov = MakeLinearProvider(1.0e-6, { 2 }, { 1.0e-10 });

  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers["MODE1"][miam::AerosolProperty::EffectiveRadius] = r_eff_prov;
  providers["MODE1"][miam::AerosolProperty::NumberConcentration] = N_prov;
  providers["MODE1"][miam::AerosolProperty::PhaseVolumeFraction] = phi_prov;

  MatrixPolicy params(1, 2);
  params[0][0] = HLC_ref;
  params[0][1] = 298.15;

  MatrixPolicy vars(1, 5);
  vars[0][0] = 1.0e-3;
  vars[0][1] = 1.0e-5;
  vars[0][2] = 55000.0;
  vars[0][3] = 10.0;
  vars[0][4] = 5.0;

  CheckFiniteDifferenceJacobian(process, phase_prefixes, spi, svi, providers, params, vars);
}

// ---------------------------------------------------------------------------
// Test: Linear providers — multi-cell FD check
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDLinearProvidersMultiCell)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["MODE1.AQUEOUS.CO2_aq"] = 1;
  svi["MODE1.AQUEOUS.H2O"] = 2;
  svi["MODE1.AQUEOUS.EXTRA1"] = 3;
  svi["MODE1.AQUEOUS.EXTRA2"] = 4;

  auto r_eff_prov = MakeLinearProvider(1.0e-6, { 3 }, { 1.0e-8 });
  auto N_prov = MakeLinearProvider(1.0e8, { 4 }, { 1.0e6 });
  auto phi_prov = MakeLinearProvider(1.0e-6, { 2 }, { 1.0e-10 });

  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers["MODE1"][miam::AerosolProperty::EffectiveRadius] = r_eff_prov;
  providers["MODE1"][miam::AerosolProperty::NumberConcentration] = N_prov;
  providers["MODE1"][miam::AerosolProperty::PhaseVolumeFraction] = phi_prov;

  std::size_t nc = 3;
  MatrixPolicy params(nc, 2);
  MatrixPolicy vars(nc, 5);
  for (std::size_t c = 0; c < nc; ++c)
  {
    params[c][0] = HLC_ref;
    params[c][1] = 298.15;
    vars[c][0] = 1.0e-3 * (c + 1);
    vars[c][1] = 1.0e-5 * (c + 1);
    vars[c][2] = 55000.0 - 5000.0 * c;
    vars[c][3] = 10.0 + 5.0 * c;
    vars[c][4] = 5.0 + 3.0 * c;
  }

  CheckFiniteDifferenceJacobian(process, phase_prefixes, spi, svi, providers, params, vars);
}

// ---------------------------------------------------------------------------
// Test: Accumulation — non-zero initial forcing and Jacobian
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, ForcingAccumulates)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["MODE1.AQUEOUS.CO2_aq"] = 1;
  svi["MODE1.AQUEOUS.H2O"] = 2;

  auto providers = MakeTestProviders("MODE1", 2.0e-6, 5.0e8, 1.0e-4);

  auto ff = process.ForcingFunction<MatrixPolicy>(phase_prefixes, spi, svi, providers);

  MatrixPolicy params(1, 2);
  params[0][0] = HLC_ref; params[0][1] = 298.15;

  MatrixPolicy vars(1, 3);
  vars[0][0] = 1.0e-3; vars[0][1] = 1.0e-5; vars[0][2] = 55000.0;

  MatrixPolicy forcing(1, 3, 0.0);
  ff(params, vars, forcing);
  double f0_gas = forcing[0][0];
  double f0_aq = forcing[0][1];

  // Second call should accumulate
  ff(params, vars, forcing);
  EXPECT_NEAR(forcing[0][0], 2.0 * f0_gas, std::abs(f0_gas) * 1e-10);
  EXPECT_NEAR(forcing[0][1], 2.0 * f0_aq, std::abs(f0_aq) * 1e-10);
}

TEST(HenryLawPhaseTransfer, JacobianAccumulates)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["MODE1.AQUEOUS.CO2_aq"] = 1;
  svi["MODE1.AQUEOUS.H2O"] = 2;

  auto providers = MakeTestProviders("MODE1", 2.0e-6, 5.0e8, 1.0e-4);

  auto jacobian = BuildJacobian(process, phase_prefixes, svi, providers, 1);
  auto jf = process.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, spi, svi, jacobian, providers);

  MatrixPolicy params(1, 2);
  params[0][0] = HLC_ref; params[0][1] = 298.15;

  MatrixPolicy vars(1, 3);
  vars[0][0] = 1.0e-3; vars[0][1] = 1.0e-5; vars[0][2] = 55000.0;

  jf(params, vars, jacobian);
  double j_gg_once = jacobian[0][0][0];

  jf(params, vars, jacobian);
  EXPECT_NEAR(jacobian[0][0][0], 2.0 * j_gg_once, std::abs(j_gg_once) * 1e-10);
}

// ---------------------------------------------------------------------------
// Test: Linear providers with multiple dependencies per provider
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDMultipleDepsPerProvider)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["MODE1.AQUEOUS.CO2_aq"] = 1;
  svi["MODE1.AQUEOUS.H2O"] = 2;
  svi["MODE1.AQUEOUS.NaCl"] = 3;
  svi["MODE1.AQUEOUS.NH4"] = 4;
  svi["MODE1.AQUEOUS.SO4"] = 5;

  // r_eff depends on NaCl AND NH4
  auto r_eff_prov = MakeLinearProvider(1.0e-6, { 3, 4 }, { 5.0e-9, 3.0e-9 });
  // N depends on SO4
  auto N_prov = MakeLinearProvider(1.0e8, { 5 }, { 2.0e6 });
  // phi depends on H2O and NaCl
  auto phi_prov = MakeLinearProvider(1.0e-6, { 2, 3 }, { 1.0e-10, 5.0e-8 });

  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers["MODE1"][miam::AerosolProperty::EffectiveRadius] = r_eff_prov;
  providers["MODE1"][miam::AerosolProperty::NumberConcentration] = N_prov;
  providers["MODE1"][miam::AerosolProperty::PhaseVolumeFraction] = phi_prov;

  MatrixPolicy params(1, 2);
  params[0][0] = HLC_ref;
  params[0][1] = 298.15;

  MatrixPolicy vars(1, 6);
  vars[0][0] = 1.0e-3;
  vars[0][1] = 1.0e-5;
  vars[0][2] = 55000.0;
  vars[0][3] = 0.5;
  vars[0][4] = 0.3;
  vars[0][5] = 0.1;

  CheckFiniteDifferenceJacobian(process, phase_prefixes, spi, svi, providers, params, vars);
}

// ---------------------------------------------------------------------------
// Test: Multiple instances + linear providers + multi-cell — the kitchen sink
// ---------------------------------------------------------------------------
TEST(HenryLawPhaseTransfer, JacobianFDKitchenSink)
{
  auto process = MakeTestProcess();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> spi;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".hlc"] = 0;
  spi["MODE1.AQUEOUS." + process.uuid_ + ".temperature"] = 1;
  spi["MODE2.AQUEOUS." + process.uuid_ + ".hlc"] = 2;
  spi["MODE2.AQUEOUS." + process.uuid_ + ".temperature"] = 3;

  std::unordered_map<std::string, std::size_t> svi;
  svi["CO2_g"] = 0;
  svi["MODE1.AQUEOUS.CO2_aq"] = 1;
  svi["MODE1.AQUEOUS.H2O"] = 2;
  svi["MODE1.AQUEOUS.NaCl"] = 3;
  svi["MODE2.AQUEOUS.CO2_aq"] = 4;
  svi["MODE2.AQUEOUS.H2O"] = 5;
  svi["MODE2.AQUEOUS.NaCl"] = 6;

  // MODE1: r_eff depends on NaCl(3), N constant, phi depends on H2O(2)
  auto r_eff_prov1 = MakeLinearProvider(1.0e-6, { 3 }, { 5.0e-9 });
  auto N_prov1 = MakeConstantProvider<MatrixPolicy>(1.0e8);
  auto phi_prov1 = MakeLinearProvider(1.0e-6, { 2 }, { 1.0e-10 });

  // MODE2: r_eff constant, N depends on NaCl(6), phi constant
  auto r_eff_prov2 = MakeConstantProvider<MatrixPolicy>(5.0e-6);
  auto N_prov2 = MakeLinearProvider(1.0e7, { 6 }, { 1.0e5 });
  auto phi_prov2 = MakeConstantProvider<MatrixPolicy>(1.0e-4);

  std::map<std::string, std::map<miam::AerosolProperty, miam::AerosolPropertyProvider<MatrixPolicy>>> providers;
  providers["MODE1"][miam::AerosolProperty::EffectiveRadius] = r_eff_prov1;
  providers["MODE1"][miam::AerosolProperty::NumberConcentration] = N_prov1;
  providers["MODE1"][miam::AerosolProperty::PhaseVolumeFraction] = phi_prov1;
  providers["MODE2"][miam::AerosolProperty::EffectiveRadius] = r_eff_prov2;
  providers["MODE2"][miam::AerosolProperty::NumberConcentration] = N_prov2;
  providers["MODE2"][miam::AerosolProperty::PhaseVolumeFraction] = phi_prov2;

  std::size_t nc = 2;
  MatrixPolicy params(nc, 4);
  MatrixPolicy vars(nc, 7);
  for (std::size_t c = 0; c < nc; ++c)
  {
    params[c][0] = HLC_ref; params[c][1] = 298.15;
    params[c][2] = HLC_ref; params[c][3] = 298.15;
    vars[c][0] = 1.0e-3 * (c + 1);
    vars[c][1] = 1.0e-5 * (c + 1); vars[c][2] = 55000.0 - 5000.0 * c;
    vars[c][3] = 0.5 + 0.1 * c;
    vars[c][4] = 2.0e-5 * (c + 1); vars[c][5] = 45000.0 + 3000.0 * c;
    vars[c][6] = 0.3 + 0.2 * c;
  }

  CheckFiniteDifferenceJacobian(process, phase_prefixes, spi, svi, providers, params, vars);
}

// ======================== Builder ========================

TEST(HenryLawPhaseTransferBuilder, BuildSuccess)
{
  miam::process::constant::HenrysLawConstantParameters hlc_params;
  hlc_params.HLC_ref_ = HLC_ref;
  hlc_params.C_ = 2400.0;
  hlc_params.T0_ = 298.15;
  miam::process::constant::HenrysLawConstant hlc(hlc_params);

  auto process = HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(MakeAqueousPhase())
      .SetGasSpecies(MakeGasSpecies())
      .SetCondensedSpecies(MakeCondensedSpecies())
      .SetSolvent(MakeSolvent())
      .SetHenrysLawConstant(hlc)
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
      .Build();

  EXPECT_DOUBLE_EQ(process.D_g_, D_g);
  EXPECT_DOUBLE_EQ(process.alpha_, alpha);
  EXPECT_DOUBLE_EQ(process.Mw_gas_, Mw_gas);
  EXPECT_DOUBLE_EQ(process.Mw_solvent_, Mw_solvent);
  EXPECT_DOUBLE_EQ(process.rho_solvent_, rho_solvent);
  EXPECT_EQ(process.gas_species_.name_, "CO2_g");
  EXPECT_FALSE(process.uuid_.empty());
}

TEST(HenryLawPhaseTransferBuilder, MissingCondensedPhaseThrows)
{
  miam::process::constant::HenrysLawConstantParameters hlc_params;
  hlc_params.HLC_ref_ = HLC_ref;
  miam::process::constant::HenrysLawConstant hlc(hlc_params);

  EXPECT_THROW(
      HenryLawPhaseTransferBuilder()
          .SetGasSpecies(MakeGasSpecies())
          .SetCondensedSpecies(MakeCondensedSpecies())
          .SetSolvent(MakeSolvent())
          .SetHenrysLawConstant(hlc)
          .SetDiffusionCoefficient(D_g)
          .SetAccommodationCoefficient(alpha)
          .Build(),
      std::runtime_error);
}

TEST(HenryLawPhaseTransferBuilder, MissingGasSpeciesThrows)
{
  miam::process::constant::HenrysLawConstantParameters hlc_params;
  hlc_params.HLC_ref_ = HLC_ref;
  miam::process::constant::HenrysLawConstant hlc(hlc_params);

  EXPECT_THROW(
      HenryLawPhaseTransferBuilder()
          .SetCondensedPhase(MakeAqueousPhase())
          .SetCondensedSpecies(MakeCondensedSpecies())
          .SetSolvent(MakeSolvent())
          .SetHenrysLawConstant(hlc)
          .SetDiffusionCoefficient(D_g)
          .SetAccommodationCoefficient(alpha)
          .Build(),
      std::runtime_error);
}

TEST(HenryLawPhaseTransferBuilder, MissingHLCThrows)
{
  EXPECT_THROW(
      HenryLawPhaseTransferBuilder()
          .SetCondensedPhase(MakeAqueousPhase())
          .SetGasSpecies(MakeGasSpecies())
          .SetCondensedSpecies(MakeCondensedSpecies())
          .SetSolvent(MakeSolvent())
          .SetDiffusionCoefficient(D_g)
          .SetAccommodationCoefficient(alpha)
          .Build(),
      std::runtime_error);
}

TEST(HenryLawPhaseTransferBuilder, MissingDiffusionCoefficientThrows)
{
  miam::process::constant::HenrysLawConstantParameters hlc_params;
  hlc_params.HLC_ref_ = HLC_ref;
  miam::process::constant::HenrysLawConstant hlc(hlc_params);

  EXPECT_THROW(
      HenryLawPhaseTransferBuilder()
          .SetCondensedPhase(MakeAqueousPhase())
          .SetGasSpecies(MakeGasSpecies())
          .SetCondensedSpecies(MakeCondensedSpecies())
          .SetSolvent(MakeSolvent())
          .SetHenrysLawConstant(hlc)
          .SetAccommodationCoefficient(alpha)
          .Build(),
      std::runtime_error);
}

TEST(HenryLawPhaseTransferBuilder, MissingAccommodationCoefficientThrows)
{
  miam::process::constant::HenrysLawConstantParameters hlc_params;
  hlc_params.HLC_ref_ = HLC_ref;
  miam::process::constant::HenrysLawConstant hlc(hlc_params);

  EXPECT_THROW(
      HenryLawPhaseTransferBuilder()
          .SetCondensedPhase(MakeAqueousPhase())
          .SetGasSpecies(MakeGasSpecies())
          .SetCondensedSpecies(MakeCondensedSpecies())
          .SetSolvent(MakeSolvent())
          .SetHenrysLawConstant(hlc)
          .SetDiffusionCoefficient(D_g)
          .Build(),
      std::runtime_error);
}

TEST(HenryLawPhaseTransferBuilder, BuiltProcessHLCWorks)
{
  miam::process::constant::HenrysLawConstantParameters hlc_params;
  hlc_params.HLC_ref_ = HLC_ref;
  hlc_params.C_ = 2400.0;
  hlc_params.T0_ = 298.15;
  miam::process::constant::HenrysLawConstant hlc(hlc_params);

  auto process = HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(MakeAqueousPhase())
      .SetGasSpecies(MakeGasSpecies())
      .SetCondensedSpecies(MakeCondensedSpecies())
      .SetSolvent(MakeSolvent())
      .SetHenrysLawConstant(hlc)
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
      .Build();

  micm::Conditions cond;
  cond.temperature_ = 298.15;
  double result = process.henry_law_constant_(cond);
  EXPECT_NEAR(result, HLC_ref, 1e-10);

  cond.temperature_ = 280.0;
  result = process.henry_law_constant_(cond);
  double expected = HLC_ref * std::exp(2400.0 * (1.0 / 280.0 - 1.0 / 298.15));
  EXPECT_NEAR(result, expected, expected * 1e-10);
}
