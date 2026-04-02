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
