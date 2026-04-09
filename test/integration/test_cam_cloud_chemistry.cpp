// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Integration test: CAM cloud aqueous-phase sulfate chemistry in MIAM.
//
// Built up incrementally to validate each layer of the mechanism:
//
//   Test 1: Single HLC (SO2 gas ⇌ aqueous) with mass conservation
//   Test 2: HLC + SO2 first dissociation (Ka1) + mass conservation
//   Test 3: All 3 HLCs + all dissociations + charge balance + mass conservation
//   Test 4: Full system with kinetic reactions
//   Test 5: FD Jacobian verification for each subsystem
//
// UNIT CONVENTIONS:
//   MIAM state variables:     mol/m³
//   Henry's Law constants:    mol/(m³·Pa)
//   Equilibrium K_miam:       K_lit / c_H2O  (for 1→2 dissociation)
//   Kinetic k_miam:           k_lit × c_H2O  (compensates solvent normalization)

#include <miam/miam.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>
#include <micm/util/jacobian_verification.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

using namespace micm;
using namespace miam;

using DenseMatrix = micm::Matrix<double>;
using SparseMatrixFD = micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>;

namespace
{
  constexpr double R_gas = miam::util::R_gas;

  // M/atm → mol m⁻³ Pa⁻¹
  constexpr double M_ATM_TO_MOL_M3_PA = 1000.0 / 101325.0;

  // Pure water concentration
  constexpr double c_H2O_M = 55.556;           // mol/L
  constexpr double C_H2O   = 55556.0;          // mol/m³ (state variable value for f_v=1)

  constexpr double T0 = 298.15;

  // Common solver helper: integrate and return convergence status
  template<typename SolverT, typename StateT>
  bool IntegrateDAE(SolverT& solver, StateT& state, double target_time,
                    double dt0, bool verbose = false)
  {
    double total_time = 0.0;
    double dt = dt0;
    while (total_time < target_time - 1.0e-10)
    {
      double step = std::min(dt, target_time - total_time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(step, state);
      if (result.state_ != SolverState::Converged)
      {
        if (verbose)
          std::cerr << "Solver failed at t=" << total_time << " s" << std::endl;
        return false;
      }
      total_time += step;
      // Ramp timestep
      if (total_time > 0.1 && dt < 0.1) dt = 0.1;
      if (total_time > 1.0 && dt < 1.0)  dt = 1.0;
      if (total_time > 10.0 && dt < 10.0)  dt = 10.0;
      if (total_time > 100.0 && dt < 100.0) dt = 100.0;
    }
    return true;
  }

  // FD Jacobian verification helpers (reused from test_jacobian_verification.cpp)
  struct IndexMaps
  {
    std::unordered_map<std::string, std::size_t> variable_indices;
    std::unordered_map<std::string, std::size_t> parameter_indices;
    std::size_t num_variables;
    std::size_t num_parameters;
  };

  IndexMaps BuildIndexMaps(const Model& model)
  {
    IndexMaps result;
    auto var_names = model.StateVariableNames();
    auto species_used = model.SpeciesUsed();
    for (const auto& name : species_used)
      var_names.insert(name);
    auto param_names = model.StateParameterNames();
    std::size_t idx = 0;
    for (const auto& name : var_names)
      result.variable_indices[name] = idx++;
    result.num_variables = idx;
    idx = 0;
    for (const auto& name : param_names)
      result.parameter_indices[name] = idx++;
    result.num_parameters = idx;
    return result;
  }

  void VerifyProcessJacobian(
      const Model& model, const IndexMaps& maps,
      const DenseMatrix& variables, const DenseMatrix& parameters,
      const std::vector<Conditions>& conditions,
      double atol = 1.0e-5, double rtol = 1.0e-4)
  {
    const std::size_t num_species = maps.num_variables;
    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    DenseMatrix params_copy(parameters);
    update_fn(conditions, params_copy);
    auto nz_elements = model.NonZeroJacobianElements(maps.variable_indices);
    auto builder = SparseMatrixFD::Create(num_species)
                       .SetNumberOfBlocks(variables.NumRows()).InitialValue(0.0);
    for (const auto& elem : nz_elements)
      builder = builder.WithElement(elem.first, elem.second);
    SparseMatrixFD analytical_jac(builder);
    auto jacobian_fn =
        model.JacobianFunction<DenseMatrix, SparseMatrixFD>(maps.parameter_indices, maps.variable_indices, analytical_jac);
    jacobian_fn(params_copy, variables, analytical_jac);
    auto forcing_fn = model.ForcingFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
    auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
    { forcing_fn(params_copy, vars, forcing); };
    auto fd_jac = FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);
    auto comparison = CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(
        analytical_jac, fd_jac, num_species, atol, rtol);
    EXPECT_TRUE(comparison.passed) << "Process Jacobian mismatch: row=" << comparison.worst_row
                                   << " col=" << comparison.worst_col
                                   << " analytical=" << comparison.worst_analytical
                                   << " fd=" << comparison.worst_fd;
    auto sparsity = CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrixFD>(
        analytical_jac, fd_jac, num_species);
    EXPECT_TRUE(sparsity.passed) << "Missing sparsity at row=" << sparsity.worst_row
                                 << " col=" << sparsity.worst_col
                                 << " fd_value=" << sparsity.worst_fd;
  }

  void VerifyConstraintJacobian(
      const Model& model, const IndexMaps& maps,
      const DenseMatrix& variables, const DenseMatrix& parameters,
      const std::vector<Conditions>& conditions,
      double atol = 1.0e-5, double rtol = 1.0e-4)
  {
    const std::size_t num_species = maps.num_variables;
    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    DenseMatrix params_copy(parameters);
    update_fn(conditions, params_copy);
    auto nz_elements = model.NonZeroConstraintJacobianElements(maps.variable_indices);
    auto builder = SparseMatrixFD::Create(num_species)
                       .SetNumberOfBlocks(variables.NumRows()).InitialValue(0.0);
    for (const auto& elem : nz_elements)
      builder = builder.WithElement(elem.first, elem.second);
    SparseMatrixFD analytical_jac(builder);
    auto jac_fn =
        model.ConstraintJacobianFunction<DenseMatrix, SparseMatrixFD>(maps.variable_indices, analytical_jac);
    jac_fn(variables, analytical_jac);
    auto residual_fn = model.ConstraintResidualFunction<DenseMatrix>(maps.variable_indices);
    auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
    { residual_fn(vars, forcing); };
    auto fd_jac = FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);
    auto comparison = CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(
        analytical_jac, fd_jac, num_species, atol, rtol);
    EXPECT_TRUE(comparison.passed) << "Constraint Jacobian mismatch: row=" << comparison.worst_row
                                   << " col=" << comparison.worst_col
                                   << " analytical=" << comparison.worst_analytical
                                   << " fd=" << comparison.worst_fd;
    auto sparsity = CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrixFD>(
        analytical_jac, fd_jac, num_species);
    EXPECT_TRUE(sparsity.passed) << "Missing constraint sparsity at row=" << sparsity.worst_row
                                 << " col=" << sparsity.worst_col
                                 << " fd_value=" << sparsity.worst_fd;
  }

  // Helper for variable index lookup
  template<typename StateT>
  std::size_t FindIdx(const StateT& state, const std::string& name)
  {
    auto it = std::find(state.variable_names_.begin(),
                        state.variable_names_.end(), name);
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  }
}  // namespace

// ════════════════════════════════════════════════════════════════════════
// TEST 1: Single HLC with mass conservation
//
// SO2(g) ⇌ SO2(aq)  with  [SO2_g] + [SO2_aq] = total
//
// Verifies that: (a) HLC equilibrium is reached, (b) mass is conserved,
// (c) gas depletes properly via LinearConstraint.
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step1_SingleHLC)
{
  double T = 280.0;
  double HLC_ref = 1.23 * M_ATM_TO_MOL_M3_PA;
  double C_hlc = 3120.0;

  auto so2_g  = Species{ "SO2" };
  auto so2_aq = Species{ "SO2_aq" };
  auto h2o    = Species{ "H2O",
      {{ "molecular weight [kg mol-1]", 0.018 },
       { "density [kg m-3]", 1000.0 }} };

  Phase gas_phase{ "GAS", { so2_g } };
  Phase aqueous_phase{ "AQUEOUS", { so2_aq, h2o } };
  auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

  auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(so2_g)
      .SetCondensedSpecies(so2_aq)
      .SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = HLC_ref, .C_ = C_hlc }))
      .SetMwSolvent(0.018)
      .SetRhoSolvent(1000.0)
      .Build();

  // Mass conservation: [SO2_g] + [SO2_aq] = total  (SO2_g algebraic)
  double gas0 = 3.01e-8;   // ~ 1 ppb
  double total_so2 = gas0;  // aq0 = 0

  auto mass_so2 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, so2_g)
      .AddTerm(gas_phase, so2_g, 1.0)
      .AddTerm(aqueous_phase, so2_aq, 1.0)
      .SetConstant(total_so2)
      .Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddConstraints(hl_so2, mass_so2);

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 70000.0;
  state.conditions_[0].CalculateIdealAirDensity();

  auto i_g  = FindIdx(state, "SO2");
  auto i_aq = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
  auto i_w  = FindIdx(state, "CLOUD.AQUEOUS.H2O");

  state.variables_[0][i_g]  = gas0;
  state.variables_[0][i_aq] = 0.0;
  state.variables_[0][i_w]  = C_H2O;
  cloud.SetDefaultParameters(state);

  // Analytical equilibrium: aq = α*g, g + aq = total → g = total/(1+α)
  double hlc_T = HLC_ref * std::exp(C_hlc * (1.0/T - 1.0/T0));
  double alpha = hlc_T * R_gas * T;  // f_v = 1
  double g_eq = total_so2 / (1.0 + alpha);
  double aq_eq = total_so2 * alpha / (1.0 + alpha);

  bool ok1 = IntegrateDAE(solver, state, 10.0, 0.01);
  ASSERT_TRUE(ok1);

  double g_f  = state.variables_[0][i_g];
  double aq_f = state.variables_[0][i_aq];

  // Mass conservation
  EXPECT_NEAR(g_f + aq_f, total_so2, 1e-15 * total_so2)
      << "Mass conservation violated";

  // HLC equilibrium
  EXPECT_NEAR(g_f, g_eq, 1e-6 * total_so2) << "SO2_g equilibrium";
  EXPECT_NEAR(aq_f, aq_eq, 1e-6 * total_so2) << "SO2_aq equilibrium";

  std::cout << "\n=== Step 1 PASSED ===" << std::endl;
  std::cout << "alpha=" << alpha << " g=" << g_f << " aq=" << aq_f
            << " total=" << g_f + aq_f << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// TEST 1b: Kw dissociation alone — pure water equilibrium
//
// H2O ⇌ H+ + OH-    (Kw, OHm algebraic)
// [H+] = [OH-]       (charge balance, Hp algebraic)
//
// Expected: [H+] = [OH-] = sqrt(Kw_miam * S²)
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step1b_KwOnly)
{
  double T = 280.0;

  auto hp   = Species{ "Hp" };
  auto ohm  = Species{ "OHm" };
  auto h2o  = Species{ "H2O",
      {{ "molecular weight [kg mol-1]", 0.018 },
       { "density [kg m-3]", 1000.0 }} };

  Phase gas_phase{ "GAS", {} };
  Phase aqueous_phase{ "AQUEOUS", { hp, ohm, h2o } };
  auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

  double Kw_miam = 1.0e-14 / (c_H2O_M * c_H2O_M);

  auto eq_kw = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ h2o })
      .SetProducts({ hp, ohm })
      .SetAlgebraicSpecies(ohm)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = Kw_miam }))
      .Build();

  // Charge balance: H+ = OH-
  auto charge = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, hp)
      .AddTerm(aqueous_phase, hp, 1.0)
      .AddTerm(aqueous_phase, ohm, -1.0)
      .SetConstant(0.0)
      .Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddConstraints(eq_kw, charge);

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 70000.0;

  auto i_hp  = FindIdx(state, "CLOUD.AQUEOUS.Hp");
  auto i_oh  = FindIdx(state, "CLOUD.AQUEOUS.OHm");
  auto i_w   = FindIdx(state, "CLOUD.AQUEOUS.H2O");

  // Expected: [H+] = sqrt(Kw_miam * S^2) = sqrt(1e-14) * 1000 = 1e-4 mol/m³
  // (converting from mol/L: pH=7 → [H+]=1e-7 M = 1e-4 mol/m³)
  double expected_hp = std::sqrt(Kw_miam) * C_H2O;

  state.variables_[0][i_hp] = 1e-4;
  state.variables_[0][i_oh] = 1e-4;
  state.variables_[0][i_w]  = C_H2O;
  cloud.SetDefaultParameters(state);

  std::cout << "\nKw_miam=" << Kw_miam << " expected [H+]=" << expected_hp << std::endl;

  bool ok = IntegrateDAE(solver, state, 1.0, 0.01);
  ASSERT_TRUE(ok) << "Solver failed for Kw-only system";

  double hp_f = state.variables_[0][i_hp];
  double oh_f = state.variables_[0][i_oh];

  std::cout << "H+=" << hp_f << " OH-=" << oh_f << " expected=" << expected_hp << std::endl;
  double pH = -std::log10(hp_f / 1000.0);
  std::cout << "pH=" << pH << std::endl;

  EXPECT_NEAR(hp_f, expected_hp, 0.01 * expected_hp) << "H+ equilibrium";
  EXPECT_NEAR(oh_f, expected_hp, 0.01 * expected_hp) << "OH- equilibrium";
  EXPECT_NEAR(hp_f, oh_f, 1e-10) << "Charge balance: H+ != OH-";

  std::cout << "=== Step 1b PASSED ===" << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// TEST 2: HLC + first dissociation + mass conservation + charge balance
//
// SO2(g) ⇌ SO2(aq)        [HLC]
// SO2(aq) ⇌ HSO3⁻ + H⁺   [Ka1]
// H2O ⇌ H⁺ + OH⁻         [Kw]
// [SO2_g] + [SO2_aq] + [HSO3⁻] = total_S   [mass conservation, SO2_g algebraic]
// [H⁺] = [OH⁻] + [HSO3⁻]                   [charge balance, H⁺ algebraic]
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step2_HLC_Plus_Dissociation)
{
  double T = 280.0;

  auto so2_g  = Species{ "SO2" };
  auto so2_aq = Species{ "SO2_aq" };
  auto hso3m  = Species{ "HSO3m" };
  auto hp     = Species{ "Hp" };
  auto ohm    = Species{ "OHm" };
  auto h2o    = Species{ "H2O",
      {{ "molecular weight [kg mol-1]", 0.018 },
       { "density [kg m-3]", 1000.0 }} };

  Phase gas_phase{ "GAS", { so2_g } };
  Phase aqueous_phase{ "AQUEOUS", { so2_aq, hso3m, hp, ohm, h2o } };
  auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

  double HLC_ref = 1.23 * M_ATM_TO_MOL_M3_PA;
  double C_hlc = 3120.0;

  auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(so2_g)
      .SetCondensedSpecies(so2_aq)
      .SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = HLC_ref, .C_ = C_hlc }))
      .SetMwSolvent(0.018)
      .SetRhoSolvent(1000.0)
      .Build();

  // Kw: H2O ⇌ H+ + OH-
  auto eq_kw = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ h2o })
      .SetProducts({ hp, ohm })
      .SetAlgebraicSpecies(ohm)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M) }))
      .Build();

  // Ka1: SO2_aq ⇌ HSO3⁻ + H⁺
  auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ so2_aq })
      .SetProducts({ hso3m, hp })
      .SetAlgebraicSpecies(hso3m)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
      .Build();

  // S-budget: [SO2_g] + [SO2_aq] + [HSO3−] = total_S  (SO2_g algebraic)
  double gas0 = 3.01e-8;
  double total_S = gas0;

  auto mass_S = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, so2_g)
      .AddTerm(gas_phase, so2_g, 1.0)
      .AddTerm(aqueous_phase, so2_aq, 1.0)
      .AddTerm(aqueous_phase, hso3m, 1.0)
      .SetConstant(total_S)
      .Build();

  // Charge balance: [H+] = [OH-] + [HSO3-]  (H+ algebraic)
  auto charge = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, hp)
      .AddTerm(aqueous_phase, hp,   1.0)
      .AddTerm(aqueous_phase, ohm, -1.0)
      .AddTerm(aqueous_phase, hso3m, -1.0)
      .SetConstant(0.0)
      .Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddConstraints(hl_so2, eq_kw, eq_ka1, mass_S, charge);

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 70000.0;
  state.conditions_[0].CalculateIdealAirDensity();

  auto i_g  = FindIdx(state, "SO2");
  auto i_aq = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
  auto i_hs = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
  auto i_hp = FindIdx(state, "CLOUD.AQUEOUS.Hp");
  auto i_oh = FindIdx(state, "CLOUD.AQUEOUS.OHm");
  auto i_w  = FindIdx(state, "CLOUD.AQUEOUS.H2O");

  // Initial conditions: compute self-consistent equilibrium values
  // by damped fixed-point iteration on [H+].
  //
  // Equations:
  //   HLC:    SO2_aq = α · SO2_g
  //   Ka1:    HSO3⁻ = Ka1 · SO2_aq · S / H⁺
  //   Kw:     OH⁻   = Kw · S² / H⁺
  //   Mass:   SO2_g + SO2_aq + HSO3⁻ = total_S
  //   Charge: H⁺ = OH⁻ + HSO3⁻
  //
  // Substituting everything into the charge balance gives a fixed-point
  // on H⁺ that converges with simple damping.
  double hlc_T = (1.23 * M_ATM_TO_MOL_M3_PA) * std::exp(3120.0 * (1.0/T - 1.0/T0));
  double alpha = hlc_T * R_gas * T;  // f_v = 1
  double Ka1_T = (1.7e-2 / c_H2O_M) * std::exp(2090.0 * (1.0/T0 - 1.0/T));
  double Kw_T  = 1.0e-14 / (c_H2O_M * c_H2O_M);

  double hp_iter = std::sqrt(Kw_T) * C_H2O;  // start at pH 7
  double so2_g_ic = 0, so2_aq_ic = 0, hso3_ic = 0, oh_ic = 0;
  for (int it = 0; it < 50; ++it)
  {
    oh_ic     = Kw_T * C_H2O * C_H2O / hp_iter;
    so2_g_ic  = total_S / (1.0 + alpha + Ka1_T * alpha * C_H2O / hp_iter);
    so2_aq_ic = alpha * so2_g_ic;
    hso3_ic   = Ka1_T * so2_aq_ic * C_H2O / hp_iter;
    double hp_new = oh_ic + hso3_ic;
    if (std::abs(hp_new - hp_iter) < 1e-15 * hp_iter) break;
    hp_iter = 0.5 * (hp_iter + hp_new);  // damped to avoid oscillation
  }

  state.variables_[0][i_g]  = so2_g_ic;
  state.variables_[0][i_aq] = so2_aq_ic;
  state.variables_[0][i_hs] = hso3_ic;
  state.variables_[0][i_hp] = hp_iter;
  state.variables_[0][i_oh] = oh_ic;
  state.variables_[0][i_w]  = C_H2O;
  cloud.SetDefaultParameters(state);

  std::cout << "Step 2 initial: SO2_g=" << so2_g_ic << " SO2_aq=" << so2_aq_ic
            << " HSO3-=" << hso3_ic << " H+=" << hp_iter
            << " OH-=" << oh_ic << std::endl;

  // FD Jacobian check for this system
  {
    auto maps = BuildIndexMaps(model);
    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("SO2")] = so2_g_ic;
    variables[0][maps.variable_indices.at("CLOUD.AQUEOUS.SO2_aq")] = so2_aq_ic;
    variables[0][maps.variable_indices.at("CLOUD.AQUEOUS.HSO3m")] = hso3_ic;
    variables[0][maps.variable_indices.at("CLOUD.AQUEOUS.Hp")] = hp_iter;
    variables[0][maps.variable_indices.at("CLOUD.AQUEOUS.OHm")] = oh_ic;
    variables[0][maps.variable_indices.at("CLOUD.AQUEOUS.H2O")] = C_H2O;
    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = T;
    conditions[0].pressure_ = 70000.0;
    std::cout << "Checking Step 2 constraint Jacobian..." << std::endl;
    VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
    std::cout << "Step 2 constraint Jacobian OK" << std::endl;
  }

  // Try a single step first
  solver.UpdateStateParameters(state);
  auto result_dbg = solver.Solve(0.001, state);
  std::cout << "After single 0.001s step: converged=" << (result_dbg.state_ == SolverState::Converged)
            << " SO2_g=" << state.variables_[0][i_g] << " SO2_aq=" << state.variables_[0][i_aq]
            << " HSO3-=" << state.variables_[0][i_hs] << " H+=" << state.variables_[0][i_hp]
            << " OH-=" << state.variables_[0][i_oh] << " H2O=" << state.variables_[0][i_w]
            << std::endl;

  bool ok2 = IntegrateDAE(solver, state, 10.0, 0.01);
  ASSERT_TRUE(ok2);

  double g_f  = state.variables_[0][i_g];
  double aq_f = state.variables_[0][i_aq];
  double hs_f = state.variables_[0][i_hs];
  double hp_f = state.variables_[0][i_hp];
  double oh_f = state.variables_[0][i_oh];

  // Sulfur mass conservation
  double total_S_f = g_f + aq_f + hs_f;
  EXPECT_NEAR(total_S_f, total_S, 1e-10 * total_S) << "S budget violated";

  // Charge balance
  double cb_err = std::abs(hp_f - oh_f - hs_f);
  EXPECT_LT(cb_err, 1e-8 * hp_f) << "Charge balance violated: H+=" << hp_f
      << " OH-=" << oh_f << " HSO3-=" << hs_f;

  // All positive
  EXPECT_GT(g_f, 0) << "SO2_g negative";
  EXPECT_GT(aq_f, 0) << "SO2_aq negative";
  EXPECT_GT(hs_f, 0) << "HSO3m negative";
  EXPECT_GT(hp_f, 0) << "Hp negative";
  EXPECT_GT(oh_f, 0) << "OHm negative";

  double pH = -std::log10(hp_f / 1000.0);  // convert mol/m³ to mol/L
  std::cout << "\n=== Step 2 PASSED ===" << std::endl;
  std::cout << "pH=" << pH << " SO2_g=" << g_f << " SO2_aq=" << aq_f
            << " HSO3-=" << hs_f << " total_S=" << total_S_f << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// TEST 3: Full equilibrium system (3 HLCs + 3 dissociations + charge
//         balance + 3 mass conservations) - NO kinetic reactions yet
//
// Algebraic variables:
//   SO2_g    (from S mass conservation)
//   H2O2_g   (from H2O2 mass conservation)
//   O3_g     (from O3 mass conservation)
//   SO2_aq   (from HLC_SO2)
//   H2O2_aq  (from HLC_H2O2)
//   O3_aq    (from HLC_O3)
//   OHm      (from Kw)
//   HSO3m    (from Ka1)
//   SO3mm    (from Ka2)
//   Hp       (from charge balance)
//
// Differential variables:
//   SO4mm    (only differential species - no source yet, stays at initial)
//   H2O      (solvent, effectively constant)
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step3_FullEquilibrium)
{
  double T = 280.0;

  auto so2_g   = Species{ "SO2" };
  auto h2o2_g  = Species{ "H2O2" };
  auto o3_g    = Species{ "O3" };
  auto so2_aq  = Species{ "SO2_aq" };
  auto h2o2_aq = Species{ "H2O2_aq" };
  auto o3_aq   = Species{ "O3_aq" };
  auto hp      = Species{ "Hp" };
  auto ohm     = Species{ "OHm" };
  auto hso3m   = Species{ "HSO3m" };
  auto so3mm   = Species{ "SO3mm" };
  auto so4mm   = Species{ "SO4mm" };
  auto h2o     = Species{ "H2O",
      {{ "molecular weight [kg mol-1]", 0.018 },
       { "density [kg m-3]", 1000.0 }} };

  Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g } };
  Phase aqueous_phase{ "AQUEOUS", {
      h2o, so2_aq, h2o2_aq, o3_aq,
      hp, ohm, hso3m, so3mm, so4mm } };
  auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

  // --- Henry's Law constraints ---
  auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(so2_g).SetCondensedSpecies(so2_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto hl_h2o2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(h2o2_g).SetCondensedSpecies(h2o2_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 7.4e4 * M_ATM_TO_MOL_M3_PA, .C_ = 6621.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto hl_o3 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(o3_g).SetCondensedSpecies(o3_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 1.15e-2 * M_ATM_TO_MOL_M3_PA, .C_ = 2560.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  // --- Dissociation equilibria ---
  auto eq_kw = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ h2o }).SetProducts({ hp, ohm })
      .SetAlgebraicSpecies(ohm).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M) }))
      .Build();

  auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ so2_aq }).SetProducts({ hso3m, hp })
      .SetAlgebraicSpecies(hso3m).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
      .Build();

  auto eq_ka2 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ hso3m }).SetProducts({ so3mm, hp })
      .SetAlgebraicSpecies(so3mm).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 6.0e-8 / c_H2O_M, .C_ = 1120.0 }))
      .Build();

  // --- Mass conservation (gas species become algebraic) ---
  // S budget: [SO2_g] + [SO2_aq] + [HSO3-] + [SO3--] = total_S
  double gas0_so2  = 3.01e-8;
  double gas0_h2o2 = 3.01e-8;
  double gas0_o3   = 1.50e-6;

  auto mass_S = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, so2_g)
      .AddTerm(gas_phase, so2_g, 1.0)
      .AddTerm(aqueous_phase, so2_aq, 1.0)
      .AddTerm(aqueous_phase, hso3m, 1.0)
      .AddTerm(aqueous_phase, so3mm, 1.0)
      .SetConstant(gas0_so2)
      .Build();

  auto mass_H2O2 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, h2o2_g)
      .AddTerm(gas_phase, h2o2_g, 1.0)
      .AddTerm(aqueous_phase, h2o2_aq, 1.0)
      .SetConstant(gas0_h2o2)
      .Build();

  auto mass_O3 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, o3_g)
      .AddTerm(gas_phase, o3_g, 1.0)
      .AddTerm(aqueous_phase, o3_aq, 1.0)
      .SetConstant(gas0_o3)
      .Build();

  // --- Charge balance: [H+] = [OH-] + [HSO3-] + 2[SO3--] + 2[SO4--] ---
  auto charge = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, hp)
      .AddTerm(aqueous_phase, hp,     1.0)
      .AddTerm(aqueous_phase, ohm,   -1.0)
      .AddTerm(aqueous_phase, hso3m, -1.0)
      .AddTerm(aqueous_phase, so3mm, -2.0)
      .AddTerm(aqueous_phase, so4mm, -2.0)
      .SetConstant(0.0)
      .Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddConstraints(hl_so2, hl_h2o2, hl_o3,
                       eq_kw, eq_ka1, eq_ka2,
                       mass_S, mass_H2O2, mass_O3, charge);

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 70000.0;
  state.conditions_[0].CalculateIdealAirDensity();

  auto i_so2_g_   = FindIdx(state, "SO2");
  auto i_h2o2_g_  = FindIdx(state, "H2O2");
  auto i_o3_g_    = FindIdx(state, "O3");
  auto i_h2o_     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so2_aq_  = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
  auto i_h2o2_aq_ = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
  auto i_o3_aq_   = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
  auto i_hp_      = FindIdx(state, "CLOUD.AQUEOUS.Hp");
  auto i_ohm_     = FindIdx(state, "CLOUD.AQUEOUS.OHm");
  auto i_hso3m_   = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
  auto i_so3mm_   = FindIdx(state, "CLOUD.AQUEOUS.SO3mm");
  auto i_so4mm_   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");

  // Initial conditions: compute self-consistent equilibrium values
  // by damped fixed-point iteration on [H+].
  double so4mm0 = 1.0;  // 1 mM in droplet (background sulfate)

  double hlc_SO2_T  = (1.23 * M_ATM_TO_MOL_M3_PA)  * std::exp(3120.0 * (1.0/T - 1.0/T0));
  double hlc_H2O2_T = (7.4e4 * M_ATM_TO_MOL_M3_PA) * std::exp(6621.0 * (1.0/T - 1.0/T0));
  double hlc_O3_T   = (1.15e-2 * M_ATM_TO_MOL_M3_PA) * std::exp(2560.0 * (1.0/T - 1.0/T0));
  double alpha_SO2  = hlc_SO2_T * R_gas * T;
  double alpha_H2O2 = hlc_H2O2_T * R_gas * T;
  double alpha_O3   = hlc_O3_T * R_gas * T;
  double Ka1_T = (1.7e-2 / c_H2O_M) * std::exp(2090.0 * (1.0/T0 - 1.0/T));
  double Ka2_T = (6.0e-8 / c_H2O_M) * std::exp(1120.0 * (1.0/T0 - 1.0/T));
  double Kw_T  = 1.0e-14 / (c_H2O_M * c_H2O_M);

  // H2O2 and O3 have no dissociation — simple HLC split
  double ic_h2o2_g  = gas0_h2o2 / (1.0 + alpha_H2O2);
  double ic_h2o2_aq = alpha_H2O2 * ic_h2o2_g;
  double ic_o3_g    = gas0_o3 / (1.0 + alpha_O3);
  double ic_o3_aq   = alpha_O3 * ic_o3_g;

  // Iterate on [H+] for SO2 equilibria + charge balance
  double hp_ic = 2.0 * so4mm0;  // charge balance dominated by SO4
  double ic_so2_g = 0, ic_so2_aq = 0, ic_hso3m = 0, ic_so3mm = 0, ic_ohm = 0;
  for (int it = 0; it < 100; ++it)
  {
    ic_ohm    = Kw_T * C_H2O * C_H2O / hp_ic;
    double f  = 1.0 + alpha_SO2
              + Ka1_T * alpha_SO2 * C_H2O / hp_ic
              + Ka2_T * Ka1_T * alpha_SO2 * C_H2O * C_H2O / (hp_ic * hp_ic);
    ic_so2_g  = gas0_so2 / f;
    ic_so2_aq = alpha_SO2 * ic_so2_g;
    ic_hso3m  = Ka1_T * ic_so2_aq * C_H2O / hp_ic;
    ic_so3mm  = Ka2_T * ic_hso3m * C_H2O / hp_ic;
    double hp_new = ic_ohm + ic_hso3m + 2.0 * ic_so3mm + 2.0 * so4mm0;
    if (std::abs(hp_new - hp_ic) < 1e-15 * hp_ic) break;
    hp_ic = 0.5 * (hp_ic + hp_new);
  }

  state.variables_[0][i_so2_g_]   = ic_so2_g;
  state.variables_[0][i_h2o2_g_]  = ic_h2o2_g;
  state.variables_[0][i_o3_g_]    = ic_o3_g;
  state.variables_[0][i_h2o_]     = C_H2O;
  state.variables_[0][i_so2_aq_]  = ic_so2_aq;
  state.variables_[0][i_h2o2_aq_] = ic_h2o2_aq;
  state.variables_[0][i_o3_aq_]   = ic_o3_aq;
  state.variables_[0][i_hp_]      = hp_ic;
  state.variables_[0][i_ohm_]     = ic_ohm;
  state.variables_[0][i_hso3m_]   = ic_hso3m;
  state.variables_[0][i_so3mm_]   = ic_so3mm;
  state.variables_[0][i_so4mm_]   = so4mm0;
  cloud.SetDefaultParameters(state);

  bool ok3 = IntegrateDAE(solver, state, 10.0, 0.01);
  ASSERT_TRUE(ok3);

  double so2_g_f   = state.variables_[0][i_so2_g_];
  double h2o2_g_f  = state.variables_[0][i_h2o2_g_];
  double o3_g_f    = state.variables_[0][i_o3_g_];
  double so2_aq_f  = state.variables_[0][i_so2_aq_];
  double h2o2_aq_f = state.variables_[0][i_h2o2_aq_];
  double o3_aq_f   = state.variables_[0][i_o3_aq_];
  double hp_f      = state.variables_[0][i_hp_];
  double ohm_f     = state.variables_[0][i_ohm_];
  double hso3m_f   = state.variables_[0][i_hso3m_];
  double so3mm_f   = state.variables_[0][i_so3mm_];
  double so4mm_f   = state.variables_[0][i_so4mm_];

  // Mass conservation checks
  EXPECT_NEAR(so2_g_f + so2_aq_f + hso3m_f + so3mm_f, gas0_so2,
              1e-10 * gas0_so2) << "S budget violated";
  EXPECT_NEAR(h2o2_g_f + h2o2_aq_f, gas0_h2o2,
              1e-10 * gas0_h2o2) << "H2O2 budget violated";
  EXPECT_NEAR(o3_g_f + o3_aq_f, gas0_o3,
              1e-10 * gas0_o3) << "O3 budget violated";

  // Charge balance
  double cb = hp_f - ohm_f - hso3m_f - 2*so3mm_f - 2*so4mm_f;
  EXPECT_NEAR(cb, 0.0, 1e-8 * hp_f) << "Charge balance violated";

  // SO4 should be unchanged (no kinetics)
  EXPECT_NEAR(so4mm_f, 1.0, 1e-6) << "SO4 should not change without kinetics";

  // All positive
  EXPECT_GT(so2_g_f, 0);
  EXPECT_GT(hp_f, 0);
  EXPECT_GT(hso3m_f, 0);

  double pH = -std::log10(hp_f / 1000.0);
  std::cout << "\n=== Step 3 PASSED ===" << std::endl;
  std::cout << "pH=" << pH << std::endl;
  std::cout << "SO2: g=" << so2_g_f << " aq=" << so2_aq_f << " HSO3-=" << hso3m_f
            << " SO3--=" << so3mm_f << " total=" << so2_g_f + so2_aq_f + hso3m_f + so3mm_f
            << std::endl;
  std::cout << "H2O2: g=" << h2o2_g_f << " aq=" << h2o2_aq_f << std::endl;
  std::cout << "O3: g=" << o3_g_f << " aq=" << o3_aq_f << std::endl;
  std::cout << "SO4--=" << so4mm_f << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// TEST 4: Full system with kinetic reactions
//
// Everything from Step 3, plus three kinetic S(IV)→S(VI) oxidation
// reactions. The mass conservation for S now includes SO4:
//   [SO2_g] + [SO2_aq] + [HSO3-] + [SO3--] + [SO4--] = total_S
//
// H2O2 is consumed by R1 and O3 by R2/R3, so their mass conservation
// totals will decrease — tracked by including the product in the budget.
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step4_FullSystemWithKinetics)
{
  double T = 280.0;

  auto so2_g   = Species{ "SO2" };
  auto h2o2_g  = Species{ "H2O2" };
  auto o3_g    = Species{ "O3" };
  auto so2_aq  = Species{ "SO2_aq" };
  auto h2o2_aq = Species{ "H2O2_aq" };
  auto o3_aq   = Species{ "O3_aq" };
  auto hp      = Species{ "Hp" };
  auto ohm     = Species{ "OHm" };
  auto hso3m   = Species{ "HSO3m" };
  auto so3mm   = Species{ "SO3mm" };
  auto so4mm   = Species{ "SO4mm" };
  auto h2o     = Species{ "H2O",
      {{ "molecular weight [kg mol-1]", 0.018 },
       { "density [kg m-3]", 1000.0 }} };

  Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g } };
  Phase aqueous_phase{ "AQUEOUS", {
      h2o, so2_aq, h2o2_aq, o3_aq,
      hp, ohm, hso3m, so3mm, so4mm } };
  auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

  // --- Henry's Law constraints ---
  auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(so2_g).SetCondensedSpecies(so2_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto hl_h2o2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(h2o2_g).SetCondensedSpecies(h2o2_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 7.4e4 * M_ATM_TO_MOL_M3_PA, .C_ = 6621.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto hl_o3 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(o3_g).SetCondensedSpecies(o3_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 1.15e-2 * M_ATM_TO_MOL_M3_PA, .C_ = 2560.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  // --- Dissociation equilibria ---
  auto eq_kw = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ h2o }).SetProducts({ hp, ohm })
      .SetAlgebraicSpecies(ohm).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M) }))
      .Build();

  auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ so2_aq }).SetProducts({ hso3m, hp })
      .SetAlgebraicSpecies(hso3m).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
      .Build();

  auto eq_ka2 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ hso3m }).SetProducts({ so3mm, hp })
      .SetAlgebraicSpecies(so3mm).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 6.0e-8 / c_H2O_M, .C_ = 1120.0 }))
      .Build();

  // --- Mass conservation ---
  // S budget: total S is conserved (all oxidation states tracked)
  double gas0_so2  = 3.01e-8;
  double gas0_h2o2 = 3.01e-8;
  double gas0_o3   = 1.50e-6;
  double so4mm0    = 1.0;  // 1 μM = 1.0 mol/m³
  double total_S   = gas0_so2 + so4mm0;

  auto mass_S = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, so2_g)
      .AddTerm(gas_phase, so2_g, 1.0)
      .AddTerm(aqueous_phase, so2_aq, 1.0)
      .AddTerm(aqueous_phase, hso3m, 1.0)
      .AddTerm(aqueous_phase, so3mm, 1.0)
      .AddTerm(aqueous_phase, so4mm, 1.0)
      .SetConstant(total_S)
      .Build();

  auto mass_H2O2 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, h2o2_g)
      .AddTerm(gas_phase, h2o2_g, 1.0)
      .AddTerm(aqueous_phase, h2o2_aq, 1.0)
      .SetConstant(gas0_h2o2)
      .Build();

  auto mass_O3 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, o3_g)
      .AddTerm(gas_phase, o3_g, 1.0)
      .AddTerm(aqueous_phase, o3_aq, 1.0)
      .SetConstant(gas0_o3)
      .Build();

  // Charge balance
  auto charge = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, hp)
      .AddTerm(aqueous_phase, hp,     1.0)
      .AddTerm(aqueous_phase, ohm,   -1.0)
      .AddTerm(aqueous_phase, hso3m, -1.0)
      .AddTerm(aqueous_phase, so3mm, -2.0)
      .AddTerm(aqueous_phase, so4mm, -2.0)
      .SetConstant(0.0)
      .Build();

  // --- Kinetic reactions ---
  // R1: HSO3⁻ + H2O2(aq) → SO4²⁻ + H2O + H⁺
  auto rxn1 = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ hso3m, h2o2_aq })
      .SetProducts({ so4mm, h2o, hp })
      .SetSolvent(h2o)
      .SetRateConstant([](const Conditions& c) -> double {
          return c_H2O_M * 7.45e7 *
                 std::exp(-4430.0 * (1.0/c.temperature_ - 1.0/298.0));
      })
      .Build();

  // R2: HSO3⁻ + O3(aq) → SO4²⁻ + H⁺
  auto rxn2 = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ hso3m, o3_aq })
      .SetProducts({ so4mm, hp })
      .SetSolvent(h2o)
      .SetRateConstant([](const Conditions& c) -> double {
          return c_H2O_M * 3.75e5 *
                 std::exp(-5530.0 * (1.0/c.temperature_ - 1.0/298.0));
      })
      .Build();

  // R3: SO3²⁻ + O3(aq) → SO4²⁻
  auto rxn3 = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ so3mm, o3_aq })
      .SetProducts({ so4mm })
      .SetSolvent(h2o)
      .SetRateConstant([](const Conditions& c) -> double {
          return c_H2O_M * 1.59e9 *
                 std::exp(-5280.0 * (1.0/c.temperature_ - 1.0/298.0));
      })
      .Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddProcesses({ rxn1, rxn2, rxn3 });
  model.AddConstraints(hl_so2, hl_h2o2, hl_o3,
                       eq_kw, eq_ka1, eq_ka2,
                       mass_S, mass_H2O2, mass_O3, charge);

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 70000.0;
  state.conditions_[0].CalculateIdealAirDensity();

  auto i_so2_g_   = FindIdx(state, "SO2");
  auto i_h2o2_g_  = FindIdx(state, "H2O2");
  auto i_o3_g_    = FindIdx(state, "O3");
  auto i_h2o_     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so2_aq_  = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
  auto i_h2o2_aq_ = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
  auto i_o3_aq_   = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
  auto i_hp_      = FindIdx(state, "CLOUD.AQUEOUS.Hp");
  auto i_ohm_     = FindIdx(state, "CLOUD.AQUEOUS.OHm");
  auto i_hso3m_   = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
  auto i_so3mm_   = FindIdx(state, "CLOUD.AQUEOUS.SO3mm");
  auto i_so4mm_   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");

  // Initial conditions: compute self-consistent equilibrium values
  // Same iteration as Step 3, but with SO4mm in the S-budget
  double hlc_SO2_T  = (1.23 * M_ATM_TO_MOL_M3_PA)  * std::exp(3120.0 * (1.0/T - 1.0/T0));
  double hlc_H2O2_T = (7.4e4 * M_ATM_TO_MOL_M3_PA) * std::exp(6621.0 * (1.0/T - 1.0/T0));
  double hlc_O3_T   = (1.15e-2 * M_ATM_TO_MOL_M3_PA) * std::exp(2560.0 * (1.0/T - 1.0/T0));
  double alpha_SO2  = hlc_SO2_T * R_gas * T;
  double alpha_H2O2 = hlc_H2O2_T * R_gas * T;
  double alpha_O3   = hlc_O3_T * R_gas * T;
  double Ka1_T = (1.7e-2 / c_H2O_M) * std::exp(2090.0 * (1.0/T0 - 1.0/T));
  double Ka2_T = (6.0e-8 / c_H2O_M) * std::exp(1120.0 * (1.0/T0 - 1.0/T));
  double Kw_T  = 1.0e-14 / (c_H2O_M * c_H2O_M);

  double ic_h2o2_g  = gas0_h2o2 / (1.0 + alpha_H2O2);
  double ic_h2o2_aq = alpha_H2O2 * ic_h2o2_g;
  double ic_o3_g    = gas0_o3 / (1.0 + alpha_O3);
  double ic_o3_aq   = alpha_O3 * ic_o3_g;

  // Note: In Step 4 the S-budget includes SO4: total_S = gas0 + so4mm0
  // The "non-SO4" sulfur is total_S - so4mm0 = gas0_so2
  double s4_budget = gas0_so2;  // S(IV) portion = total_S - so4mm0
  double hp_ic = 2.0 * so4mm0;  // charge balance dominated by SO4
  double ic_so2_g = 0, ic_so2_aq = 0, ic_hso3m = 0, ic_so3mm = 0, ic_ohm = 0;
  for (int it = 0; it < 100; ++it)
  {
    ic_ohm    = Kw_T * C_H2O * C_H2O / hp_ic;
    double f  = 1.0 + alpha_SO2
              + Ka1_T * alpha_SO2 * C_H2O / hp_ic
              + Ka2_T * Ka1_T * alpha_SO2 * C_H2O * C_H2O / (hp_ic * hp_ic);
    ic_so2_g  = s4_budget / f;
    ic_so2_aq = alpha_SO2 * ic_so2_g;
    ic_hso3m  = Ka1_T * ic_so2_aq * C_H2O / hp_ic;
    ic_so3mm  = Ka2_T * ic_hso3m * C_H2O / hp_ic;
    double hp_new = ic_ohm + ic_hso3m + 2.0 * ic_so3mm + 2.0 * so4mm0;
    if (std::abs(hp_new - hp_ic) < 1e-15 * hp_ic) break;
    hp_ic = 0.5 * (hp_ic + hp_new);
  }

  state.variables_[0][i_so2_g_]   = ic_so2_g;
  state.variables_[0][i_h2o2_g_]  = ic_h2o2_g;
  state.variables_[0][i_o3_g_]    = ic_o3_g;
  state.variables_[0][i_h2o_]     = C_H2O;
  state.variables_[0][i_so2_aq_]  = ic_so2_aq;
  state.variables_[0][i_h2o2_aq_] = ic_h2o2_aq;
  state.variables_[0][i_o3_aq_]   = ic_o3_aq;
  state.variables_[0][i_hp_]      = hp_ic;
  state.variables_[0][i_ohm_]     = ic_ohm;
  state.variables_[0][i_hso3m_]   = ic_hso3m;
  state.variables_[0][i_so3mm_]   = ic_so3mm;
  state.variables_[0][i_so4mm_]   = so4mm0;
  cloud.SetDefaultParameters(state);

  state.PrintHeader();
  state.PrintState(0);

  double total_time = 0.0;
  double dt = 0.001;
  bool converged = true;

  while (total_time < 1800.0 - 1.0e-10)
  {
    double step = std::min(dt, 1800.0 - total_time);
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(step, state);
    if (result.state_ != SolverState::Converged)
    {
      std::cerr << "Solver failed at t=" << total_time << " s" << std::endl;
      converged = false;
      break;
    }
    total_time += step;
    if (total_time > 0.01 && dt < 0.01) dt = 0.01;
    if (total_time > 0.1 && dt < 0.1) dt = 0.1;
    if (total_time > 1.0 && dt < 1.0)  dt = 1.0;
    if (total_time > 10.0 && dt < 10.0)  dt = 10.0;
    if (total_time > 100.0 && dt < 100.0) dt = 100.0;
  }

  state.PrintState(static_cast<int>(total_time));
  ASSERT_TRUE(converged) << "DAE solver failed to converge for full system";

  double so2_g_f   = state.variables_[0][i_so2_g_];
  double h2o2_g_f  = state.variables_[0][i_h2o2_g_];
  double o3_g_f    = state.variables_[0][i_o3_g_];
  double so2_aq_f  = state.variables_[0][i_so2_aq_];
  double h2o2_aq_f = state.variables_[0][i_h2o2_aq_];
  double o3_aq_f   = state.variables_[0][i_o3_aq_];
  double hp_f      = state.variables_[0][i_hp_];
  double ohm_f     = state.variables_[0][i_ohm_];
  double hso3m_f   = state.variables_[0][i_hso3m_];
  double so3mm_f   = state.variables_[0][i_so3mm_];
  double so4mm_f   = state.variables_[0][i_so4mm_];

  // Total S conservation (SO2_g + SO2_aq + HSO3- + SO3-- + SO4-- = total_S)
  double total_S_f = so2_g_f + so2_aq_f + hso3m_f + so3mm_f + so4mm_f;
  EXPECT_NEAR(total_S_f, total_S, 1e-6 * total_S) << "S budget violated";

  // Charge balance
  double cb = hp_f - ohm_f - hso3m_f - 2*so3mm_f - 2*so4mm_f;
  EXPECT_NEAR(cb, 0.0, 0.01 * hp_f) << "Charge balance violated";

  // SO4 should increase
  EXPECT_GT(so4mm_f, so4mm0) << "SO4 must increase from S(IV) oxidation";

  // Non-negative concentrations
  for (std::size_t v = 0; v < state.variables_.NumColumns(); ++v)
  {
    EXPECT_GE(state.variables_[0][v], -1.0e-10)
        << state.variable_names_[v] << " = " << state.variables_[0][v];
  }

  double pH = (hp_f > 0) ? -std::log10(hp_f / 1000.0) : 99.0;
  double air = state.conditions_[0].air_density_;
  std::cout << "\n=== Step 4 PASSED ===" << std::endl;
  std::cout << "pH=" << pH << std::endl;
  std::cout << "SO2(g)=" << so2_g_f/air*1e9 << " ppb" << std::endl;
  std::cout << "H2O2(g)=" << h2o2_g_f/air*1e9 << " ppb" << std::endl;
  std::cout << "O3(g)=" << o3_g_f/air*1e9 << " ppb" << std::endl;
  std::cout << "SO4=" << so4mm_f << " mol/m³ (" << so4mm_f/1000*1e6 << " μM)" << std::endl;
  std::cout << "S budget: " << total_S_f << " (expected " << total_S << ")" << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// TEST 5: Finite-difference Jacobian verification for the full system
//
// Verifies both process (kinetic) and constraint Jacobians at a
// physically realistic state point.
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step5_JacobianVerification)
{
  double T = 280.0;

  auto so2_g   = Species{ "SO2" };
  auto h2o2_g  = Species{ "H2O2" };
  auto o3_g    = Species{ "O3" };
  auto so2_aq  = Species{ "SO2_aq" };
  auto h2o2_aq = Species{ "H2O2_aq" };
  auto o3_aq   = Species{ "O3_aq" };
  auto hp      = Species{ "Hp" };
  auto ohm     = Species{ "OHm" };
  auto hso3m   = Species{ "HSO3m" };
  auto so3mm   = Species{ "SO3mm" };
  auto so4mm   = Species{ "SO4mm" };
  auto h2o     = Species{ "H2O",
      {{ "molecular weight [kg mol-1]", 0.018 },
       { "density [kg m-3]", 1000.0 }} };

  Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g } };
  Phase aqueous_phase{ "AQUEOUS", {
      h2o, so2_aq, h2o2_aq, o3_aq,
      hp, ohm, hso3m, so3mm, so4mm } };
  auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

  // Same constraints as Step 4
  auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(so2_g).SetCondensedSpecies(so2_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto hl_h2o2 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(h2o2_g).SetCondensedSpecies(h2o2_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 7.4e4 * M_ATM_TO_MOL_M3_PA, .C_ = 6621.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto hl_o3 = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(o3_g).SetCondensedSpecies(o3_aq).SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant({
          .HLC_ref_ = 1.15e-2 * M_ATM_TO_MOL_M3_PA, .C_ = 2560.0 }))
      .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

  auto eq_kw = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ h2o }).SetProducts({ hp, ohm })
      .SetAlgebraicSpecies(ohm).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M) }))
      .Build();

  auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ so2_aq }).SetProducts({ hso3m, hp })
      .SetAlgebraicSpecies(hso3m).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
      .Build();

  auto eq_ka2 = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase).SetReactants({ hso3m }).SetProducts({ so3mm, hp })
      .SetAlgebraicSpecies(so3mm).SetSolvent(h2o)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant({
          .A_ = 6.0e-8 / c_H2O_M, .C_ = 1120.0 }))
      .Build();

  double total_S = 3.01e-8 + 1.0;
  auto mass_S = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, so2_g)
      .AddTerm(gas_phase, so2_g, 1.0)
      .AddTerm(aqueous_phase, so2_aq, 1.0)
      .AddTerm(aqueous_phase, hso3m, 1.0)
      .AddTerm(aqueous_phase, so3mm, 1.0)
      .AddTerm(aqueous_phase, so4mm, 1.0)
      .SetConstant(total_S)
      .Build();

  auto mass_H2O2 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, h2o2_g)
      .AddTerm(gas_phase, h2o2_g, 1.0)
      .AddTerm(aqueous_phase, h2o2_aq, 1.0)
      .SetConstant(3.01e-8)
      .Build();

  auto mass_O3 = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, o3_g)
      .AddTerm(gas_phase, o3_g, 1.0)
      .AddTerm(aqueous_phase, o3_aq, 1.0)
      .SetConstant(1.50e-6)
      .Build();

  auto charge = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, hp)
      .AddTerm(aqueous_phase, hp,     1.0)
      .AddTerm(aqueous_phase, ohm,   -1.0)
      .AddTerm(aqueous_phase, hso3m, -1.0)
      .AddTerm(aqueous_phase, so3mm, -2.0)
      .AddTerm(aqueous_phase, so4mm, -2.0)
      .SetConstant(0.0)
      .Build();

  // Kinetic reactions
  auto rxn1 = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase).SetReactants({ hso3m, h2o2_aq })
      .SetProducts({ so4mm, h2o, hp }).SetSolvent(h2o)
      .SetRateConstant([](const Conditions& c) -> double {
          return c_H2O_M * 7.45e7 * std::exp(-4430.0 * (1.0/c.temperature_ - 1.0/298.0));
      }).Build();

  auto rxn2 = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase).SetReactants({ hso3m, o3_aq })
      .SetProducts({ so4mm, hp }).SetSolvent(h2o)
      .SetRateConstant([](const Conditions& c) -> double {
          return c_H2O_M * 3.75e5 * std::exp(-5530.0 * (1.0/c.temperature_ - 1.0/298.0));
      }).Build();

  auto rxn3 = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase).SetReactants({ so3mm, o3_aq })
      .SetProducts({ so4mm }).SetSolvent(h2o)
      .SetRateConstant([](const Conditions& c) -> double {
          return c_H2O_M * 1.59e9 * std::exp(-5280.0 * (1.0/c.temperature_ - 1.0/298.0));
      }).Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddProcesses({ rxn1, rxn2, rxn3 });
  model.AddConstraints(hl_so2, hl_h2o2, hl_o3,
                       eq_kw, eq_ka1, eq_ka2,
                       mass_S, mass_H2O2, mass_O3, charge);

  auto maps = BuildIndexMaps(model);

  // Test at two physically reasonable state points
  DenseMatrix variables(2, maps.num_variables, 0.0);
  auto set_var = [&](int block, const std::string& name, double val) {
    variables[block][maps.variable_indices.at(name)] = val;
  };

  // State point 1: near-initial conditions
  set_var(0, "SO2", 2.5e-8);
  set_var(0, "H2O2", 2.0e-8);
  set_var(0, "O3", 1.5e-6);
  set_var(0, "CLOUD.AQUEOUS.H2O", C_H2O);
  set_var(0, "CLOUD.AQUEOUS.SO2_aq", 1.0e-6);
  set_var(0, "CLOUD.AQUEOUS.H2O2_aq", 0.1);
  set_var(0, "CLOUD.AQUEOUS.O3_aq", 5.0e-7);
  set_var(0, "CLOUD.AQUEOUS.Hp", 5.0e-3);
  set_var(0, "CLOUD.AQUEOUS.OHm", 2.0e-6);
  set_var(0, "CLOUD.AQUEOUS.HSO3m", 3.0e-3);
  set_var(0, "CLOUD.AQUEOUS.SO3mm", 3.0e-5);
  set_var(0, "CLOUD.AQUEOUS.SO4mm", 1.0);

  // State point 2: somewhat evolved
  set_var(1, "SO2", 1.0e-8);
  set_var(1, "H2O2", 1.0e-8);
  set_var(1, "O3", 1.4e-6);
  set_var(1, "CLOUD.AQUEOUS.H2O", C_H2O);
  set_var(1, "CLOUD.AQUEOUS.SO2_aq", 5.0e-7);
  set_var(1, "CLOUD.AQUEOUS.H2O2_aq", 0.05);
  set_var(1, "CLOUD.AQUEOUS.O3_aq", 3.0e-7);
  set_var(1, "CLOUD.AQUEOUS.Hp", 8.0e-3);
  set_var(1, "CLOUD.AQUEOUS.OHm", 1.0e-6);
  set_var(1, "CLOUD.AQUEOUS.HSO3m", 1.5e-3);
  set_var(1, "CLOUD.AQUEOUS.SO3mm", 1.0e-5);
  set_var(1, "CLOUD.AQUEOUS.SO4mm", 1.5);

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = T;
  conditions[0].pressure_ = 70000.0;
  conditions[1].temperature_ = 290.0;
  conditions[1].pressure_ = 80000.0;

  std::cout << "\n=== Step 5: FD Jacobian Verification ===" << std::endl;

  // Verify process (kinetic) Jacobian
  std::cout << "Checking process Jacobian..." << std::endl;
  VerifyProcessJacobian(model, maps, variables, parameters, conditions);

  // Verify constraint Jacobian
  std::cout << "Checking constraint Jacobian..." << std::endl;
  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);

  std::cout << "=== Step 5 PASSED ===" << std::endl;
}
