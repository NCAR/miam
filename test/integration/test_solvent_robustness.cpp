// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Solvent-robustness stress tests for 3D-model readiness.
//
// Characterizes solver behavior as cloud water (solvent) approaches zero.
// All grid cells always run chemistry (cannot skip per-cell); the only
// mitigation for solvent→0 issues will be rate damping.
//
// Phase A: Individual process/constraint sweeps vs. solvent concentration
// Phase B: Full CAM cloud mechanism with vanishing cloud water
// Phase C: Multi-cell vertical column with mixed cloudy/clear cells
//
// Every test logs solver statistics (steps, rejections) to support
// future solver optimization (4- vs 6-stage DAE, tolerance tuning).

#include <miam/miam.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace micm;
using namespace miam;

using DenseMatrix = micm::Matrix<double>;

namespace
{
  constexpr double R_gas = miam::util::R_gas;
  constexpr double M_ATM_TO_MOL_M3_PA = 1000.0 / 101325.0;
  constexpr double c_H2O_M = 55.556;
  constexpr double Mw_water = 0.018;
  constexpr double rho_water = 1000.0;
  constexpr double T0 = 298.15;

  // Normal cloud water content
  constexpr double C_H2O_NORMAL = 0.017;  // mol/m³ air (~ 0.3 g/m³)

  // Solvent sweep levels (mol/m³ air)
  const std::vector<double> SOLVENT_SWEEP = {
      0.017, 1e-3, 1e-6, 1e-9, 1e-12, 0.0
  };

  // Extended sweep without the hard-zero for tests that can't handle exactly 0
  const std::vector<double> SOLVENT_SWEEP_NONZERO = {
      0.017, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15
  };

  // ── Accumulated solver statistics ──
  struct AccumulatedStats
  {
    uint64_t total_steps{};
    uint64_t total_accepted{};
    uint64_t total_rejected{};
    uint64_t total_function_calls{};
    uint64_t total_decompositions{};
    uint64_t outer_calls{};
    bool converged{ true };
    SolverState last_failure_state{ SolverState::Converged };
  };

  // ── Integrate with stats accumulation ──
  template <typename SolverT, typename StateT>
  AccumulatedStats IntegrateWithStats(SolverT& solver, StateT& state,
                                      double target_time, double dt0)
  {
    AccumulatedStats acc;
    double total_time = 0.0;
    double dt = dt0;
    while (total_time < target_time - 1.0e-10)
    {
      double step = std::min(dt, target_time - total_time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(step, state);
      acc.total_steps += result.stats_.number_of_steps_;
      acc.total_accepted += result.stats_.accepted_;
      acc.total_rejected += result.stats_.rejected_;
      acc.total_function_calls += result.stats_.function_calls_;
      acc.total_decompositions += result.stats_.decompositions_;
      acc.outer_calls++;
      if (result.state_ != SolverState::Converged)
      {
        acc.converged = false;
        acc.last_failure_state = result.state_;
        break;
      }
      total_time += step;
      if (total_time > 0.01 && dt < 0.01) dt = 0.01;
      if (total_time > 0.1 && dt < 0.1) dt = 0.1;
      if (total_time > 1.0 && dt < 1.0)  dt = 1.0;
      if (total_time > 10.0 && dt < 10.0)  dt = 10.0;
      if (total_time > 100.0 && dt < 100.0) dt = 100.0;
    }
    return acc;
  }

  // ── Print sweep table header ──
  void PrintSweepHeader(const std::string& test_name)
  {
    std::cout << "\n┌─── " << test_name << " ───" << std::endl;
    std::cout << std::setw(14) << "H2O [mol/m³]"
              << std::setw(12) << "converged"
              << std::setw(10) << "steps"
              << std::setw(10) << "rejected"
              << std::setw(10) << "f_calls"
              << std::setw(14) << "failure_mode"
              << std::endl;
    std::cout << std::string(70, '-') << std::endl;
  }

  // ── Print one sweep row ──
  void PrintSweepRow(double h2o, const AccumulatedStats& s)
  {
    std::cout << std::scientific << std::setprecision(1) << std::setw(14) << h2o
              << std::setw(12) << (s.converged ? "YES" : "NO")
              << std::dec << std::setw(10) << s.total_steps
              << std::setw(10) << s.total_rejected
              << std::setw(10) << s.total_function_calls
              << std::setw(14)
              << (s.converged ? "-" : SolverStateToString(s.last_failure_state))
              << std::endl;
  }

  // ── Variable index lookup ──
  template <typename StateT>
  std::size_t FindIdx(const StateT& state, const std::string& name)
  {
    auto it = std::find(state.variable_names_.begin(),
                        state.variable_names_.end(), name);
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  }

  // ── Check if any variable is NaN or Inf ──
  template <typename StateT>
  bool HasNanOrInf(const StateT& state, std::size_t cell = 0)
  {
    for (std::size_t v = 0; v < state.variables_.NumColumns(); ++v)
    {
      double val = state.variables_[cell][v];
      if (std::isnan(val) || std::isinf(val)) return true;
    }
    return false;
  }

  // ══════════════════════════════════════════════════════════
  // Builder helpers for the full CAM cloud mechanism
  // ══════════════════════════════════════════════════════════

  struct CamCloudSpecies
  {
    Species so2_g{ "SO2" };
    Species h2o2_g{ "H2O2" };
    Species o3_g{ "O3" };
    Species so2_aq{ "SO2_aq" };
    Species h2o2_aq{ "H2O2_aq" };
    Species o3_aq{ "O3_aq" };
    Species hp{ "Hp" };
    Species ohm{ "OHm" };
    Species hso3m{ "HSO3m" };
    Species so3mm{ "SO3mm" };
    Species so4mm{ "SO4mm" };
    Species so2oohm{ "SO2OOHm" };
    Species h2o{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g } };
    Phase aqueous_phase{ "AQUEOUS", {
        h2o, so2_aq, h2o2_aq, o3_aq,
        hp, ohm, hso3m, so3mm, so4mm, so2oohm } };
  };

  // Build the full model (HLCs + equilibria + mass/charge + kinetics)
  auto BuildCamCloudModel(CamCloudSpecies& sp)
  {
    auto cloud = representation::UniformSection{ "CLOUD", { sp.aqueous_phase } };

    auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
        .SetGasSpecies(sp.so2_g).SetCondensedSpecies(sp.so2_aq).SetSolvent(sp.h2o)
        .SetCondensedPhase(sp.aqueous_phase)
        .SetHenryLawConstant(process::constant::HenrysLawConstant({
            .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
        .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

    auto hl_h2o2 = constraint::HenryLawEquilibriumConstraintBuilder()
        .SetGasSpecies(sp.h2o2_g).SetCondensedSpecies(sp.h2o2_aq).SetSolvent(sp.h2o)
        .SetCondensedPhase(sp.aqueous_phase)
        .SetHenryLawConstant(process::constant::HenrysLawConstant({
            .HLC_ref_ = 7.4e4 * M_ATM_TO_MOL_M3_PA, .C_ = 6621.0 }))
        .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

    auto hl_o3 = constraint::HenryLawEquilibriumConstraintBuilder()
        .SetGasSpecies(sp.o3_g).SetCondensedSpecies(sp.o3_aq).SetSolvent(sp.h2o)
        .SetCondensedPhase(sp.aqueous_phase)
        .SetHenryLawConstant(process::constant::HenrysLawConstant({
            .HLC_ref_ = 1.15e-2 * M_ATM_TO_MOL_M3_PA, .C_ = 2560.0 }))
        .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

    auto eq_kw = constraint::DissolvedEquilibriumConstraintBuilder()
        .SetPhase(sp.aqueous_phase).SetReactants({ sp.h2o }).SetProducts({ sp.hp, sp.ohm })
        .SetAlgebraicSpecies(sp.ohm).SetSolvent(sp.h2o)
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M), .C_ = 6710.0 }))
        .Build();

    auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
        .SetPhase(sp.aqueous_phase).SetReactants({ sp.so2_aq }).SetProducts({ sp.hso3m, sp.hp })
        .SetAlgebraicSpecies(sp.hso3m).SetSolvent(sp.h2o)
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
        .Build();

    auto eq_ka2 = constraint::DissolvedEquilibriumConstraintBuilder()
        .SetPhase(sp.aqueous_phase).SetReactants({ sp.hso3m }).SetProducts({ sp.so3mm, sp.hp })
        .SetAlgebraicSpecies(sp.so3mm).SetSolvent(sp.h2o)
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 6.0e-8 / c_H2O_M, .C_ = 1120.0 }))
        .Build();

    auto mass_S = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(sp.gas_phase, sp.so2_g)
        .AddTerm(sp.gas_phase, sp.so2_g, 1.0)
        .AddTerm(sp.aqueous_phase, sp.so2_aq, 1.0)
        .AddTerm(sp.aqueous_phase, sp.hso3m, 1.0)
        .AddTerm(sp.aqueous_phase, sp.so3mm, 1.0)
        .AddTerm(sp.aqueous_phase, sp.so4mm, 1.0)
        .AddTerm(sp.aqueous_phase, sp.so2oohm, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto mass_H2O2 = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(sp.gas_phase, sp.h2o2_g)
        .AddTerm(sp.gas_phase, sp.h2o2_g, 1.0)
        .AddTerm(sp.aqueous_phase, sp.h2o2_aq, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto mass_O3 = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(sp.gas_phase, sp.o3_g)
        .AddTerm(sp.gas_phase, sp.o3_g, 1.0)
        .AddTerm(sp.aqueous_phase, sp.o3_aq, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto charge = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(sp.aqueous_phase, sp.hp)
        .AddTerm(sp.aqueous_phase, sp.hp, 1.0)
        .AddTerm(sp.aqueous_phase, sp.ohm, -1.0)
        .AddTerm(sp.aqueous_phase, sp.hso3m, -1.0)
        .AddTerm(sp.aqueous_phase, sp.so3mm, -2.0)
        .AddTerm(sp.aqueous_phase, sp.so4mm, -2.0)
        .AddTerm(sp.aqueous_phase, sp.so2oohm, -1.0)
        .SetConstant(0.0)
        .Build();

    auto rxn1a = process::DissolvedReversibleReactionBuilder()
        .SetPhase(sp.aqueous_phase)
        .SetReactants({ sp.hso3m, sp.h2o2_aq })
        .SetProducts({ sp.so2oohm, sp.h2o })
        .SetSolvent(sp.h2o)
        .SetForwardRateConstant(process::constant::EquilibriumConstant({
            .A_ = c_H2O_M * (7.45e7 / 13.0), .C_ = 4430.0 }))
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 1725.0 }))
        .Build();

    auto rxn1b = process::DissolvedReactionBuilder()
        .SetPhase(sp.aqueous_phase)
        .SetReactants({ sp.so2oohm, sp.hp })
        .SetProducts({ sp.so4mm })
        .SetSolvent(sp.h2o)
        .SetRateConstant([](const Conditions& c) -> double {
            return c_H2O_M * 2.4e6 *
                   std::exp(-4430.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
        })
        .Build();

    auto rxn2 = process::DissolvedReactionBuilder()
        .SetPhase(sp.aqueous_phase)
        .SetReactants({ sp.hso3m, sp.o3_aq })
        .SetProducts({ sp.so4mm, sp.hp })
        .SetSolvent(sp.h2o)
        .SetRateConstant([](const Conditions& c) -> double {
            return c_H2O_M * 3.75e5 *
                   std::exp(-5530.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
        })
        .Build();

    auto rxn3 = process::DissolvedReactionBuilder()
        .SetPhase(sp.aqueous_phase)
        .SetReactants({ sp.so3mm, sp.o3_aq })
        .SetProducts({ sp.so4mm })
        .SetSolvent(sp.h2o)
        .SetRateConstant([](const Conditions& c) -> double {
            return c_H2O_M * 1.59e9 *
                   std::exp(-5280.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
        })
        .Build();

    auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
    model.AddProcesses(rxn1a, rxn1b, rxn2, rxn3);
    model.AddConstraints(hl_so2, hl_h2o2, hl_o3,
                         eq_kw, eq_ka1, eq_ka2,
                         mass_S, mass_H2O2, mass_O3, charge);
    return std::make_pair(std::move(model), std::move(cloud));
  }
}  // namespace

// ════════════════════════════════════════════════════════════════════════
// PHASE A: Individual process/constraint solvent sweeps
// ════════════════════════════════════════════════════════════════════════

// ────────────────────────────────────────────────────────────────────────
// A1: DissolvedReaction — HSO3⁻ + O3(aq) → SO4²⁻ + H⁺
//     Rate ∝ 1/[S]^(n_r - 1), n_r = 2 → divides by [S]
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, A1_DissolvedReaction_SolventSweep)
{
  PrintSweepHeader("A1: DissolvedReaction");

  for (double h2o_level : SOLVENT_SWEEP)
  {
    auto o3_aq  = Species{ "O3_aq" };
    auto hso3m  = Species{ "HSO3m" };
    auto so4mm  = Species{ "SO4mm" };
    auto hp     = Species{ "Hp" };
    auto h2o    = Species{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase aqueous_phase{ "AQUEOUS", { h2o, o3_aq, hso3m, so4mm, hp } };
    auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

    auto rxn = process::DissolvedReactionBuilder()
        .SetPhase(aqueous_phase)
        .SetReactants({ hso3m, o3_aq })
        .SetProducts({ so4mm, hp })
        .SetSolvent(h2o)
        .SetRateConstant([](const Conditions& c) -> double {
            return c_H2O_M * 3.75e5 *
                   std::exp(-5530.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
        })
        .Build();

    // Mass-S: HSO3m + SO4mm = const
    auto mass_S = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, hso3m)
        .AddTerm(aqueous_phase, hso3m, 1.0)
        .AddTerm(aqueous_phase, so4mm, 1.0)
        .SetConstant(1.0e-3)
        .Build();

    // Charge: Hp - HSO3m - 2*SO4mm = 0
    auto charge = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, hp)
        .AddTerm(aqueous_phase, hp, 1.0)
        .AddTerm(aqueous_phase, hso3m, -1.0)
        .AddTerm(aqueous_phase, so4mm, -2.0)
        .SetConstant(0.0)
        .Build();

    auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
    model.AddProcesses({ rxn });
    model.AddConstraints(mass_S, charge);

    Phase dummy_gas{ "GAS", {} };
    auto system = System(dummy_gas, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_h2o   = FindIdx(state, "CLOUD.AQUEOUS.H2O");
    auto i_o3    = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
    auto i_hso3m = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
    auto i_so4mm = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
    auto i_hp    = FindIdx(state, "CLOUD.AQUEOUS.Hp");

    state.variables_[0][i_h2o]   = h2o_level;
    state.variables_[0][i_o3]    = 5.0e-7;
    state.variables_[0][i_hso3m] = 1.0e-3;
    state.variables_[0][i_so4mm] = 0.0;
    state.variables_[0][i_hp]    = 1.0e-3;
    cloud.SetDefaultParameters(state);

    auto stats = IntegrateWithStats(solver, state, 60.0, 0.001);
    PrintSweepRow(h2o_level, stats);

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf in final state at H2O=" << h2o_level;
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// A2: DissolvedReversibleReaction — HSO3⁻ + H2O2(aq) ⇌ SO2OOH⁻ + H2O
//     Forward ∝ 1/[S]^(n_r-1), Reverse ∝ 1/[S]^(n_p-1); both n=2
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, A2_DissolvedReversibleReaction_SolventSweep)
{
  PrintSweepHeader("A2: DissolvedReversibleReaction");

  for (double h2o_level : SOLVENT_SWEEP)
  {
    auto hso3m   = Species{ "HSO3m" };
    auto h2o2_aq = Species{ "H2O2_aq" };
    auto so2oohm = Species{ "SO2OOHm" };
    auto h2o     = Species{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase aqueous_phase{ "AQUEOUS", { h2o, hso3m, h2o2_aq, so2oohm } };
    auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

    auto rxn = process::DissolvedReversibleReactionBuilder()
        .SetPhase(aqueous_phase)
        .SetReactants({ hso3m, h2o2_aq })
        .SetProducts({ so2oohm, h2o })
        .SetSolvent(h2o)
        .SetForwardRateConstant(process::constant::EquilibriumConstant({
            .A_ = c_H2O_M * (7.45e7 / 13.0), .C_ = 4430.0 }))
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 1725.0 }))
        .Build();

    // Mass-S: HSO3m + SO2OOHm = const
    auto mass = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, hso3m)
        .AddTerm(aqueous_phase, hso3m, 1.0)
        .AddTerm(aqueous_phase, so2oohm, 1.0)
        .SetConstant(1.0e-3)
        .Build();

    auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
    model.AddProcesses({ rxn });
    model.AddConstraints(mass);

    Phase dummy_gas{ "GAS", {} };
    auto system = System(dummy_gas, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_h2o   = FindIdx(state, "CLOUD.AQUEOUS.H2O");
    auto i_hso3m = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
    auto i_h2o2  = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
    auto i_inter = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");

    state.variables_[0][i_h2o]   = h2o_level;
    state.variables_[0][i_hso3m] = 1.0e-3;
    state.variables_[0][i_h2o2]  = 1.0e-3;
    state.variables_[0][i_inter] = 0.0;
    cloud.SetDefaultParameters(state);

    auto stats = IntegrateWithStats(solver, state, 60.0, 0.001);
    PrintSweepRow(h2o_level, stats);

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf in final state at H2O=" << h2o_level;
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// A3: DissolvedEquilibriumConstraint — SO2(aq) ⇌ HSO3⁻ + H⁺ (Ka1)
//     Constraint residual divides by [S]^(n_p - 1) = [S]
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, A3_DissolvedEquilibriumConstraint_SolventSweep)
{
  PrintSweepHeader("A3: DissolvedEquilibriumConstraint (Ka1)");

  for (double h2o_level : SOLVENT_SWEEP)
  {
    auto so2_aq = Species{ "SO2_aq" };
    auto hso3m  = Species{ "HSO3m" };
    auto hp     = Species{ "Hp" };
    auto h2o    = Species{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase aqueous_phase{ "AQUEOUS", { h2o, so2_aq, hso3m, hp } };
    auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

    auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
        .SetPhase(aqueous_phase).SetReactants({ so2_aq }).SetProducts({ hso3m, hp })
        .SetAlgebraicSpecies(hso3m).SetSolvent(h2o)
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
        .Build();

    // Charge: Hp - HSO3m = 0
    auto charge = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, hp)
        .AddTerm(aqueous_phase, hp, 1.0)
        .AddTerm(aqueous_phase, hso3m, -1.0)
        .SetConstant(0.0)
        .Build();

    auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
    model.AddConstraints(eq_ka1, charge);

    Phase dummy_gas{ "GAS", {} };
    auto system = System(dummy_gas, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_h2o   = FindIdx(state, "CLOUD.AQUEOUS.H2O");
    auto i_so2   = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
    auto i_hso3m = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
    auto i_hp    = FindIdx(state, "CLOUD.AQUEOUS.Hp");

    state.variables_[0][i_h2o]   = h2o_level;
    state.variables_[0][i_so2]   = 1.0e-4;
    state.variables_[0][i_hso3m] = 0.0;
    state.variables_[0][i_hp]    = 1.0e-4;
    cloud.SetDefaultParameters(state);

    // Pure constraint system — just integrate a short time
    auto stats = IntegrateWithStats(solver, state, 10.0, 0.01);
    PrintSweepRow(h2o_level, stats);

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf at H2O=" << h2o_level;
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// A4: HenryLawEquilibriumConstraint — SO2(g) ⇌ SO2(aq)
//     Multiplies by [S] (safe) — should work at all solvent levels.
//     When [S]=0: constraint becomes -[SO2_aq] = 0 (trivially satisfied)
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, A4_HenryLawConstraint_SolventSweep)
{
  PrintSweepHeader("A4: HenryLawEquilibriumConstraint");

  for (double h2o_level : SOLVENT_SWEEP)
  {
    auto so2_g  = Species{ "SO2" };
    auto so2_aq = Species{ "SO2_aq" };
    auto h2o    = Species{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase gas_phase{ "GAS", { so2_g } };
    Phase aqueous_phase{ "AQUEOUS", { so2_aq, h2o } };
    auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

    auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
        .SetGasSpecies(so2_g).SetCondensedSpecies(so2_aq).SetSolvent(h2o)
        .SetCondensedPhase(aqueous_phase)
        .SetHenryLawConstant(process::constant::HenrysLawConstant({
            .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
        .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

    auto mass = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(gas_phase, so2_g)
        .AddTerm(gas_phase, so2_g, 1.0)
        .AddTerm(aqueous_phase, so2_aq, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
    model.AddConstraints(hl_so2, mass);

    auto system = System(gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_g  = FindIdx(state, "SO2");
    auto i_aq = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
    auto i_w  = FindIdx(state, "CLOUD.AQUEOUS.H2O");

    state.variables_[0][i_g]  = 3.01e-8;
    state.variables_[0][i_aq] = 0.0;
    state.variables_[0][i_w]  = h2o_level;
    cloud.SetDefaultParameters(state);

    auto stats = IntegrateWithStats(solver, state, 10.0, 0.01);
    PrintSweepRow(h2o_level, stats);

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf at H2O=" << h2o_level;
      double g_f  = state.variables_[0][i_g];
      double aq_f = state.variables_[0][i_aq];
      EXPECT_NEAR(g_f + aq_f, 3.01e-8, 1e-12) << "Mass violated at H2O=" << h2o_level;

      // When H2O ≈ 0, all SO2 should stay in gas phase
      if (h2o_level < 1e-9)
      {
        EXPECT_NEAR(g_f, 3.01e-8, 1e-10 * 3.01e-8)
            << "SO2 should be entirely in gas phase at H2O=" << h2o_level;
      }
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// A5: Combined HLC + Dissociation (HLC safe, dissociation divides by [S])
//     SO2(g) ⇌ SO2(aq) → HSO3⁻ + H⁺  with mass balance
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, A5_HLC_Plus_Dissociation_SolventSweep)
{
  PrintSweepHeader("A5: HLC + Dissociation");

  for (double h2o_level : SOLVENT_SWEEP)
  {
    auto so2_g  = Species{ "SO2" };
    auto so2_aq = Species{ "SO2_aq" };
    auto hso3m  = Species{ "HSO3m" };
    auto hp     = Species{ "Hp" };
    auto h2o    = Species{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase gas_phase{ "GAS", { so2_g } };
    Phase aqueous_phase{ "AQUEOUS", { h2o, so2_aq, hso3m, hp } };
    auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

    auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
        .SetGasSpecies(so2_g).SetCondensedSpecies(so2_aq).SetSolvent(h2o)
        .SetCondensedPhase(aqueous_phase)
        .SetHenryLawConstant(process::constant::HenrysLawConstant({
            .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
        .SetMwSolvent(0.018).SetRhoSolvent(1000.0).Build();

    auto eq_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
        .SetPhase(aqueous_phase).SetReactants({ so2_aq }).SetProducts({ hso3m, hp })
        .SetAlgebraicSpecies(hso3m).SetSolvent(h2o)
        .SetEquilibriumConstant(process::constant::EquilibriumConstant({
            .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
        .Build();

    auto mass_S = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(gas_phase, so2_g)
        .AddTerm(gas_phase, so2_g, 1.0)
        .AddTerm(aqueous_phase, so2_aq, 1.0)
        .AddTerm(aqueous_phase, hso3m, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto charge = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, hp)
        .AddTerm(aqueous_phase, hp, 1.0)
        .AddTerm(aqueous_phase, hso3m, -1.0)
        .SetConstant(0.0)
        .Build();

    auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
    model.AddConstraints(hl_so2, eq_ka1, mass_S, charge);

    auto system = System(gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_g     = FindIdx(state, "SO2");
    auto i_aq    = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
    auto i_hso3m = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
    auto i_hp    = FindIdx(state, "CLOUD.AQUEOUS.Hp");
    auto i_w     = FindIdx(state, "CLOUD.AQUEOUS.H2O");

    state.variables_[0][i_g]     = 3.01e-8;
    state.variables_[0][i_aq]    = 0.0;
    state.variables_[0][i_hso3m] = 0.0;
    state.variables_[0][i_hp]    = 1.0e-4;  // arbitrary guess
    state.variables_[0][i_w]     = h2o_level;
    cloud.SetDefaultParameters(state);

    auto stats = IntegrateWithStats(solver, state, 10.0, 0.01);
    PrintSweepRow(h2o_level, stats);

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf at H2O=" << h2o_level;
      double total = state.variables_[0][i_g] + state.variables_[0][i_aq] +
                     state.variables_[0][i_hso3m];
      EXPECT_NEAR(total, 3.01e-8, 1e-12) << "S budget violated at H2O=" << h2o_level;
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// A6: Full equilibrium system (no kinetics) — Step 3 configuration
//     All 3 HLCs + Kw + Ka1 + Ka2 + mass/charge constraints
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, A6_FullEquilibrium_SolventSweep)
{
  PrintSweepHeader("A6: Full Equilibrium (no kinetics)");

  for (double h2o_level : SOLVENT_SWEEP)
  {
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
    auto h2o     = Species{ "H2O",
        {{ "molecular weight [kg mol-1]", 0.018 },
         { "density [kg m-3]", 1000.0 }} };

    Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g } };
    Phase aqueous_phase{ "AQUEOUS", {
        h2o, so2_aq, h2o2_aq, o3_aq, hp, ohm, hso3m, so3mm } };
    auto cloud = representation::UniformSection{ "CLOUD", { aqueous_phase } };

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
            .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M), .C_ = 6710.0 }))
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

    auto mass_S = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(gas_phase, so2_g)
        .AddTerm(gas_phase, so2_g, 1.0)
        .AddTerm(aqueous_phase, so2_aq, 1.0)
        .AddTerm(aqueous_phase, hso3m, 1.0)
        .AddTerm(aqueous_phase, so3mm, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto mass_H2O2 = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(gas_phase, h2o2_g)
        .AddTerm(gas_phase, h2o2_g, 1.0)
        .AddTerm(aqueous_phase, h2o2_aq, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto mass_O3 = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(gas_phase, o3_g)
        .AddTerm(gas_phase, o3_g, 1.0)
        .AddTerm(aqueous_phase, o3_aq, 1.0)
        .DiagnoseConstantFromState()
        .Build();

    auto charge = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, hp)
        .AddTerm(aqueous_phase, hp, 1.0)
        .AddTerm(aqueous_phase, ohm, -1.0)
        .AddTerm(aqueous_phase, hso3m, -1.0)
        .AddTerm(aqueous_phase, so3mm, -2.0)
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
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_so2_g   = FindIdx(state, "SO2");
    auto i_h2o2_g  = FindIdx(state, "H2O2");
    auto i_o3_g    = FindIdx(state, "O3");
    auto i_w       = FindIdx(state, "CLOUD.AQUEOUS.H2O");
    auto i_so2_aq  = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
    auto i_h2o2_aq = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
    auto i_o3_aq   = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
    auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");
    auto i_ohm     = FindIdx(state, "CLOUD.AQUEOUS.OHm");
    auto i_hso3m   = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
    auto i_so3mm   = FindIdx(state, "CLOUD.AQUEOUS.SO3mm");

    state.variables_[0][i_so2_g]   = 3.01e-8;
    state.variables_[0][i_h2o2_g]  = 3.01e-8;
    state.variables_[0][i_o3_g]    = 1.50e-6;
    state.variables_[0][i_w]       = h2o_level;
    state.variables_[0][i_so2_aq]  = 0.0;
    state.variables_[0][i_h2o2_aq] = 0.0;
    state.variables_[0][i_o3_aq]   = 0.0;
    state.variables_[0][i_hp]      = 1.0;  // arbitrary guess
    state.variables_[0][i_ohm]     = 0.0;
    state.variables_[0][i_hso3m]   = 0.0;
    state.variables_[0][i_so3mm]   = 0.0;
    cloud.SetDefaultParameters(state);

    auto stats = IntegrateWithStats(solver, state, 10.0, 0.01);
    PrintSweepRow(h2o_level, stats);

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf at H2O=" << h2o_level;
    }
  }
  std::cout << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// PHASE B: Full CAM cloud mechanism with vanishing cloud water
// ════════════════════════════════════════════════════════════════════════

// ────────────────────────────────────────────────────────────────────────
// B1: Static sweep — run full mechanism at each solvent level for 1800s
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, B1_FullMechanism_StaticSweep)
{
  PrintSweepHeader("B1: Full CAM Cloud Mechanism (1800s)");

  // Extended header for this test
  std::cout << std::setw(14) << "H2O [mol/m³]"
            << std::setw(12) << "converged"
            << std::setw(10) << "steps"
            << std::setw(10) << "rejected"
            << std::setw(10) << "f_calls"
            << std::setw(12) << "SO4 [mol/m³]"
            << std::setw(14) << "failure_mode"
            << std::endl;
  std::cout << std::string(82, '-') << std::endl;

  const std::vector<double> sweep = { 0.017, 1e-3, 1e-5, 1e-7, 1e-9, 1e-12 };

  for (double h2o_level : sweep)
  {
    CamCloudSpecies sp;
    auto [model, cloud] = BuildCamCloudModel(sp);

    auto system = System(sp.gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();
    state.conditions_[0].temperature_ = 280.0;
    state.conditions_[0].pressure_ = 70000.0;
    state.conditions_[0].CalculateIdealAirDensity();

    auto i_so2_g   = FindIdx(state, "SO2");
    auto i_h2o2_g  = FindIdx(state, "H2O2");
    auto i_o3_g    = FindIdx(state, "O3");
    auto i_h2o     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
    auto i_so2_aq  = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
    auto i_h2o2_aq = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
    auto i_o3_aq   = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
    auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");
    auto i_ohm     = FindIdx(state, "CLOUD.AQUEOUS.OHm");
    auto i_hso3m   = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
    auto i_so3mm   = FindIdx(state, "CLOUD.AQUEOUS.SO3mm");
    auto i_so4mm   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
    auto i_so2oohm = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");

    // Naive ICs
    state.variables_[0][i_so2_g]   = 3.01e-8;
    state.variables_[0][i_h2o2_g]  = 3.01e-8;
    state.variables_[0][i_o3_g]    = 1.50e-6;
    state.variables_[0][i_h2o]     = h2o_level;
    state.variables_[0][i_so2_aq]  = 0.0;
    state.variables_[0][i_h2o2_aq] = 0.0;
    state.variables_[0][i_o3_aq]   = 0.0;
    state.variables_[0][i_hp]      = 1.0;
    state.variables_[0][i_ohm]     = 0.0;
    state.variables_[0][i_hso3m]   = 0.0;
    state.variables_[0][i_so3mm]   = 0.0;
    state.variables_[0][i_so4mm]   = 1.0;
    state.variables_[0][i_so2oohm] = 0.0;
    cloud.SetDefaultParameters(state);

    auto stats = IntegrateWithStats(solver, state, 1800.0, 0.001);

    double so4_f = stats.converged ? state.variables_[0][i_so4mm] : -1.0;
    std::cout << std::scientific << std::setprecision(1) << std::setw(14) << h2o_level
              << std::setw(12) << (stats.converged ? "YES" : "NO")
              << std::dec << std::setw(10) << stats.total_steps
              << std::setw(10) << stats.total_rejected
              << std::setw(10) << stats.total_function_calls
              << std::scientific << std::setprecision(3) << std::setw(12) << so4_f
              << std::setw(14)
              << (stats.converged ? "-" : SolverStateToString(stats.last_failure_state))
              << std::endl;

    if (stats.converged)
    {
      EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf at H2O=" << h2o_level;

      // S conservation
      double total_S = state.variables_[0][i_so2_g] + state.variables_[0][i_so2_aq] +
                        state.variables_[0][i_hso3m] + state.variables_[0][i_so3mm] +
                        state.variables_[0][i_so4mm] + state.variables_[0][i_so2oohm];
      EXPECT_NEAR(total_S, 1.0 + 3.01e-8, 1e-4 * total_S)
          << "S budget violated at H2O=" << h2o_level;
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// B2: Dynamic evaporation — cloud water decreases during integration
//     Piecewise-constant H2O, halving every 200s from 0.017 to ~1e-12
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, B2_DynamicEvaporation)
{
  std::cout << "\n┌─── B2: Dynamic Evaporation ───" << std::endl;

  CamCloudSpecies sp;
  auto [model, cloud] = BuildCamCloudModel(sp);

  auto system = System(sp.gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = 280.0;
  state.conditions_[0].pressure_ = 70000.0;
  state.conditions_[0].CalculateIdealAirDensity();

  auto i_so2_g   = FindIdx(state, "SO2");
  auto i_h2o2_g  = FindIdx(state, "H2O2");
  auto i_o3_g    = FindIdx(state, "O3");
  auto i_h2o     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so4mm   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
  auto i_so2oohm = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");
  auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");

  // Start with normal cloud water
  state.variables_[0][i_so2_g]   = 3.01e-8;
  state.variables_[0][i_h2o2_g]  = 3.01e-8;
  state.variables_[0][i_o3_g]    = 1.50e-6;
  state.variables_[0][i_h2o]     = C_H2O_NORMAL;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.SO2_aq")]  = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq")] = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.O3_aq")]   = 0.0;
  state.variables_[0][i_hp]      = 1.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.OHm")]     = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.HSO3m")]   = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.SO3mm")]   = 0.0;
  state.variables_[0][i_so4mm]   = 1.0;
  state.variables_[0][i_so2oohm] = 0.0;
  cloud.SetDefaultParameters(state);

  std::cout << std::setw(10) << "t_start"
            << std::setw(10) << "t_end"
            << std::setw(14) << "H2O [mol/m³]"
            << std::setw(12) << "converged"
            << std::setw(10) << "steps"
            << std::setw(10) << "rejected"
            << std::setw(12) << "SO4"
            << std::endl;
  std::cout << std::string(78, '-') << std::endl;

  double h2o = C_H2O_NORMAL;
  double t_start = 0.0;
  double segment = 200.0;
  bool failed = false;

  for (int seg = 0; seg < 10 && !failed; ++seg)
  {
    double t_end = t_start + segment;
    state.variables_[0][i_h2o] = h2o;

    auto stats = IntegrateWithStats(solver, state, segment, 0.001);

    std::cout << std::fixed << std::setprecision(0)
              << std::setw(10) << t_start
              << std::setw(10) << t_end
              << std::scientific << std::setprecision(1) << std::setw(14) << h2o
              << std::setw(12) << (stats.converged ? "YES" : "NO")
              << std::dec << std::setw(10) << stats.total_steps
              << std::setw(10) << stats.total_rejected
              << std::scientific << std::setprecision(3) << std::setw(12) << state.variables_[0][i_so4mm]
              << std::endl;

    if (!stats.converged)
    {
      std::cout << "  *** Failed at H2O=" << h2o << ": "
                << SolverStateToString(stats.last_failure_state) << std::endl;
      failed = true;
    }

    t_start = t_end;
    h2o *= 0.1;  // decrease by 10× each segment
  }
  std::cout << std::endl;

  // Not asserting full convergence — this is characterization
  // But we do check whether NaN/Inf crept into the state
  if (!failed)
  {
    EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf in final state after evaporation";
  }
}

// ────────────────────────────────────────────────────────────────────────
// B3: Cloud formation — start with no cloud, then add water
//     Tests constraint initialization when cloud appears mid-simulation
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, B3_CloudFormation)
{
  std::cout << "\n┌─── B3: Cloud Formation (clear→cloudy) ───" << std::endl;

  CamCloudSpecies sp;
  auto [model, cloud] = BuildCamCloudModel(sp);

  auto system = System(sp.gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();
  state.conditions_[0].temperature_ = 280.0;
  state.conditions_[0].pressure_ = 70000.0;
  state.conditions_[0].CalculateIdealAirDensity();

  auto i_so2_g   = FindIdx(state, "SO2");
  auto i_h2o2_g  = FindIdx(state, "H2O2");
  auto i_o3_g    = FindIdx(state, "O3");
  auto i_h2o     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so4mm   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
  auto i_so2oohm = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");
  auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");

  // Phase 1: Very small cloud water (near-clear)
  state.variables_[0][i_so2_g]   = 3.01e-8;
  state.variables_[0][i_h2o2_g]  = 3.01e-8;
  state.variables_[0][i_o3_g]    = 1.50e-6;
  state.variables_[0][i_h2o]     = 1e-15;  // nearly zero
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.SO2_aq")]  = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq")] = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.O3_aq")]   = 0.0;
  state.variables_[0][i_hp]      = 1.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.OHm")]     = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.HSO3m")]   = 0.0;
  state.variables_[0][FindIdx(state, "CLOUD.AQUEOUS.SO3mm")]   = 0.0;
  state.variables_[0][i_so4mm]   = 1.0;
  state.variables_[0][i_so2oohm] = 0.0;
  cloud.SetDefaultParameters(state);

  std::cout << "Phase 1: near-clear (H2O=1e-15) for 100s..." << std::endl;
  auto stats1 = IntegrateWithStats(solver, state, 100.0, 0.001);
  std::cout << "  converged=" << (stats1.converged ? "YES" : "NO")
            << " steps=" << stats1.total_steps
            << " rejected=" << stats1.total_rejected << std::endl;

  // Phase 2: Cloud formation — set H2O to normal
  state.variables_[0][i_h2o] = C_H2O_NORMAL;

  std::cout << "Phase 2: cloud formed (H2O=0.017) for 1800s..." << std::endl;
  auto stats2 = IntegrateWithStats(solver, state, 1800.0, 0.001);
  std::cout << "  converged=" << (stats2.converged ? "YES" : "NO")
            << " steps=" << stats2.total_steps
            << " rejected=" << stats2.total_rejected
            << " SO4=" << state.variables_[0][i_so4mm] << std::endl;

  if (stats2.converged)
  {
    EXPECT_FALSE(HasNanOrInf(state)) << "NaN/Inf after cloud formation";
    EXPECT_GT(state.variables_[0][i_so4mm], 1.0) << "SO4 should increase after cloud forms";
  }
  std::cout << std::endl;
}

// ════════════════════════════════════════════════════════════════════════
// PHASE C: Multi-cell vertical column tests
// ════════════════════════════════════════════════════════════════════════

// ────────────────────────────────────────────────────────────────────────
// C1: Mixed cloudy/clear column — 5 cells
//     [0, 1e-4, 0.017, 1e-4, 0]  with vertical T/P profile
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, C1_MixedColumn)
{
  std::cout << "\n┌─── C1: Mixed Cloudy/Clear Column (5 cells) ───" << std::endl;

  constexpr std::size_t NUM_CELLS = 5;
  const double h2o_profile[NUM_CELLS]  = { 0.0, 1e-4, 0.017, 1e-4, 0.0 };
  const double T_profile[NUM_CELLS]    = { 280.0, 275.0, 270.0, 260.0, 245.0 };
  const double P_profile[NUM_CELLS]    = { 70000.0, 65000.0, 55000.0, 45000.0, 35000.0 };

  CamCloudSpecies sp;
  auto [model, cloud] = BuildCamCloudModel(sp);

  auto system = System(sp.gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState(NUM_CELLS);

  // Get indices (same for all cells)
  auto i_so2_g   = FindIdx(state, "SO2");
  auto i_h2o2_g  = FindIdx(state, "H2O2");
  auto i_o3_g    = FindIdx(state, "O3");
  auto i_h2o     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so2_aq  = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
  auto i_h2o2_aq = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
  auto i_o3_aq   = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
  auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");
  auto i_ohm     = FindIdx(state, "CLOUD.AQUEOUS.OHm");
  auto i_hso3m   = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
  auto i_so3mm   = FindIdx(state, "CLOUD.AQUEOUS.SO3mm");
  auto i_so4mm   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
  auto i_so2oohm = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");

  for (std::size_t cell = 0; cell < NUM_CELLS; ++cell)
  {
    state.conditions_[cell].temperature_ = T_profile[cell];
    state.conditions_[cell].pressure_ = P_profile[cell];
    state.conditions_[cell].CalculateIdealAirDensity();

    state.variables_[cell][i_so2_g]   = 3.01e-8;
    state.variables_[cell][i_h2o2_g]  = 3.01e-8;
    state.variables_[cell][i_o3_g]    = 1.50e-6;
    state.variables_[cell][i_h2o]     = h2o_profile[cell];
    state.variables_[cell][i_so2_aq]  = 0.0;
    state.variables_[cell][i_h2o2_aq] = 0.0;
    state.variables_[cell][i_o3_aq]   = 0.0;
    state.variables_[cell][i_hp]      = 1.0;
    state.variables_[cell][i_ohm]     = 0.0;
    state.variables_[cell][i_hso3m]   = 0.0;
    state.variables_[cell][i_so3mm]   = 0.0;
    state.variables_[cell][i_so4mm]   = 1.0;
    state.variables_[cell][i_so2oohm] = 0.0;
  }
  cloud.SetDefaultParameters(state);

  auto stats = IntegrateWithStats(solver, state, 1800.0, 0.001);

  std::cout << "Overall: converged=" << (stats.converged ? "YES" : "NO")
            << " steps=" << stats.total_steps
            << " rejected=" << stats.total_rejected;
  if (!stats.converged)
    std::cout << " failure=" << SolverStateToString(stats.last_failure_state);
  std::cout << std::endl;

  // Per-cell results
  std::cout << std::setw(6) << "cell"
            << std::setw(14) << "H2O"
            << std::setw(8) << "T"
            << std::setw(12) << "SO4"
            << std::setw(12) << "SO2_g"
            << std::setw(10) << "NaN/Inf"
            << std::endl;
  std::cout << std::string(62, '-') << std::endl;

  for (std::size_t cell = 0; cell < NUM_CELLS; ++cell)
  {
    bool nan_inf = HasNanOrInf(state, cell);
    double so4 = state.variables_[cell][i_so4mm];
    double so2_g = state.variables_[cell][i_so2_g];

    std::cout << std::setw(6) << cell
              << std::scientific << std::setprecision(1)
              << std::setw(14) << h2o_profile[cell]
              << std::fixed << std::setprecision(0)
              << std::setw(8) << T_profile[cell]
              << std::scientific << std::setprecision(3)
              << std::setw(12) << so4
              << std::setw(12) << so2_g
              << std::setw(10) << (nan_inf ? "YES" : "no")
              << std::endl;

    EXPECT_FALSE(nan_inf) << "Cell " << cell << " has NaN/Inf values";

    // Clear cells (H2O=0) should have all aqueous species at zero
    // and gas species unchanged
    if (h2o_profile[cell] == 0.0 && stats.converged)
    {
      EXPECT_NEAR(state.variables_[cell][i_so2_aq], 0.0, 1e-15)
          << "Clear cell " << cell << " has non-zero SO2_aq";
      EXPECT_NEAR(state.variables_[cell][i_h2o2_aq], 0.0, 1e-15)
          << "Clear cell " << cell << " has non-zero H2O2_aq";
    }
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// C2: Uniform cloud column — 5 identical cells as baseline
//     All H2O = 0.017; verifies multi-cell ≈ 5× single-cell
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, C2_UniformColumn)
{
  std::cout << "\n┌─── C2: Uniform Cloud Column (5 cells) ───" << std::endl;

  constexpr std::size_t NUM_CELLS = 5;

  CamCloudSpecies sp;
  auto [model, cloud] = BuildCamCloudModel(sp);

  auto system = System(sp.gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState(NUM_CELLS);

  auto i_so2_g   = FindIdx(state, "SO2");
  auto i_h2o2_g  = FindIdx(state, "H2O2");
  auto i_o3_g    = FindIdx(state, "O3");
  auto i_h2o     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so4mm   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
  auto i_so2oohm = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");
  auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");

  for (std::size_t cell = 0; cell < NUM_CELLS; ++cell)
  {
    state.conditions_[cell].temperature_ = 280.0;
    state.conditions_[cell].pressure_ = 70000.0;
    state.conditions_[cell].CalculateIdealAirDensity();

    state.variables_[cell][i_so2_g]   = 3.01e-8;
    state.variables_[cell][i_h2o2_g]  = 3.01e-8;
    state.variables_[cell][i_o3_g]    = 1.50e-6;
    state.variables_[cell][i_h2o]     = C_H2O_NORMAL;
    state.variables_[cell][FindIdx(state, "CLOUD.AQUEOUS.SO2_aq")]  = 0.0;
    state.variables_[cell][FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq")] = 0.0;
    state.variables_[cell][FindIdx(state, "CLOUD.AQUEOUS.O3_aq")]   = 0.0;
    state.variables_[cell][i_hp]      = 1.0;
    state.variables_[cell][FindIdx(state, "CLOUD.AQUEOUS.OHm")]     = 0.0;
    state.variables_[cell][FindIdx(state, "CLOUD.AQUEOUS.HSO3m")]   = 0.0;
    state.variables_[cell][FindIdx(state, "CLOUD.AQUEOUS.SO3mm")]   = 0.0;
    state.variables_[cell][i_so4mm]   = 1.0;
    state.variables_[cell][i_so2oohm] = 0.0;
  }
  cloud.SetDefaultParameters(state);

  auto stats = IntegrateWithStats(solver, state, 1800.0, 0.001);

  std::cout << "Converged=" << (stats.converged ? "YES" : "NO")
            << " steps=" << stats.total_steps
            << " rejected=" << stats.total_rejected
            << " f_calls=" << stats.total_function_calls
            << std::endl;

  EXPECT_TRUE(stats.converged) << "Uniform column should converge";

  if (stats.converged)
  {
    // All cells should produce identical SO4
    double so4_ref = state.variables_[0][i_so4mm];
    for (std::size_t cell = 0; cell < NUM_CELLS; ++cell)
    {
      EXPECT_FALSE(HasNanOrInf(state, cell)) << "Cell " << cell << " has NaN/Inf";
      double so4 = state.variables_[cell][i_so4mm];
      EXPECT_NEAR(so4, so4_ref, 1e-10 * so4_ref)
          << "Cell " << cell << " SO4=" << so4 << " differs from cell 0 SO4=" << so4_ref;
    }
    std::cout << "SO4 (all cells) = " << so4_ref << " mol/m³" << std::endl;
  }
  std::cout << std::endl;
}

// ────────────────────────────────────────────────────────────────────────
// C3: Extreme contrast — [0, 0, 0.017, 0, 0]
//     Single cloudy cell surrounded by clear; most representative of 3D
//     Tests whether one cell's Jacobian singularity causes solver-wide
//     step-size collapse (all cells share the same time step)
// ────────────────────────────────────────────────────────────────────────
TEST(SolventRobustness, C3_ExtremeContrast)
{
  std::cout << "\n┌─── C3: Extreme Contrast Column [0,0,0.017,0,0] ───" << std::endl;

  constexpr std::size_t NUM_CELLS = 5;
  const double h2o_profile[NUM_CELLS] = { 0.0, 0.0, 0.017, 0.0, 0.0 };

  CamCloudSpecies sp;
  auto [model, cloud] = BuildCamCloudModel(sp);

  auto system = System(sp.gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState(NUM_CELLS);

  auto i_so2_g   = FindIdx(state, "SO2");
  auto i_h2o2_g  = FindIdx(state, "H2O2");
  auto i_o3_g    = FindIdx(state, "O3");
  auto i_h2o     = FindIdx(state, "CLOUD.AQUEOUS.H2O");
  auto i_so2_aq  = FindIdx(state, "CLOUD.AQUEOUS.SO2_aq");
  auto i_h2o2_aq = FindIdx(state, "CLOUD.AQUEOUS.H2O2_aq");
  auto i_o3_aq   = FindIdx(state, "CLOUD.AQUEOUS.O3_aq");
  auto i_hp      = FindIdx(state, "CLOUD.AQUEOUS.Hp");
  auto i_ohm     = FindIdx(state, "CLOUD.AQUEOUS.OHm");
  auto i_hso3m   = FindIdx(state, "CLOUD.AQUEOUS.HSO3m");
  auto i_so3mm   = FindIdx(state, "CLOUD.AQUEOUS.SO3mm");
  auto i_so4mm   = FindIdx(state, "CLOUD.AQUEOUS.SO4mm");
  auto i_so2oohm = FindIdx(state, "CLOUD.AQUEOUS.SO2OOHm");

  for (std::size_t cell = 0; cell < NUM_CELLS; ++cell)
  {
    state.conditions_[cell].temperature_ = 280.0;
    state.conditions_[cell].pressure_ = 70000.0;
    state.conditions_[cell].CalculateIdealAirDensity();

    state.variables_[cell][i_so2_g]   = 3.01e-8;
    state.variables_[cell][i_h2o2_g]  = 3.01e-8;
    state.variables_[cell][i_o3_g]    = 1.50e-6;
    state.variables_[cell][i_h2o]     = h2o_profile[cell];
    state.variables_[cell][i_so2_aq]  = 0.0;
    state.variables_[cell][i_h2o2_aq] = 0.0;
    state.variables_[cell][i_o3_aq]   = 0.0;
    state.variables_[cell][i_hp]      = 1.0;
    state.variables_[cell][i_ohm]     = 0.0;
    state.variables_[cell][i_hso3m]   = 0.0;
    state.variables_[cell][i_so3mm]   = 0.0;
    state.variables_[cell][i_so4mm]   = 1.0;
    state.variables_[cell][i_so2oohm] = 0.0;
  }
  cloud.SetDefaultParameters(state);

  auto stats = IntegrateWithStats(solver, state, 1800.0, 0.001);

  std::cout << "Overall: converged=" << (stats.converged ? "YES" : "NO")
            << " steps=" << stats.total_steps
            << " rejected=" << stats.total_rejected;
  if (!stats.converged)
    std::cout << " failure=" << SolverStateToString(stats.last_failure_state);
  std::cout << std::endl;

  // Per-cell results
  std::cout << std::setw(6) << "cell"
            << std::setw(14) << "H2O"
            << std::setw(12) << "SO4"
            << std::setw(12) << "SO2_g"
            << std::setw(12) << "SO2_aq"
            << std::setw(10) << "NaN/Inf"
            << std::endl;
  std::cout << std::string(66, '-') << std::endl;

  for (std::size_t cell = 0; cell < NUM_CELLS; ++cell)
  {
    bool nan_inf = HasNanOrInf(state, cell);
    std::cout << std::setw(6) << cell
              << std::scientific << std::setprecision(1)
              << std::setw(14) << h2o_profile[cell]
              << std::setprecision(3)
              << std::setw(12) << state.variables_[cell][i_so4mm]
              << std::setw(12) << state.variables_[cell][i_so2_g]
              << std::setw(12) << state.variables_[cell][i_so2_aq]
              << std::setw(10) << (nan_inf ? "YES" : "no")
              << std::endl;

    EXPECT_FALSE(nan_inf) << "Cell " << cell << " has NaN/Inf values";

    // Clear cells should have zero aqueous concentrations
    if (h2o_profile[cell] == 0.0 && stats.converged)
    {
      EXPECT_NEAR(state.variables_[cell][i_so2_aq], 0.0, 1e-15)
          << "Clear cell " << cell << " has non-zero SO2_aq";
      EXPECT_NEAR(state.variables_[cell][i_hso3m], 0.0, 1e-15)
          << "Clear cell " << cell << " has non-zero HSO3m";
    }

    // Cloudy cell should show sulfate production
    if (h2o_profile[cell] > 0.01 && stats.converged)
    {
      EXPECT_GT(state.variables_[cell][i_so4mm], 1.0)
          << "Cloudy cell " << cell << " should have SO4 increase";
    }
  }
  std::cout << std::endl;
}
