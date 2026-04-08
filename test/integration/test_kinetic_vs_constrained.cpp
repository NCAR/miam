// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Integration tests comparing fully kinetic ODE systems with algebraically
// constrained DAE systems. At steady state, both should produce equivalent results.

#include <miam/miam.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <tuple>

using namespace micm;
using namespace miam;

namespace
{
  constexpr double R_gas = miam::util::R_gas;

  // Solve a kinetic ODE system to steady state and return final concentrations
  template<typename FindIdx>
  std::tuple<double, double, double> SolveKineticSystem(
      double k, double k_f, double k_r, double A0, FindIdx find_idx)
  {
    auto A = Species{ "A" };
    auto B = Species{ "B" };
    auto C = Species{ "C" };
    auto S = Species{ "S" };

    auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };
    auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

    auto rate = [k](const Conditions& conditions) { return k; };
    auto reaction = process::DissolvedReaction{
      rate, { A }, { B }, S, aqueous_phase
    };

    auto forward_rate = [k_f](const Conditions& conditions) { return k_f; };
    auto reverse_rate = [k_r](const Conditions& conditions) { return k_r; };
    auto reversible = process::DissolvedReversibleReaction{
      forward_rate, reverse_rate,
      { B }, { C }, S, aqueous_phase
    };

    auto model = Model{
      .name_ = "AEROSOL",
      .representations_ = { droplet }
    };
    model.AddProcesses({ reaction });
    model.AddProcesses({ reversible });

    Phase gas_phase{ "GAS", {} };
    auto system = System(gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModelProcesses(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();

    std::size_t i_A = find_idx(state, "DROPLET.AQUEOUS.A");
    std::size_t i_B = find_idx(state, "DROPLET.AQUEOUS.B");
    std::size_t i_C = find_idx(state, "DROPLET.AQUEOUS.C");
    std::size_t i_S = find_idx(state, "DROPLET.AQUEOUS.S");

    state.variables_[0][i_A] = A0;
    state.variables_[0][i_B] = 0.0;
    state.variables_[0][i_C] = 0.0;
    state.variables_[0][i_S] = 1.0e-4;

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);

    // Integrate to steady state
    for (int step = 0; step < 5000; ++step)
    {
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(0.01, state);
      if (result.state_ != SolverState::Converged)
        break;
    }

    return { state.variables_[0][i_A],
             state.variables_[0][i_B],
             state.variables_[0][i_C] };
  }
}  // namespace

// ============================================================================
// Test 1: Kinetic dissolved reversible reaction vs equilibrium constraint
//
// Kinetic system (ODE):
//   A → B with rate k (irreversible)
//   B <-> C with rates k_f, k_r (reversible, K_eq = k_f/k_r)
//
// Constrained system (DAE):
//   A → B with rate k (irreversible)
//   B <-> C with K_eq (equilibrium constraint, C algebraic)
//   [A]+[B]+[C] = total (mass conservation, B algebraic)
//
// At steady state, A → 0, and both systems should agree:
//   [B] = total / (1 + K_eq)
//   [C] = total * K_eq / (1 + K_eq)
// ============================================================================
TEST(KineticVsConstrained, DissolvedReversibleVsEquilibriumConstraint)
{
  double k = 0.1;       // rate for A → B
  double K_eq = 5.0;    // equilibrium constant for B <-> C
  double k_f = 100.0;   // forward rate for B → C (fast)
  double k_r = k_f / K_eq;  // reverse rate for C → B
  double A0 = 1.0;

  auto find_idx = [](const auto& state, const std::string& name) -> std::size_t {
    auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  // --- Kinetic system ---
  auto [A_kin, B_kin, C_kin] = [&]() {
    auto A = Species{ "A" };
    auto B = Species{ "B" };
    auto C = Species{ "C" };
    auto S = Species{ "S" };

    auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };
    auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

    auto rate = [k](const Conditions& conditions) { return k; };
    auto reaction = process::DissolvedReaction{
      rate, { A }, { B }, S, aqueous_phase
    };

    auto forward_rate = [k_f](const Conditions& conditions) { return k_f; };
    auto reverse_rate = [k_r](const Conditions& conditions) { return k_r; };
    auto reversible = process::DissolvedReversibleReaction{
      forward_rate, reverse_rate,
      { B }, { C }, S, aqueous_phase
    };

    auto model = Model{
      .name_ = "AEROSOL",
      .representations_ = { droplet }
    };
    model.AddProcesses({ reaction });
    model.AddProcesses({ reversible });

    Phase gas_phase{ "GAS", {} };
    auto system = System(gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModelProcesses(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();

    std::size_t i_A = find_idx(state, "DROPLET.AQUEOUS.A");
    std::size_t i_B = find_idx(state, "DROPLET.AQUEOUS.B");
    std::size_t i_C = find_idx(state, "DROPLET.AQUEOUS.C");
    std::size_t i_S = find_idx(state, "DROPLET.AQUEOUS.S");

    state.variables_[0][i_A] = A0;
    state.variables_[0][i_B] = 0.0;
    state.variables_[0][i_C] = 0.0;
    state.variables_[0][i_S] = 1.0e-4;

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);

    // Integrate well past all timescales
    for (int step = 0; step < 5000; ++step)
    {
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(0.01, state);
      EXPECT_EQ(result.state_, SolverState::Converged) << "Kinetic solver failed at step " << step;
    }

    return std::make_tuple(
        state.variables_[0][i_A],
        state.variables_[0][i_B],
        state.variables_[0][i_C]);
  }();

  // --- Constrained system (DAE) ---
  auto [A_dae, B_dae, C_dae] = [&]() {
    auto A = Species{ "A" };
    auto B = Species{ "B" };
    auto C = Species{ "C" };
    auto S = Species{ "S" };

    auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };
    auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

    auto rate = [k](const Conditions& conditions) { return k; };
    auto reaction = process::DissolvedReaction{
      rate, { A }, { B }, S, aqueous_phase
    };

    auto equil = constraint::DissolvedEquilibriumConstraintBuilder()
        .SetPhase(aqueous_phase)
        .SetReactants({ B })
        .SetProducts({ C })
        .SetAlgebraicSpecies(C)
        .SetSolvent(S)
        .SetEquilibriumConstant(process::constant::EquilibriumConstant(
            process::constant::EquilibriumConstantParameters{ .A_ = K_eq }))
        .Build();

    auto mass_cons = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(aqueous_phase, B)
        .AddTerm(aqueous_phase, A, 1.0)
        .AddTerm(aqueous_phase, B, 1.0)
        .AddTerm(aqueous_phase, C, 1.0)
        .SetConstant(A0)
        .Build();

    auto model = Model{
      .name_ = "AEROSOL",
      .representations_ = { droplet }
    };
    model.AddProcesses({ reaction });
    model.AddConstraints(equil, mass_cons);

    Phase gas_phase{ "GAS", {} };
    auto system = System(gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();

    std::size_t i_A = find_idx(state, "DROPLET.AQUEOUS.A");
    std::size_t i_B = find_idx(state, "DROPLET.AQUEOUS.B");
    std::size_t i_C = find_idx(state, "DROPLET.AQUEOUS.C");
    std::size_t i_S = find_idx(state, "DROPLET.AQUEOUS.S");

    state.variables_[0][i_A] = A0;
    state.variables_[0][i_B] = 0.0;
    state.variables_[0][i_C] = 0.0;
    state.variables_[0][i_S] = 1.0e-4;

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);

    // DAE system reaches equilibrium instantly for the constrained variables;
    // just need enough time for A to decay
    for (int step = 0; step < 500; ++step)
    {
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(0.1, state);
      EXPECT_EQ(result.state_, SolverState::Converged) << "DAE solver failed at step " << step;
    }

    return std::make_tuple(
        state.variables_[0][i_A],
        state.variables_[0][i_B],
        state.variables_[0][i_C]);
  }();

  // Both systems should agree at steady state
  double total = A0;
  double B_expected = total / (1.0 + K_eq);
  double C_expected = total * K_eq / (1.0 + K_eq);
  const double tolerance = 1.0e-2;

  // Verify kinetic system reached expected steady state
  EXPECT_NEAR(A_kin, 0.0, tolerance) << "Kinetic: A should be ~0";
  EXPECT_NEAR(B_kin, B_expected, tolerance) << "Kinetic: B mismatch";
  EXPECT_NEAR(C_kin, C_expected, tolerance) << "Kinetic: C mismatch";

  // Verify DAE system reached expected steady state
  EXPECT_NEAR(A_dae, 0.0, tolerance) << "DAE: A should be ~0";
  EXPECT_NEAR(B_dae, B_expected, tolerance) << "DAE: B mismatch";
  EXPECT_NEAR(C_dae, C_expected, tolerance) << "DAE: C mismatch";

  // Direct comparison: kinetic vs constrained
  EXPECT_NEAR(A_kin, A_dae, tolerance) << "Kinetic vs DAE: A mismatch";
  EXPECT_NEAR(B_kin, B_dae, tolerance) << "Kinetic vs DAE: B mismatch";
  EXPECT_NEAR(C_kin, C_dae, tolerance) << "Kinetic vs DAE: C mismatch";

  // Verify mass conservation in both systems
  double mass_kin = A_kin + B_kin + C_kin;
  double mass_dae = A_dae + B_dae + C_dae;
  EXPECT_NEAR(mass_kin, total, 1.0e-3) << "Kinetic mass conservation violation";
  EXPECT_NEAR(mass_dae, total, 1.0e-3) << "DAE mass conservation violation";
}

// ============================================================================
// Test 2: Kinetic HL phase transfer vs HL equilibrium constraint
//
// Kinetic system (ODE):
//   A_g <-> A_aq via HenryLawPhaseTransfer (kinetic rate-based mass transfer)
//
// Constrained system (DAE):
//   A_g <-> A_aq via HenryLawEquilibriumConstraint (algebraic equilibrium)
//   [A_g] + [A_aq] = total (mass conservation, A_g algebraic)
//
// At steady state, both should give:
//   [A_aq] / ([A_g] * R * T) = HLC * f_v  (adjusted Henry's Law)
// ============================================================================
TEST(KineticVsConstrained, HenryLawPhaseTransferVsEquilibriumConstraint)
{
  double T = 298.15;
  double HLC = 3.4e-2;         // mol m⁻³ Pa⁻¹
  double Mw_gas = 0.044;       // kg mol⁻¹
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double D_g = 1.5e-5;         // m² s⁻¹
  double accommodation = 0.05;
  double H2O_conc = rho_solvent / Mw_solvent;  // 55555.6 mol/m³
  double f_v = H2O_conc * Mw_solvent / rho_solvent;
  double alpha = HLC * R_gas * T * f_v;

  auto A_g = Species{ "A_g",
      { { "molecular weight [kg mol-1]", Mw_gas } } };
  auto A_aq = Species{ "A_aq",
      { { "molecular weight [kg mol-1]", Mw_gas },
        { "density [kg m-3]", 1800.0 } } };
  auto H2O = Species{ "H2O",
      { { "molecular weight [kg mol-1]", Mw_solvent },
        { "density [kg m-3]", rho_solvent } } };

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  double gas0 = 1.0e-3;  // mol/m³ initial gas concentration
  double aq0 = 0.0;
  double total = gas0 + aq0;

  // --- Kinetic system (HL phase transfer) ---
  double kinetic_A_g, kinetic_A_aq;
  {
    auto droplet = representation::SingleMomentMode{
      "DROPLET", { aqueous_phase }, 5.0e-6, 1.2
    };

    auto transfer = process::HenryLawPhaseTransferBuilder()
        .SetCondensedPhase(aqueous_phase)
        .SetGasSpecies(A_g)
        .SetCondensedSpecies(A_aq)
        .SetSolvent(H2O)
        .SetHenrysLawConstant(process::constant::HenrysLawConstant(
            process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC }))
        .SetDiffusionCoefficient(D_g)
        .SetAccommodationCoefficient(accommodation)
        .Build();

    auto model = Model{
      .name_ = "AEROSOL",
      .representations_ = { droplet }
    };
    model.AddProcesses({ transfer });

    auto system = System(gas_phase, model);
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                      .SetSystem(system)
                      .AddExternalModelProcesses(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();

    auto find_idx = [&](const std::string& name) {
      auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
      return static_cast<std::size_t>(it - state.variable_names_.begin());
    };

    state.variables_[0][find_idx("A_g")] = gas0;
    state.variables_[0][find_idx("DROPLET.AQUEOUS.A_aq")] = aq0;
    state.variables_[0][find_idx("DROPLET.AQUEOUS.H2O")] = H2O_conc;

    state.conditions_[0].temperature_ = T;
    state.conditions_[0].pressure_ = 101325.0;
    state.conditions_[0].CalculateIdealAirDensity();

    droplet.SetDefaultParameters(state);

    // Integrate to equilibrium (long enough for mass transfer)
    for (int step = 0; step < 10000; ++step)
    {
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(0.01, state);
      if (result.state_ != SolverState::Converged)
        break;
    }

    kinetic_A_g = state.variables_[0][find_idx("A_g")];
    kinetic_A_aq = state.variables_[0][find_idx("DROPLET.AQUEOUS.A_aq")];
  }

  // --- Constrained system (HL equilibrium + mass conservation) ---
  double dae_A_g, dae_A_aq;
  {
    // Use UniformSection for the constrained system — no kinetic mass transfer
    auto droplet = representation::UniformSection{
      "DROPLET", { aqueous_phase }
    };

    auto hl_constraint = constraint::HenryLawEquilibriumConstraintBuilder()
        .SetGasSpecies(A_g)
        .SetCondensedSpecies(A_aq)
        .SetSolvent(H2O)
        .SetCondensedPhase(aqueous_phase)
        .SetHenryLawConstant(process::constant::HenrysLawConstant(
            process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC }))
        .SetMwSolvent(Mw_solvent)
        .SetRhoSolvent(rho_solvent)
        .Build();

    auto mass_cons = constraint::LinearConstraintBuilder()
        .SetAlgebraicSpecies(gas_phase, A_g)
        .AddTerm(gas_phase, A_g, 1.0)
        .AddTerm(aqueous_phase, A_aq, 1.0)
        .SetConstant(total)
        .Build();

    auto model = Model{
      .name_ = "AEROSOL",
      .representations_ = { droplet }
    };
    model.AddConstraints(hl_constraint, mass_cons);

    auto system = System(gas_phase, model);

    // Need at least one differential variable — A_g and A_aq are both algebraic.
    // Create a dummy inert species to have a differential equation.
    // Actually, A_g is algebraic (mass conservation) and A_aq is algebraic (HL constraint).
    // With ONLY algebraic variables and no differential equations, the DAE solver might
    // not converge well. Let's add a dummy inert gas species that doesn't participate.
    auto Inert = Species{ "Inert" };
    Phase gas_phase_with_inert{ "GAS", { { A_g }, { Inert } } };
    auto system_with_inert = System(gas_phase_with_inert, model);

    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                      RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                      .SetSystem(system_with_inert)
                      .AddExternalModel(model)
                      .SetIgnoreUnusedSpecies(true)
                      .Build();

    State state = solver.GetState();

    auto find_idx = [&](const std::string& name) {
      auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
      return static_cast<std::size_t>(it - state.variable_names_.begin());
    };

    state.variables_[0][find_idx("A_g")] = gas0;
    state.variables_[0][find_idx("DROPLET.AQUEOUS.A_aq")] = aq0;
    state.variables_[0][find_idx("DROPLET.AQUEOUS.H2O")] = H2O_conc;
    state.variables_[0][find_idx("Inert")] = 1.0;  // dummy species

    state.conditions_[0].temperature_ = T;
    state.conditions_[0].pressure_ = 101325.0;
    state.conditions_[0].CalculateIdealAirDensity();

    droplet.SetDefaultParameters(state);

    // DAE solver should reach equilibrium quickly
    for (int step = 0; step < 100; ++step)
    {
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(1.0, state);
      EXPECT_EQ(result.state_, SolverState::Converged) << "DAE solver failed at step " << step;
    }

    dae_A_g = state.variables_[0][find_idx("A_g")];
    dae_A_aq = state.variables_[0][find_idx("DROPLET.AQUEOUS.A_aq")];
  }

  // Expected equilibrium values (analytical)
  double gas_expected = total / (1.0 + alpha);
  double aq_expected = total * alpha / (1.0 + alpha);

  const double tolerance = 5.0e-3;

  // Verify the constrained (DAE) system matches the analytical solution
  EXPECT_NEAR(dae_A_g, gas_expected, tolerance * total)
    << "DAE: A_g mismatch vs analytical";
  EXPECT_NEAR(dae_A_aq, aq_expected, tolerance * total)
    << "DAE: A_aq mismatch vs analytical";

  // Verify the kinetic system reached the same equilibrium
  EXPECT_NEAR(kinetic_A_g, gas_expected, tolerance * total)
    << "Kinetic: A_g mismatch vs analytical";
  EXPECT_NEAR(kinetic_A_aq, aq_expected, tolerance * total)
    << "Kinetic: A_aq mismatch vs analytical";

  // Direct comparison: kinetic vs constrained
  EXPECT_NEAR(kinetic_A_g, dae_A_g, tolerance * total)
    << "Kinetic vs DAE: A_g mismatch";
  EXPECT_NEAR(kinetic_A_aq, dae_A_aq, tolerance * total)
    << "Kinetic vs DAE: A_aq mismatch";

  // Mass conservation in both
  EXPECT_NEAR(kinetic_A_g + kinetic_A_aq, total, 1.0e-6)
    << "Kinetic: mass conservation";
  EXPECT_NEAR(dae_A_g + dae_A_aq, total, 1.0e-6)
    << "DAE: mass conservation";
}
