// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Integration tests for DissolvedEquilibriumConstraint and LinearConstraint working
// together with the MICM DAE solver.

#include <miam/miam.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;
using namespace miam;

// ============================================================================
// Test 1: DissolvedEquilibriumConstraint + LinearConstraint with kinetic driver
//
// System:
//   A → B  (irreversible dissolved reaction, rate k)
//   B <-> C (equilibrium constraint, K_eq), C is algebraic
//   [A] + [B] + [C] = total (mass conservation linear constraint), B is algebraic
//
// A is differential: dA/dt = -k*[A]
// B is algebraic: [A]+[B]+[C] = total
// C is algebraic: K_eq*[B] - [C] = 0
//
// Analytical solution:
//   A(t)  = A0 * exp(-k*t)
//   B(t)  = (total - A(t)) / (1 + K_eq)
//   C(t)  = K_eq * (total - A(t)) / (1 + K_eq)
// ============================================================================
TEST(EquilibriumConstraintsIntegration, DissolvedEquilibriumWithKineticDriver)
{
  // Species
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto S = Species{ "S" };  // solvent

  // Phase and representation
  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };
  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  // Parameters
  double k = 0.1;     // s⁻¹ (rate for A → B)
  double K_eq = 2.0;  // dimensionless (equilibrium constant for B <-> C)
  double A0 = 1.0;    // mol/m³ initial [A]
  double total = A0;  // total mass (B0 = C0 = 0)

  // Kinetic process: A → B (dissolved reaction)
  auto rate = [k](const Conditions& conditions) { return k; };
  auto reaction = process::DissolvedReaction{
    rate,
    { A },   // reactants
    { B },   // products
    S,       // solvent
    aqueous_phase
  };

  // Equilibrium constraint: B <-> C, K_eq, algebraic species = C
  // Residual: G = K_eq * [B] - [C] = 0
  auto equil = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(S)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant(
          process::constant::EquilibriumConstantParameters{ .A_ = K_eq }))
      .Build();

  // Mass conservation constraint: [A] + [B] + [C] = total, algebraic species = B
  // Residual: G = [A] + [B] + [C] - total = 0
  auto mass_cons = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, B)
      .AddTerm(aqueous_phase, A, 1.0)
      .AddTerm(aqueous_phase, B, 1.0)
      .AddTerm(aqueous_phase, C, 1.0)
      .SetConstant(total)
      .Build();

  // Build model
  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });
  model.AddConstraints(equil, mass_cons);

  // Build DAE solver
  Phase gas_phase{ "GAS", {} };
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  // Find variable indices
  auto find_idx = [&](const std::string& name) {
    auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
    EXPECT_NE(it, state.variable_names_.end()) << "Species " << name << " not found";
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  std::size_t i_A = find_idx("DROPLET.AQUEOUS.A");
  std::size_t i_B = find_idx("DROPLET.AQUEOUS.B");
  std::size_t i_C = find_idx("DROPLET.AQUEOUS.C");
  std::size_t i_S = find_idx("DROPLET.AQUEOUS.S");

  // Set initial conditions (consistent with constraints: B0=C0=0, total=A0)
  state.variables_[0][i_A] = A0;
  state.variables_[0][i_B] = 0.0;
  state.variables_[0][i_C] = 0.0;
  state.variables_[0][i_S] = 1.0e-4;  // solvent (constant)

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  // Time integration (1/k = 10s, integrate to 100s for 10*tau)
  std::vector<double> test_times = { 1.0, 5.0, 10.0, 100.0 };
  double time = 0.0;
  const double tolerance = 5.0e-3;  // Allow some tolerance for DAE solver

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(0.1, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s with dt = " << dt;
      time += dt;
    }

    double A_num = state.variables_[0][i_A];
    double B_num = state.variables_[0][i_B];
    double C_num = state.variables_[0][i_C];

    // Analytical solution
    double A_an = A0 * std::exp(-k * time);
    double B_an = (total - A_an) / (1.0 + K_eq);
    double C_an = K_eq * (total - A_an) / (1.0 + K_eq);

    EXPECT_NEAR(A_num, A_an, tolerance)
      << "Species A mismatch at t = " << time << " s";
    EXPECT_NEAR(B_num, B_an, tolerance)
      << "Species B mismatch at t = " << time << " s";
    EXPECT_NEAR(C_num, C_an, tolerance)
      << "Species C mismatch at t = " << time << " s";

    // Verify mass conservation
    double mass = A_num + B_num + C_num;
    EXPECT_NEAR(mass, total, 1.0e-4)
      << "Mass conservation violated at t = " << time << " s";

    // Verify equilibrium relation
    if (B_num > 1.0e-10)
    {
      EXPECT_NEAR(C_num / B_num, K_eq, 1.0e-3)
        << "Equilibrium ratio violated at t = " << time << " s";
    }
  }

  // Steady state check (at t=100, exp(-0.1*100) = exp(-10) ~ 4.5e-5)
  EXPECT_NEAR(state.variables_[0][i_A], 0.0, 1.0e-3)
    << "A should have decayed to zero by t = " << time;
  EXPECT_NEAR(state.variables_[0][i_B], total / (1.0 + K_eq), tolerance)
    << "B should be total/(1+K_eq) at steady state";
  EXPECT_NEAR(state.variables_[0][i_C], total * K_eq / (1.0 + K_eq), tolerance)
    << "C should be total*K_eq/(1+K_eq) at steady state";
}

// ============================================================================
// Test 2: Per-instance equilibrium constraint with two representations
//
// Same chemistry as Test 1 (A → B kinetics, B <-> C equilibrium), but with
// two differently-sized droplet representations (SMALL, LARGE) with different
// initial concentrations. Only the equilibrium constraint is used (no mass
// conservation constraint, since a single constant cannot serve instances with
// different totals). The kinetics naturally conserve A+B, and C is purely
// algebraic via the equilibrium relation C = K_eq * B.
//
// Analytical solution per instance (A0 = initial [A] for that instance):
//   A(t) = A0 * exp(-k*t)
//   B(t) = A0 * (1 - exp(-k*t))        [from kinetic conservation: A+B = A0]
//   C(t) = K_eq * B(t)                  [algebraic equilibrium]
//
// Steady state: A → 0, B → A0, C → K_eq * A0
// ============================================================================
TEST(EquilibriumConstraintsIntegration, PerInstanceEquilibrium)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto S = Species{ "S" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };

  auto small_drop = representation::UniformSection{
    "SMALL",
    { aqueous_phase }
  };
  auto large_drop = representation::UniformSection{
    "LARGE",
    { aqueous_phase }
  };

  double k = 0.05;
  double K_eq = 3.0;
  double A0_small = 1.0;
  double A0_large = 2.0;

  auto rate = [k](const Conditions& conditions) { return k; };
  auto reaction = process::DissolvedReaction{
    rate,
    { A },
    { B },
    S,
    aqueous_phase
  };

  // Equilibrium constraint: C = K_eq * B (C is algebraic)
  auto equil = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(S)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant(
          process::constant::EquilibriumConstantParameters{ .A_ = K_eq }))
      .Build();

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { small_drop, large_drop }
  };
  model.AddProcesses({ reaction });
  model.AddConstraints(equil);

  Phase gas_phase{ "GAS", {} };
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  auto find_idx = [&](const std::string& name) {
    auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
    EXPECT_NE(it, state.variable_names_.end()) << "Species " << name << " not found";
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  std::size_t i_A_small = find_idx("SMALL.AQUEOUS.A");
  std::size_t i_B_small = find_idx("SMALL.AQUEOUS.B");
  std::size_t i_C_small = find_idx("SMALL.AQUEOUS.C");
  std::size_t i_S_small = find_idx("SMALL.AQUEOUS.S");
  std::size_t i_A_large = find_idx("LARGE.AQUEOUS.A");
  std::size_t i_B_large = find_idx("LARGE.AQUEOUS.B");
  std::size_t i_C_large = find_idx("LARGE.AQUEOUS.C");
  std::size_t i_S_large = find_idx("LARGE.AQUEOUS.S");

  // Different starting concentrations for each instance
  state.variables_[0][i_A_small] = A0_small;
  state.variables_[0][i_B_small] = 0.0;
  state.variables_[0][i_C_small] = 0.0;
  state.variables_[0][i_S_small] = 1.0e-4;

  state.variables_[0][i_A_large] = A0_large;
  state.variables_[0][i_B_large] = 0.0;
  state.variables_[0][i_C_large] = 0.0;
  state.variables_[0][i_S_large] = 1.0e-4;

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  small_drop.SetDefaultParameters(state);
  large_drop.SetDefaultParameters(state);

  // Integrate to steady state (1/k = 20s; 200s = 10*tau for thorough decay)
  double time = 0.0;
  while (time < 200.0)
  {
    double dt = 0.1;
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged)
      << "Solver failed at t = " << time;
    time += dt;
  }

  double tolerance = 5.0e-3;

  // Verify SMALL mode: A → 0, B → A0_small, C → K_eq * A0_small
  EXPECT_NEAR(state.variables_[0][i_A_small], 0.0, tolerance);
  EXPECT_NEAR(state.variables_[0][i_B_small], A0_small, tolerance);
  EXPECT_NEAR(state.variables_[0][i_C_small], K_eq * A0_small, tolerance);

  // Verify LARGE mode: A → 0, B → A0_large, C → K_eq * A0_large
  EXPECT_NEAR(state.variables_[0][i_A_large], 0.0, tolerance);
  EXPECT_NEAR(state.variables_[0][i_B_large], A0_large, tolerance);
  EXPECT_NEAR(state.variables_[0][i_C_large], K_eq * A0_large, tolerance);

  // Verify kinetic mass conservation: A + B = A0 for each instance
  double AB_small = state.variables_[0][i_A_small] + state.variables_[0][i_B_small];
  double AB_large = state.variables_[0][i_A_large] + state.variables_[0][i_B_large];
  EXPECT_NEAR(AB_small, A0_small, 1.0e-4);
  EXPECT_NEAR(AB_large, A0_large, 1.0e-4);

  // Verify the two instances are independent (different steady states)
  EXPECT_GT(std::abs(state.variables_[0][i_B_large] - state.variables_[0][i_B_small]), 0.5)
    << "LARGE and SMALL should have distinct B values";
  EXPECT_GT(std::abs(state.variables_[0][i_C_large] - state.variables_[0][i_C_small]), 1.0)
    << "LARGE and SMALL should have distinct C values";

  // Verify equilibrium relation holds for both instances
  EXPECT_NEAR(state.variables_[0][i_C_small] / state.variables_[0][i_B_small], K_eq, 1.0e-3);
  EXPECT_NEAR(state.variables_[0][i_C_large] / state.variables_[0][i_B_large], K_eq, 1.0e-3);
}
