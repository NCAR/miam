// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/miam.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;
using namespace miam;

// ============================================================================
// Test 1: Simple first-order decay A -> B
// Analytical: [A](t) = A0 * exp(-k*t)
//             [B](t) = A0 * (1 - exp(-k*t))
// ============================================================================
TEST(DissolvedReactionIntegration, SimpleFirstOrderDecay)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };

  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  double k = 0.1;  // s^-1
  auto rate = [k](const Conditions& conditions) { return k; };

  // A -> B with solvent C
  auto reaction = process::DissolvedReaction{
    rate,
    { A },   // reactants
    { B },   // products
    C,       // solvent
    aqueous_phase
  };

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });

  Phase gas_phase{ "GAS", {} };
  double A0 = 1.0;  // mol/m^3

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

  State state = solver.GetState();

  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
  ASSERT_LT(i_A, state.variable_names_.size());
  ASSERT_LT(i_B, state.variable_names_.size());
  ASSERT_LT(i_C, state.variable_names_.size());

  state.variables_[0][i_A] = A0;
  state.variables_[0][i_B] = 0.0;
  state.variables_[0][i_C] = 1.0e-4;  // solvent (constant)

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  std::vector<double> test_times = { 5.0, 10.0, 20.0, 100.0 };
  double time = 0.0;
  const double tolerance = 2.0e-4;

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(0.01, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s";
      time += dt;
    }

    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];

    // Analytical: [A] = A0 * exp(-k*t), [B] = A0 * (1 - exp(-k*t))
    double A_analytic = A0 * std::exp(-k * time);
    double B_analytic = A0 * (1.0 - std::exp(-k * time));

    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Species A mismatch at t = " << time << " s";
    EXPECT_NEAR(B_numeric, B_analytic, tolerance)
      << "Species B mismatch at t = " << time << " s";

    // Mass conservation: A + B = A0
    EXPECT_NEAR(A_numeric + B_numeric, A0, tolerance)
      << "Mass conservation violated at t = " << time << " s";

    // Solvent should remain constant
    EXPECT_NEAR(C_numeric, 1.0e-4, tolerance)
      << "Solvent changed at t = " << time << " s";
  }

  // At t=100, A should be essentially zero: exp(-0.1*100) = exp(-10) ≈ 4.5e-5
  double A_final = state.variables_[0][i_A];
  EXPECT_NEAR(A_final, 0.0, tolerance)
    << "A should be depleted at t = 100 s";
}

// ============================================================================
// Test 2: Solvent as Reactant  A + C -> B (C is solvent)
// With [C] >> [A], rate ≈ k * [A] (pseudo-first-order)
// ============================================================================
TEST(DissolvedReactionIntegration, SolventAsReactant)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent AND reactant

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };

  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  double k = 1.0e-3;
  auto rate = [k](const Conditions& conditions) { return k; };

  // A + C -> B with solvent C
  // rate = k / [C]^(2-1) * [A] * [C] = k * [A]
  auto reaction = process::DissolvedReaction{
    rate,
    { A, C },  // reactants
    { B },     // products
    C,         // solvent
    aqueous_phase
  };

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });

  Phase gas_phase{ "GAS", {} };

  double A0 = 1.0e-5;
  double C0 = 1.0e-4;

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

  State state = solver.GetState();

  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();

  state.variables_[0][i_A] = A0;
  state.variables_[0][i_B] = 0.0;
  state.variables_[0][i_C] = C0;

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  // rate = k * [A], so [A](t) = A0 * exp(-k * t), pseudo-first-order
  std::vector<double> test_times = { 100.0, 500.0, 1000.0, 5000.0 };
  double time = 0.0;
  const double tolerance = 1.0e-8;

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(1.0, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s";
      time += dt;
    }

    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];

    // Pseudo-first-order: [A](t) = A0 * exp(-k*t)
    double A_analytic = A0 * std::exp(-k * time);
    double B_analytic = A0 - A_analytic;

    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Species A mismatch at t = " << time << " s";
    EXPECT_NEAR(B_numeric, B_analytic, tolerance)
      << "Species B mismatch at t = " << time << " s";

    // Stoichiometry: C consumed = A consumed
    double A_consumed = A0 - A_numeric;
    double C_consumed = C0 - C_numeric;
    EXPECT_NEAR(C_consumed, A_consumed, tolerance)
      << "Stoichiometry violated at t = " << time << " s";

    // Mass balance: A + B = A0
    EXPECT_NEAR(A_numeric + B_numeric, A0, tolerance)
      << "Mass conservation violated at t = " << time << " s";
  }
}

// ============================================================================
// Test 3: Solvent as Product  A -> B + C (C is solvent)
// rate = k / [C]^(1-1) * [A] = k * [A] (first-order)
// ============================================================================
TEST(DissolvedReactionIntegration, SolventAsProduct)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent AND product

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };

  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  double k = 1.0e-3;
  auto rate = [k](const Conditions& conditions) { return k; };

  // A -> B + C with solvent C
  // rate = k * [A] (1 reactant, no solvent normalization)
  auto reaction = process::DissolvedReaction{
    rate,
    { A },      // reactants
    { B, C },   // products (solvent is a product)
    C,          // solvent
    aqueous_phase
  };

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });

  Phase gas_phase{ "GAS", {} };

  double A0 = 0.1;
  double B0 = 0.0;
  double C0 = 1.0e-4;

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

  State state = solver.GetState();

  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();

  state.variables_[0][i_A] = A0;
  state.variables_[0][i_B] = B0;
  state.variables_[0][i_C] = C0;

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  std::vector<double> test_times = { 100.0, 500.0, 1000.0, 5000.0 };
  double time = 0.0;
  const double tolerance = 1.0e-4;

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(1.0, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s";
      time += dt;
    }

    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];

    // [A](t) = A0 * exp(-k*t)
    double A_analytic = A0 * std::exp(-k * time);
    double B_analytic = A0 - A_analytic;
    double C_analytic = C0 + B_analytic;  // C produced at same rate as B

    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Species A mismatch at t = " << time << " s";
    EXPECT_NEAR(B_numeric, B_analytic, tolerance)
      << "Species B mismatch at t = " << time << " s";
    EXPECT_NEAR(C_numeric, C_analytic, tolerance)
      << "Species C mismatch at t = " << time << " s";

    // Stoichiometry: ΔB = ΔC = -ΔA
    double delta_B = B_numeric - B0;
    double delta_C = C_numeric - C0;
    EXPECT_NEAR(delta_B, delta_C, tolerance)
      << "Stoichiometry violated at t = " << time << " s";

    // Mass balance: A + B = A0
    EXPECT_NEAR(A_numeric + B_numeric, A0, tolerance)
      << "Mass conservation violated at t = " << time << " s";
  }
}

// ============================================================================
// Test 4: Multi-Phase Instances (two droplet populations, different rates)
// ============================================================================
TEST(DissolvedReactionIntegration, MultiPhaseInstances)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent

  auto small_aqueous = Phase{ "AQUEOUS_SMALL", { { A }, { B }, { C } } };
  auto large_aqueous = Phase{ "AQUEOUS_LARGE", { { A }, { B }, { C } } };

  auto small_droplet = representation::UniformSection{
    "DROPLET_SMALL",
    { small_aqueous }
  };

  auto large_droplet = representation::UniformSection{
    "DROPLET_LARGE",
    { large_aqueous }
  };

  double k_small = 0.1;
  double k_large = 0.2;

  auto rate_small = [k_small](const Conditions& conditions) { return k_small; };
  auto rate_large = [k_large](const Conditions& conditions) { return k_large; };

  auto reaction_small = process::DissolvedReaction{
    rate_small,
    { A },
    { B },
    C,
    small_aqueous
  };

  auto reaction_large = process::DissolvedReaction{
    rate_large,
    { A },
    { B },
    C,
    large_aqueous
  };

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { small_droplet, large_droplet }
  };
  model.AddProcesses({ reaction_small, reaction_large });

  Phase gas_phase{ "GAS", {} };

  double A0_small = 1.0;
  double A0_large = 0.5;

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

  State state = solver.GetState();

  std::size_t i_A_small = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_SMALL.AQUEOUS_SMALL.A") - state.variable_names_.begin();
  std::size_t i_B_small = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_SMALL.AQUEOUS_SMALL.B") - state.variable_names_.begin();
  std::size_t i_C_small = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_SMALL.AQUEOUS_SMALL.C") - state.variable_names_.begin();

  std::size_t i_A_large = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_LARGE.AQUEOUS_LARGE.A") - state.variable_names_.begin();
  std::size_t i_B_large = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_LARGE.AQUEOUS_LARGE.B") - state.variable_names_.begin();
  std::size_t i_C_large = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_LARGE.AQUEOUS_LARGE.C") - state.variable_names_.begin();

  state.variables_[0][i_A_small] = A0_small;
  state.variables_[0][i_B_small] = 0.0;
  state.variables_[0][i_C_small] = 1.0e-4;

  state.variables_[0][i_A_large] = A0_large;
  state.variables_[0][i_B_large] = 0.0;
  state.variables_[0][i_C_large] = 1.0e-4;

  small_droplet.SetDefaultParameters(state);
  large_droplet.SetDefaultParameters(state);

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  std::vector<double> test_times = { 5.0, 10.0, 20.0, 100.0 };
  double time = 0.0;
  const double tolerance = 2.0e-4;

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(0.01, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s";
      time += dt;
    }

    double A_num_sm = state.variables_[0][i_A_small];
    double B_num_sm = state.variables_[0][i_B_small];
    double A_num_lg = state.variables_[0][i_A_large];
    double B_num_lg = state.variables_[0][i_B_large];

    double A_ana_sm = A0_small * std::exp(-k_small * time);
    double B_ana_sm = A0_small * (1.0 - std::exp(-k_small * time));

    double A_ana_lg = A0_large * std::exp(-k_large * time);
    double B_ana_lg = A0_large * (1.0 - std::exp(-k_large * time));

    EXPECT_NEAR(A_num_sm, A_ana_sm, tolerance)
      << "Small droplet A mismatch at t = " << time << " s";
    EXPECT_NEAR(B_num_sm, B_ana_sm, tolerance)
      << "Small droplet B mismatch at t = " << time << " s";
    EXPECT_NEAR(A_num_sm + B_num_sm, A0_small, tolerance)
      << "Small droplet mass not conserved at t = " << time << " s";

    EXPECT_NEAR(A_num_lg, A_ana_lg, tolerance)
      << "Large droplet A mismatch at t = " << time << " s";
    EXPECT_NEAR(B_num_lg, B_ana_lg, tolerance)
      << "Large droplet B mismatch at t = " << time << " s";
    EXPECT_NEAR(A_num_lg + B_num_lg, A0_large, tolerance)
      << "Large droplet mass not conserved at t = " << time << " s";

    // Solvent should be unchanged
    EXPECT_NEAR(state.variables_[0][i_C_small], 1.0e-4, tolerance);
    EXPECT_NEAR(state.variables_[0][i_C_large], 1.0e-4, tolerance);
  }

  // Large droplets should decay faster (k_large > k_small)
  double A_final_small = state.variables_[0][i_A_small];
  double A_final_large = state.variables_[0][i_A_large];
  EXPECT_LT(A_final_large / A0_large, A_final_small / A0_small)
    << "Large droplets should have decayed more";
}

// ============================================================================
// Test 5: Second-order reaction A + B -> C
// Analytical for equal initial concentrations [A]0 = [B]0 = c0:
//   [A](t) = [B](t) = c0 / (1 + k_eff * c0 * t)
//   where k_eff = k / [S]
//   [C](t) = c0 - [A](t)
// ============================================================================
TEST(DissolvedReactionIntegration, SecondOrderTwoReactants)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto S = Species{ "S" };  // solvent

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };

  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  double k = 1.0;
  auto rate = [k](const Conditions& conditions) { return k; };

  // A + B -> C with solvent S
  // rate = k / [S] * [A] * [B]
  auto reaction = process::DissolvedReaction{
    rate,
    { A, B },  // reactants
    { C },     // products
    S,         // solvent
    aqueous_phase
  };

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });

  Phase gas_phase{ "GAS", {} };

  double c0 = 0.01;    // mol/m^3 (equal initial concentrations)
  double S0 = 50.0;    // mol/m^3 (solvent, large and approximately constant)
  double k_eff = k / S0;  // effective second-order rate constant

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

  State state = solver.GetState();

  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
  std::size_t i_S = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.S") - state.variable_names_.begin();
  ASSERT_LT(i_A, state.variable_names_.size());
  ASSERT_LT(i_B, state.variable_names_.size());
  ASSERT_LT(i_C, state.variable_names_.size());
  ASSERT_LT(i_S, state.variable_names_.size());

  state.variables_[0][i_A] = c0;
  state.variables_[0][i_B] = c0;
  state.variables_[0][i_C] = 0.0;
  state.variables_[0][i_S] = S0;

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  // For [A]0 = [B]0 = c0 and d[A]/dt = -k_eff*[A]*[B] = -k_eff*[A]^2:
  // [A](t) = c0 / (1 + k_eff * c0 * t)
  // half-life: t_half = 1/(k_eff*c0)
  double t_half = 1.0 / (k_eff * c0);
  std::vector<double> test_times = { t_half * 0.1, t_half * 0.5, t_half, t_half * 5.0 };
  double time = 0.0;
  const double tolerance = 1.0e-5;

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(t_half * 0.001, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s";
      time += dt;
    }

    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];
    double S_numeric = state.variables_[0][i_S];

    // Analytical: [A](t) = c0 / (1 + k_eff * c0 * t)
    // Note: solvent [S] changes slightly as it is not a reactant/product,
    // but the analytical solution assumes [S] ≈ S0 (valid for S0 >> c0)
    double A_analytic = c0 / (1.0 + k_eff * c0 * time);
    double C_analytic = c0 - A_analytic;

    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Species A mismatch at t = " << time << " s"
      << " (t/t_half = " << time / t_half << ")";
    // A and B should remain equal (symmetric)
    EXPECT_NEAR(A_numeric, B_numeric, tolerance)
      << "A and B should be equal at t = " << time << " s";
    EXPECT_NEAR(C_numeric, C_analytic, tolerance)
      << "Species C mismatch at t = " << time << " s";

    // Solvent should remain constant (it's not a reactant or product)
    EXPECT_NEAR(S_numeric, S0, tolerance)
      << "Solvent changed at t = " << time << " s";

    // Mass balance: A + C = c0 (and B + C = c0)
    EXPECT_NEAR(A_numeric + C_numeric, c0, tolerance)
      << "Mass conservation violated at t = " << time << " s";
  }
}

// ============================================================================
// Test: max_halflife cap with zero-concentration reactant
// Verifies that the tanh rate cap doesn't produce NaN when a reactant is zero.
// Regression test for 0/0 NaN in ForcingFunctionCapped and JacobianFunctionCapped.
// ============================================================================
TEST(DissolvedReactionIntegration, MaxHalflifeZeroReactant)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // product
  auto S = Species{ "S" };  // solvent

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };

  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  double k = 1.0;  // s^-1
  auto rate = [k](const Conditions& conditions) { return k; };

  double t_half = 10.0;  // seconds

  // A + B -> C with max_halflife cap
  auto reaction = process::DissolvedReaction{
    rate,
    { A, B },  // reactants: two species so soft-min is tested
    { C },     // products
    S,         // solvent
    aqueous_phase,
    1.0e-20,   // solvent_damping_epsilon (default)
    t_half     // max_halflife — triggers the capped code path
  };

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });

  Phase gas_phase{ "GAS", {} };

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

  // --- Sub-test 1: One reactant exactly zero ---
  {
    State state = solver.GetState();

    std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
    std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
    std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
    std::size_t i_S = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.S") - state.variable_names_.begin();
    ASSERT_LT(i_A, state.variable_names_.size());
    ASSERT_LT(i_B, state.variable_names_.size());
    ASSERT_LT(i_C, state.variable_names_.size());
    ASSERT_LT(i_S, state.variable_names_.size());

    state.variables_[0][i_A] = 1.0;    // nonzero
    state.variables_[0][i_B] = 0.0;    // exactly zero — triggers the NaN bug
    state.variables_[0][i_C] = 0.0;
    state.variables_[0][i_S] = 1.0e-4; // solvent

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);
    solver.UpdateStateParameters(state);

    auto result = solver.Solve(1.0, state);
    EXPECT_EQ(result.state_, SolverState::Converged)
      << "Solver should converge when one reactant is exactly zero";

    // With B=0, no reaction should occur — concentrations should be unchanged
    EXPECT_NEAR(state.variables_[0][i_A], 1.0, 1.0e-10)
      << "A should remain unchanged when B=0";
    EXPECT_NEAR(state.variables_[0][i_B], 0.0, 1.0e-10)
      << "B should remain zero";
    EXPECT_NEAR(state.variables_[0][i_C], 0.0, 1.0e-10)
      << "C should remain zero when no reaction occurs";
  }

  // --- Sub-test 2: Both reactants zero ---
  {
    State state = solver.GetState();

    std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
    std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
    std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
    std::size_t i_S = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.S") - state.variable_names_.begin();

    state.variables_[0][i_A] = 0.0;    // both zero
    state.variables_[0][i_B] = 0.0;
    state.variables_[0][i_C] = 0.0;
    state.variables_[0][i_S] = 1.0e-4;

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);
    solver.UpdateStateParameters(state);

    auto result = solver.Solve(1.0, state);
    EXPECT_EQ(result.state_, SolverState::Converged)
      << "Solver should converge when both reactants are zero";
  }

  // --- Sub-test 3: Near-zero reactant (1e-30) — should converge smoothly ---
  {
    State state = solver.GetState();

    std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
    std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
    std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
    std::size_t i_S = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.S") - state.variable_names_.begin();

    state.variables_[0][i_A] = 1.0;
    state.variables_[0][i_B] = 1.0e-30;  // near-zero
    state.variables_[0][i_C] = 0.0;
    state.variables_[0][i_S] = 1.0e-4;

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);
    solver.UpdateStateParameters(state);

    auto result = solver.Solve(1.0, state);
    EXPECT_EQ(result.state_, SolverState::Converged)
      << "Solver should converge with near-zero reactant concentration";

    // Concentrations should remain non-negative
    EXPECT_GE(state.variables_[0][i_A], 0.0) << "A should be non-negative";
    EXPECT_GE(state.variables_[0][i_B], 0.0) << "B should be non-negative";
    EXPECT_GE(state.variables_[0][i_C], 0.0) << "C should be non-negative";
  }
}
