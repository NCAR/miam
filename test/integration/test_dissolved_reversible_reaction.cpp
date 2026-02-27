// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/miam.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;
using namespace miam;

// Analytical solution for A <-> B with constant rate constants
// K_eq = k_f / k_r
// [A]_eq = A0 / (1 + K_eq)
// [B]_eq = A0 * K_eq / (1 + K_eq)
// k_total = k_f + k_r
// [A](t) = [A]_eq + ([A]0 - [A]_eq) * exp(-k_total * t)
// [B](t) = [B]_eq - ([A]0 - [A]_eq) * exp(-k_total * t)
void analytical_solution(double A0, double k_f, double k_r, double t,
                        double& A_t, double& B_t)
{
  double K_eq = k_f / k_r;
  double A_eq = A0 / (1.0 + K_eq);
  double B_eq = A0 * K_eq / (1.0 + K_eq);
  double k_total = k_f + k_r;
  
  double exp_term = std::exp(-k_total * t);
  A_t = A_eq + (A0 - A_eq) * exp_term;
  B_t = B_eq - (A0 - A_eq) * exp_term;
}

TEST(DissolvedReversibleReactionIntegration, SimpleAtoB)
{
  // ========================================================================
  // Test 1: Simple A <-> B (No Solvent Dependence)
  // ========================================================================
  
  // Define species
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent
  
  // Define aqueous phase
  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  
  // Define representation (single uniform section)
  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };
  
  // Rate constants
  double k_forward = 0.1;  // s^-1
  double k_reverse = 0.05; // s^-1
  
  // Create rate constant functions (temperature-independent)
  auto forward_rate = [k_forward](const Conditions& conditions) { return k_forward; };
  auto reverse_rate = [k_reverse](const Conditions& conditions) { return k_reverse; };
  
  // Create reversible reaction: A <-> B with solvent C
  auto reaction = process::DissolvedReversibleReaction{
    forward_rate,
    reverse_rate,
    { A },  // reactants
    { B },  // products
    C,      // solvent
    aqueous_phase
  };
  
  // Create model
  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });
  
  // Create gas phase (empty for this test)
  Phase gas_phase{ "GAS", {} };
  
  // Test parameters
  double A0 = 1.0;  // mol/m^3
  double K_eq = k_forward / k_reverse;
  double A_eq = A0 / (1.0 + K_eq);
  double B_eq = A0 * K_eq / (1.0 + K_eq);
  
  // Build solver
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModelProcesses(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();
  
  // Set up state
  State state = solver.GetState();
  
  // Find variable indices
  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
  
  // Set initial conditions
  state.variables_[0][i_A] = A0;     // [A] = 1.0 mol/m^3
  state.variables_[0][i_B] = 0.0;    // [B] = 0.0 mol/m^3
  state.variables_[0][i_C] = 1.0e-4; // [C] = 1.0e-4 mol/m^3 (solvent, constant)
  
  // Set conditions
  state.conditions_[0].temperature_ = 298.15;  // K
  state.conditions_[0].pressure_ = 101325.0;   // Pa
  
  // Set default parameters for the representation
  droplet.SetDefaultParameters(state);
  
  // Time points to test (in seconds)
  // tau = 6.67s, so 5*tau = 33s gets to 99% of equilibrium, using 100s for good margin
  std::vector<double> test_times = { 5.0, 10.0, 20.0, 100.0 };
  
  double time = 0.0;
  const double tolerance = 1.0e-4;  // mol/m^3 (reasonable for numerical integration)
  
  // Calculate rate constants once before starting
  solver.CalculateRateConstants(state);
  
  for (double target_time : test_times)
  {
    // Integrate to target time (with small tolerance to avoid floating point issues)
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(0.01, target_time - time);  // Smaller time step
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged) 
        << "Solver failed to converge at time " << time << " s with dt = " << dt << " s";
      time += dt;
    }
    
    // Get numerical solution
    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];
    
    // Get analytical solution
    double A_analytic, B_analytic;
    analytical_solution(A0, k_forward, k_reverse, time, A_analytic, B_analytic);
    
    // Check numerical solution matches analytical solution
    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Species A concentration mismatch at t = " << time << " s";
    EXPECT_NEAR(B_numeric, B_analytic, tolerance)
      << "Species B concentration mismatch at t = " << time << " s";
    
    // Check mass conservation
    double mass_sum = A_numeric + B_numeric;
    EXPECT_NEAR(mass_sum, A0, tolerance)
      << "Mass conservation violated at t = " << time << " s";
    
    // Check solvent remains constant
    EXPECT_NEAR(C_numeric, 1.0e-4, tolerance)
      << "Solvent changed at t = " << time << " s";
  }
  
  // Check equilibrium is reached
  double A_final = state.variables_[0][i_A];
  double B_final = state.variables_[0][i_B];
  
  EXPECT_NEAR(A_final, A_eq, tolerance)
    << "Species A did not reach equilibrium at t = " << time << " s";
  EXPECT_NEAR(B_final, B_eq, tolerance)
    << "Species B did not reach equilibrium at t = " << time << " s";
}
TEST(DissolvedReversibleReactionIntegration, SolventAsReactant)
{
  // Test 2: A + C ⇌ B (Solvent is Reactant)
  
  // Define species
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent AND reactant
  
  // Define aqueous phase
  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  
  // Define representation
  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };
  
  // Rate constants
  double k_forward = 1.0e-3;   // s^-1
  double k_reverse = 1.0e-11;  // s^-1
  
  auto forward_rate = [k_forward](const Conditions& conditions) { return k_forward; };
  auto reverse_rate = [k_reverse](const Conditions& conditions) { return k_reverse; };
  
  // Create reversible reaction: A + C <-> B with solvent C
  auto reaction = process::DissolvedReversibleReaction{
    forward_rate,
    reverse_rate,
    { A, C },  // reactants: A + C
    { B },     // products: B
    C,         // solvent
    aqueous_phase
  };
  
  // Create model
  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });
  
  Phase gas_phase{ "GAS", {} };
  
  // Test parameters
  double A0 = 1.0e-5;   // mol/m^3
  double C0 = 1.0e-4;   // mol/m^3 (10x the amount needed)
  double K_eq = k_forward / k_reverse;  // 1e8
  double A_eq = A0 / (1.0 + K_eq);
  double B_eq = A0 * K_eq / (1.0 + K_eq);
  double C_eq = C0 - B_eq;
  
  // Build solver
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModelProcesses(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();
  
  State state = solver.GetState();
  
  // Find variable indices
  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(), "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
  
  // Set initial conditions
  state.variables_[0][i_A] = A0;
  state.variables_[0][i_B] = 0.0;
  state.variables_[0][i_C] = C0;
  
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  
  droplet.SetDefaultParameters(state);
  
  // Time points to test (tau ≈ 1000s)
  std::vector<double> test_times = { 1000.0, 2000.0, 4000.0, 10000.0 };
  
  double time = 0.0;
  const double tolerance = 1.0e-8;  // mol/m^3
  
  solver.CalculateRateConstants(state);
  
  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(1.0, target_time - time);
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged) 
        << "Solver failed to converge at time " << time << " s with dt = " << dt << " s";
      time += dt;
    }
    
    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];
    
    // Analytical solution
    double k_total = k_forward + k_reverse;
    double exp_term = std::exp(-k_total * time);
    double A_analytic = A_eq + (A0 - A_eq) * exp_term;
    double B_analytic = B_eq - (A0 - A_eq) * exp_term;
    
    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Species A concentration mismatch at t = " << time << " s";
    EXPECT_NEAR(B_numeric, B_analytic, tolerance)
      << "Species B concentration mismatch at t = " << time << " s";
    
    // Check stoichiometry: C consumption equals A consumption
    double A_consumed = A0 - A_numeric;
    double C_consumed = C0 - C_numeric;
    EXPECT_NEAR(C_consumed, A_consumed, tolerance)
      << "Stoichiometry violated at t = " << time << " s";
    
    // Check mass balance
    double mass_sum = A_numeric + B_numeric;
    EXPECT_NEAR(mass_sum, A0, tolerance)
      << "Mass conservation violated at t = " << time << " s";
  }
  
  // Check equilibrium
  double A_final = state.variables_[0][i_A];
  double B_final = state.variables_[0][i_B];
  double C_final = state.variables_[0][i_C];
  
  EXPECT_NEAR(A_final, A_eq, tolerance)
    << "Species A did not reach equilibrium";
  EXPECT_NEAR(B_final, B_eq, tolerance)
    << "Species B did not reach equilibrium";
  EXPECT_NEAR(C_final, C_eq, tolerance)
    << "Species C did not reach equilibrium";
  
  // Verify nearly complete conversion (K_eq = 1e8)
  EXPECT_GT(B_final / A0, 0.99)
    << "Expected nearly complete conversion to B";
}

// Test 3: A ⇌ B + C (Solvent is Product)
// Tests case where solvent C is both a product and the solvent species
TEST(DissolvedReversibleReactionIntegration, SolventAsProduct) {
  // Species: A (reactant), B (product), C (product AND solvent)
  // Reaction: A ⇌ B + C where C is both product and solvent
  // K_eq = 1e-8 (heavily favors A, minimal conversion)
  
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent AND product
  
  // Define aqueous phase
  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  
  // Define representation
  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };
  
  // Rate constants
  double k_forward = 1.0e-11;  // s^-1 (very slow forward)
  double k_reverse = 1.0e-3;   // s^-1 (fast reverse)
  
  auto forward_rate = [k_forward](const Conditions& conditions) { return k_forward; };
  auto reverse_rate = [k_reverse](const Conditions& conditions) { return k_reverse; };
  
  // Create reversible reaction: A ⇌ B + C with solvent C
  // Forward: k_f/[C]^0 × [A] = k_f × [A] (no solvent dependence)
  // Reverse: k_r/[C]^1 × [B] × [C] = k_r × [B] (C cancels)
  auto reaction = process::DissolvedReversibleReaction{
    forward_rate,
    reverse_rate,
    { A },     // reactants: A
    { B, C },  // products: B + C
    C,         // solvent
    aqueous_phase
  };
  
  // Create model
  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddProcesses({ reaction });
  
  Phase gas_phase{ "GAS", {} };
  
  // Test parameters
  double A0 = 0.1;      // mol/m^3
  double B0 = 0.0;      // mol/m^3
  double C0 = 1.0e-4;   // mol/m^3 (solvent, very small)
  double K_eq = k_forward / k_reverse;  // 1e-8 (favors A)
  
  // Since C is very small compared to A and K_eq << 1, very little conversion
  // The analytical solution for A ⇌ B + C with minimal conversion:
  // d[A]/dt = -k_f·[A] + k_r·[B] (approximately, since [C] changes slowly)
  double k_total = k_forward + k_reverse;
  double A_eq = k_reverse / (k_forward + k_reverse) * A0;  // ≈ A0
  double B_eq = k_forward / (k_forward + k_reverse) * A0;  // ≈ 0
  
  // Build solver
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModelProcesses(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();
  
  // Initialize state
  State state = solver.GetState();
  std::size_t i_A = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                               "DROPLET.AQUEOUS.A") - state.variable_names_.begin();
  std::size_t i_B = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                               "DROPLET.AQUEOUS.B") - state.variable_names_.begin();
  std::size_t i_C = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                               "DROPLET.AQUEOUS.C") - state.variable_names_.begin();
  
  state.variables_[0][i_A] = A0;
  state.variables_[0][i_B] = B0;
  state.variables_[0][i_C] = C0;
  droplet.SetDefaultParameters(state);
  
  // Tolerance for numerical comparisons
  const double tolerance = 1.0e-4;  // mol/m³
  
  // Test integration at multiple times
  std::vector<double> test_times = {1000.0, 2000.0, 4000.0, 10000.0};
  
  solver.CalculateRateConstants(state);
  double time = 0.0;
  
  for (double target_time : test_times) {
    // Integrate to target time
    while (time < target_time - 1.0e-10) {
      double dt = std::min(1.0, target_time - time);  // 1.0s timestep
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged);
      time += dt;
    }
    
    // Get numerical solution
    double A_numeric = state.variables_[0][i_A];
    double B_numeric = state.variables_[0][i_B];
    double C_numeric = state.variables_[0][i_C];
    
    // Analytical solution (assuming C ≈ constant, which is valid for small conversion)
    double exp_term = std::exp(-k_total * time);
    double A_analytic = A_eq + (A0 - A_eq) * exp_term;
    double B_analytic = B_eq + (B0 - B_eq) * exp_term;
    double C_analytic = C0 + B_analytic;  // C increases with B production
    
    // Check analytical vs numerical
    EXPECT_NEAR(A_numeric, A_analytic, tolerance)
      << "Mismatch at t=" << time << "s";
    EXPECT_NEAR(B_numeric, B_analytic, tolerance)
      << "Mismatch at t=" << time << "s";
    EXPECT_NEAR(C_numeric, C_analytic, tolerance)
      << "Mismatch at t=" << time << "s";
    
    // Check stoichiometry: ΔC = ΔB (both produced together)
    double delta_B = B_numeric - B0;
    double delta_C = C_numeric - C0;
    EXPECT_NEAR(delta_C, delta_B, tolerance)
      << "Stoichiometry violated at t=" << time << "s";
    
    // Check mass conservation: A + B = A0
    EXPECT_NEAR(A_numeric + B_numeric, A0, tolerance)
      << "Mass not conserved at t=" << time << "s";
  }
  
  // Check equilibrium at final time
  double A_final = state.variables_[0][i_A];
  double B_final = state.variables_[0][i_B];
  double C_final = state.variables_[0][i_C];
  
  // At equilibrium: K_eq = [B]·[C] / [A]
  // Since K_eq = 1e-8 << 1, we expect B ≈ 0 and A ≈ A0
  EXPECT_NEAR(B_final, 0.0, tolerance)
    << "System should be near equilibrium with B ≈ 0";
  EXPECT_NEAR(A_final, A0, tolerance)
    << "System should be near equilibrium with A ≈ A0";
}

// Test 4: Multi-Phase Instances
// Tests two independent droplet populations with same reaction but different parameters
TEST(DissolvedReversibleReactionIntegration, MultiPhaseInstances) {
  // Two independent droplet populations:
  // - Small droplets: faster reaction
  // - Large droplets: even faster reaction (same K_eq but faster kinetics)
  
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };  // solvent
  
  // Define aqueous phases for both droplet sizes
  auto small_aqueous = Phase{ "AQUEOUS_SMALL", { { A }, { B }, { C } } };
  auto large_aqueous = Phase{ "AQUEOUS_LARGE", { { A }, { B }, { C } } };
  
  // Define representations
  auto small_droplet = representation::UniformSection{
    "DROPLET_SMALL",
    { small_aqueous }
  };
  
  auto large_droplet = representation::UniformSection{
    "DROPLET_LARGE",
    { large_aqueous }
  };
  
  // Rate constants for small droplets
  double k_forward_small = 0.1;   // s^-1
  double k_reverse_small = 0.05;  // s^-1
  
  auto forward_rate_small = [k_forward_small](const Conditions& conditions) { return k_forward_small; };
  auto reverse_rate_small = [k_reverse_small](const Conditions& conditions) { return k_reverse_small; };
  
  // Create reaction for small droplets: A ⇌ B
  auto reaction_small = process::DissolvedReversibleReaction{
    forward_rate_small,
    reverse_rate_small,
    { A },  // reactants
    { B },  // products
    C,      // solvent
    small_aqueous
  };
  
  // Rate constants for large droplets (2x faster, same K_eq)
  double k_forward_large = 0.2;   // s^-1
  double k_reverse_large = 0.1;   // s^-1
  
  auto forward_rate_large = [k_forward_large](const Conditions& conditions) { return k_forward_large; };
  auto reverse_rate_large = [k_reverse_large](const Conditions& conditions) { return k_reverse_large; };
  
  // Create reaction for large droplets: A ⇌ B
  auto reaction_large = process::DissolvedReversibleReaction{
    forward_rate_large,
    reverse_rate_large,
    { A },  // reactants
    { B },  // products
    C,      // solvent
    large_aqueous
  };
  
  // Create model with both representations
  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { small_droplet, large_droplet }
  };
  model.AddProcesses({ reaction_small, reaction_large });
  
  Phase gas_phase{ "GAS", {} };
  
  // Test parameters for small droplets
  double A0_small = 1.0;  // mol/m^3
  double K_eq_small = k_forward_small / k_reverse_small;  // 2.0
  double A_eq_small = A0_small / (1.0 + K_eq_small);
  double B_eq_small = A0_small * K_eq_small / (1.0 + K_eq_small);
  double k_total_small = k_forward_small + k_reverse_small;  // 0.15 s^-1
  double tau_small = 1.0 / k_total_small;  // 6.67 s
  
  // Test parameters for large droplets
  double A0_large = 0.5;  // mol/m^3
  double K_eq_large = k_forward_large / k_reverse_large;  // 2.0 (same as small)
  double A_eq_large = A0_large / (1.0 + K_eq_large);
  double B_eq_large = A0_large * K_eq_large / (1.0 + K_eq_large);
  double k_total_large = k_forward_large + k_reverse_large;  // 0.3 s^-1
  double tau_large = 1.0 / k_total_large;  // 3.33 s
  
  // Build solver
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModelProcesses(model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();
  
  // Initialize state
  State state = solver.GetState();
  
  // Find variable indices for small droplets
  std::size_t i_A_small = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_SMALL.AQUEOUS_SMALL.A") - state.variable_names_.begin();
  std::size_t i_B_small = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_SMALL.AQUEOUS_SMALL.B") - state.variable_names_.begin();
  std::size_t i_C_small = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_SMALL.AQUEOUS_SMALL.C") - state.variable_names_.begin();
  
  // Find variable indices for large droplets
  std::size_t i_A_large = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_LARGE.AQUEOUS_LARGE.A") - state.variable_names_.begin();
  std::size_t i_B_large = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_LARGE.AQUEOUS_LARGE.B") - state.variable_names_.begin();
  std::size_t i_C_large = std::find(state.variable_names_.begin(), state.variable_names_.end(),
                                     "DROPLET_LARGE.AQUEOUS_LARGE.C") - state.variable_names_.begin();
  
  // Set initial conditions
  state.variables_[0][i_A_small] = A0_small;
  state.variables_[0][i_B_small] = 0.0;
  state.variables_[0][i_C_small] = 1.0e-4;  // solvent
  
  state.variables_[0][i_A_large] = A0_large;
  state.variables_[0][i_B_large] = 0.0;
  state.variables_[0][i_C_large] = 1.0e-4;  // solvent
  
  small_droplet.SetDefaultParameters(state);
  large_droplet.SetDefaultParameters(state);
  
  // Tolerance for numerical comparisons
  const double tolerance = 1.0e-4;  // mol/m³
  
  // Test integration at multiple times
  // Use 100s to reach equilibrium (about 15×tau_small)
  std::vector<double> test_times = { 5.0, 10.0, 20.0, 100.0 };
  
  solver.CalculateRateConstants(state);
  double time = 0.0;
  
  for (double target_time : test_times) {
    // Integrate to target time
    while (time < target_time - 1.0e-10) {
      double dt = std::min(0.01, target_time - time);  // 0.01s timestep
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged);
      time += dt;
    }
    
    // Get numerical solution for small droplets
    double A_numeric_small = state.variables_[0][i_A_small];
    double B_numeric_small = state.variables_[0][i_B_small];
    double C_numeric_small = state.variables_[0][i_C_small];
    
    // Get numerical solution for large droplets
    double A_numeric_large = state.variables_[0][i_A_large];
    double B_numeric_large = state.variables_[0][i_B_large];
    double C_numeric_large = state.variables_[0][i_C_large];
    
    // Analytical solutions for small droplets
    double exp_term_small = std::exp(-k_total_small * time);
    double A_analytic_small = A_eq_small + (A0_small - A_eq_small) * exp_term_small;
    double B_analytic_small = B_eq_small + (0.0 - B_eq_small) * exp_term_small;
    
    // Analytical solutions for large droplets
    double exp_term_large = std::exp(-k_total_large * time);
    double A_analytic_large = A_eq_large + (A0_large - A_eq_large) * exp_term_large;
    double B_analytic_large = B_eq_large + (0.0 - B_eq_large) * exp_term_large;
    
    // Check small droplets
    EXPECT_NEAR(A_numeric_small, A_analytic_small, tolerance)
      << "Small droplet A mismatch at t=" << time << "s";
    EXPECT_NEAR(B_numeric_small, B_analytic_small, tolerance)
      << "Small droplet B mismatch at t=" << time << "s";
    EXPECT_NEAR(A_numeric_small + B_numeric_small, A0_small, tolerance)
      << "Small droplet mass not conserved at t=" << time << "s";
    EXPECT_NEAR(C_numeric_small, 1.0e-4, tolerance)
      << "Small droplet solvent changed at t=" << time << "s";
    
    // Check large droplets
    EXPECT_NEAR(A_numeric_large, A_analytic_large, tolerance)
      << "Large droplet A mismatch at t=" << time << "s";
    EXPECT_NEAR(B_numeric_large, B_analytic_large, tolerance)
      << "Large droplet B mismatch at t=" << time << "s";
    EXPECT_NEAR(A_numeric_large + B_numeric_large, A0_large, tolerance)
      << "Large droplet mass not conserved at t=" << time << "s";
    EXPECT_NEAR(C_numeric_large, 1.0e-4, tolerance)
      << "Large droplet solvent changed at t=" << time << "s";
    
    // Check that both reach same equilibrium fraction (K_eq is same)
    double frac_B_small = B_numeric_small / (A_numeric_small + B_numeric_small);
    double frac_B_large = B_numeric_large / (A_numeric_large + B_numeric_large);
    double expected_frac = K_eq_small / (1.0 + K_eq_small);  // 2/3
    
    if (time >= 50.0) {  // Check at equilibrium
      EXPECT_NEAR(frac_B_small, expected_frac, 0.01)
        << "Small droplet equilibrium fraction wrong at t=" << time << "s";
      EXPECT_NEAR(frac_B_large, expected_frac, 0.01)
        << "Large droplet equilibrium fraction wrong at t=" << time << "s";
    }
  }
  
  // Verify no cross-coupling: large droplets should equilibrate faster
  // At t=10s (3×tau_large but only 1.5×tau_small):
  // Large droplets: exp(-0.3×10) = exp(-3) ≈ 0.05 (95% to equilibrium)
  // Small droplets: exp(-0.15×10) = exp(-1.5) ≈ 0.22 (78% to equilibrium)
  // So large should be closer to equilibrium
}
