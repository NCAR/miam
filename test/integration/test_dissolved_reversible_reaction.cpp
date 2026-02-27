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