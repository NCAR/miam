// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Integration tests for HenryLawEquilibriumConstraint with the MICM DAE solver.

#include <miam/miam.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;
using namespace miam;

namespace
{
  constexpr double R_gas = miam::util::R_gas;  // 8.314462618 J mol⁻¹ K⁻¹
}

// ============================================================================
// Test 1: Henry's Law equilibrium constraint with gas-phase driver
//
// System:
//   Precursor → A_g       (MICM gas-phase reaction, rate k)
//   A_g <-> A_aq           (HL equilibrium constraint, A_aq algebraic)
//   [Precursor] + [A_g] + [A_aq] = total  (mass conservation, A_g algebraic)
//
// Differential: Precursor with dP/dt = -k*[P]
// Algebraic:    A_g (mass conservation), A_aq (HL equilibrium)
//
// Let  α = HLC * R * T * f_v  (dimensionless equilibrium ratio)
//      f_v = [H2O] * Mw_H2O / rho_H2O
//
// Analytical:
//   P(t)     = P0 * exp(-k*t)
//   A_g(t)   = (total - P(t)) / (1 + α)
//   A_aq(t)  = α * (total - P(t)) / (1 + α)
// ============================================================================
TEST(HenryLawEquilibriumConstraintIntegration, GasPhaseDriverSingleInstance)
{
  // Physical parameters
  double T = 298.15;              // K
  double HLC = 4.0e-4;            // mol m⁻³ Pa⁻¹ (constant, no T dependence)
  double Mw_solvent = 0.018;      // kg mol⁻¹ (water)
  double rho_solvent = 1000.0;    // kg m⁻³ (water)
  double H2O_conc = rho_solvent / Mw_solvent;  // 55555.6 mol/m³ → f_v ≈ 1.0
  double f_v = H2O_conc * Mw_solvent / rho_solvent;  // = 1.0
  double alpha = HLC * R_gas * T * f_v;               // dimensionless

  // Species
  auto Precursor = Species{ "Precursor" };
  auto A_g = Species{ "A_g" };
  auto A_aq = Species{ "A_aq" };
  auto H2O = Species{ "H2O",
      { { "molecular weight [kg mol-1]", Mw_solvent },
        { "density [kg m-3]", rho_solvent } } };

  // Phases
  Phase gas_phase{ "GAS", { { Precursor }, { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  // Representation
  auto droplet = representation::UniformSection{
    "DROPLET",
    { aqueous_phase }
  };

  // MICM gas-phase reaction: Precursor → A_g with rate k
  double k = 0.05;  // s⁻¹
  Process gas_rxn = ChemicalReactionBuilder()
      .SetReactants({ Precursor })
      .SetProducts({ StoichSpecies(A_g, 1.0) })
      .SetRateConstant(ArrheniusRateConstant({ .A_ = k }))
      .SetPhase(gas_phase)
      .Build();

  // HL equilibrium constraint: A_g <-> A_aq (A_aq algebraic, per instance)
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

  // Mass conservation: [Precursor] + [A_g] + [A_aq] = total, A_g algebraic (global)
  double P0 = 1.0;   // mol/m³ initial [Precursor]
  double total = P0;  // total mass (A_g0 = A_aq0 = 0)

  auto mass_cons = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, A_g)
      .AddTerm(gas_phase, Precursor, 1.0)
      .AddTerm(gas_phase, A_g, 1.0)
      .AddTerm(aqueous_phase, A_aq, 1.0)
      .SetConstant(total)
      .Build();

  // Build MIAM model (no MIAM kinetic processes, only constraints)
  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { droplet }
  };
  model.AddConstraints(hl_constraint, mass_cons);

  // Build DAE solver with MICM gas-phase reaction + MIAM constraints
  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .SetReactions({ gas_rxn })
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  auto find_idx = [&](const std::string& name) {
    auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
    EXPECT_NE(it, state.variable_names_.end()) << "Species " << name << " not found";
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  std::size_t i_P = find_idx("Precursor");
  std::size_t i_gas = find_idx("A_g");
  std::size_t i_aq = find_idx("DROPLET.AQUEOUS.A_aq");
  std::size_t i_h2o = find_idx("DROPLET.AQUEOUS.H2O");

  // Set initial conditions
  state.variables_[0][i_P] = P0;
  state.variables_[0][i_gas] = 0.0;
  state.variables_[0][i_aq] = 0.0;
  state.variables_[0][i_h2o] = H2O_conc;

  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  // Time integration (1/k = 20s, 5*tau = 100s)
  std::vector<double> test_times = { 1.0, 5.0, 10.0, 120.0 };
  double time = 0.0;
  const double tolerance = 5.0e-3;

  solver.UpdateStateParameters(state);

  for (double target_time : test_times)
  {
    while (time < target_time - 1.0e-10)
    {
      double dt = std::min(0.1, target_time - time);
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time << " s";
      time += dt;
    }

    double P_num = state.variables_[0][i_P];
    double gas_num = state.variables_[0][i_gas];
    double aq_num = state.variables_[0][i_aq];

    // Analytical solution
    double P_an = P0 * std::exp(-k * time);
    double gas_an = (total - P_an) / (1.0 + alpha);
    double aq_an = alpha * (total - P_an) / (1.0 + alpha);

    EXPECT_NEAR(P_num, P_an, tolerance)
      << "Precursor mismatch at t = " << time;
    EXPECT_NEAR(gas_num, gas_an, tolerance)
      << "A_g mismatch at t = " << time;
    EXPECT_NEAR(aq_num, aq_an, tolerance)
      << "A_aq mismatch at t = " << time;

    // Verify mass conservation
    double mass = P_num + gas_num + aq_num;
    EXPECT_NEAR(mass, total, 1.0e-6)
      << "Mass conservation violated at t = " << time;

    // Verify HL equilibrium relation: A_aq = α * A_g
    if (gas_num > 1.0e-10)
    {
      EXPECT_NEAR(aq_num / gas_num, alpha, 1.0e-2)
        << "HL equilibrium ratio violated at t = " << time;
    }
  }

  // Steady state: Precursor is depleted (at t=120s, exp(-0.05*120) ≈ 0.0025)
  EXPECT_NEAR(state.variables_[0][i_P], 0.0, 5.0e-3);
  EXPECT_NEAR(state.variables_[0][i_gas], total / (1.0 + alpha), tolerance);
  EXPECT_NEAR(state.variables_[0][i_aq], total * alpha / (1.0 + alpha), tolerance);
}

// ============================================================================
// Test 2: HL equilibrium with multiple condensed phase instances
//
// One gas species A_g shared by two droplet modes (SMALL, LARGE).
// Each mode has its own A_aq and H2O.
// The HL constraint generates one algebraic equation per instance.
// ============================================================================
TEST(HenryLawEquilibriumConstraintIntegration, MultipleInstances)
{
  double T = 298.15;
  double HLC = 4.0e-4;
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double H2O_conc = rho_solvent / Mw_solvent;
  double f_v = H2O_conc * Mw_solvent / rho_solvent;
  double alpha = HLC * R_gas * T * f_v;

  auto Precursor = Species{ "Precursor" };
  auto A_g = Species{ "A_g" };
  auto A_aq = Species{ "A_aq" };
  auto H2O = Species{ "H2O",
      { { "molecular weight [kg mol-1]", Mw_solvent },
        { "density [kg m-3]", rho_solvent } } };

  Phase gas_phase{ "GAS", { { Precursor }, { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto small_drop = representation::UniformSection{ "SMALL", { aqueous_phase } };
  auto large_drop = representation::UniformSection{ "LARGE", { aqueous_phase } };

  double k = 0.05;
  Process gas_rxn = ChemicalReactionBuilder()
      .SetReactants({ Precursor })
      .SetProducts({ StoichSpecies(A_g, 1.0) })
      .SetRateConstant(ArrheniusRateConstant({ .A_ = k }))
      .SetPhase(gas_phase)
      .Build();

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

  double P0 = 1.0;
  // With 2 instances: total = [P] + [A_g] + [A_aq_SMALL] + [A_aq_LARGE]
  // At equilibrium: [A_aq_i] = α * [A_g] for each i
  // total = [P] + [A_g] + 2*α*[A_g] = [P] + [A_g]*(1 + 2*α)
  double total = P0;

  auto mass_cons = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, A_g)
      .AddTerm(gas_phase, Precursor, 1.0)
      .AddTerm(gas_phase, A_g, 1.0)
      .AddTerm(aqueous_phase, A_aq, 1.0)  // expands to SMALL + LARGE
      .SetConstant(total)
      .Build();

  auto model = Model{
    .name_ = "AEROSOL",
    .representations_ = { small_drop, large_drop }
  };
  model.AddConstraints(hl_constraint, mass_cons);

  auto system = System(gas_phase, model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                    .SetSystem(system)
                    .SetReactions({ gas_rxn })
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  auto find_idx = [&](const std::string& name) {
    auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
    EXPECT_NE(it, state.variable_names_.end()) << "Species " << name << " not found";
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  std::size_t i_P = find_idx("Precursor");
  std::size_t i_gas = find_idx("A_g");
  std::size_t i_aq_small = find_idx("SMALL.AQUEOUS.A_aq");
  std::size_t i_h2o_small = find_idx("SMALL.AQUEOUS.H2O");
  std::size_t i_aq_large = find_idx("LARGE.AQUEOUS.A_aq");
  std::size_t i_h2o_large = find_idx("LARGE.AQUEOUS.H2O");

  state.variables_[0][i_P] = P0;
  state.variables_[0][i_gas] = 0.0;
  state.variables_[0][i_aq_small] = 0.0;
  state.variables_[0][i_h2o_small] = H2O_conc;
  state.variables_[0][i_aq_large] = 0.0;
  state.variables_[0][i_h2o_large] = H2O_conc;

  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 101325.0;

  small_drop.SetDefaultParameters(state);
  large_drop.SetDefaultParameters(state);

  // Integrate to steady state (1/k = 20s, 5*tau = 100s)
  double time = 0.0;
  while (time < 120.0)
  {
    double dt = 0.1;
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged)
      << "Solver failed at t = " << time;
    time += dt;
  }

  const double tolerance = 5.0e-3;

  // Steady state: P → 0, A_g = total / (1 + 2*α), A_aq_i = α * A_g
  double gas_expected = total / (1.0 + 2.0 * alpha);
  double aq_expected = alpha * gas_expected;

  EXPECT_NEAR(state.variables_[0][i_P], 0.0, 5.0e-3);
  EXPECT_NEAR(state.variables_[0][i_gas], gas_expected, tolerance);
  EXPECT_NEAR(state.variables_[0][i_aq_small], aq_expected, tolerance);
  EXPECT_NEAR(state.variables_[0][i_aq_large], aq_expected, tolerance);

  // Both instances should have the same A_aq (same α, same A_g)
  EXPECT_NEAR(state.variables_[0][i_aq_small], state.variables_[0][i_aq_large], 1.0e-6);

  // Mass conservation
  double mass = state.variables_[0][i_P] + state.variables_[0][i_gas]
              + state.variables_[0][i_aq_small] + state.variables_[0][i_aq_large];
  EXPECT_NEAR(mass, total, 1.0e-4);
}
