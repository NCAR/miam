// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Integration test: Aqueous carbonic acid system (based on CarbonicAcidDAE.pdf)
//
// Chemistry:
//   CO2(g) <-> CO2(aq)                    [Henry's Law, K_H = 3.4e-2 M/atm]
//   CO2(aq) + H2O <-> H+ + HCO3-         [K1 = 4.3e-7 M]
//   HCO3-   <-> H+ + CO3--               [K2 = 4.7e-11 M]
//   H2O     <-> H+ + OH-                  [Kw  = 1.0e-14 M2]
//
// Test 1 - KineticODE:
//   All five aqueous variables (CO2_aq, HCO3-, CO3--, H+, OH-) are
//   differential.  Forward/reverse kinetics drive the system to equilibrium.
//   Solver: ThreeStageRosenbrock (ODE).
//
// Test 2 - DAEConstraints:
//   CO2(aq) is algebraic via HenryLawEquilibriumConstraint (instantly in
//   equilibrium with CO2(g)).  H+ is algebraic via charge balance
//   (LinearConstraint).  OH- is algebraic via Kw equilibrium
//   (DissolvedEquilibriumConstraint).  HCO3- and CO3-- are differential,
//   driven by rxn1 and a slow rxn2.
//   Solver: FourStageDifferentialAlgebraicRosenbrock (DAE).
//
// Both tests verify that equilibrium is reached at 400 ppm atmospheric CO2,
// T = 298.15 K, P = 101325 Pa.
//
// UNIT CONVENTIONS (mol m-3 air, same as test_cam_cloud_chemistry.cpp):
//   C_H2O = 0.017 mol m-3 air  (cloud LWC ~0.3 g m-3)
//   f_v   = C_H2O * Mw/rho     (volume fraction of liquid water ~3.06e-7)
//   K_miam = K_lit[M] / c_H2O^(n_p - n_r),  c_H2O = 55.556 mol/L
//   K_miam_water = Kw_lit / c_H2O^2  (H2O appears as explicit reactant AND solvent)

#include <miam/miam.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <miam/util/constants.hpp>

#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <format>
#include <iostream>
#include <string>
#include <utility>

using namespace micm;
using namespace miam;

namespace
{
  // -- Physical constants ------------------------------------------------------
  constexpr double T = 298.15;          // K
  constexpr double P = 101325.0;        // Pa
  constexpr double Mw_water = 0.018;    // kg mol-1
  constexpr double rho_water = 1000.0;  // kg m-3

  // -- Cloud liquid water content ----------------------------------------------
  constexpr double C_H2O = 0.017;     // mol m-3 air (LWC ~0.3 g m-3)
  constexpr double c_H2O_M = 55.556;  // mol/L (pure-water molar conc.)

  // -- Henry's Law constant (CO2) -----------------------------------------------
  // K_H = 3.4e-2 M/atm = 3.4e-2 * 1000/101325 mol m-3_liq Pa-1
  constexpr double K_H_lit = 3.4e-2;                   // M/atm
  constexpr double K_H = K_H_lit * 1000.0 / 101325.0;  // mol m-3_liq Pa-1

  // -- Equilibrium constants (MIAM solvent-normalised, dimensionless) ----------
  // K_miam = K_lit[M] / c_H2O^(n_p - n_r)
  //   K1: CO2_aq <-> H+ + HCO3-   (n_r=1, n_p=2) -> / c_H2O^1
  constexpr double K1_miam = 4.3e-7 / c_H2O_M;  // ~7.74e-9
  //   K2: HCO3- <-> H+ + CO3--    (n_r=1, n_p=2) -> / c_H2O^1
  constexpr double K2_miam = 4.7e-11 / c_H2O_M;  // ~8.46e-13
  //   Kw: H2O(reactant) <-> H+ + OH-  (n_r=1 for H2O, n_p=2) -> / c_H2O^2
  //   (H2O appears as explicit reactant AND solvent; net exponent is 2)
  constexpr double Kw_miam = 1.0e-14 / (c_H2O_M * c_H2O_M);  // ~3.24e-18

  // -- CO2(g) at 400 ppm: ideal-gas concentration in mol m-3 air ---------------
  constexpr double ppm_CO2 = 400.0e-6;
  constexpr double P_CO2 = ppm_CO2 * P;                   // Pa
  constexpr double CO2g_eq = P_CO2 / (GAS_CONSTANT * T);  // ~0.01635 mol m-3 air

  // -- Kinetic rate constants (reasonable defaults; user should tune) ----------
  // CO2_aq + H2O <-> H+ + HCO3- :  k1_f/k1_r = K1_miam
  constexpr double k1_f = 0.1;             // s-1
  constexpr double k1_r = k1_f / K1_miam;  // ~1.29e7 s-1
  // HCO3- <-> H+ + CO3-- :          k2_f/k2_r = K2_miam
  constexpr double k2_f = 5.0;             // s-1
  constexpr double k2_r = k2_f / K2_miam;  // ~5.91e12 s-1
  // H2O <-> H+ + OH- :              kw_f/kw_r = Kw_miam  (ODE test only)
  constexpr double kw_f = 1.0e-6;          // s-1 (effective, via solvent)
  constexpr double kw_r = kw_f / Kw_miam;  // ~3.09e11 s-1

  // -- Phase-transfer parameters (CO2) ----------------------------------------
  constexpr double D_CO2 = 1.5e-5;  // m2 s-1  (gas-phase diffusivity)
  constexpr double alpha = 0.05;    // accommodation coefficient

  // -- Shared species, phases, representation, and rxn1 -----------------------
  // Both tests use the same 7 species, the same two phases, the same droplet
  // representation, and the same first reversible reaction (rxn1).
  struct CarbonicAcidSystem
  {
    Species CO2_g, CO2_aq, HCO3m, CO3mm, Hp, OHm, H2O;
    Phase gas_phase, aqueous_phase;
    SingleMomentMode droplet;
    DissolvedReversibleReaction rxn1;
  };

  CarbonicAcidSystem MakeCarbonicAcidSystem()
  {
    // All condensed-phase species carry MW and density for SingleMomentMode.
    Species CO2_g{ "CO2_g", { { "molecular weight [kg mol-1]", 0.044 } } };
    Species CO2_aq{ "CO2_aq", { { "molecular weight [kg mol-1]", 0.044 }, { "density [kg m-3]", 1000.0 } } };
    Species HCO3m{ "HCO3-", { { "molecular weight [kg mol-1]", 0.061 }, { "density [kg m-3]", 1000.0 } } };
    Species CO3mm{ "CO3--", { { "molecular weight [kg mol-1]", 0.060 }, { "density [kg m-3]", 1000.0 } } };
    Species Hp{ "H+", { { "molecular weight [kg mol-1]", 0.001 }, { "density [kg m-3]", 1000.0 } } };
    Species OHm{ "OH-", { { "molecular weight [kg mol-1]", 0.017 }, { "density [kg m-3]", 1000.0 } } };
    Species H2O{ "H2O", { { "molecular weight [kg mol-1]", Mw_water }, { "density [kg m-3]", rho_water } } };

    Phase gas_phase{ "GAS", { { CO2_g } } };
    Phase aqueous_phase{ "AQUEOUS", { { CO2_aq }, { HCO3m }, { CO3mm }, { Hp }, { OHm }, { H2O } } };
    auto droplet = SingleMomentMode{ "DROPLET", { aqueous_phase }, 5.0e-6, 1.2 };

    // CO2(aq) + H2O <-> H+ + HCO3-  (same rate constants in both tests)
    auto rxn1 = DissolvedReversibleReaction{ [](const Conditions&) { return k1_f; },
                                             [](const Conditions&) { return k1_r; },
                                             { CO2_aq },
                                             { Hp, HCO3m },
                                             H2O,
                                             aqueous_phase };

    return CarbonicAcidSystem{
      .CO2_g = CO2_g,
      .CO2_aq = CO2_aq,
      .HCO3m = HCO3m,
      .CO3mm = CO3mm,
      .Hp = Hp,
      .OHm = OHm,
      .H2O = H2O,
      .gas_phase = gas_phase,
      .aqueous_phase = aqueous_phase,
      .droplet = std::move(droplet),
      .rxn1 = std::move(rxn1),
    };
  }

  // -- Diagnostic output -------------------------------------------------------
  void PrintConditions(
      const std::string& label,
      double CO2g,
      double CO2_aq,
      double HCO3,
      double CO3,
      double H,
      double OH,
      double H2O)
  {
    std::cout << std::format("{} (mol m-3 air):\n", label) << std::format("  CO2_g  = {:10.3e}\n", CO2g)
              << std::format("  CO2_aq = {:10.3e}\n", CO2_aq) << std::format("  HCO3-  = {:10.3e}\n", HCO3)
              << std::format("  CO3--  = {:10.3e}\n", CO3) << std::format("  H+     = {:10.3e}\n", H)
              << std::format("  OH-    = {:10.3e}\n", OH) << std::format("  H2O    = {:10.3e}\n", H2O);
  }

  void PrintEquilibriumStatus(double CO2_aq, double HCO3, double CO3, double H, double OH, double H2O)
  {
    std::cout << "Equilibrium status (% deviation from 1.0):\n";
    if (CO2_aq > 1e-30 && H2O > 0 && H > 0 && HCO3 > 0)
    {
      double r = (H * HCO3) / (CO2_aq * H2O * K1_miam);
      std::cout << std::format("  K1  [H+][HCO3-] / ([CO2_aq][H2O] K1_miam) = {:.4f}  ({:+.2f}%)\n", r, (r - 1.0) * 100.0);
    }
    if (HCO3 > 1e-30 && H2O > 0 && H > 0 && CO3 > 0)
    {
      double r = (H * CO3) / (HCO3 * H2O * K2_miam);
      std::cout << std::format("  K2  [H+][CO3--] / ([HCO3-][H2O]  K2_miam) = {:.4f}  ({:+.2f}%)\n", r, (r - 1.0) * 100.0);
    }
    if (H > 1e-30 && OH > 1e-30 && H2O > 0)
    {
      double r = (H * OH) / (H2O * H2O * Kw_miam);
      std::cout << std::format("  Kw  [H+][OH-]   / ([H2O]^2      Kw_miam)  = {:.4f}  ({:+.2f}%)\n", r, (r - 1.0) * 100.0);
    }
    if (H > 1e-30)
    {
      double r = H / (OH + HCO3 + 2.0 * CO3);
      std::cout << std::format("  CB  [H+] / ([OH-] + [HCO3-] + 2[CO3--])   = {:.4f}  ({:+.2f}%)\n", r, (r - 1.0) * 100.0);
    }
  }

  // -- Equilibrium check -------------------------------------------------------
  // Verify dimensionless equilibrium ratios and charge balance are within rtol.
  // All concentrations in mol m-3 air.
  void CheckEquilibrium(double CO2_aq, double HCO3, double CO3, double H, double OH, double H2O, double rtol = 0.05)
  {
    // K1: K1_miam = [H+][HCO3-] / ([CO2_aq][H2O])
    if (CO2_aq > 1e-30 && H2O > 0 && H > 0 && HCO3 > 0)
    {
      double ratio = (H * HCO3) / (CO2_aq * H2O * K1_miam);
      EXPECT_NEAR(ratio, 1.0, rtol) << "K1 equilibrium not satisfied";
    }

    // K2: K2_miam = [H+][CO3--] / ([HCO3-][H2O])
    if (HCO3 > 1e-30 && H2O > 0 && H > 0 && CO3 > 0)
    {
      double ratio = (H * CO3) / (HCO3 * H2O * K2_miam);
      EXPECT_NEAR(ratio, 1.0, rtol) << "K2 equilibrium not satisfied";
    }

    // Kw: Kw_miam = [H+][OH-] / [H2O]^2
    if (H > 1e-30 && OH > 1e-30 && H2O > 0)
    {
      double ratio = (H * OH) / (H2O * H2O * Kw_miam);
      EXPECT_NEAR(ratio, 1.0, rtol) << "Kw equilibrium not satisfied";
    }

    // Charge balance: [H+] = [OH-] + [HCO3-] + 2[CO3--]
    if (H > 1e-30)
    {
      double rhs = OH + HCO3 + 2.0 * CO3;
      double cb_ratio = H / rhs;
      EXPECT_NEAR(cb_ratio, 1.0, rtol) << "Charge balance not satisfied";
    }
  }
}  // namespace

// ============================================================================
// Test 1: Fully kinetic ODE system
//
// All aqueous species (CO2_aq, HCO3-, CO3--, H+, OH-) are differential.
// HenryLawPhaseTransfer drives CO2 dissolution; three DissolvedReversibleReactions
// drive the aqueous equilibria.  The ThreeStageRosenbrock ODE solver is used.
//
// Verification: at t = 3000 s (>> tau_1 ~10 s), all equilibrium ratios (K1, K2,
// Kw) and the charge balance must be satisfied to within 5 % relative error.
// ============================================================================
TEST(AqueousCarbonicAcid, KineticODE)
{
  auto sys = MakeCarbonicAcidSystem();

  // -- ODE-specific processes --
  // Henry's Law kinetic phase transfer: CO2(g) <-> CO2(aq)
  auto transfer = HenryLawPhaseTransferBuilder()
                      .SetCondensedPhase(sys.aqueous_phase)
                      .SetGasSpecies(sys.CO2_g)
                      .SetCondensedSpecies(sys.CO2_aq)
                      .SetSolvent(sys.H2O)
                      .SetHenrysLawConstant(HenrysLawConstant(HenrysLawConstantParameters{ .HLC_ref_ = K_H }))
                      .SetDiffusionCoefficient(D_CO2)
                      .SetAccommodationCoefficient(alpha)
                      .Build();

  // HCO3- <-> H+ + CO3--  (full stiff rate; stable for pure-ODE Rosenbrock)
  auto rxn2 = DissolvedReversibleReaction{ [](const Conditions&) { return k2_f; },
                                           [](const Conditions&) { return k2_r; },
                                           { sys.HCO3m },
                                           { sys.Hp, sys.CO3mm },
                                           sys.H2O,
                                           sys.aqueous_phase };

  // H2O <-> H+ + OH-  (H2O as reactant AND solvent gives correct Kw_miam)
  auto rxn_w = DissolvedReversibleReaction{ [](const Conditions&) { return kw_f; },
                                            [](const Conditions&) { return kw_r; },
                                            { sys.H2O },
                                            { sys.Hp, sys.OHm },
                                            sys.H2O,
                                            sys.aqueous_phase };

  // -- Model --
  auto model = Model{ .name_ = "CARBONIC_ACID_ODE", .representations_ = { sys.droplet } };
  model.AddProcesses({ transfer });
  model.AddProcesses({ sys.rxn1, rxn2, rxn_w });

  // -- ODE Solver --
  auto system = System(sys.gas_phase);
  auto ode_params = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  ode_params.h_start_ = 1.0e-4;  // small initial step; Rosenbrock will grow it adaptively
  ode_params.h_max_ = 5.0;       // cap step size: k2_r*h_max ~3e13, acceptable for rtol=0.05
  ode_params.max_number_of_steps_ = 1000000;
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(ode_params)
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  // -- Variable indices --
  std::size_t i_CO2g = state.variable_map_.at("CO2_g");
  std::size_t i_CO2aq = state.variable_map_.at("DROPLET.AQUEOUS.CO2_aq");
  std::size_t i_HCO3 = state.variable_map_.at("DROPLET.AQUEOUS.HCO3-");
  std::size_t i_CO3 = state.variable_map_.at("DROPLET.AQUEOUS.CO3--");
  std::size_t i_H = state.variable_map_.at("DROPLET.AQUEOUS.H+");
  std::size_t i_OH = state.variable_map_.at("DROPLET.AQUEOUS.OH-");
  std::size_t i_H2O = state.variable_map_.at("DROPLET.AQUEOUS.H2O");

  // -- Initial conditions: start far from equilibrium --
  state.variables_[0][i_CO2g] = CO2g_eq;  // initialize CO2(g) to 400 ppm
  state.variables_[0][i_CO2aq] = 0.0;
  state.variables_[0][i_HCO3] = 0.0;
  state.variables_[0][i_CO3] = 0.0;
  state.variables_[0][i_H] = 1.0e-12;  // small seed to trigger H2O autodissociation
  state.variables_[0][i_OH] = 1.0e-12;
  state.variables_[0][i_H2O] = C_H2O;

  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = P;
  sys.droplet.SetDefaultParameters(state);

  std::cout << "\n=== KineticODE ===\n";
  PrintConditions(
      "Initial",
      state.variables_[0][i_CO2g],
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);

  // -- Integrate to equilibrium --
  // 3000 s >> tau_1 = 1/k1_f = 10 s; sufficient for all modes to equilibrate.
  solver.UpdateStateParameters(state);
  constexpr double t_end_ode = 3000.0;
  auto ode_result = solver.Solve(t_end_ode, state);
  ASSERT_EQ(ode_result.state_, SolverState::Converged) << "ODE solver failed";
  std::cout << "Solver steps: " << ode_result.stats_.number_of_steps_ << "\n";

  PrintConditions(
      "Final",
      state.variables_[0][i_CO2g],
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);
  PrintEquilibriumStatus(
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);

  // -- Verify equilibrium --
  const double co2_g_final = state.variables_[0][i_CO2g];
  const double co2_aq_final = state.variables_[0][i_CO2aq];
  const double hco3_final = state.variables_[0][i_HCO3];
  const double co3_final = state.variables_[0][i_CO3];
  const double h_final = state.variables_[0][i_H];
  const double oh_final = state.variables_[0][i_OH];
  const double h2o_final = state.variables_[0][i_H2O];

  // Non-negativity: all aqueous species must be >= 0 after solve
  EXPECT_GE(co2_aq_final, 0.0) << "CO2_aq went negative";
  EXPECT_GE(hco3_final, 0.0) << "HCO3- went negative";
  EXPECT_GE(co3_final, 0.0) << "CO3-- went negative";
  EXPECT_GE(h_final, 0.0) << "H+ went negative";
  EXPECT_GE(oh_final, 0.0) << "OH- went negative";
  EXPECT_GE(h2o_final, 0.0) << "H2O went negative";

  // Assertion guards: ensure ratio checks inside CheckEquilibrium actually run
  ASSERT_GT(hco3_final, 0.0) << "HCO3- is zero — K1/K2 ratio checks would be silently skipped";
  ASSERT_GT(co3_final, 0.0) << "CO3-- is zero — K2 ratio check would be silently skipped";
  ASSERT_GT(h_final, 0.0) << "H+ is zero — equilibrium ratio checks would be silently skipped";
  ASSERT_GT(oh_final, 0.0) << "OH- is zero — Kw ratio check would be silently skipped";

  // pH check: 400 ppm CO2 should give pH ≈ 5.6, i.e. [H+] ≈ 7.4e-10 mol m-3 air
  EXPECT_NEAR(h_final / 7.4e-10, 1.0, 0.1) << "[H+] not near expected 7.4e-10 mol m-3 air (pH 5.6)";

  // Henry's Law
  const double f_v_final = h2o_final * Mw_water / rho_water;
  const double henry_expected = K_H * GAS_CONSTANT * T * f_v_final * co2_g_final;
  ASSERT_GT(henry_expected, 0.0) << "Expected positive Henry-law CO2(aq) concentration";
  EXPECT_NEAR(co2_aq_final, henry_expected, 0.05 * henry_expected)
      << std::format("Henry's-law equilibrium mismatch: CO2_aq={} expected={}", co2_aq_final, henry_expected);

  // Carbon mass conservation: total C (gas + dissolved) must be invariant
  // Initial: all C is CO2_g; aq species start at 0
  const double C_total_init = CO2g_eq;
  const double C_total_final = co2_g_final + co2_aq_final + hco3_final + co3_final;
  EXPECT_NEAR(C_total_final / C_total_init, 1.0, 1.0e-5) << "Carbon mass not conserved";

  CheckEquilibrium(co2_aq_final, hco3_final, co3_final, h_final, oh_final, h2o_final);
}

// ============================================================================
// Test 2: DAE system - CO2(aq), CO3--, OH-, and H+ as algebraic variables
//
// CO2(g) is held at 400 ppm (no gas-phase processes).
// CO2(aq) is algebraic via HenryLawEquilibriumConstraint (instant HL equilibrium).
// CO3-- is algebraic via DissolvedEquilibriumConstraint (K2): HCO3- <-> H+ + CO3--.
// OH- is algebraic via DissolvedEquilibriumConstraint (Kw): H2O <-> H+ + OH-.
// H+ is algebraic via a LinearConstraint encoding the charge balance:
//   [H+] - [OH-] - [HCO3-] - 2[CO3--] = 0.
// HCO3- is the only differential aqueous species, driven by rxn1.
// Absolute tolerances are tightened for all trace species (H+, OH-, CO3--, CO2_aq,
// HCO3-) because the default 1e-3 mol m-3 >> their equilibrium values
// (~1e-9 to ~1e-14 mol m-3), which would allow the solver to accept steps that
// drive them negative.
// Solver: FourStageDifferentialAlgebraicRosenbrock (DAE).
//
// Verification: same equilibrium checks as Test 1, plus the algebraic
// constraints (K2, Kw, charge balance, and Henry's law) must be satisfied to
// <= 1e-4 relative error (tight, because they are enforced exactly by the DAE solver).
// ============================================================================
TEST(AqueousCarbonicAcid, DAEConstraints)
{
  auto sys = MakeCarbonicAcidSystem();

  // -- Algebraic constraints --

  // CO3--: K2 equilibrium  HCO3- <-> H+ + CO3--
  //   G = K2_miam * [HCO3-] * [H2O] / ([H2O]+δ) - [H+] * [CO3--] * [H2O] / ([H2O]+δ)^2 = 0
  auto k2_constraint = DissolvedEquilibriumConstraintBuilder()
                           .SetPhase(sys.aqueous_phase)
                           .SetReactants({ sys.HCO3m })
                           .SetProducts({ sys.Hp, sys.CO3mm })
                           .SetAlgebraicSpecies(sys.CO3mm)
                           .SetSolvent(sys.H2O)
                           .SetEquilibriumConstant(EquilibriumConstant(EquilibriumConstantParameters{ .A_ = K2_miam }))
                           .Build();

  // CO2(aq): Henry's Law equilibrium with CO2(g)
  //   G = K_H * R * T * f_v * [CO2_g] - [CO2_aq] = 0
  auto hl_constraint = HenryLawEquilibriumConstraintBuilder()
                           .SetCondensedPhase(sys.aqueous_phase)
                           .SetGasSpecies(sys.CO2_g)
                           .SetCondensedSpecies(sys.CO2_aq)
                           .SetSolvent(sys.H2O)
                           .SetHenryLawConstant(HenrysLawConstant(HenrysLawConstantParameters{ .HLC_ref_ = K_H }))
                           .Build();

  // OH-: water autodissociation equilibrium
  //   H2O(reactant) <-> H+ + OH-,  K_miam = Kw_miam = Kw_lit / c_H2O^2
  //   (H2O as explicit reactant AND solvent gives the correct Kw_miam)
  auto kw_constraint = DissolvedEquilibriumConstraintBuilder()
                           .SetPhase(sys.aqueous_phase)
                           .SetReactants({ sys.H2O })
                           .SetProducts({ sys.Hp, sys.OHm })
                           .SetAlgebraicSpecies(sys.OHm)
                           .SetSolvent(sys.H2O)
                           .SetEquilibriumConstant(EquilibriumConstant(EquilibriumConstantParameters{ .A_ = Kw_miam }))
                           .Build();

  // H+: charge balance  [H+] - [OH-] - [HCO3-] - 2[CO3--] = 0
  auto charge_balance = LinearConstraintBuilder()
                            .SetAlgebraicSpecies(sys.aqueous_phase, sys.Hp)
                            .AddTerm(sys.aqueous_phase, sys.Hp, 1.0)
                            .AddTerm(sys.aqueous_phase, sys.OHm, -1.0)
                            .AddTerm(sys.aqueous_phase, sys.HCO3m, -1.0)
                            .AddTerm(sys.aqueous_phase, sys.CO3mm, -2.0)
                            .SetConstant(0.0)
                            .Build();

  // -- Model --
  auto model = Model{ .name_ = "CARBONIC_ACID_DAE", .representations_ = { sys.droplet } };
  model.AddProcesses({ sys.rxn1 });
  model.AddConstraints(hl_constraint, k2_constraint, kw_constraint, charge_balance);

  // -- DAE Solver --
  auto system = System(sys.gas_phase);
  auto params = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  params.h_start_ = 1.0e-4;               // small initial step; Rosenbrock will grow it adaptively
  params.h_max_ = 1.0;                    // cap step size: k1_r*h_max ~1.3e7, well-conditioned in debug+release
  params.max_number_of_steps_ = 1000000;  // allow enough internal steps to cover 3000 s
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(params)
                    .SetSystem(system)
                    .AddExternalModel(model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  // -- Variable indices --
  std::size_t i_CO2g = state.variable_map_.at("CO2_g");
  std::size_t i_CO2aq = state.variable_map_.at("DROPLET.AQUEOUS.CO2_aq");
  std::size_t i_HCO3 = state.variable_map_.at("DROPLET.AQUEOUS.HCO3-");
  std::size_t i_CO3 = state.variable_map_.at("DROPLET.AQUEOUS.CO3--");
  std::size_t i_H = state.variable_map_.at("DROPLET.AQUEOUS.H+");
  std::size_t i_OH = state.variable_map_.at("DROPLET.AQUEOUS.OH-");
  std::size_t i_H2O = state.variable_map_.at("DROPLET.AQUEOUS.H2O");

  // -- Tight absolute tolerances for trace aqueous species --
  // The default atol (1e-3 mol m-3 air) >> equilibrium values of H+ (~7e-10),
  // OH- (~1e-12), CO3-- (~1e-14), CO2_aq (~4e-9), and HCO3- (~7e-10).
  // Without tighter values, the error norm accepts steps that drive these
  // species negative.
  state.absolute_tolerance_[i_H] = 1e-12;
  state.absolute_tolerance_[i_OH] = 1e-14;
  state.absolute_tolerance_[i_CO3] = 1e-16;
  state.absolute_tolerance_[i_CO2aq] = 1e-11;
  state.absolute_tolerance_[i_HCO3] = 1e-12;

  // -- Initial conditions (consistent with algebraic constraints at HCO3- = CO3-- = 0) --
  // Initialize algebraic variables to pure-water equilibrium so that all
  // constraint residuals vanish at t = 0, giving the Newton projection a
  // reliable starting point.
  const double f_v_0 = C_H2O * Mw_water / rho_water;
  const double CO2aq_0 = K_H * GAS_CONSTANT * T * f_v_0 * CO2g_eq;  // Henry's Law
  const double H_0 = std::sqrt(Kw_miam) * C_H2O;                    // pure-water [H+]

  state.variables_[0][i_CO2g] = CO2g_eq;   // 400 ppm
  state.variables_[0][i_CO2aq] = CO2aq_0;  // algebraic -- consistent with HL
  state.variables_[0][i_HCO3] = 0.0;
  state.variables_[0][i_CO3] = 0.0;
  state.variables_[0][i_H] = H_0;   // algebraic -- consistent with Kw + CB
  state.variables_[0][i_OH] = H_0;  // algebraic -- consistent with Kw
  state.variables_[0][i_H2O] = C_H2O;

  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = P;
  sys.droplet.SetDefaultParameters(state);

  std::cout << "\n=== DAEConstraints ===\n";
  PrintConditions(
      "Initial",
      state.variables_[0][i_CO2g],
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);

  // -- Integrate to equilibrium --
  solver.UpdateStateParameters(state);
  constexpr double t_end = 3000.0;
  auto result = solver.Solve(t_end, state);
  ASSERT_EQ(result.state_, SolverState::Converged) << "DAE solver failed";
  std::cout << "Solver steps: " << result.stats_.number_of_steps_ << "\n";

  PrintConditions(
      "Final",
      state.variables_[0][i_CO2g],
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);
  PrintEquilibriumStatus(
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);

  // -- Verify general equilibrium --

  // Non-negativity: all aqueous species must be >= 0 after solve
  EXPECT_GE(state.variables_[0][i_CO2aq], 0.0) << "CO2_aq went negative";
  EXPECT_GE(state.variables_[0][i_HCO3], 0.0) << "HCO3- went negative";
  EXPECT_GE(state.variables_[0][i_CO3], 0.0) << "CO3-- went negative";
  EXPECT_GE(state.variables_[0][i_H], 0.0) << "H+ went negative";
  EXPECT_GE(state.variables_[0][i_OH], 0.0) << "OH- went negative";
  EXPECT_GE(state.variables_[0][i_H2O], 0.0) << "H2O went negative";

  // Assertion guards: ensure ratio checks inside CheckEquilibrium actually run
  ASSERT_GT(state.variables_[0][i_HCO3], 0.0) << "HCO3- is zero — K1/K2 ratio checks would be silently skipped";
  ASSERT_GT(state.variables_[0][i_CO3], 0.0) << "CO3-- is zero — K2 ratio check would be silently skipped";
  ASSERT_GT(state.variables_[0][i_H], 0.0) << "H+ is zero — equilibrium ratio checks would be silently skipped";
  ASSERT_GT(state.variables_[0][i_OH], 0.0) << "OH- is zero — Kw ratio check would be silently skipped";

  // pH check: 400 ppm CO2 should give pH ≈ 5.6, i.e. [H+] ≈ 7.4e-10 mol m-3 air
  EXPECT_NEAR(state.variables_[0][i_H] / 7.4e-10, 1.0, 0.1) << "[H+] not near expected 7.4e-10 mol m-3 air (pH 5.6)";

  CheckEquilibrium(
      state.variables_[0][i_CO2aq],
      state.variables_[0][i_HCO3],
      state.variables_[0][i_CO3],
      state.variables_[0][i_H],
      state.variables_[0][i_OH],
      state.variables_[0][i_H2O]);

  // -- Verify algebraic constraints are satisfied to high precision --
  // The DAE solver enforces these exactly at each step, so residuals should be
  // negligible compared to the typical species concentration.
  double H_val = state.variables_[0][i_H];
  double OH_val = state.variables_[0][i_OH];
  double HCO3_val = state.variables_[0][i_HCO3];
  double CO3_val = state.variables_[0][i_CO3];
  double H2O_val = state.variables_[0][i_H2O];
  double CO2aq_val = state.variables_[0][i_CO2aq];
  double CO2g_val = state.variables_[0][i_CO2g];

  // Henry's Law constraint residual: K_H * R * T * f_v * [CO2_g] - [CO2_aq] = 0
  double f_v = H2O_val * Mw_water / rho_water;
  double hl_expected = K_H * GAS_CONSTANT * T * f_v * CO2g_val;
  if (hl_expected > 1.0e-30)
    EXPECT_NEAR(CO2aq_val / hl_expected, 1.0, 1.0e-4) << "Henry's-law constraint residual too large";

  // Kw constraint residual: Kw_miam * [H2O]^2 - [H+][OH-] = 0
  double kw_res = Kw_miam * H2O_val * H2O_val - H_val * OH_val;
  double kw_scale = H_val * OH_val;
  if (kw_scale > 1.0e-60)
    EXPECT_NEAR(kw_res / kw_scale, 0.0, 1.0e-4) << "Kw constraint residual too large";

  // Charge balance residual
  double cb = H_val - OH_val - HCO3_val - 2.0 * CO3_val;
  if (H_val > 1.0e-30)
    EXPECT_NEAR(cb / H_val, 0.0, 1.0e-4) << "Charge balance residual too large";
}
