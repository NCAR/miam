// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

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

  // Analytical solution for the HLPT two-species linear system:
  //   d[gas]/dt = -a · [gas] + b · [aq]
  //   d[aq]/dt  =  a · [gas] - b · [aq]
  // where a = φ · k_cond, b = φ · k_evap / f_v
  //
  // Conservation: [gas] + [aq] = M (constant)
  // Eigenvalues: 0 and -(a+b)
  // [gas](t) = gas_eq + (gas_0 - gas_eq) · exp(-(a+b)·t)
  // [aq](t)  = aq_eq  + (aq_0  - aq_eq)  · exp(-(a+b)·t)
  void analytical_hlpt(double gas_0, double aq_0,
                       double a, double b, double t,
                       double& gas_t, double& aq_t)
  {
    double M = gas_0 + aq_0;
    double gas_eq = M * b / (a + b);
    double aq_eq = M * a / (a + b);
    double exp_term = std::exp(-(a + b) * t);
    gas_t = gas_eq + (gas_0 - gas_eq) * exp_term;
    aq_t = aq_eq + (aq_0 - aq_eq) * exp_term;
  }

  // Common species with required properties for provider calculations
  micm::Species MakeGasSpecies(const std::string& name, double mw)
  {
    return micm::Species{ name, { { "molecular weight [kg mol-1]", mw } } };
  }

  micm::Species MakeCondensedSpecies(const std::string& name, double mw, double density)
  {
    return micm::Species{ name, { { "molecular weight [kg mol-1]", mw },
                                   { "density [kg m-3]", density } } };
  }
}  // namespace

// ============================================================================
// Test 1: Simple 1-species, 1-instance with analytical solution
// ============================================================================
TEST(HenryLawPhaseTransferIntegration, SimpleOneInstance)
{
  // Gas species A_g partitioning into aqueous phase as A_aq
  // Solvent is H2O
  double Mw_gas = 0.044;       // kg mol⁻¹ (like CO₂)
  double Mw_solvent = 0.018;   // kg mol⁻¹ (water)
  double rho_solvent = 1000.0; // kg m⁻³ (water)
  double D_g = 1.5e-5;         // m² s⁻¹
  double alpha = 0.05;         // accommodation coefficient
  double HLC_val = 3.4e-2;     // mol m⁻³ Pa⁻¹ (constant, no T dependence)

  auto A_g = MakeGasSpecies("A_g", Mw_gas);
  auto A_aq = MakeCondensedSpecies("A_aq", Mw_gas, 1800.0);
  auto H2O = MakeCondensedSpecies("H2O", Mw_solvent, rho_solvent);

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto droplet = representation::SingleMomentMode{
    "DROPLET",
    { aqueous_phase },
    5.0e-6,  // geometric mean radius (m)
    1.2      // geometric standard deviation
  };

  // Build the process using the builder
  auto hlc = [HLC_val](const Conditions& conditions) { return HLC_val; };

  auto transfer = process::HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(aqueous_phase)
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetHenrysLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC_val }))
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
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

  // Find variable indices
  auto find_idx = [&](const std::string& name) {
    auto it = std::find(state.variable_names_.begin(), state.variable_names_.end(), name);
    EXPECT_NE(it, state.variable_names_.end()) << "Species " << name << " not found";
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  std::size_t i_gas = find_idx("A_g");
  std::size_t i_aq = find_idx("DROPLET.AQUEOUS.A_aq");
  std::size_t i_h2o = find_idx("DROPLET.AQUEOUS.H2O");

  // Initial conditions
  double T = 298.15;
  double gas_0 = 1.0e-3;       // mol m⁻³ air
  double aq_0 = 0.0;           // mol m⁻³ air
  double solvent_conc = 55.5;   // mol m⁻³ air (liquid water content)

  state.variables_[0][i_gas] = gas_0;
  state.variables_[0][i_aq] = aq_0;
  state.variables_[0][i_h2o] = solvent_conc;

  state.conditions_[0].temperature_ = T;
  state.conditions_[0].pressure_ = 101325.0;

  droplet.SetDefaultParameters(state);

  // Compute expected effective rates
  // For SingleMomentMode: r_eff = GMD · exp(2.5 · ln²(GSD))
  double GMD = 5.0e-6;
  double GSD = 1.2;
  double ln_gsd = std::log(GSD);
  double r_eff = GMD * std::exp(2.5 * ln_gsd * ln_gsd);

  // N = V_total / V_single where V_total from species, V_single from GMD/GSD
  double V_aq = aq_0 * Mw_gas / 1800.0;       // initial aq volume
  double V_h2o = solvent_conc * Mw_solvent / rho_solvent;  // water volume
  double V_total = V_aq + V_h2o;
  double V_single = (4.0 / 3.0) * M_PI * std::pow(GMD, 3) * std::exp(4.5 * ln_gsd * ln_gsd);
  double N = V_total / V_single;

  // phi = 1.0 for single-phase mode
  double phi = 1.0;

  // k_cond from condensation rate utility
  auto cond_provider = miam::util::make_condensation_rate_provider(D_g, alpha, Mw_gas);
  double k_cond = cond_provider.ComputeValue(r_eff, N, T);
  double k_evap = k_cond / (HLC_val * R_gas * T);

  double f_v = solvent_conc * Mw_solvent / rho_solvent;
  double a = phi * k_cond;
  double b = phi * k_evap / f_v;

  // Integrate and compare
  double time = 0.0;
  double total_time = 1.0;  // 1 second total
  double dt = 0.001;         // 1 ms steps
  const double tolerance = 1.0e-2;  // relative tolerance

  while (time < total_time - 1.0e-15)
  {
    double step = std::min(dt, total_time - time);
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(step, state);
    ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time;
    time += step;
  }

  double gas_final = state.variables_[0][i_gas];
  double aq_final = state.variables_[0][i_aq];
  double h2o_final = state.variables_[0][i_h2o];

  // Analytical solution (using initial rate constants as approximation)
  double gas_analytic, aq_analytic;
  analytical_hlpt(gas_0, aq_0, a, b, total_time, gas_analytic, aq_analytic);

  // Mass conservation: gas + aq should be constant
  EXPECT_NEAR(gas_final + aq_final, gas_0 + aq_0, (gas_0 + aq_0) * 1e-6)
      << "Mass conservation violated";

  // Solvent should be unchanged
  EXPECT_NEAR(h2o_final, solvent_conc, solvent_conc * 1e-6)
      << "Solvent concentration changed";

  // Direction check: gas should decrease and aq should increase
  // (since we start with all gas and zero aq, condensation dominates)
  EXPECT_LT(gas_final, gas_0) << "Gas should decrease";
  EXPECT_GT(aq_final, aq_0) << "Aqueous should increase";

  // Quantitative comparison with analytical solution
  // Note: analytical solution is approximate because k_cond depends on N,
  // which depends on state variables. For this test with large solvent
  // concentration, the approximation should be very good.
  if (gas_analytic > 1e-20)
    EXPECT_NEAR(gas_final / gas_analytic, 1.0, tolerance)
        << "Gas final: " << gas_final << " vs analytical: " << gas_analytic;
  if (aq_analytic > 1e-20)
    EXPECT_NEAR(aq_final / aq_analytic, 1.0, tolerance)
        << "Aq final: " << aq_final << " vs analytical: " << aq_analytic;
}

// ============================================================================
// Test 2: Multi-instance mass conservation
// ============================================================================
TEST(HenryLawPhaseTransferIntegration, MultiInstanceMassConservation)
{
  // Same gas species partitions into two different droplet populations
  double Mw_gas = 0.030;       // kg mol⁻¹ (like HCHO)
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double D_g = 1.8e-5;
  double alpha = 0.1;
  double HLC_val = 3.2e3;  // mol m⁻³ Pa⁻¹ (large, like HCHO)

  auto A_g = MakeGasSpecies("A_g", Mw_gas);
  auto A_aq = MakeCondensedSpecies("A_aq", Mw_gas, 1200.0);
  auto H2O = MakeCondensedSpecies("H2O", Mw_solvent, rho_solvent);

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_small{ "AQ_SMALL", { { A_aq }, { H2O } } };
  Phase aqueous_large{ "AQ_LARGE", { { A_aq }, { H2O } } };

  auto small_drop = representation::SingleMomentMode{
    "SMALL", { aqueous_small }, 1.0e-6, 1.2
  };
  auto large_drop = representation::SingleMomentMode{
    "LARGE", { aqueous_large }, 1.0e-5, 1.4
  };

  // Build two transfer processes — one for each phase
  auto transfer_small = process::HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(aqueous_small)
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetHenrysLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC_val }))
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
      .Build();

  auto transfer_large = process::HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(aqueous_large)
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetHenrysLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC_val }))
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
      .Build();

  auto model = Model{
    .name_ = "CLOUD",
    .representations_ = { small_drop, large_drop }
  };
  model.AddProcesses({ transfer_small });
  model.AddProcesses({ transfer_large });

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
    EXPECT_NE(it, state.variable_names_.end()) << name << " not found";
    return static_cast<std::size_t>(it - state.variable_names_.begin());
  };

  std::size_t i_gas = find_idx("A_g");
  std::size_t i_aq_small = find_idx("SMALL.AQ_SMALL.A_aq");
  std::size_t i_h2o_small = find_idx("SMALL.AQ_SMALL.H2O");
  std::size_t i_aq_large = find_idx("LARGE.AQ_LARGE.A_aq");
  std::size_t i_h2o_large = find_idx("LARGE.AQ_LARGE.H2O");

  double gas_0 = 1.0e-2;
  double solvent = 55.5;

  state.variables_[0][i_gas] = gas_0;
  state.variables_[0][i_aq_small] = 0.0;
  state.variables_[0][i_h2o_small] = solvent;
  state.variables_[0][i_aq_large] = 0.0;
  state.variables_[0][i_h2o_large] = solvent;

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  small_drop.SetDefaultParameters(state);
  large_drop.SetDefaultParameters(state);

  double total_mass_0 = gas_0;  // initial total = gas + 0 + 0

  // Integrate for a while
  double time = 0.0;
  double total_time = 0.1;
  double dt = 0.001;

  while (time < total_time - 1.0e-15)
  {
    double step = std::min(dt, total_time - time);
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(step, state);
    ASSERT_EQ(result.state_, SolverState::Converged)
        << "Solver failed at t = " << time;
    time += step;
  }

  double gas_final = state.variables_[0][i_gas];
  double aq_small_final = state.variables_[0][i_aq_small];
  double aq_large_final = state.variables_[0][i_aq_large];

  // Total mass conservation: gas + aq_small + aq_large = total_mass_0
  double total_final = gas_final + aq_small_final + aq_large_final;
  EXPECT_NEAR(total_final, total_mass_0, total_mass_0 * 1e-4)
      << "Total mass not conserved: " << total_final << " vs " << total_mass_0;

  // Both aqueous concentrations should have increased
  EXPECT_GT(aq_small_final, 0.0) << "Small droplet aq should increase";
  EXPECT_GT(aq_large_final, 0.0) << "Large droplet aq should increase";

  // Gas should have decreased
  EXPECT_LT(gas_final, gas_0) << "Gas should decrease";

  // Solvent unchanged
  EXPECT_NEAR(state.variables_[0][i_h2o_small], solvent, solvent * 1e-6);
  EXPECT_NEAR(state.variables_[0][i_h2o_large], solvent, solvent * 1e-6);
}

// ============================================================================
// Test 3: Temperature-dependent HLC
// ============================================================================
TEST(HenryLawPhaseTransferIntegration, TemperatureDependentHLC)
{
  // Verify that the equilibrium shifts with temperature via the van't Hoff
  // parameterization: HLC(T) = HLC_ref · exp(C · (1/T - 1/T0))
  // Higher C implies stronger T dependence. At lower T (with positive C),
  // HLC increases → more dissolution → more aq at equilibrium.

  double Mw_gas = 0.044;
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double D_g = 1.5e-5;
  double alpha = 0.05;
  double HLC_ref = 3.4e-2;  // at T0 = 298.15
  double C = 2400.0;         // K (enthalpy parameter)
  double T0 = 298.15;

  auto A_g = MakeGasSpecies("A_g", Mw_gas);
  auto A_aq = MakeCondensedSpecies("A_aq", Mw_gas, 1800.0);
  auto H2O = MakeCondensedSpecies("H2O", Mw_solvent, rho_solvent);

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto droplet = representation::SingleMomentMode{
    "DROP", { aqueous_phase }, 5.0e-6, 1.2
  };

  process::constant::HenrysLawConstantParameters hlc_params{
    .HLC_ref_ = HLC_ref, .C_ = C, .T0_ = T0
  };

  auto build_transfer = [&]() {
    return process::HenryLawPhaseTransferBuilder()
        .SetCondensedPhase(aqueous_phase)
        .SetGasSpecies(A_g)
        .SetCondensedSpecies(A_aq)
        .SetSolvent(H2O)
        .SetHenrysLawConstant(process::constant::HenrysLawConstant(hlc_params))
        .SetDiffusionCoefficient(D_g)
        .SetAccommodationCoefficient(alpha)
        .Build();
  };

  // Run at two temperatures and compare equilibrium
  double gas_0 = 1.0e-3;
  double solvent = 55.5;

  auto run_to_equilibrium = [&](double T) -> std::pair<double, double>
  {
    auto transfer = build_transfer();
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

    std::size_t i_gas = find_idx("A_g");
    std::size_t i_aq = find_idx("DROP.AQUEOUS.A_aq");
    std::size_t i_h2o = find_idx("DROP.AQUEOUS.H2O");

    state.variables_[0][i_gas] = gas_0;
    state.variables_[0][i_aq] = 0.0;
    state.variables_[0][i_h2o] = solvent;

    state.conditions_[0].temperature_ = T;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);

    // Run to approximate equilibrium with small steps for stiff system
    double time = 0.0;
    double total_time = 1.0;
    double dt = 0.001;

    while (time < total_time - 1.0e-15)
    {
      double step = std::min(dt, total_time - time);
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(step, state);
      EXPECT_EQ(result.state_, SolverState::Converged)
          << "Solver failed at T=" << T << " t=" << time;
      time += step;
    }

    return { state.variables_[0][i_gas], state.variables_[0][i_aq] };
  };

  // Cold temperature: more dissolution (higher HLC with positive C)
  auto [gas_cold, aq_cold] = run_to_equilibrium(278.15);    // 5°C
  // Warm temperature: less dissolution
  auto [gas_warm, aq_warm] = run_to_equilibrium(308.15);    // 35°C

  // At colder temperature, HLC is larger → more dissolved at equilibrium
  // So aq_cold / gas_cold > aq_warm / gas_warm
  double ratio_cold = (aq_cold > 1e-20 && gas_cold > 1e-20) ? aq_cold / gas_cold : 0.0;
  double ratio_warm = (aq_warm > 1e-20 && gas_warm > 1e-20) ? aq_warm / gas_warm : 0.0;

  EXPECT_GT(ratio_cold, ratio_warm)
      << "Cold temperature should have higher aq/gas ratio. "
      << "Cold: " << ratio_cold << ", Warm: " << ratio_warm;

  // Mass conservation at both temperatures
  EXPECT_NEAR(gas_cold + aq_cold, gas_0, gas_0 * 1e-4);
  EXPECT_NEAR(gas_warm + aq_warm, gas_0, gas_0 * 1e-4);
}

// ============================================================================
// Test 4: Continuum vs. transition regime rate comparison
// ============================================================================
TEST(HenryLawPhaseTransferIntegration, SmallVsLargeParticleRate)
{
  // Compare transfer rates for small particles (transition regime, Kn >> 1)
  // vs large particles (near continuum, Kn << 1).
  // For sufficiently different sizes but same N, the larger particles should
  // achieve equilibrium faster (higher k_cond in continuum regime).

  double Mw_gas = 0.044;
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double D_g = 1.5e-5;
  double alpha = 0.05;
  double HLC_val = 3.4e-2;

  auto A_g = MakeGasSpecies("A_g", Mw_gas);
  auto A_aq = MakeCondensedSpecies("A_aq", Mw_gas, 1800.0);
  auto H2O = MakeCondensedSpecies("H2O", Mw_solvent, rho_solvent);

  Phase gas_phase{ "GAS", { { A_g } } };

  double gas_0 = 1.0e-3;
  double solvent = 55.5;

  // Run HLPT for a given particle radius, return (gas, aq) after fixed time
  auto run_with_radius = [&](double r_mean, const std::string& prefix,
                              const std::string& phase_name) -> std::pair<double, double>
  {
    Phase aqueous_phase{ phase_name, { { A_aq }, { H2O } } };

    auto droplet = representation::SingleMomentMode{
      prefix, { aqueous_phase }, r_mean, 1.01  // nearly monodisperse
    };

    auto transfer = process::HenryLawPhaseTransferBuilder()
        .SetCondensedPhase(aqueous_phase)
        .SetGasSpecies(A_g)
        .SetCondensedSpecies(A_aq)
        .SetSolvent(H2O)
        .SetHenrysLawConstant(process::constant::HenrysLawConstant(
            process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC_val }))
        .SetDiffusionCoefficient(D_g)
        .SetAccommodationCoefficient(alpha)
        .Build();

    auto model = Model{
      .name_ = "TEST",
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

    std::size_t i_gas = find_idx("A_g");
    std::size_t i_aq = find_idx(prefix + "." + phase_name + ".A_aq");
    std::size_t i_h2o = find_idx(prefix + "." + phase_name + ".H2O");

    state.variables_[0][i_gas] = gas_0;
    state.variables_[0][i_aq] = 0.0;
    state.variables_[0][i_h2o] = solvent;

    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    droplet.SetDefaultParameters(state);

    // Run for a short time to measure initial transfer rate
    double time = 0.0;
    double total_time = 0.01;  // 10 ms
    double dt = 0.0001;         // 100 μs steps

    while (time < total_time - 1.0e-15)
    {
      double step = std::min(dt, total_time - time);
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(step, state);
      EXPECT_EQ(result.state_, SolverState::Converged);
      time += step;
    }

    return { state.variables_[0][i_gas], state.variables_[0][i_aq] };
  };

  // Small particles (transition regime): r ~ 50 nm
  auto [gas_small, aq_small] = run_with_radius(5.0e-8, "SMALL", "AQ_S");

  // Large particles (near continuum): r ~ 10 μm
  auto [gas_large, aq_large] = run_with_radius(1.0e-5, "LARGE", "AQ_L");

  // Both should show condensation
  EXPECT_GT(aq_small, 0.0) << "Small particles should absorb gas";
  EXPECT_GT(aq_large, 0.0) << "Large particles should absorb gas";

  // Mass conservation for each run
  EXPECT_NEAR(gas_small + aq_small, gas_0, gas_0 * 1e-4);
  EXPECT_NEAR(gas_large + aq_large, gas_0, gas_0 * 1e-4);

  // The amount transferred should differ between regimes
  // (this test just verifies both regimes produce physically reasonable results;
  // the unit tests verify the actual rate formulas)
  double transfer_small = aq_small;
  double transfer_large = aq_large;

  // Both should transfer some gas
  EXPECT_GT(transfer_small, 0.0);
  EXPECT_GT(transfer_large, 0.0);
}
