// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Test that the README example compiles, runs, and produces expected output.

#include <miam/miam.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>

#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace micm;
using namespace miam;

TEST(ReadmeExample, HenryLawPhaseTransfer)
{
  // Define species with physical properties required for mass transfer
  auto co2_g  = Species{ "CO2_g",
      { { "molecular weight [kg mol-1]", 0.044 } } };
  auto co2_aq = Species{ "CO2_aq",
      { { "molecular weight [kg mol-1]", 0.044 },
        { "density [kg m-3]", 1800.0 } } };
  auto h2o    = Species{ "H2O",
      { { "molecular weight [kg mol-1]", 0.018 },
        { "density [kg m-3]", 1000.0 } } };

  // Define phases
  Phase gas_phase{ "GAS", { { co2_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { co2_aq }, { h2o } } };

  // Cloud droplets with a single-moment log-normal distribution
  auto cloud = representation::SingleMomentMode{
    "CLOUD",
    { aqueous_phase },
    5.0e-6,  // geometric mean radius [m]
    1.2      // geometric standard deviation
  };

  // Henry's Law phase transfer: CO2(g) <-> CO2(aq)
  auto co2_transfer = process::HenryLawPhaseTransferBuilder()
    .SetCondensedPhase(aqueous_phase)
    .SetGasSpecies(co2_g)
    .SetGasSpeciesName("CO2_g")
    .SetCondensedSpecies(co2_aq)
    .SetSolvent(h2o)
    .SetHenrysLawConstant(process::constant::HenrysLawConstant(
        { .HLC_ref_ = 3.4e-2 }))       // mol m-3 Pa-1 at 298 K
    .SetDiffusionCoefficient(1.5e-5)    // m2 s-1
    .SetAccommodationCoefficient(5.0e-6)
    .Build();

  // Create model and add processes
  auto cloud_model = Model{
    .name_ = "CLOUD",
    .representations_ = { cloud }
  };
  cloud_model.AddProcesses({ co2_transfer });

  // Build solver (MICM Rosenbrock)
  auto system = System(gas_phase, cloud_model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModelProcesses(cloud_model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  // Initial conditions
  state.conditions_[0].temperature_ = 298.15;   // K
  state.conditions_[0].pressure_ = 101325.0;    // Pa
  state.conditions_[0].CalculateIdealAirDensity();

  state[co2_g] = 1.0e-3;                               // mol m-3 air
  state[cloud.Species(aqueous_phase, h2o)] = 300.0;      // mol m-3 (liquid water content)
  cloud.SetDefaultParameters(state);

  // Capture output
  std::stringstream buffer;
  auto* old_buf = std::cout.rdbuf(buffer.rdbuf());

  state.PrintHeader();
  state.PrintState(0);
  for (int i = 1; i <= 10; ++i)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(0.1, state);  // 100 ms steps
    state.PrintState(i * 100);
  }

  std::cout.rdbuf(old_buf);

  std::string output = buffer.str();

  // Print the output so it's visible in test logs
  std::cout << output;

  // Verify the output is not empty and contains expected header columns
  EXPECT_FALSE(output.empty());
  EXPECT_NE(output.find("time"), std::string::npos);
  EXPECT_NE(output.find("CO2_g"), std::string::npos);

  // Verify physical constraints after integration:
  // 1. All concentrations should be non-negative
  for (std::size_t var = 0; var < state.variables_.NumColumns(); ++var)
  {
    EXPECT_GE(state.variables_[0][var], 0.0)
        << "Negative concentration for variable " << var;
  }

  // 2. Gas CO2 should decrease (dissolving into aqueous phase)
  double co2_g_final = state.variables_[0][state.variable_map_["CO2_g"]];
  EXPECT_LT(co2_g_final, 1.0e-3) << "Gas CO2 should decrease from initial value";

  // 3. Aqueous CO2 should increase from zero
  double co2_aq_final = state.variables_[0][state.variable_map_["CLOUD.AQUEOUS.CO2_aq"]];
  EXPECT_GT(co2_aq_final, 0.0) << "Aqueous CO2 should increase from zero";

  // 4. Mass conservation: gas + aqueous should equal initial total
  EXPECT_NEAR(co2_g_final + co2_aq_final, 1.0e-3, 1.0e-6)
      << "Mass conservation violated";

  // 5. Water (solvent) should be essentially unchanged
  double h2o_final = state.variables_[0][state.variable_map_["CLOUD.AQUEOUS.H2O"]];
  EXPECT_NEAR(h2o_final, 300.0, 1.0e-3) << "Solvent should be approximately conserved";
}
