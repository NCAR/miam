#include <miam/model/aerosol_model.hpp>
#include <miam/model/aerosol_scheme.hpp>
#include <miam/model/mode.hpp>
#include <miam/model/section.hpp>
#include <miam/model/gas_model.hpp>
#include <miam/util/solver_utils.hpp>

#include <micm/process/transfer_coefficient/phase_transfer_coefficient.hpp>
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/Process.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

#include <iomanip>
#include <iostream>

using namespace micm;
using namespace miam;

int main()
{
  auto co2    = Species{ "CO2" };
  auto h2o    = Species{ "H2O" };
  auto ohm    = Species{ "OH-" };
  auto hp     = Species{ "H+" };
  auto hco3m  = Species{ "HCO3-" };
  auto co32m  = Species{ "CO32-" };
  auto hexane = Species{ "C6H14" };

  // MICM
  Phase gas_phase{ "GAS", { { co2, 31.2 } } };
  Phase aqueous_phase{ "AQUEOUS", { { co2 } , { h2o } , { ohm } , { hp } , { hco3m } , { co32m } } };
  Phase organic_phase{ "ORGANIC", { { co2, 16.2 }, { hexane } } };

  // Gas 
  auto gas = GasModel{ gas_phase };

  // Cloud
  auto small_drop = Mode{
    "SMALL_DROP",
    { aqueous_phase },
    DistributionType::SingleMoment, // tracks total mass in state; fixed radius; number calculated
    5.0e-6,                         // Geometric mean diameter
    1.6,                            // Geometric standard deviation
  };
  
  auto large_drop = Mode{
    "LARGE_DROP",
    { aqueous_phase },
    DistributionType::SingleMoment, // tracks total mass in state; fixed radius; number calculated
    1.0e-5,                         // Geometric mean diameter
    1.8,                            // Geometric standard deviation
  };

  AerosolModel cloud{
    "CLOUD",
    { small_drop, large_drop },
    {} // no sectional model
  };

  // Aerosol model with 2 modes and 1 section
  auto aitken = Mode{
    "AITKEN",
    { aqueous_phase },           // Multiple phases
    DistributionType::TwoMoment, // tracks total mass and number concentration in state; radius calculated
    1.0e-7,                      // Geometric mean diameter
    1.6,                         // Geometric standard deviation
  };
  
  auto accumulation = Mode{
    "ACCUMULATION",
    { aqueous_phase, organic_phase }, // Multiple phases
    DistributionType::TwoMoment,      // tracks total mass and number concentration in state; radius calculated
    1.0e-6,                           // Geometric mean diameter
    1.6,                              // Geometric standard deviation
  };

  auto dust = Section{
    "DUST",
    { organic_phase }, 
    DistributionType::TwoMoment, // tracks total mass and number concentration in state; radius calculated
    1.0e-6,                      // Minimum diameter
    0.003                        // Maximum diameter
  };

  AerosolModel aerosol{
    "AEROSOL",
    { aitken, accumulation } ,
    { dust }
  };

  System chemical_system = ConfigureSystem(gas_phase, { cloud, aerosol });

  // State array should contain
  //  1) (GAS.) CO2
  //  2) (CLOUD.) SMALL_DROP.AQUEOUS.CO2
  //  3) (CLOUD.) SMALL_DROP.AQUEOUS.H2O
  //  4) (CLOUD.) SMALL_DROP.AQUEOUS.OH-
  //  5) (CLOUD.) SMALL_DROP.AQUEOUS.H+
  //  6) (CLOUD.) SMALL_DROP.AQUEOUS.HCO3-
  //  7) (CLOUD.) SMALL_DROP.AQUEOUS.CO32-
  //  8) (CLOUD.) LARGE_DROP.AQUEOUS.CO2
  //  9) (CLOUD.) LARGE_DROP.AQUEOUS.H2O
  // 10) (CLOUD.) LARGE_DROP.AQUEOUS.OH-
  // 11) (CLOUD.) LARGE_DROP.AQUEOUS.H+
  // 12) (CLOUD.) LARGE_DROP.AQUEOUS.HCO3-
  // 13) (CLOUD.) LARGE_DROP.AQUEOUS.CO32-
  // 14) (AEROSOL.) AITKEN.AQUEOUS.CO2
  // 15) (AEROSOL.) AITKEN.AQUEOUS.HEXANE
  // 16) (AEROSOL.) ACCUMULATION.AQUEOUS.CO2
  // 17) (AEROSOL.) ACCUMULATION.AQUEOUS.H2O
  // 18) (AEROSOL.) ACCUMULATION.AQUEOUS.OH-
  // 19) (AEROSOL.) ACCUMULATION.AQUEOUS.H+
  // 20) (AEROSOL.) ACCUMULATION.AQUEOUS.HCO3-
  // 21) (AEROSOL.) ACCUMULATION.AQUEOUS.CO32-
  // 22) (AEROSOL.) ACCUMULATION.ORGANIC.CO2
  // 23) (AEROSOL.) ACCUMULATION.ORGANIC.HEXANE
  // 24) (AEROSOL.) DUST.ORGANIC.CO2
  // 25) (AEROSOL.) DUST.ORGANIC.HEXANE

  Process co2_photo = ChemicalReactionBuilder()
                      .SetReactants({ co2 })
                      .SetRateConstant(ArrheniusRateConstant({ .A_ = 1.0e-3 }))
                      .SetPhase(gas_phase)
                      .Build();

  // Build an Effective Henry's Law PhaseTransferProcess
  // (will be for a diprotic acid: CO2(g) <-> H2CO3(aq) <-> HCO3- + H+ <-> CO32- + 2H+)
  // K_a1 = first acid dissociation constant
  // K_a2 = second acid dissociation constant
  // K_H = A * exp(C / T) = Henry's Law constant
  auto hlc_co2_parms = HenrysLawCoefficientParameters{ .A_ = 32.4, .C_ = -3.2e-3, .K_a1_ = 1.0e-5, .K_a2_ = 2.0e-5 };

  Process co2_phase_transfer = PhaseTransferProcessBuilder()
                               .SetGasSpecies(gas_phase, co2 )
                               .SetCondensedSpecies(aqueous_phase, { Yield(hp, 2.0), Yield(co32m) })
                               .SetSolvent(aqueous_phase, h2o )
                               .SetTransferCoefficient( HenrysLawCoefficient(hlc_co2_parms) )
                               .Build();

  // // Condensed phase reversible reaction
  // // K_eq = A * exp(C / T) = Equilibrium constant
  // // k_r = reverse rate constant
  // // (k_f = K_eq * k_r = forward rate constant)
  Process h2o_dissociation = ChemicalReactionBuilder()
                             .SetAerosolScope(accumulation.GetScope(), aqueous_phase)
                             .SetReactants({ h2o })
                             .SetProducts({ ohm, hp })
                             .SetRateConstant(ReversibleRateConstant({ .A_ = 1.14e-2, .C_ = 2300.0, .k_r_ = 0.32 }))
                             .Build();

  std::vector<Process> reactions{ co2_photo, co2_phase_transfer, h2o_dissociation };

  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(chemical_system)
                  .SetReactions(reactions)
                  .Build();
  
  State state = solver.GetState();
  
  // Initialize state indices map
  gas.SetStateIndices(state);
  small_drop.SetStateIndices(state);
  large_drop.SetStateIndices(state);
  aitken.SetStateIndices(state);
  accumulation.SetStateIndices(state);

  // environmental conditions
  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.conditions_[0].CalculateIdealAirDensity();

  // gas
  gas.SetConcentration(state, co2, 20.0);  // mol m-3
  
  // cloud
  small_drop.SetConcentration(state, aqueous_phase, h2o, 0.3); // mol m-3
  large_drop.SetConcentration(state, aqueous_phase, h2o, 0.3); // mol m-3
  
  // aerosol modes
  aitken.SetConcentration(state, aqueous_phase, hco3m, 0.1);        // mol m-3
  accumulation.SetConcentration(state, aqueous_phase, h2o, 0.3);    // mol m-3
  accumulation.SetConcentration(state, organic_phase, hexane, 0.1); // mol m-3

  // aerosol section 
  dust.SetConcentration(state, organic_phase, co2, 0.1); // mol m-3

  aitken.SetNumberConcentration(state, 1.0e4);       // m-3
  accumulation.SetNumberConcentration(state, 1.0e3); // m-3

  small_drop.SetRadius(state);   // m
  large_drop.SetRadius(state);   // m
  aitken.SetRadius(state);       // m
  accumulation.SetRadius(state); // m
  accumulation.SetRadius(state); // m

  state.PrintHeader();
  for (int i = 0; i < 10; ++i)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(500.0, state);
    state.PrintState(i * 500);
  }
}