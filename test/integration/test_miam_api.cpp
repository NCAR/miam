#include <miam/miam.hpp>
#include <micm/CPU.hpp>

#include <iomanip>
#include <iostream>
#include <cstdlib>

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
  auto h2co3  = Species{ "H2CO3" };
  auto hexane = Species{ "C6H14" };

  Phase gas_phase{ "GAS", { { co2, 31.2 } } };
  Phase aqueous_phase{ "AQUEOUS", { { co2 } , { h2o } , { ohm } , { hp } , { hco3m } , { co32m } } };
  Phase organic_phase{ "ORGANIC", { { co2, 16.2 }, { hexane } } };

  // Cloud droplets modeled with 1-moment log-normal distributions
  auto small_drop = representation::SingleMomentMode{
    "SMALL_DROP",
    { aqueous_phase },
    1.0e-7, // geometric mean radius (m)
    1.1     // geometric standard deviation
  };
  
  // Larger cloud droplets modeled with 1-moment log-normal distributions
  auto large_drop = representation::SingleMomentMode{
    "LARGE_DROP",
    { aqueous_phase },
    1.0e-6, // geometric mean radius (m)
    1.4     // geometric standard deviation
  };

  // Aerosol model with 2-moment distribution
  auto aitken = representation::TwoMomentMode{
    "AITKEN",
    { aqueous_phase },
    1.2     // geometric standard deviation
  };
  
  // Another aerosol mode with 2-moment distribution
  auto accumulation = representation::TwoMomentMode{
    "ACCUMULATION",
    { aqueous_phase, organic_phase },  // Multiple phases
    1.4     // geometric standard deviation
  };

  // Dust particle section
  auto dust = representation::UniformSection{
    "DUST",
    { organic_phase },
    1.0e-7, // min radius (m)
    1.0e-6  // max radius (m)
  };

  auto aerosol_model = Model{
    .name_ = "AEROSOL",
    .representations_ = { aitken, accumulation, dust }
  };

  auto cloud_model = Model{
    .name_ = "CLOUD",
    .representations_ = { small_drop, large_drop }
  };

  auto system = System(gas_phase, aerosol_model, cloud_model);

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
  // 15) (AEROSOL.) AITKEN.AQUEOUS.H2O
  // 16) (AEROSOL.) AITKEN.AQUEOUS.OH-
  // 17) (AEROSOL.) AITKEN.AQUEOUS.H+
  // 18) (AEROSOL.) AITKEN.AQUEOUS.HCO3-
  // 19) (AEROSOL.) AITKEN.AQUEOUS.CO32-
  // 20) (AEROSOL.) AITKEN.number_concentration
  // 21) (AEROSOL.) ACCUMULATION.AQUEOUS.CO2
  // 22) (AEROSOL.) ACCUMULATION.AQUEOUS.H2O
  // 23) (AEROSOL.) ACCUMULATION.AQUEOUS.OH-
  // 24) (AEROSOL.) ACCUMULATION.AQUEOUS.H+
  // 25) (AEROSOL.) ACCUMULATION.AQUEOUS.HCO3-
  // 26) (AEROSOL.) ACCUMULATION.AQUEOUS.CO32-
  // 27) (AEROSOL.) ACCUMULATION.ORGANIC.CO2
  // 28) (AEROSOL.) ACCUMULATION.ORGANIC.CH6H14
  // 29) (AEROSOL.) ACCUMULATION.number_concentration
  // 30) (AEROSOL.) DUST.ORGANIC.CO2
  // 31) (AEROSOL.) DUST.ORGANIC.CH6H14
  // 32) (AEROSOL.) DUST.number_concentration

#if 0
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
                               .SetTransferCoefficient(HenrysLawCoefficient(hlc_co2_parms))
                               .Build();
#endif

  // // Condensed phase reversible reaction
  // // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // // k_r = reverse rate constant
  // // (k_f = K_eq * k_r = forward rate constant)
  auto h2o_dissociation = process::DissolvedReversibleReactionBuilder()
                             .SetPhase(aqueous_phase)
                             .SetReactants({ h2o })
                             .SetProducts({ ohm, hp })
                             .SetSolvent(h2o)
                             .SetEquilibriumConstant(process::constant::EquilibriumConstant(process::constant::EquilibriumConstantParameters{ .A_ = 1.14e-2, .C_ = 2300.0, .T0_ = 298.15 }))
                             .SetReverseRateConstant(ArrheniusRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.4e11, .C_ = 5.1e4 }))
                             .Build();

  // Condensed phase reversible reaction: CO2 hydration
  // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  auto co2_hydration = process::DissolvedReversibleReactionBuilder()
                          .SetPhase(aqueous_phase)
                          .SetReactants({ co2, h2o })
                          .SetProducts({ h2co3 })
                          .SetSolvent(h2o)
                          .SetEquilibriumConstant(process::constant::EquilibriumConstant(process::constant::EquilibriumConstantParameters{ .A_ = 1.70e3, .C_ = 2400.0, .T0_ = 298.15 }))
                          .SetReverseRateConstant(ArrheniusRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.4e11, .C_ = 5.1e4 }))
                          .Build();

  // Condensed phase reversible reaction: H2CO3 dissociation
  // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  auto h2co3_dissociation = process::DissolvedReversibleReactionBuilder()
                               .SetPhase(aqueous_phase)
                               .SetReactants({ h2co3 })
                               .SetProducts({ hco3m, hp })
                               .SetSolvent(h2o)
                               .SetEquilibriumConstant(process::constant::EquilibriumConstant(process::constant::EquilibriumConstantParameters{ .A_ = 4.27e2, .C_ = 2300.0, .T0_ = 298.15 }))
                               .SetReverseRateConstant(ArrheniusRateConstant(ArrheniusRateConstantParameters{ .A_ = 2.5e10, .C_ = 4.0e4 }))
                               .Build();

  // Condensed phase reversible reaction: HCO3- dissociation
  // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  auto hco3m_dissociation = process::DissolvedReversibleReactionBuilder()
                               .SetPhase(aqueous_phase)
                               .SetReactants({ hco3m })
                               .SetProducts({ co32m, hp })
                               .SetSolvent(h2o)
                               .SetEquilibriumConstant(process::constant::EquilibriumConstant(process::constant::EquilibriumConstantParameters{ .A_ = 1.70e1, .C_ = 2300.0, .T0_ = 298.15 }))
                               .SetReverseRateConstant(ArrheniusRateConstant(ArrheniusRateConstantParameters{ .A_ = 6.4e9, .C_ = 3.1e4 }))
                               .Build();
  
  std::vector<process::DissolvedReversibleReaction> reactions{
//    co2_photo,
//    co2_phase_transfer,
    h2o_dissociation,
    co2_hydration,
    h2co3_dissociation,
    hco3m_dissociation
  };

  aerosol_model.AddProcesses(reactions);
  cloud_model.AddProcesses(reactions);

  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModelProcesses(aerosol_model)
                  .AddExternalModelProcesses(cloud_model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();
  
  State state = solver.GetState();
  
  // Start the state in some arbitrary initial condition
  auto& state_vec = state.variables_.AsVector();
  std::generate(state_vec.begin(), state_vec.end(), [&]() { 
    return 1.0e-6 + static_cast<double>(std::rand()) / RAND_MAX * 1.0e-6; 
  });

  // environmental conditions
  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.conditions_[0].CalculateIdealAirDensity();

  // gas
  state[co2] = 0.2;

  // cloud
  state[small_drop.Species(aqueous_phase, h2o)] = 0.3;      // mol m-3
  state[large_drop.Species(aqueous_phase, h2o)] = 0.3;      // mol m-3
  
  // aerosol mode
  state[aitken.Species(aqueous_phase, hco3m)] = 0.1;        // mol m-3
  state[accumulation.Species(aqueous_phase, h2o)] = 0.3;    // mol m-3
  state[accumulation.Species(organic_phase, hexane)] = 0.1; // mol m-3
  
  // aerosol section
  state[dust.Species(organic_phase, co2)] = 0.1;            // mol m-3

  // set aerosol state parameters
  state[aitken.NumberConcentration()] = 1.0e8; // m-3
  state[accumulation.NumberConcentration()] = 1.0e7; // m-3

  // Set aerosol parameters
  small_drop.SetDefaultParameters(state);
  large_drop.SetDefaultParameters(state);
  aitken.SetDefaultParameters(state);
  accumulation.SetDefaultParameters(state);
  dust.SetDefaultParameters(state);

  // state.PrintHeader();
  for (int i = 0; i < 10; ++i)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(500.0, state);
    // state.PrintState(i * 500);
  }

  // Print with species names if variable_names_ is available
  std::cout << "\nState Variables by Species:" << std::endl;
  for (std::size_t cell = 0; cell < state.variables_.NumRows(); ++cell) 
  {
    std::cout << "Cell " << cell << ":" << std::endl;
    for (std::size_t var = 0; var < state.variables_.NumColumns(); ++var) 
    {
      std::string species_name = (var < state.variable_names_.size()) 
                              ? state.variable_names_[var] 
                              : "Var" + std::to_string(var);
      std::cout << "  " << std::setw(12) << species_name << " = " 
                << std::scientific << state.variables_[cell][var] << std::endl;
    }
    std::cout << std::endl;
  }
}