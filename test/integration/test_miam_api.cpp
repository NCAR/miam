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
  auto hexane = Species{ "C6H14" };

  Phase gas_phase{ "GAS", { { co2, 31.2 } } };
  Phase aqueous_phase{ "AQUEOUS", { { co2 } , { h2o } , { ohm } , { hp } , { hco3m } , { co32m } } };
  Phase organic_phase{ "ORGANIC", { { co2, 16.2 }, { hexane } } };

  // Cloud droplets modeled with 1-moment log-normal distributions
  auto small_drop = Distribution<shape::LogNormal, moment::Single>{
    "SMALL_DROP",
    { aqueous_phase }
  };
  
  // Larger cloud droplets modeled with 1-moment log-normal distributions
  auto large_drop = Distribution<shape::LogNormal, moment::Single>{
    "LARGE_DROP",
    { aqueous_phase }
  };

  // Aerosol model with 2-moment distribution
  auto aitken = Distribution<shape::LogNormal, moment::Two>{
    "AITKEN",
    { aqueous_phase }
  };
  
  // Another aerosol mode with 2-moment distribution
  auto accumulation = Distribution<shape::LogNormal, moment::Two>{
    "ACCUMULATION",
    { aqueous_phase, organic_phase }  // Multiple phases
  };

  // Dust particle section with 2-moment delta-function distribution
  auto dust = Distribution<shape::DeltaFunction, moment::Two>{
    "DUST",
    { organic_phase }
  };

  auto system = System({
    .gas_phase_ = gas_phase,
    .external_models_ = { small_drop, large_drop, aitken, accumulation, dust }
  });

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

  // // Condensed phase reversible reaction
  // // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // // k_r = reverse rate constant
  // // (k_f = K_eq * k_r = forward rate constant)
  Process h2o_dissociation = DissolvedReversibleReactionBuilder()
                             .SetPhase(aqueous_phase)
                             .SetReactants({ h2o })
                             .SetProducts({ ohm, hp })
                             .SetSolvent(h2o)
                             .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 1.14e-2, .C_ = 2300.0, .T0_ = 298.15 }))
                             .SetReverseRateConstant(ArrheniusRateConstant({ .A_ = 1.4e11, .E_a_ = 5.1e4 }))
                             .Build();

  // Condensed phase reversible reaction: CO2 hydration
  // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  Process co2_hydration = DissolvedReversibleReactionBuilder()
                          .SetPhase(aqueous_phase)
                          .SetReactants({ co2, h2o })
                          .SetProducts({ h2co3 })
                          .SetSolvent(h2o)
                          .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 1.70e3, .C_ = 2400.0, .T0_ = 298.15 }))
                          .SetReverseRateConstant(ArrheniusRateConstant({ .A_ = 1.4e11, .E_a_ = 5.1e4 }))
                          .Build();

  // Condensed phase reversible reaction: H2CO3 dissociation
  // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  Process h2co3_dissociation = DissolvedReversibleReactionBuilder()
                               .SetPhase(aqueous_phase)
                               .SetReactants({ h2co3 })
                               .SetProducts({ hco3m, hp })
                               .SetSolvent(h2o)
                               .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 4.27e2, .C_ = 2300.0, .T0_ = 298.15 }))
                               .SetReverseRateConstant(ArrheniusRateConstant({ .A_ = 2.5e10, .E_a_ = 4.0e4 }))
                               .Build();

  // Condensed phase reversible reaction: HCO3- dissociation
  // K_eq = A * exp( C ( 1 / T0 - 1 / T ) ) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  Process hco3m_dissociation = DissolvedReversibleReactionBuilder()
                               .SetPhase(aqueous_phase)
                               .SetReactants({ hco3m })
                               .SetProducts({ co32m, hp })
                               .SetSolvent(h2o)
                               .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 1.70e1, .C_ = 2300.0, .T0_ = 298.15 }))
                               .SetReverseRateConstant(ArrheniusRateConstant({ .A_ = 6.4e9, .E_a_ = 3.1e4 }))
                               .Build();
#endif
  
  std::vector<Process> reactions{
//    co2_photo,
//    co2_phase_transfer,
//    h2o_dissociation,
//    co2_hydration,
//    h2co3_dissociation,
//    hco3m_dissociation
  };


  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .SetReactions(reactions)
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

  // shape parameters: single moment log-normal distributions
  small_drop.SetParameter(state, small_drop.Shape().GeometricMeanRadius(), 5.0e-6); // m
  large_drop.SetParameter(state, large_drop.Shape().GeometricMeanRadius(), 1.0e-5); // m
  small_drop.SetParameter(state, small_drop.Shape().GeometricStandardDeviation(), 1.6); // unitless
  large_drop.SetParameter(state, large_drop.Shape().GeometricStandardDeviation(), 1.8); // unitless

  // shape parameters: two moment log-normal distributions
  state[aitken.Shape().NumberConcentration()] = 1.0e8; // m-3
  state[accumulation.Shape().NumberConcentration()] = 1.0e7; // m-3
  aitken.SetParameter(state, aitken.Shape().GeometricStandardDeviation(), 1.6); // unitless
  accumulation.SetParameter(state, accumulation.Shape().GeometricStandardDeviation(), 1.8); // unitless

  // shape parameters: two moment delta-function distribution
  state[dust.Shape().NumberConcentration()] = 1.0e6; // m-3
  
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