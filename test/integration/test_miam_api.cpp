#include <miam/model/aerosol_scheme.hpp>
#include <miam/model/mode.hpp>
#include <miam/model/section.hpp>
#include <miam/util/solver_utils.hpp>

#include <micm/process/transfer_coefficient/henrys_law_constant.hpp>
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/Process.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

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
  auto h2co3 = Species{ "H2CO3" };

  Phase gas_phase{ "GAS", { { co2, 31.2 } } };
  Phase aqueous_phase{ "AQUEOUS", { { co2 } , { h2o } , { ohm } , { hp } , { hco3m } , { co32m }, { h2co3 }} };
  Phase organic_phase{ "ORGANIC", { { co2, 16.2 }, { hexane } } };

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

  // Model aerosol with 2 modes and 1 section
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

  System chemical_system = ConfigureSystem(gas_phase, { small_drop, large_drop, aitken, accumulation }, { dust });

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
  // 14) (AEROSOL.) AITKEN.AQUEOUS.CO2           ( Mode )
  // 15) (AEROSOL.) AITKEN.AQUEOUS.HEXANE        ( Mode )
  // 16) (AEROSOL.) ACCUMULATION.AQUEOUS.CO2     ( Mode )
  // 17) (AEROSOL.) ACCUMULATION.AQUEOUS.H2O     ( Mode )
  // 18) (AEROSOL.) ACCUMULATION.AQUEOUS.OH-     ( Mode )
  // 19) (AEROSOL.) ACCUMULATION.AQUEOUS.H+      ( Mode )
  // 20) (AEROSOL.) ACCUMULATION.AQUEOUS.HCO3-   ( Mode )
  // 21) (AEROSOL.) ACCUMULATION.AQUEOUS.CO32-   ( Mode )
  // 22) (AEROSOL.) ACCUMULATION.ORGANIC.CO2     ( Mode )
  // 23) (AEROSOL.) ACCUMULATION.ORGANIC.HEXANE  ( Mode )
  // 24) (AEROSOL.) DUST.ORGANIC.CO2             ( Section )
  // 25) (AEROSOL.) DUST.ORGANIC.HEXANE          ( Section )

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
  Process co2_phase_transfer = PhaseTransferProcessBuilder()
                               .SetGasSpecies(gas_phase, co2 )
                               .SetCondensedSpecies(aqueous_phase, h2co3)
                               .SetSolvent(aqueous_phase, h2o )
                               .SetTransferCoefficient(
                                  HenrysLawConstant{{
                                    .H_ref_ = 1.5e-3,
                                    .enthalpy_ = -10000.0,
                                    .temperature_ref_ = 298.15 }})
                               .Build();

  // First reversible acid-dissociation reaction `H2CO3(aq) <-> HCO3-(aq) + H+(aq)`
  // auto h2co3_hco3m_hp = DissolvedReversibleProcessBuilder()
  Process h2co3_hco3m_hp = ChemicalReactionBuilder()
                            .SetPhase(aqueous_phase)
                            .SetReactants({ h2co3 })
                            .SetProducts({ { hco3m, 1.0 }, { hp, 1.0 } })
                            .SetRateConstant(ArrheniusRateConstant({ .A_ = 1.0e-3 }))
                            .Build();

  // Second reversible acid-dissociation reaction `HCO3-(aq) <-> CO3--(aq) + H+(aq)`
  Process hco3m_co32m_hp = ChemicalReactionBuilder()
                            .SetPhase(aqueous_phase)
                            .SetReactants({ hco3m })
                            .SetProducts({ { co32m, 1.0 }, {  hp, 1.0 } })
                            .SetRateConstant(ArrheniusRateConstant({ .A_ = 1.0e-3 }))
                            .Build();

  std::vector<Process> reactions{ co2_photo, co2_phase_transfer, h2co3_hco3m_hp, hco3m_co32m_hp };

  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(chemical_system)
                  .SetReactions(reactions)
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
  state.SetConcentration(co2, 0.2);
  
  // cloud
  state.SetConcentration(small_drop.Species(aqueous_phase, h2o), 0.3);      // mol m-3
  state.SetConcentration(large_drop.Species(aqueous_phase, h2o), 0.3);      // mol m-3
  
  // aerosol mode
  state.SetConcentration(aitken.Species(aqueous_phase, hco3m), 0.1);        // mol m-3
  state.SetConcentration(accumulation.Species(aqueous_phase, h2o), 0.3);    // mol m-3
  state.SetConcentration(accumulation.Species(organic_phase, hexane), 0.1); // mol m-3
  
  // aerosol section
  state.SetConcentration(dust.Species(organic_phase, co2), 0.1);            // mol m-3

  // number concentration
  state.SetConcentration(aitken.NumberConcentration(), 1.0e4);        // m-3
  state.SetConcentration(accumulation.NumberConcentration(), 1.0e3);  // m-3

  // density
  state.SetConcentration(aitken.Density(), 2.0e2);        // m-3
  state.SetConcentration(accumulation.Density(), 3.0e2);  // m-3

  // radius
  double radius_small_drop = small_drop.GetRadius(state);      // m
  double radius_large_drop = large_drop.GetRadius(state);      // m
  double radius_aitken = aitken.GetRadius(state);              // m
  double radius_accumulation = accumulation.GetRadius(state);  // m

  state.SetConcentration(small_drop.Radius(), radius_small_drop);      // m
  state.SetConcentration(large_drop.Radius(), radius_large_drop);      // m
  state.SetConcentration(aitken.Radius(), radius_aitken);              // m
  state.SetConcentration(accumulation.Radius(), radius_accumulation);  // m

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