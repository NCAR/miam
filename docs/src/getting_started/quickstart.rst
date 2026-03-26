==========
Quickstart
==========

This example sets up a minimal aerosol system: a gas-phase species (A)
partitioning into cloud droplets via Henry's Law, then solves the coupled
system using MICM's Rosenbrock solver.

Complete example
================

Save the following as ``henry_law_example.cpp``:

.. code-block:: c++

   #include <miam/miam.hpp>
   #include <miam/processes/constants/henrys_law_constant.hpp>
   #include <micm/CPU.hpp>

   #include <iomanip>
   #include <iostream>

   using namespace micm;
   using namespace miam;

   int main()
   {
     // 1. Define species with physical properties
     auto A_gas = Species{ "A_g",
         { { "molecular weight [kg mol-1]", 0.044 } } };
     auto A_aq = Species{ "A_aq",
         { { "molecular weight [kg mol-1]", 0.044 },
           { "density [kg m-3]", 1800.0 } } };
     auto H2O = Species{ "H2O",
         { { "molecular weight [kg mol-1]", 0.018 },
           { "density [kg m-3]", 1000.0 } } };

     // 2. Define phases
     Phase gas_phase{ "GAS", { { A_gas } } };
     Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

     // 3. Create a particle representation (cloud droplets)
     auto droplets = representation::SingleMomentMode{
       "DROPLET",
       { aqueous_phase },
       5.0e-6,  // geometric mean radius [m]
       1.2      // geometric standard deviation
     };

     // 4. Define the Henry's Law phase transfer process
     auto transfer = process::HenryLawPhaseTransferBuilder()
       .SetCondensedPhase(aqueous_phase)
       .SetGasSpecies(A_gas)
       .SetGasSpeciesName("A_g")
       .SetCondensedSpecies(A_aq)
       .SetSolvent(H2O)
       .SetHenrysLawConstant(process::constant::HenrysLawConstant(
           { .HLC_ref_ = 3.4e-2 }))  // mol m-3 Pa-1
       .SetDiffusionCoefficient(1.5e-5)   // m2 s-1
       .SetAccommodationCoefficient(0.05)
       .Build();

     // 5. Assemble model
     auto model = Model{
       .name_ = "CLOUD",
       .representations_ = { droplets }
     };
     model.AddProcesses({ transfer });

     // 6. Build solver
     auto system = System(gas_phase, model);
     auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
         RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
         .SetSystem(system)
         .AddExternalModelProcesses(model)
         .SetIgnoreUnusedSpecies(true)
         .Build();

     // 7. Set initial conditions
     State state = solver.GetState();
     state.conditions_[0].temperature_ = 298.15;
     state.conditions_[0].pressure_ = 101325.0;
     state.conditions_[0].CalculateIdealAirDensity();

     state[A_gas] = 1.0e-3;   // mol m-3 air
     state[droplets.Species(aqueous_phase, H2O)] = 55.5;
     droplets.SetDefaultParameters(state);

     // 8. Integrate
     state.PrintHeader();
     for (int i = 0; i < 20; ++i)
     {
       solver.CalculateRateConstants(state);
       solver.Solve(0.05, state);  // 50 ms steps
       state.PrintState(i * 50);
     }
   }

Build and run:

.. code-block:: bash

   g++ -o henry_law_example henry_law_example.cpp \
       -I/usr/local/miam/include -I/usr/local/micm/include -std=c++20
   ./henry_law_example

What happens
============

1. **Species and phases** are declared with physical properties (molecular
   weight, density) that the aerosol property providers need.

2. A **SingleMomentMode** representation describes a log-normal droplet
   population with fixed geometric mean radius and standard deviation.
   The model derives effective radius, number concentration, and phase
   volume fraction from these parameters and the species concentrations.

3. A **HenryLawPhaseTransfer** process is configured with a Henry's Law
   constant, gas-phase diffusion coefficient, and mass accommodation
   coefficient. The builder extracts molecular weights and densities from
   the species objects.

4. The **Model** bundles representations and processes. When passed to
   MICM's ``CpuSolverBuilder``, it exposes its state variables and
   parameters through MICM's ``ExternalModelSystem`` interface.

5. The **Rosenbrock solver** integrates the coupled ODE system, calling
   MIAM's forcing and Jacobian functions at each solver iteration.

Next steps
==========

- :doc:`../user_guide/concepts` — understand the architecture
- :doc:`../user_guide/representations` — all three representation types
- :doc:`../user_guide/processes` — all process types and their builders
- :doc:`../science_guide/henry_law_phase_transfer` — the full equation set
