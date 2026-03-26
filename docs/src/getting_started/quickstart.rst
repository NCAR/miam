==========
Quickstart
==========

This example sets up a cloud droplet system and models CO₂ dissolving from
the gas phase into aqueous droplets via Henry's Law phase transfer.

Complete example
================

Save the following as ``cloud_chem.cpp``:

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
     // Define species with physical properties required for mass transfer
     auto co2 = Species{ "CO2",
         { { "molecular weight [kg mol-1]", 0.044 },
           { "density [kg m-3]", 1800.0 } } };
     auto h2o = Species{ "H2O",
         { { "molecular weight [kg mol-1]", 0.018 },
           { "density [kg m-3]", 1000.0 } } };

     // Define phases
     Phase gas_phase{ "GAS", { { co2 } } };
     Phase aqueous_phase{ "AQUEOUS", { { co2 }, { h2o } } };

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
       .SetGasSpecies(co2)
       .SetCondensedSpecies(co2)
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

     state[co2] = 1.0e-3;                                  // mol m-3 air
     state[cloud.Species(aqueous_phase, h2o)] = 300.0;      // mol m-3 (liquid water content)
     cloud.SetDefaultParameters(state);

     // Integrate
     state.PrintHeader();
     state.PrintState(0);
     for (int i = 1; i <= 10; ++i)
     {
       solver.CalculateRateConstants(state);
       auto result = solver.Solve(0.1, state);  // 100 ms steps
       state.PrintState(i * 100);
     }
   }

Build and run:

.. code-block:: bash

   g++ -o cloud_chem cloud_chem.cpp \
       -I/usr/local/miam/include -I/usr/local/micm/include -std=c++20
   ./cloud_chem

Expected output — gas-phase CO₂ dissolves into the cloud droplets, approaching
Henry's Law equilibrium:

.. code-block:: text

     time,                CO2,  CLOUD.AQUEOUS.CO2,  CLOUD.AQUEOUS.H2O
        0,           1.00e-03,           0.00e+00,           3.00e+02
      100,           9.12e-04,           8.85e-05,           3.00e+02
      200,           8.48e-04,           1.52e-04,           3.00e+02
      300,           8.03e-04,           1.97e-04,           3.00e+02
      400,           7.47e-04,           2.53e-04,           3.00e+02
      500,           7.30e-04,           2.70e-04,           3.00e+02
      600,           7.18e-04,           2.82e-04,           3.00e+02
      700,           7.09e-04,           2.91e-04,           3.00e+02
      800,           7.03e-04,           2.97e-04,           3.00e+02
      900,           6.98e-04,           3.02e-04,           3.00e+02
     1000,           6.98e-04,           3.02e-04,           3.00e+02

What happens
============

1. **Species and phases** are declared with physical properties (molecular
   weight, density) that the aerosol property providers need. The same
   ``CO2`` species appears in both the gas and aqueous phases.

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
