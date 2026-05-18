================
Building a Model
================

This guide shows how to assemble representations and processes into a
complete MIAM Model and wire it into MICM's solver.

Basic Setup
===========

.. code-block:: c++

   #include <miam/miam.hpp>
   #include <micm/CPU.hpp>
   using namespace micm;
   using namespace miam;

A Model is constructed with a name and a list of representations:

.. code-block:: c++

   auto model = Model{
     .name_ = "AEROSOL",
     .representations_ = { aitken, accumulation, dust }
   };

Processes are added separately:

.. code-block:: c++

   model.AddProcesses({ h2o_dissociation, co2_hydration });

Multi-Model Systems
===================

Multiple Models can represent different particle categories (e.g., cloud
droplets vs. background aerosol). Each Model is registered independently
with MICM:

.. code-block:: c++

   auto cloud_model = Model{
     .name_ = "CLOUD",
     .representations_ = { small_drop, large_drop }
   };
   cloud_model.AddProcesses({ reactions });

   auto aerosol_model = Model{
     .name_ = "AEROSOL",
     .representations_ = { aitken, accumulation }
   };
   aerosol_model.AddProcesses({ reactions });

MICM Integration
================

Register Models as external model systems when building the solver:

.. code-block:: c++

   auto system = System(gas_phase, cloud_model, aerosol_model);

   auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
       RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
       .SetSystem(system)
       .AddExternalModel(cloud_model)
       .AddExternalModel(aerosol_model)
       .SetIgnoreUnusedSpecies(true)
       .Build();

``SetSystem`` registers all state variables (gas-phase species from MICM
plus MIAM's species and aerosol parameters). ``AddExternalModel``
registers each Model's forcing and Jacobian functions with the solver.

State Initialization
====================

After building the solver, initialize the state:

.. code-block:: c++

   State state = solver.GetState();

   // Environmental conditions
   state.conditions_[0].temperature_ = 287.45;  // K
   state.conditions_[0].pressure_ = 101319.9;   // Pa
   state.conditions_[0].CalculateIdealAirDensity();

   // Gas-phase species (short name)
   state[co2] = 0.2;  // mol m-3

   // Condensed-phase species (fully qualified name)
   state[droplets.Species(aqueous_phase, h2o)] = 300.0;  // mol m-3 (liquid water content)

   // Number concentration (TwoMomentMode only)
   state[aitken.NumberConcentration()] = 1.0e8;  // m-3

   // Representation parameters (GMD, GSD, section bounds)
   droplets.SetDefaultParameters(state);
   aitken.SetDefaultParameters(state);

Time Integration
================

.. code-block:: c++

   for (int i = 0; i < 100; ++i)
   {
     solver.CalculateRateConstants(state);
     auto result = solver.Solve(10.0, state);  // 10-second steps
   }

``CalculateRateConstants`` calls MIAM's ``UpdateStateParametersFunction``,
which evaluates temperature-dependent constants (HLC, K_eq, k_f, k_r).
The Rosenbrock solver then calls MIAM's forcing and Jacobian functions
internally at each stage.
