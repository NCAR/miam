=========
Processes
=========

Processes define the chemical and physical transformations of species
concentrations. Each process provides forcing (ODE right-hand side) and
Jacobian contributions to the solver.

DissolvedReversibleReaction
===========================

Equilibrium reactions within a condensed phase:

.. math::

   \text{Reactants} \rightleftharpoons \text{Products}

with :math:`K_\text{eq} = k_f / k_r`.

.. code-block:: c++

   #include <miam/miam.hpp>

   auto h2o_dissociation = process::DissolvedReversibleReactionBuilder()
     .SetPhase(aqueous_phase)
     .SetReactants({ h2o })
     .SetProducts({ ohm, hp })
     .SetSolvent(h2o)
     .SetEquilibriumConstant(process::constant::EquilibriumConstant(
         { .A_ = 1.14e-2, .C_ = 2300.0, .T0_ = 298.15 }))
     .SetReverseRateConstant(ArrheniusRateConstant(
         { .A_ = 1.4e11, .C_ = 5.1e4 }))
     .Build();

Builder methods
---------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``SetPhase(phase)``
     - The condensed phase where the reaction occurs
   * - ``SetReactants(species_list)``
     - Reactant species
   * - ``SetProducts(species_list)``
     - Product species
   * - ``SetSolvent(species)``
     - Solvent species (concentration used for activity)
   * - ``SetEquilibriumConstant(obj)``
     - Temperature-dependent :math:`K_\text{eq}(T)`
   * - ``SetForwardRateConstant(obj)``
     - Explicit :math:`k_f(T)` (alternative to equilibrium constant)
   * - ``SetReverseRateConstant(obj)``
     - Explicit :math:`k_r(T)`

If both ``SetEquilibriumConstant`` and ``SetReverseRateConstant`` are
provided, the forward rate constant is computed as
:math:`k_f = K_\text{eq} \cdot k_r`.

This process does **not** require aerosol properties.

HenryLawPhaseTransfer
=====================

Gas-condensed phase mass transfer governed by Henry's Law:

.. math::

   \text{A(gas)} \rightleftharpoons \text{A(condensed)}

.. code-block:: c++

   #include <miam/miam.hpp>
   #include <miam/processes/constants/henrys_law_constant.hpp>

   auto transfer = process::HenryLawPhaseTransferBuilder()
     .SetCondensedPhase(aqueous_phase)
     .SetGasSpecies(A_gas)
     .SetGasSpeciesName("A_g")
     .SetCondensedSpecies(A_aq)
     .SetSolvent(H2O)
     .SetHenrysLawConstant(process::constant::HenrysLawConstant(
         { .HLC_ref_ = 3.4e-2, .C_ = 2400.0, .T0_ = 298.15 }))
     .SetDiffusionCoefficient(1.5e-5)    // m2 s-1
     .SetAccommodationCoefficient(0.05)  // dimensionless
     .Build();

Builder methods
---------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``SetCondensedPhase(phase)``
     - Target condensed phase
   * - ``SetGasSpecies(species)``
     - Gas-phase species (must carry ``molecular weight``)
   * - ``SetGasSpeciesName(name)``
     - Name used for the gas species in the state variable map
   * - ``SetCondensedSpecies(species)``
     - Condensed-phase solute species (must carry ``molecular weight``
       and ``density``)
   * - ``SetSolvent(species)``
     - Solvent species (must carry ``molecular weight`` and ``density``)
   * - ``SetHenrysLawConstant(obj)``
     - Temperature-dependent HLC: :math:`H(T) = H_\text{ref} \exp(C(1/T - 1/T_0))`
   * - ``SetDiffusionCoefficient(D_g)``
     - Gas-phase diffusion coefficient [m² s⁻¹]
   * - ``SetAccommodationCoefficient(α)``
     - Mass accommodation coefficient [0–1]

This process **requires** three aerosol properties from the
representation: effective radius, number concentration, and phase volume
fraction. These are provided automatically by the Model.

See :doc:`../science_guide/henry_law_phase_transfer` for the full equation
set and Jacobian derivations.

Rate Constants
==============

EquilibriumConstant
-------------------

Temperature-dependent equilibrium constant:

.. math::

   K_\text{eq}(T) = A \cdot \exp\!\left( C \left( \frac{1}{T_0} - \frac{1}{T} \right) \right)

.. code-block:: c++

   process::constant::EquilibriumConstant keq(
       { .A_ = 1.14e-2, .C_ = 2300.0, .T0_ = 298.15 });

   double value = keq.Calculate(conditions);

HenrysLawConstant
------------------

Temperature-dependent Henry's Law constant:

.. math::

   H(T) = H_\text{ref} \cdot \exp\!\left( C \left( \frac{1}{T} - \frac{1}{T_0} \right) \right)

.. code-block:: c++

   process::constant::HenrysLawConstant hlc(
       { .HLC_ref_ = 3.4e-2, .C_ = 2400.0, .T0_ = 298.15 });

   double value = hlc.Calculate(conditions);  // mol m-3 Pa-1

Adding Processes to a Model
===========================

.. code-block:: c++

   auto model = Model{
     .name_ = "AEROSOL",
     .representations_ = { droplets, dust }
   };

   // Add multiple processes at once
   model.AddProcesses({ h2o_dissociation, co2_hydration, transfer });

``AddProcesses`` assigns each process a unique UUID (via ``CopyWithNewUuid``)
so that the same process definition can be reused across modes without
state parameter name collisions.
