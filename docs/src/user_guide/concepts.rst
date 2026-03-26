========
Concepts
========

MIAM (Model-Independent Aerosol Module) is organized around three core
abstractions: **Representations**, **Processes**, and the **Model** that
orchestrates them.

Architecture Overview
=====================

.. code-block:: text

   ┌─────────────────────────────────────┐
   │           MICM Solver               │
   │  (Rosenbrock, BDF, etc.)            │
   ├─────────────────────────────────────┤
   │  ExternalModelProcessSet interface  │
   │    ForcingFunction(params, vars)     │
   │    JacobianFunction(params, vars)    │
   ├─────────────────────────────────────┤
   │            MIAM Model               │
   │  ┌─────────────┐ ┌──────────────┐   │
   │  │Representations│ │  Processes    │  │
   │  │             │ │              │   │
   │  │ - Modes     │ │ - Dissolved  │   │
   │  │ - Sections  │ │   Reactions  │   │
   │  │             │ │ - Phase      │   │
   │  │ Providers:  │ │   Transfer   │   │
   │  │  r_eff, N,  │◄┤              │   │
   │  │  φ_p        │ │              │   │
   │  └─────────────┘ └──────────────┘   │
   └─────────────────────────────────────┘

Representations
---------------

A **representation** describes a particle size distribution. It owns a set
of **phases** (e.g., aqueous, organic), each containing **species** whose
concentrations are state variables. Representations provide:

- **State variables** — species concentrations plus (for TwoMomentMode) a
  number concentration variable.
- **State parameters** — fixed distribution parameters (geometric mean
  radius, geometric standard deviation, section bounds).
- **Aerosol property providers** — functions that compute effective radius,
  number concentration, and phase volume fraction from the state, along
  with their partial derivatives for the Jacobian.

MIAM provides three representation types:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``SingleMomentMode``
     - Log-normal distribution with fixed size parameters (GMD, GSD).
       Number concentration is diagnosed from total volume.
   * - ``TwoMomentMode``
     - Log-normal distribution with prognostic number concentration.
       Effective radius depends on both volume and number.
   * - ``UniformSection``
     - Sectional bin with fixed radius bounds. Number concentration is
       diagnosed from total volume.

Processes
---------

A **process** describes a transformation of species concentrations. Each
process implements a common interface that provides:

- **Forcing function** — the ODE right-hand side contribution.
- **Jacobian function** — the partial derivatives of the forcing with
  respect to state variables (stored as −J per MICM convention).
- **State parameter management** — temperature-dependent rate constants
  are updated each time step before the solver iterates.

Processes can **require aerosol properties** from representations. When
they do, the Model automatically builds and routes the appropriate
providers.

MIAM currently provides two process types:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Type
     - Description
   * - ``DissolvedReversibleReaction``
     - Equilibrium reactions within a condensed phase
       (e.g., H₂O ⇌ OH⁻ + H⁺). Does not require aerosol properties.
   * - ``HenryLawPhaseTransfer``
     - Gas ⇌ condensed-phase mass transfer governed by Henry's Law.
       Requires effective radius, number concentration, and phase volume
       fraction from the representation.

Model
-----

The **Model** is the top-level container. It:

1. Stores representations and processes.
2. Collects all state variable and parameter names for the solver.
3. Builds aerosol property providers by matching process requirements to
   the representations that own each phase prefix.
4. Returns forcing and Jacobian functions that iterate over all processes
   in a single loop.

A Model is registered with MICM as an ``ExternalModelSystem`` and
``ExternalModelProcessSet``, making MIAM's state seamlessly part of the
coupled solver state.

Data Flow
=========

Each solver time step follows this sequence:

1. **Conditions update** — the host model sets temperature and pressure
   on the MICM ``State``.

2. **Rate constant calculation** — ``solver.CalculateRateConstants(state)``
   calls MIAM's ``UpdateStateParametersFunction``, which evaluates
   temperature-dependent constants (HLC, K_eq, k_f, k_r) and writes them
   into state parameter columns.

3. **Solver iteration** — the Rosenbrock solver calls MIAM's
   ``ForcingFunction`` and ``JacobianFunction`` at each internal stage.
   These functions:

   a. Query aerosol property providers for r_eff, N, φ_p (and their
      partial derivatives, for the Jacobian).
   b. Compute condensation/evaporation rates.
   c. Accumulate forcing terms and Jacobian entries.

4. **State update** — the solver writes the updated concentrations back
   into the ``State`` object.

Species Naming Convention
=========================

State variable names follow a hierarchical pattern:

.. code-block:: text

   <mode_name>.<phase_name>.<species_name>

For example, a model named ``"CLOUD"`` with a ``SingleMomentMode`` named
``"DROPLET"`` containing an aqueous phase produces variables like:

- ``DROPLET.AQUEOUS.CO2``
- ``DROPLET.AQUEOUS.H2O``

Gas-phase species are registered directly by MICM with their short name
(e.g., ``CO2``). When passed to MICM's ``System`` constructor, the model
name is **not** prepended — only the representation and phase names form
the prefix.

Two-moment modes also introduce a number concentration variable:

- ``AITKEN.number_concentration``
