===============
Representations
===============

Representations describe how particle populations are distributed in size
space. Each representation owns one or more phases, manages the
corresponding state variables and parameters, and provides aerosol property
calculations.

SingleMomentMode
================

A log-normal size distribution with **fixed** geometric mean diameter (GMD)
and geometric standard deviation (GSD). The number concentration is
**diagnosed** from the total species volume.

.. code-block:: c++

   auto cloud = representation::SingleMomentMode{
     "CLOUD",                // mode name
     { aqueous_phase },      // phases
     1.0e-6,                 // geometric mean radius [m]
     1.4                     // geometric standard deviation
   };

**State variables**: species concentrations for each phase (no number
concentration variable).

**State parameters** (2): ``geometric_mean_radius``, ``geometric_standard_deviation``.

**Aerosol properties**:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Property
     - Formula
     - State dependencies
   * - Effective radius
     - :math:`r_\text{eff} = \text{GMD} \cdot \exp(2.5 \ln^2 \text{GSD})`
     - None (parameterized)
   * - Number concentration
     - :math:`N = V_\text{total} / V_\text{single}`
     - All species (all phases)
   * - Phase volume fraction
     - :math:`\varphi_p = V_\text{phase} / V_\text{total}`
     - All species (all phases)

TwoMomentMode
=============

A log-normal distribution with **prognostic** number concentration and
fixed GSD. Effective radius is derived from the total volume and number
concentration.

.. code-block:: c++

   auto aitken = representation::TwoMomentMode{
     "AITKEN",               // mode name
     { aqueous_phase },      // phases
     1.2                     // geometric standard deviation
   };

   // Set the number concentration in the solver state:
   state[aitken.NumberConcentration()] = 1.0e8;  // m-3

**State variables**: species concentrations + one ``number_concentration``
variable per mode.

**State parameters** (1): ``geometric_standard_deviation``.

**Aerosol properties**:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Property
     - Formula
     - State dependencies
   * - Effective radius
     - :math:`r_\text{eff} = (3 V / 4\pi N)^{1/3} \cdot \exp(2.5 \ln^2 \text{GSD})`
     - All species + N
   * - Number concentration
     - :math:`N = [\text{N\_var}]`
     - Number concentration variable
   * - Phase volume fraction
     - :math:`\varphi_p = V_\text{phase} / V_\text{total}`
     - All species (all phases)

UniformSection
==============

A sectional bin with **fixed** minimum and maximum radius bounds. Number
concentration is diagnosed from total species volume and average particle
size.

.. code-block:: c++

   auto dust = representation::UniformSection{
     "DUST",                 // section name
     { organic_phase },      // phases
     1.0e-7,                 // min radius [m]
     1.0e-6                  // max radius [m]
   };

**State variables**: species concentrations for each phase.

**State parameters** (2): ``min_radius``, ``max_radius``.

**Aerosol properties**:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Property
     - Formula
     - State dependencies
   * - Effective radius
     - :math:`r_\text{eff} = (r_\min + r_\max) / 2`
     - None (parameterized)
   * - Number concentration
     - :math:`N = V_\text{total} / V_\text{single}`
     - All species (all phases)
   * - Phase volume fraction
     - :math:`\varphi_p = V_\text{phase} / V_\text{total}`
     - All species (all phases)

Multi-Phase Modes
=================

A representation can host multiple phases. For example, an accumulation
mode with both aqueous and organic phases:

.. code-block:: c++

   auto accumulation = representation::TwoMomentMode{
     "ACCUMULATION",
     { aqueous_phase, organic_phase },  // two phases
     1.4
   };

All species across all phases contribute to the total particle volume used
for size calculations. The phase volume fraction φ_p partitions the
particle surface area among phases (see
:doc:`../science_guide/henry_law_phase_transfer` § Phase volume fraction).

Species Property Requirements
=============================

Species used in representations must carry physical properties for volume
calculations:

.. code-block:: c++

   auto nacl = Species{ "NaCl",
       { { "molecular weight [kg mol-1]", 0.058 },
         { "density [kg m-3]", 2165.0 } } };

These are accessed at provider creation time (once, not on the hot path).

Accessing State Variables
=========================

Representations provide helper methods for constructing fully qualified
state variable names:

.. code-block:: c++

   // Returns "DROPLET.AQUEOUS.CO2"
   std::string name = droplet.Species(aqueous_phase, co2);

   // For TwoMomentMode: returns "AITKEN.number_concentration"
   std::string nc_name = aitken.NumberConcentration();

Use these with the MICM ``State`` object to set initial conditions:

.. code-block:: c++

   state[droplet.Species(aqueous_phase, co2)] = 1.0e-5;  // mol m-3
