===================
Aerosol Properties
===================

Aerosol properties bridge representations and processes. Processes declare
which properties they need; the Model builds providers by querying the
representations.

AerosolProperty Enum
====================

MIAM defines three aerosol properties:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Property
     - Units
     - Description
   * - ``EffectiveRadius``
     - m
     - The effective radius of particles in the mode or section
   * - ``NumberConcentration``
     - # m⁻³
     - The total number concentration of particles
   * - ``PhaseVolumeFraction``
     - dimensionless
     - Fraction of total particle volume belonging to a specific phase

Provider Pattern
================

An ``AerosolPropertyProvider`` is a lightweight callable object created at
setup time by a representation. It captures all necessary column indices
so that runtime calls are free of string lookups or map searches.

Each provider exposes:

``dependent_variable_indices``
   A vector of state variable indices that the property depends on. Used
   for Jacobian sparsity analysis and partial derivative mapping.

``ComputeValue(params, vars, result)``
   Computes the property value for all grid cells. Called inside the
   forcing function loop.

``ComputeValueAndDerivatives(params, vars, result, partials)``
   Computes the property value **and** partial derivatives with respect to
   each dependent variable. Called inside the Jacobian function loop. The
   ``partials`` matrix has one column per dependent variable, in the same
   order as ``dependent_variable_indices``.

How Providers Are Built
=======================

When you register processes and representations with a Model, the Model's
``BuildProviders`` method:

1. Collects ``RequiredAerosolProperties()`` from each process.
2. De-duplicates properties per phase.
3. Finds the representation that owns each phase prefix.
4. Calls ``GetPropertyProvider(property, param_indices, var_indices,
   phase_name)`` on the representation.
5. Stores the resulting providers in a map keyed by phase prefix.

This happens once at solver setup time. The providers are then passed to
the process ``ForcingFunction`` and ``JacobianFunction`` closures.

Derivative Chain
================

For the Jacobian, MIAM implements a chain rule through aerosol properties.
When a process depends on an aerosol property that itself depends on state
variables, the total Jacobian contribution is:

.. math::

   \frac{\partial f_i}{\partial x_j} =
     \underbrace{\frac{\partial f_i}{\partial P}}_{\text{process}}
     \cdot
     \underbrace{\frac{\partial P}{\partial x_j}}_{\text{provider}}

where :math:`P` is an aerosol property (e.g., effective radius) and
:math:`x_j` is a state variable (e.g., species concentration).

The process computes :math:`\partial f_i / \partial P` (e.g., how forcing
depends on effective radius), and the provider supplies
:math:`\partial P / \partial x_j` via the ``partials`` matrix. The process
then assembles the complete Jacobian entries by multiplying these factors.

Per-Representation Derivatives
==============================

The partial derivatives :math:`\partial P / \partial x_j` vary by
representation type. For details, see
:doc:`../science_guide/henry_law_phase_transfer` § Aerosol property
derivatives.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Representation
     - :math:`\partial r_\text{eff} / \partial x_j`
     - Notes
   * - SingleMomentMode
     - 0
     - Radius is fixed by GMD/GSD parameters
   * - TwoMomentMode
     - Non-zero (depends on V, N)
     - Chain rule through total volume and number concentration
   * - UniformSection
     - 0
     - Radius is fixed by section bounds
