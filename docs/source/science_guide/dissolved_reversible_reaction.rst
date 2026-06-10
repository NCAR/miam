================================
Dissolved Reversible Reactions
================================

Equilibrium reactions within a condensed phase, following the form:

.. math::

   \text{Reactants} \rightleftharpoons \text{Products}

Equilibrium Constant
====================

Temperature-dependent:

.. math::

   K_\text{eq}(T) = A \cdot \exp\!\left(C \left(\frac{1}{T_0} - \frac{1}{T}\right)\right)

At each time step, the forward rate constant is derived from:

.. math::

   k_f = K_\text{eq} \cdot k_r

where :math:`k_r` is specified directly (e.g., via an Arrhenius rate
constant). Alternatively, :math:`k_f` can be provided explicitly.

Rate Expression
===============

For a generic reaction with reactants :math:`R_1, R_2, \ldots` and
products :math:`P_1, P_2, \ldots`:

.. math::

   r_f = k_f \prod_i [R_i], \qquad r_r = k_r \prod_j [P_j]

Special treatment of solvents
-----------------------------

When a species is designated as the **solvent**, its concentration is *not*
included as a separate multiplier in the rate expression. Instead, the
solvent's activity is implicitly set to unity (dilute solution
approximation). This means that reactions like:

.. math::

   \text{H₂O} \rightleftharpoons \text{OH⁻} + \text{H⁺}

use :math:`r_f = k_f` (not :math:`r_f = k_f [\text{H₂O}]`), because
water is the overwhelmingly dominant species and its activity ≈ 1.

Net Rate and Forcing
====================

.. math::

   R_\text{net} = r_f - r_r = k_f \prod_i [R_i] - k_r \prod_j [P_j]

For each reactant species:

.. math::

   \frac{d[R_i]}{dt} = -R_\text{net}

For each product species:

.. math::

   \frac{d[P_j]}{dt} = +R_\text{net}

Jacobian
========

The MICM solver stores **−J**. For each reactant–variable pair:

.. math::

   \text{stored } {-J[R_i, x]} = -\frac{\partial (-R_\text{net})}{\partial x}

For a first-order forward reaction (:math:`R_\text{net} = k_f - k_r [P]`):

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Entry
     - Mathematical :math:`J`
     - Stored :math:`-J`
   * - :math:`J[P, P]`
     - :math:`-k_r`
     - :math:`+k_r`

For a second-order reverse reaction (:math:`r_r = k_r [P_1][P_2]`):

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Entry
     - Mathematical :math:`J`
     - Stored :math:`-J`
   * - :math:`J[R, P_1]`
     - :math:`+k_r [P_2]`
     - :math:`-k_r [P_2]`
   * - :math:`J[R, P_2]`
     - :math:`+k_r [P_1]`
     - :math:`-k_r [P_1]`
   * - :math:`J[P_1, P_1]`
     - :math:`-k_r [P_2]`
     - :math:`+k_r [P_2]`
   * - :math:`J[P_1, P_2]`
     - :math:`-k_r [P_1]`
     - :math:`+k_r [P_1]`

The general pattern follows standard chemical kinetics Jacobian entries.

No Aerosol Properties Required
==============================

Unlike Henry's Law phase transfer, dissolved reversible reactions are
purely chemical and do not depend on particle size or number. They produce
no indirect Jacobian entries through aerosol properties.

Multiple Instances
==================

When the same reaction occurs in multiple condensed-phase instances
(e.g., the same aqueous chemistry in several aerosol modes), each instance
operates independently on its own set of fully-qualified state variables
(e.g., ``AITKEN.AQUEOUS.H2O`` vs ``ACCUMULATION.AQUEOUS.H2O``).
