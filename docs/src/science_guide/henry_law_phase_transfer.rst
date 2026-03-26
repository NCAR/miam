===========================
Henry's Law Phase Transfer
===========================

Gas-condensed phase mass transfer governed by Henry's Law equilibrium and
Fuchs-Sutugin transition-regime kinetics.

Physics
=======

For gas species A partitioning into a condensed phase:

.. math::

   \text{A(gas)} \rightleftharpoons \text{A(condensed)}

The condensation rate describes how fast gas molecules transfer into
particles, the evaporation rate describes the reverse, and Henry's Law
defines the equilibrium ratio.

Henry's Law Constant
====================

Temperature-dependent:

.. math::

   H(T) = H_\text{ref} \cdot \exp\!\left(C \left(\frac{1}{T} - \frac{1}{T_0}\right)\right)

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Symbol
     - Description
     - Units
   * - :math:`H_\text{ref}`
     - Reference HLC at :math:`T_0`
     - mol m⁻³ Pa⁻¹
   * - :math:`C`
     - Temperature dependence parameter
     - K
   * - :math:`T_0`
     - Reference temperature (default 298.15)
     - K

HLC is evaluated once per time step and stored as a state parameter.

Condensation Rate
=================

Mean molecular speed
--------------------

.. math::

   \bar{c} = \sqrt{\frac{8 R T}{\pi M_w}}

Mean free path
--------------

.. math::

   \lambda = \frac{3 D_g}{\bar{c}}

Knudsen number
--------------

.. math::

   \text{Kn} = \frac{\lambda}{r_\text{eff}}

Characterizes the gas-particle interaction regime (Kn ≪ 1 → continuum,
Kn ≫ 1 → free molecular).

Fuchs-Sutugin correction
-------------------------

.. math::

   f(\text{Kn}) = \frac{1 + \text{Kn}}{1 + \frac{2 \text{Kn}(1 + \text{Kn})}{\alpha}}

Interpolates between the continuum (:math:`f \to 1`) and free-molecular
(:math:`f \to \alpha / (2 \text{Kn})`) limits.

**Derivative with respect to Kn**:

.. math::

   \frac{df}{d\text{Kn}} = \frac{\alpha - 2\text{Kn}^2 - 2\text{Kn}}{\alpha \cdot D^2}

where :math:`D = 1 + 2\text{Kn}(1 + \text{Kn})/\alpha`.

Condensation rate constant
--------------------------

.. math::

   k_c = 4\pi \, r_\text{eff} \, N \, D_g \, f(\text{Kn})

Units: s⁻¹. First-order rate constant for gas-to-condensed transfer.

**Partial derivatives**:

.. math::

   \frac{\partial k_c}{\partial r_\text{eff}} =
       4\pi N D_g \left[ f(\text{Kn}) - \text{Kn} \cdot \frac{df}{d\text{Kn}} \right]

.. math::

   \frac{\partial k_c}{\partial N} = \frac{k_c}{N}

Evaporation Rate
================

.. math::

   k_e = \frac{k_c}{H \cdot R \cdot T}

**Partial derivatives** — obtained by dividing :math:`k_c` partials by the
same constant factor:

.. math::

   \frac{\partial k_e}{\partial r_\text{eff}} =
       \frac{1}{HRT} \frac{\partial k_c}{\partial r_\text{eff}},\qquad
   \frac{\partial k_e}{\partial N} = \frac{k_e}{N}

Phase Volume Fraction
=====================

For multi-phase modes, the exposed surface area of phase :math:`p` is
proportional to its volume fraction:

.. math::

   \varphi_p = \frac{V_\text{phase}}{V_\text{total}}

where

.. math::

   V_\text{phase} = \sum_i [\text{species}_{p,i}] \cdot \frac{M_{w,i}}{\rho_i}, \qquad
   V_\text{total} = \sum_q \sum_i [\text{species}_{q,i}] \cdot \frac{M_{w,i}}{\rho_i}

Solvent Volume Fraction
=======================

Converts between concentration bases:

.. math::

   f_v = [\text{solvent}] \cdot \frac{M_{w,\text{solvent}}}{\rho_\text{solvent}}

Net Transfer Rate
=================

.. math::

   R_\text{net} = \varphi_p \, k_c \, [\text{A}]_\text{gas}
     - \varphi_p \, k_e \, \frac{[\text{A}]_\text{aq}}{f_v}

The ODE right-hand sides:

.. math::

   \frac{d[\text{A}]_\text{gas}}{dt} = -R_\text{net}, \qquad
   \frac{d[\text{A}]_\text{aq}}{dt} = +R_\text{net}

When multiple condensed-phase instances exist (e.g., multiple aerosol
modes containing the same phase), each contributes its own
:math:`R_\text{net}` and contributions are summed.

Jacobian
========

The MICM Rosenbrock solver stores **−J** (negative Jacobian). All values
below show the mathematical derivative and the stored value.

Direct entries
--------------

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Entry
     - Mathematical :math:`J`
     - Stored :math:`-J`
   * - :math:`J[\text{gas}, \text{gas}]`
     - :math:`-\varphi_p k_c`
     - :math:`+\varphi_p k_c`
   * - :math:`J[\text{gas}, \text{aq}]`
     - :math:`+\varphi_p k_e / f_v`
     - :math:`-\varphi_p k_e / f_v`
   * - :math:`J[\text{gas}, \text{solvent}]`
     - :math:`-\varphi_p k_e [\text{A}]_\text{aq} / (f_v \cdot [\text{solvent}])`
     - :math:`+\varphi_p k_e [\text{A}]_\text{aq} / (f_v \cdot [\text{solvent}])`
   * - :math:`J[\text{aq}, x]`
     - :math:`-J[\text{gas}, x]`
     - :math:`-(-J[\text{gas}, x])`

The antisymmetry :math:`J[\text{aq}, x] = -J[\text{gas}, x]` is a direct
consequence of mass conservation.

Indirect entries (chain rule)
-----------------------------

When an aerosol property :math:`P` depends on state variable :math:`y_j`,
additional terms appear:

**Through** :math:`r_\text{eff}`:

.. math::

   \text{stored } {-J[\text{gas}, y_j]} \mathrel{+}=
     \varphi_p \left(
       \frac{\partial k_c}{\partial r} [\text{A}]_\text{gas}
     - \frac{\partial k_e}{\partial r} \frac{[\text{A}]_\text{aq}}{f_v}
     \right) \frac{\partial r}{\partial y_j}

**Through** :math:`N`:

.. math::

   \text{stored } {-J[\text{gas}, y_j]} \mathrel{+}=
     \varphi_p \left(
       \frac{\partial k_c}{\partial N} [\text{A}]_\text{gas}
     - \frac{\partial k_e}{\partial N} \frac{[\text{A}]_\text{aq}}{f_v}
     \right) \frac{\partial N}{\partial y_j}

**Through** :math:`\varphi_p`:

.. math::

   \text{stored } {-J[\text{gas}, y_j]} \mathrel{+}=
     R \cdot \frac{\partial \varphi_p}{\partial y_j}

where :math:`R = k_c [\text{A}]_\text{gas} - k_e [\text{A}]_\text{aq}/f_v`
is the un-scaled net rate.

Aerosol Property Derivatives
==============================

Effective radius
----------------

- **SingleMomentMode / UniformSection**: :math:`\partial r / \partial y_j = 0`
  (parameterized).
- **TwoMomentMode**:

  .. math::

     \frac{\partial r_\text{eff}}{\partial [\text{species}_{p,i}]} =
       \frac{r_\text{eff} \cdot M_{w,i}}{3 \, V_\text{total} \, \rho_i},
     \qquad
     \frac{\partial r_\text{eff}}{\partial N} =
       -\frac{r_\text{eff}}{3 N}

Number concentration
--------------------

- **SingleMomentMode / UniformSection**: :math:`N = V_\text{total}/V_\text{single}`, so

  .. math::

     \frac{\partial N}{\partial [\text{species}_{p,i}]} =
       \frac{M_{w,i}}{\rho_i \cdot V_\text{single}}

- **TwoMomentMode**: :math:`N` is a prognostic variable, so
  :math:`\partial N / \partial N_\text{var} = 1`.

Phase volume fraction
---------------------

For a mode/section with a single phase, :math:`\varphi_p = 1` and all
partials are zero. For multi-phase:

.. math::

   \frac{\partial \varphi_p}{\partial [\text{species}_{p,i}]} =
     \frac{M_{w,i}}{\rho_i} \cdot \frac{1 - \varphi_p}{V_\text{total}}

.. math::

   \frac{\partial \varphi_p}{\partial [\text{species}_{q,i}]} =
     -\frac{M_{w,i}}{\rho_i} \cdot \frac{\varphi_p}{V_\text{total}}
     \qquad (q \neq p)
