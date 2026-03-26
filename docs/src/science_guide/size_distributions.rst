==================
Size Distributions
==================

MIAM supports two families of aerosol size distributions: **modal**
(log-normal) and **sectional** (uniform bins).

Log-Normal (Modal)
==================

A log-normal distribution is parameterized by the geometric mean radius
:math:`r_g` and geometric standard deviation :math:`\sigma_g`:

.. math::

   \frac{dN}{d\ln r} = \frac{N}{\sqrt{2\pi}\,\ln\sigma_g}
   \exp\!\left(-\frac{(\ln r - \ln r_g)^2}{2\ln^2\sigma_g}\right)

Effective radius
----------------

For use in condensation rate calculations, the effective radius accounts
for the surface-area weighting of the distribution:

.. math::

   r_\text{eff} = r_g \cdot \exp(2.5\,\ln^2\sigma_g)

This moment relation applies to both ``SingleMomentMode``
(:math:`r_g` fixed) and ``TwoMomentMode`` (:math:`r_g` derived from
total volume and number).

Single-particle volume
----------------------

For ``SingleMomentMode``, the volume of a single "average" particle:

.. math::

   V_\text{single} = \frac{4}{3}\pi\,r_g^3\,\exp(4.5\,\ln^2\sigma_g)

Number concentration is diagnosed as :math:`N = V_\text{total}/V_\text{single}`.

Two-moment mean radius
-----------------------

For ``TwoMomentMode``, the geometric mean radius is derived from volume
and number:

.. math::

   r_\text{mean} = \left(\frac{3\,V_\text{total}}{4\pi\,N}\right)^{1/3}

.. math::

   r_\text{eff} = r_\text{mean}\cdot\exp(2.5\,\ln^2\sigma_g)

Sectional (UniformSection)
==========================

A sectional bin assumes particles are uniformly distributed between
:math:`r_\min` and :math:`r_\max`:

.. math::

   r_\text{eff} = \frac{r_\min + r_\max}{2}

.. math::

   V_\text{single} = \frac{4}{3}\pi\,r_\text{eff}^3

Number concentration is diagnosed from total volume:
:math:`N = V_\text{total}/V_\text{single}`.

Volume Calculations
===================

All representations compute total volume from species concentrations:

.. math::

   V_\text{total} = \sum_p \sum_i [\text{species}_{p,i}]\cdot\frac{M_{w,i}}{\rho_i}

Species must carry ``molecular weight [kg mol-1]`` and
``density [kg m-3]`` properties for this calculation.
