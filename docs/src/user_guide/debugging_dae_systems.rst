=======================================
Diagnosing Solver Failures in MIAM
=======================================

MIAM uses numerical ODE and DAE solvers (via MICM's Rosenbrock
implementation) to integrate chemical systems forward in time.  When
these solvers fail — or worse, silently produce wrong answers — the
root cause is not always obvious.  This guide covers general strategies
for diagnosing solver problems, followed by a detailed case study from
MIAM's development.

.. contents:: On this page
   :local:
   :depth: 2

Common Causes of Solver Failure
===============================

Before diving into diagnostics, it helps to know the usual suspects.

Incorrect Jacobian
-------------------

If the analytical Jacobian has a bug (wrong derivative, missing term,
wrong sign), the solver's Newton-like iteration converges to the wrong
point or diverges entirely.  Symptoms:

- ``SolverState::Converged`` is never reached, even for tiny time steps
- The solver converges but drifts away from analytical solutions over
  time
- Reducing the time step doesn't help (or makes things worse)

An incorrect Jacobian is one of the most common bugs when implementing
new process types.  Always verify it with finite differences before
debugging anything else.

Stiff systems with too-large time steps
-----------------------------------------

Rosenbrock methods are A-stable, so they handle stiffness well in
principle.  But if the initial time step is too large relative to the
fastest timescale, the first few Newton iterations can diverge before
the step-size controller has a chance to adapt.  Symptoms:

- The solver fails on the very first step
- Reducing the initial ``h_start`` or ``dt`` lets it proceed
- Works fine after the first few steps once the step-size controller
  kicks in

**Fix**: Start with a small time step and ramp up.  The solver's
adaptive step-size controller will find the right step size, but it
needs a reasonable starting point.

Inconsistent DAE initial conditions
-------------------------------------

For DAE systems (models with algebraic constraints), every constraint
equation :math:`G(y) = 0` must be satisfied at the initial state.
If :math:`G(y_0)` is large, the solver's first step includes a huge
Newton correction that overshoots wildly due to nonlinearity.  Symptoms:

- Concentrations jump by many orders of magnitude on the first step
- The solver reports convergence even though values are nonsensical
  (negative concentrations, violated mass balance)
- Each sub-system works in isolation but the combined system fails

.. note::

   **MICM now includes automatic constraint initialization.**  As of
   commit ``9d91673``, the Rosenbrock DAE solver performs Newton
   iteration on the algebraic variables before time-stepping to project
   them onto the constraint manifold.  This means you no longer need to
   compute perfectly consistent initial conditions by hand — a rough
   guess is sufficient.  For systems with very small equilibrium
   constants (e.g., :math:`K_w \approx 10^{-14}`), you may need to
   tighten the initialization parameters:

   .. code-block:: c++

      auto params = RosenbrockSolverParameters::
          FourStageDifferentialAlgebraicRosenbrockParameters();
      params.constraint_init_max_iterations_ = 50;
      params.constraint_init_tolerance_ = 1e-14;

   See the :ref:`case study <case-study-inconsistent-ics>` below for
   the full story.

Negative concentrations
------------------------

Chemical concentrations cannot be negative, but neither Rosenbrock nor
BDF solvers enforce this.  Large, fast reactions can drive a species
below zero within a single solver stage.  Symptoms:

- A species that should asymptotically approach zero goes slightly
  negative
- Subsequent steps amplify the negative value through rate expressions
  (e.g., ``k * [A] * [B]`` where ``[B] < 0`` flips the sign)

**Fix**: Clamp to zero at the physics boundary, or reduce the time
step.  MIAM's ``State::PrintState()`` is useful for spotting this.

Unit-conversion errors
-----------------------

MIAM state variables are in mol/m³, but literature equilibrium
constants are usually in mol/L.  Mixing these up introduces a factor of
1000 per concentration term.  For a dissociation like
A ⇌ B + C, this is a factor of 10⁶ error in the equilibrium residual.
Symptoms:

- Solutions are qualitatively right but quantitatively off by powers
  of 10
- The equilibrium ratio ``[products]/[reactants]`` is consistently wrong
  by the same factor

See the :doc:`../science_guide/henry_law_phase_transfer` guide for
conversion details.

Diagnostic Strategies
=====================

Strategy 1: Verify the Jacobian with finite differences
---------------------------------------------------------

This should always be your **first** diagnostic step.  MIAM and MICM
provide finite-difference Jacobian verification utilities that compare
your analytical Jacobian element-by-element against numerical
derivatives.

For **kinetic process** Jacobians:

.. code-block:: c++

   #include <micm/util/jacobian_verification.hpp>

   auto maps = BuildIndexMaps(model);
   DenseMatrix variables(1, maps.num_variables, 0.0);
   // ... set state point ...

   auto nz = model.NonZeroJacobianElements(maps.variable_indices);
   // ... build sparse matrix from nz elements ...
   auto jac_fn = model.JacobianFunction<DenseMatrix, SparseMatrix>(
       maps.parameter_indices, maps.variable_indices, analytical_jac);
   jac_fn(parameters, variables, analytical_jac);

   auto forcing_fn = model.ForcingFunction<DenseMatrix>(
       maps.parameter_indices, maps.variable_indices);
   auto fd_jac = FiniteDifferenceJacobian<DenseMatrix>(
       [&](const auto& v, auto& f) { forcing_fn(parameters, v, f); },
       variables, num_species);

   auto result = CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrix>(
       analytical_jac, fd_jac, num_species, /*atol=*/1e-5, /*rtol=*/1e-4);
   EXPECT_TRUE(result.passed);

For **constraint** Jacobians, the pattern is similar but uses
``ConstraintResidualFunction`` and ``ConstraintJacobianFunction``.  See
``test_jacobian_verification.cpp`` and ``test_cam_cloud_chemistry.cpp``
for complete working examples.

Also check **sparsity completeness** — a missing entry in the sparsity
pattern silently drops that derivative to zero:

.. code-block:: c++

   auto sparsity = CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrix>(
       analytical_jac, fd_jac, num_species);
   EXPECT_TRUE(sparsity.passed);

Strategy 2: Test sub-systems in isolation
------------------------------------------

When a complex system fails, decompose it into the smallest sub-system
that still reproduces the problem.  This is the single most effective
debugging strategy for chemical mechanisms.

Concretely:

1. Start with **one process** (e.g., one Henry's law phase transfer)
   and verify it against an analytical solution.
2. Add **one constraint** at a time (e.g., mass conservation, then
   equilibrium, then charge balance).
3. At each step, verify convergence and check that the solution makes
   physical sense.
4. When a step fails, the bug is in the most recently added piece or in
   the interaction between the new piece and an existing one.

MIAM's integration tests are structured this way (see
``test_cam_cloud_chemistry.cpp`` Steps 1 through 4).

Strategy 3: Compare kinetic and constrained formulations
---------------------------------------------------------

If you suspect the constraint implementation, compare the DAE
(constrained) system against a fully kinetic ODE system that should
give the same answer at steady state.

For example, a dissociation equilibrium B ⇌ C with :math:`K_{eq} = 5`
can be modeled two ways:

- **Kinetic ODE**: B → C with rate :math:`k_f`, C → B with rate
  :math:`k_r = k_f / K_{eq}`.  Uses the standard Rosenbrock solver.
- **Constrained DAE**: C is algebraic with :math:`[C] = K_{eq} [B]`.
  Uses the DAE Rosenbrock solver with mass-matrix formulation.

At steady state, both systems should produce identical concentrations.
If they disagree, the constraint math is wrong.  If they agree but the
DAE system produces transient garbage, the issue is initial conditions
(see below).

MIAM's ``test_kinetic_vs_constrained.cpp`` implements exactly this
comparison for both dissolved equilibrium and Henry's law systems.

Strategy 4: Print and inspect the state
-----------------------------------------

Sometimes the most effective diagnostic is simply printing the state
before and after each solver step:

.. code-block:: c++

   state.PrintHeader();
   state.PrintState(0);  // before

   auto result = solver.Solve(dt, state);

   state.PrintState(1);  // after

Look for:

- **Negative concentrations** — a sign of overshoot or Jacobian error
- **Mass conservation violation** — sum the species that should be
  conserved and compare before vs. after
- **Unreasonably large values** — indicates a diverging Newton
  iteration
- **NaN or Inf** — typically a division by zero in a rate expression

Strategy 5: Evaluate constraint residuals directly
---------------------------------------------------

For DAE systems, you can compute the constraint residual :math:`G(y)`
at any state point to check whether the algebraic equations are
satisfied:

.. code-block:: c++

   auto residual_fn = model.ConstraintResidualFunction<DenseMatrix>(
       maps.variable_indices);
   DenseMatrix forcing(1, maps.num_variables, 0.0);
   residual_fn(variables, forcing);
   // forcing[0][i_algebraic] should be near zero

If the residuals are large at your initial state, the solver will
struggle on the first step.

.. _case-study-inconsistent-ics:

Case Study: Inconsistent Initial Conditions in CAM Cloud Chemistry
===================================================================

This section walks through a real debugging session from MIAM's
development that motivated the constraint initialization feature in
MICM.

The system
-----------

A sulfur dioxide dissolution mechanism with five constraints:

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Constraint
     - Type
     - Algebraic variable
   * - SO₂(g) ⇌ SO₂(aq)
     - Henry's law
     - SO₂(aq)
   * - SO₂(aq) ⇌ HSO₃⁻ + H⁺
     - Dissociation equilibrium
     - HSO₃⁻
   * - H₂O ⇌ H⁺ + OH⁻
     - Dissociation equilibrium
     - OH⁻
   * - SO₂(g) + SO₂(aq) + HSO₃⁻ = total
     - Mass conservation
     - SO₂(g)
   * - H⁺ = OH⁻ + HSO₃⁻
     - Charge balance
     - H⁺

The symptom
------------

After a single 0.001 s solve step, the solver reported
``SolverState::Converged``, yet the state variables were wildly wrong:

.. code-block:: text

   Before:  SO2_g = 5.0e-13   HSO3- = 3.0e-08   H+ = 1.0e-06
   After:   SO2_g = 6.7e+12   HSO3- = -3.8e+14   H+ = -1.3e+17

Concentrations jumped by 25 orders of magnitude, went negative, and
mass conservation was completely violated (total S went from
3 × 10⁻⁸ to 64).

Critically, each sub-system worked in isolation:

- Henry's law + mass conservation alone → passed
- Water dissociation + charge balance alone → passed
- Combined system → catastrophic failure

Diagnosis
----------

**Step 1: Verify the Jacobian.**  Finite-difference comparison showed
the analytical constraint Jacobian was correct.  This ruled out
derivative bugs.

**Step 2: Evaluate the constraint residuals.**  Computing :math:`G(y_0)`
at the initial conditions revealed the problem:

.. code-block:: python

   # The Kw constraint was off by 4 orders of magnitude
   Hp_times_OHm = Hp * OHm          # = 1.03e-12
   expected      = Kw * S * S        # = 1.00e-08

The initial conditions did not lie on the constraint manifold.

**Step 3: Trace the source.**  The original code computed species
concentrations from a single pH guess without iteration:

.. code-block:: c++

   // WRONG — single-pass, not self-consistent
   double hp_guess = 1e-2;
   double so2_g  = total / (1 + alpha*(1 + Ka1*S/hp_guess));
   double so2_aq = alpha * so2_g;
   double hso3   = Ka1 * so2_aq * S / hp_guess;
   double oh     = Kw * S * S / hp_guess;
   double hp     = oh + hso3;   // overrides hp_guess!

``hp_guess`` (0.01) was used to compute ``oh`` and ``hso3``, then
``hp`` was recalculated from the charge balance as ~1 × 10⁻⁶.  No
constraint equation was satisfied.

Why the solver was fooled
--------------------------

MIAM's Rosenbrock DAE solver forms a linear system at each stage:

.. math::

   (\alpha M - J)\,k = f(y)

where :math:`M` is a diagonal mass matrix (:math:`M_{ii} = 0` for
algebraic variables, :math:`M_{ii} = 1` for differential).  For
algebraic rows this reduces to :math:`-J\,k = G(y)` — essentially one
Newton step.

When :math:`G(y_0)` is enormous, a single Newton step overshoots
because the constraint equations are nonlinear.  The Rosenbrock error
estimator was fooled because it assumes corrections are small relative
to the solution.

The fix (user-side)
--------------------

The immediate fix was a damped fixed-point iteration to compute
self-consistent initial conditions:

.. code-block:: c++

   double hp = initial_guess;
   for (int it = 0; it < 100; ++it)
   {
     oh     = Kw * S * S / hp;
     so2_g  = total_S / (1 + alpha
              + Ka1*alpha*S/hp + Ka2*Ka1*alpha*S*S/(hp*hp));
     so2_aq = alpha * so2_g;
     hso3   = Ka1 * so2_aq * S / hp;
     so3    = Ka2 * hso3 * S / hp;

     double hp_new = oh + hso3 + 2*so3 + 2*so4;
     if (std::abs(hp_new - hp) < 1e-15 * hp) break;
     hp = 0.5 * (hp + hp_new);  // 50% damping
   }

This converges in 2–3 iterations.

The fix (solver-side)
----------------------

While user-side initialization works, it requires domain-specific
knowledge of the chemistry to derive the fixed-point iteration.  To
eliminate this burden, MICM's Rosenbrock solver was updated (commit
``9d91673``) to automatically perform Newton iteration on the algebraic
variables before time-stepping begins.  The solver now:

1. Evaluates the constraint residuals :math:`G(y_0)`
2. Computes the constraint Jacobian :math:`\partial G / \partial y`
3. Solves :math:`-J \, \Delta y = G(y)` to get a Newton correction
4. Applies the correction to algebraic variables only
5. Repeats until :math:`\|G\|_\infty <` tolerance

This means users can now supply rough initial guesses (e.g., all
aqueous species at zero, gas at budget totals) and the solver handles
the projection onto the constraint manifold.  MIAM's integration tests
``Step3b_NaiveInitialConditions`` and ``Step4b_NaiveInitialConditions``
verify this works for the full CAM cloud chemistry system.

High-Level Takeaways
====================

1. **Start with the Jacobian.** If the solver produces nonsense, verify
   the analytical Jacobian against finite differences before anything
   else.  A wrong derivative is the most common implementation bug.

2. **Build incrementally.** Start with one process, one constraint.  Add
   complexity one piece at a time.  When a step fails, the bug is in
   the last piece you added.

3. **Compare formulations.** A kinetic ODE and an equivalent DAE
   constraint should agree at steady state.  Disagreement points to
   constraint math; transient disagreement points to initialization.

4. **The solver can report "Converged" with garbage output.** Always
   check mass conservation, charge balance, and sign of concentrations.

5. **Unit conversions matter enormously.** A factor of 1000 in mol/m³
   vs. mol/L propagates as 10⁶ in a two-product equilibrium expression.
   See :doc:`../science_guide/henry_law_phase_transfer`.

6. **For DAE systems, provide reasonable initial guesses.** MICM's
   constraint initialization will refine them, but starting in the
   right ballpark (correct sign, correct order of magnitude) helps
   convergence.  For systems with very small equilibrium constants,
   consider tightening ``constraint_init_tolerance_``.
