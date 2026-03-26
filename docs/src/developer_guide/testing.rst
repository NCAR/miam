=======
Testing
=======

MIAM uses `Google Test <https://github.com/google/googletest>`_ for all
tests, fetched automatically via CMake FetchContent.

Test Structure
==============

.. code-block:: text

   test/
   ├── CMakeLists.txt           # top-level test configuration
   ├── unit/
   │   ├── CMakeLists.txt
   │   ├── model.cpp            # Model tests
   │   ├── representation_policy.hpp   # shared test helpers
   │   ├── representations/     # representation unit tests
   │   │   ├── test_single_moment_mode.cpp
   │   │   ├── test_two_moment_mode.cpp
   │   │   └── test_uniform_section.cpp
   │   └── processes/           # process unit tests
   │       ├── test_dissolved_reversible_reaction.cpp
   │       ├── test_henry_law_phase_transfer.cpp
   │       └── test_aerosol_property_providers.cpp
   └── integration/
       ├── CMakeLists.txt
       ├── test_dissolved_reversible_reaction.cpp
       └── test_miam_api.cpp

Running Tests
=============

.. code-block:: bash

   cd build
   cmake .. -DMIAM_ENABLE_TESTS=ON
   make -j$(nproc)
   ctest --output-on-failure

To run a specific test:

.. code-block:: bash

   ctest -R test_henry_law_phase_transfer --output-on-failure

Unit Tests
==========

Unit tests verify individual components in isolation. Key patterns:

**Process tests** should verify:

- Forcing function values against known analytical solutions
- Jacobian entries against finite-difference approximations
- State parameter updates (rate constants at different temperatures)
- Sparsity patterns (NonZeroJacobianElements)

**Representation tests** should verify:

- State variable and parameter names
- Provider values (ComputeValue)
- Provider derivatives (ComputeValueAndDerivatives)

**Provider tests** verify the chain rule through aerosol properties. See
[test/unit/processes/test_aerosol_property_providers.cpp](test/unit/processes/test_aerosol_property_providers.cpp)
for examples covering all three representation types.

Integration Tests
=================

Integration tests wire up a complete system with MICM's solver and verify
physical behavior:

- **test_dissolved_reversible_reaction.cpp** — verifies that a reversible
  reaction approaches equilibrium.
- **test_miam_api.cpp** — full multi-mode, multi-process tests covering
  all three representation types with Henry's Law phase transfer.

These tests use MICM's Rosenbrock solver and verify that:

1. Mass is conserved (gas + all condensed instances).
2. Henry's Law equilibrium is approached.
3. The solver converges without excessive iterations.

Test Utilities
==============

``representation_policy.hpp`` provides a shared test header with common
type aliases:

.. code-block:: c++

   #include "representation_policy.hpp"  // provides DenseMatrixPolicy, etc.

Adding a New Test
=================

1. Create the test file in the appropriate directory (``unit/`` or
   ``integration/``).
2. Add it to the corresponding ``CMakeLists.txt``:

   .. code-block:: cmake

      create_standard_test(
          NAME test_my_new_process
          SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test_my_new_process.cpp)

3. Run ``ctest`` to verify the new test passes.
