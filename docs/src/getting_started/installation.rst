============
Installation
============

Requirements
============

- A C++20-capable compiler (GCC 11+, Clang 14+, MSVC 19.30+)
- CMake 3.21 or later

MIAM is a header-only library. Its dependencies —
`MICM <https://github.com/NCAR/micm>`_ and
`Google Test <https://github.com/google/googletest>`_ — are fetched
automatically via CMake's ``FetchContent`` mechanism.

Building from Source
====================

.. code-block:: bash

   git clone https://github.com/NCAR/miam.git
   cd miam
   mkdir build && cd build
   cmake ..
   make -j8

Running Tests
-------------

From the ``build/`` directory:

.. code-block:: bash

   make test

Or with verbose output:

.. code-block:: bash

   ctest --output-on-failure

Uninstalling
------------

.. code-block:: bash

   sudo make uninstall

CMake Options
=============

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Option
     - Default
     - Description
   * - ``MIAM_ENABLE_TESTS``
     - ``ON``
     - Build unit and integration tests
   * - ``MIAM_BUILD_DOCS``
     - ``OFF``
     - Build Sphinx/Doxygen documentation (requires Doxygen and Sphinx)

Docker
======

Build and run the MIAM container:

.. code-block:: bash

   git clone https://github.com/NCAR/miam.git
   cd miam
   docker build -t miam -f docker/Dockerfile .
   docker run -it miam bash

Inside the container:

.. code-block:: bash

   cd /build/
   make test
