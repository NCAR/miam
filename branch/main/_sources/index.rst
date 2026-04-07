.. MIAM documentation HTML titles
..
.. # (over and under) for module headings
.. = for sections
.. - for subsections
.. ^ for subsubsections
.. ~ for subsubsubsections
.. " for paragraphs

====
MIAM
====

**Model-Independent Aerosol Module** — a C++ library for representing
aerosol and cloud particle populations and their chemical and physical
transformations, designed to plug into the
`MICM <https://github.com/NCAR/micm>`_ chemistry solver.

MIAM is part of the `MUSICA <https://github.com/NCAR/musica>`_ family of
atmospheric chemistry tools developed at
`NSF NCAR <https://ncar.ucar.edu/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/installation
   getting_started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/concepts
   user_guide/representations
   user_guide/processes
   user_guide/building_a_model
   user_guide/aerosol_properties

.. toctree::
   :maxdepth: 2
   :caption: Science Guide

   science_guide/henry_law_phase_transfer
   science_guide/dissolved_reversible_reaction
   science_guide/size_distributions

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api_reference/index

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide

   developer_guide/adding_a_process
   developer_guide/adding_a_representation
   developer_guide/testing

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
