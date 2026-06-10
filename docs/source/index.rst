.. MIAM documentation HTML titles
..
.. # (over and under) for module headings
.. = for sections
.. - for subsections
.. ^ for subsubsections
.. ~ for subsubsubsections
.. " for paragraphs

================================
Welcome to MIAM's documentation!
================================

**Model-Independent Aerosol Module** — a C++ library for representing
aerosol and cloud particle populations and their chemical and physical
transformations, designed to plug into the
`MICM <https://github.com/NCAR/micm>`_ chemistry solver.

MIAM is part of the `MUSICA <https://github.com/NCAR/musica>`_ family of
atmospheric chemistry tools developed at
`NSF NCAR <https://ncar.ucar.edu/>`_.

.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card:: Getting started
        :img-top: _static/index_getting_started.svg
        :link: getting_started/installation
        :link-type: doc

        Check out the getting started guide to build and install MIAM.

    .. grid-item-card:: User guide
        :img-top: _static/index_user_guide.svg
        :link: user_guide/concepts
        :link-type: doc

        Learn how to build aerosol models and configure MIAM here!

    .. grid-item-card:: API reference
        :img-top: _static/index_api.svg
        :link: api_reference/index
        :link-type: doc

        The source code for MIAM is heavily documented. This reference will help you extend MIAM for your needs.

    .. grid-item-card:: Developer guide
        :img-top: _static/index_contribute.svg
        :link: developer_guide/adding_a_process
        :link-type: doc

        If you'd like to contribute new science code or update the docs,
        check out the developer guide!

    .. grid-item-card:: Science guide
        :img-top: _static/index_science.svg
        :link: science_guide/henry_law_phase_transfer
        :link-type: doc

        Dig into the science behind the processes and representations MIAM implements.


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
   user_guide/debugging_dae_systems

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
