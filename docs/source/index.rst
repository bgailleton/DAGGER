DAGGER's documentation
===================================

**DAGGER** is a ``c++`` header-only topographic graph library specifically designed for geomorphological applications (e.g. topographic analysis, LEMs, surface flow, ...). One of the key aspect of ``DAGGER`` is its modularity: its components are designed from scratch to follow given **standadised interfaces** allowing easy development of alternative module (e.g. replacing a regular grid ``connector`` with an irregular grid ``connector``). We provide and maintain extensive bindings for ``python``, but provide examples on how to use it with, ``Julia``, ``R`` (WIP), ``JS/Wasm`` (WIP) and MATLAB (limited). ``DAG - ger`` stands for **Directed Acyclic Graph**, the base numerical principle to process flow on numerical topography where the *vertices* are the discretised elevations and the *arcs* direction is dictated by topographic gradient. The different modules are:

- the ``connector`` module deals with the grid topology (e.g. location of nodes, area of cells, neighbouring patterns, boundary conditions, receivers/donors computation ...),
- the ``graph`` module connects the nodes beyond immediate vicinity (e.g. topological ordering, local minima resolution, watershed labelling, flow accumulation, multiple vs single flow, ...),
- the ``wrap_helper`` module manages ``I/O`` operations with the targetted platform (e.g. on ``python`` it contains the conversion functions from ``std:vector`` to  ``numpy`` arrays),
- the ``algorithms`` module(s), implementing examples of standalone algorithms taking ``graph`` and/or ``connectors`` as input.

.. note::

   This project is under active development.



Where to go from here ?
-----------------------

- Want to learn by example, see what's available or just try stuff? Go to :ref:`quickstart`
- Want to know the concepts and details behind ``DAGGER``'s structure's? Go to :ref:`modules`
- Want to develop a new module/binding? Go to :ref:`developer`
- Want details about functions, types, classes and to which language they are available? Go to :ref:`api`


Contents
----------

.. toctree::
   
   quickstart
   installation
   modules
   developer
   api

