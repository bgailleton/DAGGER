DAGGER's documentation
===================================

**DAGGER** is a `c++` header-only topographic graph library specifically designed for geomorphological applications (e.g. topographic analysis, LEMs, surface flow, ...), with bindings available for higher-level languages: `python` (full bindings), `Julia` (good coverage), `R` (WIP), JS/Wasm (WIP) and MATLAB (limited). The goal of this library is to provide a generic and versatile engine for manipulating numerical surface topographic data as a **Directed Acyclic Graph**. The code has several **modules** connected *via* **standardised interfaces** allowing easy extendability:

- the `connector` module deals wiith the grid topology (e.g. location of nodes, area of cells, neighbouring patterns, ...)
- the `graph` module connects the nodes beyond immediate vicinity (e.g. topological ordering, local minima sorting, watershed labelling, receivers/donors computation)
- the `wrap_helper` module manages `I/O` operations with the targetted platform (e.g. on python it contains the conversion functions from `std:vector` to  `numpy` arrays)
- the `algorithms` module(s) which contains standalone functions and classes calling the above modules for specific applications



.. note::

   This project is under active development.

Contents
--------

.. toctree::
   
   quickstart
   installation
   usage
   api
