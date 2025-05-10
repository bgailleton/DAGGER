# DAGGER - Directed Acyclic Graph for digital topoGrahic manipulations Eventually cRoss-platform


`DAGGER` is a library to process digital topography. It is designed as a flexible backend engine for Topographic Analysis and Landscape Evolution Models, providing a range of method to calculate flow routing (Single/Multiple flow, receivers,donors,links,...), resolve local minima (fill, carve, reroute, ...) and a lot of related problems.

Originally coded as a header-only `c++` (most of it can still be used as such), now acts as the `c++` backend of [`scabbard`](https://github.com/bgailleton/scabbard/tree/main), a `python` package for topographic analysis, landscape evolution modelling and hydrodynamics/morphodynamics modelling.

`dagger` development started parallel to `fastescapelib` [see here](https://fastscapelib.readthedocs.io/en/latest/) at the GFZ (Potsdam) as a follow up to [CHONK](https://gmd.copernicus.org/articles/17/71/2024/). `scabbard` and `DAGGER` are now moving to mutualise development efforts and gradually migrate to use `fastscapelib` and `libtopotoolbox` behind the scene. 

## What to do from here?

You want to use `graphflood`, `CHONK` or `trackscape/popscape`? See [`scabbard`](https://github.com/bgailleton/scabbard/tree/main) for using them from `python`. If you are interested by just using the `c++` header-only version, contact me (I'll eventualy add documentation and example, although a lot of the code is being offset to `scabbard` where I develop new tools with `taichi` (GPU), `numba` (JIT-CPU), `libtopotoolbox` and `fastscapelib`).


## Authors

Boris Gailleton - boris.gailleton@univ-rennes.fr




































<!--  -->
