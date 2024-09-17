---
title: 'DAGGER: A versatile graph library for Topographic Analysis and Landscape Evolution Models'
tags:
  - C++
  - Python
  - R
  - Julia
  - Geomorphology
  - LEM
  - Topography
authors:
  - name: Boris Gailleton
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)

affiliations:
 - name: Université Rennes 1, Rennes, France
   index: 1
 - name: GFZ, Potsdam, Germany
   index: 2

date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Processing numerical topography is a keystone in geomorphological studies. On the first hand, topographic analysis on existing Digital Elevation Models provide _plethora_ of insights on tectonics, climate, nature of the substrate or even on the hazard potentials of a given landscape. On the other hand Lansdcapes Evolution Models (LEMs) simulate the evolution of topography and subsequent hydrologic and sediment fluxes through space and time and across the scales. While the field of applications is extremely varied, they all rely on a common set of base tools and algorithm processing structuring topographic as a graph where elevation dictates the direction of the links connecting nodes together. With this graph structure comes a wide range of efficient algorithms commonly used for processing topography, for example traversing the landscapes from the most upstream to the most downstream nodes or calculating the maximum flow distance from outlets.

However, building and processing topographic graphs is not straightforward and rely on a number of application-dependent user choices. These significantly impact the complexity and efficiency of the subsequent algorithms. Three main aspects complicate the process. (i) Each node can have one or multiple receivers - single flow direction _vs_ multiple flow direction. (ii) Many different grid topologies can link nodes together, _e.g._ regular D8, regular D4 or irregular voronoi grids. Boundary conditions are also a prime issue and manage the existence (or not) of links around the edges of the grid (_e.g._ periodicity, cyclicity) or even within the grids (_e.g._ isolating a single watershed). (iii) Local minimas, whether they be actual data (_e.g._ lakes, endhoreic basins) or noise, break the relationship between absolute elevation and receivers. They complicate the notion of upstream/downstream and the existing solutions imply significantly different implications.

<!-- It is also crucial to consider computational constrains. Some analysis only require single execution of specific algorithms (e.g. filling topographic depression once to get rid of local minimas) while other analysis, like the extraction of river networks, rely on complex successions of algorithms and can benefit from an integrated framework optimised to facilitate repeated operations (_e.g._ getting all upstream nodes of a given one, labelling watersheds, computation of topological sorting). Finally, LEMs require optimised updates of the topographic graph each time the surface evolve, while topographic analysis require only one, more sophisticated built. -->

We introduce `DAGGER`, a library providing a versatile numerical structure specifically designed for surface topographic manipulation. It provides different ways of processing topography that can be optimised for different use cases. It aims to make the access to graph-oriented geomorphological algorithms as generic as possible without sticking to a static topology, boundary condition, local minima method or computational style. The different aspects of the graph computation are decoupled making easy the switch from a specific method to another (e.g. seamlessly changing the local minima solver or the flow direction). Written in template-heavy `C++`, `DAGGER` is inherently designed to be easily wrapped by higher-level languages if needed - we provide extended bindings for `python` and we demonstrate how it can easily be extended to other languages (_e.g._ `Julia`, `R`, `MATLAB` or even `Webassembly`) with minimal code additions. The documentation details the data structure and the algorithms with indications on optimal and pessimal use cases.

# Statement of need

Multiple mature and well-designed frameworks are already dedicated to geomorphological applications. They all offer high-quality graph-based sets of algorithms for geomorphological analysis [TTB, LSD, richDEM], designing LEMs [fastscapelib, Landlab] and others even linking the two [MuddPILE, TTLEM]. Stand-alone codes also provide embedded methods associated with publications [openLEM, CAESAR, CIDRE, cordonnier]. However, they remain tied to their primary applications and/or specific languages. It can be tedious to adapt the algorithms in different contexts without significant code overhead. For example both TTB, LSD and richDEM are designed primarily for regular grids and require significant adaptations when used with different topologies or boundary conditions. LANDLAB provides a very flexible set of grid topologies, but is tightly binded with `python` and `cython` making the use or expansion of this tool difficult from other languages/frameworks or simply for different use cases. Stand-alone codes are by nature designed for a single use and can lack documentation to integrate them for different purposes. Initiatives like `BMI` address this issue by standardising the coding interfaces between models written in different languages, but remain quite specific the modelling world.

The field of computer science also has a lot of general graph libraries in many languages [lemon, boost, ...]. However, they are not dedicated to surface topography and require adaptations to make them ingest and build topographic graphs, with associated performance costs. The entry step is also very high as the technical jargons of each fields are very different.

`DAGGER` aims to fill these gaps by providing a general purpose topography-oriented graph library designed from scratch to be usable from different target languages and for multiple use cases.

# Design and extendability

`DAGGER` is separated into 4 decoupled elements. The `connector` element manages the connectivity between a node and its immediate neighbours. This include management of boundary conditions (e.g. nodes near grid edges; no data; flow can enter, out or not) but also immediate spatial relationships (_e.g._ distance between two nodes, representative area). The `graph` element contain the code related to algorithms dedicated to anything related node connectivity beyond a single node and its immediate vicinity. This includes for example topological sorting, accessing to nodes upstream or downstream of a given one, solving local minimas using [cordonnier]'s method, watersheds labelling,... The third elements are stand-alone algorithms and data-structure taking `graph` and/or `connector` as input and returning independant elements (_e.g._ [Priority flood, lindsay carving, hillshading, ...]). Finally the fourth elements are the `wrap_helper` and manages `I/O` operations depending on the wrapper language. Functions in the three other elements that ingest or return data outside of `DAGGER` are processed by standardised formatting functions ensuring smooth and optimised behind-the-scenes data transfers. Both are defined in a `wrap_helper` dedicated to one wrapping language which converts in place objects from the wrap language (_e.g_ `numpy.ndarray` using `pybind11`) making its values accessible and modifiable through the `operator[]` in `c++`.

Breaking the codes into these elements ensure easy extendability of the code. For example, a new grid topology can be defined by simply creating a corresponding `connector`, or a new language can be wrapped by simply creating the corresponding `wrap_helper`. New stand-alone algorithms or graph functionalities only need to ensure they are written in a `connector`-agnostic way. This design philosophy is conceptually similar to `BMI` and needs well-defined standards for the `connector` and `wrappers` which should be clearly versioned in case of compatibility breaking change.


# Acknowledgements

This work has beneficed from funding from the Helmotz Fundation and from the European Research Council through the project FEASIBLe. I acknowledge Benoit Bovy, Jean Braun, Philippe Steer, Wolfgang Schwanghart and Guillaume Cordonnier for insightful discussions around the topic and Bar Oryan, Emma Graf and Marian Ruiz Sanchez-Oro for testing early version of the software.

# References




<!-- # Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->
