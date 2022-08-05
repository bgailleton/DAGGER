---
title: 'DAGGER: A versatile graph library for computing digital topography in geomorphological applications'
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

Numerically manipulating surface topography and bathymetry is a keystone of Geomorphological studies. Topographic analysis on existing Digital Elevation Models provide _plethora_ of insights on the tectonic and climatic history; on the nature of the substrate underlying hillslopes, rivers, beaches or the ocean floor; or even on the hazard potentials of a given landscape. Lansdcapes Evolution Models (LEMs) generate and/or simulate the evolution of topography through space and time applying a wide range of process equations. While the applications, needs, algorithms etc. are extremely varied, they all have a common point. Whether explicitely or implicitely, they process topography as a graph where each discretised location (_i.e._ nodes) is linked to its immediate neighbours. Most applications go further and need a Directed Acyclic Graph (DAG), where the links between nodes have an upstream (_i.e._ donor) and downstream (_i.e_ receiver) direction, leveraging the associated graph theory algorithms. 

Building and processing topographic graph can be far from straightforward. The potential complexity of topographic graph is very application-dependent. Three main points need specific care: (i) the nodes' topology, each node can be linked to one or multiple receivers (with associated differences in the donor direction) - _i.e._ Single Flow _vs_ Multiple flows graphs - and on different kind of grids - _e.g._ regular, irregular, D8, D4; (ii) local minimas (_e.g._ lakes, noise, endhoreic basins) break the relationship between absolute elevation and the notion of upstream/downstream and many different methods exists to bypass this issue with significantly different implications; (iii) boundary conditions determine the links' behaviours or existence with associated geometrical complexities (_e.g._ cyclic boundaries, processing individual watersheds, procedural separation of oceanic and continental domains). To these process-linked specificities, it is also crucial to consider computational constrains. Some analysis only require single execution of specific algorithms (e.g. filling topographic depression once and for all) and optimisation can focus on a single task. Other analysis, like the extraction of river networks, rely on complex successions of algorithms and processes and can benefit from an integrated framework optimised to facilitate repeated operations (_e.g._ getting all upstream nodes of a given one, labelling watersheds, computation of topological sorting). Finally, LEMs require the re-computation of the topographic graph each time the surface evolve, with the need of optimising graph updates.

We introduce `DAGGER`, a library providing a versatile graph structure specifically designed for surface topographic manipulation. It provides different ways of manipulating topography that can be optimised for different use cases from single one-off operations, to detailed and integrated graph structure. It aims to make the access to graph-oriented geomorphologic algorithms and manipulations as generic as possible without sticking to a static topology, boundary condition, local minima method or computational style. It does not aim to be specifically used by one application, but more to provide a versatile base tool to build dedicated frameworks around it. Written in template-heavy `C++`, `DAGGER` is inherently designed to be easily wrapped by higher-level languages if needed - we provide example for `python`, `R`, `Julia` and even `javascript/wasm`. The documentation details the different objects and functions with indications on optimal and pessimal use cases.

# Statement of need

Multiple mature and well-designed frameworks are already dedicated to geomorphological applications. They all offer high-quality graph-based sets of algorithms for general geomorphological analysis [TTB, LSD, richDEM], designing LEMs [fastscapelib, Landlab] and others even linking the two [MuddPILE, TTLEM]. Stand-alone codes also provide embedded methods associated with publications [openLEM, CAESAR, CIDRE, cordonnier]. However, they remain tied to their primary applications and/or specific languages. It can be hard to reuse their algorithms in different contexts (for example from topographic analysi to LEMs) without significant code overhead. For example both TTB, LSD and richDEM are optimised for regular grids and would require significant adaptation to be used with different topologies or boundary conditions. LANDLAB provide a very flexible set of grid topologies, but is tightly binded with `python` and `cython` making the use or expansion of this tool tedious from other languages. Stand-alone codes are by nature designed for a single use and can lack documentation to integrate them for different purposes. Initiatives like `BMI` address this issue by standardising the coding interfaces between models written in different languages, but remain quite specific the the modelling world. 

The field of computer science also comports a lot of general graph libraries in many languages [lemon, boost, ...]. However, they are not dedicated to surface topography which is not the main application of graph theory and require adaptations to make them ingest and build topographic graphs. The entry step is also high as technical jargon is very far from physical geography. 

`DAGGER` aims to fill these gaps by providing a general purpose topographic-oriented graph library designed from scratch to be usable from different target languages and for multiple uses.

# Design and extendability

`DAGGER` is separated into 4 types of elements. The `connector` element manages the connectivity between a node and its immediate neighbours. This include management of boundary conditions (e.g. nodes near grid edges, no data, flow can enter, out or not) but also spatial relationships like distance between two nodes or representative area. The `graph` element contain the code related to algorithms dedicated to graph theory or anything related node connectivity beyond a single node and its immediate vicinity. For example this includes topological sorting, accessing to nodes upstream or downstream of a given one, solving local minimas using [cordonnier]'s method, labelling watersheds, propagating signal upstream/downstream,... `graph` algorithm takes `connector` as arguments of most function, meaning the `graph` is independent from boundary conditions or grid topology. The third elements are stand-alone algorithms and data-structure taking `graph` and/or `connector` as input and returning independant elements (_e.g._ [Priority flood, lindsay carving, hillshading, ...]). Finally the fourth element is the `wrap_helper` and manages `I/O` operation depending on the wrapper language. Each function in the three other elements that need external data are processed by a `format_input` function and returned through a `format_output` function. Both are defined in a `wrap_helper` dedicated to the wrapping language which convert into an object which only require to be accessible though  `operator[]` in `c++`.

Breaking the codes into these elements ensure easy extendability of the code. For example, a new grid topology can be defined by simply creating a corresponding `connector`, or a new language can be wrapped by simply creating the corresponding `wrap_helper`. New stand-alone algorithms or graph functionalities only need to ensure they are written in a `connector`-agnostic way. This design philosophy is conceptually similar to `BMI` and needs well-defined standards for the `connector` and `wrappers` which should be clearly versioned in case of compatibility breaking change.


# Acknowledgements

This work has beneficed from funding from the Helmotz Fundation and from the European Research Council through the project FEASIBLe. I acknowledge Benoit Bovy, Jean Braun and Guillaume Cordonnier for insightful discussions around the topic and Bar Oryan, Emma Graf and Marian Ruiz Sanchez-Oro for testing early version of the software.

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

