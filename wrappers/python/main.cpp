#include <pybind11/pybind11.h>
#include "D8connector.hpp"
#include "D4connector.hpp"
#include "graph.hpp"
#include "wrap_helper.hpp"

#include "hillshading.hpp"
#include "popscape.hpp"
#include "popscape_utils.hpp"
#include "fastflood.hpp"
#include "graphflood.hpp"
#include "trackscape.hpp"
#include "utils.hpp"
#include "simple_depression_solver.hpp"

using namespace DAGGER;



template<typename CONNECTOR_T>
void declare_graph(py::module &m, std::string typestr)
{

  py::class_<graph<double, CONNECTOR_T > >(m, typestr.c_str(), R"pdoc(
Full Graph module, to plug on a connector to unlock non-local topological operations.

Description:
------------

The graph manages all non local connectivity and ensures the graph is acyclic. 
It contains all the routines for topological ordering, local minima resolving 
(multiple methods), or accumulating values in the upstream/downstream directions


Authors:
--------
B.G.)pdoc")

    .def(py::init<CONNECTOR_T&>())
    .def(
      "init_graph",
       &graph<double,CONNECTOR_T>::init_graph,
       R"pdoc(Initialise the data structure and allocate memory (mostly used internally).)pdoc"
       )
    .def(
      "set_opt_stst_rerouting",
       &graph<double,CONNECTOR_T>:: set_opt_stst_rerouting,
       py::arg("onoff"),
       R"pdoc(Activate (true) or deactivate (false) an optimiser. Most of the time does not make a difference but can _eventually_ approximate a few link a bit more precisely when rerouting local minimas.)pdoc"
       )
    .def(
      "compute_graph",
       &graph<double,CONNECTOR_T>::template compute_graph<py::array_t<double,1>, py::array >,
        py::arg("topography"),py::arg("no_MFD"),py::arg("quicksort_on"),
        R"pdoc(
Full computation of the graph (connector updates of links included).

Description
-------------

COmpute the graph by (i) updating receivers/donors, (ii) resolving local minimas
with the method and (iii) computing topological sortings. Can be set to SFD only
to save time.

Parameters
-----------

- topography (1D array): topographic field (of node size)
- no_MFD: if true, no MFD info are computed (especially the toposort MFD and the
recomputing after local minima solver that can be time consuming 
if done repeteadly - e.g. for LEMs)
- quicksort: if true, the toposort MFD is done by sorting "filled" topography by 
elevation using the std::sort algorithm (aka quicksort), otherwise it uses an 
homemade  topological sorting algorithm. Quicksort is O(nlogn) and toposort 
O(n+l) so the most efficient is case dependent.


returns:
--------

PP_topography (1D array): a preprocessed topography where LM have been filled or 
carved or processed in the wanted way.

Authors:
--------
B.G.

)pdoc"

       )
    .def(
      "is_Sstack_full",
       &graph<double,CONNECTOR_T>:: is_Sstack_full,
       R"pdoc(Debugging function, to ignore)pdoc"
       )
    .def(
      "activate_opti_sparse_border_cordonnier",
       &graph<double,CONNECTOR_T>:: activate_opti_sparse_border_cordonnier,
       R"pdoc(Debugging function, to ignore)pdoc"
       )
    .def(
      "get_all_nodes_upstream_of",
       &graph<double,CONNECTOR_T>::template get_all_nodes_upstream_of< py::array_t<int,1> >,
       py::arg("node"),py::arg("only_SFD"),
       R"pdoc(
Fecth all the nodes upstream of a given one.

Description
-------------

Use the SFD stack structure from Braun and Willett (2013) to gather all the 
nodes upstream of a given one. Effectively labels a watershed from a custom 
starting point. Can be extended to MFD using FIFO queues (slower). This is not
the fastest way to label all the watersheds and should be reserved for one-off
operations.

Parameters
-----------

- node (int): the starting node index
- only_SFD (bool): only label Steepest Descent watersheds (faster, less nodes)

returns:
--------

nodes_upstream (1D array: int): returns a 1D array containing all the node 
indices of upstream locations.

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_all_nodes_downstream_of",
       &graph<double,CONNECTOR_T>::template get_all_nodes_downstream_of< py::array_t<int,1> >,
       py::arg("node"), py::arg("only_SFD"),
       R"pdoc(
Fecth all the nodes downstream of a given one.

Description
-------------

Use the SFD stack structure from Braun and Willett (2013) to gather all the 
nodes downstream of a given one. Effectively labels a water paths from a custom 
starting point. Can be extended to MFD using FIFO queues (slower). This is not
the fastest way to label all the flow paths and should be reserved for one-off
operations.

Parameters
-----------

- node (int): the starting node index
- only_SFD (bool): only label Steepest Descent watersheds (faster, less nodes)

returns:
--------

nodes_downstream (1D array: int): returns a 1D array containing all the node 
indices of downstream locations.

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_SFD_stack",
      &graph<double,CONNECTOR_T>::template get_SFD_stack<py::array_t<size_t,1>>,
      R"pdoc(
Returns the single flow stack (Braun and Willett. 2013) in "stack order".

Description
-------------

Returns the single flow direction stack (sensu Braun and Willett (2013). 
Warning: it corresponds to the last computation (``compute_graph``) and will 
return an empty array if not computed at all yet. 

The stack is in "stack order", **i.e.** from the most dowstream nodes to the 
most upstream ones.

This can also be call the topologicaly sorted list of nodes. be careful though,
the array is of node size but contains ordered node index.

returns:
--------

SFD stack array.

Authors:
--------
B.G.

)pdoc"
      )
    .def(
      "get_MFD_stack",
      &graph<double,CONNECTOR_T>::template get_MFD_stack<py::array_t<size_t,1>>,
      R"pdoc(
Returns the Multiple flow stack in "stack order".

Description
-------------

Returns the multiple flow direction stack.

Warning: it corresponds to the last computation (``compute_graph`` with 
only_SFD = false) and will return an empty array if not computed at all yet. 

The stack is in "stack order", **i.e.** from the most dowstream nodes to the 
most upstream ones.

This can also be call the topologicaly sorted list of nodes. be careful though,
the array is of node size but contains ordered node index. It does not follow 
the same topology than Braun and Willett (2013) and is more difficult to compute
. Only use if MFD is truly needed. 

returns:
--------

MFD stack array.

Authors:
--------
B.G.

)pdoc"
      )

    .def(
      "accumulate_constant_downstream_SFD",
       &graph<double,CONNECTOR_T>::template accumulate_constant_downstream_SFD< py::array_t<double, 1> >,
       py::arg("constant_value"),
       R"pdoc(
Accumulates (integrate) a constant downstream in the SFD direction.

Description
-------------

Sums a constant downstream in the single flow direction (each node adds the 
constant value and give it to its receivers following the inverse stack order).
For example, can be the area of a regular grid cell to compute drainage area.
The operation is very efficient (needs the computed graph).

Parameters
-----------

- constant_value (float): the value to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "accumulate_variable_downstream_SFD",
       &graph<double,CONNECTOR_T>::template accumulate_variable_downstream_SFD< py::array_t<double, 1>, py::array_t<double, 1> >,
       py::arg("values"),
       R"pdoc(
Accumulates (integrate) a variable downstream in the SFD direction.

Description
-------------

Sums a variable downstream in the single flow direction (each node adds the 
variable value and give it to its receivers following the inverse stack order).
For example, can be the cell area times variable precipitation rates of a 
regular grid cell to compute weighted drainage area (~runoff).

The operation is very efficient (needs the computed graph).

Parameters
-----------

- values (1D array of node size, float): the value to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "accumulate_constant_downstream_MFD",
       &graph<double,CONNECTOR_T>::template accumulate_constant_downstream_MFD< py::array_t<double, 1>, py::array_t<double, 1> >,
       py::arg("weights"), py::arg("constant_value"),
       R"pdoc(
Accumulates (integrate) a constant downstream in the MFD direction.

Description
-------------

Sums a constant downstream in the multiple flow direction (each node adds the 
variable value and give it to its receivers following the inverse stack order).

It needs to ingest partitionning weights, an array of n links size telling each
nodes how to partition their value to their multiple receivers. The latter can 
be calculated from a ``Connector`` and is most of the time expressed function of
topographic gradient.

For example, can be the area of a regular grid cell to compute drainage area.

Parameters
-----------

- weigths (1D array of link size, float): the weights
- constant_value (float): the value to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "accumulate_variable_downstream_MFD",
       &graph<double,CONNECTOR_T>::template accumulate_variable_downstream_MFD< py::array_t<double, 1>, py::array_t<double, 1> >,
       py::arg("weights"), py::arg("values"),
       R"pdoc(
Accumulates (integrate) a variable downstream in the MFD direction.

Description
-------------

Sums a variable value of node size downstream in the multiple flow direction 
(each node adds the variable value and give it to its receivers following the 
inverse stack order).

It needs to ingest partitionning weights, an array of n links size telling each
nodes how to partition their value to their multiple receivers. The latter can 
be calculated from a ``Connector`` and is most of the time expressed function of
topographic gradient.

For example, can be the area of a regular grid cell times a variable 
precipitation rate to compute runoff.

Parameters
-----------

- weigths (1D array of link size, float): the weights
- values (1D array of link size, float): the values to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "set_LMR_method",
       &graph<double,CONNECTOR_T>:: set_LMR_method,
       py::arg("LMR_method"),
       R"pdoc(
Sets the Local Minima Resolver method.

Description
-------------

Select which LMR (Local Minima Resolver) to use. It needs to be a LMR enum value
and will dictacte if/how the local minima (internal pits) are rerouted (or not).

Can be one of the following:

- dagger.LMR.cordonnier_fill: approximate filling from Cordonnier et al. (2019)
- dagger.LMR.cordonnier_carve: approximate carving from Cordonnier et al. (2019)
- dagger.LMR.priority_flood: Fill with Barnes et al., 2014
- dagger.LMR.none: Ignore pits

There is no better solution than another, performances as well as accuracy are
extremely case dependent.


Parameters
-----------

- LMR_method (LMR enum)


Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "set_minimum_slope_for_LMR",
       &graph<double,CONNECTOR_T>:: set_minimum_slope_for_LMR,
       py::arg("slope"),
       R"pdoc(LMR solvers impose a numerical topographic gradient to avoid 0 slopes. Default is 1e-5.)pdoc"
       )
    .def(
      "set_slope_randomness_for_LMR",
       &graph<double,CONNECTOR_T>:: set_slope_randomness_for_LMR,
       py::arg("magnitude"),
       R"pdoc(Avoid falt surfaces by imposing a very small randomness when processing local minimas. Must be an order of magitude smaller than the  minimal slope.)pdoc"
       )

    // Distance functions
    .def(
      "get_SFD_distance_from_outlets",
       &graph<double,CONNECTOR_T>::template get_SFD_distance_from_outlets< py::array_t<double,1> >,
       R"pdoc(
Calculates the distance from the outlets following the SFD.

Description
-------------

Uses the SFD stack order to integrate the distance from the outlets to the
sources. Useful for long profiles.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"

       )
    .def(
      "get_SFD_min_distance_from_sources",
       &graph<double,CONNECTOR_T>::template get_SFD_min_distance_from_sources< py::array_t<double,1> >,
       R"pdoc(
Calculates the minimum distance from the sources following the SFD.

Description
-------------

Uses the SFD stack in inverse order to integrate the distance from the sources 
toward the outlets. It keeps the minimum distance. Useful to identify the 
closest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_SFD_max_distance_from_sources",
       &graph<double,CONNECTOR_T>::template get_SFD_max_distance_from_sources< py::array_t<double,1> >,
       R"pdoc(
Calculates the maximum distance from the sources following the SFD.

Description
-------------

Uses the SFD stack in inverse order to integrate the distance from the sources 
toward the outlets. It keeps the maximum distance. Useful to identify the 
furthest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_MFD_max_distance_from_sources",
       &graph<double,CONNECTOR_T>::template get_MFD_max_distance_from_sources< py::array_t<double,1> >,
       R"pdoc(
Calculates the maximum distance from the sources following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the sources 
toward the outlets. It keeps the maximum distance. Useful to identify the 
furthest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_MFD_min_distance_from_sources",
       &graph<double,CONNECTOR_T>::template get_MFD_min_distance_from_sources< py::array_t<double,1> >,
       R"pdoc(
Calculates the minimum distance from the sources following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the sources 
toward the outlets. It keeps the minimum distance. Useful to identify the 
closest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_MFD_max_distance_from_outlets",
       &graph<double,CONNECTOR_T>::template get_MFD_max_distance_from_outlets< py::array_t<double,1> >,
       R"pdoc(
Calculates the maximum distance from the outlet following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the outlet 
toward the outlets. It keeps the maximum distance. Useful to identify the 
furthest outlet.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
       )
    .def(
      "get_MFD_min_distance_from_outlets",
       &graph<double,CONNECTOR_T>::template get_MFD_min_distance_from_outlets< py::array_t<double,1> >,
       R"pdoc(
Calculates the minimum distance from the outlet following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the outlet 
toward the outlets. It keeps the minimum distance. Useful to identify the 
closest outlet.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
       )
    
    // Watershed labelling
    .def(
      "get_SFD_basin_labels",
      &graph<double,CONNECTOR_T>::template get_SFD_basin_labels< py::array_t<int,1> >,
      R"pdoc(
Labels SFD watersheds with unique ID.

Description
-------------

Uses the SFD stack order to label very efficiently basins.

returns:
--------

1D array of integer labels

Authors:
--------
B.G.

)pdoc"
    )


    .def(
      "get_drainage_area_SFD",
      &graph<double,CONNECTOR_T>::template get_drainage_area_SFD< py::array_t<double,1> >,
      R"pdoc(
Labels SFD watersheds with unique ID.

Description
-------------

Uses the SFD stack order to label very efficiently basins.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
    )

    .def(
      "get_drainage_area_MFD",
      &graph<double,CONNECTOR_T>::template get_drainage_area_MFD< py::array_t<double,1>, py::array_t<double,1> >,
      R"pdoc(
Labels MFD watersheds with unique ID.

Description
-------------

Uses the MFD stack order to label very efficiently basins.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"
    )
    

    .def(
      "get_n_pits",
       &graph<double,CONNECTOR_T>:: get_n_pits,
       R"pdoc(Return the number of internal pits (prior solving).)pdoc"
    )

    .def(
      "get_debug_mask",
       &graph<double,CONNECTOR_T>:: get_debug_mask,
       R"pdoc(Ignore. Internal debugging mask process. Changes purpose and is mostly deactivated.)pdoc"
    )

    .def(
      "get_debug_int",
       &graph<double,CONNECTOR_T>:: get_debug_int,
       R"pdoc(Ignore. Internal debugging int process. Changes purpose and is mostly deactivated.)pdoc"
    )
    
  ;
}

template<typename CONNECTOR_T>
void declare_popscape_old(py::module &m, std::string typestr)
{
  py::class_<popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T > >(m, typestr.c_str())
    .def(py::init<RANDNOISE,int,int,double,double>())
    // .def_readwrite("graph",  &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::graph)
    // .def_readwrite("connector",  &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::connector)
    .def("solve_generic", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::solve_generic)
    .def("get_topo", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_topo<py::array>)
    .def("get_QA", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_QA<py::array>)
    .def("compute_graph", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::compute_graph)
    .def("compute_DA_SFD", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::compute_DA_SFD)
    .def("apply_uplift", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::apply_uplift)
    .def("apply_variable_uplift", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template apply_variable_uplift<py::array_t<double,1> >)
    .def("solve_SFD_SPL_imp", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::solve_SFD_SPL_imp)
    .def("hydraulic_erosion_v0", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::hydraulic_erosion_v0)
    .def("normalise_topography", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::normalise_topography)
    // .def("run_SFD_exp_latmag", &popscape_old<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::run_SFD_exp_latmag)
    
  ;
}

template<typename CONNECTOR_T>
void declare_popscape(py::module &m, std::string typestr)
{
  py::class_<popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T > >(m, typestr.c_str())
    .def(py::init<DAGGER::graph<double, CONNECTOR_T>&, CONNECTOR_T&>())
    .def("StSt", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::StSt)
    .def("restriction", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::restriction)
    .def("interpolation", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::interpolation)
    .def("smooth", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::smooth)
    .def("set_topo", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_topo<py::array_t<double,1> >)
    .def("get_topo", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_topo<py::array_t<double,1> >)
    .def("get_QA", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_QA<py::array_t<double,1> >)
    .def("get_chistar", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_chistar<py::array_t<double,1> >)
    .def("simple_Kfchi", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::simple_Kfchi)
    .def("simple_Kfz", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::simple_Kfz)
  ;
}

template<typename CONNECTOR_T>
void declare_trackscape(py::module &m, std::string typestr)
{
  py::class_<trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T > >(m, typestr.c_str())
    .def(py::init<>())
    .def_readwrite("graph", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::graph)
    .def_readwrite("connector", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::connector)
    .def("init_random", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::init_random)
    .def("get_topo", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_topo<py::array>)
    .def("get_hillshade", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_hillshade<py::array>)
    .def("get_h_sed", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_h_sed<py::array>)
    .def("get_Qw", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_Qw<py::array>)
    .def("get_precipitations", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_precipitations<py::array>)
    .def("run_SFD", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::run_SFD)
    .def("block_uplift", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::block_uplift)
    .def("external_uplift", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template external_uplift<py::array_t<double,1>& >)    
    .def("init_TSP_module", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template init_TSP_module<py::array_t<double,1>& >)
    .def("update_TSP_source",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template update_TSP_source<py::array_t<double,1>& >)
    .def("sample_carrot_TSP", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template sample_carrot_TSP<py::array_t<double,1> >)    
    .def("sample_carrot_Ch_MTSI", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template sample_carrot_Ch_MTSI<py::array_t<double,1> >)
    .def("get_transect_Ch_MTSI", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_transect_Ch_MTSI<py::array_t<double,1> >)
    .def("get_transect_TSP", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_transect_TSP<py::array_t<double,1> >)
    .def("set_single_Ks", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_Ks)
    .def("set_single_Kr", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_Kr)
    .def("set_single_depcoeff", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_depcoeff)
    .def("set_single_precipitations", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_precipitations)
    .def("set_single_kappa_s", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_kappa_s)
    .def("set_single_kappa_r", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_kappa_r)
    .def("set_single_Sc", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_Sc)
    .def("set_single_Sc_M", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_Sc_M)
    .def("set_single_lambda", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_lambda)
    .def("set_single_sea_level", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_single_sea_level)
    .def("hillslopes_on", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::hillslopes_on)
    .def("hillslopes_off", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::hillslopes_off)
    .def("fluvial_on", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::fluvial_on)
    .def("fluvial_off", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::fluvial_off)
    .def("marine_on", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::marine_on)
    .def("marine_off", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::marine_off)
    .def("fill_up", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::fill_up)
    .def("init_Ch_MTSI", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::init_Ch_MTSI)
    .def("rise_boundary_by", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::rise_boundary_by)
    .def("get_TSP_surface_concentrations", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_TSP_surface_concentrations<py::array>)
    .def("get_Ch_MTIS_surface_age", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_Ch_MTIS_surface_age<py::array>)
    .def("set_variable_precipitations", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_precipitations<py::array_t<double,1>& >)    
    .def("set_variable_Kr", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_Kr<py::array_t<double,1>& >)    
    .def("set_variable_Ks",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_Ks<py::array_t<double,1>& >)
    .def("set_variable_depcoeff",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_depcoeff<py::array_t<double,1>& >)
    .def("set_variable_kappa_s",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_kappa_s<py::array_t<double,1>& >)
    .def("set_variable_kappa_r",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_kappa_r<py::array_t<double,1>& >)
    .def("set_variable_Sc",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_Sc<py::array_t<double,1>& >)
    .def("set_variable_Sc_M",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_Sc_M<py::array_t<double,1>& >)
    .def("set_variable_Ke",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_Ke<py::array_t<double,1>& >)
    .def("set_variable_lambda",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_lambda<py::array_t<double,1>& >)
    .def("set_variable_sea_level",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template set_variable_sea_level<py::array_t<double,1>& >)
    .def("feed_topo",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template feed_topo<py::array_t<double,1>& >)
    .def("set_m",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_m)
    .def("set_n",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::set_n)
    .def("run_SFD_implicit",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::run_SFD_implicit)
    .def("lithify",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::lithify)
    .def("strip_sediment",&trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::strip_sediment)

    
    

    
    
  ;
}

template<typename CONNECTOR_T>
void declare_ff(py::module &m, std::string typestr)
{
  py::class_<fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> > >(m, typestr.c_str())
    .def(py::init<DAGGER::graph<double, CONNECTOR_T>&, CONNECTOR_T&, py::array_t<double,1>&,py::array_t<double,1>& >())
    .def_readwrite("rec", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::rec)
    // .def("run_SFD", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_SFD)
    // .def("run_SFD_with_erosion", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_SFD_with_erosion)
    // .def("run_MFD_erosion", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD_erosion)
    // .def("run_MFD_erosion_B", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD_erosion_B)
    .def("run_MFD_static", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD_static)
    .def("run_MFD_static_SPL", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD_static_SPL)
    // .def("run_MFD", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD)
    // .def("run_MFD_dynamic", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD_dynamic)
    // .def("run_MFD_exp", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::run_MFD_exp)
    .def("get_hw", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_hw<py::array >)
    // .def("get_spatial_dts", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_spatial_dts<py::array >)
    .def("get_Qwin", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_Qwin<py::array >)
    .def("get_Qwout", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_Qwout<py::array >)
    .def("get_Qs", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_Qs<py::array >)
    .def("get_topography", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_topography<py::array >)
    .def("add_to_hw", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::add_to_hw)
    .def("set_Qbase", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template set_Qbase<py::array_t<double,1> >)
    .def("set_Qs_entry_points", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template set_Qs_entry_points<py::array_t<double,1>,py::array_t<int,1> >)
    .def("increment_hw_from_Qbase", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::increment_hw_from_Qbase    )
    .def("caesar_lisflood", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::caesar_lisflood    )
    .def("set_topological_number", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_topological_number    )
    // .def("basicFloodos", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::basicFloodos)
    // .def("basicFloodos_v2", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::basicFloodos_v2)
    // .def("basicFloodos_v3", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::basicFloodos_v3)
    // .def("fill_up", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::fill_up)
    .def("set_manning", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_mannings)
    // .def("testDebugWalk", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::testDebugWalk)
    .def("set_parting_coeff", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_parting_coeff)
    // .def("check_SD_val", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::check_SD_val)
    .def("set_out_boundaries_to_permissive", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_out_boundaries_to_permissive)
    // .def("set_edges_to_0", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_edges_to_0)
    .def("get_a_eff", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_a_eff<py::array >)
    .def("get_w_eff", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_w_eff<py::array >)
    .def("get_hydraulic_slope_D8", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template get_hydraulic_slope_D8<py::array >)
    // .def("spatial_dt",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template spatial_dt<py::array_t<double,1> >)
    .def("set_dt",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_dt)
    // .def("enable_Afdt",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::enable_Afdt)
    // .def("disable_Afdt",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::disable_Afdt)
    // .def("config_Afdt",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::config_Afdt)
    .def("enable_hflow",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::enable_hflow)
    .def("disable_hflow",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::disable_hflow)
    .def("set_sensibility_to_flowdepth",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_sensibility_to_flowdepth)
    .def("get_sensibility_to_flowdepth",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::get_sensibility_to_flowdepth)
    // .def("fill_topo",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::fill_topo)
    .def("set_stochaslope",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_stochaslope)
    .def("out_boundary_match_donors",&fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::out_boundary_match_donors)
    .def("set_boundary_slope", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::set_boundary_slope)







#ifdef OPENMP_YOLO  
    .def("check_devices", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template check_devices, py::call_guard<py::gil_scoped_release>()    )
    .def("caesar_lisflood_OMP", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,  DAGGER::numvec<double> >::template caesar_lisflood_OMP, py::call_guard<py::gil_scoped_release>()    )
#endif

  ;
}


template<typename fT, typename GRAPH_T, typename CONNECTOR_T>
void declare_graphflood(py::module &m, std::string typestr)
{
  py::class_< graphflood<fT, GRAPH_T, CONNECTOR_T> >(m, typestr.c_str())
    .def(py::init<GRAPH_T&, CONNECTOR_T&>())
    .def("run", &graphflood<fT, GRAPH_T, CONNECTOR_T>::run, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("set_topo", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template set_topo <py::array_t<double,1> >, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("set_hw", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template set_hw <py::array_t<double,1> >, R"pdoc(Main function running the model from all the input params)pdoc")
    
    .def("get_hw", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_hw <py::array_t<double,1> >, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("get_surface_topo", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_surface_topo <py::array_t<double,1> >, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("get_bedrock_topo", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_bedrock_topo <py::array_t<double,1> >, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("get_Qwin", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Qwin <py::array_t<double,1> >, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("get_SSTACKDEBUG", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_SSTACKDEBUG <py::array_t<size_t,1> >, R"pdoc(Main function running the model from all the input params)pdoc")

    .def("enable_MFD", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_MFD, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("enable_SFD", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_SFD, R"pdoc(Main function running the model from all the input params)pdoc")
    .def("set_dt_hydro",&graphflood<fT, GRAPH_T, CONNECTOR_T>::set_dt_hydro)
    .def("fill_minima",&graphflood<fT, GRAPH_T, CONNECTOR_T>::fill_minima)
    .def("reroute_minima",&graphflood<fT, GRAPH_T, CONNECTOR_T>::reroute_minima)
    .def("ignore_minima",&graphflood<fT, GRAPH_T, CONNECTOR_T>::ignore_minima)
    .def("enable_morpho",&graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_morpho)
    .def("disable_morpho",&graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_morpho)
    .def("set_dt_morpho_multiplier",&graphflood<fT, GRAPH_T, CONNECTOR_T>::set_dt_morpho_multiplier)


    .def("set_dt_morpho", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_dt_morpho)
    .def("set_single_aexp", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_single_aexp)
    .def("set_single_ke", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_single_ke)
    .def("set_single_ke_lateral", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_single_ke_lateral)
    .def("set_single_kd", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_single_kd)
    .def("set_single_kd_lateral", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_single_kd_lateral)
    .def("set_single_tau_c", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_single_kd_lateral)
    .def("set_variable_ke", &graphflood<fT,GRAPH_T,CONNECTOR_T>::template set_variable_ke<py::array_t<double,1> >)
    


    .def("compute_tuqQ", &graphflood<fT,GRAPH_T,CONNECTOR_T>::template compute_tuqQ<py::array_t<double,1> >)
    .def("compute_elemental_transfer", &graphflood<fT,GRAPH_T,CONNECTOR_T>::template compute_elemental_transfer<py::array_t<double,1>, py::array_t<double,1> >)




    .def("set_water_input_by_entry_points", &graphflood<fT,GRAPH_T,CONNECTOR_T>::template set_water_input_by_entry_points<py::array_t<double,1> ,py::array_t<int,1> >)
    .def("set_water_input_by_constant_precipitation_rate", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_water_input_by_constant_precipitation_rate)
    .def("set_water_input_by_variable_precipitation_rate", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template set_water_input_by_variable_precipitation_rate<py::array_t<double,1> >)
    
    .def("set_sed_input_by_entry_points", &graphflood<fT,GRAPH_T,CONNECTOR_T>::template set_sed_input_by_entry_points<py::array_t<double,1> ,py::array_t<int,1> >)

    .def("enable_Qwout_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_Qwout_recording)
    .def("disable_Qwout_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_Qwout_recording)
    .def("get_Qwout_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Qwout_recording<py::array_t<double,1> >)

    .def("enable_Sw_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_Sw_recording)
    .def("disable_Sw_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_Sw_recording)
    .def("get_Sw_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Sw_recording<py::array_t<double,1> >)

    .def("enable_dhw_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_dhw_recording)
    .def("disable_dhw_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_dhw_recording)
    .def("get_dhw_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_dhw_recording<py::array_t<double,1> >)

    .def("enable_filling_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_filling_recording)
    .def("disable_filling_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_filling_recording)
    .def("get_filling_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_filling_recording<py::array_t<double,1> >)

    .def("enable_edot_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_edot_recording)
    .def("disable_edot_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_edot_recording)
    .def("get_edot_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_edot_recording<py::array_t<double,1> >)

    .def("enable_flowvec_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_flowvec_recording)
    .def("disable_flowvec_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_flowvec_recording)
    .def("get_flowvec_recording", &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_flowvec_recording<py::array_t<double,1> >)

  
    .def("get_tot_Qw_input", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qw_input)
    .def("get_tot_Qw_output", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qw_output)
    .def("get_tot_Qwin_output", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qwin_output)
    .def("get_tot_Qs_output", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qs_output)


    .def("set_stochaslope", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_stochaslope)
    .def("disable_stochaslope", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_stochaslope)
    .def("set_fixed_hw_at_boundaries", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_fixed_hw_at_boundaries)
    .def("set_fixed_slope_at_boundaries", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_fixed_slope_at_boundaries)
    .def("get_dt_hydro",&graphflood<fT, GRAPH_T, CONNECTOR_T>::get_dt_hydro)

    .def("set_partition_method", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_partition_method)
    .def("set_topological_number", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_topological_number)
    .def("get_topological_number", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_topological_number)


    .def("get_courant_dt_hydro", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_courant_dt_hydro)
    .def("set_courant_numer", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_courant_numer)
    .def("set_max_courant_dt_hydro", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_max_courant_dt_hydro)
    .def("set_min_courant_dt_hydro", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_min_courant_dt_hydro)
    .def("enable_courant_dt_hydro", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_courant_dt_hydro)
    .def("set_Qwin_crit", &graphflood<fT,GRAPH_T,CONNECTOR_T>::set_Qwin_crit)
    .def("get_nT", &graphflood<fT,GRAPH_T,CONNECTOR_T>::get_nT)


    .def("enable_hydrostationnary", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_hydrostationnary)
    .def("disable_hydrostationnary", &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_hydrostationnary)

    
    

    .def("block_uplift", &graphflood<fT,GRAPH_T,CONNECTOR_T>::block_uplift)
    .def("variable_uplift", &graphflood<fT,GRAPH_T,CONNECTOR_T>::template variable_uplift<py::array_t<double,1> >)

    .def("run_precipitions",  &graphflood<fT,GRAPH_T,CONNECTOR_T>::run_precipitions)
    .def("run_precipitions_exp",  &graphflood<fT,GRAPH_T,CONNECTOR_T>::run_precipitions_exp)
    .def("run_graphipiton",  &graphflood<fT,GRAPH_T,CONNECTOR_T>::run_graphipiton)
    .def("run_exp", &graphflood<fT,GRAPH_T,CONNECTOR_T>::run_exp)

    .def("define_precipitations_Ath", &graphflood<fT,GRAPH_T,CONNECTOR_T>::define_precipitations_Ath)

        
    

  ;
}



PYBIND11_MODULE(dagger, m) {
  m.doc() = R"pbdoc(
      DAGGER - python API
      ===================
      
      Quick API
      ---------

      .. autosummary::

          graph
          graph.init_graph
          graph.set_opt_stst_rerouting
          graph.compute_graph
          graph.is_Sstack_full
          graph.activate_opti_sparse_border_cordonnier
          graph.get_all_nodes_upstream_of
          graph.get_all_nodes_downstream_of
          graph.get_SFD_stack
          graph.get_MFD_stack
          graph.accumulate_constant_downstream_SFD
          graph.accumulate_variable_downstream_SFD
          graph.accumulate_constant_downstream_MFD
          graph.accumulate_variable_downstream_MFD
          graph.set_LMR_method
          graph.set_minimum_slope_for_LMR
          graph.set_slope_randomness_for_LMR
          graph.get_SFD_distance_from_outlets
          graph.get_SFD_min_distance_from_sources
          graph.get_SFD_max_distance_from_sources
          graph.get_MFD_max_distance_from_sources
          graph.get_MFD_min_distance_from_sources
          graph.get_MFD_max_distance_from_outlets
          graph.get_MFD_min_distance_from_outlets
          graph.get_SFD_basin_labels

          D8N
          D8N.__init__
          D8N.set_default_boundaries
          D8N.set_custom_boundaries
          D8N.print_dim
          D8N.get_HS
          D8N.get_mask_array
          D8N.set_values_at_boundaries
          D8N.set_out_boundaries_to_permissive
          D8N.get_boundary_at_node
          D8N.get_rowcol_Sreceivers
          D8N.print_receivers
          D8N.get_rec_array_size
          D8N.update_links_MFD_only
          D8N.update_links
          D8N.update_links_from_topo
          D8N.sum_at_outlets
          D8N.keep_only_at_outlets
          D8N.get_SFD_receivers
          D8N.get_SFD_dx
          D8N.get_SFD_ndonors
          D8N.get_SFD_donors_flat
          D8N.get_SFD_donors_list
          D8N.get_links
          D8N.get_linknodes_flat
          D8N.get_linknodes_list
          D8N.get_linknodes_list_oriented
          D8N.get_SFD_receivers_at_node
          D8N.get_SFD_dx_at_node
          D8N.get_SFD_ndonors_at_node
          D8N.get_SFD_donors_at_node
          D8N.get_SFD_gradient
          D8N.get_links_gradient
          D8N.get_MFD_mean_gradient
          D8N.get_MFD_weighted_gradient
          D8N.get_link_weights
          D8N.set_stochaticiy_for_SFD


      Full API
      ---------

      .. autoclass:: D8N
        :members:

      .. autoclass:: graph
        :members:

          
  )pbdoc";


  py::enum_<DEPRES>(m, "LMR")
    .value("cordonnier_fill", DEPRES::cordonnier_fill)
    .value("cordonnier_carve", DEPRES::cordonnier_carve)
    .value("cordonnier_simple", DEPRES::cordonnier_simple)
    .value("priority_flood", DEPRES::priority_flood)
    .value("priority_flood_opti", DEPRES::priority_flood_opti)
    .value("priority_full_MFD", DEPRES::priority_full_MFD)
    
    .value("none", DEPRES::none)
    .value("dagger_carve", DEPRES::dagger_carve)
    .value("dagger_fill", DEPRES::dagger_fill)
  ;

  py::enum_<MFD_PARTITIONNING>(m, "MFD_PARTITIONNING")
    .value("PROPOSLOPE", MFD_PARTITIONNING::PROPOSLOPE)
    .value("SQRTSLOPE", MFD_PARTITIONNING::SQRTSLOPE)
    .value("PROPOREC", MFD_PARTITIONNING::PROPOREC)
  ;




  py::class_<D8connector<double> >(
m, "D8N", R"pbdoc(
D8 regular grid connector D8N

Description:
------------

D8N class: the built-in connector for regular grids. D8N regular grids are the most used type of grids, like most global DEMs for example.
They are defined by a number of rows (ny) and a number of columns (nx) with associated spacing between nodes in both directions (dy,dx).
The D8N constructor defaults to ``4edges`` boundary preset and precompiles neighbouring. Neighbouring and indexing is optimised by geometrical relationships:
Node index is simply a flat index in the row major direction (idx = row_id * nx + col_id) and the link ID for the top-right, right, bottomright and bottom neighbours are respectiveley idx *4, idx*4 +1, idx * 4+2 and idx*4 + 3

The constructor only needs basics geometrical information.

Parameters:
-----------

    * nx (int): number of nodes in the X direction (columns)
    * ny (int): number of nodes in the Y direction (rows)
    * dx (float64): distance between nodes in the x directions
    * dy (float64): distance between nodes in the y directions
    * x_min (float64): X coordinates of the bottom left corner
    * y_min (float64): Y coordinates of the top left corner


Authors:
--------
B.G.

)pbdoc")

    .def(
      py::init<int,int,double,double,double,double>(),
      py::arg("nx"),py::arg("ny"),py::arg("dx"),py::arg("dy"),py::arg("x_min"), py::arg("x_max"),
      R"pbdoc(
D8 regular grid connector D8N

Description:
------------

D8N class: the built-in connector for regular grids. D8N regular grids are the most used type of grids, like most global DEMs for example.
They are defined by a number of rows (ny) and a number of columns (nx) with associated spacing between nodes in both directions (dy,dx).
The D8N constructor defaults to ``4edges`` boundary preset and precompiles neighbouring. Neighbouring and indexing is optimised by geometrical relationships:
Node index is simply a flat index in the row major direction (idx = row_id * nx + col_id) and the link ID for the top-right, right, bottomright and bottom neighbours are respectiveley idx *4, idx*4 +1, idx * 4+2 and idx*4 + 3

The constructor only needs basics geometrical information.

Parameters:
-----------

    * nx (int): number of nodes in the X direction (columns)
    * ny (int): number of nodes in the Y direction (rows)
    * dx (float64): distance between nodes in the x directions
    * dy (float64): distance between nodes in the y directions
    * x_min (float64): X coordinates of the bottom left corner
    * y_min (float64): Y coordinates of the top left corner


Authors:
--------
B.G.

)pbdoc"
      )
    .def(
      "set_default_boundaries", &D8connector<double>::set_default_boundaries,
      py::arg("boundary_preset"),
      R"pbdoc(

Automatically sets the default boundary system based on predefined presets.

Description:
------------

Boundary conditions can be tricky to set. Luckily, most use cases needs 
classical opened bounary at the edge of the dem, or periodic (sometimes called 
cyclic) bounaries at EW or NS edges. This functions automate these (i.e. set the
right boundary codes and recompute the linknode array) 

Parameters:
-----------

    * boundary_preset (str): the preset. Can be "4edges", "periodic_EW" or 
      "periodic_NS"


Authors:
--------
B.G.

)pbdoc"
      )

    .def(
      "set_custom_boundaries", 
      &D8connector<double>::set_custom_boundaries<py::array_t<int,1> >,
      py::arg("boundary_codes"),
      R"pbdoc(

Manually sets the boundary codes as uint8 1D array (see doc for options).

Description:
------------

Function to manually specify boundary conditions, for cases where the classic 
presets are not providing the desired option. For example if one needs specific
outlets/inlets location or deactivate some nodes. This is a good function to 
automate the ingestion of nodata, isolating a watershed, or defining entry
points to the landscape.

Parameters:
-----------

    * boundary_codes (1D Array): the uint8 array of node size containing the
      boundary codes (see section :ref:`boundary`)


Authors:
--------
B.G.

)pbdoc"


      )
    .def("print_dim", &D8connector<double>::print_dim, R"pdoc(Debugging function)pdoc")
    .def("get_HS", &D8connector<double>::get_HS<std::vector<double>, py::array >, R"pdoc(Deprecated, kept for legacy (see hillshade function outside of ``connector``))pdoc")
    .def("get_mask_array",&D8connector<double>::get_mask_array, R"pdoc(Returns a 1D array of bool where false are nodata)pdoc")

    .def(
      "set_values_at_boundaries", 
      &D8connector<double>::set_values_at_boundaries<py::array_t<double,1> >,
      py::arg("array_to_modify"), py::arg("value"),
      R"pdoc(
Modify an array in place with given value where connector can out flux

Description:
------------

Every nodes on the input array where the connector satisfy the coundary condition 
``can_out`` will be set to a given value. It modify the array in place (i.e. 
does not return anything but modify the input).

Useful to impose boundary condition in LEMs in a versatile way.

Parameters:
-----------

    * array_to_modify (1D Array): flot64 array to be modified in place
    * value (float64): value to impose


Authors:
--------
B.G.

)pdoc"
      )

    .def(
      "set_out_boundaries_to_permissive",
      &D8connector<double>::set_out_boundaries_to_permissive,
      R"pdoc(
Converts OUT boundaries to CAN_OUT.

Description:
------------

Converts OUT boundaries, where flow entering this cell has to out the model to 
CAN_OUT, where flow leaves the model if the cell has no downslope receivers.

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_boundary_at_node", &D8connector<double>::get_boundary_at_node,py::arg("node_idx"), R"pdoc(Returns the boundary code at node index)pdoc"
    )
    
    .def("get_rowcol_Sreceivers",&D8connector<double>:: get_rowcol_Sreceivers,py::arg("row_index"), py::arg("col_index"), R"pdoc(Debug function to get the receiver (node) indices of a node from its row and column index)pdoc")
    
    .def("print_receivers", &D8connector<double>::template print_receivers<std::vector<double> >,py::arg("node_index"), py::arg("topography"), R"pdoc(Debuggin function printing to the terminal the receivers of a node index and their topography (post graph computation! so the topographic field may not be the one used for the receivers/LM computations))pdoc")
    
    .def("get_rec_array_size",&D8connector<double>::get_rec_array_size, R"pdoc(Debug function - ignore)pdoc")
    
    .def(
      "update_links_MFD_only", &D8connector<double>::template update_links_MFD_only<std::vector<double> >,
      py::arg("topography"),
      R"pdoc(
Updates all the link directionalities - but not the SFD receiver/donors.

Description:
------------

Updates all the link directionalities based on a given topography. Note that it
does not process the SFD receiver/donors and is (paradoxally) faster. However,
the use of this function is reserved to experienced user who seek tuning perform
ances as many routines are speed up by SFD info behind the scene -  even in MFD.

Parameters
----------

- topography (1D array): the flat topography


Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "update_links_from_topo", 
      &D8connector<double>::template update_links_from_topo<std::vector<double> >,
      py::arg("topography"),
      R"pdoc(
Updates all the link directionalities - but not the SFD receiver/donors.

Description:
------------

Updates all the link directionalities based on a given topography. This is the 
full updating function computing MFD/SFD links/receivers/donors/... .

Parameters
----------

- topography (1D array): the flat topography


Authors:
--------
B.G.

)pdoc"
    )
    
    .def(
      "sum_at_outlets", 
      &D8connector<double>::template sum_at_outlets<py::array_t<double,1>, double >,
      py::arg("array"),py::arg("include_pits"),
      R"pdoc(
Sum the values contains in the input array where flux out the model.

Description:
------------

Sum the values contains in the input array where flux out the model. It can also
include the internal pit if needed. Internal pits are only referring to nodes
where local minimas have not been resolved

DEPRECATION WARNING: Will be detached to a standalone algorithm in a future
update.

Parameters
----------

- array (1D array): the array to sum (of node size)
- include_pits (bool): if true, internal pits (receiverless nodes) are sumed too

Returns
----------

The sum.

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "keep_only_at_outlets", 
      &D8connector<double>::template keep_only_at_outlets<py::array_t<double,1>, py::array >,
      py::arg("array"),py::arg("include_pits"),
      R"pdoc(
return a copy of the input array where all the non-outting nodes are set to 0.

Description:
------------

returns a copy of the input array where all the outting nodes are set to 0. 
It can also include the internal pit if needed. Internal pits are only referring 
to nodes where local minimas have not been resolved.

DEPRECATION WARNING: Will be detached to a standalone algorithm in a future
update.

Parameters
----------

- array (1D array): the array to sum (of node size)
- include_pits (bool): if true, internal pits (receiverless nodes) are sumed too

Returns
----------

The modified array

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_SFD_receivers",
      &D8connector<double>::template get_SFD_receivers<py::array_t<int,1>>,
      R"pdoc(
returns the array of SFD receivers.

Description:
------------

returns the array of SFD receivers according to the current state of the \
connector. These can be obtained after computing links from a topography, but 
can also be modified when a graph resolves local minimas.

Returns
----------

array of node size containing the node index of the SFD receivers

Authors:
--------
B.G.

)pdoc"
    )
    
    .def(
      "get_SFD_dx",
      &D8connector<double>::template get_SFD_dx<py::array_t<double,1>>,
      R"pdoc(
returns the array of SFD distance to receivers.

Description:
------------

returns the array of SFD dist. to receivers according to the current state of
the connector. These can be obtained after computing links from a topography,  
but can also be modified when a graph resolves local minimas.

Returns
----------

array of floating point distance to the SFD receiver

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_SFD_ndonors",
      &D8connector<double>::template get_SFD_ndonors<py::array_t<int,1>>,
      R"pdoc(
returns the array of SFD number of donors.

Description:
------------

returns the array of number of SFD donors for every nodes.

Returns
----------

array of integer of node size with the number of SFD donors for each nodes

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_SFD_donors_flat",
      &D8connector<double>::template get_SFD_donors_flat<py::array_t<int,1>>,
      R"pdoc(
returns a flat array of SFD donors(read description for indexing!).

Description:
------------

returns a flat array of SFD donors. It is a sparse array for the sake of simplicity:
the Donors for a node are every >=0 nodes from the index node_index * 8 to 
node_index * 8 + nSdonors[node_index].

Returns
----------

The flat array of donor indices

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_SFD_donors_list",
      &D8connector<double>::template get_SFD_donors_list<std::vector<std::vector<int> > >,
      R"pdoc(
returns a list (not an array!) of irregular size with donor indices.

Description:
------------

returns a list (not an array!) of irregular size with donor indices. The first 
index is the node index and it points to a list of ndonors size with all the
SFD donors of the given node.

Returns
----------

THe 2D irregular list of node size

Authors:
--------
B.G.

)pdoc"
    )
    
    .def(
      "get_links",
      &D8connector<double>::template get_links<std::vector<std::uint8_t> >,
      R"pdoc(
returns an array of link size with link type.

Description:
------------

returns an array of link size with link type. See the documentation for the 
relationships between nodes and links. the node indices for a link index can be
found in the linknode array, at link_index*2 and link_index*2 + 1. The link type
is an uint8 number: 0 for inverse (node 2 give to node 1), 1 for normal (node 1
fives to node 2) or 3 invalid/inactive link.

Returns
----------

the array of link size with link code (uint8)

Authors:
--------
B.G.

)pdoc"
    )
    
    .def(
      "get_linknodes_flat",
      &D8connector<double>::template get_linknodes_flat<py::array_t<int,1>>,
      R"pdoc(
returns a flat array of linknodes (node indices pair for each links).

Description:
------------

1D flat array of linknodes. The node corresponding to link index li can be found 
at linknodes[li*2] and linknodes[li*2 + 1]. Non-existing or invalid links have 
node indices of -1.

Returns
----------

the array of 2 * link size with node indicies for each links

Authors:
--------
B.G.

)pdoc"
    )
    
    
    .def(
      "get_linknodes_list",
      &D8connector<double>::template get_linknodes_list<std::vector<std::vector<int> > >,
      R"pdoc(
returns a list (not an array!) of link nodes.

Description:
------------

2D link of link nodes. index is link index and it provides a list of nodes
composing the link.

Returns
----------

List of link size.

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_linknodes_list_oriented",
      &D8connector<double>::template get_linknodes_list_oriented<std::vector<std::vector<int> > >,
      R"pdoc(
returns a list (not an array!) of link nodes, donor first, rec second.

Description:
------------

2D link of link nodes. index is link index and it provides a list of nodes
composing the link. This version orients the list with the donor always first 
(and the receiver always second).

Returns
----------

List of link size.

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_SFD_receivers_at_node", 
      &D8connector<double>:: get_SFD_receivers_at_node,
      R"pdoc(Returns the node index of the SFD receivers for a given node index)pdoc"
      )
    
    .def(
      "get_SFD_dx_at_node", 
      &D8connector<double>:: get_SFD_dx_at_node,
      R"pdoc(Returns the distance to the SFD receivers for a given node index)pdoc"
      )
    
    .def(
      "get_SFD_ndonors_at_node", 
      &D8connector<double>:: get_SFD_ndonors_at_node

      )
    
    .def(
      "get_SFD_donors_at_node", 
      &D8connector<double>::template get_SFD_donors_at_node<std::vector<int> >,
      R"pdoc(Returns a list of SFD donors for a given node index)pdoc"
      )
    
    .def(
      "get_SFD_gradient", 
      &D8connector<double>::template get_SFD_gradient<py::array_t<double,1>, py::array_t<double,1> >,
      py::arg("topography"),
      R"pdoc(
returns an array of node size with the topographic gradient

Description:
------------

Takes a topography and returns the steepest descent gradient using precomputed
SFD receivers informations. Note that if you feed the function with a different 
topography than the one the graph has been computed with - or if local minimas 
have been resolved, you may obtain local negative values.

Parameters
-----------

- Topography (array 1D): the topography to use for the gradient computation

Returns
----------

Array of node size fo topographic gradient

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_links_gradient", 
      &D8connector<double>::template get_links_gradient< py::array_t<double,1>, py::array_t<double,1> >,
      py::arg("topography"), py::arg("minimum_slope"),
      R"pdoc(
returns an array of link size with the topographic gradient for each of them.

Description:
------------

Takes a topography and returns the gradient for each topographic link. It can 
recast teh negative/small gradients to a minimum value if needed.

Parameters
-----------

- Topography (array 1D): the topography to use for the gradient computation
- min_slope (float): the minimum slope - set to very negative value if not 
needed.

Returns
----------

Array of link size of topographic gradient

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_MFD_mean_gradient", 
      &D8connector<double>::template get_MFD_mean_gradient< py::array_t<double,1>, py::array_t<double,1> >,
      py::arg("topography"),
      R"pdoc(
returns an array of node size with the mean topographic gradient for each nodes.

Description:
------------

returns an array of node size with the mean topographic gradient for each nodes.
The gradient is computed for all receivers and averaged.

Parameters
-----------

- Topography (array 1D): the topography to use for the gradient computation


Returns
----------

Array of node size of mean topographic gradient

Authors:
--------
B.G.

)pdoc"

      )
    
    .def(
      "get_MFD_weighted_gradient", 
      &D8connector<double>::template get_MFD_weighted_gradient< py::array_t<double,1>, py::array_t<double,1> >,
      py::arg("topography"), py::arg("weights"),
      R"pdoc(
returns an array of node size with the weighted mean gradient for each nodes.

Description:
------------

returns an array of node size with the weighted mean gradient for each nodes.
It requires a weight for each links which can be obtained using the 
``get_link_weights`` function.

Parameters
-----------

- Topography (array 1D): the topography to use for the gradient computation
- weight (array 1D): link size array obtained with ``get_link_weights`` function


Returns
----------

Array of node size of weighted topographic gradient

Authors:
--------
B.G.

)pdoc"
      )
    
    .def(
      "get_link_weights", 
      &D8connector<double>::template get_link_weights< py::array_t<double,1>, py::array_t<double,1> >,
      py::arg("gradients"),py::arg("exponent"),
      R"pdoc(
Computes partition weights for each link function of rec slopes per node basis.

Description:
------------

Assign a [0,1] weight to each link to partition flow in MFD from every nodes to
their receivers. The partition uses an exponent to set the sensitivity to slope
differences. in a general manner, the higher the value, the more weight would 
go toward the steepest receiver. There are three particular cases (numerically 
and conceptually):

- exp = 0: equally parted towards all receivers regardless of their slopes
- exp = 0.5: proportional to the squareroot of the slopes
- exp = 1: perfectly proportional to the slopes

Parameters
-----------

- link_gradients (1D array): gradients per links (``get_link_gradients``)
- exponent (float): a positive exponent setting the sensitivity to slope

Returns
----------

Array of link size of partition weights

Authors:
--------
B.G.)pdoc"
      )
    
    .def(
      "set_stochaticiy_for_SFD", 
      &D8connector<double>::set_stochaticiy_for_SFD,
      py::arg("magnitude"),
      R"pdoc(EXPERIMENTAL: adds stochasticity to the SFD receivers calculation. Best to ignore.)pdoc"
      )
  ;

  py::class_<D4connector<double> >(m, "D4N", R"pdoc(DEPRECATED - will be back at some points, keeping for legacy)pdoc")
    .def(py::init<int,int,double,double,double,double>())
    .def("set_default_boundaries", &D4connector<double>::set_default_boundaries)
    .def("set_custom_boundaries", &D4connector<double>::set_custom_boundaries<py::array_t<int,1> >)
    .def("print_dim", &D4connector<double>::print_dim)
    .def("get_HS", &D4connector<double>::get_HS<std::vector<double>, py::array >)
    // .def("fill_barne_2014", &D4connector<double>::fill_barne_2014<std::vector<double> >)
    .def("get_mask_array",&D4connector<double>::get_mask_array)
    .def("set_values_at_boundaries", &D4connector<double>::set_values_at_boundaries<py::array_t<double,1> >)
    .def("set_out_boundaries_to_permissive", &D4connector<double>::set_out_boundaries_to_permissive)
    .def("get_boundary_at_node", &D4connector<double>::get_boundary_at_node)
  ;

  py::class_<numvec<double> >(m,"numvecf64")
    .def(py::init<py::array_t<double,1>&>())
    .def("get", &numvec<double>::get)
    .def("set", &numvec<double>::set)
  ;

  declare_graph<D8connector<double> >(m,"graph");
  // declare_graph<D4connector<double> >(m,"graphD4");

//=============================================================================================
//=============================================================================================
//===================== Standalone Algorithms =================================================
//=============================================================================================
//=============================================================================================

  m.def(
    "hillshade",
    &hillshade<D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1>, double>,
    py::arg("connector"),py::arg("topography"),
     R"pbdoc(
Hillshading function for visualisation

Description:
------------

Returns a [0,1] hillshade for a regular grid. Negative/nodata are set to 0

Parameters:
-----------

    * D8connector
    * Flat topography of node size (1D array)

Returns:
--------

    * Flat hillshade of node size (1D array)

Authors:
--------
B.G.

)pbdoc"
     );


  m.def(
    "rayshade",
    &rayshade< DAGGER::graph<double, DAGGER::D8connector<double> >, D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1>, double> 
    );

  m.def(
    "label_depressions_PQ",
    &label_depressions_PQ< py::array_t<double, 1>, py::array_t<int, 1>, D8connector<double> > 
    );

  m.def(
    "label_ocean",
    &label_ocean< py::array_t<double, 1>, py::array_t<int, 1>, D8connector<double> > 
    );


  m.def(
    "standalone_priority_flood",
    &standalone_priority_flood<D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1>, double >,
   py::arg("topography"), py::arg("connector")
  );

  m.def(
    "standalone_priority_flood_opti",
    &standalone_priority_flood_opti<D8connector<double>,  DAGGER::graph<double, DAGGER::D8connector<double> >, py::array_t<double, 1>, py::array_t<double, 1>, double >,
   py::arg("topography"), py::arg("connector"), py::arg("graph")
  );


  declare_popscape_old<DAGGER::D8connector<double> >(m,"popscape_old");
  declare_popscape<DAGGER::D8connector<double> >(m,"popscape");
  // // declare_popscape_old<DAGGER::D4connector<double> >(m,"popscape_oldD4");
  declare_trackscape<DAGGER::D8connector<double> >(m,"trackscape");
  // // declare_trackscape<DAGGER::D4connector<double> >(m,"trackscapeD4");

  py::enum_<RANDNOISE>(m, "NOISE")
    .value("WHITE", RANDNOISE::WHITE)
    .value("RED", RANDNOISE::RED)
    .value("PERLIN", RANDNOISE::PERLIN)
  ;


  declare_graphflood< double, DAGGER::graph<double, DAGGER::D8connector<double> >, DAGGER::D8connector<double> >(m, "graphflood");


  m.def("generate_perlin_noise_2D", &generate_perlin_noise_2D<py::array_t<double,1>, double, D8connector<double> >);

  




};
;





































// end of file




