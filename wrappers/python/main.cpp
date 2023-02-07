#include <pybind11/pybind11.h>
#include "D8connector.hpp"
#include "D4connector.hpp"
#include "graph.hpp"
#include "wrap_helper.hpp"

#include "hillshading.hpp"
#include "popscape.hpp"
#include "popscape_utils.hpp"
#include "fastflood.hpp"
#include "trackscape.hpp"
#include "utils.hpp"

using namespace DAGGER;



template<typename CONNECTOR_T>
void declare_graph(py::module &m, std::string typestr)
{

  py::class_<graph<double, CONNECTOR_T > >(m, typestr.c_str())
    .def(py::init<CONNECTOR_T&>())
    .def("init_graph", &graph<double,CONNECTOR_T>::init_graph)
    .def("set_opt_stst_rerouting", &graph<double,CONNECTOR_T>:: set_opt_stst_rerouting)
    .def("compute_graph", &graph<double,CONNECTOR_T>::template compute_graph<py::array_t<double,1>, py::array >)
    // .def("compute_graph_timer", &graph<double,CONNECTOR_T>::template compute_graph_timer<CONNECTOR_T, py::array_t<double,1>, py::array >)
    .def("is_Sstack_full", &graph<double,CONNECTOR_T>:: is_Sstack_full)
    .def("activate_opti_sparse_border_cordonnier", &graph<double,CONNECTOR_T>:: activate_opti_sparse_border_cordonnier)
    .def("get_all_nodes_upstream_of", &graph<double,CONNECTOR_T>::template get_all_nodes_upstream_of< py::array_t<int,1> > )
    .def("get_all_nodes_downstream_of", &graph<double,CONNECTOR_T>::template get_all_nodes_downstream_of< py::array_t<int,1> > )
    .def("get_SFD_stack",&graph<double,CONNECTOR_T>::template get_SFD_stack<py::array_t<size_t,1>>)
    .def("get_MFD_stack",&graph<double,CONNECTOR_T>::template get_MFD_stack<py::array_t<size_t,1>>)

    .def("accumulate_constant_downstream_SFD", &graph<double,CONNECTOR_T>::template accumulate_constant_downstream_SFD< py::array_t<double, 1> > )
    .def("accumulate_variable_downstream_SFD", &graph<double,CONNECTOR_T>::template accumulate_variable_downstream_SFD< py::array_t<double, 1>, py::array_t<double, 1> > )
    .def("accumulate_constant_downstream_MFD", &graph<double,CONNECTOR_T>::template accumulate_constant_downstream_MFD< py::array_t<double, 1>, py::array_t<double, 1> > )
    .def("accumulate_variable_downstream_MFD", &graph<double,CONNECTOR_T>::template accumulate_variable_downstream_MFD< py::array_t<double, 1>, py::array_t<double, 1> > )
    .def("set_LMR_method", &graph<double,CONNECTOR_T>:: set_LMR_method)
    .def("set_minimum_slope_for_LMR", &graph<double,CONNECTOR_T>:: set_minimum_slope_for_LMR)
    .def("set_slope_randomness_for_LMR", &graph<double,CONNECTOR_T>:: set_slope_randomness_for_LMR)

    // Distance functions
    .def("get_SFD_distance_from_outlets", &graph<double,CONNECTOR_T>::template get_SFD_distance_from_outlets< py::array_t<double,1> >)
    .def("get_SFD_min_distance_from_sources", &graph<double,CONNECTOR_T>::template get_SFD_min_distance_from_sources< py::array_t<double,1> >)
    .def("get_SFD_max_distance_from_sources", &graph<double,CONNECTOR_T>::template get_SFD_max_distance_from_sources< py::array_t<double,1> >)
    .def("get_MFD_max_distance_from_sources", &graph<double,CONNECTOR_T>::template get_MFD_max_distance_from_sources< py::array_t<double,1> >)
    .def("get_MFD_min_distance_from_sources", &graph<double,CONNECTOR_T>::template get_MFD_min_distance_from_sources< py::array_t<double,1> >)
    .def("get_MFD_max_distance_from_outlets", &graph<double,CONNECTOR_T>::template get_MFD_max_distance_from_outlets< py::array_t<double,1> >)
    .def("get_MFD_min_distance_from_outlets", &graph<double,CONNECTOR_T>::template get_MFD_min_distance_from_outlets< py::array_t<double,1> >)
    
    // Watershed labelling
    .def("get_SFD_basin_labels",  &graph<double,CONNECTOR_T>::template get_MFD_min_distance_from_outlets< py::array_t<int,1> >)
    
  ;
}

template<typename CONNECTOR_T>
void declare_popscape(py::module &m, std::string typestr)
{
  py::class_<popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T > >(m, typestr.c_str())
    .def(py::init<RANDNOISE,int,int,float_t,float_t>())
    // .def_readwrite("graph",  &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::graph)
    // .def_readwrite("connector",  &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::connector)
    .def("solve_generic", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::solve_generic)
    .def("get_topo", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_topo<py::array>)
    .def("get_QA", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_QA<py::array>)
    .def("compute_graph", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::compute_graph)
    .def("compute_DA_SFD", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::compute_DA_SFD)
    .def("apply_uplift", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::apply_uplift)
    .def("apply_variable_uplift", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template apply_variable_uplift<py::array_t<double,1> >)
    .def("solve_SFD_SPL_imp", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::solve_SFD_SPL_imp)
    .def("hydraulic_erosion_v0", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::hydraulic_erosion_v0)
    .def("normalise_topography", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::normalise_topography)
    .def("run_SFD_exp_latmag", &popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::run_SFD_exp_latmag)
    
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

PYBIND11_MODULE(dagger, m) {
  m.doc() = R"pbdoc(
      DAGGER - python API
      -----------------------
      .. rubric:: Classes

      .. autoclass:: D8N
        :members:

          


      .. rubric:: D8N Members

      .. autosummary::

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
          
  )pbdoc";


  py::enum_<DEPRES>(m, "LMR")
    .value("cordonnier_fill", DEPRES::cordonnier_fill)
    .value("cordonnier_carve", DEPRES::cordonnier_carve)
    .value("priority_flood", DEPRES::priority_flood)
    .value("none", DEPRES::none)
  ;

  py::class_<D8connector<double> >(m, "D8N", R"pbdoc(Connector for regular grids)pbdoc")
    .def(
      py::init<int,int,double,double,double,double>(),
      py::arg("nx"),py::arg("ny"),py::arg("dx"),py::arg("dy"),py::arg("x_min"), py::arg("x_max"),
      R"pbdoc(
Constructor for the D8 regular grid connector D8N

Description:
------------

Instanciate the D8N class: the built-in connector for regular grids. D8N regular grids are the most used type of grids, like most global DEMs for example.
They are defined by a number of rows (ny) and a number of columns (nx) with associated spacing between nodes in both directions (dy,dx).
The D8N constructor defaults to ``4edges`` boundary preset and precompiles neighbouring.

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
B.G. (08/2022)

)pbdoc"
      )
    .def(
      "set_default_boundaries", &D8connector<double>::set_default_boundaries,
      R"pbdoc(
Constructor for the D8 regular grid connector D8N

Description:
------------

Instanciate the D8N class: the built-in connector for regular grids. D8N regular grids are the most used type of grids, like most global DEMs for example.
They are defined by a number of rows (ny) and a number of columns (nx) with associated spacing between nodes in both directions (dy,dx).
The D8N constructor defaults to ``4edges`` boundary preset and precompiles neighbouring.

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
B.G. (08/2022)

)pbdoc"
      )
    .def("set_custom_boundaries", &D8connector<double>::set_custom_boundaries<py::array_t<int,1> >)
    .def("print_dim", &D8connector<double>::print_dim)
    .def("get_HS", &D8connector<double>::get_HS<std::vector<double>, py::array >)
    .def("get_mask_array",&D8connector<double>::get_mask_array)
    .def("set_values_at_boundaries", &D8connector<double>::set_values_at_boundaries<py::array_t<double,1> >)
    .def("set_out_boundaries_to_permissive", &D8connector<double>::set_out_boundaries_to_permissive)
    .def("get_boundary_at_node", &D8connector<double>::get_boundary_at_node)
    .def("get_rowcol_Sreceivers",&D8connector<double>:: get_rowcol_Sreceivers)
    .def("print_receivers", &D8connector<double>::template print_receivers<std::vector<double> >)
    .def("get_rec_array_size",&D8connector<double>::get_rec_array_size)
    .def("update_links_MFD_only", &D8connector<double>::template update_links_MFD_only<std::vector<double> >)
    .def("update_links", &D8connector<double>::template update_links<std::vector<double> >)
    .def("update_links_from_topo", &D8connector<double>::template update_links_from_topo<py::array_t<double,1> >)
    .def("sum_at_outlets", &D8connector<double>::template sum_at_outlets<py::array_t<double,1>, double >)
    .def("keep_only_at_outlets", &D8connector<double>::template keep_only_at_outlets<py::array_t<double,1>, py::array >)
    .def("get_SFD_receivers",&D8connector<double>::template get_SFD_receivers<py::array_t<int,1>>)
    .def("get_SFD_dx",&D8connector<double>::template get_SFD_dx<py::array_t<double,1>>)
    .def("get_SFD_ndonors",&D8connector<double>::template get_SFD_ndonors<py::array_t<int,1>>)
    .def("get_SFD_donors_flat",&D8connector<double>::template get_SFD_donors_flat<py::array_t<int,1>>)
    .def("get_SFD_donors_list",&D8connector<double>::template get_SFD_donors_list<std::vector<std::vector<int> > >)
    .def("get_links",&D8connector<double>::template get_links<std::vector<std::uint8_t> >)
    .def("get_linknodes_flat",&D8connector<double>::template get_linknodes_flat<py::array_t<int,1>>)
    .def("get_linknodes_list",&D8connector<double>::template get_linknodes_list<std::vector<std::vector<int> > >)
    .def("get_linknodes_list_oriented",&D8connector<double>::template get_linknodes_list_oriented<std::vector<std::vector<int> > >)
    .def("get_SFD_receivers_at_node", &D8connector<double>:: get_SFD_receivers_at_node)
    .def("get_SFD_dx_at_node", &D8connector<double>:: get_SFD_dx_at_node)
    .def("get_SFD_ndonors_at_node", &D8connector<double>:: get_SFD_ndonors_at_node)
    .def("get_SFD_donors_at_node", &D8connector<double>::template get_SFD_donors_at_node<std::vector<int> >)
    .def("get_SFD_gradient", &D8connector<double>::template get_SFD_gradient<py::array_t<double,1>, py::array_t<double,1> >)
    .def("get_links_gradient", &D8connector<double>::template get_links_gradient< py::array_t<double,1>, py::array_t<double,1> >)
    .def("get_MFD_mean_gradient", &D8connector<double>::template get_MFD_mean_gradient< py::array_t<double,1>, py::array_t<double,1> >)
    .def("get_MFD_weighted_gradient", &D8connector<double>::template get_MFD_weighted_gradient< py::array_t<double,1>, py::array_t<double,1> >)
    .def("get_link_weights", &D8connector<double>::template get_link_weights< py::array_t<double,1>, py::array_t<double,1> >)
    .def("set_stochaticiy_for_SFD", &D8connector<double>::set_stochaticiy_for_SFD)
  ;

  py::class_<D4connector<double> >(m, "D4N")
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
B.G. (08/2022)

)pbdoc"
     );


  declare_popscape<DAGGER::D8connector<double> >(m,"popscape");
  // // declare_popscape<DAGGER::D4connector<double> >(m,"popscapeD4");
  declare_trackscape<DAGGER::D8connector<double> >(m,"trackscape");
  // // declare_trackscape<DAGGER::D4connector<double> >(m,"trackscapeD4");

  py::enum_<RANDNOISE>(m, "NOISE")
    .value("WHITE", RANDNOISE::WHITE)
    .value("RED", RANDNOISE::RED)
    .value("PERLIN", RANDNOISE::PERLIN)
  ;

  py::class_<fastflood_recorder<double> >(m, "fastflood_recorder")
    .def(py::init<>())
    .def("enable_edot_recording",&fastflood_recorder<double>::enable_edot_recording)
    .def("disable_edot_recording",&fastflood_recorder<double>::disable_edot_recording)
    .def("get_edot",&fastflood_recorder<double>::get_edot<py::array_t<double,1> >)
    .def("enable_ddot_recording",&fastflood_recorder<double>::enable_ddot_recording)
    .def("disable_ddot_recording",&fastflood_recorder<double>::disable_ddot_recording)
    .def("get_ddot",&fastflood_recorder<double>::get_ddot<py::array_t<double,1> >)
    .def("enable_lateral_edot_recording",&fastflood_recorder<double>::enable_lateral_edot_recording)
    .def("disable_lateral_edot_recording",&fastflood_recorder<double>::disable_lateral_edot_recording)
    .def("get_lateral_edot",&fastflood_recorder<double>::get_lateral_edot<py::array_t<double,1> >)
    .def("enable_lateral_ddot_recording",&fastflood_recorder<double>::enable_lateral_ddot_recording)
    .def("disable_lateral_ddot_recording",&fastflood_recorder<double>::disable_lateral_ddot_recording)
    .def("get_lateral_ddot",&fastflood_recorder<double>::get_lateral_ddot<py::array_t<double,1> >)
    .def("enable_qs_recording",&fastflood_recorder<double>::enable_qs_recording)
    .def("disable_qs_recording",&fastflood_recorder<double>::disable_qs_recording)
    .def("get_qs",&fastflood_recorder<double>::get_qs<py::array_t<double,1> >)
    .def("enable_dhw_recording",&fastflood_recorder<double>::enable_dhw_recording)
    .def("disable_dhw_recording",&fastflood_recorder<double>::disable_dhw_recording)
    .def("get_dhw",&fastflood_recorder<double>::get_dhw<py::array_t<double,1> >)
    .def("enable_tau_recording",&fastflood_recorder<double>::enable_tau_recording)
    .def("disable_tau_recording",&fastflood_recorder<double>::disable_tau_recording)
    .def("get_tau",&fastflood_recorder<double>::get_tau<py::array_t<double,1> >)
    .def("enable_vmot_recording",&fastflood_recorder<double>::enable_vmot_recording)
    .def("disable_vmot_recording",&fastflood_recorder<double>::disable_vmot_recording)
    .def("get_vmot",&fastflood_recorder<double>::get_vmot<py::array_t<double,1> >)
  ;

  declare_ff<DAGGER::D8connector<double> >(m,"FF");
  // // declare_ff<DAGGER::D4connector<double> >(m,"FFD4");


  // m.def("generate_perlin_noise_2D", &generate_perlin_noise_2D<py::array_t<double,1>, double, D8connector<double> >);

  




};
;





































// end of file




