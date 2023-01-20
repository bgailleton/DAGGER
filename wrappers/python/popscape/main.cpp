#include <pybind11/pybind11.h>
#include "wrap_helper_python.hpp"
#include "D8connector.hpp"
#include "D4connector.hpp"
#include "graph.hpp"
#include "hillshading.hpp"
#include "popscape.hpp"
#include "trackscape.hpp"

#ifdef OPENMP_YOLO
#include<omp.h>
#endif


namespace py = pybind11;

namespace popscape
{

template<typename CONNECTOR_T>
void declare_popscape(py::module &m, std::string typestr)
{
	py::class_<popscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T > >(m, typestr.c_str())
		.def(py::init<RANDNOISE,int,int,float_t,float_t>())
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
		.def("run_SFD", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::run_SFD)
		.def("block_uplift", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::block_uplift)
		.def("external_uplift", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template external_uplift<py::array_t<double,1>& >)		
		.def("init_TSP_module", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template init_TSP_module<py::array_t<double,1>& >)		
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
		.def("hillslopes_on", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::hillslopes_on)
		.def("hillslopes_off", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::hillslopes_off)
		.def("fluvial_on", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::fluvial_on)
		.def("fluvial_off", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::fluvial_off)
		.def("fill_up", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::fill_up)
		.def("init_Ch_MTSI", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::init_Ch_MTSI)
		.def("rise_boundary_by", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::rise_boundary_by)
		.def("get_TSP_surface_concentrations", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_TSP_surface_concentrations<py::array>)
		.def("get_Ch_MTIS_surface_age", &trackscape<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T >::template get_Ch_MTIS_surface_age<py::array>)

		

		
		
	;
}

PYBIND11_MODULE(popscape_dagger, m) {
  m.doc() = R"pbdoc(
      popscape-dagger example plugin
      -----------------------
      .. currentmodule:: popscape-dagger
      .. autosummary::
         :toctree: _generate

  )pbdoc";

	declare_popscape<DAGGER::D8connector<double> >(m,"popscape");
	declare_popscape<DAGGER::D4connector<double> >(m,"popscapeD4");
	declare_trackscape<DAGGER::D8connector<double> >(m,"trackscape");
	declare_trackscape<DAGGER::D4connector<double> >(m,"trackscapeD4");

	py::enum_<RANDNOISE>(m, "NOISE")
    .value("WHITE", RANDNOISE::WHITE)
    .value("RED", RANDNOISE::RED)
    .value("PERLIN", RANDNOISE::PERLIN)
  ;


  


};


};	




















