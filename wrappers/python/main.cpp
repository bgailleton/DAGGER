#include <pybind11/pybind11.h>
#include "D8connector.hpp"
#include "graph.hpp"
#include "wrap_helper.hpp"

#include "hillshading.hpp"


namespace DAGGER
{

PYBIND11_MODULE(dagger, m) {
  m.doc() = R"pbdoc(
      dagger example plugin
      -----------------------
      .. currentmodule:: dagger
      .. autosummary::
         :toctree: _generate

  )pbdoc";


  py::enum_<DEPRES>(m, "LMR")
    .value("cordonnier_fill", DEPRES::cordonnier_fill)
    .value("cordonnier_carve", DEPRES::cordonnier_carve)
    .value("priority_flood", DEPRES::priority_flood)
    .value("none", DEPRES::none)
  ;

	py::class_<D8connector<double> >(m, "D8N")
		.def(py::init<int,int,double,double,double,double>())
		.def("set_default_boundaries", &D8connector<double>::set_default_boundaries)
		.def("set_custom_boundaries", &D8connector<double>::set_custom_boundaries<py::array_t<int,1> >)
		.def("print_dim", &D8connector<double>::print_dim)
		.def("get_HS", &D8connector<double>::get_HS<std::vector<double>, py::array >)
		// .def("fill_barne_2014", &D8connector<double>::fill_barne_2014<std::vector<double> >)
		.def("get_mask_array",&D8connector<double>::get_mask_array)
		.def("set_values_at_boundaries", &D8connector<double>::set_values_at_boundaries<py::array_t<double,1> >)
	;

	py::class_<graph<double> >(m, "graph")
		.def(py::init<int,int>())
		.def("get_DA_proposlope", &graph<double>::get_DA_proposlope<D8connector<double>, std::vector<double>, py::array >)
		.def("get_DA_SS", &graph<double>::get_DA_SS<D8connector<double>, std::vector<double>, py::array >)
		.def("get_rowcol_Sreceivers",&graph<double>::get_rowcol_Sreceivers<D8connector<double> >)
		.def("print_receivers", &graph<double>::print_receivers<D8connector<double>, std::vector<double> >)
		.def("get_rec_array_size",&graph<double>::get_rec_array_size)
		.def("init_graph", &graph<double>::init_graph<D8connector<double> >)
		.def("compute_graph", &graph<double>::compute_graph<D8connector<double>, py::array_t<double,1>, py::array >)
		.def("test_Srecs", &graph<double>::test_Srecs<py::array >)
		.def("sum_at_outlets", &graph<double>::sum_at_outlets<D8connector<double>, py::array_t<double,1>, double >)
		.def("keep_only_at_outlets", &graph<double>::keep_only_at_outlets<D8connector<double>, py::array_t<double,1>, py::array >)
		.def("is_Sstack_full", &graph<double>::is_Sstack_full)
		.def("has_Srecs", &graph<double>::has_Srecs)
		.def("get_all_nodes_upstream_of", &graph<double>::get_all_nodes_upstream_of< D8connector<double>, py::array_t<int,1> > )
		.def("get_all_nodes_downstream_of", &graph<double>::get_all_nodes_downstream_of< D8connector<double>, py::array_t<int,1> > )
		.def("get_SFD_receivers",&graph<double>::get_SFD_receivers<py::array_t<int,1>>)
		.def("get_SFD_dx",&graph<double>::get_SFD_dx<py::array_t<double,1>>)
		.def("get_SFD_ndonors",&graph<double>::get_SFD_ndonors<py::array_t<int,1>>)
		.def("get_SFD_donors_flat",&graph<double>::get_SFD_donors_flat<py::array_t<int,1>>)
		.def("get_SFD_donors_list",&graph<double>::get_SFD_donors_list<std::vector<std::vector<int> > >)
		.def("get_SFD_stack",&graph<double>::get_SFD_stack<py::array_t<size_t,1>>)
		.def("get_MFD_stack",&graph<double>::get_MFD_stack<py::array_t<size_t,1>>)
		.def("get_links",&graph<double>::get_links<std::vector<bool> >)
		.def("get_linknodes_flat",&graph<double>::get_linknodes_flat<py::array_t<int,1>>)
		.def("get_linknodes_list",&graph<double>::get_linknodes_list<std::vector<std::vector<int> > >)
		.def("get_linknodes_list_oriented",&graph<double>::get_linknodes_list_oriented<std::vector<std::vector<int> > >)
		.def("get_SFD_receivers_at_node", &graph<double>::get_SFD_receivers_at_node)
		.def("get_SFD_dx_at_node", &graph<double>::get_SFD_dx_at_node)
		.def("get_SFD_ndonors_at_node", &graph<double>::get_SFD_ndonors_at_node)
		.def("get_SFD_donors_at_node", &graph<double>::get_SFD_donors_at_node<std::vector<int> >)
		.def("get_SFD_gradient", &graph<double>::get_SFD_gradient<py::array_t<double,1>, py::array_t<double,1> >)
		.def("get_links_gradient", &graph<double>::get_links_gradient< D8connector<double>, py::array_t<double,1>, py::array_t<double,1> >)
		.def("get_MFD_mean_gradient", &graph<double>::get_MFD_mean_gradient< D8connector<double>, py::array_t<double,1>, py::array_t<double,1> >)
		.def("get_MFD_weighted_gradient", &graph<double>::get_MFD_weighted_gradient< D8connector<double>, py::array_t<double,1>, py::array_t<double,1> >)
		.def("get_link_weights", &graph<double>::get_link_weights< py::array_t<double,1>, py::array_t<double,1> >)
		.def("accumulate_constant_downstream_SFD", &graph<double>::accumulate_constant_downstream_SFD< D8connector<double>, py::array_t<double, 1> > )
		.def("accumulate_variable_downstream_SFD", &graph<double>::accumulate_variable_downstream_SFD< D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1> > )
		.def("accumulate_constant_downstream_MFD", &graph<double>::accumulate_constant_downstream_MFD< D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1> > )
		.def("accumulate_variable_downstream_MFD", &graph<double>::accumulate_variable_downstream_MFD< D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1> > )
		.def("speed_test_links",&graph<double>::speed_test_links<D8connector<double> >)
		.def("set_LMR_method", &graph<double>::set_LMR_method)
		.def("set_minimum_slope_for_LMR", &graph<double>::set_minimum_slope_for_LMR)
		.def("set_slope_randomness_for_LMR", &graph<double>::set_slope_randomness_for_LMR)
	;

	py::class_<numvec<double> >(m,"numvecf64")
		.def(py::init<py::array_t<double,1>&>())
		.def("get", &numvec<double>::get)
		.def("set", &numvec<double>::set)
	;

	m.def("hillshade",&hillshade<D8connector<double>, py::array_t<double, 1>, py::array_t<double, 1>, double>);

};
;




// end of namespace DAGGER
};



































// end of file




