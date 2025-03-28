#include "D4connector.hpp"
#include "D8connector.hpp"
#include "graph.hpp"
#include "wrap_helper.hpp"
#include <pybind11/pybind11.h>

#include "hillshading.hpp"

namespace DAGGER {

template<typename CONNECTOR_T>
void
declare_graph(py::module& m, std::string typestr)
{

	py::class_<graph<double, CONNECTOR_T>>(m, typestr.c_str())
		.def(py::init<int, int>())
		.def(
			"get_DA_proposlope",
			&graph<double, CONNECTOR_T>::
				template get_DA_proposlope<CONNECTOR_T, std::vector<double>, py::array>)
		.def("get_DA_SS",
				 &graph<double, CONNECTOR_T>::
					 template get_DA_SS<CONNECTOR_T, std::vector<double>, py::array>)
		.def(
			"get_rowcol_Sreceivers",
			&graph<double, CONNECTOR_T>::template get_rowcol_Sreceivers<CONNECTOR_T>)
		.def("print_receivers",
				 &graph<double, CONNECTOR_T>::
					 template print_receivers<CONNECTOR_T, std::vector<double>>)
		.def("get_rec_array_size", &graph<double, CONNECTOR_T>::get_rec_array_size)
		.def("init_graph",
				 &graph<double, CONNECTOR_T>::template init_graph<CONNECTOR_T>)
		.def(
			"update_Mrecs",
			&graph<double, CONNECTOR_T>::template update_Mrecs<CONNECTOR_T,
																												 std::vector<double>>)
		.def("update_recs",
				 &graph<double, CONNECTOR_T>::template update_recs<CONNECTOR_T,
																													 std::vector<double>>)
		.def("set_opt_stst_rerouting",
				 &graph<double, CONNECTOR_T>::set_opt_stst_rerouting)
		.def(
			"compute_graph",
			&graph<double, CONNECTOR_T>::
				template compute_graph<CONNECTOR_T, py::array_t<double, 1>, py::array>)
		.def(
			"compute_graph_timer",
			&graph<double,
						 CONNECTOR_T>::template compute_graph_timer<CONNECTOR_T,
																												py::array_t<double, 1>,
																												py::array>)
		.def("test_Srecs",
				 &graph<double, CONNECTOR_T>::template test_Srecs<py::array>)
		.def("sum_at_outlets",
				 &graph<double, CONNECTOR_T>::
					 template sum_at_outlets<CONNECTOR_T, py::array_t<double, 1>, double>)
		.def(
			"keep_only_at_outlets",
			&graph<double,
						 CONNECTOR_T>::template keep_only_at_outlets<CONNECTOR_T,
																												 py::array_t<double, 1>,
																												 py::array>)
		.def("is_Sstack_full", &graph<double, CONNECTOR_T>::is_Sstack_full)
		.def("has_Srecs", &graph<double, CONNECTOR_T>::has_Srecs)
		.def("activate_opti_sparse_border_cordonnier",
				 &graph<double, CONNECTOR_T>::activate_opti_sparse_border_cordonnier)
		.def("get_all_nodes_upstream_of",
				 &graph<double, CONNECTOR_T>::
					 template get_all_nodes_upstream_of<CONNECTOR_T, py::array_t<int, 1>>)
		.def(
			"get_all_nodes_downstream_of",
			&graph<double, CONNECTOR_T>::
				template get_all_nodes_downstream_of<CONNECTOR_T, py::array_t<int, 1>>)
		.def("get_SFD_receivers",
				 &graph<double,
								CONNECTOR_T>::template get_SFD_receivers<py::array_t<int, 1>>)
		.def(
			"get_SFD_dx",
			&graph<double, CONNECTOR_T>::template get_SFD_dx<py::array_t<double, 1>>)
		.def("get_SFD_ndonors",
				 &graph<double,
								CONNECTOR_T>::template get_SFD_ndonors<py::array_t<int, 1>>)
		.def("get_SFD_donors_flat",
				 &graph<double,
								CONNECTOR_T>::template get_SFD_donors_flat<py::array_t<int, 1>>)
		.def("get_SFD_donors_list",
				 &graph<double, CONNECTOR_T>::template get_SFD_donors_list<
					 std::vector<std::vector<int>>>)
		.def("get_SFD_stack",
				 &graph<double,
								CONNECTOR_T>::template get_SFD_stack<py::array_t<size_t, 1>>)
		.def("get_MFD_stack",
				 &graph<double,
								CONNECTOR_T>::template get_MFD_stack<py::array_t<size_t, 1>>)
		.def("get_links",
				 &graph<double, CONNECTOR_T>::template get_links<std::vector<bool>>)
		.def("get_linknodes_flat",
				 &graph<double,
								CONNECTOR_T>::template get_linknodes_flat<py::array_t<int, 1>>)
		.def(
			"get_linknodes_flat_D4",
			&graph<double,
						 CONNECTOR_T>::template get_linknodes_flat_D4<py::array_t<int, 1>>)
		.def("get_linkdx_flat_D4",
				 &graph<double, CONNECTOR_T>::
					 template get_linkdx_flat_D4<py::array_t<double, 1>, CONNECTOR_T>)

		.def("get_linknodes_list",
				 &graph<double, CONNECTOR_T>::template get_linknodes_list<
					 std::vector<std::vector<int>>>)
		.def("get_linknodes_list_oriented",
				 &graph<double, CONNECTOR_T>::template get_linknodes_list_oriented<
					 std::vector<std::vector<int>>>)
		.def("get_SFD_receivers_at_node",
				 &graph<double, CONNECTOR_T>::get_SFD_receivers_at_node)
		.def("get_SFD_dx_at_node", &graph<double, CONNECTOR_T>::get_SFD_dx_at_node)
		.def("get_SFD_ndonors_at_node",
				 &graph<double, CONNECTOR_T>::get_SFD_ndonors_at_node)
		.def("get_SFD_donors_at_node",
				 &graph<double,
								CONNECTOR_T>::template get_SFD_donors_at_node<std::vector<int>>)
		.def("get_SFD_gradient",
				 &graph<double,
								CONNECTOR_T>::template get_SFD_gradient<py::array_t<double, 1>,
																												py::array_t<double, 1>>)
		.def(
			"get_links_gradient",
			&graph<double,
						 CONNECTOR_T>::template get_links_gradient<CONNECTOR_T,
																											 py::array_t<double, 1>,
																											 py::array_t<double, 1>>)
		.def("get_MFD_mean_gradient",
				 &graph<double, CONNECTOR_T>::template get_MFD_mean_gradient<
					 CONNECTOR_T,
					 py::array_t<double, 1>,
					 py::array_t<double, 1>>)
		.def("get_MFD_weighted_gradient",
				 &graph<double, CONNECTOR_T>::template get_MFD_weighted_gradient<
					 CONNECTOR_T,
					 py::array_t<double, 1>,
					 py::array_t<double, 1>>)
		.def("get_link_weights",
				 &graph<double,
								CONNECTOR_T>::template get_link_weights<py::array_t<double, 1>,
																												py::array_t<double, 1>>)
		.def(
			"accumulate_constant_downstream_SFD",
			&graph<double, CONNECTOR_T>::template accumulate_constant_downstream_SFD<
				CONNECTOR_T,
				py::array_t<double, 1>>)
		.def(
			"accumulate_variable_downstream_SFD",
			&graph<double, CONNECTOR_T>::template accumulate_variable_downstream_SFD<
				CONNECTOR_T,
				py::array_t<double, 1>,
				py::array_t<double, 1>>)
		.def(
			"accumulate_constant_downstream_MFD",
			&graph<double, CONNECTOR_T>::template accumulate_constant_downstream_MFD<
				CONNECTOR_T,
				py::array_t<double, 1>,
				py::array_t<double, 1>>)
		.def(
			"accumulate_variable_downstream_MFD",
			&graph<double, CONNECTOR_T>::template accumulate_variable_downstream_MFD<
				CONNECTOR_T,
				py::array_t<double, 1>,
				py::array_t<double, 1>>)
		.def("set_LMR_method", &graph<double, CONNECTOR_T>::set_LMR_method)
		.def("set_minimum_slope_for_LMR",
				 &graph<double, CONNECTOR_T>::set_minimum_slope_for_LMR)
		.def("set_slope_randomness_for_LMR",
				 &graph<double, CONNECTOR_T>::set_slope_randomness_for_LMR)

		// Distance functions
		.def("get_SFD_distance_from_outlets",
				 &graph<double, CONNECTOR_T>::template get_SFD_distance_from_outlets<
					 CONNECTOR_T,
					 py::array_t<double, 1>>)
		.def(
			"get_SFD_min_distance_from_sources",
			&graph<double, CONNECTOR_T>::template get_SFD_min_distance_from_sources<
				CONNECTOR_T,
				py::array_t<double, 1>>)
		.def(
			"get_SFD_max_distance_from_sources",
			&graph<double, CONNECTOR_T>::template get_SFD_max_distance_from_sources<
				CONNECTOR_T,
				py::array_t<double, 1>>)
		.def(
			"get_MFD_max_distance_from_sources",
			&graph<double, CONNECTOR_T>::template get_MFD_max_distance_from_sources<
				CONNECTOR_T,
				py::array_t<double, 1>>)
		.def(
			"get_MFD_min_distance_from_sources",
			&graph<double, CONNECTOR_T>::template get_MFD_min_distance_from_sources<
				CONNECTOR_T,
				py::array_t<double, 1>>)
		.def(
			"get_MFD_max_distance_from_outlets",
			&graph<double, CONNECTOR_T>::template get_MFD_max_distance_from_outlets<
				CONNECTOR_T,
				py::array_t<double, 1>>)
		.def(
			"get_MFD_min_distance_from_outlets",
			&graph<double, CONNECTOR_T>::template get_MFD_min_distance_from_outlets<
				CONNECTOR_T,
				py::array_t<double, 1>>)

		// Watershed labelling
		.def(
			"get_SFD_basin_labels",
			&graph<double, CONNECTOR_T>::template get_MFD_min_distance_from_outlets<
				CONNECTOR_T,
				py::array_t<int, 1>>)

		;
};

PYBIND11_MODULE(dagger, m)
{
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
		.value("none", DEPRES::none);

	py::class_<D8connector<double>>(m, "D8N")
		.def(py::init<int, int, double, double, double, double>())
		.def("set_default_boundaries", &D8connector<double>::set_default_boundaries)
		.def("set_custom_boundaries",
				 &D8connector<double>::set_custom_boundaries<py::array_t<int, 1>>)
		.def("print_dim", &D8connector<double>::print_dim)
		.def("get_HS", &D8connector<double>::get_HS<std::vector<double>, py::array>)
		// .def("fill_barne_2014",
		// &D8connector<double>::fill_barne_2014<std::vector<double> >)
		.def("get_mask_array", &D8connector<double>::get_mask_array)
		.def("set_values_at_boundaries",
				 &D8connector<double>::set_values_at_boundaries<py::array_t<double, 1>>)
		.def("set_out_boundaries_to_permissive",
				 &D8connector<double>::set_out_boundaries_to_permissive)
		.def("get_boundary_at_node", &D8connector<double>::get_boundary_at_node);

	py::class_<D4connector<double>>(m, "D4N")
		.def(py::init<int, int, double, double, double, double>())
		.def("set_default_boundaries", &D4connector<double>::set_default_boundaries)
		.def("set_custom_boundaries",
				 &D4connector<double>::set_custom_boundaries<py::array_t<int, 1>>)
		.def("print_dim", &D4connector<double>::print_dim)
		.def("get_HS", &D4connector<double>::get_HS<std::vector<double>, py::array>)
		// .def("fill_barne_2014",
		// &D4connector<double>::fill_barne_2014<std::vector<double> >)
		.def("get_mask_array", &D4connector<double>::get_mask_array)
		.def("set_values_at_boundaries",
				 &D4connector<double>::set_values_at_boundaries<py::array_t<double, 1>>)
		.def("set_out_boundaries_to_permissive",
				 &D4connector<double>::set_out_boundaries_to_permissive)
		.def("get_boundary_at_node", &D4connector<double>::get_boundary_at_node);

	py::class_<numvec<double>>(m, "numvecf64")
		.def(py::init<py::array_t<double, 1>&>())
		.def("get", &numvec<double>::get)
		.def("set", &numvec<double>::set);

	m.def("hillshade",
				&hillshade<D8connector<double>,
									 py::array_t<double, 1>,
									 py::array_t<double, 1>,
									 double>);
	declare_graph<D8connector<double>>(m, "graph");
	declare_graph<D4connector<double>>(m, "graphD4");
};
;

// end of namespace DAGGER
}; // namespace DAGGER

// end of file
