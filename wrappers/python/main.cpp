#include "declare_connectors.hpp"
#include "declare_enums.hpp"
#include "declare_graphfloods.hpp"
#include "declare_graphs.hpp"
#include "declare_includes.hpp"
#include "declare_popscapes.hpp"
#include "declare_trackscapes.hpp"

using namespace DAGGER;

PYBIND11_MODULE(dagger, m)
{
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

	delclare_enums(m);
	declare_D8connector(m, "D8N");
	declare_graph<D8connector<FLOATING_POINT_DAGGER>>(m, "graph");

	//=============================================================================================
	//=============================================================================================
	//===================== Standalone Algorithms
	//=================================================
	//=============================================================================================
	//=============================================================================================

	m.def("hillshade",
				&hillshade<D8connector<FLOATING_POINT_DAGGER>,
									 py::array_t<FLOATING_POINT_DAGGER, 1>,
									 py::array_t<FLOATING_POINT_DAGGER, 1>,
									 FLOATING_POINT_DAGGER>,
				py::arg("connector"),
				py::arg("topography"),
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

)pbdoc");

	m.def("rayshade",
				&rayshade<DAGGER::graph<FLOATING_POINT_DAGGER,
																DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
									D8connector<FLOATING_POINT_DAGGER>,
									py::array_t<FLOATING_POINT_DAGGER, 1>,
									py::array_t<FLOATING_POINT_DAGGER, 1>,
									FLOATING_POINT_DAGGER>);

	m.def("set_BC_to_remove_seas",
				&set_BC_to_remove_seas<D8connector<FLOATING_POINT_DAGGER>,
															 py::array_t<FLOATING_POINT_DAGGER, 1>,
															 FLOATING_POINT_DAGGER>);

	m.def("label_depressions_PQ",
				&label_depressions_PQ<py::array_t<FLOATING_POINT_DAGGER, 1>,
															py::array_t<int, 1>,
															D8connector<FLOATING_POINT_DAGGER>>);

	m.def("label_ocean",
				&label_ocean<py::array_t<FLOATING_POINT_DAGGER, 1>,
										 py::array_t<int, 1>,
										 D8connector<FLOATING_POINT_DAGGER>>);

	m.def("standalone_priority_flood",
				&standalone_priority_flood<D8connector<FLOATING_POINT_DAGGER>,
																	 py::array_t<FLOATING_POINT_DAGGER, 1>,
																	 py::array_t<FLOATING_POINT_DAGGER, 1>,
																	 FLOATING_POINT_DAGGER>,
				py::arg("topography"),
				py::arg("connector"));

	m.def("standalone_priority_flood_opti",
				&standalone_priority_flood_opti<
					D8connector<FLOATING_POINT_DAGGER>,
					DAGGER::graph<FLOATING_POINT_DAGGER,
												DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
					py::array_t<FLOATING_POINT_DAGGER, 1>,
					py::array_t<FLOATING_POINT_DAGGER, 1>,
					FLOATING_POINT_DAGGER>,
				py::arg("topography"),
				py::arg("connector"),
				py::arg("graph"));

	m.def(
		"RiverNetwork",
		RiverNetwork<FLOATING_POINT_DAGGER,
								 DAGGER::D8connector<FLOATING_POINT_DAGGER>,
								 DAGGER::graph<FLOATING_POINT_DAGGER,
															 DAGGER::D8connector<FLOATING_POINT_DAGGER>>>);
	m.def(
		"DrainageDivides",
		DrainageDivides<FLOATING_POINT_DAGGER,
										DAGGER::D8connector<FLOATING_POINT_DAGGER>,
										DAGGER::graph<FLOATING_POINT_DAGGER,
																	DAGGER::D8connector<FLOATING_POINT_DAGGER>>>);

	// m.def(
	//   "check_connector_template",
	//   &check_connector_template< D8connector<FLOATING_POINT_DAGGER>,
	//   FLOATING_POINT_DAGGER >
	// );

	declare_popscape_old<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(
		m, "popscape_old");
	declare_popscape<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m, "popscape");
	// // declare_popscape_old<DAGGER::D4connector<FLOATING_POINT_DAGGER>
	// >(m,"popscape_oldD4");
	declare_trackscape<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m,
																																 "trackscape");
	// // declare_trackscape<DAGGER::D4connector<FLOATING_POINT_DAGGER>
	// >(m,"trackscapeD4");

	py::enum_<RANDNOISE>(m, "NOISE")
		.value("WHITE", RANDNOISE::WHITE)
		.value("RED", RANDNOISE::RED)
		.value("PERLIN", RANDNOISE::PERLIN);

	declare_graphflood<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER,
																	 DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
										 DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m,
																																 "graphflood");

	m.def("generate_perlin_noise_2D",
				&generate_perlin_noise_2D<py::array_t<FLOATING_POINT_DAGGER, 1>,
																	FLOATING_POINT_DAGGER,
																	D8connector<FLOATING_POINT_DAGGER>>);

	m.def("quick_fluvial_topo",
				&quick_fluvial_topo<float, py::array_t<float, 1>>);
};
;

// end of file
