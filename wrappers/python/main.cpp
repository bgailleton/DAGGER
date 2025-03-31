
#include "declare_algo.hpp"
#include "declare_connectors.hpp"
#include "declare_dbag.hpp"
#include "declare_enums.hpp"
#include "declare_graphfloods.hpp"
#include "declare_graphs.hpp"
#include "declare_includes.hpp"
#include "declare_parambag.hpp"
#include "declare_popscapes.hpp"
#include "declare_rivnets.hpp"
#include "declare_trackscapes.hpp"

#include "rd_graph.hpp"
#include "rd_neighbourer.hpp"
#include "rd_select_watershed.hpp"

#include "pybind11/pybind11.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

using namespace DAGGER;

PYBIND11_MODULE(dagger, m)
{

	xt::import_numpy();
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
	declare_dbag(m);
	declare_param(m);
	declare_D8connector(m, "D8N");
	declare_graph<D8connector<FLOATING_POINT_DAGGER>>(m, "graph");

	//=============================================================================================
	//=============================================================================================
	//===================== Standalone Algorithms
	//=================================================
	//=============================================================================================
	//=============================================================================================

	declare_algos(m);

	// m.def(
	//   "check_connector_template",
	//   &check_connector_template< D8connector<FLOATING_POINT_DAGGER>,
	//   FLOATING_POINT_DAGGER >
	// );

	declare_popscape_old<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(
		m, "popscape_old");
	declare_popscape<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m, "popscape");
	declare_trackscape<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m,
																																 "trackscape");

	declare_graphflood<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER,
																	 DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
										 DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m,
																																 "graphflood");
	declare_rivnet(m);

	// Xtensor-python backend

	py::enum_<boundaries>(m, "boundaries")
		.value("normal", boundaries::normal)
		.value("periodicEW", boundaries::periodicEW)
		.value("periodicNS", boundaries::periodicNS)
		.value("customs", boundaries::customs);

	py::class_<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>>(m,
																																 "GridCPP_f32")
		.def(py::init<int, int, float, float, std::uint8_t>());

	py::class_<GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>>(m,
																																	"GridCPP_f64")
		.def(py::init<int, int, double, double, std::uint8_t>());

	m.def("_PriorityFlood_D4_f64",
				&_PriorityFlood_D4<xt::pytensor<double, 2>,
													 double,
													 GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
													 xt::pytensor<std::uint8_t, 2>>);
	m.def("_PriorityFlood_D4_f32",
				&_PriorityFlood_D4<xt::pytensor<float, 2>,
													 float,
													 GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
													 xt::pytensor<std::uint8_t, 2>>);

	//   template<class GRID_T, class ARR_INT, class ARR_FT, class ARR_BCs>
	m.def(
		"compute_SF_stack_D4_full_f32",
		&compute_full_SF_graph<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
													 float>);

	m.def(
		"compute_SF_stack_D4_full_f64",
		&compute_full_SF_graph<GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
													 double>);

	m.def(
		"BCs_to_mask_f64",
		&BCs_to_mask<GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>, double>);
	m.def(
		"BCs_to_mask_f32",
		&BCs_to_mask<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>, float>);
	m.def(
		"mask_to_BCs_f64",
		&mask_to_BCs<GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>, double>);
	m.def(
		"mask_to_BCs_f32",
		&mask_to_BCs<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>, float>);
	m.def("mask_watersheds_above_elevations_f64",
				&mask_watersheds_above_elevations<
					GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
					double>);
	m.def("mask_watersheds_above_elevations_f32",
				&mask_watersheds_above_elevations<
					GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
					float>);
	m.def(
		"label_watersheds_mask_f32",
		&label_watersheds_mask<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
													 float>);
	m.def(
		"label_watersheds_mask_f64",
		&label_watersheds_mask<GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
													 double>);
	m.def("mask_watersheds_min_area_f32",
				&mask_watersheds_min_area<
					GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
					float>);
	m.def("mask_watersheds_min_area_f64",
				&mask_watersheds_min_area<
					GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
					double>);
	m.def("correct_Sreceivers_from_mask_f32",
				&correct_Sreceivers_from_mask<
					GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
					float>);
	m.def("correct_Sreceivers_from_mask_f64",
				&correct_Sreceivers_from_mask<
					GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
					double>);
	m.def(
		"bounding_box_from_label_f32",
		&bounding_box_from_label<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
														 float>);
	m.def("bounding_box_from_label_f64",
				&bounding_box_from_label<
					GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
					double>);
	m.def("mask_upstream_MFD_f32",
				&mask_upstream_MFD<GridCPP<int, float, xt::pytensor<std::uint8_t, 2>>,
													 float>);
	m.def("mask_upstream_MFD_f64",
				&mask_upstream_MFD<GridCPP<int, double, xt::pytensor<std::uint8_t, 2>>,
													 double>);
};
;

// end of file
