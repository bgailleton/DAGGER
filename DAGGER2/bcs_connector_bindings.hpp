#pragma once

#include "dg2_BCs.hpp"
#include "dg2_BCs_helper.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace dagger2 {

// ==============================================
// ARRAY BINDINGS
// ==============================================

// Pybind11 binding helper
template<typename T>
void
bind_array_types(py::module& m, const std::string& suffix)
{
	using AR = ArrayRef<T>;
	using G2D = Grid2D<T>;
	using G3D = Grid3D<T>;

	py::class_<AR, std::shared_ptr<AR>>(m, ("ArrayRef" + suffix).c_str())
		.def(py::init<py::array_t<T>>())
		.def("__getitem__", [](AR& a, size_t i) { return a[i]; })
		.def("__setitem__", [](AR& a, size_t i, T val) { a[i] = val; })
		.def("__len__", &AR::size)
		.def_property_readonly("size", &AR::size)
		.def("as_numpy", &AR::as_numpy);

	py::class_<G2D, std::shared_ptr<G2D>>(m, ("Grid2D" + suffix).c_str())
		.def(py::init<py::array_t<T>, size_t, size_t>(),
				 py::arg("array"),
				 py::arg("rows") = 0,
				 py::arg("cols") = 0)
		.def(
			"__call__",
			[](G2D& g, size_t r, size_t c) -> T& { return g(r, c); },
			py::return_value_policy::reference_internal)
		.def("__getitem__", [](G2D& g, size_t i) { return g[i]; })
		.def("__setitem__", [](G2D& g, size_t i, T val) { g[i] = val; })
		.def_property_readonly("rows", &G2D::rows)
		.def_property_readonly("cols", &G2D::cols)
		.def_property_readonly("size", &G2D::size)
		.def("as_numpy", &G2D::as_numpy_2d)
		.def("get_arref", &G2D::get_arref)
		// Copy operations
		.def("copy", &G2D::copy, "Create a deep copy of this Grid2D")
		.def("copy_as_f32", &G2D::template copy_as<float>, "Copy as float Grid2D")
		.def("copy_as_f64", &G2D::template copy_as<double>, "Copy as double Grid2D")
		.def("copy_as_i32", &G2D::template copy_as<int32_t>, "Copy as int32 Grid2D")
		.def("copy_as_i64", &G2D::template copy_as<int64_t>, "Copy as int64 Grid2D")
		.def("copy_from",
				 static_cast<void (G2D::*)(const G2D&)>(&G2D::copy_from),
				 "Copy data from another Grid2D of same type");

	py::class_<G3D, std::shared_ptr<G3D>>(m, ("Grid3D" + suffix).c_str())
		.def(py::init<py::array_t<T>, size_t, size_t, size_t>(),
				 py::arg("array"),
				 py::arg("depth") = 0,
				 py::arg("rows") = 0,
				 py::arg("cols") = 0)
		.def(
			"__call__",
			[](G3D& g, size_t d, size_t r, size_t c) -> T& { return g(d, r, c); },
			py::return_value_policy::reference_internal)
		.def("__getitem__", [](G3D& g, size_t i) { return g[i]; })
		.def("__setitem__", [](G3D& g, size_t i, T val) { g[i] = val; })
		.def_property_readonly("depth", &G3D::depth)
		.def_property_readonly("rows", &G3D::rows)
		.def_property_readonly("cols", &G3D::cols)
		.def_property_readonly("size", &G3D::size)
		.def("as_numpy", &G3D::as_numpy_3d)
		.def("get_slice", &G3D::get_slice, "Get 2D slice at specific depth")
		.def("set_slice", &G3D::set_slice, "Set 2D slice at specific depth");

	// Bind utility functions
	m.def(
		("compute_sum_" + suffix).c_str(), &compute_sum<T>, "Compute sum of array");
	m.def(("compute_grid_average_" + suffix).c_str(),
				&compute_grid_average<T>,
				"Compute average of grid");
}

// ==============================================
// CORE ENUM BINDINGS
// ==============================================

inline void
bind_bc_enums(py::module& m)
{
	// NodeType enum
	py::enum_<NodeType>(m, "NodeType")
		.value("NO_DATA", NodeType::NO_DATA, "Invalid or uninitialized node")
		.value("NORMAL", NodeType::NORMAL, "Standard internal node")
		.value("HAS_TO_OUT", NodeType::HAS_TO_OUT, "Forced outlet node")
		.value("CAN_OUT", NodeType::CAN_OUT, "Optional outlet node")
		.value("IN", NodeType::IN, "Inlet/source node")
		.value("PERIODIC", NodeType::PERIODIC, "Periodic boundary condition")
		.value("REFLECT", NodeType::REFLECT, "Reflective boundary condition")
		.export_values();

	// Direction enum
	py::enum_<Direction>(m, "Direction")
		.value("NORTH", Direction::NORTH)
		.value("EAST", Direction::EAST)
		.value("SOUTH", Direction::SOUTH)
		.value("WEST", Direction::WEST)
		.value("NORTHEAST", Direction::NORTHEAST)
		.value("SOUTHEAST", Direction::SOUTHEAST)
		.value("SOUTHWEST", Direction::SOUTHWEST)
		.value("NORTHWEST", Direction::NORTHWEST)
		.value("CENTER", Direction::CENTER)
		.value("INVALID", Direction::INVALID)
		.export_values();

	// ConnectivityType enum
	py::enum_<ConnectivityType>(m, "ConnectivityType")
		.value("D4", ConnectivityType::D4, "4-connected (cardinal directions only)")
		.value("D8", ConnectivityType::D8, "8-connected (cardinal + diagonal)")
		.export_values();
}

// ==============================================
// NODETYPE UTILITIES BINDINGS
// ==============================================

inline void
bind_node_type_utils(py::module& m)
{
	py::class_<NodeTypeUtils>(m, "NodeTypeUtils")
		.def_static(
			"to_string", &NodeTypeUtils::to_string, "Convert NodeType to string")
		.def_static(
			"from_string", &NodeTypeUtils::from_string, "Create NodeType from string")
		.def_static(
			"is_active", &NodeTypeUtils::is_active, "Check if node type is active")
		.def_static("is_boundary",
								&NodeTypeUtils::is_boundary,
								"Check if node type is boundary")
		.def_static("allows_outflow",
								&NodeTypeUtils::allows_outflow,
								"Check if allows outflow")
		.def_static(
			"allows_inflow", &NodeTypeUtils::allows_inflow, "Check if allows inflow")
		.def_static("get_description",
								&NodeTypeUtils::get_description,
								"Get detailed description");
}
// ==============================================
// NEIGHBOR STRUCT BINDING
// ==============================================

inline void
bind_neighbor(py::module& m)
{
	py::class_<Neighbor>(m, "Neighbor")
		.def(py::init<>())
		.def(py::init<size_t, size_t, size_t, Direction, double, NodeType, bool>(),
				 py::arg("index"),
				 py::arg("row"),
				 py::arg("col"),
				 py::arg("direction"),
				 py::arg("distance"),
				 py::arg("boundary_type"),
				 py::arg("is_valid"))
		.def_readwrite("index", &Neighbor::index)
		.def_readwrite("row", &Neighbor::row)
		.def_readwrite("col", &Neighbor::col)
		.def_readwrite("direction", &Neighbor::direction)
		.def_readwrite("distance", &Neighbor::distance)
		.def_readwrite("boundary_type", &Neighbor::boundary_type)
		.def_readwrite("is_valid", &Neighbor::is_valid);
}

// ==============================================
// BCBUILDER BINDINGS
// ==============================================

template<typename T>
void
bind_bc_builder(py::module& m, const std::string& suffix)
{
	using BCBuilderType = BCBuilder<T>;

	py::class_<BCBuilderType>(m, ("BCBuilder" + suffix).c_str())
		// Constructors
		.def(py::init<size_t, size_t>(), py::arg("rows"), py::arg("cols"))

		// Basic properties
		.def("rows", &BCBuilderType::rows)
		.def("cols", &BCBuilderType::cols)
		.def("size", &BCBuilderType::size)

		// Grid access
		.def("get",
				 py::overload_cast<size_t, size_t>(&BCBuilderType::get, py::const_))
		.def("set",
				 py::overload_cast<size_t, size_t, NodeType>(&BCBuilderType::set))
		.def("get_grid", &BCBuilderType::get_grid)

		// Basic patterns - method chaining
		.def("fill", &BCBuilderType::fill, py::arg("type"))
		.def("to_numpy_grid", &BCBuilderType::to_numpy_grid)

		.def("set_borders", &BCBuilderType::set_borders, py::arg("type"))

		// Border-specific setters
		.def("set_north_border", &BCBuilderType::set_north_border, py::arg("type"))
		.def("set_south_border", &BCBuilderType::set_south_border, py::arg("type"))
		.def("set_east_border", &BCBuilderType::set_east_border, py::arg("type"))
		.def("set_west_border", &BCBuilderType::set_west_border, py::arg("type"))

		// Preset patterns
		.def("open_borders", &BCBuilderType::open_borders)
		.def("closed_borders", &BCBuilderType::closed_borders)
		.def("periodic_borders", &BCBuilderType::periodic_borders)
		.def("ns_periodic_ew_closed", &BCBuilderType::ns_periodic_ew_closed)
		.def("ns_periodic_ew_open", &BCBuilderType::ns_periodic_ew_open)
		.def("ew_periodic_ns_closed", &BCBuilderType::ew_periodic_ns_closed)
		.def("ew_periodic_ns_open", &BCBuilderType::ew_periodic_ns_open)
		.def("inlet_north_outlet_south", &BCBuilderType::inlet_north_outlet_south)
		.def("inlet_west_outlet_east", &BCBuilderType::inlet_west_outlet_east)
		.def("corner_outlets", &BCBuilderType::corner_outlets)
		.def("center_outlet", &BCBuilderType::center_outlet)
		.def("random_outlets",
				 &BCBuilderType::random_outlets,
				 py::arg("num_outlets"),
				 py::arg("seed") = 0)

		// Elevation-based methods
		.def("lowest_elevation_outlets",
				 &BCBuilderType::lowest_elevation_outlets,
				 py::arg("elevation"),
				 py::arg("num_outlets") = 1)
		.def("highest_elevation_inlets",
				 &BCBuilderType::highest_elevation_inlets,
				 py::arg("elevation"),
				 py::arg("num_inlets") = 1)
		.def("set_sea_level",
				 &BCBuilderType::set_sea_level,
				 py::arg("elevation"),
				 py::arg("sea_level"))

		// Advanced patterns
		.def("expand_boundary_type",
				 &BCBuilderType::expand_boundary_type,
				 py::arg("source_type"),
				 py::arg("target_type"),
				 py::arg("expansion_cells"))

		// Update methods
		.def("update_add_outlets",
				 &BCBuilderType::update_add_outlets,
				 py::arg("elevation"),
				 py::arg("additional_outlets"))
		.def("update_add_inlets",
				 &BCBuilderType::update_add_inlets,
				 py::arg("elevation"),
				 py::arg("additional_inlets"))
		.def("update_selective",
				 &BCBuilderType::update_selective,
				 py::arg("target_type"),
				 py::arg("new_type"))
		.def("update_by_elevation",
				 &BCBuilderType::update_by_elevation,
				 py::arg("elevation"),
				 py::arg("min_elev"),
				 py::arg("max_elev"),
				 py::arg("new_type"))
		.def("update_conditional",
				 &BCBuilderType::update_conditional,
				 py::arg("predicate"),
				 py::arg("new_type"))
		.def("update_by_distance",
				 &BCBuilderType::update_by_distance,
				 py::arg("reference_points"),
				 py::arg("max_distance"),
				 py::arg("new_type"))

		// Analysis methods
		.def("count_boundary_types", &BCBuilderType::count_boundary_types)
		.def("validate", &BCBuilderType::validate)
		.def("get_summary", &BCBuilderType::get_summary)

		// Connector creation
		.def("create_connector",
				 &BCBuilderType::create_connector,
				 py::arg("connectivity") = ConnectivityType::D8)
		.def("create_connectors", &BCBuilderType::create_connectors)

		// Serialization
		.def("export_as_uint8", &BCBuilderType::export_as_uint8)
		.def(
			"import_from_uint8", &BCBuilderType::import_from_uint8, py::arg("data"))
		.def("clone", &BCBuilderType::clone);
}

// ==============================================
// QUICKBC BINDINGS
// ==============================================

template<typename T>
void
bind_quick_bc(py::module& m, const std::string& suffix)
{
	using QuickBCType = QuickBC<T>;

	py::class_<QuickBCType>(m, ("QuickBC" + suffix).c_str())
		.def_static("open_domain",
								&QuickBCType::open_domain,
								py::arg("rows"),
								py::arg("cols"))
		.def_static("closed_domain",
								&QuickBCType::closed_domain,
								py::arg("rows"),
								py::arg("cols"))
		.def_static("periodic_ew_domain",
								&QuickBCType::periodic_ew_domain,
								py::arg("rows"),
								py::arg("cols"))
		.def_static("periodic_ns_domain",
								&QuickBCType::periodic_ns_domain,
								py::arg("rows"),
								py::arg("cols"))
		.def_static("sea_level_domain",
								&QuickBCType::sea_level_domain,
								py::arg("rows"),
								py::arg("cols"),
								py::arg("elevation"),
								py::arg("sea_level"));
}

// ==============================================
// CONNECTOR BINDINGS
// ==============================================

template<typename T>
void
bind_connector(py::module& m, const std::string& suffix)
{
	using ConnectorType = Connector<T>;
	using FlowTransferType = typename ConnectorType::FlowTransfer;

	// FlowTransfer struct
	py::class_<FlowTransferType>(m, ("FlowTransfer" + suffix).c_str())
		.def(py::init<>())
		.def(py::init<size_t, T, bool, bool, bool>(),
				 py::arg("dest"),
				 py::arg("mult"),
				 py::arg("periodic") = false,
				 py::arg("reflect") = false,
				 py::arg("exits") = false)
		.def_readwrite("destination_index", &FlowTransferType::destination_index)
		.def_readwrite("flux_multiplier", &FlowTransferType::flux_multiplier)
		.def_readwrite("is_periodic_wrap", &FlowTransferType::is_periodic_wrap)
		.def_readwrite("is_reflection", &FlowTransferType::is_reflection)
		.def_readwrite("exits_domain", &FlowTransferType::exits_domain);

	// Main Connector class with shared_ptr holder
	py::class_<ConnectorType, std::shared_ptr<ConnectorType>>(
		m, ("Connector" + suffix).c_str())
		// Constructors
		.def(py::init<size_t, size_t, ConnectivityType>(),
				 py::arg("rows"),
				 py::arg("cols"),
				 py::arg("connectivity") = ConnectivityType::D8)
		.def(py::init<size_t,
									size_t,
									std::shared_ptr<Grid2D<NodeType>>,
									ConnectivityType>(),
				 py::arg("rows"),
				 py::arg("cols"),
				 py::arg("boundary_grid"),
				 py::arg("connectivity") = ConnectivityType::D8)

		// Basic properties
		.def("rows", &ConnectorType::rows)
		.def("cols", &ConnectorType::cols)
		.def("size", &ConnectorType::size)
		.def("connectivity_type", &ConnectorType::connectivity_type)
		.def("num_directions", &ConnectorType::num_directions)

		// Index conversion
		.def("to_1d", &ConnectorType::to_1d)
		.def("to_2d", &ConnectorType::to_2d)
		.def("is_valid_coord",
				 py::overload_cast<size_t, size_t>(&ConnectorType::is_valid_coord,
																					 py::const_))
		.def(
			"is_valid_coord",
			py::overload_cast<int, int>(&ConnectorType::is_valid_coord, py::const_))
		.def("is_valid_index", &ConnectorType::is_valid_index)

		// Boundary conditions
		.def("get_boundary_type",
				 py::overload_cast<size_t, size_t>(&ConnectorType::get_boundary_type,
																					 py::const_))
		.def(
			"get_boundary_type",
			py::overload_cast<size_t>(&ConnectorType::get_boundary_type, py::const_))
		.def("set_boundary_type",
				 py::overload_cast<size_t, size_t, NodeType>(
					 &ConnectorType::set_boundary_type))
		.def("set_boundary_type",
				 py::overload_cast<size_t, NodeType>(&ConnectorType::set_boundary_type))
		.def("is_active_node",
				 py::overload_cast<size_t, size_t>(&ConnectorType::is_active_node,
																					 py::const_))
		.def("is_active_node",
				 py::overload_cast<size_t>(&ConnectorType::is_active_node, py::const_))
		.def("is_boundary_node",
				 py::overload_cast<size_t, size_t>(&ConnectorType::is_boundary_node,
																					 py::const_))
		.def(
			"is_boundary_node",
			py::overload_cast<size_t>(&ConnectorType::is_boundary_node, py::const_))

		// Direction utilities
		.def("get_opposite_direction", &ConnectorType::get_opposite_direction)
		.def("get_direction_distance", &ConnectorType::get_direction_distance)
		.def("is_diagonal_direction", &ConnectorType::is_diagonal_direction)
		.def("is_cardinal_direction", &ConnectorType::is_cardinal_direction)

		// Single neighbor access (unchecked - fastest)
		.def("north_unchecked", &ConnectorType::north_unchecked)
		.def("east_unchecked", &ConnectorType::east_unchecked)
		.def("south_unchecked", &ConnectorType::south_unchecked)
		.def("west_unchecked", &ConnectorType::west_unchecked)
		.def("northeast_unchecked", &ConnectorType::northeast_unchecked)
		.def("southeast_unchecked", &ConnectorType::southeast_unchecked)
		.def("southwest_unchecked", &ConnectorType::southwest_unchecked)
		.def("northwest_unchecked", &ConnectorType::northwest_unchecked)

		// Single neighbor access (checked)
		.def("get_neighbor",
				 py::overload_cast<size_t, size_t, Direction>(
					 &ConnectorType::get_neighbor, py::const_))
		.def("get_neighbor",
				 py::overload_cast<size_t, Direction>(&ConnectorType::get_neighbor,
																							py::const_))
		.def("north", py::overload_cast<size_t>(&ConnectorType::north, py::const_))
		.def("east", py::overload_cast<size_t>(&ConnectorType::east, py::const_))
		.def("south", py::overload_cast<size_t>(&ConnectorType::south, py::const_))
		.def("west", py::overload_cast<size_t>(&ConnectorType::west, py::const_))
		.def("northeast",
				 py::overload_cast<size_t>(&ConnectorType::northeast, py::const_))
		.def("southeast",
				 py::overload_cast<size_t>(&ConnectorType::southeast, py::const_))
		.def("southwest",
				 py::overload_cast<size_t>(&ConnectorType::southwest, py::const_))
		.def("northwest",
				 py::overload_cast<size_t>(&ConnectorType::northwest, py::const_))
		.def("north",
				 py::overload_cast<size_t, size_t>(&ConnectorType::north, py::const_))
		.def("east",
				 py::overload_cast<size_t, size_t>(&ConnectorType::east, py::const_))
		.def("south",
				 py::overload_cast<size_t, size_t>(&ConnectorType::south, py::const_))
		.def("west",
				 py::overload_cast<size_t, size_t>(&ConnectorType::west, py::const_))
		.def(
			"northeast",
			py::overload_cast<size_t, size_t>(&ConnectorType::northeast, py::const_))
		.def(
			"southeast",
			py::overload_cast<size_t, size_t>(&ConnectorType::southeast, py::const_))
		.def(
			"southwest",
			py::overload_cast<size_t, size_t>(&ConnectorType::southwest, py::const_))
		.def(
			"northwest",
			py::overload_cast<size_t, size_t>(&ConnectorType::northwest, py::const_))

		// All neighbors access
		.def("get_all_neighbors",
				 py::overload_cast<size_t, size_t>(&ConnectorType::get_all_neighbors,
																					 py::const_))
		.def(
			"get_all_neighbors",
			py::overload_cast<size_t>(&ConnectorType::get_all_neighbors, py::const_))
		.def("get_valid_neighbors",
				 py::overload_cast<size_t, size_t>(&ConnectorType::get_valid_neighbors,
																					 py::const_))
		.def("get_valid_neighbors",
				 py::overload_cast<size_t>(&ConnectorType::get_valid_neighbors,
																	 py::const_))

		// Vectorized operations
		.def("get_neighbor_indices", &ConnectorType::get_neighbor_indices)
		.def("get_valid_neighbor_indices",
				 &ConnectorType::get_valid_neighbor_indices)

		// Boundary condition specific methods
		.def("apply_periodic_bc", &ConnectorType::apply_periodic_bc)
		.def("get_boundary_nodes", &ConnectorType::get_boundary_nodes)
		.def("get_all_boundary_nodes", &ConnectorType::get_all_boundary_nodes)
		.def("set_border_boundary", &ConnectorType::set_border_boundary)
		.def("set_periodic_boundaries", &ConnectorType::set_periodic_boundaries)
		.def("set_reflective_boundaries", &ConnectorType::set_reflective_boundaries)

		// Effective neighbor access (with BC handling)
		.def("get_effective_neighbor", &ConnectorType::get_effective_neighbor)
		.def("get_effective_neighbors",
				 py::overload_cast<size_t, size_t>(
					 &ConnectorType::get_effective_neighbors, py::const_))
		.def("get_effective_neighbors",
				 py::overload_cast<size_t>(&ConnectorType::get_effective_neighbors,
																	 py::const_))
		.def("get_effective_valid_neighbors",
				 py::overload_cast<size_t, size_t>(
					 &ConnectorType::get_effective_valid_neighbors, py::const_))
		.def("get_effective_valid_neighbors",
				 py::overload_cast<size_t>(
					 &ConnectorType::get_effective_valid_neighbors, py::const_))
		.def("has_boundary_handling", &ConnectorType::has_boundary_handling)

		// Advanced flow handling
		.def("transfer_flow", &ConnectorType::transfer_flow);
}

// ==============================================
// CONVENIENCE BINDING FUNCTIONS
// ==============================================

template<typename T>
void
bind_all_bc_types(py::module& m, const std::string& suffix)
{
	bind_bc_builder<T>(m, suffix);
	bind_quick_bc<T>(m, suffix);
	bind_connector<T>(m, suffix);
	// bind_all_fast_connector_types<T>(m, suffix);
}

// Main binding function for all BC functionality
inline void
bind_arrbcconn(py::module& m)
{
	// Bind core enums and utilities
	bind_bc_enums(m);
	bind_node_type_utils(m);
	bind_neighbor(m);

	// Bind the template classes for different types
	bind_array_types<float>(m, "F32");
	bind_array_types<double>(m, "F64");
	bind_array_types<int32_t>(m, "I32");
	bind_array_types<int64_t>(m, "I64");
	bind_array_types<size_t>(m, "U64");
	bind_array_types<uint8_t>(m, "U8");

	// Add this in bind_arrbcconn function, after binding the enums:
	py::class_<Grid2D<NodeType>, std::shared_ptr<Grid2D<NodeType>>>(
		m, "Grid2DNodeType")
		.def(py::init<py::array_t<uint8_t>, size_t, size_t>())
		.def(
			"__call__",
			[](Grid2D<NodeType>& g, size_t r, size_t c) -> NodeType& {
				return g(r, c);
			},
			py::return_value_policy::reference_internal)
		.def("__getitem__",
				 [](Grid2D<NodeType>& g, size_t i) -> NodeType { return g[i]; })
		.def("__setitem__",
				 [](Grid2D<NodeType>& g, size_t i, NodeType val) { g[i] = val; })
		.def_property_readonly("rows", &Grid2D<NodeType>::rows)
		.def_property_readonly("cols", &Grid2D<NodeType>::cols)
		.def_property_readonly("size", &Grid2D<NodeType>::size);

	// Bind utility functions
	m.def("compute_sum_f32", &compute_sum<float>, "Compute sum of float array");
	m.def("compute_sum_f64", &compute_sum<double>, "Compute sum of double array");
	m.def("compute_sum_i32", &compute_sum<int32_t>, "Compute sum of int32 array");
	m.def("compute_sum_i64", &compute_sum<int64_t>, "Compute sum of int64 array");

	m.def("compute_grid_average_f32",
				&compute_grid_average<float>,
				"Compute average of float grid");
	m.def("compute_grid_average_f64",
				&compute_grid_average<double>,
				"Compute average of double grid");

	m.def("create_test_grid",
				&create_test_grid,
				"Create a test grid with random values",
				py::arg("rows"),
				py::arg("cols"));

	m.def("apply_smoothing_filter",
				&apply_smoothing_filter,
				"Apply 3x3 smoothing filter to grid",
				py::arg("grid"));

	// Bind templated classes for common types
	bind_all_bc_types<float>(m, "F32");
	bind_all_bc_types<double>(m, "F64");
	bind_all_bc_types<int32_t>(m, "I32");
	bind_all_bc_types<int64_t>(m, "I64");
}

} // namespace dagger2
