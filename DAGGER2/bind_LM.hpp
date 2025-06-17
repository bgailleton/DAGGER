#pragma once

#include "dg2_fastconnector.hpp"
#include "dg2_priority_flood.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace dagger2 {

// ==============================================
// PRIORITY FLOOD CONFIG BINDINGS
// ==============================================

template<typename T>
void
bind_priority_flood_config(py::module& m, const std::string& suffix)
{
	using ConfigType = typename PriorityFlood<T>::Config;

	py::class_<ConfigType>(m, ("PriorityFloodConfig" + suffix).c_str())
		.def(py::init<>())
		.def_readwrite("epsilon",
									 &ConfigType::epsilon,
									 "Minimum elevation increment for filling")
		.def_readwrite("preserve_flat_areas",
									 &ConfigType::preserve_flat_areas,
									 "Keep flat areas flat vs add tiny gradient")
		.def_readwrite("validate_result",
									 &ConfigType::validate_result,
									 "Run post-processing validation")
		.def_readwrite("generate_statistics",
									 &ConfigType::generate_statistics,
									 "Collect detailed statistics")
		.def_readwrite("max_fill_depth",
									 &ConfigType::max_fill_depth,
									 "Maximum allowed fill depth")
		.def_readwrite("max_iterations",
									 &ConfigType::max_iterations,
									 "Safety limit on iterations");
}

// ==============================================
// PRIORITY FLOOD RESULTS BINDINGS
// ==============================================

template<typename T>
void
bind_priority_flood_results(py::module& m, const std::string& suffix)
{
	using ResultsType = typename PriorityFlood<T>::Results;

	py::class_<ResultsType>(m, ("PriorityFloodResults" + suffix).c_str())
		.def(py::init<>())
		.def_readonly(
			"nodes_processed", &ResultsType::nodes_processed, "Total nodes processed")
		.def_readonly("nodes_filled",
									&ResultsType::nodes_filled,
									"Nodes that had elevation increased")
		.def_readonly("total_fill_volume",
									&ResultsType::total_fill_volume,
									"Total volume of material added")
		.def_readonly("max_fill_depth",
									&ResultsType::max_fill_depth,
									"Maximum fill depth encountered")
		.def_readonly(
			"iterations", &ResultsType::iterations, "Number of algorithm iterations")
		.def_readonly("convergence_reached",
									&ResultsType::convergence_reached,
									"Whether algorithm converged normally")
		.def_readonly("warnings", &ResultsType::warnings, "Any warnings generated")
		.def("__repr__", [](const ResultsType& r) {
			return "<PriorityFloodResults processed=" +
						 std::to_string(r.nodes_processed) +
						 " filled=" + std::to_string(r.nodes_filled) +
						 " volume=" + std::to_string(r.total_fill_volume) + ">";
		});
}

// ==============================================
// PRIORITY FLOOD MAIN CLASS BINDINGS
// ==============================================

template<typename T>
void
bind_priority_flood(py::module& m, const std::string& suffix)
{
	using PriorityFloodType = PriorityFlood<T>;
	using ConfigType = typename PriorityFloodType::Config;
	using ResultsType = typename PriorityFloodType::Results;

	py::class_<PriorityFloodType>(m, ("PriorityFlood" + suffix).c_str())
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation, const Connector<T>& connector) {
				return PriorityFloodType::fill_depressions(
					elevation, connector, typename PriorityFloodType::Config{});
			},
			"Fill depressions with default config",
			py::arg("elevation"),
			py::arg("connector"))
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation, const Connector<T>& connector, T epsilon) {
				return PriorityFloodType::fill_depressions(
					elevation, connector, epsilon);
			},
			"Fill depressions with custom epsilon",
			py::arg("elevation"),
			py::arg("connector"),
			py::arg("epsilon"))
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation,
				 const Connector<T>& connector,
				 const ConfigType& config) {
				return PriorityFloodType::fill_depressions(
					elevation, connector, config);
			},
			"Fill depressions with full config",
			py::arg("elevation"),
			py::arg("connector"),
			py::arg("config"))
		.def_static("backup_elevation",
								&PriorityFloodType::backup_elevation,
								"Create backup of elevation grid",
								py::arg("elevation"))
		.def_static("restore_elevation",
								&PriorityFloodType::restore_elevation,
								"Restore elevation from backup",
								py::arg("elevation"),
								py::arg("backup"))
		.def_static("generate_report",
								&PriorityFloodType::generate_report,
								"Generate detailed results report",
								py::arg("results"));
}

// ==============================================
// PRIORITY FLOOD FASTCONNECTOR BINDINGS
// ==============================================

template<typename T>
void
bind_priority_flood_fast(py::module& m, const std::string& suffix)
{
	using PriorityFloodType = PriorityFlood<T>;
	using ConfigType = typename PriorityFloodType::Config;
	using FastConnectorD4 = FastConnector<T, ConnectivityType::D4>;
	using FastConnectorD8 = FastConnector<T, ConnectivityType::D8>;

	// Create dummy classes that don't conflict
	struct PriorityFloodFastD4
	{};
	struct PriorityFloodFastD8
	{};

	// FastConnector D4 bindings
	py::class_<PriorityFloodFastD4>(m, ("PriorityFloodFastD4" + suffix).c_str())
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation, const FastConnectorD4& connector) {
				return PriorityFloodType::template fill_depressions<FastConnectorD4>(
					elevation, connector, typename PriorityFloodType::Config{});
			},
			"Fill depressions with FastConnector D4 (default config)",
			py::arg("elevation"),
			py::arg("connector"))
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation, const FastConnectorD4& connector, T epsilon) {
				return PriorityFloodType::template fill_depressions<FastConnectorD4>(
					elevation, connector, epsilon);
			},
			"Fill depressions with FastConnector D4 (custom epsilon)",
			py::arg("elevation"),
			py::arg("connector"),
			py::arg("epsilon"))
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation,
				 const FastConnectorD4& connector,
				 const ConfigType& config) {
				return PriorityFloodType::template fill_depressions<FastConnectorD4>(
					elevation, connector, config);
			},
			"Fill depressions with FastConnector D4 (full config)",
			py::arg("elevation"),
			py::arg("connector"),
			py::arg("config"));

	// FastConnector D8 bindings
	py::class_<PriorityFloodFastD8>(m, ("PriorityFloodFastD8" + suffix).c_str())
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation, const FastConnectorD8& connector) {
				return PriorityFloodType::template fill_depressions<FastConnectorD8>(
					elevation, connector, typename PriorityFloodType::Config{});
			},
			"Fill depressions with FastConnector D8 (default config)",
			py::arg("elevation"),
			py::arg("connector"))
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation, const FastConnectorD8& connector, T epsilon) {
				return PriorityFloodType::template fill_depressions<FastConnectorD8>(
					elevation, connector, epsilon);
			},
			"Fill depressions with FastConnector D8 (custom epsilon)",
			py::arg("elevation"),
			py::arg("connector"),
			py::arg("epsilon"))
		.def_static(
			"fill_depressions",
			[](Grid2D<T>& elevation,
				 const FastConnectorD8& connector,
				 const ConfigType& config) {
				return PriorityFloodType::template fill_depressions<FastConnectorD8>(
					elevation, connector, config);
			},
			"Fill depressions with FastConnector D8 (full config)",
			py::arg("elevation"),
			py::arg("connector"),
			py::arg("config"));
}

// ==============================================
// CONVENIENCE BINDING FUNCTIONS
// ==============================================

template<typename T>
void
bind_all_priority_flood_types(py::module& m, const std::string& suffix)
{
	bind_priority_flood_config<T>(m, suffix);
	bind_priority_flood_results<T>(m, suffix);
	bind_priority_flood<T>(m, suffix);
	bind_priority_flood_fast<T>(m, suffix);
}

// ==============================================
// MAIN BINDING FUNCTION
// ==============================================

inline void
bind_priority_flood_module(py::module& m)
{
	// Bind template classes for different types
	bind_all_priority_flood_types<float>(m, "F32");
	bind_all_priority_flood_types<double>(m, "F64");

	// Module documentation
	m.doc() = R"pbdoc(
        🏔️ DAGGER2 PRIORITY FLOOD MODULE 🏔️
        ===================================

        Advanced depression filling using Priority Flood + Epsilon algorithm.
        Handles all boundary conditions and provides comprehensive statistics.

        Supports both regular Connector and optimized FastConnector interfaces.
    )pbdoc";
}

} // namespace dagger2
