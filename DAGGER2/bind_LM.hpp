#pragma once

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
									 "Safety limit on iterations")
		.def_readwrite("use_optimized",
									 &ConfigType::use_optimized,
									 "Use ultra-optimized raw index version");
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
		.def_readonly("depressions_filled",
									&ResultsType::depressions_filled,
									"Number of separate depressions identified")
		.def_readonly("total_fill_volume",
									&ResultsType::total_fill_volume,
									"Total volume of material added")
		.def_readonly("max_fill_depth",
									&ResultsType::max_fill_depth,
									"Maximum fill depth encountered")
		.def_readonly("min_elevation",
									&ResultsType::min_elevation,
									"Minimum elevation after filling")
		.def_readonly("max_elevation",
									&ResultsType::max_elevation,
									"Maximum elevation after filling")
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
		.def_static("fill_depressions",
								py::overload_cast<Grid2D<T>&, const Connector<T>&>(
									&PriorityFloodType::fill_depressions),
								"Fill depressions with default config",
								py::arg("elevation"),
								py::arg("connector"))
		.def_static("fill_depressions",
								py::overload_cast<Grid2D<T>&, const Connector<T>&, T>(
									&PriorityFloodType::fill_depressions),
								"Fill depressions with custom epsilon",
								py::arg("elevation"),
								py::arg("connector"),
								py::arg("epsilon"))
		.def_static("fill_depressions",
								py::overload_cast<Grid2D<T>&, const Connector<T>&, bool>(
									&PriorityFloodType::fill_depressions),
								"Fill depressions with optimization flag",
								py::arg("elevation"),
								py::arg("connector"),
								py::arg("use_optimized"))
		.def_static(
			"fill_depressions",
			py::overload_cast<Grid2D<T>&, const Connector<T>&, const ConfigType&>(
				&PriorityFloodType::fill_depressions),
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
		.def_static("has_proper_drainage",
								&PriorityFloodType::has_proper_drainage,
								"Check if grid has proper drainage",
								py::arg("elevation"),
								py::arg("connector"),
								py::arg("use_optimized") = true)
		.def_static("find_sinks",
								&PriorityFloodType::find_sinks,
								"Find remaining sink nodes",
								py::arg("elevation"),
								py::arg("connector"),
								py::arg("use_optimized") = true)
		.def_static("estimate_epsilon",
								&PriorityFloodType::estimate_epsilon,
								"Estimate appropriate epsilon value",
								py::arg("elevation"),
								py::arg("connector"),
								py::arg("use_optimized") = true)
		.def_static("generate_report",
								&PriorityFloodType::generate_report,
								"Generate detailed results report",
								py::arg("results"));
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

        Features:
        - Ultra-optimized raw index version (use_optimized=True)
        - Original Neighbor structure version (use_optimized=False)
        - Full boundary condition support (PERIODIC, REFLECT, etc.)
        - Detailed statistics and validation
    )pbdoc";
}

} // namespace dagger2
