#pragma once

#include "dg2_cordonnier_flow_rerouter.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace dagger2 {

// ==============================================
// ENUM BINDINGS
// ==============================================

template<typename T>
void
bind_cordonnier_enums(py::module& m, const std::string& suffix)
{
	// Enforcement Strategy enum
	using yolo_T1 = typename CordonnierFlowRerouter<T>::EnforcementStrategy;
	py::enum_<yolo_T1>(m,
										 ("EnforcementStrategy" + suffix).c_str(),
										 "Flow enforcement strategy within depressions")
		.value("SIMPLE_CORRECTION",
					 yolo_T1::SIMPLE_CORRECTION,
					 "Only update local minima receivers (minimal change)")
		.value("DEPRESSION_CARVING",
					 yolo_T1::DEPRESSION_CARVING,
					 "Create narrow channels through depressions")
		.value("DEPRESSION_FILLING",
					 yolo_T1::DEPRESSION_FILLING,
					 "Route flow as if depressions were filled with water")
		.export_values();

	// MST Algorithm enum
	using yolo_T2 = typename CordonnierFlowRerouter<T>::MSTAlgorithm;
	py::enum_<yolo_T2>(m,
										 ("MSTAlgorithm" + suffix).c_str(),
										 "Minimum spanning tree algorithm choice")
		.value("KRUSKAL_LOGN",
					 yolo_T2::KRUSKAL_LOGN,
					 "O(n log n) complexity using Kruskal's algorithm")
		.value("PLANAR_LINEAR",
					 yolo_T2::PLANAR_LINEAR,
					 "O(n) complexity leveraging planar graph properties")
		.export_values();
}

// ==============================================
// CONFIG BINDINGS
// ==============================================

template<typename T>
void
bind_cordonnier_config(py::module& m, const std::string& suffix)
{
	using ConfigType = typename CordonnierFlowRerouter<T>::Config;
	using EnforcementStrategy =
		typename CordonnierFlowRerouter<T>::EnforcementStrategy;
	using MSTAlgorithm = typename CordonnierFlowRerouter<T>::MSTAlgorithm;

	py::class_<ConfigType>(m, ("CordonnierConfig" + suffix).c_str())
		.def(py::init<>())
		.def_readwrite("strategy",
									 &ConfigType::strategy,
									 "Flow enforcement strategy within depressions")
		.def_readwrite("mst_algorithm",
									 &ConfigType::mst_algorithm,
									 "Minimum spanning tree algorithm choice")
		.def_readwrite("min_gradient_threshold",
									 &ConfigType::min_gradient_threshold,
									 "Minimum gradient threshold for flow detection")
		.def_readwrite("preserve_original_receivers",
									 &ConfigType::preserve_original_receivers,
									 "Keep backup of original receivers")
		.def_readwrite("generate_statistics",
									 &ConfigType::generate_statistics,
									 "Collect detailed computation statistics")
		.def_readwrite("validate_result",
									 &ConfigType::validate_result,
									 "Run post-processing validation checks")
		.def("__repr__", [](const ConfigType& config) {
			std::ostringstream oss;
			oss << "CordonnierConfig(strategy=";
			switch (config.strategy) {
				case EnforcementStrategy::SIMPLE_CORRECTION:
					oss << "SIMPLE_CORRECTION";
					break;
				case EnforcementStrategy::DEPRESSION_CARVING:
					oss << "DEPRESSION_CARVING";
					break;
				case EnforcementStrategy::DEPRESSION_FILLING:
					oss << "DEPRESSION_FILLING";
					break;
			}
			oss << ", mst_algorithm=";
			switch (config.mst_algorithm) {
				case MSTAlgorithm::KRUSKAL_LOGN:
					oss << "KRUSKAL_LOGN";
					break;
				case MSTAlgorithm::PLANAR_LINEAR:
					oss << "PLANAR_LINEAR";
					break;
			}
			oss << ", min_gradient_threshold=" << config.min_gradient_threshold
					<< ", validate_result=" << (config.validate_result ? "True" : "False")
					<< ")";
			return oss.str();
		});
}

// ==============================================
// RESULTS BINDINGS
// ==============================================

template<typename T>
void
bind_cordonnier_results(py::module& m, const std::string& suffix)
{
	using ResultsType = typename CordonnierFlowRerouter<T>::Results;

	py::class_<ResultsType>(m, ("CordonnierResults" + suffix).c_str())
		.def(py::init<>())
		.def_readonly("num_basins_processed",
									&ResultsType::num_basins_processed,
									"Total number of basins processed")
		.def_readonly("num_local_minima",
									&ResultsType::num_local_minima,
									"Number of local minima found")
		.def_readonly("num_inner_basins",
									&ResultsType::num_inner_basins,
									"Number of inner (non-boundary) basins")
		.def_readonly("num_boundary_basins",
									&ResultsType::num_boundary_basins,
									"Number of boundary basins")
		.def_readonly("num_links_created",
									&ResultsType::num_links_created,
									"Number of basin links created")
		.def_readonly("num_receivers_updated",
									&ResultsType::num_receivers_updated,
									"Number of flow receivers updated")
		.def_readonly("total_energy_minimized",
									&ResultsType::total_energy_minimized,
									"Total energy minimized in MST")
		.def_readonly("computation_time_ms",
									&ResultsType::computation_time_ms,
									"Total computation time in milliseconds")
		.def_readonly("convergence_achieved",
									&ResultsType::convergence_achieved,
									"Whether algorithm converged successfully")
		.def_readonly(
			"warnings", &ResultsType::warnings, "List of warning messages")
		.def("__repr__",
				 [](const ResultsType& results) {
					 std::ostringstream oss;
					 oss << "CordonnierResults(\n"
							 << "  basins_processed=" << results.num_basins_processed << ",\n"
							 << "  local_minima=" << results.num_local_minima << ",\n"
							 << "  receivers_updated=" << results.num_receivers_updated
							 << ",\n"
							 << "  computation_time=" << results.computation_time_ms
							 << "ms,\n"
							 << "  convergence="
							 << (results.convergence_achieved ? "True" : "False");
					 if (!results.warnings.empty()) {
						 oss << ",\n  warnings=" << results.warnings.size();
					 }
					 oss << "\n)";
					 return oss.str();
				 })
		.def(
			"summary",
			[](const ResultsType& results) {
				std::ostringstream oss;
				oss << "🌊 CORDONNIER FLOW REROUTER RESULTS 🌊\n"
						<< "======================================\n"
						<< "Basins processed: " << results.num_basins_processed << "\n"
						<< "Local minima found: " << results.num_local_minima << "\n"
						<< "Inner basins: " << results.num_inner_basins << "\n"
						<< "Boundary basins: " << results.num_boundary_basins << "\n"
						<< "Basin links created: " << results.num_links_created << "\n"
						<< "Receivers updated: " << results.num_receivers_updated << "\n"
						<< "Total energy minimized: " << results.total_energy_minimized
						<< "\n"
						<< "Computation time: " << results.computation_time_ms << " ms\n"
						<< "Convergence: "
						<< (results.convergence_achieved ? "✅ Success" : "❌ Failed")
						<< "\n";

				if (!results.warnings.empty()) {
					oss << "\nWarnings ⚠️:\n";
					for (const auto& warning : results.warnings) {
						oss << "  • " << warning << "\n";
					}
				}

				return oss.str();
			},
			"Generate detailed summary report");
}

// ==============================================
// MAIN CLASS BINDINGS
// ==============================================

template<typename T>
void
bind_cordonnier_flow_rerouter(py::module& m, const std::string& suffix)
{
	using CordonnierType = CordonnierFlowRerouter<T>;
	using ConfigType = typename CordonnierType::Config;
	using ResultsType = typename CordonnierType::Results;

	py::class_<CordonnierType>(m, ("CordonnierFlowRerouter" + suffix).c_str())
		// Static methods for flow rerouting
		.def_static("reroute_flow",
								&CordonnierType::reroute_flow,
								R"pbdoc(
                    Execute flow rerouting using Cordonnier et al. 2019 algorithm

                    This method implements the linear complexity algorithm for flow routing
                    in topographies with depressions. It explicitly computes flow paths
                    both within and across depressions through basin graph construction.

                    Parameters
                    ----------
                    elevation : Grid2D
                        Input topographic elevation data
                    receivers : ArrayRef
                        Flow receiver array to be updated (modified in-place)
                    connector : Connector
                        Grid connectivity and boundary condition handler
                    config : CordonnierConfig, optional
                        Configuration parameters for the algorithm

                    Returns
                    -------
                    CordonnierResults
                        Detailed results and statistics from the computation

                    Notes
                    -----
                    - Time complexity: O(n) for planar graphs, O(n log n) for Kruskal's MST
                    - Space complexity: O(n) where n is the number of grid nodes
                    - Handles all boundary conditions through the Connector
                    - Multiple enforcement strategies available for different use cases

                    Examples
                    --------
                    >>> import dg2
                    >>> elevation = dg2.Grid2DF64(elevation_data, rows, cols)
                    >>> receivers = dg2.ArrayRefU64(np.zeros(elevation.size, dtype=np.uint64))
                    >>> connector = dg2.ConnectorF64(rows, cols)
                    >>>
                    >>> config = dg2.CordonnierConfigF64()
                    >>> config.strategy = dg2.EnforcementStrategy.DEPRESSION_FILLING
                    >>> config.mst_algorithm = dg2.MSTAlgorithm.PLANAR_LINEAR
                    >>>
                    >>> results = dg2.CordonnierFlowRerouterF64.reroute_flow(
                    ...     elevation, receivers, connector, config)
                    >>> print(results.summary())
                    )pbdoc",
								py::arg("elevation"),
								py::arg("receivers"),
								py::arg("connector"),
								py::arg("config") = ConfigType())

		// Convenience methods with different parameter combinations
		.def_static(
			"reroute_flow_simple",
			[](Grid2D<T>& elevation,
				 ArrayRef<size_t>& receivers,
				 const Connector<T>& connector,
				 typename CordonnierType::EnforcementStrategy strategy) {
				ConfigType config;
				config.strategy = strategy;
				return CordonnierType::reroute_flow(
					elevation, receivers, connector, config);
			},
			"Reroute flow with specified enforcement strategy",
			py::arg("elevation"),
			py::arg("receivers"),
			py::arg("connector"),
			py::arg("strategy"))

		.def_static(
			"reroute_flow_fast",
			[](Grid2D<T>& elevation,
				 ArrayRef<size_t>& receivers,
				 const Connector<T>& connector) {
				ConfigType config;
				config.strategy =
					CordonnierType::EnforcementStrategy::DEPRESSION_FILLING;
				config.mst_algorithm = CordonnierType::MSTAlgorithm::PLANAR_LINEAR;
				config.validate_result = false;
				return CordonnierType::reroute_flow(
					elevation, receivers, connector, config);
			},
			"Fast rerouting with optimal settings (no validation)",
			py::arg("elevation"),
			py::arg("receivers"),
			py::arg("connector"))

		.def_static(
			"reroute_flow_detailed",
			[](Grid2D<T>& elevation,
				 ArrayRef<size_t>& receivers,
				 const Connector<T>& connector) {
				ConfigType config;
				config.strategy =
					CordonnierType::EnforcementStrategy::DEPRESSION_FILLING;
				config.mst_algorithm = CordonnierType::MSTAlgorithm::PLANAR_LINEAR;
				config.generate_statistics = true;
				config.validate_result = true;
				return CordonnierType::reroute_flow(
					elevation, receivers, connector, config);
			},
			"Detailed rerouting with full statistics and validation",
			py::arg("elevation"),
			py::arg("receivers"),
			py::arg("connector"));
}

// ==============================================
// CONVENIENCE BINDING FUNCTIONS
// ==============================================

template<typename T>
void
bind_all_cordonnier_types(py::module& m, const std::string& suffix)
{
	bind_cordonnier_config<T>(m, suffix);
	bind_cordonnier_results<T>(m, suffix);
	bind_cordonnier_flow_rerouter<T>(m, suffix);
}

// ==============================================
// UTILITY FUNCTIONS
// ==============================================

template<typename T>
void
bind_cordonnier_utilities(py::module& m, const std::string& suffix)
{
	// Helper functions for working with flow networks
	m.def(("create_simple_receivers" + suffix).c_str(),
				[](const Grid2D<T>& elevation, const Connector<T>& connector) {
					std::vector<size_t> receivers(elevation.size());

					// Initialize with steepest descent receivers
					for (size_t i = 0; i < elevation.size(); ++i) {
						if (!connector.is_active_node(i)) {
							receivers[i] = SIZE_MAX;
							continue;
						}

						auto [row, col] = connector.to_2d(i);
						auto neighbors = connector.get_effective_valid_neighbors(row, col);

						size_t best_receiver = SIZE_MAX;
						T steepest_gradient = 0;

						for (const auto& neighbor : neighbors) {
							T gradient =
								(elevation[i] - elevation[neighbor.index]) / neighbor.distance;
							if (gradient > steepest_gradient) {
								steepest_gradient = gradient;
								best_receiver = neighbor.index;
							}
						}

						receivers[i] = best_receiver;
					}

					return std::make_shared<ArrayRef<size_t>>(receivers);
				},
				"Create initial receivers using steepest descent",
				py::arg("elevation"),
				py::arg("connector"));

	m.def(("validate_flow_network" + suffix).c_str(),
				[](const ArrayRef<size_t>& receivers, const Connector<T>& connector) {
					std::vector<std::string> issues;

					// Check for basic validity
					for (size_t i = 0; i < receivers.size(); ++i) {
						if (!connector.is_active_node(i))
							continue;

						size_t receiver = receivers[i];
						if (receiver != SIZE_MAX) {
							if (receiver >= receivers.size()) {
								issues.push_back("Node " + std::to_string(i) +
																 " has out-of-bounds receiver");
							} else if (!connector.is_active_node(receiver)) {
								issues.push_back("Node " + std::to_string(i) +
																 " flows to inactive node");
							}
						}
					}

					// Simple cycle detection
					for (size_t i = 0; i < receivers.size(); ++i) {
						if (!connector.is_active_node(i))
							continue;

						std::unordered_set<size_t> visited;
						size_t current = i;

						while (current != SIZE_MAX &&
									 visited.find(current) == visited.end()) {
							visited.insert(current);
							current = receivers[current];

							if (visited.size() > receivers.size()) {
								issues.push_back("Potential infinite loop from node " +
																 std::to_string(i));
								break;
							}
						}
					}

					return issues;
				},
				"Validate flow receiver network for common issues",
				py::arg("receivers"),
				py::arg("connector"));
}

// ==============================================
// MAIN BINDING FUNCTION
// ==============================================

inline void
bind_cordonnier_flow_rerouter_module(py::module& m)
{
	// Bind enums first
	bind_cordonnier_enums<float>(m, "F32");
	bind_cordonnier_enums<double>(m, "F64");

	// Bind template classes for different types
	bind_all_cordonnier_types<float>(m, "F32");
	bind_all_cordonnier_types<double>(m, "F64");

	// Bind utility functions
	bind_cordonnier_utilities<float>(m, "F32");
	bind_cordonnier_utilities<double>(m, "F64");

	// Module documentation
	m.doc() = R"pbdoc(
        🌊 DAGGER2 CORDONNIER FLOW REROUTER MODULE 🌊
        ============================================

        Advanced flow routing for topographies with depressions using the
        Cordonnier et al. 2019 linear complexity algorithm.

        This module provides explicit flow path enforcement through basin graph
        construction and minimum spanning tree optimization, offering superior
        performance compared to traditional priority flood approaches.

        Key Features:
        - Linear O(n) time complexity for planar graphs
        - Multiple enforcement strategies (filling, carving, correction)
        - Comprehensive boundary condition support
        - Detailed statistics and validation
        - Basin graph reusability for other algorithms

        Reference:
        Cordonnier, G., Bovy, B., & Braun, J. (2019). A versatile, linear
        complexity algorithm for flow routing in topographies with depressions.
        Earth Surface Dynamics, 7(2), 549-562.
    )pbdoc";
}

} // namespace dagger2
