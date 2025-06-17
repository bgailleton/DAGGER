#pragma once

#include "dg2_graph.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace dagger2 {

// ==============================================
// ENUM BINDINGS
// ==============================================

inline void
bind_graph_enums(py::module& m)
{
	// FlowMode enumeration
	py::enum_<FlowMode>(m, "FlowMode")
		.value(
			"SINGLE_FLOW", FlowMode::SINGLE_FLOW, "Single receiver per donor (D4/D8)")
		.value("MULTIPLE_FLOW",
					 FlowMode::MULTIPLE_FLOW,
					 "Multiple receivers per donor (MFD)")
		.export_values();

	// SingleFlowMethod enumeration
	py::enum_<SingleFlowMethod>(m, "SingleFlowMethod")
		.value("STEEPEST_DESCENT",
					 SingleFlowMethod::STEEPEST_DESCENT,
					 "Always choose steepest downslope neighbor")
		.value("RANDOM_UNIFORM",
					 SingleFlowMethod::RANDOM_UNIFORM,
					 "Random selection among downslope neighbors")
		.value("RANDOM_WEIGHTED",
					 SingleFlowMethod::RANDOM_WEIGHTED,
					 "Gradient-weighted random selection")
		.export_values();

	// MultipleFlowMethod enumeration
	py::enum_<MultipleFlowMethod>(m, "MultipleFlowMethod")
		.value("FREEMAN_1991",
					 MultipleFlowMethod::FREEMAN_1991,
					 "Freeman (1991) - proportional to slope")
		.value("QUINN_1991",
					 MultipleFlowMethod::QUINN_1991,
					 "Quinn et al. (1991) - proportional to slope^alpha")
		.value("SEIBERT_MCGLYNN_2007",
					 MultipleFlowMethod::SEIBERT_MCGLYNN_2007,
					 "Seibert & McGlynn (2007) - adaptive")
		.value("TARBOTON_1997",
					 MultipleFlowMethod::TARBOTON_1997,
					 "Tarboton (1997) - infinite direction (D∞)")
		.value("HOLMGREN_1994",
					 MultipleFlowMethod::HOLMGREN_1994,
					 "Holmgren (1994) - with exponent parameter")
		.export_values();
}

// Note: Vector2 is already bound in the Perlin noise module, so we skip it here
// to avoid "type already registered" conflicts

// ==============================================
// FLOW LINK BINDINGS
// ==============================================

template<typename T>
void
bind_flow_link(py::module& m, const std::string& suffix)
{
	using FlowLinkType = FlowLink<T>;

	py::class_<FlowLinkType>(m, ("FlowLink" + suffix).c_str())
		.def(py::init<>())
		.def(py::init<size_t, size_t, Direction, T, T, T, bool>(),
				 py::arg("from"),
				 py::arg("to"),
				 py::arg("direction"),
				 py::arg("weight"),
				 py::arg("gradient"),
				 py::arg("distance"),
				 py::arg("boundary_exit") = false)
		.def_readwrite("from_index", &FlowLinkType::from_index)
		.def_readwrite("to_index", &FlowLinkType::to_index)
		.def_readwrite("direction", &FlowLinkType::direction)
		.def_readwrite("weight", &FlowLinkType::weight)
		.def_readwrite("gradient", &FlowLinkType::gradient)
		.def_readwrite("distance", &FlowLinkType::distance)
		.def_readwrite("is_boundary_exit", &FlowLinkType::is_boundary_exit)
		.def("is_valid", &FlowLinkType::is_valid);
}

// ==============================================
// FLOW NODE BINDINGS
// ==============================================

template<typename T>
void
bind_flow_node(py::module& m, const std::string& suffix)
{
	using FlowNodeType = FlowNode<T>;

	py::class_<FlowNodeType>(m, ("FlowNode" + suffix).c_str())
		.def(py::init<>())
		.def_readwrite("index", &FlowNodeType::index)
		.def_readwrite("row", &FlowNodeType::row)
		.def_readwrite("col", &FlowNodeType::col)
		.def_readwrite("elevation", &FlowNodeType::elevation)
		.def_readwrite("boundary_type", &FlowNodeType::boundary_type)
		.def_readwrite("topological_order", &FlowNodeType::topological_order)
		.def_readwrite("num_donors", &FlowNodeType::num_donors)
		.def_readwrite("num_receivers", &FlowNodeType::num_receivers)
		.def_readwrite("donor_links", &FlowNodeType::donor_links)
		.def_readwrite("receiver_links", &FlowNodeType::receiver_links)
		.def_readwrite("is_source", &FlowNodeType::is_source)
		.def_readwrite("is_sink", &FlowNodeType::is_sink)
		.def_readwrite("is_outlet", &FlowNodeType::is_outlet)
		.def_readwrite("is_pit", &FlowNodeType::is_pit);
}

// ==============================================
// FLOW GRAPH BINDINGS
// ==============================================

// ==============================================
// FLOW GRAPH BINDINGS (Complete with Accumulation Functions) - FIXED
// ==============================================

template<typename T>
void
bind_flow_graph(py::module& m, const std::string& suffix)
{
	using FlowGraphType = FlowGraph<T>;
	using ConnectorType = Connector<T>;
	using ArrayRefType = ArrayRef<T>;
	using ArrayRefSizeType = ArrayRef<size_t>;

	py::class_<FlowGraphType>(m, ("FlowGraph" + suffix).c_str())
		// Fixed constructor - uses std::shared_ptr<ArrayRef<T>>
		.def(
			py::init<std::shared_ptr<ConnectorType>, std::shared_ptr<ArrayRefType>>(),
			py::arg("connector"),
			py::arg("elevation"))
		.def(py::init<std::shared_ptr<ConnectorType>,
									std::shared_ptr<ArrayRefType>,
									FlowMode>(),
				 py::arg("connector"),
				 py::arg("elevation"),
				 py::arg("flow_mode"))
		.def(py::init<std::shared_ptr<ConnectorType>,
									std::shared_ptr<ArrayRefType>,
									FlowMode,
									SingleFlowMethod>(),
				 py::arg("connector"),
				 py::arg("elevation"),
				 py::arg("flow_mode"),
				 py::arg("single_method"))
		.def(py::init<std::shared_ptr<ConnectorType>,
									std::shared_ptr<ArrayRefType>,
									FlowMode,
									SingleFlowMethod,
									MultipleFlowMethod>(),
				 py::arg("connector"),
				 py::arg("elevation"),
				 py::arg("flow_mode"),
				 py::arg("single_method"),
				 py::arg("multiple_method"))

		// Basic properties
		.def_property_readonly("size", &FlowGraphType::size)
		.def_property_readonly("num_nodes", &FlowGraphType::num_nodes)
		.def_property_readonly("num_links", &FlowGraphType::num_links)
		.def_property_readonly("num_sources", &FlowGraphType::num_sources)
		.def_property_readonly("num_sinks", &FlowGraphType::num_sinks)
		.def_property_readonly("num_outlets", &FlowGraphType::num_outlets)
		.def_property_readonly("is_built", &FlowGraphType::is_built)

		// Configuration getters
		.def("get_flow_mode", &FlowGraphType::get_flow_mode)
		.def("get_single_method", &FlowGraphType::get_single_method)
		.def("get_multiple_method", &FlowGraphType::get_multiple_method)
		.def("get_flow_exponent", &FlowGraphType::get_flow_exponent)
		.def("get_min_gradient", &FlowGraphType::get_min_gradient)
		.def("get_random_seed", &FlowGraphType::get_random_seed)

		// Configuration setters
		.def("set_flow_mode", &FlowGraphType::set_flow_mode)
		.def("set_single_method", &FlowGraphType::set_single_method)
		.def("set_multiple_method", &FlowGraphType::set_multiple_method)
		.def("set_flow_exponent", &FlowGraphType::set_flow_exponent)
		.def("set_min_gradient", &FlowGraphType::set_min_gradient)
		.def("set_random_seed", &FlowGraphType::set_random_seed)

		// Graph building
		.def("build", &FlowGraphType::build)

		// Node and link access
		.def("get_node",
				 py::overload_cast<size_t>(&FlowGraphType::get_node, py::const_),
				 py::return_value_policy::reference_internal)
		.def(
			"get_node",
			py::overload_cast<size_t, size_t>(&FlowGraphType::get_node, py::const_),
			py::return_value_policy::reference_internal)
		.def("get_link",
				 &FlowGraphType::get_link,
				 py::return_value_policy::reference_internal)

		// Connectivity queries
		.def("get_donors",
				 py::overload_cast<size_t>(&FlowGraphType::get_donors, py::const_))
		.def("get_receivers",
				 py::overload_cast<size_t>(&FlowGraphType::get_receivers, py::const_))
		.def("get_sources", &FlowGraphType::get_sources)
		.def("get_sinks", &FlowGraphType::get_sinks)
		.def("get_outlets", &FlowGraphType::get_outlets)

		// ==============================================
		// ACCUMULATION FUNCTIONS
		// ==============================================

		// Generic accumulation with custom function
		.def(
			"accumulate",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& input,
				 py::function accumulate_func,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					input,
					[accumulate_func](T upstream, T local, T weight) -> T {
						return accumulate_func(upstream, local, weight).template cast<T>();
					},
					include_self);
			},
			py::arg("result"),
			py::arg("input"),
			py::arg("accumulate_func"),
			py::arg("include_self") = true,
			"Generic accumulation with custom function")

		// // Drainage area accumulation
		// .def("accumulate_drainage_area",
		// &FlowGraphType::accumulate_drainage_area,
		//      py::arg("drainage_area"), py::arg("cell_area"),
		//      "Accumulate drainage area using cell areas")

		// Drainage area accumulation
		.def("compute_drainage_area",
				 &FlowGraphType::template compute_drainage_area<T>,
				 py::arg("cell_area"),
				 "Accumulate drainage area using cell areas and returns a Grid2D")

		// Flow accumulation (cell count)
		.def(
			"accumulate_flow",
			[](FlowGraphType& self, ArrayRefSizeType& flow_accumulation) {
				self.accumulate_flow(flow_accumulation);
			},
			py::arg("flow_accumulation"),
			"Accumulate flow (count of upstream cells)")

		// Specialized accumulation variants with lambda wrappers - Fixed template
		// calls
		.def(
			"accumulate_sum",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& input,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					input,
					[](T upstream, T local, T weight) -> T {
						return upstream; // Simple sum
					},
					include_self);
			},
			py::arg("result"),
			py::arg("input"),
			py::arg("include_self") = true,
			"Simple sum accumulation")

		.def(
			"accumulate_weighted_sum",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& input,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					input,
					[](T upstream, T local, T weight) -> T {
						return upstream * weight; // Weighted by flow proportion
					},
					include_self);
			},
			py::arg("result"),
			py::arg("input"),
			py::arg("include_self") = true,
			"Weighted sum accumulation using flow weights")

		.def(
			"accumulate_with_loss",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& input,
				 T loss_rate,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					input,
					[loss_rate](T upstream, T local, T weight) -> T {
						// Apply loss rate during transport
						return upstream * (T(1.0) - loss_rate);
					},
					include_self);
			},
			py::arg("result"),
			py::arg("input"),
			py::arg("loss_rate"),
			py::arg("include_self") = true,
			"Accumulation with loss during transport")

		.def(
			"accumulate_discharge",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& precipitation,
				 T runoff_coefficient,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					precipitation,
					[runoff_coefficient](T upstream, T local, T weight) -> T {
						return upstream +
									 (local * runoff_coefficient); // Add local contribution
					},
					include_self);
			},
			py::arg("result"),
			py::arg("precipitation"),
			py::arg("runoff_coefficient"),
			py::arg("include_self") = true,
			"Discharge accumulation with runoff coefficient")

		.def(
			"accumulate_sediment",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& erosion_rate,
				 T transport_capacity,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					erosion_rate,
					[transport_capacity](T upstream, T local, T weight) -> T {
						// Simplified sediment transport with capacity limitation
						T total = upstream + local;
						return std::min(total, transport_capacity);
					},
					include_self);
			},
			py::arg("result"),
			py::arg("erosion_rate"),
			py::arg("transport_capacity"),
			py::arg("include_self") = true,
			"Sediment accumulation with transport capacity")

		.def(
			"accumulate_max",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& input,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					input,
					[](T upstream, T local, T weight) -> T {
						return std::max(upstream, local); // Maximum value
					},
					include_self);
			},
			py::arg("result"),
			py::arg("input"),
			py::arg("include_self") = true,
			"Maximum value accumulation")

		.def(
			"accumulate_min",
			[](FlowGraphType& self,
				 ArrayRefType& result,
				 const ArrayRefType& input,
				 bool include_self) {
				self.template accumulate<T>(
					result,
					input,
					[](T upstream, T local, T weight) -> T {
						return std::min(upstream, local); // Minimum value
					},
					include_self);
			},
			py::arg("result"),
			py::arg("input"),
			py::arg("include_self") = true,
			"Minimum value accumulation")

		// ==============================================
		// GRAPH ANALYSIS
		// ==============================================

		.def("get_watershed",
				 &FlowGraphType::get_watershed,
				 py::arg("outlet_index"),
				 "Get all nodes in watershed contributing to outlet")

		.def("find_strongly_connected_components",
				 &FlowGraphType::find_strongly_connected_components,
				 "Find strongly connected components for cycle detection")

		// Statistics
		.def("compute_statistics", &FlowGraphType::compute_statistics)
		.def("get_connector", &FlowGraphType::get_connector);
}

// ==============================================
// GRAPH ALGORITHMS BINDINGS (Basic only)
// ==============================================

template<typename T>
void
bind_graph_algorithms(py::module& m, const std::string& suffix)
{
	using AlgorithmsType = GraphAlgorithms<T>;
	using EdgeType = typename AlgorithmsType::Edge;

	// Edge structure
	py::class_<EdgeType>(m, ("Edge" + suffix).c_str())
		.def(py::init<size_t, size_t, T, Direction>(),
				 py::arg("from"),
				 py::arg("to"),
				 py::arg("weight"),
				 py::arg("direction"))
		.def_readwrite("from", &EdgeType::from)
		.def_readwrite("to", &EdgeType::to)
		.def_readwrite("weight", &EdgeType::weight)
		.def_readwrite("direction", &EdgeType::direction);

	// Main algorithms class (only bind methods that actually exist)
	py::class_<AlgorithmsType>(m, ("GraphAlgorithms" + suffix).c_str())
		.def(py::init<std::shared_ptr<FlowGraph<T>>>())

		// Basic algorithms that exist
		.def("kruskal_mst", &AlgorithmsType::kruskal_mst)
		.def("prim_mst", &AlgorithmsType::prim_mst)
		.def("dijkstra", &AlgorithmsType::dijkstra, py::arg("source"))
		.def("betweenness_centrality", &AlgorithmsType::betweenness_centrality)
		.def("pagerank",
				 &AlgorithmsType::pagerank,
				 py::arg("damping_factor") = 0.85,
				 py::arg("tolerance") = 1e-8,
				 py::arg("max_iterations") = 100)
		.def("louvain_communities", &AlgorithmsType::louvain_communities);
}

// ==============================================
// QUANTUM GRAPH OPERATIONS BINDINGS (Basic only)
// ==============================================

template<typename T>
void
bind_quantum_graph_ops(py::module& m, const std::string& suffix)
{
	using QuantumOpsType = QuantumGraphOps<T>;

	py::class_<QuantumOpsType>(m, ("QuantumGraphOps" + suffix).c_str())
		.def(py::init<std::shared_ptr<FlowGraph<T>>>())

		// Only bind methods that actually exist
		.def("quantum_walk",
				 &QuantumOpsType::quantum_walk,
				 py::arg("start_node"),
				 py::arg("steps"),
				 py::arg("dt") = 0.1);
}

// ==============================================
// GRAPH NEURAL NETWORKS BINDINGS (Basic only)
// ==============================================

template<typename T>
void
bind_graph_neural_net(py::module& m, const std::string& suffix)
{
	using NeuralNetType = GraphNeuralNet<T>;

	py::class_<NeuralNetType>(m, ("GraphNeuralNet" + suffix).c_str())
		.def(py::init<std::shared_ptr<FlowGraph<T>>>())

		// Only bind methods that actually exist
		.def("gcn_forward", &NeuralNetType::gcn_forward);
}

// ==============================================
// OPTIMIZATION ALGORITHMS BINDINGS (Basic only)
// ==============================================

template<typename T>
void
bind_optimization_algorithms(py::module& m, const std::string& suffix)
{
	using AntColonyType = AntColonyGraphOpt<T>;
	using EvolutionaryType = EvolutionaryGraphOpt<T>;

	// Ant Colony Optimization (only bind methods that exist)
	py::class_<AntColonyType>(m, ("AntColonyGraphOpt" + suffix).c_str())
		.def(py::init<std::shared_ptr<FlowGraph<T>>>())
		.def("solve_tsp", &AntColonyType::solve_tsp);

	// Evolutionary Optimization (only bind methods that exist)
	py::class_<EvolutionaryType>(m, ("EvolutionaryGraphOpt" + suffix).c_str())
		.def(py::init<std::shared_ptr<FlowGraph<T>>>());
}

// ==============================================
// ULTIMATE GRAPH ENGINE BINDINGS
// ==============================================

template<typename T>
void
bind_ultimate_graph_engine(py::module& m, const std::string& suffix)
{
	using EngineType = UltimateGraphEngine<T>;
	using AnalysisType = typename EngineType::ComprehensiveAnalysis;

	// Comprehensive analysis structure
	py::class_<AnalysisType>(m, ("ComprehensiveAnalysis" + suffix).c_str())
		.def_readonly("basic_stats", &AnalysisType::basic_stats)
		.def_readonly("centrality_measures", &AnalysisType::centrality_measures)
		.def_readonly("communities", &AnalysisType::communities)
		.def_readonly("quantum_amplitudes", &AnalysisType::quantum_amplitudes)
		.def_readonly("degree_entropy", &AnalysisType::degree_entropy)
		.def_readonly("box_counting_dimension",
									&AnalysisType::box_counting_dimension)
		.def_readonly("summary_report", &AnalysisType::summary_report);

	// Ultimate Graph Engine
	py::class_<EngineType>(m, ("UltimateGraphEngine" + suffix).c_str())
		.def(py::init<std::shared_ptr<FlowGraph<T>>>())

		// Access to sub-modules
		.def("get_flow_graph", &EngineType::get_flow_graph)
		.def("get_algorithms",
				 &EngineType::get_algorithms,
				 py::return_value_policy::reference_internal)
		.def("get_quantum_ops",
				 &EngineType::get_quantum_ops,
				 py::return_value_policy::reference_internal)
		.def("get_neural_net",
				 &EngineType::get_neural_net,
				 py::return_value_policy::reference_internal)
		.def("get_evolutionary_opt",
				 &EngineType::get_evolutionary_opt,
				 py::return_value_policy::reference_internal)
		.def("get_ant_colony",
				 &EngineType::get_ant_colony,
				 py::return_value_policy::reference_internal)

		// Comprehensive analysis
		.def("analyze_everything", &EngineType::analyze_everything)
		.def("benchmark_all_algorithms", &EngineType::benchmark_all_algorithms);
}

// ==============================================
// CONVENIENCE FACTORY FUNCTIONS
// ==============================================

template<typename T>
void
bind_graph_factory_functions(py::module& m, const std::string& suffix)
{
	// Factory functions
	m.def(("create_ultimate_graph_engine" + suffix).c_str(),
				&create_ultimate_graph_engine<T>,
				"Create ultimate graph engine with basic parameters",
				py::arg("rows"),
				py::arg("cols"),
				py::arg("elevation"),
				py::arg("flow_mode") = FlowMode::SINGLE_FLOW,
				py::arg("connectivity") = ConnectivityType::D8);

	m.def(("create_landscape_graph_engine" + suffix).c_str(),
				&create_landscape_graph_engine<T>,
				"Create graph engine optimized for landscape modeling",
				py::arg("elevation_data"),
				py::arg("rows"),
				py::arg("cols"));
}

// ==============================================
// COMPREHENSIVE BINDING FUNCTION
// ==============================================

template<typename T>
void
bind_all_graph_types(py::module& m, const std::string& suffix)
{
	// Basic types (skip Vector2 as it's already bound in Perlin noise module)
	bind_flow_link<T>(m, suffix);
	bind_flow_node<T>(m, suffix);

	// Core graph classes
	bind_flow_graph<T>(m, suffix);
	bind_graph_algorithms<T>(m, suffix);

	// Advanced modules (only what exists)
	bind_quantum_graph_ops<T>(m, suffix);
	bind_graph_neural_net<T>(m, suffix);
	bind_optimization_algorithms<T>(m, suffix);

	// Ultimate engine
	bind_ultimate_graph_engine<T>(m, suffix);

	// Factory functions
	bind_graph_factory_functions<T>(m, suffix);
}

// ==============================================
// MAIN BINDING FUNCTION
// ==============================================

inline void
bind_graph_module(py::module& m)
{
	// Bind enums first
	bind_graph_enums(m);

	// Bind for different numeric types
	bind_all_graph_types<float>(m, "F32");
	bind_all_graph_types<double>(m, "F64");

	// Module documentation
	m.doc() = R"pbdoc(
        🌟 DAGGER2 GRAPH ANALYSIS ENGINE 🌟
        ===================================

        Advanced graph analysis capabilities for flow routing and network analysis!

        📚 FEATURES:
        • Flow graph construction with single/multiple flow directions
        • Comprehensive node and link management
        • Topological ordering and level sets
        • Graph algorithms (shortest paths, MST, centrality)
        • Basic quantum-inspired operations
        • Machine learning on graphs (GCN)
        • Optimization algorithms (ant colony, evolutionary)

        🎯 APPLICATIONS:
        • Landscape evolution modeling
        • Hydrology and drainage analysis
        • Flow accumulation and routing
        • Watershed delineation
        • Network analysis

        🚀 USAGE:
        ```python
        import dagger2

        # Create flow graph
        connector = dagger2.ConnectorF64(rows, cols)
        elevation = dagger2.ArrayRefF64(elevation_data)
        graph = dagger2.FlowGraphF64(connector, elevation, dagger2.FlowMode.SINGLE_FLOW)
        graph.build()

        # Run algorithms
        algorithms = dagger2.GraphAlgorithmsF64(graph)
        centrality = algorithms.betweenness_centrality()
        communities = algorithms.louvain_communities()

        # Create ultimate engine
        engine = dagger2.create_ultimate_graph_engineF64(rows, cols, elevation)
        analysis = engine.analyze_everything()
        ```

        ✨ Perfect for advanced landscape modeling and network analysis! ✨
    )pbdoc";
}

} // namespace dagger2
