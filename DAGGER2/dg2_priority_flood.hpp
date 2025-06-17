#pragma once

#include "dg2_BCs.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include "dg2_fastconnector.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace dagger2 {

// ===================================================================
// CONNECTOR ADAPTER TRAITS - Unified interface for both connectors
// ===================================================================

/**
 * Connector traits to abstract differences between Connector and FastConnector
 */
template<typename ConnectorType>
struct ConnectorTraits;

// Specialization for regular Connector<T>
template<typename T>
struct ConnectorTraits<Connector<T>>
{
	using value_type = T;
	using neighbor_type = Neighbor; // Global Neighbor struct

	static constexpr bool is_fast_connector = false;

	static inline size_t to_1d(const Connector<T>& conn, size_t row, size_t col)
	{
		return conn.to_1d(row, col);
	}

	static inline std::pair<size_t, size_t> to_2d(const Connector<T>& conn,
																								size_t idx)
	{
		return conn.to_2d(idx);
	}

	static inline NodeType get_boundary_type(const Connector<T>& conn,
																					 size_t row,
																					 size_t col)
	{
		return conn.get_boundary_type(row, col);
	}

	static inline auto get_effective_valid_neighbors(const Connector<T>& conn,
																									 size_t row,
																									 size_t col)
	{
		return conn.get_effective_valid_neighbors(row, col);
	}
};

/**
 * Specialization for FastConnector types
 */
template<typename T, ConnectivityType CONN>
struct ConnectorTraits<FastConnector<T, CONN>>
{
	using value_type = T;
	using neighbor_type = typename FastConnector<T, CONN>::Neighbor;

	static constexpr bool is_fast_connector = true;

	static inline size_t to_1d(const FastConnector<T, CONN>& conn,
														 size_t row,
														 size_t col)
	{
		return conn.to_1d(row, col);
	}

	static inline std::pair<size_t, size_t> to_2d(
		const FastConnector<T, CONN>& conn,
		size_t idx)
	{
		return conn.to_2d(idx);
	}

	static inline NodeType get_boundary_type(const FastConnector<T, CONN>& conn,
																					 size_t row,
																					 size_t col)
	{
		return conn.get_boundary_type(row, col);
	}

	static inline auto get_effective_valid_neighbors(
		const FastConnector<T, CONN>& conn,
		size_t row,
		size_t col)
	{
		return conn.get_effective_valid_neighbors(row, col);
	}
};

/**
 * Priority Flood Plus Epsilon: Advanced depression filling for topographic data
 *
 * This templated implementation works with both Connector<T> and
 * FastConnector<T, CONN> types through the ConnectorTraits abstraction,
 * maintaining full backward compatibility.
 */
template<typename T>
class PriorityFlood
{
public:
	using PriorityQueue = std::priority_queue<std::pair<T, size_t>,
																						std::vector<std::pair<T, size_t>>,
																						std::greater<>>;

	/**
	 * Configuration structure for the algorithm
	 */
	struct Config
	{
		T epsilon = static_cast<T>(1e-6);
		T max_fill_depth = std::numeric_limits<T>::max();
		size_t max_iterations = SIZE_MAX;
		bool preserve_flat_areas = true;
		bool validate_result = false;
		bool generate_statistics = false;
		bool verbose = false;
	};

	/**
	 * Results structure containing algorithm outputs and statistics
	 */
	struct Results
	{
		size_t iterations = 0;
		size_t nodes_processed = 0;
		size_t nodes_filled = 0;
		T total_fill_volume = 0;
		T max_fill_depth = 0;
		bool convergence_reached = false;
		std::vector<std::string> warnings;

		// Additional statistics
		size_t boundary_seeds = 0;
		double processing_time = 0.0;
		size_t memory_usage = 0;
	};

	/**
	 * Main depression filling method - works with any connector type
	 */
	template<typename ConnectorType>
	static Results fill_depressions(Grid2D<T>& elevation,
																	const ConnectorType& connector,
																	const Config& config = Config())
	{
		using Traits = ConnectorTraits<ConnectorType>;

		Results results;
		auto start_time = std::chrono::high_resolution_clock::now();

		// Validate inputs
		if (elevation.rows() != connector.rows() ||
				elevation.cols() != connector.cols()) {
			throw std::invalid_argument(
				"Elevation grid and connector dimensions mismatch");
		}

		// Initialize working configuration
		Config working_config = config;
		if (working_config.epsilon <= 0) {
			working_config.epsilon = static_cast<T>(1e-6);
		}

		// Create backup if statistics are requested
		std::unique_ptr<Grid2D<T>> original_elevation;
		if (working_config.generate_statistics) {
			original_elevation = backup_elevation(elevation);
		}

		// Initialize data structures
		PriorityQueue pq;
		std::vector<bool> processed(elevation.size(), false);

		// Phase 1: Initialize boundary seeds (templated for both connector types)
		initialize_boundary_seeds_templated(
			elevation, connector, pq, processed, results);

		if (pq.empty()) {
			results.warnings.push_back(
				"No valid boundary seeds found - skipping flood fill");
			return results;
		}

		// Phase 2: Priority flood algorithm (templated for both connector types)
		execute_priority_flood_templated(
			elevation, connector, pq, processed, working_config, results);

		// Phase 3: Post-processing and validation
		if (working_config.validate_result) {
			validate_result_templated(elevation, connector, results);
		}

		if (working_config.generate_statistics && original_elevation) {
			compute_statistics(elevation, *original_elevation, results);
		}

		// Record timing
		auto end_time = std::chrono::high_resolution_clock::now();
		results.processing_time =
			std::chrono::duration<double>(end_time - start_time).count();

		return results;
	}

	/**
	 * Simplified interface with default configuration - backward compatible
	 */
	template<typename ConnectorType>
	static Results fill_depressions(Grid2D<T>& elevation,
																	const ConnectorType& connector)
	{
		return fill_depressions(elevation, connector, Config());
	}

	/**
	 * Fill depressions with custom epsilon value - backward compatible
	 */
	template<typename ConnectorType>
	static Results fill_depressions(Grid2D<T>& elevation,
																	const ConnectorType& connector,
																	T epsilon)
	{
		Config config;
		config.epsilon = epsilon;
		return fill_depressions(elevation, connector, config);
	}

	/**
	 * Create a backup copy of elevation before filling
	 */
	static std::unique_ptr<Grid2D<T>> backup_elevation(const Grid2D<T>& elevation)
	{
		std::vector<T> backup_data(elevation.begin(), elevation.end());
		return std::make_unique<Grid2D<T>>(
			backup_data, elevation.rows(), elevation.cols());
	}

	/**
	 * Restore elevation from backup
	 */
	static void restore_elevation(Grid2D<T>& elevation, const Grid2D<T>& backup)
	{
		if (elevation.size() != backup.size()) {
			throw std::invalid_argument("Elevation and backup size mismatch");
		}
		std::copy(backup.begin(), backup.end(), elevation.begin());
	}

private:
	/**
	 * Templated boundary seed initialization for both connector types
	 */
	template<typename ConnectorType>
	static void initialize_boundary_seeds_templated(
		const Grid2D<T>& elevation,
		const ConnectorType& connector,
		PriorityQueue& pq,
		std::vector<bool>& processed,
		Results& results)
	{
		using Traits = ConnectorTraits<ConnectorType>;

		const size_t rows = elevation.rows();
		const size_t cols = elevation.cols();

		// Pre-allocate boundary seeds vector
		std::vector<size_t> boundary_candidates;
		boundary_candidates.reserve(2 * (rows + cols));

		// Only check actual perimeter nodes
		for (size_t r = 0; r < rows; ++r) {
			if (r == 0 || r == rows - 1) {
				// Top and bottom rows - check all columns
				for (size_t c = 0; c < cols; ++c) {
					boundary_candidates.push_back(Traits::to_1d(connector, r, c));
				}
			} else {
				// Middle rows - only check left and right edges
				boundary_candidates.push_back(Traits::to_1d(connector, r, 0));
				boundary_candidates.push_back(Traits::to_1d(connector, r, cols - 1));
			}
		}

		// Process boundary candidates
		for (size_t idx : boundary_candidates) {
			auto [r, c] = Traits::to_2d(connector, idx);
			NodeType node_type = Traits::get_boundary_type(connector, r, c);

			// Skip NO_DATA nodes
			if (node_type == NodeType::NO_DATA) {
				processed[idx] = true;
				continue;
			}

			// Only add nodes that allow outflow
			if (NodeTypeUtils::allows_outflow(node_type)) {
				T elev = elevation(r, c);
				pq.emplace(elev, idx);
				processed[idx] = true;
				results.boundary_seeds++;
			}
		}
	}

	/**
	 * Templated priority flood execution for both connector types
	 */
	template<typename ConnectorType>
	static void execute_priority_flood_templated(Grid2D<T>& elevation,
																							 const ConnectorType& connector,
																							 PriorityQueue& pq,
																							 std::vector<bool>& processed,
																							 const Config& config,
																							 Results& results)
	{
		using Traits = ConnectorTraits<ConnectorType>;

		results.iterations = 0;

		while (!pq.empty() && results.iterations < config.max_iterations) {
			auto [current_elev, current_idx] = pq.top();
			pq.pop();

			results.iterations++;
			results.nodes_processed++;

			// Update elevation if needed
			T original_elev = elevation[current_idx];
			if (current_elev > original_elev) {
				T fill_depth = current_elev - original_elev;

				// Check maximum fill depth constraint
				if (fill_depth > config.max_fill_depth) {
					results.warnings.push_back("Maximum fill depth exceeded at node " +
																		 std::to_string(current_idx));
					continue;
				}

				elevation[current_idx] = current_elev;
				results.nodes_filled++;
				results.total_fill_volume += fill_depth;
				results.max_fill_depth = std::max(results.max_fill_depth, fill_depth);
			}

			// Process neighbors using connector abstraction
			auto [row, col] = Traits::to_2d(connector, current_idx);
			auto neighbors =
				Traits::get_effective_valid_neighbors(connector, row, col);

			for (const auto& neighbor : neighbors) {
				if (processed[neighbor.index])
					continue;

				// Skip NO_DATA neighbors
				if (neighbor.boundary_type == NodeType::NO_DATA)
					continue;

				T neighbor_elev = elevation[neighbor.index];
				T required_fill_elevation = current_elev;

				// Apply epsilon increment if needed
				if (neighbor_elev <= current_elev) {
					if (config.preserve_flat_areas &&
							std::abs(neighbor_elev - current_elev) < config.epsilon) {
						required_fill_elevation = current_elev;
					} else {
						required_fill_elevation = current_elev + config.epsilon;
					}
				} else {
					required_fill_elevation = neighbor_elev;
				}

				// Add to priority queue
				pq.emplace(required_fill_elevation, neighbor.index);
				processed[neighbor.index] = true;
			}
		}

		results.convergence_reached = pq.empty();

		if (results.iterations >= config.max_iterations) {
			results.warnings.push_back(
				"Maximum iterations reached - algorithm may not have converged");
		}
	}

	/**
	 * Templated result validation for both connector types
	 */
	template<typename ConnectorType>
	static void validate_result_templated(const Grid2D<T>& elevation,
																				const ConnectorType& connector,
																				Results& results)
	{
		using Traits = ConnectorTraits<ConnectorType>;

		const size_t rows = elevation.rows();
		const size_t cols = elevation.cols();
		size_t violations = 0;

		for (size_t r = 0; r < rows; ++r) {
			for (size_t c = 0; c < cols; ++c) {
				size_t current_idx = Traits::to_1d(connector, r, c);

				// Skip NO_DATA nodes
				if (Traits::get_boundary_type(connector, r, c) == NodeType::NO_DATA) {
					continue;
				}

				T current_elev = elevation[current_idx];
				auto neighbors = Traits::get_effective_valid_neighbors(connector, r, c);

				for (const auto& neighbor : neighbors) {
					if (neighbor.boundary_type == NodeType::NO_DATA)
						continue;

					T neighbor_elev = elevation[neighbor.index];
					if (current_elev > neighbor_elev + static_cast<T>(1e-10)) {
						violations++;
					}
				}
			}
		}

		if (violations > 0) {
			results.warnings.push_back("Validation found " +
																 std::to_string(violations) +
																 " drainage violations");
		}
	}

	/**
	 * Compute algorithm statistics
	 */
	static void compute_statistics(const Grid2D<T>& elevation,
																 const Grid2D<T>& original_elevation,
																 Results& results)
	{
		T total_change = 0;
		T max_change = 0;
		size_t changed_nodes = 0;

		for (size_t i = 0; i < elevation.size(); ++i) {
			T change = elevation[i] - original_elevation[i];
			if (change > static_cast<T>(1e-10)) {
				total_change += change;
				max_change = std::max(max_change, change);
				changed_nodes++;
			}
		}

		results.total_fill_volume = total_change;
		results.max_fill_depth = max_change;
		results.nodes_filled = changed_nodes;
	}

public:
	/**
	 * Generate a detailed report of the results
	 */
	static std::string generate_report(const Results& results)
	{
		std::ostringstream report;
		report << "=== Priority Flood Algorithm Results ===\n";
		report << "Iterations: " << results.iterations << "\n";
		report << "Nodes processed: " << results.nodes_processed << "\n";
		report << "Nodes filled: " << results.nodes_filled << "\n";
		report << "Boundary seeds: " << results.boundary_seeds << "\n";
		report << "Total fill volume: " << results.total_fill_volume << "\n";
		report << "Maximum fill depth: " << results.max_fill_depth << "\n";
		report << "Processing time: " << results.processing_time << " seconds\n";
		report << "Convergence reached: "
					 << (results.convergence_reached ? "Yes" : "No") << "\n";

		if (!results.warnings.empty()) {
			report << "\nWarnings:\n";
			for (const auto& warning : results.warnings) {
				report << "  - " << warning << "\n";
			}
		}

		return report.str();
	}
};

// Convenience type aliases
using PriorityFloodF32 = PriorityFlood<float>;
using PriorityFloodF64 = PriorityFlood<double>;

} // namespace dagger2

/*
The version bellow is more optimsed but does not manages the boundary conditions
periodic and all

*/

// #pragma once

// #include "dg2_BCs.hpp"
// #include "dg2_array.hpp"
// #include "dg2_connector.hpp"
// #include <algorithm>
// #include <cmath>
// #include <limits>
// #include <memory>
// #include <queue>
// #include <unordered_set>
// #include <vector>

// namespace dagger2 {

// /**
//  * Priority Flood Plus Epsilon: Advanced depression filling for topographic
//  data
//  *
//  * This implementation provides comprehensive depression filling using the
//  * Priority Flood algorithm with epsilon increment, properly handling all
//  * boundary conditions and no-data values. The algorithm ensures proper
//  drainage
//  * by filling local minima while preserving the overall topographic
//  structure.
//  *
//  * Key Features:
//  * - Handles all NodeType boundary conditions through Connector
//  * - In-place elevation modification with optional backup
//  * - Configurable epsilon increment for depression depth control
//  * - Robust handling of NO_DATA regions
//  * - Performance optimized with priority queue and visited tracking
//  * - Detailed statistics and validation
//  */
// template<typename T = double>
// class PriorityFlood
// {
// public:
//     /**
//      * Configuration structure for priority flood parameters
//      */
//     struct Config
//     {
//         T epsilon; // Minimum elevation increment for filling (default: 1e-6)
//         bool preserve_flat_areas; // Keep flat areas flat vs add tiny
//         gradient
//                                                             // (default:
//                                                             false)
//         bool validate_result;           // Run post-processing validation
//         (default: true) bool generate_statistics; // Collect detailed
//         statistics (default: false) T max_fill_depth;                   //
//         Maximum allowed fill depth (default: infinity) size_t
//             max_iterations; // Safety limit on iterations (default: grid_size
//             * 2)

//         Config()
//             : epsilon(static_cast<T>(1e-6))
//             , preserve_flat_areas(false)
//             , validate_result(false)
//             , generate_statistics(false)
//             , max_fill_depth(std::numeric_limits<T>::infinity())
//             , max_iterations(0) // Will be set to grid_size * 2 if 0
//         {
//         }
//     };

//     /**
//      * Results structure containing flood fill statistics and validation info
//      */
//     struct Results
//     {
//         size_t nodes_processed;      // Total nodes processed
//         size_t nodes_filled;             // Nodes that had elevation
//         increased size_t depressions_filled; // Number of separate
//         depressions identified T total_fill_volume;             // Total
//         volume of material added T max_fill_depth;                    //
//         Maximum fill depth encountered T min_elevation; // Minimum elevation
//         after filling T max_elevation;                     // Maximum
//         elevation after filling size_t iterations;               // Number of
//         algorithm iterations bool convergence_reached;    // Whether
//         algorithm converged normally std::vector<std::string> warnings; //
//         Any warnings generated

//         Results()
//             : nodes_processed(0)
//             , nodes_filled(0)
//             , depressions_filled(0)
//             , total_fill_volume(static_cast<T>(0))
//             , max_fill_depth(static_cast<T>(0))
//             , min_elevation(std::numeric_limits<T>::infinity())
//             , max_elevation(-std::numeric_limits<T>::infinity())
//             , iterations(0)
//             , convergence_reached(false)
//         {
//         }
//     };

// private:
//     // OPTIMIZATION 1: Replace FloodNode with simple pair
//     using PQElement = std::pair<T, size_t>; // elevation, index
//     using PriorityQueue = std::priority_queue<PQElement,
//     std::vector<PQElement>, std::greater<PQElement>>;

// public:
//     /**
//      * Main priority flood function - fills depressions in elevation grid
//      *
//      * @param elevation Grid2D containing elevation data (modified in place)
//      * @param connector Connector for neighbor operations and boundary
//      handling
//      * @param config Configuration parameters for the algorithm
//      * @return Results structure with statistics and validation info
//      */
//     static Results fill_depressions(Grid2D<T>& elevation,
//                                                                     const
//                                                                     Connector<T>&
//                                                                     connector,
//                                                                     const
//                                                                     Config&
//                                                                     config =
//                                                                     Config())
//     {
//         Results results;

//         // Validate inputs
//         if (elevation.rows() != connector.rows() ||
//                 elevation.cols() != connector.cols()) {
//             results.warnings.push_back("Grid and connector dimension
//             mismatch"); return results;
//         }

//         if (elevation.size() == 0) {
//             results.warnings.push_back("Empty elevation grid");
//             return results;
//         }

//         // Set up algorithm parameters
//         Config working_config = config;
//         if (working_config.max_iterations == 0) {
//             working_config.max_iterations = elevation.size() * 2;
//         }

//         // OPTIMIZATION 4: Pre-allocate vectors
//         PriorityQueue pq;
//         std::vector<bool> processed(elevation.size(), false);
//         std::vector<T> original_elevation;

//         if (working_config.generate_statistics) {
//             original_elevation.assign(elevation.begin(), elevation.end());
//         }

//         // Phase 1: Initialize boundary seeds (OPTIMIZATION 2)
//         initialize_boundary_seeds_optimized(elevation, connector, pq,
//         processed, results);

//         if (pq.empty()) {
//             results.warnings.push_back(
//                 "No valid boundary seeds found - cannot perform flood fill");
//             return results;
//         }

//         // Phase 2: Priority flood algorithm (OPTIMIZATION 3)
//         execute_priority_flood_optimized(
//             elevation, connector, pq, processed, working_config, results);

//         // Phase 3: Post-processing and validation
//         if (working_config.validate_result) {
//             validate_result(elevation, connector, results);
//         }

//         if (working_config.generate_statistics) {
//             compute_statistics(elevation, original_elevation, results);
//         }

//         return results;
//     }

//     /**
//      * Simplified interface with default configuration
//      */
//     static Results fill_depressions(Grid2D<T>& elevation,
//                                                                     const
//                                                                     Connector<T>&
//                                                                     connector)
//     {
//         return fill_depressions(elevation, connector, Config());
//     }

//     /**
//      * Fill depressions with custom epsilon value
//      */
//     static Results fill_depressions(Grid2D<T>& elevation,
//                                                                     const
//                                                                     Connector<T>&
//                                                                     connector,
//                                                                     T
//                                                                     epsilon)
//     {
//         Config config;
//         config.epsilon = epsilon;
//         return fill_depressions(elevation, connector, config);
//     }

//     /**
//      * Create a backup copy of elevation before filling
//      */
//     static std::unique_ptr<Grid2D<T>> backup_elevation(const Grid2D<T>&
//     elevation)
//     {
//         std::vector<T> backup_data(elevation.begin(), elevation.end());
//         return std::make_unique<Grid2D<T>>(
//             backup_data, elevation.rows(), elevation.cols());
//     }

//     /**
//      * Restore elevation from backup
//      */
//     static void restore_elevation(Grid2D<T>& elevation, const Grid2D<T>&
//     backup)
//     {
//         if (elevation.size() != backup.size()) {
//             throw std::invalid_argument("Elevation and backup size
//             mismatch");
//         }

//         std::copy(backup.begin(), backup.end(), elevation.begin());
//     }

// private:
//     /**
//      * OPTIMIZATION 2: Only check perimeter nodes for boundary seeds
//      */
//     static void initialize_boundary_seeds_optimized(const Grid2D<T>&
//     elevation,
//                                                                                                     const Connector<T>& connector,
//                                                                                                     PriorityQueue& pq,
//                                                                                                     std::vector<bool>& processed,
//                                                                                                     Results& results)
//     {
//         const size_t rows = elevation.rows();
//         const size_t cols = elevation.cols();

//         // Pre-allocate boundary seeds vector
//         std::vector<size_t> boundary_candidates;
//         boundary_candidates.reserve(2 * (rows + cols));

//         // Only check actual perimeter nodes
//         for (size_t r = 0; r < rows; ++r) {
//             if (r == 0 || r == rows - 1) {
//                 // Top and bottom rows - check all columns
//                 for (size_t c = 0; c < cols; ++c) {
//                     boundary_candidates.push_back(connector.to_1d(r, c));
//                 }
//             } else {
//                 // Middle rows - only check left and right edges
//                 boundary_candidates.push_back(connector.to_1d(r, 0));
//                 boundary_candidates.push_back(connector.to_1d(r, cols - 1));
//             }
//         }

//         // Process boundary candidates
//         for (size_t idx : boundary_candidates) {
//             auto [r, c] = connector.to_2d(idx);
//             NodeType node_type = connector.get_boundary_type(r, c);

//             // Skip NO_DATA nodes
//             if (node_type == NodeType::NO_DATA) {
//                 processed[idx] = true;
//                 continue;
//             }

//             // Only add nodes that allow outflow
//             if (NodeTypeUtils::allows_outflow(node_type)) {
//                 T elev = elevation(r, c);
//                 pq.emplace(elev, idx);
//                 processed[idx] = true;
//             }
//         }
//     }

//     /**
//      * OPTIMIZATION 3: Direct neighbor access with unrolled loops
//      */
//     static void execute_priority_flood_optimized(Grid2D<T>& elevation,
//                                                                                                 const Connector<T>& connector,
//                                                                                                 PriorityQueue& pq,
//                                                                                                 std::vector<bool>& processed,
//                                                                                                 const Config& config,
//                                                                                                 Results& results)
//     {
//         results.iterations = 0;
//         const size_t rows = elevation.rows();
//         const size_t cols = elevation.cols();

//         while (!pq.empty() && results.iterations < config.max_iterations) {
//             auto [current_elev, current_idx] = pq.top();
//             pq.pop();

//             results.iterations++;
//             results.nodes_processed++;

//             // Update elevation if needed
//             T original_elev = elevation[current_idx];
//             if (current_elev > original_elev) {
//                 T fill_depth = current_elev - original_elev;

//                 // Check maximum fill depth constraint
//                 if (fill_depth > config.max_fill_depth) {
//                     results.warnings.push_back("Maximum fill depth exceeded
//                     at node " +
//                                                                          std::to_string(current_idx));
//                     continue;
//                 }

//                 elevation[current_idx] = current_elev;
//                 results.nodes_filled++;
//                 results.total_fill_volume += fill_depth;
//                 results.max_fill_depth = std::max(results.max_fill_depth,
//                 fill_depth);
//             }

//             // OPTIMIZATION 3: Direct neighbor processing with manual bounds
//             checking auto [row, col] = connector.to_2d(current_idx);

//             // Unrolled D8 neighbor access for maximum performance
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row - 1, col, rows, cols, current_elev, config); // North
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row, col + 1, rows, cols, current_elev, config); // East
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row + 1, col, rows, cols, current_elev, config); // South
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row, col - 1, rows, cols, current_elev, config); // West
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row - 1, col + 1, rows, cols, current_elev, config); // NE
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row + 1, col + 1, rows, cols, current_elev, config); // SE
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row + 1, col - 1, rows, cols, current_elev, config); // SW
//             process_neighbor_if_valid(elevation, connector, pq, processed,
//             row - 1, col - 1, rows, cols, current_elev, config); // NW
//         }

//         results.convergence_reached = pq.empty();

//         if (results.iterations >= config.max_iterations) {
//             results.warnings.push_back(
//                 "Maximum iterations reached - algorithm may not have
//                 converged");
//         }
//     }

//     /**
//      * OPTIMIZATION 3: Inline neighbor processing
//      */
//     static inline void process_neighbor_if_valid(Grid2D<T>& elevation,
//                                                                                              const Connector<T>& connector,
//                                                                                              PriorityQueue& pq,
//                                                                                              std::vector<bool>& processed,
//                                                                                              size_t nr, size_t nc,
//                                                                                              size_t rows, size_t cols,
//                                                                                              T current_elev,
//                                                                                              const Config& config)
//     {
//         // Bounds check
//         if (nr >= rows || nc >= cols) return;

//         size_t neighbor_idx = connector.to_1d(nr, nc);
//         if (processed[neighbor_idx]) return;

//         // Check boundary type
//         NodeType boundary_type = connector.get_boundary_type(nr, nc);
//         if (boundary_type == NodeType::NO_DATA) return;

//         T neighbor_elev = elevation[neighbor_idx];
//         T required_fill_elevation = current_elev;

//         // Apply epsilon increment if needed
//         if (neighbor_elev <= current_elev) {
//             if (config.preserve_flat_areas &&
//                     std::abs(neighbor_elev - current_elev) < config.epsilon)
//                     {
//                 required_fill_elevation = current_elev;
//             } else {
//                 required_fill_elevation = current_elev + config.epsilon;
//             }
//         } else {
//             required_fill_elevation = neighbor_elev;
//         }

//         // Add to priority queue
//         pq.emplace(required_fill_elevation, neighbor_idx);
//         processed[neighbor_idx] = true;
//     }

//     /**
//      * Validate the flood fill result
//      */
//     static void validate_result(const Grid2D<T>& elevation,
//                                                             const
//                                                             Connector<T>&
//                                                             connector,
//                                                             Results& results)
//     {
//         size_t invalid_drainage = 0;
//         const size_t rows = elevation.rows();
//         const size_t cols = elevation.cols();

//         for (size_t r = 0; r < rows; ++r) {
//             for (size_t c = 0; c < cols; ++c) {
//                 NodeType node_type = connector.get_boundary_type(r, c);

//                 // Skip NO_DATA nodes
//                 if (node_type == NodeType::NO_DATA) {
//                     continue;
//                 }

//                 T center_elev = elevation(r, c);
//                 auto neighbors = connector.get_effective_valid_neighbors(r,
//                 c);

//                 bool has_downhill_path = false;
//                 bool is_outlet = NodeTypeUtils::allows_outflow(node_type);

//                 for (const auto& neighbor : neighbors) {
//                     if (neighbor.boundary_type != NodeType::NO_DATA) {
//                         T neighbor_elev = elevation[neighbor.index];
//                         if (neighbor_elev < center_elev) {
//                             has_downhill_path = true;
//                             break;
//                         }
//                     }
//                 }

//                 // Nodes should either have a downhill path or be outlets
//                 if (!has_downhill_path && !is_outlet) {
//                     invalid_drainage++;
//                 }

//                 // Update elevation bounds
//                 results.min_elevation = std::min(results.min_elevation,
//                 center_elev); results.max_elevation =
//                 std::max(results.max_elevation, center_elev);
//             }
//         }

//         if (invalid_drainage > 0) {
//             results.warnings.push_back(
//                 "Found " + std::to_string(invalid_drainage) +
//                 " nodes without proper drainage after filling");
//         }
//     }

//     /**
//      * Compute detailed statistics about the flood fill operation
//      */
//     static void compute_statistics(const Grid2D<T>& elevation,
//                                                                  const
//                                                                  std::vector<T>&
//                                                                  original_elevation,
//                                                                  Results&
//                                                                  results)
//     {
//         if (original_elevation.size() != elevation.size()) {
//             results.warnings.push_back(
//                 "Cannot compute statistics - original elevation data
//                 missing");
//             return;
//         }

//         // Count separate depressions (connected components of filled areas)
//         std::vector<bool> visited(elevation.size(), false);
//         results.depressions_filled = 0;

//         for (size_t i = 0; i < elevation.size(); ++i) {
//             if (!visited[i] && elevation[i] > original_elevation[i]) {
//                 // Found start of a new depression - trace it
//                 results.depressions_filled++;
//                 trace_depression(i, elevation, original_elevation, visited);
//             }
//         }
//     }

//     /**
//      * Trace a single depression for statistics
//      */
//     static void trace_depression(size_t start_index,
//                                                              const Grid2D<T>&
//                                                              elevation, const
//                                                              std::vector<T>&
//                                                              original_elevation,
//                                                              std::vector<bool>&
//                                                              visited)
//     {
//         std::queue<size_t> queue;
//         queue.push(start_index);
//         visited[start_index] = true;

//         const size_t rows = elevation.rows();
//         const size_t cols = elevation.cols();

//         while (!queue.empty()) {
//             size_t current = queue.front();
//             queue.pop();

//             auto [r, c] = std::make_pair(current / cols, current % cols);

//             // Check 4-connected neighbors for depression extent
//             static const std::array<std::pair<int, int>, 4> directions = {
//                 { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } }
//             };

//             for (const auto& [dr, dc] : directions) {
//                 int nr = static_cast<int>(r) + dr;
//                 int nc = static_cast<int>(c) + dc;

//                 if (nr >= 0 && nr < static_cast<int>(rows) && nc >= 0 &&
//                         nc < static_cast<int>(cols)) {

//                     size_t neighbor_idx =
//                         static_cast<size_t>(nr) * cols +
//                         static_cast<size_t>(nc);

//                     if (!visited[neighbor_idx] &&
//                             elevation[neighbor_idx] >
//                             original_elevation[neighbor_idx]) {
//                         visited[neighbor_idx] = true;
//                         queue.push(neighbor_idx);
//                     }
//                 }
//             }
//         }
//     }

// public:
//     /**
//      * Utility functions for working with flood fill results
//      */

//     /**
//      * Check if elevation grid has proper drainage (no internal sinks)
//      */
//     static bool has_proper_drainage(const Grid2D<T>& elevation,
//                                                                     const
//                                                                     Connector<T>&
//                                                                     connector)
//     {
//         Results dummy_results;
//         validate_result(elevation, connector, dummy_results);
//         return dummy_results.warnings.empty();
//     }

//     /**
//      * Find all remaining sinks in elevation data
//      */
//     static std::vector<size_t> find_sinks(const Grid2D<T>& elevation,
//                                                                                 const Connector<T>& connector)
//     {
//         std::vector<size_t> sinks;
//         const size_t rows = elevation.rows();
//         const size_t cols = elevation.cols();

//         for (size_t r = 0; r < rows; ++r) {
//             for (size_t c = 0; c < cols; ++c) {
//                 NodeType node_type = connector.get_boundary_type(r, c);

//                 if (node_type == NodeType::NO_DATA) {
//                     continue;
//                 }

//                 size_t idx = connector.to_1d(r, c);
//                 T center_elev = elevation(r, c);
//                 auto neighbors = connector.get_effective_valid_neighbors(r,
//                 c);

//                 bool is_sink = true;
//                 bool is_outlet = NodeTypeUtils::allows_outflow(node_type);

//                 if (is_outlet) {
//                     is_sink = false;
//                 } else {
//                     for (const auto& neighbor : neighbors) {
//                         if (neighbor.boundary_type != NodeType::NO_DATA) {
//                             T neighbor_elev = elevation[neighbor.index];
//                             if (neighbor_elev < center_elev) {
//                                 is_sink = false;
//                                 break;
//                             }
//                         }
//                     }
//                 }

//                 if (is_sink) {
//                     sinks.push_back(idx);
//                 }
//             }
//         }

//         return sinks;
//     }

//     /**
//      * Estimate required epsilon for proper drainage
//      */
//     static T estimate_epsilon(const Grid2D<T>& elevation,
//                                                         const Connector<T>&
//                                                         connector)
//     {
//         T min_elevation_diff = std::numeric_limits<T>::infinity();
//         const size_t rows = elevation.rows();
//         const size_t cols = elevation.cols();

//         for (size_t r = 0; r < rows; ++r) {
//             for (size_t c = 0; c < cols; ++c) {
//                 NodeType node_type = connector.get_boundary_type(r, c);

//                 if (node_type == NodeType::NO_DATA) {
//                     continue;
//                 }

//                 T center_elev = elevation(r, c);
//                 auto neighbors = connector.get_effective_valid_neighbors(r,
//                 c);

//                 for (const auto& neighbor : neighbors) {
//                     if (neighbor.boundary_type != NodeType::NO_DATA) {
//                         T neighbor_elev = elevation[neighbor.index];
//                         T diff = std::abs(neighbor_elev - center_elev);
//                         if (diff > std::numeric_limits<T>::epsilon()) {
//                             min_elevation_diff = std::min(min_elevation_diff,
//                             diff);
//                         }
//                     }
//                 }
//             }
//         }

//         // Return 1/1000th of minimum elevation difference, with reasonable
//         bounds T estimated_epsilon = min_elevation_diff /
//         static_cast<T>(1000); return std::max(static_cast<T>(1e-9),
//                                         std::min(estimated_epsilon,
//                                         static_cast<T>(1e-3)));
//     }

//     /**
//      * Generate detailed report of flood fill results
//      */
//     static std::string generate_report(const Results& results)
//     {
//         std::ostringstream report;

//         report << "Priority Flood Results Summary:\n";
//         report << "==============================\n";
//         report << "Nodes processed: " << results.nodes_processed << "\n";
//         report << "Nodes filled: " << results.nodes_filled << "\n";
//         report << "Depressions filled: " << results.depressions_filled <<
//         "\n"; report << "Total fill volume: " << results.total_fill_volume <<
//         "\n"; report << "Maximum fill depth: " << results.max_fill_depth <<
//         "\n"; report << "Elevation range: [" << results.min_elevation << ", "
//                      << results.max_elevation << "]\n";
//         report << "Algorithm iterations: " << results.iterations << "\n";
//         report << "Convergence: " << (results.convergence_reached ? "Yes" :
//         "No")
//                      << "\n";

//         if (!results.warnings.empty()) {
//             report << "\nWarnings:\n";
//             for (const auto& warning : results.warnings) {
//                 report << "  - " << warning << "\n";
//             }
//         }

//         return report.str();
//     }
// };

// // Convenience type aliases
// using PriorityFloodF32 = PriorityFlood<float>;
// using PriorityFloodF64 = PriorityFlood<double>;

// } // namespace dagger2
