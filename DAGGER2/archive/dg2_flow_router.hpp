#pragma once

#include "dg2_BCs.hpp"
#include "dg2_LM.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include "dg2_cordonnier_method.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace dagger2 {

/**
 * Flow routing method enumeration
 */
enum class FlowRoutingMethod : uint8_t
{
	// Single flow direction methods
	D8_STEEPEST = 0, // D8 steepest descent
	D4_STEEPEST = 1, // D4 steepest descent
	DINF = 2,				 // D-infinity (Tarboton 1997)

	// Multiple flow direction methods
	MFD_FREEMAN = 3,	// Freeman (1991) multiple flow
	MFD_QUINN = 4,		// Quinn et al. (1991)
	MFD_HOLMGREN = 5, // Holmgren (1994)
	MFD_SEIBERT = 6,	// Seibert & McGlynn (2007)

	// Advanced methods
	ADAPTIVE_MFD = 7,			 // Adaptive multiple flow
	KINEMATIC_ROUTING = 8, // Kinematic wave routing

	// Cordonnier basin graph methods
	CORDONNIER_SIMPLE = 9,	 // Simple correction
	CORDONNIER_CARVING = 10, // Carving approach
	CORDONNIER_FILLING = 11	 // Filling approach
};

/**
 * Flow routing configuration
 */
template<typename T>
struct FlowRoutingConfig
{
	FlowRoutingMethod method;
	T flow_exponent;								 // For MFD methods (typically 1.1)
	T min_gradient;									 // Minimum gradient threshold
	bool handle_depressions;				 // Whether to handle local minima
	FillingMethod depression_method; // How to handle depressions
	bool compute_accumulation;			 // Whether to compute flow accumulation
	bool compute_drainage_area;			 // Whether to compute drainage area
	T cell_area;										 // Area of each cell
	bool enforce_monotonic;					 // Enforce monotonic flow paths
	uint32_t random_seed;						 // For stochastic methods

	FlowRoutingConfig()
		: method(FlowRoutingMethod::D8_STEEPEST)
		, flow_exponent(1.1)
		, min_gradient(1e-8)
		, handle_depressions(true)
		, depression_method(FillingMethod::PRIORITY_FLOOD_EPSILON)
		, compute_accumulation(true)
		, compute_drainage_area(true)
		, cell_area(1.0)
		, enforce_monotonic(true)
		, random_seed(42)
	{
	}
};

/**
 * Flow routing results
 */
template<typename T>
struct FlowRoutingResult
{
	bool success;
	std::shared_ptr<ArrayRef<size_t>> receivers; // Single flow receivers
	std::shared_ptr<ArrayRef<std::vector<size_t>>>
		donors;																			 // Multiple donors per node
	std::shared_ptr<ArrayRef<size_t>> stack_order; // Topological order
	std::shared_ptr<ArrayRef<size_t>> flow_accumulation; // Flow accumulation
	std::shared_ptr<ArrayRef<T>> drainage_area;					 // Drainage area
	std::shared_ptr<ArrayRef<T>> flow_weights;					 // Flow weights (for MFD)
	std::shared_ptr<ArrayRef<T>> gradients;							 // Flow gradients
	FillingResult<T> depression_result; // Depression filling results
	std::chrono::milliseconds processing_time;

	// Multiple flow results (for MFD methods)
	std::shared_ptr<ArrayRef<std::vector<std::pair<size_t, T>>>>
		mfd_receivers; // (receiver, weight) pairs

	// 2D versions (when using Grid2D input)
	std::shared_ptr<Grid2D<size_t>> receivers_2d;
	std::shared_ptr<Grid2D<std::vector<size_t>>> donors_2d;
	std::shared_ptr<Grid2D<size_t>> stack_order_2d;
	std::shared_ptr<Grid2D<size_t>> flow_accumulation_2d;
	std::shared_ptr<Grid2D<T>> drainage_area_2d;
	std::shared_ptr<Grid2D<T>> flow_weights_2d;
	std::shared_ptr<Grid2D<T>> gradients_2d;
	std::shared_ptr<Grid2D<std::vector<std::pair<size_t, T>>>> mfd_receivers_2d;

	bool is_2d = false;
	size_t rows = 0;
	size_t cols = 0;

	FlowRoutingResult()
		: success(false)
		, processing_time(0)
	{
	}

	// Convenience methods for accessing data
	size_t get_receiver(size_t idx) const
	{
		return receivers ? (*receivers)[idx] : SIZE_MAX;
	}

	size_t get_receiver_2d(size_t row, size_t col) const
	{
		return receivers_2d ? (*receivers_2d)(row, col) : SIZE_MAX;
	}

	size_t get_flow_accumulation(size_t idx) const
	{
		return flow_accumulation ? (*flow_accumulation)[idx] : 1;
	}

	size_t get_flow_accumulation_2d(size_t row, size_t col) const
	{
		return flow_accumulation_2d ? (*flow_accumulation_2d)(row, col) : 1;
	}

	T get_drainage_area(size_t idx) const
	{
		return drainage_area ? (*drainage_area)[idx] : static_cast<T>(1);
	}

	T get_drainage_area_2d(size_t row, size_t col) const
	{
		return drainage_area_2d ? (*drainage_area_2d)(row, col) : static_cast<T>(1);
	}

	// Convert between 1D and 2D indices
	size_t to_1d(size_t row, size_t col) const { return row * cols + col; }

	std::pair<size_t, size_t> to_2d(size_t idx) const
	{
		return { idx / cols, idx % cols };
	}

	// Sync 1D and 2D data
	void sync_to_2d()
	{
		if (!is_2d || rows == 0 || cols == 0)
			return;

		// Create 2D versions if they don't exist
		if (receivers && !receivers_2d) {
			receivers_2d = std::make_shared<Grid2D<size_t>>(rows, cols);
			for (size_t i = 0; i < receivers->size(); ++i) {
				auto [r, c] = to_2d(i);
				(*receivers_2d)(r, c) = (*receivers)[i];
			}
		}

		if (flow_accumulation && !flow_accumulation_2d) {
			flow_accumulation_2d = std::make_shared<Grid2D<size_t>>(rows, cols);
			for (size_t i = 0; i < flow_accumulation->size(); ++i) {
				auto [r, c] = to_2d(i);
				(*flow_accumulation_2d)(r, c) = (*flow_accumulation)[i];
			}
		}

		if (drainage_area && !drainage_area_2d) {
			drainage_area_2d = std::make_shared<Grid2D<T>>(rows, cols);
			for (size_t i = 0; i < drainage_area->size(); ++i) {
				auto [r, c] = to_2d(i);
				(*drainage_area_2d)(r, c) = (*drainage_area)[i];
			}
		}

		if (gradients && !gradients_2d) {
			gradients_2d = std::make_shared<Grid2D<T>>(rows, cols);
			for (size_t i = 0; i < gradients->size(); ++i) {
				auto [r, c] = to_2d(i);
				(*gradients_2d)(r, c) = (*gradients)[i];
			}
		}

		// Add other arrays as needed...
	}
};

/**
 * Comprehensive Flow Router
 *
 * This class provides a unified interface for all flow routing methods
 * in the DAGGER2 framework, including depression handling and various
 * single/multiple flow direction algorithms.
 */
template<typename T = double>
class FlowRouter
{
private:
	// Core components
	std::shared_ptr<Connector<T>> connector_;
	std::shared_ptr<ArrayRef<T>> elevation_;
	std::shared_ptr<Grid2D<T>> elevation_2d_;
	FlowRoutingConfig<T> config_;

	// Grid properties
	size_t rows_;
	size_t cols_;
	size_t size_;
	bool is_2d_;

	// Working elevation (after depression filling)
	mutable std::shared_ptr<ArrayRef<T>> working_elevation_;

	// Depression filler and Cordonnier router
	mutable std::unique_ptr<DepressionFiller<T>> depression_filler_;
	mutable std::unique_ptr<CordonnierFlowRouter<T>> cordonnier_router_;

	// Random number generator
	mutable std::mt19937 rng_;

public:
	/**
	 * Constructor for 1D ArrayRef input
	 */
	FlowRouter(std::shared_ptr<Connector<T>> connector,
						 std::shared_ptr<ArrayRef<T>> elevation,
						 const FlowRoutingConfig<T>& config = FlowRoutingConfig<T>())
		: connector_(connector)
		, elevation_(elevation)
		, config_(config)
		, is_2d_(false)
		, rng_(config.random_seed)
	{
		if (!connector_) {
			throw std::invalid_argument("Connector cannot be null");
		}

		if (!elevation_) {
			throw std::invalid_argument("Elevation array cannot be null");
		}

		rows_ = connector_->rows();
		cols_ = connector_->cols();
		size_ = connector_->size();

		if (elevation_->size() != size_) {
			throw std::invalid_argument("Elevation array size must match grid size");
		}

		// Initialize working elevation
		working_elevation_ = std::make_shared<ArrayRef<T>>(*elevation_);
	}

	/**
	 * Constructor for 2D Grid2D input
	 */
	FlowRouter(std::shared_ptr<Connector<T>> connector,
						 std::shared_ptr<Grid2D<T>> elevation,
						 const FlowRoutingConfig<T>& config = FlowRoutingConfig<T>())
		: connector_(connector)
		, elevation_2d_(elevation)
		, config_(config)
		, is_2d_(true)
		, rng_(config.random_seed)
	{
		if (!connector_) {
			throw std::invalid_argument("Connector cannot be null");
		}

		if (!elevation) {
			throw std::invalid_argument("Elevation grid cannot be null");
		}

		rows_ = connector_->rows();
		cols_ = connector_->cols();
		size_ = connector_->size();

		if (elevation_2d_->rows() != rows_ || elevation_2d_->cols() != cols_) {
			throw std::invalid_argument(
				"Elevation grid dimensions must match connector");
		}

		// Convert Grid2D to ArrayRef for internal processing
		std::vector<T> temp_data(size_);
		for (size_t i = 0; i < rows_; ++i) {
			for (size_t j = 0; j < cols_; ++j) {
				temp_data[i * cols_ + j] = (*elevation_2d_)(i, j);
			}
		}
		elevation_ = std::make_shared<ArrayRef<T>>(temp_data);

		// Initialize working elevation
		working_elevation_ = std::make_shared<ArrayRef<T>>(*elevation_);
	}

	// ======================
	// CONFIGURATION
	// ======================

	const FlowRoutingConfig<T>& get_config() const { return config_; }
	void set_config(const FlowRoutingConfig<T>& config)
	{
		config_ = config;
		rng_.seed(config_.random_seed);
	}

	void set_method(FlowRoutingMethod method) { config_.method = method; }
	void set_depression_handling(
		bool enable,
		FillingMethod method = FillingMethod::PRIORITY_FLOOD_EPSILON)
	{
		config_.handle_depressions = enable;
		config_.depression_method = method;
	}

	bool is_2d() const { return is_2d_; }
	size_t rows() const { return rows_; }
	size_t cols() const { return cols_; }
	size_t size() const { return size_; }

	// ======================
	// MAIN FLOW ROUTING
	// ======================

	/**
	 * Compute flow routing using configured method
	 */
	FlowRoutingResult<T> compute_flow_routing() const
	{
		auto start_time = std::chrono::high_resolution_clock::now();

		FlowRoutingResult<T> result;
		result.is_2d = is_2d_;
		result.rows = rows_;
		result.cols = cols_;

		try {
			// Step 1: Handle depressions if enabled
			if (config_.handle_depressions) {
				result.depression_result = handle_depressions();
				if (!result.depression_result.success) {
					result.success = false;
					return result;
				}
			}

			// Step 2: Compute flow routing based on method
			if (is_cordonnier_method(config_.method)) {
				compute_cordonnier_flow_routing(result);
			} else if (is_mfd_method(config_.method)) {
				compute_mfd_flow_routing(result);
			} else {
				compute_sfd_flow_routing(result);
			}

			// Step 3: Compute accumulation and drainage area if requested
			if (config_.compute_accumulation) {
				compute_flow_accumulation(result);
			}

			if (config_.compute_drainage_area) {
				compute_drainage_area(result);
			}

			// Step 4: Create 2D versions if input was 2D
			if (is_2d_) {
				result.sync_to_2d();
			}

			result.success = true;

		} catch (const std::exception& e) {
			result.success = false;
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		result.processing_time =
			std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
																														start_time);

		return result;
	}

	// ======================
	// CONVENIENCE METHODS
	// ======================

	/**
	 * Quick flow accumulation computation
	 */
	std::shared_ptr<ArrayRef<size_t>> compute_flow_accumulation_quick() const
	{
		FlowRoutingConfig<T> quick_config = config_;
		quick_config.compute_accumulation = true;
		quick_config.compute_drainage_area = false;

		FlowRouter<T> router(connector_, elevation_, quick_config);
		auto result = router.compute_flow_routing();

		return result.flow_accumulation;
	}

	/**
	 * Quick drainage area computation
	 */
	std::shared_ptr<ArrayRef<T>> compute_drainage_area_quick() const
	{
		FlowRoutingConfig<T> quick_config = config_;
		quick_config.compute_accumulation = false;
		quick_config.compute_drainage_area = true;

		FlowRouter<T> router(connector_, elevation_, quick_config);
		auto result = router.compute_flow_routing();

		return result.drainage_area;
	}

	/**
	 * Quick single flow receivers
	 */
	std::shared_ptr<ArrayRef<size_t>> compute_receivers() const
	{
		FlowRoutingConfig<T> quick_config = config_;
		quick_config.compute_accumulation = false;
		quick_config.compute_drainage_area = false;

		FlowRouter<T> router(connector_, elevation_, quick_config);
		auto result = router.compute_flow_routing();

		return result.receivers;
	}

	// ======================
	// ANALYSIS UTILITIES
	// ======================

	/**
	 * Generate flow routing report
	 */
	std::string generate_report(const FlowRoutingResult<T>& result) const
	{
		std::ostringstream report;

		report << "=== FLOW ROUTING REPORT ===\n\n";

		// Method information
		report << "Method: ";
		switch (config_.method) {
			case FlowRoutingMethod::D8_STEEPEST:
				report << "D8 Steepest Descent\n";
				break;
			case FlowRoutingMethod::D4_STEEPEST:
				report << "D4 Steepest Descent\n";
				break;
			case FlowRoutingMethod::MFD_FREEMAN:
				report << "Multiple Flow Direction (Freeman 1991)\n";
				break;
			case FlowRoutingMethod::MFD_QUINN:
				report << "Multiple Flow Direction (Quinn et al. 1991)\n";
				break;
			case FlowRoutingMethod::CORDONNIER_FILLING:
				report << "Cordonnier Basin Graph (Filling)\n";
				break;
			default:
				report << "Custom Method\n";
				break;
		}

		// Configuration
		if (is_mfd_method(config_.method)) {
			report << "Flow Exponent: " << config_.flow_exponent << "\n";
		}
		report << "Minimum Gradient: " << config_.min_gradient << "\n";
		report << "Depression Handling: "
					 << (config_.handle_depressions ? "Enabled" : "Disabled") << "\n";

		// Results summary
		report << "\n=== RESULTS SUMMARY ===\n";
		report << "Processing Time: " << result.processing_time.count() << " ms\n";
		report << "Success: " << (result.success ? "Yes" : "No") << "\n";
		report << "Data Format: " << (result.is_2d ? "2D Grid" : "1D Array")
					 << "\n";

		if (result.success) {
			// Flow statistics
			size_t nodes_with_receivers = 0;
			size_t source_nodes = 0;
			size_t sink_nodes = 0;

			if (result.receivers) {
				for (size_t i = 0; i < size_; ++i) {
					if (!connector_->is_active_node(i))
						continue;

					if ((*result.receivers)[i] != SIZE_MAX) {
						nodes_with_receivers++;
					}

					if (result.donors && (*result.donors)[i].empty()) {
						source_nodes++;
					}

					if ((*result.receivers)[i] == SIZE_MAX) {
						sink_nodes++;
					}
				}
			}

			report << "Nodes with Receivers: " << nodes_with_receivers << "\n";
			report << "Source Nodes: " << source_nodes << "\n";
			report << "Sink Nodes: " << sink_nodes << "\n";

			// Accumulation statistics
			if (result.flow_accumulation) {
				auto min_it = std::min_element(result.flow_accumulation->begin(),
																			 result.flow_accumulation->end());
				auto max_it = std::max_element(result.flow_accumulation->begin(),
																			 result.flow_accumulation->end());
				report << "Flow Accumulation Range: " << *min_it << " - " << *max_it
							 << "\n";

				// Calculate total accumulated flow
				size_t total_flow = std::accumulate(result.flow_accumulation->begin(),
																						result.flow_accumulation->end(),
																						static_cast<size_t>(0));
				report << "Total Accumulated Flow: " << total_flow << "\n";
			}

			// Drainage area statistics
			if (result.drainage_area) {
				auto min_it = std::min_element(result.drainage_area->begin(),
																			 result.drainage_area->end());
				auto max_it = std::max_element(result.drainage_area->begin(),
																			 result.drainage_area->end());
				report << "Drainage Area Range: " << *min_it << " - " << *max_it
							 << "\n";

				T total_area = std::accumulate(result.drainage_area->begin(),
																			 result.drainage_area->end(),
																			 static_cast<T>(0));
				report << "Total Drainage Area: " << total_area << "\n";
			}
		}

		// Depression handling results
		if (config_.handle_depressions && result.depression_result.success) {
			report << "\n=== DEPRESSION HANDLING ===\n";
			report << "Depressions Filled: "
						 << result.depression_result.depressions_filled << "\n";
			report << "Cells Modified: " << result.depression_result.cells_modified
						 << "\n";
			report << "Volume Filled: "
						 << result.depression_result.total_volume_filled << "\n";
			report << "Volume Carved: "
						 << result.depression_result.total_volume_carved << "\n";
			report << "Max Modification: "
						 << result.depression_result.max_modification << "\n";
		}

		return report.str();
	}

	/**
	 * Validate flow routing results
	 */
	bool validate_results(const FlowRoutingResult<T>& result) const
	{
		if (!result.success)
			return false;

		// Check array sizes
		if (result.receivers && result.receivers->size() != size_)
			return false;
		if (result.donors && result.donors->size() != size_)
			return false;
		if (result.flow_accumulation && result.flow_accumulation->size() != size_)
			return false;
		if (result.drainage_area && result.drainage_area->size() != size_)
			return false;

		// Check flow consistency
		if (result.receivers && result.donors) {
			for (size_t i = 0; i < size_; ++i) {
				if (!connector_->is_active_node(i))
					continue;

				size_t receiver = (*result.receivers)[i];
				if (receiver != SIZE_MAX) {
					// Check that receiver has this node as donor
					const auto& receiver_donors = (*result.donors)[receiver];
					if (std::find(receiver_donors.begin(), receiver_donors.end(), i) ==
							receiver_donors.end()) {
						return false;
					}
				}
			}
		}

		// Check that all active nodes eventually reach a boundary
		if (result.receivers) {
			for (size_t i = 0; i < size_; ++i) {
				if (!connector_->is_active_node(i))
					continue;

				std::unordered_set<size_t> visited;
				size_t current = i;
				bool reaches_boundary = false;

				while (current != SIZE_MAX && visited.find(current) == visited.end()) {
					visited.insert(current);

					NodeType boundary_type = connector_->get_boundary_type(current);
					if (NodeTypeUtils::allows_outflow(boundary_type)) {
						reaches_boundary = true;
						break;
					}

					current = (*result.receivers)[current];
				}

				if (!reaches_boundary)
					return false;
			}
		}

		return true;
	}

	/**
	 * Export flow directions as direction map
	 */
	std::shared_ptr<ArrayRef<uint8_t>> export_flow_directions(
		const FlowRoutingResult<T>& result) const
	{
		std::vector<uint8_t> directions_vec(
			size_, static_cast<uint8_t>(Direction::INVALID));
		auto directions = std::make_shared<ArrayRef<uint8_t>>(directions_vec);

		if (!result.receivers)
			return directions;

		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i) || (*result.receivers)[i] == SIZE_MAX)
				continue;

			auto [from_row, from_col] = connector_->to_2d(i);
			auto [to_row, to_col] = connector_->to_2d((*result.receivers)[i]);

			int dr = static_cast<int>(to_row) - static_cast<int>(from_row);
			int dc = static_cast<int>(to_col) - static_cast<int>(from_col);

			Direction dir = Direction::INVALID;
			if (dr == -1 && dc == 0)
				dir = Direction::NORTH;
			else if (dr == 0 && dc == 1)
				dir = Direction::EAST;
			else if (dr == 1 && dc == 0)
				dir = Direction::SOUTH;
			else if (dr == 0 && dc == -1)
				dir = Direction::WEST;
			else if (dr == -1 && dc == 1)
				dir = Direction::NORTHEAST;
			else if (dr == 1 && dc == 1)
				dir = Direction::SOUTHEAST;
			else if (dr == 1 && dc == -1)
				dir = Direction::SOUTHWEST;
			else if (dr == -1 && dc == -1)
				dir = Direction::NORTHWEST;

			(*directions)[i] = static_cast<uint8_t>(dir);
		}

		return directions;
	}

	/**
	 * Export flow weights for MFD methods
	 */
	std::shared_ptr<ArrayRef<std::vector<T>>> export_flow_weights(
		const FlowRoutingResult<T>& result) const
	{
		std::vector<std::vector<T>> weights_vec(size_);
		auto weights = std::make_shared<ArrayRef<std::vector<T>>>(weights_vec);

		if (!result.mfd_receivers)
			return weights;

		for (size_t i = 0; i < size_; ++i) {
			(*weights)[i].reserve((*result.mfd_receivers)[i].size());
			for (const auto& [receiver, weight] : (*result.mfd_receivers)[i]) {
				(*weights)[i].push_back(weight);
			}
		}

		return weights;
	}

private:
	// ======================
	// PRIVATE IMPLEMENTATION
	// ======================

	/**
	 * Handle depressions using configured method
	 */
	FillingResult<T> handle_depressions() const
	{
		if (!depression_filler_) {
			FillingConfig<T> fill_config;
			fill_config.method = config_.depression_method;
			fill_config.epsilon =
				config_.min_gradient * 10.0; // Use slightly larger epsilon

			depression_filler_ = std::make_unique<DepressionFiller<T>>(
				connector_, elevation_, fill_config);
		}

		auto result = depression_filler_->fill_depressions();

		if (result.success) {
			// Update working elevation with filled results
			auto filled_elevation = depression_filler_->get_filled_elevation();
			working_elevation_ = std::make_shared<ArrayRef<T>>(filled_elevation);
		}

		return result;
	}

	/**
	 * Check if method is Cordonnier-based
	 */
	bool is_cordonnier_method(FlowRoutingMethod method) const
	{
		return method == FlowRoutingMethod::CORDONNIER_SIMPLE ||
					 method == FlowRoutingMethod::CORDONNIER_CARVING ||
					 method == FlowRoutingMethod::CORDONNIER_FILLING;
	}

	/**
	 * Check if method is multiple flow direction
	 */
	bool is_mfd_method(FlowRoutingMethod method) const
	{
		return method == FlowRoutingMethod::MFD_FREEMAN ||
					 method == FlowRoutingMethod::MFD_QUINN ||
					 method == FlowRoutingMethod::MFD_HOLMGREN ||
					 method == FlowRoutingMethod::MFD_SEIBERT ||
					 method == FlowRoutingMethod::ADAPTIVE_MFD;
	}

	/**
	 * Compute Cordonnier flow routing
	 */
	void compute_cordonnier_flow_routing(FlowRoutingResult<T>& result) const
	{
		if (!cordonnier_router_) {
			FlowEnforcementStrategy strategy;
			switch (config_.method) {
				case FlowRoutingMethod::CORDONNIER_SIMPLE:
					strategy = FlowEnforcementStrategy::SIMPLE_CORRECTION;
					break;
				case FlowRoutingMethod::CORDONNIER_CARVING:
					strategy = FlowEnforcementStrategy::DEPRESSION_CARVING;
					break;
				case FlowRoutingMethod::CORDONNIER_FILLING:
					strategy = FlowEnforcementStrategy::DEPRESSION_FILLING;
					break;
				default:
					strategy = FlowEnforcementStrategy::DEPRESSION_FILLING;
					break;
			}

			cordonnier_router_ = std::make_unique<CordonnierFlowRouter<T>>(
				connector_, working_elevation_, strategy);
		}

		cordonnier_router_->compute_flow_routing();

		// Convert results to ArrayRef
		auto receivers_vec = cordonnier_router_->get_receivers();
		result.receivers = std::make_shared<ArrayRef<size_t>>(receivers_vec);

		auto donors_vec = cordonnier_router_->get_donors();
		result.donors = std::make_shared<ArrayRef<std::vector<size_t>>>(donors_vec);

		auto stack_vec = cordonnier_router_->get_stack_order();
		result.stack_order = std::make_shared<ArrayRef<size_t>>(stack_vec);
	}

	/**
	 * Compute single flow direction routing
	 */
	void compute_sfd_flow_routing(FlowRoutingResult<T>& result) const
	{
		std::vector<size_t> receivers_vec(size_, SIZE_MAX);
		result.receivers = std::make_shared<ArrayRef<size_t>>(receivers_vec);

		std::vector<std::vector<size_t>> donors_vec(size_);
		result.donors = std::make_shared<ArrayRef<std::vector<size_t>>>(donors_vec);

		std::vector<T> gradients_vec(size_, static_cast<T>(0));
		result.gradients = std::make_shared<ArrayRef<T>>(gradients_vec);

		// Clear donors
		for (size_t i = 0; i < size_; ++i) {
			(*result.donors)[i].clear();
		}

		// Compute receivers
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			auto receiver_info = compute_single_receiver(i);
			(*result.receivers)[i] = receiver_info.first;
			(*result.gradients)[i] = receiver_info.second;
		}

		// Build donors from receivers
		for (size_t i = 0; i < size_; ++i) {
			if ((*result.receivers)[i] != SIZE_MAX) {
				(*result.donors)[(*result.receivers)[i]].push_back(i);
			}
		}

		// Compute topological order
		result.stack_order =
			compute_topological_order(*result.receivers, *result.donors);
	}

	/**
	 * Compute multiple flow direction routing
	 */
	void compute_mfd_flow_routing(FlowRoutingResult<T>& result) const
	{
		// Create ArrayRef objects with proper initialization
		std::vector<std::vector<std::pair<size_t, T>>> mfd_receivers_vec(size_);
		result.mfd_receivers =
			std::make_shared<ArrayRef<std::vector<std::pair<size_t, T>>>>(
				mfd_receivers_vec);

		std::vector<std::vector<size_t>> donors_vec(size_);
		result.donors = std::make_shared<ArrayRef<std::vector<size_t>>>(donors_vec);

		std::vector<T> gradients_vec(size_, static_cast<T>(0));
		result.gradients = std::make_shared<ArrayRef<T>>(gradients_vec);

		// Clear arrays
		for (size_t i = 0; i < size_; ++i) {
			(*result.mfd_receivers)[i].clear();
			(*result.donors)[i].clear();
		}

		// Compute MFD receivers
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			compute_mfd_receivers(
				i, (*result.mfd_receivers)[i], (*result.gradients)[i]);
		}

		// Build donors from MFD receivers
		for (size_t i = 0; i < size_; ++i) {
			for (const auto& [receiver, weight] : (*result.mfd_receivers)[i]) {
				(*result.donors)[receiver].push_back(i);
			}
		}

		// For compatibility, set primary receiver as steepest
		std::vector<size_t> receivers_vec(size_, SIZE_MAX);
		result.receivers = std::make_shared<ArrayRef<size_t>>(receivers_vec);

		for (size_t i = 0; i < size_; ++i) {
			if (!(*result.mfd_receivers)[i].empty()) {
				// Find receiver with highest weight
				auto max_it = std::max_element(
					(*result.mfd_receivers)[i].begin(),
					(*result.mfd_receivers)[i].end(),
					[](const auto& a, const auto& b) { return a.second < b.second; });
				(*result.receivers)[i] = max_it->first;
			}
		}

		// Compute topological order using primary receivers
		result.stack_order =
			compute_topological_order(*result.receivers, *result.donors);
	}

	/**
	 * Compute single flow receiver for a node
	 */
	std::pair<size_t, T> compute_single_receiver(size_t index) const
	{
		T center_elev = (*working_elevation_)[index];
		auto neighbors = get_valid_flow_neighbors(index);

		size_t best_receiver = SIZE_MAX;
		T max_gradient = -std::numeric_limits<T>::infinity();

		for (const auto& neighbor : neighbors) {
			T neighbor_elev = (*working_elevation_)[neighbor.index];
			T gradient = (center_elev - neighbor_elev) / neighbor.distance;

			if (gradient > max_gradient && gradient > config_.min_gradient) {
				max_gradient = gradient;
				best_receiver = neighbor.index;
			}
		}

		return { best_receiver, max_gradient };
	}

	/**
	 * Get valid flow neighbors based on routing method
	 */
	std::vector<Neighbor> get_valid_flow_neighbors(size_t index) const
	{
		auto neighbors = connector_->get_valid_neighbors(index);

		// Filter based on connectivity for certain methods
		if (config_.method == FlowRoutingMethod::D4_STEEPEST) {
			// Only cardinal directions
			neighbors.erase(std::remove_if(neighbors.begin(),
																		 neighbors.end(),
																		 [](const Neighbor& n) {
																			 return n.direction >=
																							Direction::NORTHEAST;
																		 }),
											neighbors.end());
		}

		return neighbors;
	}

	/**
	 * Compute topological order for flow processing
	 */
	std::shared_ptr<ArrayRef<size_t>> compute_topological_order(
		const ArrayRef<size_t>& receivers,
		const ArrayRef<std::vector<size_t>>& donors) const
	{
		std::vector<size_t> stack_order_vec(size_);
		auto stack_order = std::make_shared<ArrayRef<size_t>>(stack_order_vec);
		size_t stack_index = 0;

		std::vector<size_t> in_degree(size_, 0);
		std::queue<size_t> queue;

		// Calculate in-degrees
		for (size_t i = 0; i < size_; ++i) {
			in_degree[i] = donors[i].size();
		}

		// Add sources (nodes with no donors)
		for (size_t i = 0; i < size_; ++i) {
			if (connector_->is_active_node(i) && in_degree[i] == 0) {
				queue.push(i);
			}
		}

		// Topological sort
		while (!queue.empty()) {
			size_t current = queue.front();
			queue.pop();
			(*stack_order)[stack_index++] = current;

			size_t receiver = receivers[current];
			if (receiver != SIZE_MAX && connector_->is_active_node(receiver)) {
				in_degree[receiver]--;
				if (in_degree[receiver] == 0) {
					queue.push(receiver);
				}
			}
		}

		return stack_order;
	}

	/**
	 * Compute multiple flow receivers for a node
	 */
	void compute_mfd_receivers(size_t index,
														 std::vector<std::pair<size_t, T>>& receivers,
														 T& max_gradient) const
	{
		receivers.clear();
		max_gradient = 0;

		T center_elev = (*working_elevation_)[index];
		auto neighbors = get_valid_flow_neighbors(index);

		// Find all downslope neighbors
		std::vector<std::tuple<size_t, T, T>>
			downslope_neighbors; // index, gradient, distance
		T total_weight = 0;

		for (const auto& neighbor : neighbors) {
			T neighbor_elev = (*working_elevation_)[neighbor.index];
			T gradient = (center_elev - neighbor_elev) / neighbor.distance;

			if (gradient > config_.min_gradient) {
				downslope_neighbors.emplace_back(
					neighbor.index, gradient, neighbor.distance);
				max_gradient = std::max(max_gradient, gradient);
			}
		}

		if (downslope_neighbors.empty())
			return;

		// Compute flow proportions based on method
		std::vector<T> weights;
		weights.reserve(downslope_neighbors.size());

		switch (config_.method) {
			case FlowRoutingMethod::MFD_FREEMAN: {
				// Proportional to slope
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					weights.push_back(gradient);
					total_weight += gradient;
				}
				break;
			}

			case FlowRoutingMethod::MFD_QUINN: {
				// Proportional to slope^exponent
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					T weight = std::pow(gradient, config_.flow_exponent);
					weights.push_back(weight);
					total_weight += weight;
				}
				break;
			}

			case FlowRoutingMethod::MFD_HOLMGREN: {
				// Holmgren method with angle-based weighting
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					T angle = std::atan(gradient);
					T weight = std::pow(std::tan(angle), config_.flow_exponent);
					weights.push_back(weight);
					total_weight += weight;
				}
				break;
			}

			case FlowRoutingMethod::MFD_SEIBERT: {
				// Seibert & McGlynn adaptive method
				T adaptive_exp = 1.0 + 0.5 * downslope_neighbors.size();
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					T weight = std::pow(gradient, adaptive_exp);
					weights.push_back(weight);
					total_weight += weight;
				}
				break;
			}

			case FlowRoutingMethod::ADAPTIVE_MFD: {
				// Adaptive based on local slope variance
				T mean_gradient = 0;
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					mean_gradient += gradient;
				}
				mean_gradient /= downslope_neighbors.size();

				T variance = 0;
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					variance += (gradient - mean_gradient) * (gradient - mean_gradient);
				}
				variance /= downslope_neighbors.size();

				T adaptive_exp =
					1.0 + variance; // Higher variance = more concentrated flow
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					T weight = std::pow(gradient, adaptive_exp);
					weights.push_back(weight);
					total_weight += weight;
				}
				break;
			}

			default: {
				// Fall back to Freeman method
				for (const auto& [idx, gradient, dist] : downslope_neighbors) {
					weights.push_back(gradient);
					total_weight += gradient;
				}
				break;
			}
		}

		// Normalize weights and create receiver list
		if (total_weight > 0) {
			for (size_t i = 0; i < downslope_neighbors.size(); ++i) {
				T normalized_weight = weights[i] / total_weight;
				receivers.emplace_back(std::get<0>(downslope_neighbors[i]),
															 normalized_weight);
			}
		}
	}

	// ======================
	// ACCUMULATION COMPUTATION
	// ======================

	/**
	 * Compute flow accumulation
	 */
	void compute_flow_accumulation(FlowRoutingResult<T>& result) const
	{
		std::vector<size_t> accumulation_vec(size_, 1);
		result.flow_accumulation =
			std::make_shared<ArrayRef<size_t>>(accumulation_vec);

		if (is_mfd_method(config_.method)) {
			compute_mfd_flow_accumulation(result);
		} else {
			compute_sfd_flow_accumulation(result);
		}
	}

	/**
	 * Compute single flow direction accumulation
	 */
	void compute_sfd_flow_accumulation(FlowRoutingResult<T>& result) const
	{
		// Process in reverse topological order
		for (int i = static_cast<int>(size_) - 1; i >= 0; --i) {
			size_t node = (*result.stack_order)[i];
			if (!connector_->is_active_node(node))
				continue;

			size_t receiver = (*result.receivers)[node];
			if (receiver != SIZE_MAX && connector_->is_active_node(receiver)) {
				(*result.flow_accumulation)[receiver] +=
					(*result.flow_accumulation)[node];
			}
		}
	}

	/**
	 * Compute multiple flow direction accumulation
	 */
	void compute_mfd_flow_accumulation(FlowRoutingResult<T>& result) const
	{
		// Process in reverse topological order
		for (int i = static_cast<int>(size_) - 1; i >= 0; --i) {
			size_t node = (*result.stack_order)[i];
			if (!connector_->is_active_node(node))
				continue;

			// Distribute flow to all receivers based on weights
			for (const auto& [receiver, weight] : (*result.mfd_receivers)[node]) {
				if (connector_->is_active_node(receiver)) {
					(*result.flow_accumulation)[receiver] +=
						static_cast<size_t>((*result.flow_accumulation)[node] * weight);
				}
			}
		}
	}

	/**
	 * Compute drainage area
	 */
	void compute_drainage_area(FlowRoutingResult<T>& result) const
	{
		std::vector<T> drainage_area_vec(size_, config_.cell_area);
		result.drainage_area = std::make_shared<ArrayRef<T>>(drainage_area_vec);

		if (is_mfd_method(config_.method)) {
			compute_mfd_drainage_area(result);
		} else {
			compute_sfd_drainage_area(result);
		}
	}

	/**
	 * Compute single flow direction drainage area
	 */
	void compute_sfd_drainage_area(FlowRoutingResult<T>& result) const
	{
		// Process in reverse topological order
		for (int i = static_cast<int>(size_) - 1; i >= 0; --i) {
			size_t node = (*result.stack_order)[i];
			if (!connector_->is_active_node(node))
				continue;

			size_t receiver = (*result.receivers)[node];
			if (receiver != SIZE_MAX && connector_->is_active_node(receiver)) {
				(*result.drainage_area)[receiver] += (*result.drainage_area)[node];
			}
		}
	}

	/**
	 * Compute multiple flow direction drainage area
	 */
	void compute_mfd_drainage_area(FlowRoutingResult<T>& result) const
	{
		// Process in reverse topological order
		for (int i = static_cast<int>(size_) - 1; i >= 0; --i) {
			size_t node = (*result.stack_order)[i];
			if (!connector_->is_active_node(node))
				continue;

			// Distribute drainage area to all receivers based on weights
			for (const auto& [receiver, weight] : (*result.mfd_receivers)[node]) {
				if (connector_->is_active_node(receiver)) {
					(*result.drainage_area)[receiver] +=
						(*result.drainage_area)[node] * weight;
				}
			}
		}
	}
};

// ======================
// CONVENIENCE FUNCTIONS
// ======================

/**
 * Quick flow routing for common use cases
 */
template<typename T = double>
class QuickFlow
{
public:
	/**
	 * D8 flow accumulation with depression filling
	 */
	static std::shared_ptr<ArrayRef<size_t>> d8_flow_accumulation(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = FlowRoutingMethod::D8_STEEPEST;
		config.handle_depressions = true;
		config.compute_accumulation = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("D8 flow routing failed");
		}

		return result.flow_accumulation;
	}

	/**
	 * D4 flow accumulation
	 */
	static std::shared_ptr<ArrayRef<size_t>> d4_flow_accumulation(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = FlowRoutingMethod::D4_STEEPEST;
		config.handle_depressions = true;
		config.compute_accumulation = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("D4 flow routing failed");
		}

		return result.flow_accumulation;
	}

	/**
	 * Multiple flow direction (Freeman) accumulation
	 */
	static std::shared_ptr<ArrayRef<size_t>> mfd_flow_accumulation(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		T flow_exponent = 1.1)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = FlowRoutingMethod::MFD_FREEMAN;
		config.flow_exponent = flow_exponent;
		config.handle_depressions = true;
		config.compute_accumulation = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("MFD flow routing failed");
		}

		return result.flow_accumulation;
	}

	/**
	 * Drainage area computation
	 */
	static std::shared_ptr<ArrayRef<T>> drainage_area(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		T cell_area = 1.0,
		FlowRoutingMethod method = FlowRoutingMethod::D8_STEEPEST)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = method;
		config.cell_area = cell_area;
		config.handle_depressions = true;
		config.compute_drainage_area = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("Drainage area computation failed");
		}

		return result.drainage_area;
	}

	/**
	 * Flow receivers (single flow direction)
	 */
	static std::shared_ptr<ArrayRef<size_t>> flow_receivers(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		FlowRoutingMethod method = FlowRoutingMethod::D8_STEEPEST)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = method;
		config.handle_depressions = true;
		config.compute_accumulation = false;
		config.compute_drainage_area = false;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("Flow receiver computation failed");
		}

		return result.receivers;
	}

	/**
	 * Cordonnier flow routing with depression handling
	 */
	static std::shared_ptr<ArrayRef<size_t>> cordonnier_flow(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		bool use_filling = true)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = use_filling ? FlowRoutingMethod::CORDONNIER_FILLING
																: FlowRoutingMethod::CORDONNIER_CARVING;
		config.handle_depressions =
			false; // Cordonnier handles depressions internally
		config.compute_accumulation = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("Cordonnier flow routing failed");
		}

		return result.flow_accumulation;
	}

	// ======================
	// 2D GRID VERSIONS
	// ======================

	/**
	 * D8 flow accumulation with depression filling (Grid2D version)
	 */
	static std::shared_ptr<Grid2D<size_t>> d8_flow_accumulation_2d(
		std::shared_ptr<Connector<T>> connector,
		const Grid2D<T>& elevation)
	{
		auto elevation_ref = std::make_shared<Grid2D<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = FlowRoutingMethod::D8_STEEPEST;
		config.handle_depressions = true;
		config.compute_accumulation = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("D8 flow routing failed");
		}

		return result.flow_accumulation_2d;
	}

	/**
	 * Flow receivers (Grid2D version)
	 */
	static std::shared_ptr<Grid2D<size_t>> flow_receivers_2d(
		std::shared_ptr<Connector<T>> connector,
		const Grid2D<T>& elevation,
		FlowRoutingMethod method = FlowRoutingMethod::D8_STEEPEST)
	{
		auto elevation_ref = std::make_shared<Grid2D<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = method;
		config.handle_depressions = true;
		config.compute_accumulation = false;
		config.compute_drainage_area = false;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("Flow receiver computation failed");
		}

		return result.receivers_2d;
	}

	/**
	 * Drainage area computation (Grid2D version)
	 */
	static std::shared_ptr<Grid2D<T>> drainage_area_2d(
		std::shared_ptr<Connector<T>> connector,
		const Grid2D<T>& elevation,
		T cell_area = 1.0,
		FlowRoutingMethod method = FlowRoutingMethod::D8_STEEPEST)
	{
		auto elevation_ref = std::make_shared<Grid2D<T>>(elevation);

		FlowRoutingConfig<T> config;
		config.method = method;
		config.cell_area = cell_area;
		config.handle_depressions = true;
		config.compute_drainage_area = true;

		FlowRouter<T> router(connector, elevation_ref, config);
		auto result = router.compute_flow_routing();

		if (!result.success) {
			throw std::runtime_error("Drainage area computation failed");
		}

		return result.drainage_area_2d;
	}
};

// ======================
// TYPE ALIASES
// ======================

// Type aliases for common usage
using FlowRouterF32 = FlowRouter<float>;
using FlowRouterF64 = FlowRouter<double>;
using FlowRoutingResultF32 = FlowRoutingResult<float>;
using FlowRoutingResultF64 = FlowRoutingResult<double>;
using FlowRoutingConfigF32 = FlowRoutingConfig<float>;
using FlowRoutingConfigF64 = FlowRoutingConfig<double>;

} // namespace dagger2
