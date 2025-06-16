#pragma once

#include "dg2_BCs.hpp"
#include "dg2_BCs_helper.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include "dg2_cordonnier_method.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * Depression filling method enumeration
 */
enum class FillingMethod : uint8_t
{
	// Priority flood variants
	PRIORITY_FLOOD_EPSILON = 0, // Barnes et al. 2014 - Priority Flood + epsilon
	PRIORITY_FLOOD_OPTIMAL = 1, // Wang & Liu 2006 - Priority Flood optimal
	PRIORITY_FLOOD_STRICT = 2,	// Strict priority flood (no epsilon)

	// Hybrid methods (Cordonnier et al. 2019)
	HYBRID_PRIORITY_FLOODING = 3, // Combination of different strategies
	SLOPED_FACETS = 4,						// Zhou et al. 2016 - sloped facets
	IMPACT_REDUCTION = 5,					// Minimize impact on original DEM

	// Carving methods
	LEAST_COST_PATH = 6,		// Carve channels using least cost paths
	BREACH_DEPRESSIONS = 7, // Breach depression outlets
	HYBRID_FILL_BREACH = 8, // Combined fill and breach approach

	// Advanced methods
	FLOW_BREACHING = 9,				 // Lindsay 2016 - flow accumulation based
	CONSTRAINED_FLOW = 10,		 // Constrained to minimize modification
	MORPHOLOGICAL_BREACH = 11, // Morphology-aware breaching

	// Specialized methods
	DEPRESSION_HIERARCHY = 12,		// Hierarchical depression treatment
	HYDROLOGICAL_CORRECTION = 13, // Preserve hydrological features
	TERRAIN_CARVING = 14,					// Advanced terrain modification
	STOCHASTIC_DEPRESSION = 15,		// Stochastic approach for uncertainty

	// Cordonnier et al. 2019 methods
	CORDONNIER_SIMPLE = 16,	 // Simple correction (minimal receiver updates)
	CORDONNIER_CARVING = 17, // Depression carving via basin graph
	CORDONNIER_FILLING = 18	 // Depression filling via basin graph
};

/**
 * Depression filling configuration
 */
template<typename T>
struct FillingConfig
{
	FillingMethod method;
	T epsilon;								// Small increment for priority flood
	T breach_threshold;				// Maximum depth to breach vs fill
	T carving_cost_threshold; // Cost threshold for path carving
	bool preserve_ridges;			// Avoid modifying ridge lines
	bool preserve_valleys;		// Maintain valley bottoms
	bool enforce_drainage;		// Ensure all cells drain
	size_t max_iterations;		// Maximum iterations for iterative methods
	T convergence_tolerance;	// Convergence threshold
	uint32_t random_seed;			// For stochastic methods
	bool parallel_processing; // Enable parallel computation
	bool verbose;							// Output progress information

	// Advanced parameters
	T slope_threshold;					// Minimum slope to maintain
	T flow_accumulation_weight; // Weight for flow-based decisions
	T morphological_weight;			// Weight for morphological preservation
	T min_gradient;
	bool use_flow_barriers; // Respect flow barriers in BCs
	bool adaptive_epsilon;	// Adjust epsilon based on local conditions

	FillingConfig()
		: method(FillingMethod::PRIORITY_FLOOD_EPSILON)
		, epsilon(1e-4)
		, breach_threshold(5.0)
		, carving_cost_threshold(10.0)
		, preserve_ridges(false)
		, preserve_valleys(true)
		, enforce_drainage(true)
		, max_iterations(1000)
		, convergence_tolerance(1e-6)
		, random_seed(42)
		, parallel_processing(true)
		, verbose(false)
		, slope_threshold(1e-6)
		, flow_accumulation_weight(1.0)
		, morphological_weight(0.5)
		, min_gradient(static_cast<T>(1e-6))
		, use_flow_barriers(true)
		, adaptive_epsilon(false)

	{
	}
};

/**
 * Depression filling results and statistics
 */
template<typename T>
struct FillingResult
{
	bool success;							 // Whether filling completed successfully
	size_t iterations_used;		 // Number of iterations performed
	T total_volume_filled;		 // Total volume of depressions filled
	T total_volume_carved;		 // Total volume carved/breached
	T max_modification;				 // Maximum elevation change
	T mean_modification;			 // Mean elevation change
	size_t cells_modified;		 // Number of cells modified
	size_t depressions_filled; // Number of depressions identified and filled
	size_t outlets_created;		 // Number of artificial outlets created
	std::chrono::milliseconds processing_time; // Processing duration
	std::vector<std::pair<size_t, T>>
		major_modifications; // Large modifications (index, change)

	FillingResult()
		: success(false)
		, iterations_used(0)
		, total_volume_filled(0)
		, total_volume_carved(0)
		, max_modification(0)
		, mean_modification(0)
		, cells_modified(0)
		, depressions_filled(0)
		, outlets_created(0)
		, processing_time(0)
	{
	}
};

/**
 * Depression information structure
 */
template<typename T>
struct Depression
{
	size_t id;											 // Unique depression identifier
	size_t outlet_index;						 // Index of depression outlet
	T outlet_elevation;							 // Elevation of outlet
	T min_elevation;								 // Minimum elevation in depression
	T volume;												 // Volume of depression
	size_t area;										 // Number of cells in depression
	std::vector<size_t> cells;			 // All cells in this depression
	T fill_level;										 // Level to fill to
	bool should_breach;							 // Whether to breach instead of fill
	std::vector<size_t> breach_path; // Path to carve if breaching

	Depression()
		: id(SIZE_MAX)
		, outlet_index(SIZE_MAX)
		, outlet_elevation(0)
		, min_elevation(0)
		, volume(0)
		, area(0)
		, fill_level(0)
		, should_breach(false)
	{
	}
};

/**
 * Comprehensive Depression Filling Engine
 *
 * This class implements multiple state-of-the-art depression filling algorithms
 * with full support for boundary conditions and the DAGGER2 framework.
 * Based on Barnes et al. 2014, Cordonnier et al. 2019, and other recent
 * advances.
 */
template<typename T = double>
class DepressionFiller
{
private:
	// Core components
	std::shared_ptr<Connector<T>> connector_;
	std::shared_ptr<ArrayRef<T>> elevation_;
	FillingConfig<T> config_;

	// Grid properties
	size_t rows_;
	size_t cols_;
	size_t size_;

	// Working arrays
	mutable std::vector<T> working_elevation_;
	mutable std::vector<bool> processed_;
	mutable std::vector<bool> in_queue_;
	mutable std::vector<size_t> labels_;
	mutable std::vector<T> flow_accumulation_;

	// Cordonnier method instance
	mutable std::unique_ptr<CordonnierFlowRouter<T>> cordonnier_router_;

	// Depression tracking
	mutable std::vector<Depression<T>> depressions_;
	mutable std::unordered_map<size_t, size_t> cell_to_depression_;

	// Random number generator
	mutable std::mt19937 rng_;

	// Priority queue element for flood algorithms
	struct QueueElement
	{
		size_t index;
		T elevation;
		size_t priority_level; // For tie-breaking

		QueueElement(size_t idx, T elev, size_t level = 0)
			: index(idx)
			, elevation(elev)
			, priority_level(level)
		{
		}

		bool operator<(const QueueElement& other) const
		{
			if (std::abs(elevation - other.elevation) <
					std::numeric_limits<T>::epsilon()) {
				return priority_level >
							 other.priority_level; // Lower priority level = higher priority
			}
			return elevation > other.elevation; // Lower elevation = higher priority
		}
	};

public:
	/**
	 * Constructor
	 */
	DepressionFiller(std::shared_ptr<Connector<T>> connector,
									 std::shared_ptr<ArrayRef<T>> elevation,
									 const FillingConfig<T>& config = FillingConfig<T>())
		: connector_(connector)
		, elevation_(elevation)
		, config_(config)
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

		// Initialize working arrays
		working_elevation_.resize(size_);
		processed_.resize(size_);
		in_queue_.resize(size_);
		labels_.resize(size_);
		flow_accumulation_.resize(size_);

		// Copy initial elevation data
		for (size_t i = 0; i < size_; ++i) {
			working_elevation_[i] = (*elevation_)[i];
		}
	}

	// ======================
	// CONFIGURATION
	// ======================

	const FillingConfig<T>& get_config() const { return config_; }

	void set_config(const FillingConfig<T>& config)
	{
		config_ = config;
		rng_.seed(config_.random_seed);
	}

	void set_method(FillingMethod method) { config_.method = method; }
	void set_epsilon(T epsilon) { config_.epsilon = epsilon; }
	void set_breach_threshold(T threshold)
	{
		config_.breach_threshold = threshold;
	}
	void set_preserve_ridges(bool preserve)
	{
		config_.preserve_ridges = preserve;
	}
	void set_preserve_valleys(bool preserve)
	{
		config_.preserve_valleys = preserve;
	}

	// ======================
	// MAIN FILLING INTERFACE
	// ======================

	/**
	 * Fill depressions using configured method
	 */
	FillingResult<T> fill_depressions() const
	{
		auto start_time = std::chrono::high_resolution_clock::now();

		FillingResult<T> result;

		try {
			switch (config_.method) {
				case FillingMethod::PRIORITY_FLOOD_EPSILON:
					result = priority_flood_epsilon();
					break;
				case FillingMethod::PRIORITY_FLOOD_OPTIMAL:
					result = priority_flood_optimal();
					break;
				case FillingMethod::PRIORITY_FLOOD_STRICT:
					result = priority_flood_strict();
					break;
				case FillingMethod::HYBRID_PRIORITY_FLOODING:
					result = hybrid_priority_flooding();
					break;
				case FillingMethod::SLOPED_FACETS:
					result = sloped_facets_method();
					break;
				case FillingMethod::IMPACT_REDUCTION:
					result = impact_reduction_method();
					break;
				case FillingMethod::LEAST_COST_PATH:
					result = least_cost_path_carving();
					break;
				case FillingMethod::BREACH_DEPRESSIONS:
					result = breach_depressions_method();
					break;
				case FillingMethod::HYBRID_FILL_BREACH:
					result = hybrid_fill_breach();
					break;
				case FillingMethod::FLOW_BREACHING:
					result = flow_breaching_method();
					break;
				case FillingMethod::CONSTRAINED_FLOW:
					result = constrained_flow_method();
					break;
				case FillingMethod::MORPHOLOGICAL_BREACH:
					result = morphological_breach_method();
					break;
				case FillingMethod::DEPRESSION_HIERARCHY:
					result = depression_hierarchy_method();
					break;
				case FillingMethod::HYDROLOGICAL_CORRECTION:
					result = hydrological_correction_method();
					break;
				case FillingMethod::TERRAIN_CARVING:
					result = terrain_carving_method();
					break;
				case FillingMethod::STOCHASTIC_DEPRESSION:
					result = stochastic_depression_method();
					break;
				case FillingMethod::CORDONNIER_SIMPLE:
					result = cordonnier_simple_method();
					break;
				case FillingMethod::CORDONNIER_CARVING:
					result = cordonnier_carving_method();
					break;
				case FillingMethod::CORDONNIER_FILLING:
					result = cordonnier_filling_method();
					break;
				default:
					throw std::invalid_argument("Unknown filling method");
			}

			result.success = true;

		} catch (const std::exception& e) {
			result.success = false;
			if (config_.verbose) {
				std::cerr << "Depression filling failed: " << e.what() << std::endl;
			}
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		result.processing_time =
			std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
																														start_time);

		return result;
	}

	/**
	 * Get filled elevation data
	 */
	std::vector<T> get_filled_elevation() const { return working_elevation_; }

	/**
	 * Get filled elevation as ArrayRef
	 */
	ArrayRef<T> get_filled_elevation_ref() const
	{
		return ArrayRef<T>(working_elevation_);
	}

	/**
	 * Apply results to original elevation array
	 */
	void apply_to_elevation() const
	{
		for (size_t i = 0; i < size_; ++i) {
			(*elevation_)[i] = working_elevation_[i];
		}
	}

	/**
	 * Get depression information
	 */
	const std::vector<Depression<T>>& get_depressions() const
	{
		return depressions_;
	}

	// ======================
	// PRIORITY FLOOD METHODS
	// ======================

private:
	/**
	 * Priority Flood + Epsilon (Barnes et al. 2014)
	 * The gold standard for depression filling
	 */
	FillingResult<T> priority_flood_epsilon() const
	{
		FillingResult<T> result;

		// Initialize working arrays
		reset_working_arrays();

		std::priority_queue<QueueElement> open_set;

		// Add boundary cells to open set
		initialize_boundary_queue(open_set);

		T current_epsilon = config_.epsilon;
		size_t cells_processed = 0;

		while (!open_set.empty()) {
			QueueElement current = open_set.top();
			open_set.pop();

			if (processed_[current.index])
				continue;

			processed_[current.index] = true;
			cells_processed++;

			// Adaptive epsilon adjustment
			if (config_.adaptive_epsilon) {
				current_epsilon = calculate_adaptive_epsilon(current.index);
			}

			// Process neighbors
			auto neighbors = connector_->get_effective_valid_neighbors(current.index);
			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || processed_[neighbor.index])
					continue;

				T original_elevation = (*elevation_)[neighbor.index];
				T new_elevation =
					std::max(original_elevation,
									 working_elevation_[current.index] + current_epsilon);

				// Check boundary conditions
				if (!should_modify_cell(neighbor.index, new_elevation))
					continue;

				if (new_elevation > working_elevation_[neighbor.index]) {
					working_elevation_[neighbor.index] = new_elevation;

					if (new_elevation > original_elevation) {
						result.total_volume_filled += (new_elevation - original_elevation);
						result.cells_modified++;
						result.max_modification = std::max(
							result.max_modification, new_elevation - original_elevation);
					}
				}

				if (!in_queue_[neighbor.index]) {
					open_set.emplace(neighbor.index,
													 working_elevation_[neighbor.index],
													 cells_processed);
					in_queue_[neighbor.index] = true;
				}
			}
		}

		result.iterations_used = cells_processed;
		result.mean_modification =
			result.cells_modified > 0
				? result.total_volume_filled / result.cells_modified
				: 0;

		return result;
	}

	/**
	 * Optimal Priority Flood (Wang & Liu 2006)
	 * More conservative approach with better preservation of original terrain
	 */
	FillingResult<T> priority_flood_optimal() const
	{
		FillingResult<T> result;

		reset_working_arrays();

		// First pass: identify all depressions
		identify_depressions();

		// Second pass: fill each depression optimally
		for (const auto& depression : depressions_) {
			if (depression.volume < config_.breach_threshold) {
				// Small depressions: fill completely
				fill_depression_optimal(depression, result);
			} else {
				// Large depressions: consider breaching
				if (should_breach_depression(depression)) {
					breach_depression(depression, result);
				} else {
					fill_depression_optimal(depression, result);
				}
			}
		}

		return result;
	}

	/**
	 * Strict Priority Flood (no epsilon)
	 * Ensures exact drainage without artificial gradients
	 */
	FillingResult<T> priority_flood_strict() const
	{
		FillingResult<T> result;

		reset_working_arrays();

		std::priority_queue<QueueElement> open_set;
		initialize_boundary_queue(open_set);

		while (!open_set.empty()) {
			QueueElement current = open_set.top();
			open_set.pop();

			if (processed_[current.index])
				continue;
			processed_[current.index] = true;

			auto neighbors = connector_->get_effective_valid_neighbors(current.index);
			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || processed_[neighbor.index])
					continue;

				T original_elevation = (*elevation_)[neighbor.index];
				T new_elevation =
					std::max(original_elevation, working_elevation_[current.index]);

				if (!should_modify_cell(neighbor.index, new_elevation))
					continue;

				working_elevation_[neighbor.index] = new_elevation;

				if (new_elevation > original_elevation) {
					result.total_volume_filled += (new_elevation - original_elevation);
					result.cells_modified++;
					result.max_modification = std::max(
						result.max_modification, new_elevation - original_elevation);
				}

				if (!in_queue_[neighbor.index]) {
					open_set.emplace(neighbor.index, working_elevation_[neighbor.index]);
					in_queue_[neighbor.index] = true;
				}
			}
		}

		return result;
	}

	// ======================
	// HYBRID METHODS (CORDONNIER ET AL. 2019)
	// ======================

	/**
	 * Hybrid Priority Flooding
	 * Combines multiple strategies based on depression characteristics
	 */
	FillingResult<T> hybrid_priority_flooding() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();
		compute_flow_accumulation();

		for (const auto& depression : depressions_) {
			// Decision matrix for treatment method
			bool is_large = depression.volume > config_.breach_threshold;
			bool has_flow = has_significant_flow(depression);
			bool is_deep = (depression.fill_level - depression.min_elevation) >
										 config_.breach_threshold;
			bool preserves_morphology =
				config_.preserve_valleys && is_valley_depression(depression);

			if (is_large && has_flow && !preserves_morphology) {
				// Breach large, flowing depressions
				breach_depression(depression, result);
			} else if (is_deep && !preserves_morphology) {
				// Carve deep depressions
				carve_depression_path(depression, result);
			} else {
				// Fill smaller or morphologically important depressions
				if (config_.preserve_ridges && affects_ridges(depression)) {
					fill_depression_minimal_impact(depression, result);
				} else {
					fill_depression_optimal(depression, result);
				}
			}
		}

		return result;
	}

	/**
	 * Carve depression path - missing function
	 */
	void carve_depression_path(const Depression<T>& depression,
														 FillingResult<T>& result) const
	{
		if (depression.breach_path.empty()) {
			return;
		}

		// Calculate target elevations along path
		T start_elev = depression.min_elevation;
		T end_elev = working_elevation_[depression.outlet_index];

		if (depression.breach_path.size() < 2) {
			return;
		}

		// Ensure minimum gradient
		T total_distance = static_cast<T>(depression.breach_path.size() - 1);
		T min_drop = config_.min_gradient * total_distance;
		T actual_drop = start_elev - end_elev;

		if (actual_drop < min_drop) {
			end_elev = start_elev - min_drop;
		}

		// Carve along the path
		for (size_t i = 0; i < depression.breach_path.size(); ++i) {
			size_t cell_idx = depression.breach_path[i];

			// Calculate target elevation with gradient
			T progress = static_cast<T>(i) / (depression.breach_path.size() - 1);
			T target_elev = start_elev - progress * (start_elev - end_elev);

			if (target_elev < working_elevation_[cell_idx]) {
				T volume_carved = working_elevation_[cell_idx] - target_elev;
				working_elevation_[cell_idx] = target_elev;

				result.total_volume_carved += volume_carved;
				result.cells_modified++;
				result.max_modification =
					std::max(result.max_modification, volume_carved);
			}
		}
	}

	/**
	 * Sloped Facets Method (Zhou et al. 2016)
	 * Creates gradual slopes instead of flat areas
	 */
	FillingResult<T> sloped_facets_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		for (const auto& depression : depressions_) {
			create_sloped_facets(depression, result);
		}

		return result;
	}

	/**
	 * Impact Reduction Method
	 * Minimizes modification to original DEM
	 */
	FillingResult<T> impact_reduction_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		// Sort depressions by impact potential
		auto sorted_depressions = depressions_;
		std::sort(sorted_depressions.begin(),
							sorted_depressions.end(),
							[](const Depression<T>& a, const Depression<T>& b) {
								return a.volume < b.volume; // Fill smallest first
							});

		for (const auto& depression : sorted_depressions) {
			// Try multiple approaches and choose the one with least impact
			std::vector<T> fill_impact = estimate_fill_impact(depression);
			std::vector<T> breach_impact = estimate_breach_impact(depression);
			std::vector<T> carve_impact = estimate_carve_impact(depression);

			// Choose method with minimum impact
			T min_fill = *std::min_element(fill_impact.begin(), fill_impact.end());
			T min_breach =
				*std::min_element(breach_impact.begin(), breach_impact.end());
			T min_carve = *std::min_element(carve_impact.begin(), carve_impact.end());

			if (min_fill <= min_breach && min_fill <= min_carve) {
				fill_depression_minimal_impact(depression, result);
			} else if (min_breach <= min_carve) {
				breach_depression(depression, result);
			} else {
				carve_depression_path(depression, result);
			}
		}

		return result;
	}

	// ======================
	// CARVING AND BREACHING METHODS
	// ======================

	/**
	 * Least Cost Path Carving
	 * Creates drainage channels using optimal path algorithms
	 */
	FillingResult<T> least_cost_path_carving() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		for (const auto& depression : depressions_) {
			// Find least cost path from depression to nearest outlet
			auto path = find_least_cost_path(depression);
			if (!path.empty()) {
				carve_path(path, depression, result);
			} else {
				// Fall back to filling if no path found
				fill_depression_optimal(depression, result);
			}
		}

		return result;
	}

	/**
	 * Breach Depressions Method
	 * Simple breaching approach for all depressions
	 */
	FillingResult<T> breach_depressions_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		for (const auto& depression : depressions_) {
			breach_depression(depression, result);
		}

		return result;
	}

	/**
	 * Hybrid Fill-Breach Method
	 * Intelligently combines filling and breaching
	 */
	FillingResult<T> hybrid_fill_breach() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();
		compute_flow_accumulation();

		for (const auto& depression : depressions_) {
			T fill_cost = estimate_fill_cost(depression);
			T breach_cost = estimate_breach_cost(depression);

			// Include flow accumulation in decision
			T flow_weight = get_depression_flow_weight(depression);
			breach_cost *= (1.0 + flow_weight * config_.flow_accumulation_weight);

			if (breach_cost < fill_cost) {
				breach_depression(depression, result);
			} else {
				fill_depression_optimal(depression, result);
			}
		}

		return result;
	}

	// ======================
	// ADVANCED METHODS
	// ======================

	/**
	 * Flow Breaching Method (Lindsay 2016)
	 * Uses flow accumulation to guide breaching decisions
	 */
	FillingResult<T> flow_breaching_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		compute_flow_accumulation();
		identify_depressions();

		for (const auto& depression : depressions_) {
			if (has_significant_flow(depression)) {
				// Create channel following flow accumulation patterns
				auto flow_path = trace_flow_path(depression);
				if (!flow_path.empty()) {
					carve_flow_channel(flow_path, depression, result);
					continue;
				}
			}

			// Fall back to standard breach or fill
			if (depression.volume > config_.breach_threshold) {
				breach_depression(depression, result);
			} else {
				fill_depression_optimal(depression, result);
			}
		}

		return result;
	}

	/**
	 * Constrained Flow Method
	 * Minimizes modifications while ensuring drainage
	 */
	FillingResult<T> constrained_flow_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		// Create constraint map
		std::vector<bool> constraints = create_constraint_map();

		for (const auto& depression : depressions_) {
			if (violates_constraints(depression, constraints)) {
				// Use minimal modification approach
				fill_depression_constrained(depression, constraints, result);
			} else {
				// Standard approach
				fill_depression_optimal(depression, result);
			}
		}

		return result;
	}

	/**
	 * Morphological Breach Method
	 * Preserves important morphological features
	 */
	FillingResult<T> morphological_breach_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		// Identify morphological features
		auto morphological_map = identify_morphological_features();

		for (const auto& depression : depressions_) {
			if (affects_morphological_features(depression, morphological_map)) {
				// Use morphology-preserving approach
				breach_depression_morphology_aware(
					depression, morphological_map, result);
			} else {
				// Standard breach
				breach_depression(depression, result);
			}
		}

		return result;
	}

	/**
	 * Depression Hierarchy Method
	 * Treats depressions in hierarchical order
	 */
	FillingResult<T> depression_hierarchy_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		// Build depression hierarchy
		auto hierarchy = build_depression_hierarchy();

		// Process from most important to least important
		for (const auto& level : hierarchy) {
			for (size_t depression_id : level) {
				const auto& depression = depressions_[depression_id];

				// Treatment depends on hierarchy level
				if (is_primary_depression(depression, hierarchy)) {
					fill_depression_optimal(depression, result);
				} else {
					// Secondary depressions: try breaching first
					if (can_breach_safely(depression)) {
						breach_depression(depression, result);
					} else {
						fill_depression_optimal(depression, result);
					}
				}
			}
		}

		return result;
	}

	/**
	 * Hydrological Correction Method
	 * Maintains hydrological consistency
	 */
	FillingResult<T> hydrological_correction_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();
		compute_flow_accumulation();

		// Identify hydrologically important areas
		auto hydro_map = identify_hydrological_features();

		for (const auto& depression : depressions_) {
			if (is_hydrologically_significant(depression, hydro_map)) {
				// Preserve hydrological function
				create_hydrological_channel(depression, hydro_map, result);
			} else {
				// Standard treatment
				fill_depression_optimal(depression, result);
			}
		}

		return result;
	}

	/**
	 * Terrain Carving Method
	 * Advanced terrain modification with multiple criteria
	 */
	FillingResult<T> terrain_carving_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();
		compute_flow_accumulation();

		// Multi-criteria evaluation for each depression
		for (const auto& depression : depressions_) {
			T terrain_score = evaluate_terrain_modification_score(depression);

			if (terrain_score > config_.carving_cost_threshold) {
				// High score: use sophisticated carving
				carve_terrain_sophisticated(depression, result);
			} else {
				// Low score: simple treatment
				fill_depression_optimal(depression, result);
			}
		}

		return result;
	}

	/**
	 * Stochastic Depression Method
	 * Handles uncertainty in depression treatment
	 */
	FillingResult<T> stochastic_depression_method() const
	{
		FillingResult<T> result;

		reset_working_arrays();
		identify_depressions();

		std::uniform_real_distribution<T> random_dist(0.0, 1.0);

		for (const auto& depression : depressions_) {
			// Stochastic decision making
			T fill_probability = calculate_fill_probability(depression);
			T random_value = random_dist(rng_);

			if (random_value < fill_probability) {
				fill_depression_optimal(depression, result);
			} else {
				// Breach with some randomness in path selection
				breach_depression_stochastic(depression, result);
			}
		}

		return result;
	}

	// ======================
	// UTILITY FUNCTIONS
	// ======================

private:
	/**
	 * Reset all working arrays to initial state
	 */
	void reset_working_arrays() const
	{
		std::fill(processed_.begin(), processed_.end(), false);
		std::fill(in_queue_.begin(), in_queue_.end(), false);
		std::fill(labels_.begin(), labels_.end(), 0);

		for (size_t i = 0; i < size_; ++i) {
			working_elevation_[i] = (*elevation_)[i];
		}

		depressions_.clear();
		cell_to_depression_.clear();
	}

	/**
	 * Initialize boundary queue for priority flood algorithms
	 */
	void initialize_boundary_queue(std::priority_queue<QueueElement>& queue) const
	{
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			NodeType boundary_type = connector_->get_boundary_type(i);

			// Add boundary nodes and nodes adjacent to boundaries
			if (NodeTypeUtils::allows_outflow(boundary_type) ||
					NodeTypeUtils::is_boundary(boundary_type)) {
				queue.emplace(i, working_elevation_[i], 0);
				in_queue_[i] = true;
			} else {
				// Check if adjacent to boundary
				auto neighbors = connector_->get_all_neighbors(i);
				for (const auto& neighbor : neighbors) {
					if (neighbor.is_valid) {
						NodeType neighbor_type =
							connector_->get_boundary_type(neighbor.index);
						if (NodeTypeUtils::allows_outflow(neighbor_type)) {
							queue.emplace(i, working_elevation_[i], 1);
							in_queue_[i] = true;
							break;
						}
					} else if (config_.use_flow_barriers) {
						// Adjacent to grid boundary
						queue.emplace(i, working_elevation_[i], 1);
						in_queue_[i] = true;
						break;
					}
				}
			}
		}
	}

	/**
	 * Check if a cell should be modified based on boundary conditions
	 */
	bool should_modify_cell(size_t index, T new_elevation) const
	{
		NodeType boundary_type = connector_->get_boundary_type(index);

		// Never modify NO_DATA cells
		if (boundary_type == NodeType::NO_DATA)
			return false;

		// Respect flow barriers if configured
		if (config_.use_flow_barriers) {
			if (boundary_type == NodeType::REFLECT)
				return false;
			if (boundary_type == NodeType::IN && new_elevation > (*elevation_)[index])
				return false;
		}

		// Preserve ridge lines if configured
		if (config_.preserve_ridges && is_ridge_cell(index))
			return false;

		return true;
	}

	/**
	 * Calculate adaptive epsilon based on local conditions
	 */
	T calculate_adaptive_epsilon(size_t index) const
	{
		if (!config_.adaptive_epsilon)
			return config_.epsilon;

		// Base epsilon on local slope and neighborhood characteristics
		auto neighbors = connector_->get_valid_neighbors(index);
		if (neighbors.empty())
			return config_.epsilon;

		T local_slope = 0;
		T elevation_variance = 0;
		T center_elev = working_elevation_[index];

		for (const auto& neighbor : neighbors) {
			T neighbor_elev = working_elevation_[neighbor.index];
			T gradient = std::abs(center_elev - neighbor_elev) / neighbor.distance;
			local_slope += gradient;
			elevation_variance +=
				(neighbor_elev - center_elev) * (neighbor_elev - center_elev);
		}

		local_slope /= neighbors.size();
		elevation_variance /= neighbors.size();

		// Adaptive epsilon: larger in flat areas, smaller in steep areas
		T adaptive_factor =
			1.0 + std::sqrt(elevation_variance) / (local_slope + config_.epsilon);
		return config_.epsilon * std::min(adaptive_factor, static_cast<T>(10.0));
	}

	/**
	 * Identify all depressions in the current elevation data
	 */
	// void identify_depressions() const
	// {
	// 	depressions_.clear();
	// 	cell_to_depression_.clear();

	// 	std::vector<bool> visited(size_, false);
	// 	size_t depression_id = 0;

	// 	for (size_t i = 0; i < size_; ++i) {
	// 		if (!connector_->is_active_node(i) || visited[i])
	// 			continue;

	// 		// Check if this is a local minimum
	// 		if (is_local_minimum(i)) {
	// 			Depression<T> depression;
	// 			depression.id = depression_id++;
	// 			depression.min_elevation = working_elevation_[i];

	// 			// Flood fill to find entire depression
	// 			flood_fill_depression(i, depression, visited);

	// 			if (depression.area > 0) {
	// 				depressions_.push_back(depression);
	// 			}
	// 		}
	// 	}

	// 	// Calculate depression characteristics
	// 	for (auto& depression : depressions_) {
	// 		calculate_depression_properties(depression);
	// 	}
	// }

	/**
	 * Improved depression identification using priority queue approach
	 * This method properly identifies depressions and their spill elevations
	 */
	void identify_depressions() const
	{
		depressions_.clear();
		cell_to_depression_.clear();

		// Arrays for the identification process
		std::vector<bool> visited(size_, false);
		std::vector<bool> is_boundary_or_processed(size_, false);
		std::vector<T> spill_elevation(size_, std::numeric_limits<T>::max());
		std::vector<size_t> depression_id(size_, SIZE_MAX);

		// Priority queue for processing cells from low to high elevation
		std::priority_queue<QueueElement> pq;

		// Step 1: Initialize with boundary cells and cells that can drain
		initialize_boundary_queue(pq);

		// Mark boundary cells as processed
		while (!pq.empty()) {
			QueueElement current = pq.top();
			pq.pop();

			if (is_boundary_or_processed[current.index])
				continue;

			is_boundary_or_processed[current.index] = true;
			spill_elevation[current.index] = working_elevation_[current.index];

			// Add unprocessed neighbors to queue
			auto neighbors = connector_->get_valid_neighbors(current.index);
			for (const auto& neighbor : neighbors) {
				if (!is_boundary_or_processed[neighbor.index]) {
					T neighbor_elev = working_elevation_[neighbor.index];
					// Can drain to current cell if neighbor is higher or equal
					if (neighbor_elev >= working_elevation_[current.index]) {
						pq.emplace(neighbor.index, neighbor_elev);
					}
				}
			}
		}

		// Step 2: Identify remaining cells as potential depression cells
		std::vector<size_t> depression_candidates;
		for (size_t i = 0; i < size_; ++i) {
			if (connector_->is_active_node(i) && !is_boundary_or_processed[i]) {
				depression_candidates.push_back(i);
			}
		}

		if (depression_candidates.empty()) {
			return; // No depressions found
		}

		// Step 3: Process depression candidates to find spill elevations
		std::fill(visited.begin(), visited.end(), false);
		size_t current_depression_id = 0;

		for (size_t candidate_idx : depression_candidates) {
			if (visited[candidate_idx])
				continue;

			// Start a new depression from this candidate
			Depression<T> depression;
			depression.id = current_depression_id++;
			depression.min_elevation = std::numeric_limits<T>::max();
			depression.outlet_elevation = std::numeric_limits<T>::max();
			depression.outlet_index = SIZE_MAX;

			// Use flood-fill to find all cells in this depression
			std::queue<size_t> flood_queue;
			flood_queue.push(candidate_idx);
			visited[candidate_idx] = true;

			// Track potential outlets (boundary between depression and drained areas)
			std::vector<std::pair<size_t, T>> potential_outlets;

			while (!flood_queue.empty()) {
				size_t current = flood_queue.front();
				flood_queue.pop();

				T current_elev = working_elevation_[current];
				depression.cells.push_back(current);
				depression.area++;
				depression.min_elevation =
					std::min(depression.min_elevation, current_elev);
				cell_to_depression_[current] = depression.id;
				depression_id[current] = depression.id;

				auto neighbors = connector_->get_valid_neighbors(current);
				for (const auto& neighbor : neighbors) {
					T neighbor_elev = working_elevation_[neighbor.index];

					if (is_boundary_or_processed[neighbor.index]) {
						// This neighbor can drain - it's a potential outlet
						T spill_elev = std::max(current_elev, neighbor_elev);
						potential_outlets.emplace_back(neighbor.index, spill_elev);
					} else if (!visited[neighbor.index]) {
						// Check if neighbor should be part of same depression
						if (can_flow_between_cells(current, neighbor.index)) {
							visited[neighbor.index] = true;
							flood_queue.push(neighbor.index);
						} else {
							// Different depression or higher area
							T spill_elev = std::max(current_elev, neighbor_elev);
							potential_outlets.emplace_back(neighbor.index, spill_elev);
						}
					}
				}
			}

			// Step 4: Find the true outlet (lowest spill elevation)
			if (!potential_outlets.empty()) {
				auto min_outlet = *std::min_element(
					potential_outlets.begin(),
					potential_outlets.end(),
					[](const auto& a, const auto& b) { return a.second < b.second; });
				depression.outlet_index = min_outlet.first;
				depression.outlet_elevation = min_outlet.second;
				depression.fill_level = depression.outlet_elevation;
			} else {
				// No outlet found - this shouldn't happen if boundary initialization
				// worked
				depression.outlet_elevation =
					depression.min_elevation + config_.epsilon;
				depression.fill_level = depression.outlet_elevation;
			}

			// Only add non-trivial depressions
			if (depression.area > 0 &&
					depression.fill_level > depression.min_elevation + config_.epsilon) {
				depressions_.push_back(depression);
			}
		}

		// Step 5: Calculate final depression properties
		for (auto& depression : depressions_) {
			calculate_depression_properties_corrected(depression);
		}
	}

	/**
	 * Helper function to determine if flow can occur between two adjacent cells
	 */
	bool can_flow_between_cells(size_t from_idx, size_t to_idx) const
	{
		T from_elev = working_elevation_[from_idx];
		T to_elev = working_elevation_[to_idx];

		// Cells are in same depression if:
		// 1. Elevations are very similar (within epsilon)
		// 2. OR there's no significant barrier between them
		T elevation_diff = std::abs(from_elev - to_elev);

		if (elevation_diff <= config_.epsilon) {
			return true;
		}

		// Check for monotonic flow path
		if (from_elev <= to_elev) {
			return true; // Water can flow from lower to equal/higher with epsilon
		}

		return false;
	}

	/**
	 * Corrected depression properties calculation
	 */
	void calculate_depression_properties_corrected(
		Depression<T>& depression) const
	{
		depression.volume = 0;
		depression.min_elevation = std::numeric_limits<T>::max();

		// Recalculate min elevation and volume
		for (size_t cell_idx : depression.cells) {
			T cell_elev = working_elevation_[cell_idx];
			depression.min_elevation = std::min(depression.min_elevation, cell_elev);
		}

		// Calculate volume based on correct fill level
		for (size_t cell_idx : depression.cells) {
			T cell_elev = working_elevation_[cell_idx];
			if (depression.fill_level > cell_elev) {
				depression.volume += (depression.fill_level - cell_elev);
			}
		}

		// Determine breaching strategy
		depression.should_breach = should_breach_depression(depression);

		// Find breach path if needed
		if (depression.should_breach && depression.outlet_index != SIZE_MAX) {
			depression.breach_path = find_breach_path_corrected(depression);
		}
	}

	/**
	 * Improved breach path finding that ensures valid path
	 */
	std::vector<size_t> find_breach_path_corrected(
		const Depression<T>& depression) const
	{
		if (depression.outlet_index == SIZE_MAX || depression.cells.empty()) {
			return {};
		}

		// Find the depression cell closest to the outlet
		size_t start_idx = depression.cells[0];
		T min_distance =
			connector_->euclidean_distance(start_idx, depression.outlet_index);

		for (size_t cell_idx : depression.cells) {
			T distance =
				connector_->euclidean_distance(cell_idx, depression.outlet_index);
			if (distance < min_distance) {
				min_distance = distance;
				start_idx = cell_idx;
			}
		}

		// Use A* pathfinding to outlet, but ensure path goes through depression
		// boundary
		auto path = find_astar_path(start_idx, depression.outlet_index);

		// Validate path doesn't go through other depressions inappropriately
		if (!path.empty() && is_valid_breach_path(path, depression)) {
			return path;
		}

		return {}; // No valid breach path found
	}

	/**
	 * Validate that breach path is reasonable
	 */
	bool is_valid_breach_path(const std::vector<size_t>& path,
														const Depression<T>& depression) const
	{
		if (path.size() < 2) {
			return false;
		}

		// Check path doesn't go through significantly higher elevations
		T max_allowed_elevation =
			depression.outlet_elevation + config_.breach_threshold;

		for (size_t cell_idx : path) {
			if (working_elevation_[cell_idx] > max_allowed_elevation) {
				return false;
			}
		}

		// Path should start in depression and end at outlet
		bool starts_in_depression =
			std::find(depression.cells.begin(), depression.cells.end(), path[0]) !=
			depression.cells.end();

		return starts_in_depression && path.back() == depression.outlet_index;
	}

	/**
	 * Check if a cell is a local minimum
	 */
	bool is_local_minimum(size_t index) const
	{
		T center_elev = working_elevation_[index];
		auto neighbors = connector_->get_valid_neighbors(index);

		for (const auto& neighbor : neighbors) {
			if (working_elevation_[neighbor.index] < center_elev) {
				return false;
			}
		}

		return !neighbors.empty(); // Not a minimum if no neighbors
	}

	/**
	 * Flood fill to identify entire depression
	 */
	void flood_fill_depression(size_t start_index,
														 Depression<T>& depression,
														 std::vector<bool>& visited) const
	{
		std::queue<size_t> queue;
		queue.push(start_index);
		visited[start_index] = true;

		T current_level = working_elevation_[start_index];
		depression.outlet_elevation = std::numeric_limits<T>::max();
		depression.outlet_index = SIZE_MAX;

		while (!queue.empty()) {
			size_t current = queue.front();
			queue.pop();

			depression.cells.push_back(current);
			depression.area++;
			cell_to_depression_[current] = depression.id;

			auto neighbors = connector_->get_valid_neighbors(current);
			for (const auto& neighbor : neighbors) {
				T neighbor_elev = working_elevation_[neighbor.index];

				if (!visited[neighbor.index]) {
					if (neighbor_elev <= current_level + config_.epsilon) {
						// Part of the same depression
						visited[neighbor.index] = true;
						queue.push(neighbor.index);
						current_level = std::max(current_level, neighbor_elev);
					} else {
						// Potential outlet
						if (neighbor_elev < depression.outlet_elevation) {
							depression.outlet_elevation = neighbor_elev;
							depression.outlet_index = neighbor.index;
						}
					}
				}
			}
		}

		depression.fill_level = depression.outlet_elevation;
	}

	/**
	 * Calculate depression properties (volume, etc.)
	 */
	void calculate_depression_properties(Depression<T>& depression) const
	{
		depression.volume = 0;
		depression.min_elevation = std::numeric_limits<T>::max();

		for (size_t cell_idx : depression.cells) {
			T cell_elev = working_elevation_[cell_idx];
			depression.min_elevation = std::min(depression.min_elevation, cell_elev);

			if (depression.fill_level > cell_elev) {
				depression.volume += (depression.fill_level - cell_elev);
			}
		}

		// Determine if depression should be breached
		depression.should_breach = should_breach_depression(depression);
	}

	/**
	 * Determine if a depression should be breached rather than filled
	 */
	bool should_breach_depression(const Depression<T>& depression) const
	{
		// Size-based decision
		if (depression.volume > config_.breach_threshold)
			return true;

		// Flow-based decision
		if (has_significant_flow(depression))
			return true;

		// Depth-based decision
		T max_depth = depression.fill_level - depression.min_elevation;
		if (max_depth > config_.breach_threshold)
			return true;

		return false;
	}

	/**
	 * Check if depression has significant flow accumulation
	 */
	bool has_significant_flow(const Depression<T>& depression) const
	{
		if (flow_accumulation_.empty())
			return false;

		T total_flow = 0;
		for (size_t cell_idx : depression.cells) {
			total_flow += flow_accumulation_[cell_idx];
		}

		T average_flow = total_flow / depression.area;
		return average_flow >
					 (config_.flow_accumulation_weight * size_ * 0.001); // 0.1% of grid
	}

	/**
	 * Compute flow accumulation for flow-based methods
	 */
	void compute_flow_accumulation() const
	{
		std::fill(flow_accumulation_.begin(), flow_accumulation_.end(), 1.0);

		// Simple D8 flow accumulation
		std::vector<size_t> processing_order = get_processing_order();

		for (size_t idx : processing_order) {
			if (!connector_->is_active_node(idx))
				continue;

			// Find steepest descent neighbor
			size_t receiver = find_steepest_receiver(idx);
			if (receiver != SIZE_MAX && connector_->is_active_node(receiver)) {
				flow_accumulation_[receiver] += flow_accumulation_[idx];
			}
		}
	}

	/**
	 * Get processing order for flow accumulation (topological sort)
	 */
	std::vector<size_t> get_processing_order() const
	{
		std::vector<size_t> order;
		std::vector<int> in_degree(size_, 0);

		// Calculate in-degrees
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			size_t receiver = find_steepest_receiver(i);
			if (receiver != SIZE_MAX) {
				in_degree[receiver]++;
			}
		}

		// Topological sort
		std::queue<size_t> queue;
		for (size_t i = 0; i < size_; ++i) {
			if (connector_->is_active_node(i) && in_degree[i] == 0) {
				queue.push(i);
			}
		}

		while (!queue.empty()) {
			size_t current = queue.front();
			queue.pop();
			order.push_back(current);

			size_t receiver = find_steepest_receiver(current);
			if (receiver != SIZE_MAX) {
				in_degree[receiver]--;
				if (in_degree[receiver] == 0) {
					queue.push(receiver);
				}
			}
		}

		return order;
	}

	/**
	 * Find steepest descent neighbor
	 */
	size_t find_steepest_receiver(size_t index) const
	{
		T center_elev = working_elevation_[index];
		auto neighbors = connector_->get_valid_neighbors(index);

		size_t steepest_neighbor = SIZE_MAX;
		T max_gradient = -std::numeric_limits<T>::infinity();

		for (const auto& neighbor : neighbors) {
			T neighbor_elev = working_elevation_[neighbor.index];
			T gradient = (center_elev - neighbor_elev) / neighbor.distance;

			if (gradient > max_gradient && gradient > config_.min_gradient) {
				max_gradient = gradient;
				steepest_neighbor = neighbor.index;
			}
		}

		return steepest_neighbor;
	}

	/**
	 * Fill depression optimally
	 */
	void fill_depression_optimal(const Depression<T>& depression,
															 FillingResult<T>& result) const
	{
		for (size_t cell_idx : depression.cells) {
			T original_elev = (*elevation_)[cell_idx];
			T new_elev = std::max(original_elev, depression.fill_level);

			if (new_elev > working_elevation_[cell_idx] &&
					should_modify_cell(cell_idx, new_elev)) {
				T volume_added = new_elev - working_elevation_[cell_idx];
				working_elevation_[cell_idx] = new_elev;

				result.total_volume_filled += volume_added;
				result.cells_modified++;
				result.max_modification =
					std::max(result.max_modification, new_elev - original_elev);
			}
		}

		result.depressions_filled++;
	}

	/**
	 * Breach depression by carving outlet channel
	 */
	void breach_depression(const Depression<T>& depression,
												 FillingResult<T>& result) const
	{
		if (depression.outlet_index == SIZE_MAX) {
			// No outlet found, fall back to filling
			fill_depression_optimal(depression, result);
			return;
		}

		// Find path from depression to outlet
		auto path = find_breach_path(depression);
		if (path.empty()) {
			fill_depression_optimal(depression, result);
			return;
		}

		// Carve the path
		carve_path(path, depression, result);
		result.outlets_created++;
	}

	/**
	 * Find optimal breach path
	 */
	std::vector<size_t> find_breach_path(const Depression<T>& depression) const
	{
		if (depression.outlet_index == SIZE_MAX)
			return {};

		// Find lowest point in depression as starting point
		size_t start_idx = depression.cells[0];
		T min_elev = working_elevation_[start_idx];

		for (size_t cell_idx : depression.cells) {
			if (working_elevation_[cell_idx] < min_elev) {
				min_elev = working_elevation_[cell_idx];
				start_idx = cell_idx;
			}
		}

		// Use A* to find path to outlet
		return find_astar_path(start_idx, depression.outlet_index);
	}

	/**
	 * A* pathfinding for breach paths
	 */
	std::vector<size_t> find_astar_path(size_t start, size_t goal) const
	{
		struct AStarNode
		{
			size_t index;
			T g_cost; // Distance from start
			T h_cost; // Heuristic distance to goal
			T f_cost() const { return g_cost + h_cost; }
			size_t parent;

			// Default constructor - THIS WAS MISSING
			AStarNode()
				: index(SIZE_MAX)
				, g_cost(static_cast<T>(0))
				, h_cost(static_cast<T>(0))
				, parent(SIZE_MAX)
			{
			}

			// Parameterized constructor
			AStarNode(size_t idx, T g, T h, size_t p)
				: index(idx)
				, g_cost(g)
				, h_cost(h)
				, parent(p)
			{
			}
		};

		auto compare = [](const AStarNode& a, const AStarNode& b) {
			return a.f_cost() > b.f_cost();
		};

		std::priority_queue<AStarNode, std::vector<AStarNode>, decltype(compare)>
			open_set(compare);
		std::unordered_set<size_t> closed_set;
		std::unordered_map<size_t, AStarNode> nodes;

		T h_start = static_cast<T>(connector_->euclidean_distance(start, goal));
		open_set.emplace(start, static_cast<T>(0), h_start, SIZE_MAX);
		nodes.emplace(start,
									AStarNode(start, static_cast<T>(0), h_start, SIZE_MAX));

		while (!open_set.empty()) {
			AStarNode current = open_set.top();
			open_set.pop();

			if (current.index == goal) {
				// Reconstruct path
				std::vector<size_t> path;
				size_t idx = goal;
				while (idx != SIZE_MAX) {
					path.push_back(idx);
					idx = nodes[idx].parent;
				}
				std::reverse(path.begin(), path.end());
				return path;
			}

			closed_set.insert(current.index);

			auto neighbors = connector_->get_valid_neighbors(current.index);
			for (const auto& neighbor : neighbors) {
				if (closed_set.count(neighbor.index))
					continue;

				// Cost includes elevation gain and distance
				T elevation_cost = std::max(static_cast<T>(0),
																		working_elevation_[neighbor.index] -
																			working_elevation_[current.index]);
				T tentative_g = current.g_cost + neighbor.distance +
												elevation_cost * static_cast<T>(10);
				T h =
					static_cast<T>(connector_->euclidean_distance(neighbor.index, goal));

				auto it = nodes.find(neighbor.index);
				if (it == nodes.end() || tentative_g < it->second.g_cost) {
					nodes[neighbor.index] =
						AStarNode(neighbor.index, tentative_g, h, current.index);
					open_set.emplace(neighbor.index, tentative_g, h, current.index);
				}
			}
		}

		return {}; // No path found
	}

	/**
	 * Carve elevation along a path
	 */
	void carve_path(const std::vector<size_t>& path,
									const Depression<T>& depression,
									FillingResult<T>& result) const
	{
		if (path.size() < 2)
			return;

		// Calculate target elevations along path
		T start_elev = depression.min_elevation;
		T end_elev = working_elevation_[path.back()];
		T elevation_drop = start_elev - end_elev;

		// Ensure minimum gradient
		T min_drop_per_cell = config_.slope_threshold;
		if (elevation_drop < min_drop_per_cell * path.size()) {
			end_elev = start_elev - min_drop_per_cell * path.size();
		}

		for (size_t i = 0; i < path.size(); ++i) {
			size_t cell_idx = path[i];
			if (!should_modify_cell(cell_idx, working_elevation_[cell_idx]))
				continue;

			// Linear interpolation of target elevation
			T progress = static_cast<T>(i) / (path.size() - 1);
			T target_elev = start_elev - progress * (start_elev - end_elev);

			if (target_elev < working_elevation_[cell_idx]) {
				T volume_carved = working_elevation_[cell_idx] - target_elev;
				working_elevation_[cell_idx] = target_elev;

				result.total_volume_carved += volume_carved;
				result.cells_modified++;
				result.max_modification =
					std::max(result.max_modification, volume_carved);
			}
		}
	}

	/**
	 * Create sloped facets instead of flat areas
	 */
	void create_sloped_facets(const Depression<T>& depression,
														FillingResult<T>& result) const
	{
		if (depression.cells.empty())
			return;

		// Find outlets around depression
		std::vector<std::pair<size_t, T>> outlets;
		for (size_t cell_idx : depression.cells) {
			auto neighbors = connector_->get_valid_neighbors(cell_idx);
			for (const auto& neighbor : neighbors) {
				if (std::find(depression.cells.begin(),
											depression.cells.end(),
											neighbor.index) == depression.cells.end()) {
					// This neighbor is outside the depression
					outlets.emplace_back(neighbor.index,
															 working_elevation_[neighbor.index]);
				}
			}
		}

		if (outlets.empty()) {
			fill_depression_optimal(depression, result);
			return;
		}

		// Create gradual slopes toward outlets
		for (size_t cell_idx : depression.cells) {
			if (!should_modify_cell(cell_idx, working_elevation_[cell_idx]))
				continue;

			// Find distance-weighted elevation from outlets
			T weighted_elevation = 0;
			T total_weight = 0;

			for (const auto& [outlet_idx, outlet_elev] : outlets) {
				T distance =
					static_cast<T>(connector_->euclidean_distance(cell_idx, outlet_idx));
				T weight = 1.0 / (distance + config_.epsilon);

				weighted_elevation += outlet_elev * weight;
				total_weight += weight;
			}

			if (total_weight > 0) {
				T target_elev = weighted_elevation / total_weight;
				T original_elev = (*elevation_)[cell_idx];

				// Ensure we don't lower below original
				target_elev = std::max(target_elev, original_elev);

				if (target_elev > working_elevation_[cell_idx]) {
					T volume_added = target_elev - working_elevation_[cell_idx];
					working_elevation_[cell_idx] = target_elev;

					result.total_volume_filled += volume_added;
					result.cells_modified++;
					result.max_modification =
						std::max(result.max_modification, target_elev - original_elev);
				}
			}
		}

		result.depressions_filled++;
	}

	/**
	 * Check if cell is on a ridge line
	 */
	bool is_ridge_cell(size_t index) const
	{
		T center_elev = working_elevation_[index];
		auto neighbors = connector_->get_valid_neighbors(index);

		if (neighbors.size() < 3)
			return false; // Need sufficient neighbors

		size_t higher_count = 0;
		for (const auto& neighbor : neighbors) {
			if (working_elevation_[neighbor.index] < center_elev) {
				higher_count++;
			}
		}

		// Ridge if most neighbors are lower
		return higher_count >= neighbors.size() * 0.75;
	}

	/**
	 * Check if depression is in a valley
	 */
	bool is_valley_depression(const Depression<T>& depression) const
	{
		// Count surrounding high elevation areas
		std::unordered_set<size_t> depression_set(depression.cells.begin(),
																							depression.cells.end());
		size_t high_neighbors = 0;
		size_t total_neighbors = 0;

		for (size_t cell_idx : depression.cells) {
			auto neighbors = connector_->get_valid_neighbors(cell_idx);
			for (const auto& neighbor : neighbors) {
				if (depression_set.find(neighbor.index) == depression_set.end()) {
					total_neighbors++;
					if (working_elevation_[neighbor.index] > depression.fill_level) {
						high_neighbors++;
					}
				}
			}
		}

		return total_neighbors > 0 &&
					 static_cast<double>(high_neighbors) / total_neighbors > 0.6;
	}

	// Additional placeholder implementations for complex methods
	// These would be fully implemented based on specific algorithm requirements

	void fill_depression_minimal_impact(const Depression<T>& depression,
																			FillingResult<T>& result) const
	{
		// Simplified implementation - would use optimization to minimize total
		// modification
		fill_depression_optimal(depression, result);
	}

	bool affects_ridges(const Depression<T>& depression) const
	{
		for (size_t cell_idx : depression.cells) {
			if (is_ridge_cell(cell_idx))
				return true;
		}
		return false;
	}

	std::vector<T> estimate_fill_impact(const Depression<T>& depression) const
	{
		return {
			depression.volume
		}; // Simplified - would compute multiple impact metrics
	}

	std::vector<T> estimate_breach_impact(const Depression<T>& depression) const
	{
		auto path = find_breach_path(depression);
		return { static_cast<T>(path.size() * config_.epsilon) }; // Simplified
	}

	std::vector<T> estimate_carve_impact(const Depression<T>& depression) const
	{
		return { depression.volume *
						 0.5 }; // Simplified - would compute carving cost
	}

	std::vector<size_t> find_least_cost_path(
		const Depression<T>& depression) const
	{
		return find_breach_path(
			depression); // Simplified - would use more sophisticated cost function
	}

	T estimate_fill_cost(const Depression<T>& depression) const
	{
		return depression.volume;
	}

	T estimate_breach_cost(const Depression<T>& depression) const
	{
		auto path = find_breach_path(depression);
		return static_cast<T>(path.size());
	}

	T get_depression_flow_weight(const Depression<T>& depression) const
	{
		if (flow_accumulation_.empty())
			return 0;

		T total_flow = 0;
		for (size_t cell_idx : depression.cells) {
			total_flow += flow_accumulation_[cell_idx];
		}
		return total_flow / depression.area;
	}

	std::vector<size_t> trace_flow_path(const Depression<T>& depression) const
	{
		return find_breach_path(depression); // Simplified
	}

	void carve_flow_channel(const std::vector<size_t>& path,
													const Depression<T>& depression,
													FillingResult<T>& result) const
	{
		carve_path(path, depression, result); // Simplified
	}

	std::vector<bool> create_constraint_map() const
	{
		std::vector<bool> constraints(size_, false);

		// Mark ridge cells as constrained if preserve_ridges is enabled
		if (config_.preserve_ridges) {
			for (size_t i = 0; i < size_; ++i) {
				if (connector_->is_active_node(i) && is_ridge_cell(i)) {
					constraints[i] = true;
				}
			}
		}

		return constraints;
	}

	bool violates_constraints(const Depression<T>& depression,
														const std::vector<bool>& constraints) const
	{
		for (size_t cell_idx : depression.cells) {
			if (constraints[cell_idx])
				return true;
		}
		return false;
	}

	void fill_depression_constrained(const Depression<T>& depression,
																	 const std::vector<bool>& constraints,
																	 FillingResult<T>& result) const
	{
		// Simplified - would implement constrained optimization
		fill_depression_optimal(depression, result);
	}

	// Additional method stubs for completeness
	std::vector<bool> identify_morphological_features() const
	{
		return std::vector<bool>(size_, false);
	}

	bool affects_morphological_features(
		const Depression<T>& depression,
		const std::vector<bool>& morphological_map) const
	{
		return false;
	}

	void breach_depression_morphology_aware(
		const Depression<T>& depression,
		const std::vector<bool>& morphological_map,
		FillingResult<T>& result) const
	{
		breach_depression(depression, result);
	}

	std::vector<std::vector<size_t>> build_depression_hierarchy() const
	{
		std::vector<std::vector<size_t>> hierarchy;
		if (!depressions_.empty()) {
			std::vector<size_t> all_depressions;
			for (size_t i = 0; i < depressions_.size(); ++i) {
				all_depressions.push_back(i);
			}
			hierarchy.push_back(all_depressions);
		}
		return hierarchy;
	}

	bool is_primary_depression(
		const Depression<T>& depression,
		const std::vector<std::vector<size_t>>& hierarchy) const
	{
		return true; // Simplified
	}

	bool can_breach_safely(const Depression<T>& depression) const
	{
		return depression.outlet_index != SIZE_MAX;
	}

	std::vector<bool> identify_hydrological_features() const
	{
		return std::vector<bool>(size_, false);
	}

	bool is_hydrologically_significant(const Depression<T>& depression,
																		 const std::vector<bool>& hydro_map) const
	{
		return has_significant_flow(depression);
	}

	void create_hydrological_channel(const Depression<T>& depression,
																	 const std::vector<bool>& hydro_map,
																	 FillingResult<T>& result) const
	{
		breach_depression(depression, result);
	}

	T evaluate_terrain_modification_score(const Depression<T>& depression) const
	{
		return depression.volume + depression.area * config_.morphological_weight;
	}

	void carve_terrain_sophisticated(const Depression<T>& depression,
																	 FillingResult<T>& result) const
	{
		breach_depression(depression, result);
	}

	T calculate_fill_probability(const Depression<T>& depression) const
	{
		T size_factor = std::min(static_cast<T>(1.0),
														 depression.volume / config_.breach_threshold);
		return 1.0 - size_factor; // Larger depressions less likely to be filled
	}

	void breach_depression_stochastic(const Depression<T>& depression,
																		FillingResult<T>& result) const
	{
		breach_depression(depression, result); // Simplified
	}

	// ======================
	// CORDONNIER ET AL. 2019 METHODS
	// ======================

	/**
	 * Initialize Cordonnier router if needed
	 */
	void initialize_cordonnier_router() const
	{
		if (!cordonnier_router_) {
			cordonnier_router_ = std::make_unique<CordonnierFlowRouter<T>>(
				connector_, elevation_, FlowEnforcementStrategy::SIMPLE_CORRECTION);
		}
	}

	/**
	 * Cordonnier Simple Correction Method
	 */
	FillingResult<T> cordonnier_simple_method() const
	{
		FillingResult<T> result;

		initialize_cordonnier_router();
		cordonnier_router_->set_strategy(
			FlowEnforcementStrategy::SIMPLE_CORRECTION);
		cordonnier_router_->compute_flow_routing();

		// Apply the flow routing results to working elevation
		// The Cordonnier method doesn't modify elevations, only flow paths
		// So we keep the original elevations but could optionally apply flow-based
		// corrections

		result.success = true;
		result.iterations_used = 1;

		// The Cordonnier method doesn't fill depressions in the traditional sense
		// It routes flow through them without changing elevations
		result.total_volume_filled = 0;
		result.cells_modified = 0;

		auto stats = cordonnier_router_->compute_statistics();
		result.depressions_filled = stats.num_inner_basins;

		return result;
	}

	/**
	 * Cordonnier Carving Method
	 */
	FillingResult<T> cordonnier_carving_method() const
	{
		FillingResult<T> result;

		initialize_cordonnier_router();
		cordonnier_router_->set_strategy(
			FlowEnforcementStrategy::DEPRESSION_CARVING);
		cordonnier_router_->compute_flow_routing();

		// Get the flow receivers and apply carving-like elevation modifications
		const auto& receivers = cordonnier_router_->get_receivers();
		const auto& basin_labels = cordonnier_router_->get_basin_labels();
		const auto& basins = cordonnier_router_->get_basins();

		// Apply minimal carving to enforce flow paths
		for (const auto& basin : basins) {
			if (!basin.is_boundary_basin && basin.spill_node != SIZE_MAX) {
				// Trace and slightly carve the flow path
				carve_cordonnier_path(basin, receivers, result);
			}
		}

		result.success = true;
		result.iterations_used = 1;

		auto stats = cordonnier_router_->compute_statistics();
		result.depressions_filled = stats.num_inner_basins;

		return result;
	}

	/**
	 * Cordonnier Filling Method
	 */
	FillingResult<T> cordonnier_filling_method() const
	{
		FillingResult<T> result;

		initialize_cordonnier_router();
		cordonnier_router_->set_strategy(
			FlowEnforcementStrategy::DEPRESSION_FILLING);
		cordonnier_router_->compute_flow_routing();

		// Get water levels and apply them as filled elevations
		const auto& water_levels = cordonnier_router_->get_water_levels();
		const auto& basin_labels = cordonnier_router_->get_basin_labels();
		const auto& basins = cordonnier_router_->get_basins();

		// Apply water levels as filled elevations for inner basins
		for (const auto& basin : basins) {
			if (!basin.is_boundary_basin) {
				for (size_t node : basin.nodes) {
					T original_elev = (*elevation_)[node];
					T water_level = water_levels[node];

					if (water_level > original_elev &&
							should_modify_cell(node, water_level)) {
						T volume_added = water_level - working_elevation_[node];
						working_elevation_[node] = water_level;

						result.total_volume_filled += volume_added;
						result.cells_modified++;
						result.max_modification =
							std::max(result.max_modification, water_level - original_elev);
					}
				}
			}
		}

		result.success = true;
		result.iterations_used = 1;
		result.mean_modification =
			result.cells_modified > 0
				? result.total_volume_filled / result.cells_modified
				: 0;

		auto stats = cordonnier_router_->compute_statistics();
		result.depressions_filled = stats.num_inner_basins;

		return result;
	}

	/**
	 * Apply minimal carving along Cordonnier flow path
	 */
	void carve_cordonnier_path(const Basin<T>& basin,
														 const std::vector<size_t>& receivers,
														 FillingResult<T>& result) const
	{
		if (basin.nodes.empty() || basin.spill_node == SIZE_MAX)
			return;

		// Find the flow path from basin minimum to spill using Cordonnier receivers
		std::vector<size_t> path;
		size_t current = basin.singular_node;
		std::unordered_set<size_t> visited;

		while (current != SIZE_MAX && visited.find(current) == visited.end() &&
					 current != basin.spill_node) {
			path.push_back(current);
			visited.insert(current);

			// Follow the Cordonnier receivers
			current = receivers[current];
		}

		if (current == basin.spill_node) {
			path.push_back(current);
		}

		// Apply minimal carving along this path
		if (path.size() > 1) {
			T start_elev = working_elevation_[path[0]];
			T end_elev = working_elevation_[path.back()];

			// Ensure monotonic decrease
			for (size_t i = 1; i < path.size(); ++i) {
				size_t node = path[i];
				size_t prev_node = path[i - 1];

				T prev_elev = working_elevation_[prev_node];
				T current_elev = working_elevation_[node];

				// Ensure this node is lower than previous
				T target_elev = prev_elev - config_.epsilon;

				if (current_elev > target_elev &&
						should_modify_cell(node, target_elev)) {
					T volume_carved = working_elevation_[node] - target_elev;
					working_elevation_[node] = target_elev;

					result.total_volume_carved += volume_carved;
					result.cells_modified++;
					result.max_modification =
						std::max(result.max_modification, volume_carved);
				}
			}
		}
	}

public:
	// ======================
	// ANALYSIS AND VALIDATION
	// ======================

	/**
	 * Validate that all depressions have been properly handled
	 */
	bool validate_drainage() const
	{
		// Check that no local minima remain (except at boundaries)
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			NodeType boundary_type = connector_->get_boundary_type(i);
			if (NodeTypeUtils::allows_outflow(boundary_type))
				continue;

			if (is_local_minimum(i)) {
				return false; // Found remaining local minimum
			}
		}

		return true;
	}

	/**
	 * Get statistics about modifications made
	 */
	struct ModificationStatistics
	{
		size_t total_cells_modified;
		T total_volume_change;
		T max_elevation_increase;
		T max_elevation_decrease;
		T mean_elevation_change;
		T std_elevation_change;
		size_t cells_filled;
		size_t cells_carved;
		std::vector<std::pair<size_t, T>> largest_modifications;

		ModificationStatistics()
			: total_cells_modified(0)
			, total_volume_change(0)
			, max_elevation_increase(0)
			, max_elevation_decrease(0)
			, mean_elevation_change(0)
			, std_elevation_change(0)
			, cells_filled(0)
			, cells_carved(0)
		{
		}
	};

	ModificationStatistics analyze_modifications() const
	{
		ModificationStatistics stats;

		std::vector<T> changes;
		changes.reserve(size_);

		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			T original = (*elevation_)[i];
			T modified = working_elevation_[i];
			T change = modified - original;

			if (std::abs(change) > config_.epsilon) {
				stats.total_cells_modified++;
				stats.total_volume_change += std::abs(change);
				changes.push_back(change);

				if (change > 0) {
					stats.cells_filled++;
					stats.max_elevation_increase =
						std::max(stats.max_elevation_increase, change);
				} else {
					stats.cells_carved++;
					stats.max_elevation_decrease =
						std::min(stats.max_elevation_decrease, change);
				}

				// Track largest modifications
				if (std::abs(change) > config_.breach_threshold) {
					stats.largest_modifications.emplace_back(i, change);
				}
			}
		}

		// Calculate statistics
		if (!changes.empty()) {
			T sum =
				std::accumulate(changes.begin(), changes.end(), static_cast<T>(0));
			stats.mean_elevation_change = sum / changes.size();

			T sq_sum = 0;
			for (T change : changes) {
				sq_sum += (change - stats.mean_elevation_change) *
									(change - stats.mean_elevation_change);
			}
			stats.std_elevation_change = std::sqrt(sq_sum / changes.size());
		}

		// Sort largest modifications by magnitude
		std::sort(stats.largest_modifications.begin(),
							stats.largest_modifications.end(),
							[](const auto& a, const auto& b) {
								return std::abs(a.second) > std::abs(b.second);
							});

		return stats;
	}

	/**
	 * Export depression map for visualization
	 */
	std::vector<size_t> get_depression_map() const
	{
		std::vector<size_t> depression_map(size_, SIZE_MAX);

		for (const auto& [cell_idx, depression_id] : cell_to_depression_) {
			depression_map[cell_idx] = depression_id;
		}

		return depression_map;
	}

	/**
	 * Get modification map showing elevation changes
	 */
	std::vector<T> get_modification_map() const
	{
		std::vector<T> modification_map(size_);

		for (size_t i = 0; i < size_; ++i) {
			modification_map[i] = working_elevation_[i] - (*elevation_)[i];
		}

		return modification_map;
	}

	/**
	 * Generate summary report
	 */
	std::string generate_report() const
	{
		std::ostringstream report;

		report << "=== DEPRESSION FILLING REPORT ===\n\n";

		// Method information
		report << "Method: ";
		switch (config_.method) {
			case FillingMethod::PRIORITY_FLOOD_EPSILON:
				report << "Priority Flood + Epsilon (Barnes et al. 2014)";
				break;
			case FillingMethod::PRIORITY_FLOOD_OPTIMAL:
				report << "Priority Flood Optimal (Wang & Liu 2006)";
				break;
			case FillingMethod::HYBRID_PRIORITY_FLOODING:
				report << "Hybrid Priority Flooding (Cordonnier et al. 2019)";
				break;
			case FillingMethod::SLOPED_FACETS:
				report << "Sloped Facets (Zhou et al. 2016)";
				break;
			case FillingMethod::LEAST_COST_PATH:
				report << "Least Cost Path Carving";
				break;
			default:
				report << "Custom Method";
				break;
		}
		report << "\n";

		// Configuration
		report << "Epsilon: " << config_.epsilon << "\n";
		report << "Breach Threshold: " << config_.breach_threshold << "\n";
		report << "Preserve Ridges: " << (config_.preserve_ridges ? "Yes" : "No")
					 << "\n";
		report << "Preserve Valleys: " << (config_.preserve_valleys ? "Yes" : "No")
					 << "\n\n";

		// Depression statistics
		report << "=== DEPRESSION ANALYSIS ===\n";
		report << "Total Depressions Found: " << depressions_.size() << "\n";

		if (!depressions_.empty()) {
			T total_volume = 0;
			size_t total_area = 0;
			T max_volume = 0;
			T max_depth = 0;

			for (const auto& depression : depressions_) {
				total_volume += depression.volume;
				total_area += depression.area;
				max_volume = std::max(max_volume, depression.volume);
				max_depth =
					std::max(max_depth, depression.fill_level - depression.min_elevation);
			}

			report << "Total Depression Volume: " << total_volume << "\n";
			report << "Total Depression Area: " << total_area << " cells\n";
			report << "Average Depression Size: "
						 << (total_area / depressions_.size()) << " cells\n";
			report << "Largest Depression Volume: " << max_volume << "\n";
			report << "Maximum Depression Depth: " << max_depth << "\n\n";
		}

		// Modification statistics
		auto mod_stats = analyze_modifications();
		report << "=== MODIFICATION STATISTICS ===\n";
		report << "Cells Modified: " << mod_stats.total_cells_modified << " ("
					 << (100.0 * mod_stats.total_cells_modified / size_) << "%)\n";
		report << "Cells Filled: " << mod_stats.cells_filled << "\n";
		report << "Cells Carved: " << mod_stats.cells_carved << "\n";
		report << "Total Volume Change: " << mod_stats.total_volume_change << "\n";
		report << "Maximum Elevation Increase: " << mod_stats.max_elevation_increase
					 << "\n";
		report << "Maximum Elevation Decrease: " << mod_stats.max_elevation_decrease
					 << "\n";
		report << "Mean Elevation Change: " << mod_stats.mean_elevation_change
					 << "\n";
		report << "Std Elevation Change: " << mod_stats.std_elevation_change
					 << "\n\n";

		// Validation
		bool drainage_valid = validate_drainage();
		report << "=== VALIDATION ===\n";
		report << "Drainage Validation: " << (drainage_valid ? "PASSED" : "FAILED")
					 << "\n";

		if (!drainage_valid) {
			report << "Warning: Some local minima may remain in the processed DEM.\n";
		}

		// Largest modifications
		if (!mod_stats.largest_modifications.empty()) {
			report << "\n=== LARGEST MODIFICATIONS ===\n";
			size_t count = std::min(static_cast<size_t>(10),
															mod_stats.largest_modifications.size());
			for (size_t i = 0; i < count; ++i) {
				auto [idx, change] = mod_stats.largest_modifications[i];
				auto [row, col] = connector_->to_2d(idx);
				report << "Cell (" << row << ", " << col
							 << "): " << (change > 0 ? "+" : "") << change << "\n";
			}
		}

		return report.str();
	}

	// ======================
	// ADVANCED UTILITIES
	// ======================

	/**
	 * Apply custom depression treatment function
	 */
	template<typename TreatmentFunc>
	FillingResult<T> apply_custom_treatment(TreatmentFunc treatment_func) const
	{
		FillingResult<T> result;
		auto start_time = std::chrono::high_resolution_clock::now();

		reset_working_arrays();
		identify_depressions();

		for (const auto& depression : depressions_) {
			treatment_func(depression, working_elevation_, result);
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		result.processing_time =
			std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
																														start_time);
		result.success = true;

		return result;
	}

	/**
	 * Iterative refinement of depression filling
	 */
	FillingResult<T> iterative_refinement(size_t max_iterations = 5) const
	{
		FillingResult<T> result;

		for (size_t iter = 0; iter < max_iterations; ++iter) {
			auto iteration_result = fill_depressions();

			if (!iteration_result.success) {
				result = iteration_result;
				break;
			}

			// Accumulate results
			result.total_volume_filled += iteration_result.total_volume_filled;
			result.total_volume_carved += iteration_result.total_volume_carved;
			result.cells_modified += iteration_result.cells_modified;
			result.depressions_filled += iteration_result.depressions_filled;
			result.outlets_created += iteration_result.outlets_created;
			result.processing_time += iteration_result.processing_time;
			result.iterations_used = iter + 1;

			// Check convergence
			if (iteration_result.cells_modified == 0 ||
					iteration_result.total_volume_filled <
						config_.convergence_tolerance) {
				result.success = true;
				break;
			}

			// Update elevation for next iteration
			for (size_t i = 0; i < size_; ++i) {
				(*elevation_)[i] = working_elevation_[i];
			}
		}

		return result;
	}

	/**
	 * Parallel processing version (simplified implementation)
	 */
	FillingResult<T> fill_depressions_parallel() const
	{
		if (!config_.parallel_processing) {
			return fill_depressions();
		}

		// For now, fall back to sequential processing
		// In a full implementation, this would use thread pools and
		// partition the work among multiple threads
		return fill_depressions();
	}

	/**
	 * Progressive depression filling with user callback
	 */
	template<typename ProgressCallback>
	FillingResult<T> fill_depressions_progressive(
		ProgressCallback progress_callback) const
	{
		FillingResult<T> result;
		auto start_time = std::chrono::high_resolution_clock::now();

		reset_working_arrays();
		identify_depressions();

		size_t total_depressions = depressions_.size();
		size_t processed_depressions = 0;

		for (const auto& depression : depressions_) {
			// Process depression
			if (should_breach_depression(depression)) {
				breach_depression(depression, result);
			} else {
				fill_depression_optimal(depression, result);
			}

			processed_depressions++;

			// Call progress callback
			double progress =
				static_cast<double>(processed_depressions) / total_depressions;
			if (!progress_callback(
						progress, processed_depressions, total_depressions)) {
				// User requested cancellation
				result.success = false;
				break;
			}
		}

		if (processed_depressions == total_depressions) {
			result.success = true;
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		result.processing_time =
			std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
																														start_time);

		return result;
	}

	/**
	 * Create filled DEM with elevation constraints
	 */
	FillingResult<T> fill_with_constraints(const ArrayRef<T>& min_elevation,
																				 const ArrayRef<T>& max_elevation) const
	{
		if (min_elevation.size() != size_ || max_elevation.size() != size_) {
			throw std::invalid_argument("Constraint arrays must match grid size");
		}

		FillingResult<T> result = fill_depressions();

		if (!result.success)
			return result;

		// Apply constraints
		size_t constraint_violations = 0;
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			if (working_elevation_[i] < min_elevation[i]) {
				working_elevation_[i] = min_elevation[i];
				constraint_violations++;
			} else if (working_elevation_[i] > max_elevation[i]) {
				working_elevation_[i] = max_elevation[i];
				constraint_violations++;
			}
		}

		if (constraint_violations > 0 && config_.verbose) {
			std::cout << "Applied constraints to " << constraint_violations
								<< " cells\n";
		}

		return result;
	}

	/**
	 * Export results to different formats
	 */
	void export_results(const std::string& filename_prefix) const
	{
		// Export filled elevation
		auto filled_elevation = get_filled_elevation();
		// Implementation would depend on desired output format (ASCII, GeoTIFF,
		// etc.)

		// Export depression map
		auto depression_map = get_depression_map();
		// Implementation would export depression IDs

		// Export modification map
		auto modification_map = get_modification_map();
		// Implementation would export elevation changes

		// Export report
		std::string report = generate_report();
		// Implementation would write to text file
	}

	// ======================
	// FACTORY METHODS
	// ======================

	/**
	 * Create pre-configured depression filler for common scenarios
	 */
	static DepressionFiller create_for_hydrology(
		std::shared_ptr<Connector<T>> connector,
		std::shared_ptr<ArrayRef<T>> elevation)
	{
		FillingConfig<T> config;
		config.method = FillingMethod::HYBRID_PRIORITY_FLOODING;
		config.preserve_valleys = true;
		config.enforce_drainage = true;
		config.flow_accumulation_weight = 2.0;
		config.use_flow_barriers = true;

		return DepressionFiller(connector, elevation, config);
	}

	static DepressionFiller create_for_geomorphology(
		std::shared_ptr<Connector<T>> connector,
		std::shared_ptr<ArrayRef<T>> elevation)
	{
		FillingConfig<T> config;
		config.method = FillingMethod::IMPACT_REDUCTION;
		config.preserve_ridges = true;
		config.preserve_valleys = true;
		config.morphological_weight = 2.0;
		config.breach_threshold = 2.0; // Lower threshold for more breaching

		return DepressionFiller(connector, elevation, config);
	}

	static DepressionFiller create_for_minimal_modification(
		std::shared_ptr<Connector<T>> connector,
		std::shared_ptr<ArrayRef<T>> elevation)
	{
		FillingConfig<T> config;
		config.method = FillingMethod::CONSTRAINED_FLOW;
		config.epsilon = 1e-6; // Very small epsilon
		config.preserve_ridges = true;
		config.preserve_valleys = true;
		config.breach_threshold = 1.0; // Aggressive breaching

		return DepressionFiller(connector, elevation, config);
	}

	static DepressionFiller create_for_performance(
		std::shared_ptr<Connector<T>> connector,
		std::shared_ptr<ArrayRef<T>> elevation)
	{
		FillingConfig<T> config;
		config.method = FillingMethod::PRIORITY_FLOOD_EPSILON;
		config.parallel_processing = true;
		config.verbose = false;
		config.max_iterations = 100; // Limit iterations for speed

		return DepressionFiller(connector, elevation, config);
	}
};

// ======================
// CONVENIENCE FUNCTIONS
// ======================

/**
 * Quick depression filling for common use cases
 */
template<typename T = double>
class QuickFill
{
public:
	/**
	 * Simple priority flood with epsilon
	 */
	static ArrayRef<T> priority_flood(std::shared_ptr<Connector<T>> connector,
																		const ArrayRef<T>& elevation,
																		T epsilon = 1e-4)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FillingConfig<T> config;
		config.method = FillingMethod::PRIORITY_FLOOD_EPSILON;
		config.epsilon = epsilon;

		DepressionFiller<T> filler(connector, elevation_ref, config);
		auto result = filler.fill_depressions();

		if (!result.success) {
			throw std::runtime_error("Depression filling failed");
		}

		return filler.get_filled_elevation_ref();
	}

	/**
	 * Hybrid fill and breach approach
	 */
	static std::vector<T> hybrid_fill_breach(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		T breach_threshold = 5.0)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FillingConfig<T> config;
		config.method = FillingMethod::HYBRID_FILL_BREACH;
		config.breach_threshold = breach_threshold;

		DepressionFiller<T> filler(connector, elevation_ref, config);
		auto result = filler.fill_depressions();

		if (!result.success) {
			throw std::runtime_error("Depression filling failed");
		}

		return filler.get_filled_elevation();
	}

	/**
	 * Breach all depressions
	 */
	static std::vector<T> breach_all(std::shared_ptr<Connector<T>> connector,
																	 const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FillingConfig<T> config;
		config.method = FillingMethod::BREACH_DEPRESSIONS;

		DepressionFiller<T> filler(connector, elevation_ref, config);
		auto result = filler.fill_depressions();

		if (!result.success) {
			throw std::runtime_error("Depression filling failed");
		}

		return filler.get_filled_elevation();
	}

	/**
	 * Minimal impact filling
	 */
	static std::vector<T> minimal_impact(std::shared_ptr<Connector<T>> connector,
																			 const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FillingConfig<T> config;
		config.method = FillingMethod::IMPACT_REDUCTION;
		config.preserve_ridges = true;
		config.preserve_valleys = true;

		DepressionFiller<T> filler(connector, elevation_ref, config);
		auto result = filler.fill_depressions();

		if (!result.success) {
			throw std::runtime_error("Depression filling failed");
		}

		return filler.get_filled_elevation();
	}

	/**
	 * Linear complexity flow routing using basin graph (Cordonnier 2019)
	 */
	static std::vector<T> cordonnier_method(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		bool use_carving = false)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		FillingConfig<T> config;
		config.method = use_carving ? FillingMethod::CORDONNIER_CARVING
																: FillingMethod::CORDONNIER_FILLING;

		DepressionFiller<T> filler(connector, elevation_ref, config);
		auto result = filler.fill_depressions();

		if (!result.success) {
			throw std::runtime_error("Cordonnier depression filling failed");
		}

		return filler.get_filled_elevation();
	}
};

// Type aliases for common usage
using DepressionFillerF32 = DepressionFiller<float>;
using DepressionFillerF64 = DepressionFiller<double>;

} // namespace dagger2
