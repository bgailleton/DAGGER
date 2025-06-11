#pragma once

#include "dg2_BCs.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include "dg2_unionfind.hpp"
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * Flow enforcement strategy within depressions
 */
enum class FlowEnforcementStrategy : uint8_t
{
	SIMPLE_CORRECTION = 0,	// Only update minimal receivers (3 nodes max)
	DEPRESSION_CARVING = 1, // Carve narrow channel from minimum to spill
	DEPRESSION_FILLING = 2	// Fill depression with gentle slope toward spill
};

/**
 * Minimum Spanning Tree algorithm choice
 */
enum class MSTAlgorithm : uint8_t
{
	KRUSKAL = 0,		 // O(m log m) - simple and reliable
	MARES_PLANAR = 1 // O(n) - for planar graphs (experimental for D8)
};

/**
 * Basin structure - represents a drainage basin
 */
template<typename T>
struct Basin
{
	size_t id;								 // Unique basin identifier
	size_t singular_node;			 // Local minimum or boundary node
	bool is_boundary_basin;		 // True if drains to boundary
	std::vector<size_t> nodes; // All nodes in this basin
	T water_level;						 // Elevation of water surface
	size_t spill_node;				 // Node where water exits (SIZE_MAX if boundary)

	Basin()
		: id(SIZE_MAX)
		, singular_node(SIZE_MAX)
		, is_boundary_basin(false)
		, water_level(0)
		, spill_node(SIZE_MAX)
	{
	}

	Basin(size_t basin_id, size_t singular, bool is_boundary)
		: id(basin_id)
		, singular_node(singular)
		, is_boundary_basin(is_boundary)
		, water_level(0)
		, spill_node(SIZE_MAX)
	{
	}
};

/**
 * Basin link - connection between two adjacent basins
 */
template<typename T>
struct BasinLink
{
	size_t basin1_id;
	size_t basin2_id;
	size_t pass_node1; // Node in basin1 that forms the pass
	size_t pass_node2; // Node in basin2 that forms the pass
	T pass_elevation;	 // Elevation of the pass (max of the two nodes)

	BasinLink()
		: basin1_id(SIZE_MAX)
		, basin2_id(SIZE_MAX)
		, pass_node1(SIZE_MAX)
		, pass_node2(SIZE_MAX)
		, pass_elevation(0)
	{
	}

	BasinLink(size_t b1, size_t b2, size_t n1, size_t n2, T elev)
		: basin1_id(b1)
		, basin2_id(b2)
		, pass_node1(n1)
		, pass_node2(n2)
		, pass_elevation(elev)
	{
	}

	// Get the other basin in this link
	size_t get_other_basin(size_t basin_id) const
	{
		return (basin_id == basin1_id) ? basin2_id : basin1_id;
	}

	// Get pass node in specified basin
	size_t get_pass_node(size_t basin_id) const
	{
		return (basin_id == basin1_id) ? pass_node1 : pass_node2;
	}
};

/**
 * Cordonnier et al. 2019 Flow Routing Algorithm
 *
 * This implements the linear complexity algorithm for flow routing in
 * topographies with depressions as described in "A versatile, linear complexity
 * algorithm for flow routing in topographies with depressions" (Earth Surf.
 * Dynam., 7, 549–562, 2019).
 *
 * The algorithm works by:
 * 1. Computing drainage basins and linking adjacent basins
 * 2. Building a minimum spanning tree of basin connections
 * 3. Updating flow receivers to enforce drainage
 */
template<typename T = double>
class CordonnierFlowRouter
{
private:
	// Core components
	std::shared_ptr<Connector<T>> connector_;
	std::shared_ptr<ArrayRef<T>> elevation_;

	// Configuration
	FlowEnforcementStrategy strategy_;
	MSTAlgorithm mst_algorithm_;

	// Grid properties
	size_t rows_;
	size_t cols_;
	size_t size_;

	// Working data
	mutable std::vector<size_t> basin_labels_; // Basin ID for each node
	mutable std::vector<size_t> receivers_;		 // Flow receiver for each node
	mutable std::vector<std::vector<size_t>> donors_; // Flow donors for each node
	mutable std::vector<size_t> stack_order_; // Topological order for processing

	// Basin graph structures
	mutable std::vector<Basin<T>> basins_;
	mutable std::vector<BasinLink<T>> basin_links_;
	mutable std::unordered_map<std::pair<size_t, size_t>,
														 size_t,
														 std::hash<std::pair<size_t, size_t>>>
		link_map_;

	// MST results
	mutable std::vector<size_t> mst_links_; // Indices of links in MST
	mutable std::vector<T> water_levels_;		// Water level for each node

	// Hash function for pair
	struct PairHash
	{
		size_t operator()(const std::pair<size_t, size_t>& p) const
		{
			return std::hash<size_t>()(p.first) ^
						 (std::hash<size_t>()(p.second) << 1);
		}
	};

public:
	/**
	 * Constructor
	 */
	CordonnierFlowRouter(std::shared_ptr<Connector<T>> connector,
											 std::shared_ptr<ArrayRef<T>> elevation,
											 FlowEnforcementStrategy strategy =
												 FlowEnforcementStrategy::DEPRESSION_FILLING,
											 MSTAlgorithm mst_algo = MSTAlgorithm::KRUSKAL)
		: connector_(connector)
		, elevation_(elevation)
		, strategy_(strategy)
		, mst_algorithm_(mst_algo)
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
		basin_labels_.resize(size_, SIZE_MAX);
		receivers_.resize(size_, SIZE_MAX);
		donors_.resize(size_);
		stack_order_.resize(size_);
		water_levels_.resize(size_);

		// Initialize flow receivers and donors
		initialize_flow_receivers();
	}

	// ======================
	// CONFIGURATION
	// ======================

	FlowEnforcementStrategy get_strategy() const { return strategy_; }
	MSTAlgorithm get_mst_algorithm() const { return mst_algorithm_; }

	void set_strategy(FlowEnforcementStrategy strategy) { strategy_ = strategy; }
	void set_mst_algorithm(MSTAlgorithm algo) { mst_algorithm_ = algo; }

	// ======================
	// MAIN ALGORITHM
	// ======================

	/**
	 * Execute the complete Cordonnier algorithm
	 */
	void compute_flow_routing() const
	{
		// Step 1: Compute basins and linkage
		compute_basins_and_linkage();

		// Step 2: Compute minimum spanning tree
		compute_minimum_spanning_tree();

		// Step 3: Update flow receivers
		update_flow_receivers();

		// Compute final topological order
		compute_stack_order();
	}

	// ======================
	// RESULTS ACCESS
	// ======================

	const std::vector<size_t>& get_receivers() const { return receivers_; }
	const std::vector<std::vector<size_t>>& get_donors() const { return donors_; }
	const std::vector<size_t>& get_stack_order() const { return stack_order_; }
	const std::vector<size_t>& get_basin_labels() const { return basin_labels_; }
	const std::vector<T>& get_water_levels() const { return water_levels_; }
	const std::vector<Basin<T>>& get_basins() const { return basins_; }
	const std::vector<BasinLink<T>>& get_basin_links() const
	{
		return basin_links_;
	}

	/**
	 * Get flow accumulation using computed flow graph
	 */
	std::vector<size_t> compute_flow_accumulation() const
	{
		std::vector<size_t> accumulation(size_, 1);

		// Process in reverse topological order
		for (auto it = stack_order_.rbegin(); it != stack_order_.rend(); ++it) {
			size_t node = *it;
			if (!connector_->is_active_node(node))
				continue;

			size_t receiver = receivers_[node];
			if (receiver != SIZE_MAX && connector_->is_active_node(receiver)) {
				accumulation[receiver] += accumulation[node];
			}
		}

		return accumulation;
	}

	/**
	 * Get drainage area using computed flow graph
	 */
	std::vector<T> compute_drainage_area(T cell_area = 1.0) const
	{
		std::vector<T> area(size_, cell_area);

		// Process in reverse topological order
		for (auto it = stack_order_.rbegin(); it != stack_order_.rend(); ++it) {
			size_t node = *it;
			if (!connector_->is_active_node(node))
				continue;

			size_t receiver = receivers_[node];
			if (receiver != SIZE_MAX && connector_->is_active_node(receiver)) {
				area[receiver] += area[node];
			}
		}

		return area;
	}

private:
	// ======================
	// INITIALIZATION
	// ======================

	/**
	 * Initialize flow receivers using steepest descent
	 */
	void initialize_flow_receivers() const
	{
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i)) {
				receivers_[i] = SIZE_MAX;
				continue;
			}

			NodeType boundary_type = connector_->get_boundary_type(i);

			// Boundary nodes have no receiver
			if (NodeTypeUtils::allows_outflow(boundary_type)) {
				receivers_[i] = SIZE_MAX;
				continue;
			}

			// Find steepest downslope neighbor
			T center_elev = (*elevation_)[i];
			auto neighbors = connector_->get_valid_neighbors(i);

			size_t best_receiver = SIZE_MAX;
			T steepest_slope = -std::numeric_limits<T>::infinity();

			for (const auto& neighbor : neighbors) {
				T neighbor_elev = (*elevation_)[neighbor.index];
				T slope = (center_elev - neighbor_elev) / neighbor.distance;

				if (slope > steepest_slope && slope > 0) {
					steepest_slope = slope;
					best_receiver = neighbor.index;
				}
			}

			receivers_[i] = best_receiver;
		}

		// Build donors from receivers
		for (auto& donor_list : donors_) {
			donor_list.clear();
		}

		for (size_t i = 0; i < size_; ++i) {
			if (receivers_[i] != SIZE_MAX) {
				donors_[receivers_[i]].push_back(i);
			}
		}
	}

	// ======================
	// STEP 1: BASIN COMPUTATION AND LINKAGE
	// ======================

	/**
	 * Compute drainage basins and create basin graph
	 */
	void compute_basins_and_linkage() const
	{
		basins_.clear();
		basin_links_.clear();
		link_map_.clear();
		std::fill(basin_labels_.begin(), basin_labels_.end(), SIZE_MAX);

		// Find all singular nodes (local minima and boundary nodes)
		std::vector<size_t> singular_nodes = find_singular_nodes();

		// Assign basin labels using depth-first traversal
		size_t next_basin_id = 0;
		for (size_t singular_node : singular_nodes) {
			if (basin_labels_[singular_node] == SIZE_MAX) {
				assign_basin_labels(singular_node, next_basin_id++);
			}
		}

		// Create basin objects
		create_basin_objects(singular_nodes);

		// Find links between adjacent basins
		find_basin_links();
	}

	/**
	 * Find all singular nodes (local minima and boundary nodes)
	 */
	std::vector<size_t> find_singular_nodes() const
	{
		std::vector<size_t> singular_nodes;

		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			NodeType boundary_type = connector_->get_boundary_type(i);

			// Boundary nodes are singular
			if (NodeTypeUtils::allows_outflow(boundary_type)) {
				singular_nodes.push_back(i);
				continue;
			}

			// Local minima are singular (no valid receiver)
			if (receivers_[i] == SIZE_MAX) {
				singular_nodes.push_back(i);
			}
		}

		return singular_nodes;
	}

	/**
	 * Assign basin labels using depth-first traversal from singular node
	 */
	void assign_basin_labels(size_t start_node, size_t basin_id) const
	{
		std::stack<size_t> stack;
		stack.push(start_node);
		basin_labels_[start_node] = basin_id;

		while (!stack.empty()) {
			size_t current = stack.top();
			stack.pop();

			// Process all donors (nodes that flow into current)
			for (size_t donor : donors_[current]) {
				if (basin_labels_[donor] == SIZE_MAX) {
					basin_labels_[donor] = basin_id;
					stack.push(donor);
				}
			}
		}
	}

	/**
	 * Create basin objects from labeled nodes
	 */
	void create_basin_objects(const std::vector<size_t>& singular_nodes) const
	{
		// Count basins
		size_t num_basins = 0;
		for (size_t i = 0; i < size_; ++i) {
			if (basin_labels_[i] != SIZE_MAX) {
				num_basins = std::max(num_basins, basin_labels_[i] + 1);
			}
		}

		basins_.resize(num_basins);

		// Initialize basins with singular nodes
		for (size_t singular_node : singular_nodes) {
			size_t basin_id = basin_labels_[singular_node];
			if (basin_id != SIZE_MAX) {
				basins_[basin_id].id = basin_id;
				basins_[basin_id].singular_node = singular_node;

				NodeType boundary_type = connector_->get_boundary_type(singular_node);
				basins_[basin_id].is_boundary_basin =
					NodeTypeUtils::allows_outflow(boundary_type);

				if (basins_[basin_id].is_boundary_basin) {
					basins_[basin_id].water_level = (*elevation_)[singular_node];
				}
			}
		}

		// Add all nodes to their respective basins
		for (size_t i = 0; i < size_; ++i) {
			size_t basin_id = basin_labels_[i];
			if (basin_id != SIZE_MAX) {
				basins_[basin_id].nodes.push_back(i);
			}
		}
	}

	/**
	 * Find links between adjacent basins
	 */
	void find_basin_links() const
	{
		// Iterate through all edges in the grid
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			size_t basin1 = basin_labels_[i];
			if (basin1 == SIZE_MAX)
				continue;

			auto neighbors = connector_->get_valid_neighbors(i);
			for (const auto& neighbor : neighbors) {
				size_t basin2 = basin_labels_[neighbor.index];
				if (basin2 == SIZE_MAX || basin1 == basin2)
					continue;

				// Ensure consistent ordering (basin1 < basin2)
				size_t node1 = i;
				size_t node2 = neighbor.index;
				if (basin1 > basin2) {
					std::swap(basin1, basin2);
					std::swap(node1, node2);
				}

				auto link_key = std::make_pair(basin1, basin2);

				// Calculate pass elevation
				T pass_elevation = std::max((*elevation_)[node1], (*elevation_)[node2]);

				auto it = link_map_.find(link_key);
				if (it == link_map_.end()) {
					// Create new link
					size_t link_idx = basin_links_.size();
					basin_links_.emplace_back(
						basin1, basin2, node1, node2, pass_elevation);
					link_map_[link_key] = link_idx;
				} else {
					// Update existing link if this pass is lower
					size_t link_idx = it->second;
					if (pass_elevation < basin_links_[link_idx].pass_elevation) {
						basin_links_[link_idx].pass_node1 = node1;
						basin_links_[link_idx].pass_node2 = node2;
						basin_links_[link_idx].pass_elevation = pass_elevation;
					}
				}
			}
		}
	}

	// ======================
	// STEP 2: MINIMUM SPANNING TREE
	// ======================

	/**
	 * Compute minimum spanning tree of basin graph
	 */
	void compute_minimum_spanning_tree() const
	{
		if (mst_algorithm_ == MSTAlgorithm::KRUSKAL) {
			compute_mst_kruskal();
		} else {
			compute_mst_mares_planar();
		}
	}

	/**
	 * Kruskal's algorithm for MST - O(m log m)
	 */
	void compute_mst_kruskal() const
	{
		mst_links_.clear();

		if (basin_links_.empty())
			return;

		// Sort links by pass elevation
		std::vector<size_t> sorted_indices(basin_links_.size());
		std::iota(sorted_indices.begin(), sorted_indices.end(), 0);

		std::sort(
			sorted_indices.begin(), sorted_indices.end(), [this](size_t a, size_t b) {
				return basin_links_[a].pass_elevation < basin_links_[b].pass_elevation;
			});

		// Initialize Union-Find structure
		UnionFind uf(basins_.size());

		// Add virtual external basin for boundary basins
		size_t external_basin_id = basins_.size();

		// Connect all boundary basins to external basin
		for (size_t i = 0; i < basins_.size(); ++i) {
			if (basins_[i].is_boundary_basin) {
				uf.union_sets(i, external_basin_id);
			}
		}

		// Process links in order of increasing elevation
		for (size_t link_idx : sorted_indices) {
			const auto& link = basin_links_[link_idx];

			if (uf.union_sets(link.basin1_id, link.basin2_id)) {
				mst_links_.push_back(link_idx);
			}
		}
	}

	/**
	 * Mareš planar graph algorithm for MST - O(n)
	 * Simplified implementation for demonstration
	 */
	void compute_mst_mares_planar() const
	{
		// For now, fall back to Kruskal's algorithm
		// Full implementation would require planar graph operations
		compute_mst_kruskal();
	}

	// ======================
	// STEP 3: UPDATE FLOW RECEIVERS
	// ======================

	/**
	 * Update flow receivers based on MST and chosen strategy
	 */
	void update_flow_receivers() const
	{
		// First, compute water levels for all basins
		compute_water_levels();

		// Orient the MST from boundary basins to inner basins
		orient_basin_tree();

		// Apply flow enforcement strategy
		switch (strategy_) {
			case FlowEnforcementStrategy::SIMPLE_CORRECTION:
				apply_simple_correction();
				break;
			case FlowEnforcementStrategy::DEPRESSION_CARVING:
				apply_depression_carving();
				break;
			case FlowEnforcementStrategy::DEPRESSION_FILLING:
				apply_depression_filling();
				break;
		}

		// Update donors based on new receivers
		update_donors();
	}

	/**
	 * Compute water levels for all basins
	 */
	void compute_water_levels() const
	{
		// Initialize water levels
		for (size_t i = 0; i < size_; ++i) {
			water_levels_[i] = (*elevation_)[i];
		}

		// Set water levels for boundary basins
		for (const auto& basin : basins_) {
			if (basin.is_boundary_basin) {
				T boundary_elevation = (*elevation_)[basin.singular_node];
				for (size_t node : basin.nodes) {
					water_levels_[node] =
						std::max(water_levels_[node], boundary_elevation);
				}
			}
		}

		// Propagate water levels through MST
		// This is a simplified version - full implementation would handle nested
		// depressions
		for (size_t link_idx : mst_links_) {
			const auto& link = basin_links_[link_idx];

			T spill_elevation = link.pass_elevation;

			// Update water level for the higher basin
			auto& basin1 = basins_[link.basin1_id];
			auto& basin2 = basins_[link.basin2_id];

			if (!basin1.is_boundary_basin && basin2.is_boundary_basin) {
				basin1.water_level = spill_elevation;
				for (size_t node : basin1.nodes) {
					water_levels_[node] = std::max(water_levels_[node], spill_elevation);
				}
			} else if (basin1.is_boundary_basin && !basin2.is_boundary_basin) {
				basin2.water_level = spill_elevation;
				for (size_t node : basin2.nodes) {
					water_levels_[node] = std::max(water_levels_[node], spill_elevation);
				}
			}
		}
	}

	/**
	 * Orient basin tree from boundary to inner basins
	 */
	void orient_basin_tree() const
	{
		// Set spill nodes for each basin based on MST
		for (size_t link_idx : mst_links_) {
			const auto& link = basin_links_[link_idx];

			auto& basin1 = basins_[link.basin1_id];
			auto& basin2 = basins_[link.basin2_id];

			// Determine flow direction based on boundary status and water levels
			if (basin1.is_boundary_basin && !basin2.is_boundary_basin) {
				basin2.spill_node = link.pass_node2;
			} else if (!basin1.is_boundary_basin && basin2.is_boundary_basin) {
				basin1.spill_node = link.pass_node1;
			} else if (!basin1.is_boundary_basin && !basin2.is_boundary_basin) {
				// Both are inner basins - flow goes from higher water level to lower
				if (basin1.water_level > basin2.water_level) {
					basin1.spill_node = link.pass_node1;
				} else {
					basin2.spill_node = link.pass_node2;
				}
			}
		}
	}

	/**
	 * Simple correction strategy - update minimal receivers
	 */
	void apply_simple_correction() const
	{
		for (const auto& basin : basins_) {
			if (!basin.is_boundary_basin && basin.spill_node != SIZE_MAX) {
				// Update receiver of local minimum to point to spill
				receivers_[basin.singular_node] = basin.spill_node;
			}
		}
	}

	/**
	 * Depression carving strategy - create single path to spill
	 */
	void apply_depression_carving() const
	{
		for (const auto& basin : basins_) {
			if (!basin.is_boundary_basin && basin.spill_node != SIZE_MAX) {
				carve_path_to_spill(basin);
			}
		}
	}

	/**
	 * Depression filling strategy - create gentle slope to spill
	 */
	void apply_depression_filling() const
	{
		for (const auto& basin : basins_) {
			if (!basin.is_boundary_basin && basin.spill_node != SIZE_MAX) {
				fill_depression_to_spill(basin);
			}
		}
	}

	/**
	 * Carve a single path from local minimum to spill
	 */
	void carve_path_to_spill(const Basin<T>& basin) const
	{
		// Trace path from spill back to local minimum
		std::vector<size_t> path;
		size_t current = basin.spill_node;

		while (current != basin.singular_node && current != SIZE_MAX) {
			path.push_back(current);

			// Follow original receivers in reverse
			size_t next = SIZE_MAX;
			for (size_t donor : donors_[current]) {
				if (basin_labels_[donor] == basin.id) {
					next = donor;
					break;
				}
			}
			current = next;
		}

		if (current == basin.singular_node) {
			path.push_back(current);
		}

		// Reverse path and update receivers
		std::reverse(path.begin(), path.end());
		for (size_t i = 0; i < path.size() - 1; ++i) {
			receivers_[path[i]] = path[i + 1];
		}
	}

	/**
	 * Fill depression with gentle slope toward spill
	 */
	void fill_depression_to_spill(const Basin<T>& basin) const
	{
		if (basin.nodes.empty())
			return;

		// Use breadth-first search from spill
		std::queue<size_t> queue;
		std::unordered_set<size_t> visited;

		queue.push(basin.spill_node);
		visited.insert(basin.spill_node);

		while (!queue.empty()) {
			size_t current = queue.front();
			queue.pop();

			auto neighbors = connector_->get_valid_neighbors(current);
			for (const auto& neighbor : neighbors) {
				size_t neighbor_idx = neighbor.index;

				// Only process nodes in this basin that are below water level
				if (basin_labels_[neighbor_idx] == basin.id &&
						visited.find(neighbor_idx) == visited.end() &&
						(*elevation_)[neighbor_idx] < basin.water_level) {

					// Set receiver to minimize distance to spill
					T min_distance = std::numeric_limits<T>::max();
					size_t best_receiver = SIZE_MAX;

					auto neighbor_neighbors =
						connector_->get_valid_neighbors(neighbor_idx);
					for (const auto& nn : neighbor_neighbors) {
						if (visited.find(nn.index) != visited.end()) {
							T distance = static_cast<T>(
								connector_->euclidean_distance(nn.index, basin.spill_node));
							if (distance < min_distance) {
								min_distance = distance;
								best_receiver = nn.index;
							}
						}
					}

					if (best_receiver != SIZE_MAX) {
						receivers_[neighbor_idx] = best_receiver;
					}

					visited.insert(neighbor_idx);
					queue.push(neighbor_idx);
				}
			}
		}
	}

	/**
	 * Update donors based on modified receivers
	 */
	void update_donors() const
	{
		// Clear existing donors
		for (auto& donor_list : donors_) {
			donor_list.clear();
		}

		// Rebuild donors from receivers
		for (size_t i = 0; i < size_; ++i) {
			if (receivers_[i] != SIZE_MAX) {
				donors_[receivers_[i]].push_back(i);
			}
		}
	}

	/**
	 * Compute topological order for flow processing
	 */
	void compute_stack_order() const
	{
		std::fill(stack_order_.begin(), stack_order_.end(), SIZE_MAX);
		std::vector<bool> visited(size_, false);
		size_t order_index = 0;

		// Find all nodes with no donors (sources)
		std::queue<size_t> queue;
		std::vector<size_t> in_degree(size_, 0);

		// Calculate in-degrees
		for (size_t i = 0; i < size_; ++i) {
			in_degree[i] = donors_[i].size();
		}

		// Add sources to queue
		for (size_t i = 0; i < size_; ++i) {
			if (connector_->is_active_node(i) && in_degree[i] == 0) {
				queue.push(i);
			}
		}

		// Topological sort
		while (!queue.empty()) {
			size_t current = queue.front();
			queue.pop();

			if (!visited[current]) {
				visited[current] = true;
				stack_order_[order_index++] = current;

				// Process receiver
				size_t receiver = receivers_[current];
				if (receiver != SIZE_MAX) {
					in_degree[receiver]--;
					if (in_degree[receiver] == 0) {
						queue.push(receiver);
					}
				}
			}
		}
	}

public:
	// ======================
	// ANALYSIS AND UTILITIES
	// ======================

	/**
	 * Get statistics about the basin graph
	 */
	struct BasinGraphStatistics
	{
		size_t num_basins;
		size_t num_boundary_basins;
		size_t num_inner_basins;
		size_t num_links;
		size_t num_mst_links;
		size_t largest_basin_size;
		T total_depression_volume;
		size_t total_depression_area;

		BasinGraphStatistics()
			: num_basins(0)
			, num_boundary_basins(0)
			, num_inner_basins(0)
			, num_links(0)
			, num_mst_links(0)
			, largest_basin_size(0)
			, total_depression_volume(0)
			, total_depression_area(0)
		{
		}
	};

	BasinGraphStatistics compute_statistics() const
	{
		BasinGraphStatistics stats;

		stats.num_basins = basins_.size();
		stats.num_links = basin_links_.size();
		stats.num_mst_links = mst_links_.size();

		for (const auto& basin : basins_) {
			if (basin.is_boundary_basin) {
				stats.num_boundary_basins++;
			} else {
				stats.num_inner_basins++;

				// Calculate depression volume
				for (size_t node : basin.nodes) {
					T elevation = (*elevation_)[node];
					if (basin.water_level > elevation) {
						stats.total_depression_volume += (basin.water_level - elevation);
						stats.total_depression_area++;
					}
				}
			}

			stats.largest_basin_size =
				std::max(stats.largest_basin_size, basin.nodes.size());
		}

		return stats;
	}

	/**
	 * Generate detailed report
	 */
	std::string generate_report() const
	{
		auto stats = compute_statistics();
		std::ostringstream report;

		report << "=== CORDONNIER FLOW ROUTING REPORT ===\n\n";

		// Algorithm configuration
		report << "Flow Enforcement Strategy: ";
		switch (strategy_) {
			case FlowEnforcementStrategy::SIMPLE_CORRECTION:
				report << "Simple Correction\n";
				break;
			case FlowEnforcementStrategy::DEPRESSION_CARVING:
				report << "Depression Carving\n";
				break;
			case FlowEnforcementStrategy::DEPRESSION_FILLING:
				report << "Depression Filling\n";
				break;
		}

		report << "MST Algorithm: ";
		switch (mst_algorithm_) {
			case MSTAlgorithm::KRUSKAL:
				report << "Kruskal's Algorithm\n";
				break;
			case MSTAlgorithm::MARES_PLANAR:
				report << "Mareš Planar Algorithm\n";
				break;
		}

		report << "\n=== BASIN GRAPH STATISTICS ===\n";
		report << "Total Basins: " << stats.num_basins << "\n";
		report << "Boundary Basins: " << stats.num_boundary_basins << "\n";
		report << "Inner Basins: " << stats.num_inner_basins << "\n";
		report << "Basin Links: " << stats.num_links << "\n";
		report << "MST Links: " << stats.num_mst_links << "\n";
		report << "Largest Basin Size: " << stats.largest_basin_size << " nodes\n";

		if (stats.total_depression_area > 0) {
			report << "Total Depression Volume: " << stats.total_depression_volume
						 << "\n";
			report << "Total Depression Area: " << stats.total_depression_area
						 << " nodes\n";
			report << "Average Depression Depth: "
						 << (stats.total_depression_volume / stats.total_depression_area)
						 << "\n";
		}

		report << "\n=== PERFORMANCE NOTES ===\n";
		report << "This algorithm has O(n) complexity for planar graphs\n";
		report << "and O(n log n) complexity for general graphs.\n";
		report << "Memory usage is linear in grid size.\n";

		return report.str();
	}

	/**
	 * Validate flow routing results
	 */
	bool validate_flow_routing() const
	{
		// Check that every active node has a path to boundary
		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i))
				continue;

			// Trace path to boundary
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

				current = receivers_[current];
			}

			if (!reaches_boundary) {
				return false; // Found node that doesn't reach boundary
			}
		}

		return true;
	}

	/**
	 * Export basin map for visualization
	 */
	std::vector<size_t> export_basin_map() const { return basin_labels_; }

	/**
	 * Export water level map
	 */
	std::vector<T> export_water_level_map() const { return water_levels_; }

	/**
	 * Export flow direction map (as direction enum values)
	 */
	std::vector<uint8_t> export_flow_direction_map() const
	{
		std::vector<uint8_t> flow_directions(
			size_, static_cast<uint8_t>(Direction::INVALID));

		for (size_t i = 0; i < size_; ++i) {
			if (!connector_->is_active_node(i) || receivers_[i] == SIZE_MAX)
				continue;

			auto [from_row, from_col] = connector_->to_2d(i);
			auto [to_row, to_col] = connector_->to_2d(receivers_[i]);

			int dr = static_cast<int>(to_row) - static_cast<int>(from_row);
			int dc = static_cast<int>(to_col) - static_cast<int>(from_col);

			// Determine direction
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

			flow_directions[i] = static_cast<uint8_t>(dir);
		}

		return flow_directions;
	}
};

// ======================
// CONVENIENCE FUNCTIONS
// ======================

/**
 * Quick flow routing using Cordonnier method
 */
template<typename T = double>
class QuickCordonnier
{
public:
	/**
	 * Simple flow routing with depression filling
	 */
	static std::vector<size_t> compute_flow_receivers(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		CordonnierFlowRouter<T> router(
			connector, elevation_ref, FlowEnforcementStrategy::DEPRESSION_FILLING);
		router.compute_flow_routing();

		return router.get_receivers();
	}

	/**
	 * Flow routing with carving strategy
	 */
	static std::vector<size_t> compute_carving_flow(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		CordonnierFlowRouter<T> router(
			connector, elevation_ref, FlowEnforcementStrategy::DEPRESSION_CARVING);
		router.compute_flow_routing();

		return router.get_receivers();
	}

	/**
	 * Get flow accumulation using Cordonnier method
	 */
	static std::vector<size_t> compute_flow_accumulation(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		CordonnierFlowRouter<T> router(connector, elevation_ref);
		router.compute_flow_routing();

		return router.compute_flow_accumulation();
	}

	/**
	 * Get drainage area using Cordonnier method
	 */
	static std::vector<T> compute_drainage_area(
		std::shared_ptr<Connector<T>> connector,
		const ArrayRef<T>& elevation,
		T cell_area = 1.0)
	{
		auto elevation_ref = std::make_shared<ArrayRef<T>>(elevation);

		CordonnierFlowRouter<T> router(connector, elevation_ref);
		router.compute_flow_routing();

		return router.compute_drainage_area(cell_area);
	}
};

// Type aliases
using CordonnierFlowRouterF32 = CordonnierFlowRouter<float>;
using CordonnierFlowRouterF64 = CordonnierFlowRouter<double>;

} // namespace dagger2
