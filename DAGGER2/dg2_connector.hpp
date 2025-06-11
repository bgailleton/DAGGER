#pragma once

#include "dg2_BCs.hpp"
#include "dg2_array.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <random>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * Direction enumeration for D4/D8 connectivity
 */
enum class Direction : uint8_t
{
	// D4 directions (cardinal)
	NORTH = 0, // Up
	EAST = 1,	 // Right
	SOUTH = 2, // Down
	WEST = 3,	 // Left

	// D8 additional directions (diagonal)
	NORTHEAST = 4, // Up-Right
	SOUTHEAST = 5, // Down-Right
	SOUTHWEST = 6, // Down-Left
	NORTHWEST = 7, // Up-Left

	// Special values
	CENTER = 8,		// Current node (self)
	INVALID = 255 // Invalid direction
};

/**
 * Connectivity type
 */
enum class ConnectivityType : uint8_t
{
	D4 = 4, // 4-connected (cardinal directions only)
	D8 = 8	// 8-connected (cardinal + diagonal)
};

/**
 * Neighbor information structure
 */
struct Neighbor
{
	size_t index;				 // 1D index in grid
	size_t row, col;		 // 2D coordinates
	Direction direction; // Direction from center
	double distance; // Distance from center (1.0 for D4, sqrt(2) for diagonals)
	NodeType boundary_type; // Boundary condition of this neighbor
	bool is_valid;					// Whether this neighbor exists and is active

	Neighbor()
		: index(0)
		, row(0)
		, col(0)
		, direction(Direction::INVALID)
		, distance(0.0)
		, boundary_type(NodeType::NO_DATA)
		, is_valid(false)
	{
	}

	Neighbor(size_t idx,
					 size_t r,
					 size_t c,
					 Direction dir,
					 double dist,
					 NodeType bt,
					 bool valid)
		: index(idx)
		, row(r)
		, col(c)
		, direction(dir)
		, distance(dist)
		, boundary_type(bt)
		, is_valid(valid)
	{
	}
};

/**
 * Connector: High-performance neighbor management for grid-based computations
 *
 * This class provides comprehensive neighboring operations for D4/D8
 * connectivity with full boundary condition support. It's designed to be the
 * central hub for all spatial operations in landscape modeling.
 */
template<typename T = double>
class Connector
{
private:
	// Grid dimensions
	size_t rows_;
	size_t cols_;
	size_t size_;

	// Connectivity configuration
	ConnectivityType connectivity_;
	size_t num_directions_;

	// Boundary conditions grid
	std::shared_ptr<Grid2D<NodeType>> boundary_grid_;

	// Direction vectors for efficient neighbor calculation
	// Order: N, E, S, W, NE, SE, SW, NW
	static constexpr std::array<int, 8> dr_ = { -1, 0, 1, 0, -1, 1, 1, -1 };
	static constexpr std::array<int, 8> dc_ = { 0, 1, 0, -1, 1, 1, -1, -1 };

	// Distances for each direction
	static constexpr std::array<double, 8> distances_ = {
		1.0,
		1.0,
		1.0,
		1.0,								// D4 distances
		1.4142135623730951, // sqrt(2) for diagonals
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951
	};

	// Opposite directions lookup
	static constexpr std::array<Direction, 8> opposite_directions_ = {
		Direction::SOUTH,			// Opposite of NORTH
		Direction::WEST,			// Opposite of EAST
		Direction::NORTH,			// Opposite of SOUTH
		Direction::EAST,			// Opposite of WEST
		Direction::SOUTHWEST, // Opposite of NORTHEAST
		Direction::NORTHWEST, // Opposite of SOUTHEAST
		Direction::NORTHEAST, // Opposite of SOUTHWEST
		Direction::SOUTHEAST	// Opposite of NORTHWEST
	};

public:
	/**
	 * Constructor
	 */
	Connector(size_t rows,
						size_t cols,
						ConnectivityType connectivity = ConnectivityType::D8)
		: rows_(rows)
		, cols_(cols)
		, size_(rows * cols)
		, connectivity_(connectivity)
		, num_directions_(static_cast<size_t>(connectivity))
	{

		if (rows == 0 || cols == 0) {
			throw std::invalid_argument("Grid dimensions must be positive");
		}

		// Initialize boundary grid with all NORMAL nodes
		std::vector<NodeType> bc_data(size_, NodeType::NORMAL);
		boundary_grid_ = std::make_shared<Grid2D<NodeType>>(bc_data, rows_, cols_);
	}

	/**
	 * Constructor with boundary conditions
	 */
	Connector(size_t rows,
						size_t cols,
						std::shared_ptr<Grid2D<NodeType>> boundary_grid,
						ConnectivityType connectivity = ConnectivityType::D8)
		: rows_(rows)
		, cols_(cols)
		, size_(rows * cols)
		, connectivity_(connectivity)
		, num_directions_(static_cast<size_t>(connectivity))
		, boundary_grid_(boundary_grid)
	{

		if (rows == 0 || cols == 0) {
			throw std::invalid_argument("Grid dimensions must be positive");
		}

		if (!boundary_grid_ || boundary_grid_->rows() != rows_ ||
				boundary_grid_->cols() != cols_) {
			throw std::invalid_argument(
				"Boundary grid dimensions must match connector dimensions");
		}
	}

	// ======================
	// BASIC PROPERTIES
	// ======================

	size_t rows() const { return rows_; }
	size_t cols() const { return cols_; }
	size_t size() const { return size_; }
	ConnectivityType connectivity_type() const { return connectivity_; }
	size_t num_directions() const { return num_directions_; }

	// ======================
	// INDEX CONVERSION
	// ======================

	inline size_t to_1d(size_t row, size_t col) const
	{
		return row * cols_ + col;
	}

	inline std::pair<size_t, size_t> to_2d(size_t index) const
	{
		return { index / cols_, index % cols_ };
	}

	inline bool is_valid_coord(size_t row, size_t col) const
	{
		return row < rows_ && col < cols_;
	}

	inline bool is_valid_coord(int row, int col) const
	{
		return row >= 0 && col >= 0 && static_cast<size_t>(row) < rows_ &&
					 static_cast<size_t>(col) < cols_;
	}

	inline bool is_valid_index(size_t index) const { return index < size_; }

	// ======================
	// BOUNDARY CONDITIONS
	// ======================

	NodeType get_boundary_type(size_t row, size_t col) const
	{
		return (*boundary_grid_)(row, col);
	}

	NodeType get_boundary_type(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_boundary_type(row, col);
	}

	void set_boundary_type(size_t row, size_t col, NodeType type)
	{
		(*boundary_grid_)(row, col) = type;
	}

	void set_boundary_type(size_t index, NodeType type)
	{
		auto [row, col] = to_2d(index);
		set_boundary_type(row, col, type);
	}

	bool is_active_node(size_t row, size_t col) const
	{
		return NodeTypeUtils::is_active(get_boundary_type(row, col));
	}

	bool is_active_node(size_t index) const
	{
		return NodeTypeUtils::is_active(get_boundary_type(index));
	}

	bool is_boundary_node(size_t row, size_t col) const
	{
		return NodeTypeUtils::is_boundary(get_boundary_type(row, col));
	}

	bool is_boundary_node(size_t index) const
	{
		return NodeTypeUtils::is_boundary(get_boundary_type(index));
	}

	// ======================
	// DIRECTION UTILITIES
	// ======================

	Direction get_opposite_direction(Direction dir) const
	{
		uint8_t idx = static_cast<uint8_t>(dir);
		if (idx >= 8)
			return Direction::INVALID;
		return opposite_directions_[idx];
	}

	double get_direction_distance(Direction dir) const
	{
		uint8_t idx = static_cast<uint8_t>(dir);
		if (idx >= 8)
			return 0.0;
		return distances_[idx];
	}

	bool is_diagonal_direction(Direction dir) const
	{
		return dir >= Direction::NORTHEAST && dir <= Direction::NORTHWEST;
	}

	bool is_cardinal_direction(Direction dir) const
	{
		return dir >= Direction::NORTH && dir <= Direction::WEST;
	}

	// ======================
	// SINGLE NEIGHBOR ACCESS (UNCHECKED - FASTEST)
	// ======================

	inline size_t get_neighbor_index_unchecked(size_t index, Direction dir) const
	{
		auto [row, col] = to_2d(index);
		uint8_t d = static_cast<uint8_t>(dir);
		return to_1d(row + dr_[d], col + dc_[d]);
	}

	inline std::pair<size_t, size_t>
	get_neighbor_coord_unchecked(size_t row, size_t col, Direction dir) const
	{
		uint8_t d = static_cast<uint8_t>(dir);
		return { row + dr_[d], col + dc_[d] };
	}

	// Specific direction accessors (unchecked) - D8 only for diagonals
	inline size_t north_unchecked(size_t index) const
	{
		return get_neighbor_index_unchecked(index, Direction::NORTH);
	}
	inline size_t east_unchecked(size_t index) const
	{
		return get_neighbor_index_unchecked(index, Direction::EAST);
	}
	inline size_t south_unchecked(size_t index) const
	{
		return get_neighbor_index_unchecked(index, Direction::SOUTH);
	}
	inline size_t west_unchecked(size_t index) const
	{
		return get_neighbor_index_unchecked(index, Direction::WEST);
	}
	inline size_t northeast_unchecked(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor_index_unchecked(index, Direction::NORTHEAST)
						 : SIZE_MAX;
	}
	inline size_t southeast_unchecked(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor_index_unchecked(index, Direction::SOUTHEAST)
						 : SIZE_MAX;
	}
	inline size_t southwest_unchecked(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor_index_unchecked(index, Direction::SOUTHWEST)
						 : SIZE_MAX;
	}
	inline size_t northwest_unchecked(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor_index_unchecked(index, Direction::NORTHWEST)
						 : SIZE_MAX;
	}

	// ======================
	// SINGLE NEIGHBOR ACCESS (CHECKED)
	// ======================

	Neighbor get_neighbor(size_t row, size_t col, Direction dir) const
	{
		uint8_t d = static_cast<uint8_t>(dir);
		if (d >= num_directions_) {
			return Neighbor(); // Invalid direction for this connectivity
		}

		int new_row = static_cast<int>(row) + dr_[d];
		int new_col = static_cast<int>(col) + dc_[d];

		if (!is_valid_coord(new_row, new_col)) {
			return Neighbor(); // Out of bounds
		}

		size_t nr = static_cast<size_t>(new_row);
		size_t nc = static_cast<size_t>(new_col);
		size_t nidx = to_1d(nr, nc);
		NodeType bt = get_boundary_type(nr, nc);
		bool valid = NodeTypeUtils::is_active(bt);

		return Neighbor(nidx, nr, nc, dir, distances_[d], bt, valid);
	}

	Neighbor get_neighbor(size_t index, Direction dir) const
	{
		auto [row, col] = to_2d(index);
		return get_neighbor(row, col, dir);
	}

	// Specific direction accessors (checked) - D8 only for diagonals
	Neighbor north(size_t index) const
	{
		return get_neighbor(index, Direction::NORTH);
	}
	Neighbor east(size_t index) const
	{
		return get_neighbor(index, Direction::EAST);
	}
	Neighbor south(size_t index) const
	{
		return get_neighbor(index, Direction::SOUTH);
	}
	Neighbor west(size_t index) const
	{
		return get_neighbor(index, Direction::WEST);
	}
	Neighbor northeast(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(index, Direction::NORTHEAST)
						 : Neighbor();
	}
	Neighbor southeast(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(index, Direction::SOUTHEAST)
						 : Neighbor();
	}
	Neighbor southwest(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(index, Direction::SOUTHWEST)
						 : Neighbor();
	}
	Neighbor northwest(size_t index) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(index, Direction::NORTHWEST)
						 : Neighbor();
	}

	Neighbor north(size_t row, size_t col) const
	{
		return get_neighbor(row, col, Direction::NORTH);
	}
	Neighbor east(size_t row, size_t col) const
	{
		return get_neighbor(row, col, Direction::EAST);
	}
	Neighbor south(size_t row, size_t col) const
	{
		return get_neighbor(row, col, Direction::SOUTH);
	}
	Neighbor west(size_t row, size_t col) const
	{
		return get_neighbor(row, col, Direction::WEST);
	}
	Neighbor northeast(size_t row, size_t col) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(row, col, Direction::NORTHEAST)
						 : Neighbor();
	}
	Neighbor southeast(size_t row, size_t col) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(row, col, Direction::SOUTHEAST)
						 : Neighbor();
	}
	Neighbor southwest(size_t row, size_t col) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(row, col, Direction::SOUTHWEST)
						 : Neighbor();
	}
	Neighbor northwest(size_t row, size_t col) const
	{
		return (connectivity_ == ConnectivityType::D8)
						 ? get_neighbor(row, col, Direction::NORTHWEST)
						 : Neighbor();
	}
	// ALL NEIGHBORS ACCESS
	// ======================

	std::vector<Neighbor> get_all_neighbors(size_t row, size_t col) const
	{
		std::vector<Neighbor> neighbors;
		neighbors.reserve(num_directions_);

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);
			neighbors.push_back(n);
		}

		return neighbors;
	}

	std::vector<Neighbor> get_all_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_all_neighbors(row, col);
	}

	std::vector<Neighbor> get_valid_neighbors(size_t row, size_t col) const
	{
		std::vector<Neighbor> neighbors;
		neighbors.reserve(num_directions_);

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);
			if (n.is_valid) {
				neighbors.push_back(n);
			}
		}

		return neighbors;
	}

	std::vector<Neighbor> get_valid_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_valid_neighbors(row, col);
	}

	// ======================
	// VECTORIZED OPERATIONS
	// ======================

	std::vector<size_t> get_neighbor_indices(size_t index) const
	{
		std::vector<size_t> indices;
		indices.reserve(num_directions_);

		auto [row, col] = to_2d(index);
		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);
			indices.push_back(n.is_valid ? n.index
																	 : SIZE_MAX); // SIZE_MAX for invalid
		}

		return indices;
	}

	std::vector<size_t> get_valid_neighbor_indices(size_t index) const
	{
		std::vector<size_t> indices;
		indices.reserve(num_directions_);

		auto [row, col] = to_2d(index);
		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);
			if (n.is_valid) {
				indices.push_back(n.index);
			}
		}

		return indices;
	}

	// ======================
	// BULK OPERATIONS
	// ======================

	void get_all_neighbors_bulk(const std::vector<size_t>& indices,
															std::vector<std::vector<Neighbor>>& results) const
	{
		results.clear();
		results.reserve(indices.size());

		for (size_t idx : indices) {
			results.push_back(get_all_neighbors(idx));
		}
	}

	void get_valid_neighbors_bulk(
		const std::vector<size_t>& indices,
		std::vector<std::vector<Neighbor>>& results) const
	{
		results.clear();
		results.reserve(indices.size());

		for (size_t idx : indices) {
			results.push_back(get_valid_neighbors(idx));
		}
	}

	// ======================
	// FLOW DIRECTION HELPERS
	// ======================

	Direction get_steepest_descent_direction(size_t index,
																					 const ArrayRef<T>& elevation) const
	{
		auto [row, col] = to_2d(index);
		T center_elev = elevation[index];

		Direction steepest_dir = Direction::INVALID;
		T max_gradient = 0.0;

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);

			if (n.is_valid) {
				T neighbor_elev = elevation[n.index];
				T gradient = (center_elev - neighbor_elev) / n.distance;

				if (gradient > max_gradient) {
					max_gradient = gradient;
					steepest_dir = dir;
				}
			}
		}

		return steepest_dir;
	}

	std::vector<Direction> get_downhill_directions(
		size_t index,
		const ArrayRef<T>& elevation) const
	{
		auto [row, col] = to_2d(index);
		T center_elev = elevation[index];

		std::vector<Direction> downhill_dirs;
		downhill_dirs.reserve(num_directions_);

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);

			if (n.is_valid && elevation[n.index] < center_elev) {
				downhill_dirs.push_back(dir);
			}
		}

		return downhill_dirs;
	}

	// ======================
	// DISTANCE CALCULATIONS
	// ======================

	double euclidean_distance(size_t idx1, size_t idx2) const
	{
		auto [r1, c1] = to_2d(idx1);
		auto [r2, c2] = to_2d(idx2);

		double dr = static_cast<double>(r2) - static_cast<double>(r1);
		double dc = static_cast<double>(c2) - static_cast<double>(c1);

		return std::sqrt(dr * dr + dc * dc);
	}

	double manhattan_distance(size_t idx1, size_t idx2) const
	{
		auto [r1, c1] = to_2d(idx1);
		auto [r2, c2] = to_2d(idx2);

		return std::abs(static_cast<int>(r2) - static_cast<int>(r1)) +
					 std::abs(static_cast<int>(c2) - static_cast<int>(c1));
	}

	// ======================
	// SPECIAL BOUNDARY CONDITION HANDLING
	// ======================

	/**
	 * Get the actual neighbor after applying boundary condition logic
	 * This handles REFLECT and PERIODIC boundary conditions properly
	 * based on the individual node's NodeType
	 */
	Neighbor get_effective_neighbor(size_t row, size_t col, Direction dir) const
	{
		uint8_t d = static_cast<uint8_t>(dir);
		int new_row = static_cast<int>(row) + dr_[d];
		int new_col = static_cast<int>(col) + dc_[d];

		// Check if we're going out of bounds
		bool out_of_bounds = !is_valid_coord(new_row, new_col);

		if (!out_of_bounds) {
			// Normal in-bounds neighbor - check if it's active
			size_t nr = static_cast<size_t>(new_row);
			size_t nc = static_cast<size_t>(new_col);
			size_t nidx = to_1d(nr, nc);
			NodeType neighbor_bc = get_boundary_type(nr, nc);
			bool is_active = NodeTypeUtils::is_active(neighbor_bc);

			return Neighbor(nidx, nr, nc, dir, distances_[d], neighbor_bc, is_active);
		}

		// We went out of bounds - check current node's boundary condition
		NodeType current_bc = get_boundary_type(row, col);

		// Handle different boundary conditions based on current node type
		switch (current_bc) {
			case NodeType::PERIODIC: {
				// Wrap around to opposite side
				auto [wrapped_row, wrapped_col] = apply_periodic_bc(new_row, new_col);

				size_t wrapped_idx = to_1d(wrapped_row, wrapped_col);
				NodeType wrapped_bc = get_boundary_type(wrapped_row, wrapped_col);
				bool is_active = NodeTypeUtils::is_active(wrapped_bc);

				return Neighbor(wrapped_idx,
												wrapped_row,
												wrapped_col,
												dir,
												distances_[d],
												wrapped_bc,
												is_active);
			}

			case NodeType::REFLECT: {
				// Reflect back to current node with opposite direction
				return Neighbor(to_1d(row, col),
												row,
												col,
												get_opposite_direction(dir),
												0.0,
												current_bc,
												true);
			}

			case NodeType::HAS_TO_OUT:
			case NodeType::CAN_OUT: {
				// Flow exits the domain - return invalid neighbor
				// The calling code should handle this as flux leaving the system
				return Neighbor();
			}

			case NodeType::IN: {
				// Input nodes don't allow outflow - return invalid neighbor
				return Neighbor();
			}

			default:
				// For NORMAL and NO_DATA nodes at boundary, return invalid
				return Neighbor();
		}
	}

	/**
	 * Get effective neighbors with boundary condition handling
	 */
	std::vector<Neighbor> get_effective_neighbors(size_t row, size_t col) const
	{
		std::vector<Neighbor> neighbors;
		neighbors.reserve(num_directions_);

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_effective_neighbor(row, col, dir);
			neighbors.push_back(n);
		}

		return neighbors;
	}

	std::vector<Neighbor> get_effective_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_effective_neighbors(row, col);
	}

	/**
	 * Get effective valid neighbors (only those that are active)
	 */
	std::vector<Neighbor> get_effective_valid_neighbors(size_t row,
																											size_t col) const
	{
		std::vector<Neighbor> neighbors;
		neighbors.reserve(num_directions_);

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_effective_neighbor(row, col, dir);
			if (n.is_valid) {
				neighbors.push_back(n);
			}
		}

		return neighbors;
	}

	std::vector<Neighbor> get_effective_valid_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_effective_valid_neighbors(row, col);
	}

	/**
	 * Check if a flow from current node in given direction would be handled
	 * by boundary conditions (PERIODIC/REFLECT) rather than going out of domain
	 */
	bool has_boundary_handling(size_t row, size_t col, Direction dir) const
	{
		// Check if we would go out of bounds
		uint8_t d = static_cast<uint8_t>(dir);
		int new_row = static_cast<int>(row) + dr_[d];
		int new_col = static_cast<int>(col) + dc_[d];

		if (is_valid_coord(new_row, new_col)) {
			return false; // Normal in-bounds neighbor
		}

		// Out of bounds - check current node's boundary condition
		NodeType bc = get_boundary_type(row, col);
		return bc == NodeType::PERIODIC || bc == NodeType::REFLECT;
	}

	// ======================
	// ADVANCED FLOW HANDLING
	// ======================

	/**
	 * Transfer flow/flux considering boundary conditions
	 * Returns the actual destination index and modified flux
	 */
	struct FlowTransfer
	{
		size_t destination_index;
		T flux_multiplier; // For reflection (can be negative)
		bool is_periodic_wrap;
		bool is_reflection;
		bool exits_domain;

		FlowTransfer()
			: destination_index(SIZE_MAX)
			, flux_multiplier(0.0)
			, is_periodic_wrap(false)
			, is_reflection(false)
			, exits_domain(true)
		{
		}

		FlowTransfer(size_t dest,
								 T mult,
								 bool periodic = false,
								 bool reflect = false,
								 bool exits = false)
			: destination_index(dest)
			, flux_multiplier(mult)
			, is_periodic_wrap(periodic)
			, is_reflection(reflect)
			, exits_domain(exits)
		{
		}
	};

	FlowTransfer transfer_flow(size_t from_index,
														 Direction direction,
														 T flux) const
	{
		auto [row, col] = to_2d(from_index);
		NodeType from_bc = get_boundary_type(row, col);

		uint8_t d = static_cast<uint8_t>(direction);
		int target_row = static_cast<int>(row) + dr_[d];
		int target_col = static_cast<int>(col) + dc_[d];
		bool going_out_of_bounds = !is_valid_coord(target_row, target_col);

		if (!going_out_of_bounds) {
			// Normal in-bounds transfer
			size_t target_idx =
				to_1d(static_cast<size_t>(target_row), static_cast<size_t>(target_col));
			NodeType target_bc = get_boundary_type(static_cast<size_t>(target_row),
																						 static_cast<size_t>(target_col));

			if (NodeTypeUtils::is_active(target_bc)) {
				return FlowTransfer(
					target_idx, static_cast<T>(1.0), false, false, false);
			} else {
				// Target is inactive (NO_DATA) - flow is lost
				return FlowTransfer();
			}
		}

		// Going out of bounds - handle based on source node's boundary condition
		switch (from_bc) {
			case NodeType::PERIODIC: {
				// Wrap around to opposite side
				auto [wrapped_row, wrapped_col] =
					apply_periodic_bc(target_row, target_col);
				size_t wrapped_idx = to_1d(wrapped_row, wrapped_col);
				NodeType wrapped_bc = get_boundary_type(wrapped_row, wrapped_col);

				if (NodeTypeUtils::is_active(wrapped_bc)) {
					return FlowTransfer(
						wrapped_idx, static_cast<T>(1.0), true, false, false);
				} else {
					// Wrapped to inactive node - flow is lost
					return FlowTransfer();
				}
			}

			case NodeType::REFLECT: {
				// Flow bounces back to source node
				return FlowTransfer(
					from_index, static_cast<T>(-1.0), false, true, false);
			}

			case NodeType::HAS_TO_OUT: {
				// Flow must exit domain
				return FlowTransfer();
			}

			case NodeType::CAN_OUT: {
				// Flow can exit domain (decision depends on algorithm)
				return FlowTransfer();
			}

			case NodeType::IN: {
				// Input nodes don't allow outflow
				return FlowTransfer();
			}

			default: {
				// NORMAL or NO_DATA nodes at boundary - flow exits
				return FlowTransfer();
			}
		}
	}

	/**
	 * Apply periodic boundary conditions to coordinates
	 */
	std::pair<size_t, size_t> apply_periodic_bc(int row, int col) const
	{
		size_t new_row, new_col;

		// Handle row wrapping
		if (row < 0) {
			new_row = rows_ - 1; // Wrap to bottom
		} else if (row >= static_cast<int>(rows_)) {
			new_row = 0; // Wrap to top
		} else {
			new_row = static_cast<size_t>(row);
		}

		// Handle column wrapping
		if (col < 0) {
			new_col = cols_ - 1; // Wrap to right
		} else if (col >= static_cast<int>(cols_)) {
			new_col = 0; // Wrap to left
		} else {
			new_col = static_cast<size_t>(col);
		}

		return { new_row, new_col };
	}

	// ======================
	// STANDARD BOUNDARY HANDLING
	// ======================

	std::vector<size_t> get_boundary_nodes(NodeType boundary_type) const
	{
		std::vector<size_t> boundary_nodes;

		for (size_t i = 0; i < size_; ++i) {
			if (get_boundary_type(i) == boundary_type) {
				boundary_nodes.push_back(i);
			}
		}

		return boundary_nodes;
	}

	std::vector<size_t> get_all_boundary_nodes() const
	{
		std::vector<size_t> boundary_nodes;

		for (size_t i = 0; i < size_; ++i) {
			if (is_boundary_node(i)) {
				boundary_nodes.push_back(i);
			}
		}

		return boundary_nodes;
	}

	void set_border_boundary(NodeType boundary_type)
	{
		// Set top and bottom rows
		for (size_t c = 0; c < cols_; ++c) {
			set_boundary_type(0, c, boundary_type);					// Top row
			set_boundary_type(rows_ - 1, c, boundary_type); // Bottom row
		}

		// Set left and right columns
		for (size_t r = 0; r < rows_; ++r) {
			set_boundary_type(r, 0, boundary_type);					// Left column
			set_boundary_type(r, cols_ - 1, boundary_type); // Right column
		}
	}

	/**
	 * Set up periodic boundary conditions on opposite borders
	 */
	void set_periodic_boundaries() { set_border_boundary(NodeType::PERIODIC); }

	/**
	 * Set up reflective boundary conditions
	 */
	void set_reflective_boundaries() { set_border_boundary(NodeType::REFLECT); }

	// ======================
	// PERFORMANCE UTILITIES
	// ======================

	size_t count_valid_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		size_t count = 0;

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);
			if (n.is_valid)
				++count;
		}

		return count;
	}

	bool has_valid_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);

		for (size_t d = 0; d < num_directions_; ++d) {
			Direction dir = static_cast<Direction>(d);
			Neighbor n = get_neighbor(row, col, dir);
			if (n.is_valid)
				return true;
		}

		return false;
	}

	// ======================
	// ADVANCED UTILITIES & CREATIVE TOOLS
	// ======================

	/**
	 * Random number generator for consistent seeding
	 */
private:
	mutable std::mt19937 rng_{ std::random_device{}() };

public:
	void set_random_seed(uint32_t seed) const { rng_.seed(seed); }

	// ======================
	// RANDOM ACCESS UTILITIES
	// ======================

	/**
	 * Get random valid index in grid
	 */
	size_t random_index() const
	{
		std::uniform_int_distribution<size_t> dist(0, size_ - 1);
		return dist(rng_);
	}

	/**
	 * Get random valid 2D coordinate
	 */
	std::pair<size_t, size_t> random_coord() const
	{
		std::uniform_int_distribution<size_t> row_dist(0, rows_ - 1);
		std::uniform_int_distribution<size_t> col_dist(0, cols_ - 1);
		return { row_dist(rng_), col_dist(rng_) };
	}

	/**
	 * Get random active node (not NO_DATA)
	 */
	size_t random_active_index() const
	{
		size_t attempts = 0;
		const size_t max_attempts = size_ * 2; // Prevent infinite loop

		while (attempts < max_attempts) {
			size_t idx = random_index();
			if (is_active_node(idx)) {
				return idx;
			}
			++attempts;
		}

		// Fallback: find any active node
		for (size_t i = 0; i < size_; ++i) {
			if (is_active_node(i))
				return i;
		}

		throw std::runtime_error("No active nodes found in grid");
	}

	/**
	 * Get random node of specific boundary type
	 */
	size_t random_boundary_node(NodeType boundary_type) const
	{
		auto boundary_nodes = get_boundary_nodes(boundary_type);
		if (boundary_nodes.empty()) {
			throw std::runtime_error("No nodes of specified boundary type found");
		}

		std::uniform_int_distribution<size_t> dist(0, boundary_nodes.size() - 1);
		return boundary_nodes[dist(rng_)];
	}

	/**
	 * Get random neighbor of a given node
	 */
	Neighbor random_neighbor(size_t index) const
	{
		auto neighbors = get_valid_neighbors(index);
		if (neighbors.empty()) {
			return Neighbor(); // No valid neighbors
		}

		std::uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);
		return neighbors[dist(rng_)];
	}

	/**
	 * Get random direction for a node (respects connectivity)
	 */
	Direction random_direction() const
	{
		std::uniform_int_distribution<size_t> dist(0, num_directions_ - 1);
		return static_cast<Direction>(dist(rng_));
	}

	// ======================
	// TRAVERSAL UTILITIES
	// ======================

	/**
	 * Generate random traversal order (visit all nodes once in random order)
	 */
	std::vector<size_t> random_traversal() const
	{
		std::vector<size_t> indices(size_);
		std::iota(indices.begin(), indices.end(), 0);
		std::shuffle(indices.begin(), indices.end(), rng_);
		return indices;
	}

	/**
	 * Generate random traversal of only active nodes
	 */
	std::vector<size_t> random_active_traversal() const
	{
		std::vector<size_t> active_indices;
		active_indices.reserve(size_);

		for (size_t i = 0; i < size_; ++i) {
			if (is_active_node(i)) {
				active_indices.push_back(i);
			}
		}

		std::shuffle(active_indices.begin(), active_indices.end(), rng_);
		return active_indices;
	}

	/**
	 * Generate spiral traversal from center outward
	 */
	std::vector<size_t> spiral_traversal(size_t center_row,
																			 size_t center_col) const
	{
		std::vector<size_t> spiral_order;
		spiral_order.reserve(size_);

		std::vector<std::vector<bool>> visited(rows_,
																					 std::vector<bool>(cols_, false));

		// Directions for spiral: right, down, left, up
		std::array<int, 4> dr = { 0, 1, 0, -1 };
		std::array<int, 4> dc = { 1, 0, -1, 0 };

		size_t r = center_row, c = center_col;
		int dir = 0;
		int steps = 1;

		spiral_order.push_back(to_1d(r, c));
		visited[r][c] = true;

		while (spiral_order.size() < size_) {
			for (int i = 0; i < 2; ++i) { // Two sides per step count
				for (int j = 0; j < steps; ++j) {
					int new_r = static_cast<int>(r) + dr[dir];
					int new_c = static_cast<int>(c) + dc[dir];

					if (is_valid_coord(new_r, new_c) && !visited[new_r][new_c]) {
						r = static_cast<size_t>(new_r);
						c = static_cast<size_t>(new_c);
						spiral_order.push_back(to_1d(r, c));
						visited[r][c] = true;
					}
				}
				dir = (dir + 1) % 4;
			}
			++steps;
		}

		return spiral_order;
	}

	/**
	 * Generate spiral traversal from grid center
	 */
	std::vector<size_t> spiral_traversal() const
	{
		return spiral_traversal(rows_ / 2, cols_ / 2);
	}

	/**
	 * Breadth-first traversal from starting point
	 */
	std::vector<size_t> bfs_traversal(size_t start_index) const
	{
		std::vector<size_t> order;
		order.reserve(size_);

		std::vector<bool> visited(size_, false);
		std::queue<size_t> queue;

		queue.push(start_index);
		visited[start_index] = true;

		while (!queue.empty()) {
			size_t current = queue.front();
			queue.pop();
			order.push_back(current);

			auto neighbors = get_valid_neighbors(current);
			for (const auto& neighbor : neighbors) {
				if (!visited[neighbor.index]) {
					visited[neighbor.index] = true;
					queue.push(neighbor.index);
				}
			}
		}

		return order;
	}

	/**
	 * Depth-first traversal from starting point
	 */
	std::vector<size_t> dfs_traversal(size_t start_index) const
	{
		std::vector<size_t> order;
		order.reserve(size_);

		std::vector<bool> visited(size_, false);
		std::stack<size_t> stack;

		stack.push(start_index);

		while (!stack.empty()) {
			size_t current = stack.top();
			stack.pop();

			if (!visited[current]) {
				visited[current] = true;
				order.push_back(current);

				auto neighbors = get_valid_neighbors(current);
				// Reverse order for consistent left-to-right traversal
				for (auto it = neighbors.rbegin(); it != neighbors.rend(); ++it) {
					if (!visited[it->index]) {
						stack.push(it->index);
					}
				}
			}
		}

		return order;
	}

	// ======================
	// NEIGHBORHOOD ANALYSIS
	// ======================

	/**
	 * Calculate local density (number of active neighbors / total possible
	 * neighbors)
	 */
	double local_density(size_t index) const
	{
		size_t valid_count = count_valid_neighbors(index);
		return static_cast<double>(valid_count) /
					 static_cast<double>(num_directions_);
	}

	/**
	 * Get neighborhood statistics
	 */
	struct NeighborhoodStats
	{
		size_t total_neighbors;
		size_t valid_neighbors;
		size_t boundary_neighbors;
		double density;
		double avg_distance;
		std::unordered_map<NodeType, size_t> boundary_type_counts;

		NeighborhoodStats()
			: total_neighbors(0)
			, valid_neighbors(0)
			, boundary_neighbors(0)
			, density(0.0)
			, avg_distance(0.0)
		{
		}
	};

	NeighborhoodStats analyze_neighborhood(size_t index) const
	{
		NeighborhoodStats stats;
		auto neighbors = get_all_neighbors(index);

		stats.total_neighbors = neighbors.size();
		double total_distance = 0.0;

		for (const auto& neighbor : neighbors) {
			if (neighbor.is_valid) {
				++stats.valid_neighbors;
				total_distance += neighbor.distance;

				if (is_boundary_node(neighbor.index)) {
					++stats.boundary_neighbors;
				}

				++stats.boundary_type_counts[neighbor.boundary_type];
			}
		}

		stats.density = static_cast<double>(stats.valid_neighbors) /
										static_cast<double>(stats.total_neighbors);
		stats.avg_distance = (stats.valid_neighbors > 0)
													 ? total_distance / stats.valid_neighbors
													 : 0.0;

		return stats;
	}

	// ======================
	// PATHFINDING & CONNECTIVITY
	// ======================

	/**
	 * Find shortest path between two nodes (A* algorithm)
	 */
	std::vector<size_t> find_path(size_t start, size_t goal) const
	{
		if (start == goal)
			return { start };

		struct Node
		{
			size_t index;
			double g_cost; // Distance from start
			double h_cost; // Heuristic distance to goal
			double f_cost() const { return g_cost + h_cost; }
			size_t parent;

			Node(size_t idx, double g, double h, size_t p)
				: index(idx)
				, g_cost(g)
				, h_cost(h)
				, parent(p)
			{
			}
		};

		auto compare = [](const Node& a, const Node& b) {
			return a.f_cost() > b.f_cost();
		};

		std::priority_queue<Node, std::vector<Node>, decltype(compare)> open_set(
			compare);
		std::unordered_set<size_t> closed_set;
		std::unordered_map<size_t, Node> nodes;

		double h_start = euclidean_distance(start, goal);
		open_set.emplace(start, 0.0, h_start, SIZE_MAX);
		nodes.emplace(start, Node(start, 0.0, h_start, SIZE_MAX));

		while (!open_set.empty()) {
			Node current = open_set.top();
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

			auto neighbors = get_valid_neighbors(current.index);
			for (const auto& neighbor : neighbors) {
				if (closed_set.count(neighbor.index))
					continue;

				double tentative_g = current.g_cost + neighbor.distance;
				double h = euclidean_distance(neighbor.index, goal);

				auto it = nodes.find(neighbor.index);
				if (it == nodes.end() || tentative_g < it->second.g_cost) {
					nodes[neighbor.index] =
						Node(neighbor.index, tentative_g, h, current.index);
					open_set.emplace(neighbor.index, tentative_g, h, current.index);
				}
			}
		}

		return {}; // No path found
	}

	/**
	 * Check if two nodes are connected
	 */
	bool are_connected(size_t start, size_t goal) const
	{
		return !find_path(start, goal).empty();
	}

	/**
	 * Find all connected components
	 */
	std::vector<std::vector<size_t>> find_connected_components() const
	{
		std::vector<std::vector<size_t>> components;
		std::vector<bool> visited(size_, false);

		for (size_t i = 0; i < size_; ++i) {
			if (!visited[i] && is_active_node(i)) {
				std::vector<size_t> component;
				std::queue<size_t> queue;

				queue.push(i);
				visited[i] = true;

				while (!queue.empty()) {
					size_t current = queue.front();
					queue.pop();
					component.push_back(current);

					auto neighbors = get_valid_neighbors(current);
					for (const auto& neighbor : neighbors) {
						if (!visited[neighbor.index]) {
							visited[neighbor.index] = true;
							queue.push(neighbor.index);
						}
					}
				}

				components.push_back(std::move(component));
			}
		}

		return components;
	}

	// ======================
	// GEOMETRIC UTILITIES
	// ======================

	/**
	 * Get nodes within a certain radius
	 */
	std::vector<size_t> nodes_within_radius(size_t center, double radius) const
	{
		std::vector<size_t> result;
		auto [center_row, center_col] = to_2d(center);

		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				size_t idx = to_1d(r, c);
				if (is_active_node(idx)) {
					double dist = std::sqrt(
						std::pow(static_cast<double>(r) - static_cast<double>(center_row),
										 2) +
						std::pow(static_cast<double>(c) - static_cast<double>(center_col),
										 2));
					if (dist <= radius) {
						result.push_back(idx);
					}
				}
			}
		}

		return result;
	}

	/**
	 * Get nodes in a rectangular region
	 */
	std::vector<size_t> nodes_in_rectangle(size_t top_left_row,
																				 size_t top_left_col,
																				 size_t bottom_right_row,
																				 size_t bottom_right_col) const
	{
		std::vector<size_t> result;

		for (size_t r = top_left_row; r <= std::min(bottom_right_row, rows_ - 1);
				 ++r) {
			for (size_t c = top_left_col; c <= std::min(bottom_right_col, cols_ - 1);
					 ++c) {
				size_t idx = to_1d(r, c);
				if (is_active_node(idx)) {
					result.push_back(idx);
				}
			}
		}

		return result;
	}

	// ======================
	// FLOW PATTERN ANALYSIS
	// ======================

	/**
	 * Trace flow path from a starting point
	 */
	std::vector<size_t> trace_flow_path(size_t start,
																			const ArrayRef<T>& elevation,
																			size_t max_steps = 1000) const
	{
		std::vector<size_t> path;
		size_t current = start;
		std::unordered_set<size_t> visited;

		path.push_back(current);
		visited.insert(current);

		for (size_t step = 0; step < max_steps; ++step) {
			Direction steepest = get_steepest_descent_direction(current, elevation);

			if (steepest == Direction::INVALID)
				break;

			Neighbor next = get_effective_neighbor(
				to_2d(current).first, to_2d(current).second, steepest);
			if (!next.is_valid || visited.count(next.index))
				break;

			current = next.index;
			path.push_back(current);
			visited.insert(current);

			// Check if we've reached an outlet
			NodeType bc = get_boundary_type(current);
			if (bc == NodeType::HAS_TO_OUT || bc == NodeType::CAN_OUT)
				break;
		}

		return path;
	}
};

// Type aliases for common use cases
using ConnectorF32 = Connector<float>;
using ConnectorF64 = Connector<double>;
using ConnectorI32 = Connector<int32_t>;

} // namespace dagger2
