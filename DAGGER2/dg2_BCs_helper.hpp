#pragma once

#include "dg2_BCs.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * BCBuilder: Comprehensive boundary condition builder with many preset options
 *
 * This class provides a wide range of tools to create and update boundary
 * conditions for landscape modeling, including elevation-based rules, geometric
 * patterns, and complex multi-criteria setups.
 */
template<typename T = double>
class BCBuilder
{
private:
	size_t rows_;
	size_t cols_;
	size_t size_;
	std::shared_ptr<Grid2D<NodeType>> bc_grid_;

public:
	/**
	 * Constructor - creates a new boundary condition grid
	 */
	BCBuilder(size_t rows, size_t cols)
		: rows_(rows)
		, cols_(cols)
		, size_(rows * cols)
	{

		if (rows == 0 || cols == 0) {
			throw std::invalid_argument("Grid dimensions must be positive");
		}

		// Initialize with all NORMAL nodes
		std::vector<NodeType> bc_data(size_, NodeType::NORMAL);
		bc_grid_ = std::make_shared<Grid2D<NodeType>>(bc_data, rows_, cols_);
	}

	/**
	 * Constructor from existing boundary grid (for updates)
	 */
	explicit BCBuilder(std::shared_ptr<Grid2D<NodeType>> bc_grid)
		: bc_grid_(bc_grid)
	{

		if (!bc_grid_) {
			throw std::invalid_argument("Boundary grid cannot be null");
		}

		rows_ = bc_grid_->rows();
		cols_ = bc_grid_->cols();
		size_ = bc_grid_->size();
	}

	// ======================
	// BASIC ACCESSORS
	// ======================

	size_t rows() const { return rows_; }
	size_t cols() const { return cols_; }
	size_t size() const { return size_; }

	std::shared_ptr<Grid2D<NodeType>> get_grid() const { return bc_grid_; }

	NodeType get(size_t row, size_t col) const { return (*bc_grid_)(row, col); }

	void set(size_t row, size_t col, NodeType type)
	{
		(*bc_grid_)(row, col) = type;
	}

	inline size_t to_1d(size_t row, size_t col) const
	{
		return row * cols_ + col;
	}

	inline std::pair<size_t, size_t> to_2d(size_t index) const
	{
		return { index / cols_, index % cols_ };
	}

	// ======================
	// BASIC PATTERNS
	// ======================

	/**
	 * Fill entire grid with specified type
	 */
	BCBuilder& fill(NodeType type)
	{
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			set(r, c, type);
		}
		return *this;
	}

	/**
	 * Convert BC grid to Grid2D<uint8_t> for numpy compatibility
	 */
	std::shared_ptr<Grid2D<uint8_t>> to_numpy_grid() const
	{
		std::vector<uint8_t> data(size_);
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			data[i] = static_cast<uint8_t>(get(r, c));
		}
		return std::make_shared<Grid2D<uint8_t>>(data, rows_, cols_);
	}

	/**
	 * Reset to all NORMAL nodes
	 */
	BCBuilder& reset() { return fill(NodeType::NORMAL); }

	/**
	 * Set all border nodes to specified type
	 */
	BCBuilder& set_borders(NodeType type)
	{
		// Top and bottom rows
		for (size_t c = 0; c < cols_; ++c) {
			set(0, c, type);				 // Top
			set(rows_ - 1, c, type); // Bottom
		}

		// Left and right columns (avoiding corners already set)
		for (size_t r = 1; r < rows_ - 1; ++r) {
			set(r, 0, type);				 // Left
			set(r, cols_ - 1, type); // Right
		}

		return *this;
	}

	/**
	 * Set specific border edges
	 */
	BCBuilder& set_north_border(NodeType type)
	{
		for (size_t c = 0; c < cols_; ++c) {
			set(0, c, type);
		}
		return *this;
	}

	BCBuilder& set_south_border(NodeType type)
	{
		for (size_t c = 0; c < cols_; ++c) {
			set(rows_ - 1, c, type);
		}
		return *this;
	}

	BCBuilder& set_east_border(NodeType type)
	{
		for (size_t r = 0; r < rows_; ++r) {
			set(r, cols_ - 1, type);
		}
		return *this;
	}

	BCBuilder& set_west_border(NodeType type)
	{
		for (size_t r = 0; r < rows_; ++r) {
			set(r, 0, type);
		}
		return *this;
	}

	// ======================
	// PRESET GEOMETRIC LAYOUTS
	// ======================

	/**
	 * All borders open (CAN_OUT)
	 */
	BCBuilder& open_borders() { return set_borders(NodeType::CAN_OUT); }

	/**
	 * All borders closed (REFLECT)
	 */
	BCBuilder& closed_borders() { return set_borders(NodeType::REFLECT); }

	/**
	 * All borders forced outlets
	 */
	BCBuilder& outlet_borders() { return set_borders(NodeType::HAS_TO_OUT); }

	/**
	 * Periodic boundaries on all edges
	 */
	BCBuilder& periodic_borders() { return set_borders(NodeType::PERIODIC); }

	/**
	 * North-South periodic, East-West open
	 */
	BCBuilder& ns_periodic_ew_open()
	{
		set_north_border(NodeType::PERIODIC);
		set_south_border(NodeType::PERIODIC);
		set_east_border(NodeType::CAN_OUT);
		set_west_border(NodeType::CAN_OUT);
		return *this;
	}

	/**
	 * East-West periodic, North-South open
	 */
	BCBuilder& ew_periodic_ns_open()
	{
		set_east_border(NodeType::PERIODIC);
		set_west_border(NodeType::PERIODIC);
		set_north_border(NodeType::CAN_OUT);
		set_south_border(NodeType::CAN_OUT);
		return *this;
	}

	/**
	 * North-South periodic, East-West closed
	 */
	BCBuilder& ns_periodic_ew_closed()
	{
		set_north_border(NodeType::PERIODIC);
		set_south_border(NodeType::PERIODIC);
		set_east_border(NodeType::REFLECT);
		set_west_border(NodeType::REFLECT);
		return *this;
	}

	/**
	 * East-West periodic, North-South closed
	 */
	BCBuilder& ew_periodic_ns_closed()
	{
		set_east_border(NodeType::PERIODIC);
		set_west_border(NodeType::PERIODIC);
		set_north_border(NodeType::REFLECT);
		set_south_border(NodeType::REFLECT);
		return *this;
	}

	/**
	 * Mixed boundary: North inlet, South outlet, East-West closed
	 */
	BCBuilder& inlet_north_outlet_south()
	{
		set_north_border(NodeType::IN);
		set_south_border(NodeType::HAS_TO_OUT);
		set_east_border(NodeType::REFLECT);
		set_west_border(NodeType::REFLECT);
		return *this;
	}

	/**
	 * Mixed boundary: West inlet, East outlet, North-South closed
	 */
	BCBuilder& inlet_west_outlet_east()
	{
		set_west_border(NodeType::IN);
		set_east_border(NodeType::HAS_TO_OUT);
		set_north_border(NodeType::REFLECT);
		set_south_border(NodeType::REFLECT);
		return *this;
	}

	/**
	 * Corner outlets - only corners are outlets, rest closed
	 */
	BCBuilder& corner_outlets()
	{
		set_borders(NodeType::REFLECT);

		// Set corners as outlets
		set(0, 0, NodeType::HAS_TO_OUT);								 // NW
		set(0, cols_ - 1, NodeType::HAS_TO_OUT);				 // NE
		set(rows_ - 1, 0, NodeType::HAS_TO_OUT);				 // SW
		set(rows_ - 1, cols_ - 1, NodeType::HAS_TO_OUT); // SE

		return *this;
	}

	/**
	 * Single center outlet (for circular drainage)
	 */
	BCBuilder& center_outlet()
	{
		set_borders(NodeType::REFLECT);
		set(rows_ / 2, cols_ / 2, NodeType::HAS_TO_OUT);
		return *this;
	}

	/**
	 * Multiple random outlets
	 */
	BCBuilder& random_outlets(size_t num_outlets, uint32_t seed = 0)
	{
		std::mt19937 rng(seed);
		set_borders(NodeType::REFLECT);

		std::uniform_int_distribution<size_t> row_dist(0, rows_ - 1);
		std::uniform_int_distribution<size_t> col_dist(0, cols_ - 1);

		std::unordered_set<size_t> outlet_positions;

		while (outlet_positions.size() < num_outlets) {
			size_t r = row_dist(rng);
			size_t c = col_dist(rng);
			size_t idx = to_1d(r, c);

			if (outlet_positions.find(idx) == outlet_positions.end()) {
				outlet_positions.insert(idx);
				set(r, c, NodeType::HAS_TO_OUT);
			}
		}

		return *this;
	}

	// ======================
	// ELEVATION-BASED PATTERNS
	// ======================

	/**
	 * Set sea level - nodes below become NO_DATA, neighbors become CAN_OUT
	 */
	BCBuilder& set_sea_level(const ArrayRef<T>& elevation, T sea_level)
	{
		// First pass: identify sea nodes
		std::vector<bool> is_sea(size_, false);
		for (size_t i = 0; i < size_; ++i) {
			if (elevation[i] <= sea_level) {
				auto [r, c] = to_2d(i);
				set(r, c, NodeType::NO_DATA);
				is_sea[i] = true;
			}
		}

		// Second pass: set coastal nodes
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				size_t idx = to_1d(r, c);
				if (!is_sea[idx] && get(r, c) == NodeType::NORMAL) {
					// Check if this land node neighbors any sea
					bool neighbors_sea = false;

					// Check 8 neighbors
					for (int dr = -1; dr <= 1; ++dr) {
						for (int dc = -1; dc <= 1; ++dc) {
							if (dr == 0 && dc == 0)
								continue;

							int nr = static_cast<int>(r) + dr;
							int nc = static_cast<int>(c) + dc;

							if (nr >= 0 && nr < static_cast<int>(rows_) && nc >= 0 &&
									nc < static_cast<int>(cols_)) {
								size_t nidx = to_1d(nr, nc);
								if (is_sea[nidx]) {
									neighbors_sea = true;
									break;
								}
							}
						}
						if (neighbors_sea)
							break;
					}

					if (neighbors_sea) {
						set(r, c, NodeType::CAN_OUT);
					}
				}
			}
		}

		return *this;
	}

	/**
	 * Set elevation-based outlets - lowest N points become outlets
	 */
	BCBuilder& lowest_elevation_outlets(const ArrayRef<T>& elevation,
																			size_t num_outlets)
	{
		// Create vector of (elevation, index) pairs
		std::vector<std::pair<T, size_t>> elev_idx;
		elev_idx.reserve(size_);

		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			if (get(r, c) != NodeType::NO_DATA) {
				elev_idx.emplace_back(elevation[i], i);
			}
		}

		// Sort by elevation
		std::sort(elev_idx.begin(), elev_idx.end());

		// Set lowest points as outlets
		size_t count = std::min(num_outlets, elev_idx.size());
		for (size_t i = 0; i < count; ++i) {
			auto [r, c] = to_2d(elev_idx[i].second);
			set(r, c, NodeType::HAS_TO_OUT);
		}

		return *this;
	}

	/**
	 * Set elevation-based inlets - highest N points become inlets
	 */
	BCBuilder& highest_elevation_inlets(const ArrayRef<T>& elevation,
																			size_t num_inlets)
	{
		// Create vector of (elevation, index) pairs
		std::vector<std::pair<T, size_t>> elev_idx;
		elev_idx.reserve(size_);

		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			if (get(r, c) != NodeType::NO_DATA) {
				elev_idx.emplace_back(elevation[i], i);
			}
		}

		// Sort by elevation (descending)
		std::sort(elev_idx.begin(), elev_idx.end(), std::greater<>());

		// Set highest points as inlets
		size_t count = std::min(num_inlets, elev_idx.size());
		for (size_t i = 0; i < count; ++i) {
			auto [r, c] = to_2d(elev_idx[i].second);
			set(r, c, NodeType::IN);
		}

		return *this;
	}

	/**
	 * Set mountain ranges as ridges (REFLECT) based on elevation threshold
	 */
	BCBuilder& mountain_ridges(const ArrayRef<T>& elevation, T ridge_threshold)
	{
		for (size_t i = 0; i < size_; ++i) {
			if (elevation[i] >= ridge_threshold) {
				auto [r, c] = to_2d(i);
				if (get(r, c) != NodeType::NO_DATA) {
					set(r, c, NodeType::REFLECT);
				}
			}
		}
		return *this;
	}

	/**
	 * Valley bottoms as outlets (local elevation minima)
	 */
	BCBuilder& valley_outlets(const ArrayRef<T>& elevation, T tolerance = 0.01)
	{
		for (size_t r = 1; r < rows_ - 1; ++r) {
			for (size_t c = 1; c < cols_ - 1; ++c) {
				size_t idx = to_1d(r, c);
				if (get(r, c) == NodeType::NO_DATA)
					continue;

				T center_elev = elevation[idx];
				bool is_minimum = true;

				// Check all 8 neighbors
				for (int dr = -1; dr <= 1; ++dr) {
					for (int dc = -1; dc <= 1; ++dc) {
						if (dr == 0 && dc == 0)
							continue;

						size_t nr = r + dr;
						size_t nc = c + dc;
						size_t nidx = to_1d(nr, nc);

						if (get(nr, nc) != NodeType::NO_DATA &&
								elevation[nidx] < center_elev + tolerance) {
							is_minimum = false;
							break;
						}
					}
					if (!is_minimum)
						break;
				}

				if (is_minimum) {
					set(r, c, NodeType::HAS_TO_OUT);
				}
			}
		}
		return *this;
	}

	/**
	 * Elevation bands - different BCs for different elevation ranges
	 */
	BCBuilder& elevation_bands(const ArrayRef<T>& elevation,
														 const std::vector<T>& thresholds,
														 const std::vector<NodeType>& types)
	{
		if (thresholds.size() != types.size()) {
			throw std::invalid_argument("Thresholds and types must have same size");
		}

		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			if (get(r, c) == NodeType::NO_DATA)
				continue;

			T elev = elevation[i];
			for (size_t j = 0; j < thresholds.size(); ++j) {
				if (elev <= thresholds[j]) {
					set(r, c, types[j]);
					break;
				}
			}
		}

		return *this;
	}

	// ======================
	// MASK-BASED PATTERNS
	// ======================

	/**
	 * Apply mask - set nodes to NO_DATA where mask is false
	 */
	BCBuilder& apply_mask(const ArrayRef<bool>& mask)
	{
		if (mask.size() != size_) {
			throw std::invalid_argument("Mask size must match grid size");
		}

		for (size_t i = 0; i < size_; ++i) {
			if (!mask[i]) {
				auto [r, c] = to_2d(i);
				set(r, c, NodeType::NO_DATA);
			}
		}

		return *this;
	}

	/**
	 * Apply inverted mask - set nodes to NO_DATA where mask is true
	 */
	BCBuilder& apply_inverted_mask(const ArrayRef<bool>& mask)
	{
		if (mask.size() != size_) {
			throw std::invalid_argument("Mask size must match grid size");
		}

		for (size_t i = 0; i < size_; ++i) {
			if (mask[i]) {
				auto [r, c] = to_2d(i);
				set(r, c, NodeType::NO_DATA);
			}
		}

		return *this;
	}

	/**
	 * Circular mask - nodes outside radius become NO_DATA
	 */
	BCBuilder& circular_mask(size_t center_row, size_t center_col, double radius)
	{
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				double dr = static_cast<double>(r) - static_cast<double>(center_row);
				double dc = static_cast<double>(c) - static_cast<double>(center_col);
				double dist = std::sqrt(dr * dr + dc * dc);

				if (dist > radius) {
					set(r, c, NodeType::NO_DATA);
				}
			}
		}

		return *this;
	}

	/**
	 * Circular mask centered on grid
	 */
	BCBuilder& circular_mask(double radius)
	{
		return circular_mask(rows_ / 2, cols_ / 2, radius);
	}

	/**
	 * Rectangular mask
	 */
	BCBuilder& rectangular_mask(size_t top_row,
															size_t left_col,
															size_t bottom_row,
															size_t right_col)
	{
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				if (r < top_row || r > bottom_row || c < left_col || c > right_col) {
					set(r, c, NodeType::NO_DATA);
				}
			}
		}

		return *this;
	}

	/**
	 * Diamond/rhombus mask
	 */
	BCBuilder& diamond_mask(size_t center_row, size_t center_col, size_t radius)
	{
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				size_t manhattan_dist =
					std::abs(static_cast<int>(r) - static_cast<int>(center_row)) +
					std::abs(static_cast<int>(c) - static_cast<int>(center_col));

				if (manhattan_dist > radius) {
					set(r, c, NodeType::NO_DATA);
				}
			}
		}

		return *this;
	}

	// ======================
	// NEIGHBOR-BASED UPDATES
	// ======================

	/**
	 * Fix borders around NO_DATA - neighboring nodes become CAN_OUT
	 */
	BCBuilder& fix_nodata_borders()
	{
		// First, identify all NO_DATA nodes
		std::vector<bool> is_nodata(size_, false);
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			if (get(r, c) == NodeType::NO_DATA) {
				is_nodata[i] = true;
			}
		}

		// Second pass: update neighbors of NO_DATA
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				size_t idx = to_1d(r, c);
				NodeType current_type = get(r, c);

				// Only update NORMAL nodes
				if (current_type == NodeType::NORMAL) {
					bool neighbors_nodata = false;

					// Check 8 neighbors
					for (int dr = -1; dr <= 1; ++dr) {
						for (int dc = -1; dc <= 1; ++dc) {
							if (dr == 0 && dc == 0)
								continue;

							int nr = static_cast<int>(r) + dr;
							int nc = static_cast<int>(c) + dc;

							if (nr >= 0 && nr < static_cast<int>(rows_) && nc >= 0 &&
									nc < static_cast<int>(cols_)) {
								size_t nidx = to_1d(nr, nc);
								if (is_nodata[nidx]) {
									neighbors_nodata = true;
									break;
								}
							} else {
								// Also consider grid boundary as NO_DATA
								neighbors_nodata = true;
								break;
							}
						}
						if (neighbors_nodata)
							break;
					}

					if (neighbors_nodata) {
						set(r, c, NodeType::CAN_OUT);
					}
				}
			}
		}

		return *this;
	}

	/**
	 * Expand boundary type by N cells
	 */
	BCBuilder& expand_boundary_type(NodeType source_type,
																	NodeType target_type,
																	size_t expansion_cells)
	{
		for (size_t iteration = 0; iteration < expansion_cells; ++iteration) {
			std::vector<std::pair<size_t, size_t>> to_update;

			for (size_t r = 0; r < rows_; ++r) {
				for (size_t c = 0; c < cols_; ++c) {
					if (get(r, c) == NodeType::NORMAL) {
						// Check if neighbors source_type
						bool neighbors_source = false;

						for (int dr = -1; dr <= 1; ++dr) {
							for (int dc = -1; dc <= 1; ++dc) {
								if (dr == 0 && dc == 0)
									continue;

								int nr = static_cast<int>(r) + dr;
								int nc = static_cast<int>(c) + dc;

								if (nr >= 0 && nr < static_cast<int>(rows_) && nc >= 0 &&
										nc < static_cast<int>(cols_)) {
									if (get(nr, nc) == source_type) {
										neighbors_source = true;
										break;
									}
								}
							}
							if (neighbors_source)
								break;
						}

						if (neighbors_source) {
							to_update.emplace_back(r, c);
						}
					}
				}
			}

			// Apply updates
			for (const auto& [r, c] : to_update) {
				set(r, c, target_type);
			}
		}

		return *this;
	}

	/**
	 * Create buffer zones around specific boundary types
	 */
	BCBuilder& create_buffer_zone(NodeType center_type,
																NodeType buffer_type,
																size_t buffer_size)
	{
		return expand_boundary_type(center_type, buffer_type, buffer_size);
	}

	// ======================
	// COMPLEX PATTERNS
	// ======================

	/**
	 * River system - create dendritic drainage from outlets
	 */
	BCBuilder& river_system(const ArrayRef<T>& elevation,
													size_t num_outlets,
													double width_factor = 0.02)
	{
		// First set some random outlets at low elevation
		lowest_elevation_outlets(elevation, num_outlets);

		// Find all outlets
		std::vector<size_t> outlets;
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			if (get(r, c) == NodeType::HAS_TO_OUT) {
				outlets.push_back(i);
			}
		}

		// For each outlet, trace upstream and create rivers
		double river_width = std::max(1.0, width_factor * std::min(rows_, cols_));

		for (size_t outlet_idx : outlets) {
			// Simple upstream propagation based on elevation
			std::queue<size_t> river_queue;
			std::unordered_set<size_t> visited;

			river_queue.push(outlet_idx);
			visited.insert(outlet_idx);

			while (!river_queue.empty()) {
				size_t current = river_queue.front();
				river_queue.pop();

				auto [r, c] = to_2d(current);
				T current_elev = elevation[current];

				// Check neighbors for upstream flow
				for (int dr = -1; dr <= 1; ++dr) {
					for (int dc = -1; dc <= 1; ++dc) {
						if (dr == 0 && dc == 0)
							continue;

						int nr = static_cast<int>(r) + dr;
						int nc = static_cast<int>(c) + dc;

						if (nr >= 0 && nr < static_cast<int>(rows_) && nc >= 0 &&
								nc < static_cast<int>(cols_)) {
							size_t nidx = to_1d(nr, nc);

							if (visited.find(nidx) == visited.end() &&
									elevation[nidx] > current_elev &&
									get(nr, nc) == NodeType::NORMAL) {

								// Make this part of river system
								set(nr, nc, NodeType::CAN_OUT);
								visited.insert(nidx);
								river_queue.push(nidx);
							}
						}
					}
				}
			}
		}

		return *this;
	}

	/**
	 * Watershed divides based on elevation
	 */
	BCBuilder& watershed_divides(const ArrayRef<T>& elevation,
															 T threshold_gradient = 0.1)
	{
		for (size_t r = 1; r < rows_ - 1; ++r) {
			for (size_t c = 1; c < cols_ - 1; ++c) {
				size_t idx = to_1d(r, c);
				if (get(r, c) != NodeType::NORMAL)
					continue;

				T center_elev = elevation[idx];
				T max_gradient = 0.0;

				// Calculate maximum gradient to neighbors
				for (int dr = -1; dr <= 1; ++dr) {
					for (int dc = -1; dc <= 1; ++dc) {
						if (dr == 0 && dc == 0)
							continue;

						size_t nr = r + dr;
						size_t nc = c + dc;
						size_t nidx = to_1d(nr, nc);

						double distance = (dr == 0 || dc == 0) ? 1.0 : 1.4142135623730951;
						T gradient = std::abs(elevation[nidx] - center_elev) / distance;
						max_gradient = std::max(max_gradient, gradient);
					}
				}

				if (max_gradient > threshold_gradient) {
					set(r, c, NodeType::REFLECT);
				}
			}
		}

		return *this;
	}

	/**
	 * Fault lines - linear features as barriers
	 */
	BCBuilder& add_fault_line(size_t start_row,
														size_t start_col,
														size_t end_row,
														size_t end_col,
														size_t width = 1,
														NodeType fault_type = NodeType::REFLECT)
	{
		// Bresenham's line algorithm with width
		int x0 = static_cast<int>(start_col);
		int y0 = static_cast<int>(start_row);
		int x1 = static_cast<int>(end_col);
		int y1 = static_cast<int>(end_row);

		int dx = std::abs(x1 - x0);
		int dy = std::abs(y1 - y0);
		int sx = (x0 < x1) ? 1 : -1;
		int sy = (y0 < y1) ? 1 : -1;
		int err = dx - dy;

		int x = x0, y = y0;

		while (true) {
			// Set fault with width
			for (int wy = -static_cast<int>(width / 2);
					 wy <= static_cast<int>(width / 2);
					 ++wy) {
				for (int wx = -static_cast<int>(width / 2);
						 wx <= static_cast<int>(width / 2);
						 ++wx) {
					int fx = x + wx;
					int fy = y + wy;

					if (fx >= 0 && fx < static_cast<int>(cols_) && fy >= 0 &&
							fy < static_cast<int>(rows_)) {
						set(fy, fx, fault_type);
					}
				}
			}

			if (x == x1 && y == y1)
				break;

			int e2 = 2 * err;
			if (e2 > -dy) {
				err -= dy;
				x += sx;
			}
			if (e2 < dx) {
				err += dx;
				y += sy;
			}
		}

		return *this;
	}

	/**
	 * Island chains - multiple circular features
	 */
	BCBuilder& island_chain(const std::vector<std::pair<size_t, size_t>>& centers,
													const std::vector<double>& radii,
													NodeType island_type = NodeType::REFLECT)
	{
		if (centers.size() != radii.size()) {
			throw std::invalid_argument("Centers and radii must have same size");
		}

		for (size_t i = 0; i < centers.size(); ++i) {
			auto [center_r, center_c] = centers[i];
			double radius = radii[i];

			for (size_t r = 0; r < rows_; ++r) {
				for (size_t c = 0; c < cols_; ++c) {
					double dr = static_cast<double>(r) - static_cast<double>(center_r);
					double dc = static_cast<double>(c) - static_cast<double>(center_c);
					double dist = std::sqrt(dr * dr + dc * dc);

					if (dist <= radius) {
						set(r, c, island_type);
					}
				}
			}
		}

		return *this;
	}

	// ======================
	// UPDATE FUNCTIONS
	// ======================

	/**
	 * Update existing BC based on elevation and sea level
	 */
	BCBuilder& update_sea_level(const ArrayRef<T>& elevation, T sea_level)
	{
		return set_sea_level(elevation, sea_level);
	}

	/**
	 * Update BC to add more outlets at low elevations
	 */
	BCBuilder& update_add_outlets(const ArrayRef<T>& elevation,
																size_t additional_outlets)
	{
		return lowest_elevation_outlets(elevation, additional_outlets);
	}

	/**
	 * Update BC to add more inlets at high elevations
	 */
	BCBuilder& update_add_inlets(const ArrayRef<T>& elevation,
															 size_t additional_inlets)
	{
		return highest_elevation_inlets(elevation, additional_inlets);
	}

	/**
	 * Selective update - only modify nodes of specific type
	 */
	BCBuilder& update_selective(NodeType target_type, NodeType new_type)
	{
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			if (get(r, c) == target_type) {
				set(r, c, new_type);
			}
		}
		return *this;
	}

	/**
	 * Conditional update based on elevation
	 */
	BCBuilder& update_by_elevation(const ArrayRef<T>& elevation,
																 T min_elev,
																 T max_elev,
																 NodeType new_type)
	{
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			T elev = elevation[i];
			if (elev >= min_elev && elev <= max_elev &&
					get(r, c) != NodeType::NO_DATA) {
				set(r, c, new_type);
			}
		}
		return *this;
	}

	/**
	 * Update based on custom predicate function
	 */
	BCBuilder& update_conditional(
		std::function<bool(size_t, size_t, NodeType)> predicate,
		NodeType new_type)
	{
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				NodeType current = get(r, c);
				if (predicate(r, c, current)) {
					set(r, c, new_type);
				}
			}
		}
		return *this;
	}

	/**
	 * Update based on distance from specific points
	 */
	BCBuilder& update_by_distance(
		const std::vector<std::pair<size_t, size_t>>& reference_points,
		double max_distance,
		NodeType new_type)
	{
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				if (get(r, c) == NodeType::NO_DATA)
					continue;

				bool within_distance = false;
				for (const auto& [ref_r, ref_c] : reference_points) {
					double dr = static_cast<double>(r) - static_cast<double>(ref_r);
					double dc = static_cast<double>(c) - static_cast<double>(ref_c);
					double dist = std::sqrt(dr * dr + dc * dc);

					if (dist <= max_distance) {
						within_distance = true;
						break;
					}
				}

				if (within_distance) {
					set(r, c, new_type);
				}
			}
		}
		return *this;
	}

	/**
	 * Smooth boundaries - replace isolated boundary types
	 */
	BCBuilder& smooth_boundaries(size_t iterations = 1)
	{
		for (size_t iter = 0; iter < iterations; ++iter) {
			auto old_grid = *bc_grid_; // Copy current state

			for (size_t r = 1; r < rows_ - 1; ++r) {
				for (size_t c = 1; c < cols_ - 1; ++c) {
					NodeType current = old_grid(r, c);
					if (current == NodeType::NO_DATA)
						continue;

					// Count neighbor types
					std::unordered_map<NodeType, size_t> neighbor_counts;

					for (int dr = -1; dr <= 1; ++dr) {
						for (int dc = -1; dc <= 1; ++dc) {
							if (dr == 0 && dc == 0)
								continue;

							NodeType neighbor_type = old_grid(r + dr, c + dc);
							if (neighbor_type != NodeType::NO_DATA) {
								neighbor_counts[neighbor_type]++;
							}
						}
					}

					// Find most common neighbor type
					NodeType most_common = current;
					size_t max_count = 0;

					for (const auto& [type, count] : neighbor_counts) {
						if (count > max_count) {
							max_count = count;
							most_common = type;
						}
					}

					// Only change if significantly outnumbered
					size_t current_type_count = neighbor_counts[current];
					if (max_count > current_type_count + 2) {
						set(r, c, most_common);
					}
				}
			}
		}

		return *this;
	}

	// ======================
	// ANALYSIS AND VALIDATION
	// ======================

	/**
	 * Count nodes of each boundary type
	 */
	std::unordered_map<NodeType, size_t> count_boundary_types() const
	{
		std::unordered_map<NodeType, size_t> counts;

		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			NodeType type = get(r, c);
			counts[type]++;
		}

		return counts;
	}

	/**
	 * Validate boundary consistency
	 */
	std::vector<std::string> validate() const
	{
		std::vector<std::string> issues;

		auto counts = count_boundary_types();

		// Check for common issues
		if (counts[NodeType::NO_DATA] == size_) {
			issues.push_back("All nodes are NO_DATA");
		}

		if (counts[NodeType::HAS_TO_OUT] == 0 && counts[NodeType::CAN_OUT] == 0) {
			issues.push_back("No outlet nodes found - flow cannot exit domain");
		}

		// Check for orphaned inlet nodes
		if (counts[NodeType::IN] > 0) {
			bool has_flow_path = false;
			for (size_t i = 0; i < size_; ++i) {
				auto [r, c] = to_2d(i);
				if (get(r, c) == NodeType::IN) {
					// Simple check: inlet should have at least one normal neighbor
					for (int dr = -1; dr <= 1; ++dr) {
						for (int dc = -1; dc <= 1; ++dc) {
							if (dr == 0 && dc == 0)
								continue;

							int nr = static_cast<int>(r) + dr;
							int nc = static_cast<int>(c) + dc;

							if (nr >= 0 && nr < static_cast<int>(rows_) && nc >= 0 &&
									nc < static_cast<int>(cols_)) {
								NodeType neighbor_type = get(nr, nc);
								if (neighbor_type == NodeType::NORMAL ||
										neighbor_type == NodeType::CAN_OUT) {
									has_flow_path = true;
									break;
								}
							}
						}
						if (has_flow_path)
							break;
					}
					if (has_flow_path)
						break;
				}
			}

			if (!has_flow_path) {
				issues.push_back("Inlet nodes may be isolated from flow network");
			}
		}

		// Check for periodic boundary consistency
		if (counts[NodeType::PERIODIC] > 0) {
			// Basic check: should have even number on opposite boundaries
			size_t north_periodic = 0, south_periodic = 0, east_periodic = 0,
						 west_periodic = 0;

			for (size_t c = 0; c < cols_; ++c) {
				if (get(0, c) == NodeType::PERIODIC)
					north_periodic++;
				if (get(rows_ - 1, c) == NodeType::PERIODIC)
					south_periodic++;
			}

			for (size_t r = 0; r < rows_; ++r) {
				if (get(r, 0) == NodeType::PERIODIC)
					west_periodic++;
				if (get(r, cols_ - 1) == NodeType::PERIODIC)
					east_periodic++;
			}

			if (north_periodic != south_periodic) {
				issues.push_back("North-South periodic boundaries are unbalanced");
			}

			if (east_periodic != west_periodic) {
				issues.push_back("East-West periodic boundaries are unbalanced");
			}
		}

		return issues;
	}

	/**
	 * Get statistics summary
	 */
	std::string get_summary() const
	{
		auto counts = count_boundary_types();
		std::ostringstream ss;

		ss << "Boundary Condition Summary (" << rows_ << "x" << cols_
			 << " grid):\n";

		for (const auto& [type, count] : counts) {
			if (count > 0) {
				double percentage = 100.0 * count / size_;
				ss << "  " << NodeTypeUtils::to_string(type) << ": " << count << " ("
					 << std::fixed << std::setprecision(1) << percentage << "%)\n";
			}
		}

		auto issues = validate();
		if (!issues.empty()) {
			ss << "\nValidation Issues:\n";
			for (const auto& issue : issues) {
				ss << "  - " << issue << "\n";
			}
		} else {
			ss << "\nValidation: PASSED\n";
		}

		return ss.str();
	}

	// ======================
	// CONNECTOR CREATION
	// ======================

	/**
	 * Create a Connector with this boundary condition grid
	 */
	std::shared_ptr<Connector<T>> create_connector(
		ConnectivityType connectivity = ConnectivityType::D8) const
	{
		return std::make_shared<Connector<T>>(rows_, cols_, bc_grid_, connectivity);
	}

	/**
	 * Create multiple connectors with different connectivity
	 */
	std::pair<std::shared_ptr<Connector<T>>, std::shared_ptr<Connector<T>>>
	create_connectors() const
	{
		auto d4_connector = std::make_shared<Connector<T>>(
			rows_, cols_, bc_grid_, ConnectivityType::D4);
		auto d8_connector = std::make_shared<Connector<T>>(
			rows_, cols_, bc_grid_, ConnectivityType::D8);
		return { d4_connector, d8_connector };
	}

	// ======================
	// SERIALIZATION HELPERS
	// ======================

	/**
	 * Export boundary grid as integer array (for saving/loading)
	 */
	std::vector<uint8_t> export_as_uint8() const
	{
		std::vector<uint8_t> result(size_);
		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			result[i] = static_cast<uint8_t>(get(r, c));
		}
		return result;
	}

	/**
	 * Import boundary grid from integer array
	 */
	BCBuilder& import_from_uint8(const std::vector<uint8_t>& data)
	{
		if (data.size() != size_) {
			throw std::invalid_argument("Data size must match grid size");
		}

		for (size_t i = 0; i < size_; ++i) {
			auto [r, c] = to_2d(i);
			set(r, c, static_cast<NodeType>(data[i]));
		}

		return *this;
	}

	/**
	 * Clone the builder (deep copy)
	 */
	BCBuilder clone() const
	{
		BCBuilder new_builder(rows_, cols_);

		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				new_builder.set(r, c, get(r, c));
			}
		}

		return new_builder;
	}
};

// ======================
// CONVENIENCE FUNCTIONS
// ======================

/**
 * Quick builders for common scenarios
 */
template<typename T = double>
class QuickBC
{
public:
	/**
	 * Create simple open domain (all borders CAN_OUT)
	 */
	static std::shared_ptr<Grid2D<NodeType>> open_domain(size_t rows, size_t cols)
	{
		return BCBuilder<T>(rows, cols).open_borders().get_grid();
	}

	/**
	 * Create closed domain (all borders REFLECT)
	 */
	static std::shared_ptr<Grid2D<NodeType>> closed_domain(size_t rows,
																												 size_t cols)
	{
		return BCBuilder<T>(rows, cols).closed_borders().get_grid();
	}

	/**
	 * Create periodic domain
	 */
	static std::shared_ptr<Grid2D<NodeType>> periodic_ew_domain(size_t rows,
																															size_t cols)
	{
		return BCBuilder<T>(rows, cols).ew_periodic_ns_open().get_grid();
	}

	/**
	 * Create periodic domain
	 */
	static std::shared_ptr<Grid2D<NodeType>> periodic_ns_domain(size_t rows,
																															size_t cols)
	{
		return BCBuilder<T>(rows, cols).ns_periodic_ew_open().get_grid();
	}

	/**
	 * Create domain with sea level
	 */
	static std::shared_ptr<Grid2D<NodeType>> sea_level_domain(
		size_t rows,
		size_t cols,
		const ArrayRef<T>& elevation,
		T sea_level)
	{
		return BCBuilder<T>(rows, cols)
			.set_sea_level(elevation, sea_level)
			.get_grid();
	}

	/**
	 * Create circular island domain
	 */
	static std::shared_ptr<Grid2D<NodeType>>
	island_domain(size_t rows, size_t cols, double radius_fraction = 0.4)
	{
		double radius = radius_fraction * std::min(rows, cols) / 2.0;
		return BCBuilder<T>(rows, cols)
			.fill(NodeType::NO_DATA)
			.circular_mask(radius)
			.set(rows / 2, cols / 2, NodeType::NORMAL) // Ensure center is active
			.fix_nodata_borders()
			.get_grid();
	}

	/**
	 * Create simple drainage basin
	 */
	static std::shared_ptr<Grid2D<NodeType>> drainage_basin(
		size_t rows,
		size_t cols,
		const ArrayRef<T>& elevation,
		size_t num_outlets = 1)
	{
		return BCBuilder<T>(rows, cols)
			.closed_borders()
			.lowest_elevation_outlets(elevation, num_outlets)
			.get_grid();
	}

	/**
	 * Create inlet-outlet channel
	 */
	static std::shared_ptr<Grid2D<NodeType>> flow_channel(size_t rows,
																												size_t cols)
	{
		return BCBuilder<T>(rows, cols).inlet_west_outlet_east().get_grid();
	}
};

/**
 * Preset landscape scenarios
 */
template<typename T = double>
class LandscapePresets
{
public:
	/**
	 * Mountain watershed with ridges and valleys
	 */
	static BCBuilder<T> mountain_watershed(size_t rows,
																				 size_t cols,
																				 const ArrayRef<T>& elevation)
	{
		// Get elevation statistics
		T min_elev = *std::min_element(elevation.begin(), elevation.end());
		T max_elev = *std::max_element(elevation.begin(), elevation.end());
		T range = max_elev - min_elev;

		T ridge_threshold = min_elev + 0.8 * range;
		T valley_threshold = min_elev + 0.2 * range;

		return BCBuilder<T>(rows, cols)
			.closed_borders()														 // Start with closed domain
			.mountain_ridges(elevation, ridge_threshold) // High elevations as ridges
			.valley_outlets(elevation)									 // Valley bottoms as outlets
			.fix_nodata_borders();											 // Clean up boundaries
	}

	/**
	 * Coastal landscape with sea level
	 */
	static BCBuilder<T> coastal_landscape(size_t rows,
																				size_t cols,
																				const ArrayRef<T>& elevation,
																				T sea_level)
	{
		return BCBuilder<T>(rows, cols)
			.set_sea_level(elevation, sea_level) // Set sea and coast
			.fix_nodata_borders();							 // Ensure proper drainage
	}

	/**
	 * River delta system
	 */
	static BCBuilder<T> river_delta(size_t rows,
																	size_t cols,
																	const ArrayRef<T>& elevation)
	{
		return BCBuilder<T>(rows, cols)
			.inlet_north_outlet_south()				// Flow from north to south
			.river_system(elevation, 3, 0.05) // Create tributary network
			.fix_nodata_borders();
	}

	/**
	 * Island archipelago
	 */
	static BCBuilder<T> archipelago(size_t rows,
																	size_t cols,
																	size_t num_islands = 5,
																	uint32_t seed = 42)
	{
		std::mt19937 rng(seed);
		std::uniform_int_distribution<size_t> row_dist(rows / 4, 3 * rows / 4);
		std::uniform_int_distribution<size_t> col_dist(cols / 4, 3 * cols / 4);
		std::uniform_real_distribution<double> radius_dist(
			std::min(rows, cols) * 0.05, std::min(rows, cols) * 0.15);

		std::vector<std::pair<size_t, size_t>> centers;
		std::vector<double> radii;

		for (size_t i = 0; i < num_islands; ++i) {
			centers.emplace_back(row_dist(rng), col_dist(rng));
			radii.push_back(radius_dist(rng));
		}

		return BCBuilder<T>(rows, cols)
			.fill(NodeType::NO_DATA)												// Start with all water
			.island_chain(centers, radii, NodeType::NORMAL) // Add islands
			.fix_nodata_borders();													// Set coastal boundaries
	}

	/**
	 * Fault-bounded valley
	 */
	static BCBuilder<T> fault_valley(size_t rows,
																	 size_t cols,
																	 const ArrayRef<T>& elevation)
	{
		size_t mid_row = rows / 2;

		return BCBuilder<T>(rows, cols)
			.closed_borders()
			.add_fault_line(
				0, cols / 4, rows - 1, cols / 4, 2, NodeType::REFLECT) // West fault
			.add_fault_line(0,
											3 * cols / 4,
											rows - 1,
											3 * cols / 4,
											2,
											NodeType::REFLECT)			// East fault
			.lowest_elevation_outlets(elevation, 2) // Valley outlets
			.fix_nodata_borders();
	}
};

} // namespace dagger2
