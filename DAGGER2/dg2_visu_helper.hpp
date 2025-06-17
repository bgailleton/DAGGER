#pragma once

/**
 * dg2_visualization_helpers.hpp
 *
 * Efficient C++ backend functions to generate raw Grid2D data for
 * visualization. These functions use the Connector's D8 connectivity and
 * boundary conditions to create hillshade relief and other creative
 * visualization grids.
 */

#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace dagger2 {

// ======================
// HILLSHADE RELIEF
// ======================

/**
 * Generate D8 hillshade relief using Horn's algorithm with connector awareness
 */
template<typename T>
Grid2D<T>
create_d8_hillshade_relief(
	const Grid2D<T>& elevation,
	const Connector<T>& connector,
	T sun_azimuth = static_cast<T>(315.0),	// degrees (315 = NW)
	T sun_elevation = static_cast<T>(45.0), // degrees
	T z_factor = static_cast<T>(1.0),				// vertical exaggeration
	T cell_size = static_cast<T>(1.0),			// grid cell size
	T no_data_value = static_cast<T>(0.0)		// value for NO_DATA areas
)
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> hillshade_data(rows * cols);

	// Convert angles to radians
	T azimuth_rad = sun_azimuth * static_cast<T>(M_PI) / static_cast<T>(180.0);
	T elevation_rad =
		sun_elevation * static_cast<T>(M_PI) / static_cast<T>(180.0);

	// Calculate sun direction components
	T sun_x = std::sin(azimuth_rad) * std::cos(elevation_rad);
	T sun_y = std::cos(azimuth_rad) * std::cos(elevation_rad);
	T sun_z = std::sin(elevation_rad);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			// Check boundary condition
			NodeType node_type = connector.get_boundary_type(r, c);
			if (node_type == NodeType::NO_DATA) {
				hillshade_data[idx] = no_data_value;
				continue;
			}

			// Get D8 neighbors using connector
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			// Calculate gradients using Horn's algorithm adapted for D8
			T dz_dx = static_cast<T>(0.0);
			T dz_dy = static_cast<T>(0.0);
			T center_elev = elevation(r, c);

			// Use D8 directions to compute gradients
			// East-West gradient (dz/dx)
			T east_sum = static_cast<T>(0.0), west_sum = static_cast<T>(0.0);
			int east_count = 0, west_count = 0;

			// North-South gradient (dz/dy)
			T north_sum = static_cast<T>(0.0), south_sum = static_cast<T>(0.0);
			int north_count = 0, south_count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				Direction dir = neighbor.direction;

				// Accumulate by direction for Horn's algorithm
				switch (dir) {
					case Direction::EAST:
						east_sum += neighbor_elev;
						east_count++;
						break;
					case Direction::WEST:
						west_sum += neighbor_elev;
						west_count++;
						break;
					case Direction::NORTH:
						north_sum += neighbor_elev;
						north_count++;
						break;
					case Direction::SOUTH:
						south_sum += neighbor_elev;
						south_count++;
						break;
					case Direction::NORTHEAST:
						east_sum += neighbor_elev;
						north_sum += neighbor_elev;
						east_count++;
						north_count++;
						break;
					case Direction::NORTHWEST:
						west_sum += neighbor_elev;
						north_sum += neighbor_elev;
						west_count++;
						north_count++;
						break;
					case Direction::SOUTHEAST:
						east_sum += neighbor_elev;
						south_sum += neighbor_elev;
						east_count++;
						south_count++;
						break;
					case Direction::SOUTHWEST:
						west_sum += neighbor_elev;
						south_sum += neighbor_elev;
						west_count++;
						south_count++;
						break;
					default:
						break;
				}
			}

			// Calculate gradients with fallback to center elevation
			if (east_count > 0 && west_count > 0) {
				dz_dx = (east_sum / east_count - west_sum / west_count) /
								(static_cast<T>(2.0) * cell_size);
			} else if (east_count > 0) {
				dz_dx = (east_sum / east_count - center_elev) / cell_size;
			} else if (west_count > 0) {
				dz_dx = (center_elev - west_sum / west_count) / cell_size;
			}

			if (north_count > 0 && south_count > 0) {
				dz_dy = (south_sum / south_count - north_sum / north_count) /
								(static_cast<T>(2.0) * cell_size);
			} else if (north_count > 0) {
				dz_dy = (center_elev - north_sum / north_count) / cell_size;
			} else if (south_count > 0) {
				dz_dy = (south_sum / south_count - center_elev) / cell_size;
			}

			// Apply z-factor
			dz_dx *= z_factor;
			dz_dy *= z_factor;

			// Calculate surface normal
			T slope = std::sqrt(dz_dx * dz_dx + dz_dy * dz_dy);
			T aspect = std::atan2(dz_dy, -dz_dx);

			// Calculate hillshade using dot product
			T hillshade;
			if (slope == static_cast<T>(0.0)) {
				hillshade = std::sin(elevation_rad);
			} else {
				T cos_zenith_slope = std::cos(
					static_cast<T>(M_PI) / static_cast<T>(2.0) - std::atan(slope));
				T sin_zenith_slope = std::sin(
					static_cast<T>(M_PI) / static_cast<T>(2.0) - std::atan(slope));

				hillshade = cos_zenith_slope * std::sin(elevation_rad) +
										sin_zenith_slope * std::cos(elevation_rad) *
											std::cos(azimuth_rad - aspect);
			}

			// Normalize to 1-0 range (inverting and normalsing)
			hillshade_data[idx] = 1 - hillshade;
		}
	}

	return Grid2D<T>(hillshade_data, rows, cols);
}

// ======================
// SLOPE ANALYSIS
// ======================

/**
 * Calculate slope magnitude using D8 connectivity
 */
template<typename T>
Grid2D<T>
create_slope_grid(const Grid2D<T>& elevation,
									const Connector<T>& connector,
									T cell_size = static_cast<T>(1.0),
									bool return_degrees = true,
									T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> slope_data(rows * cols);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				slope_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			T max_slope = static_cast<T>(0.0);

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				T rise = std::abs(neighbor_elev - center_elev);
				T run = neighbor.distance * cell_size;
				T slope = rise / run;

				max_slope = std::max(max_slope, slope);
			}

			if (return_degrees) {
				max_slope =
					std::atan(max_slope) * static_cast<T>(180.0) / static_cast<T>(M_PI);
			}

			slope_data[idx] = max_slope;
		}
	}

	return Grid2D<T>(slope_data, rows, cols);
}

// ======================
// ASPECT ANALYSIS
// ======================

/**
 * Calculate aspect (slope direction) using D8 connectivity
 */
template<typename T>
Grid2D<T>
create_aspect_grid(const Grid2D<T>& elevation,
									 const Connector<T>& connector,
									 T cell_size = static_cast<T>(1.0),
									 T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> aspect_data(rows * cols);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				aspect_data[idx] = no_data_value;
				continue;
			}

			auto neighbors = connector.get_effective_valid_neighbors(r, c);
			T center_elev = elevation(r, c);

			// Calculate gradients for aspect
			T dz_dx = static_cast<T>(0.0);
			T dz_dy = static_cast<T>(0.0);
			int count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				T dz = neighbor_elev - center_elev;

				// Project onto cardinal directions
				switch (neighbor.direction) {
					case Direction::EAST:
						dz_dx += dz;
						count++;
						break;
					case Direction::WEST:
						dz_dx -= dz;
						count++;
						break;
					case Direction::NORTH:
						dz_dy -= dz;
						count++;
						break;
					case Direction::SOUTH:
						dz_dy += dz;
						count++;
						break;
					case Direction::NORTHEAST:
						dz_dx += dz * static_cast<T>(0.707);
						dz_dy -= dz * static_cast<T>(0.707);
						count++;
						break;
					case Direction::NORTHWEST:
						dz_dx -= dz * static_cast<T>(0.707);
						dz_dy -= dz * static_cast<T>(0.707);
						count++;
						break;
					case Direction::SOUTHEAST:
						dz_dx += dz * static_cast<T>(0.707);
						dz_dy += dz * static_cast<T>(0.707);
						count++;
						break;
					case Direction::SOUTHWEST:
						dz_dx -= dz * static_cast<T>(0.707);
						dz_dy += dz * static_cast<T>(0.707);
						count++;
						break;
					default:
						break;
				}
			}

			T aspect;
			if (count == 0 ||
					(dz_dx == static_cast<T>(0.0) && dz_dy == static_cast<T>(0.0))) {
				aspect = static_cast<T>(-1.0); // Flat area
			} else {
				aspect = std::atan2(dz_dy, -dz_dx) * static_cast<T>(180.0) /
								 static_cast<T>(M_PI);
				if (aspect < static_cast<T>(0.0)) {
					aspect += static_cast<T>(360.0);
				}
			}

			aspect_data[idx] = aspect;
		}
	}

	return Grid2D<T>(aspect_data, rows, cols);
}

// ======================
// FLOW ACCUMULATION PROXY
// ======================

/**
 * Create a simple flow accumulation proxy based on local drainage density
 */
template<typename T>
Grid2D<T>
create_drainage_density_grid(const Grid2D<T>& elevation,
														 const Connector<T>& connector,
														 T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> drainage_data(rows * cols, static_cast<T>(0.0));

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				drainage_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			T drainage_score = static_cast<T>(0.0);
			int lower_neighbors = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				if (neighbor_elev < center_elev) {
					lower_neighbors++;
					drainage_score += (center_elev - neighbor_elev) / neighbor.distance;
				}
			}

			// Normalize by number of neighbors for local drainage intensity
			if (lower_neighbors > 0) {
				drainage_score = drainage_score * lower_neighbors *
												 lower_neighbors; // Emphasize convergent areas
			}

			drainage_data[idx] = drainage_score;
		}
	}

	return Grid2D<T>(drainage_data, rows, cols);
}

// ======================
// TERRAIN ROUGHNESS
// ======================

/**
 * Calculate terrain roughness index using D8 neighborhood
 */
template<typename T>
Grid2D<T>
create_roughness_grid(const Grid2D<T>& elevation,
											const Connector<T>& connector,
											T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> roughness_data(rows * cols);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				roughness_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			T sum_squared_diff = static_cast<T>(0.0);
			int count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				T diff = neighbor_elev - center_elev;
				sum_squared_diff += diff * diff;
				count++;
			}

			T roughness =
				(count > 0) ? std::sqrt(sum_squared_diff / count) : static_cast<T>(0.0);
			roughness_data[idx] = roughness;
		}
	}

	return Grid2D<T>(roughness_data, rows, cols);
}

// ======================
// CURVATURE ANALYSIS
// ======================

/**
 * Calculate plan curvature (contour curvature) using D8 connectivity
 */
template<typename T>
Grid2D<T>
create_plan_curvature_grid(const Grid2D<T>& elevation,
													 const Connector<T>& connector,
													 T cell_size = static_cast<T>(1.0),
													 T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> curvature_data(rows * cols);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				curvature_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			// Simple curvature estimation using second derivatives
			T curvature = static_cast<T>(0.0);
			int count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				T second_deriv =
					(neighbor_elev - center_elev) /
					(neighbor.distance * neighbor.distance * cell_size * cell_size);
				curvature += second_deriv;
				count++;
			}

			if (count > 0) {
				curvature /= count;
			}

			curvature_data[idx] = curvature;
		}
	}

	return Grid2D<T>(curvature_data, rows, cols);
}

// ======================
// CREATIVE VISUALIZATION GRIDS
// ======================

/**
 * Generate procedural erosion pattern grid
 */
template<typename T>
Grid2D<T>
create_erosion_pattern_grid(const Grid2D<T>& elevation,
														const Connector<T>& connector,
														T intensity = static_cast<T>(1.0),
														uint32_t seed = 12345,
														T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> erosion_data(rows * cols);

	std::mt19937 gen(seed);
	std::uniform_real_distribution<T> noise_dist(static_cast<T>(0.0),
																							 static_cast<T>(1.0));

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				erosion_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			// Calculate erosion susceptibility based on local geometry
			T slope_factor = static_cast<T>(0.0);
			T convergence_factor = static_cast<T>(0.0);
			int count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				T local_slope =
					std::abs(neighbor_elev - center_elev) / neighbor.distance;
				slope_factor += local_slope;

				if (neighbor_elev > center_elev) {
					convergence_factor += static_cast<T>(1.0);
				}
				count++;
			}

			if (count > 0) {
				slope_factor /= count;
				convergence_factor /= count;
			}

			// Combine with noise for realistic patterns
			T noise = noise_dist(gen);
			T erosion = (slope_factor + convergence_factor) * intensity *
									(static_cast<T>(0.5) + static_cast<T>(0.5) * noise);

			erosion_data[idx] = erosion;
		}
	}

	return Grid2D<T>(erosion_data, rows, cols);
}

/**
 * Generate a visibility/exposure grid (how "exposed" each cell is)
 */
template<typename T>
Grid2D<T>
create_exposure_grid(const Grid2D<T>& elevation,
										 const Connector<T>& connector,
										 T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> exposure_data(rows * cols);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				exposure_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			T relative_height = static_cast<T>(0.0);
			int count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				relative_height += (center_elev - neighbor_elev);
				count++;
			}

			T exposure = (count > 0) ? relative_height / count : static_cast<T>(0.0);
			exposure =
				std::max(static_cast<T>(0.0), exposure); // Only positive exposure

			exposure_data[idx] = exposure;
		}
	}

	return Grid2D<T>(exposure_data, rows, cols);
}

/**
 * Generate a shelter/protection grid (inverse of exposure)
 */
template<typename T>
Grid2D<T>
create_shelter_grid(const Grid2D<T>& elevation,
										const Connector<T>& connector,
										T no_data_value = static_cast<T>(-1.0))
{
	const size_t rows = elevation.rows();
	const size_t cols = elevation.cols();
	std::vector<T> shelter_data(rows * cols);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			size_t idx = connector.to_1d(r, c);

			if (connector.get_boundary_type(r, c) == NodeType::NO_DATA) {
				shelter_data[idx] = no_data_value;
				continue;
			}

			T center_elev = elevation(r, c);
			auto neighbors = connector.get_effective_valid_neighbors(r, c);

			T protection = static_cast<T>(0.0);
			int count = 0;

			for (const auto& neighbor : neighbors) {
				if (!neighbor.is_valid || neighbor.boundary_type == NodeType::NO_DATA) {
					continue;
				}

				T neighbor_elev = elevation[neighbor.index];
				if (neighbor_elev > center_elev) {
					protection += (neighbor_elev - center_elev) / neighbor.distance;
				}
				count++;
			}

			T shelter = (count > 0) ? protection / count : static_cast<T>(0.0);
			shelter_data[idx] = shelter;
		}
	}

	return Grid2D<T>(shelter_data, rows, cols);
}

} // namespace dagger2
