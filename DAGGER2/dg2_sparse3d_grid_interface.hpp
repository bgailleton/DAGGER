#pragma once

#include "dg2_array.hpp"
#include "dg2_sparse3d_data.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace dagger2 {

/**
 * Sparse3DGridInterface: Bridge between sparse 3D data and regular
 * Grid2D/Grid3D
 *
 * This class provides conversion utilities and cross-section extraction
 * capabilities between sparse and dense grid representations.
 */
template<typename T>
class Sparse3DGridInterface
{
private:
	std::shared_ptr<Sparse3DContainer<T>> sparse_data_;

public:
	/**
	 * Constructor
	 */
	explicit Sparse3DGridInterface(
		std::shared_ptr<Sparse3DContainer<T>> sparse_data)
		: sparse_data_(sparse_data)
	{
		if (!sparse_data_) {
			throw std::invalid_argument("Sparse data container cannot be null");
		}
	}

	// ======================
	// GRID2D EXTRACTION
	// ======================

	/**
	 * Extract XY plane at specific Z coordinate
	 */
	std::shared_ptr<Grid2D<T>> extract_xy_plane(int64_t z_coord) const
	{
		auto coords_at_z = sparse_data_->get_coords_at_z(z_coord);
		if (coords_at_z.empty()) {
			// Return empty 1x1 grid with default value
			std::vector<T> data = { sparse_data_->get_default_value() };
			return std::make_shared<Grid2D<T>>(data, 1, 1);
		}

		// Find XY bounds
		int64_t min_x = coords_at_z[0].x, max_x = coords_at_z[0].x;
		int64_t min_y = coords_at_z[0].y, max_y = coords_at_z[0].y;

		for (const auto& coord : coords_at_z) {
			min_x = std::min(min_x, coord.x);
			max_x = std::max(max_x, coord.x);
			min_y = std::min(min_y, coord.y);
			max_y = std::max(max_y, coord.y);
		}

		size_t width = static_cast<size_t>(max_x - min_x + 1);
		size_t height = static_cast<size_t>(max_y - min_y + 1);

		// Fill grid data
		std::vector<T> grid_data(width * height, sparse_data_->get_default_value());

		for (const auto& coord : coords_at_z) {
			size_t grid_x = static_cast<size_t>(coord.x - min_x);
			size_t grid_y = static_cast<size_t>(coord.y - min_y);
			size_t grid_idx = grid_y * width + grid_x;
			grid_data[grid_idx] = sparse_data_->get_value(coord);
		}

		return std::make_shared<Grid2D<T>>(grid_data, height, width);
	}

	/**
	 * Extract XZ plane at specific Y coordinate
	 */
	std::shared_ptr<Grid2D<T>> extract_xz_plane(int64_t y_coord) const
	{
		auto coords_at_y = sparse_data_->get_coords_at_y(y_coord);
		if (coords_at_y.empty()) {
			std::vector<T> data = { sparse_data_->get_default_value() };
			return std::make_shared<Grid2D<T>>(data, 1, 1);
		}

		// Find XZ bounds
		int64_t min_x = coords_at_y[0].x, max_x = coords_at_y[0].x;
		int64_t min_z = coords_at_y[0].z, max_z = coords_at_y[0].z;

		for (const auto& coord : coords_at_y) {
			min_x = std::min(min_x, coord.x);
			max_x = std::max(max_x, coord.x);
			min_z = std::min(min_z, coord.z);
			max_z = std::max(max_z, coord.z);
		}

		size_t width = static_cast<size_t>(max_x - min_x + 1);
		size_t depth = static_cast<size_t>(max_z - min_z + 1);

		// Fill grid data (z as rows, x as columns)
		std::vector<T> grid_data(width * depth, sparse_data_->get_default_value());

		for (const auto& coord : coords_at_y) {
			size_t grid_x = static_cast<size_t>(coord.x - min_x);
			size_t grid_z = static_cast<size_t>(coord.z - min_z);
			size_t grid_idx = grid_z * width + grid_x;
			grid_data[grid_idx] = sparse_data_->get_value(coord);
		}

		return std::make_shared<Grid2D<T>>(grid_data, depth, width);
	}

	/**
	 * Extract YZ plane at specific X coordinate
	 */
	std::shared_ptr<Grid2D<T>> extract_yz_plane(int64_t x_coord) const
	{
		auto coords_at_x = sparse_data_->get_coords_at_x(x_coord);
		if (coords_at_x.empty()) {
			std::vector<T> data = { sparse_data_->get_default_value() };
			return std::make_shared<Grid2D<T>>(data, 1, 1);
		}

		// Find YZ bounds
		int64_t min_y = coords_at_x[0].y, max_y = coords_at_x[0].y;
		int64_t min_z = coords_at_x[0].z, max_z = coords_at_x[0].z;

		for (const auto& coord : coords_at_x) {
			min_y = std::min(min_y, coord.y);
			max_y = std::max(max_y, coord.y);
			min_z = std::min(min_z, coord.z);
			max_z = std::max(max_z, coord.z);
		}

		size_t height = static_cast<size_t>(max_y - min_y + 1);
		size_t depth = static_cast<size_t>(max_z - min_z + 1);

		// Fill grid data (z as rows, y as columns)
		std::vector<T> grid_data(height * depth, sparse_data_->get_default_value());

		for (const auto& coord : coords_at_x) {
			size_t grid_y = static_cast<size_t>(coord.y - min_y);
			size_t grid_z = static_cast<size_t>(coord.z - min_z);
			size_t grid_idx = grid_z * height + grid_y;
			grid_data[grid_idx] = sparse_data_->get_value(coord);
		}

		return std::make_shared<Grid2D<T>>(grid_data, depth, height);
	}

	/**
	 * Extract arbitrary cross-section defined by a plane
	 */
	struct PlaneDefinition
	{
		Coord3D point;	// Point on the plane
		Coord3D normal; // Normal vector (will be normalized)

		PlaneDefinition(const Coord3D& p, const Coord3D& n)
			: point(p)
			, normal(n)
		{
		}
	};

	std::shared_ptr<Grid2D<T>> extract_plane_cross_section(
		const PlaneDefinition& plane,
		size_t resolution_u = 100,
		size_t resolution_v = 100,
		double extent_u = 10.0,
		double extent_v = 10.0) const
	{

		// Create two orthogonal vectors in the plane
		Coord3D u_vec, v_vec;

		// Find a vector not parallel to normal
		if (std::abs(plane.normal.x) < 0.9) {
			u_vec = Coord3D(1, 0, 0);
		} else {
			u_vec = Coord3D(0, 1, 0);
		}

		// Cross product to get first tangent vector
		// u_vec = u_vec - (u_vec · normal) * normal
		double dot_product =
			static_cast<double>(u_vec.x * plane.normal.x + u_vec.y * plane.normal.y +
													u_vec.z * plane.normal.z);
		double normal_len_sq = static_cast<double>(plane.normal.x * plane.normal.x +
																							 plane.normal.y * plane.normal.y +
																							 plane.normal.z * plane.normal.z);

		if (normal_len_sq > 0) {
			double factor = dot_product / normal_len_sq;
			u_vec = u_vec - Coord3D(static_cast<int64_t>(factor * plane.normal.x),
															static_cast<int64_t>(factor * plane.normal.y),
															static_cast<int64_t>(factor * plane.normal.z));
		}

		// Second tangent vector: v_vec = normal × u_vec
		v_vec = Coord3D(plane.normal.y * u_vec.z - plane.normal.z * u_vec.y,
										plane.normal.z * u_vec.x - plane.normal.x * u_vec.z,
										plane.normal.x * u_vec.y - plane.normal.y * u_vec.x);

		// Sample the plane
		std::vector<T> grid_data(resolution_u * resolution_v,
														 sparse_data_->get_default_value());

		for (size_t i = 0; i < resolution_u; ++i) {
			for (size_t j = 0; j < resolution_v; ++j) {
				double u =
					(static_cast<double>(i) / (resolution_u - 1) - 0.5) * extent_u;
				double v =
					(static_cast<double>(j) / (resolution_v - 1) - 0.5) * extent_v;

				Coord3D sample_point =
					plane.point +
					Coord3D(static_cast<int64_t>(u * u_vec.x + v * v_vec.x),
									static_cast<int64_t>(u * u_vec.y + v * v_vec.y),
									static_cast<int64_t>(u * u_vec.z + v * v_vec.z));

				// Find nearest voxel (simple nearest neighbor interpolation)
				T value = sparse_data_->get_value(sample_point);
				grid_data[i * resolution_v + j] = value;
			}
		}

		return std::make_shared<Grid2D<T>>(grid_data, resolution_u, resolution_v);
	}

	// ======================
	// GRID3D CONVERSION
	// ======================

	/**
	 * Convert sparse data to dense Grid3D
	 */
	std::shared_ptr<Grid3D<T>> to_dense_grid3d() const
	{
		if (sparse_data_->empty()) {
			std::vector<T> data = { sparse_data_->get_default_value() };
			return std::make_shared<Grid3D<T>>(data, 1, 1, 1);
		}

		auto bbox = sparse_data_->get_bounding_box();
		if (!bbox.is_valid) {
			std::vector<T> data = { sparse_data_->get_default_value() };
			return std::make_shared<Grid3D<T>>(data, 1, 1, 1);
		}

		size_t depth = static_cast<size_t>(bbox.max_coord.z - bbox.min_coord.z + 1);
		size_t height =
			static_cast<size_t>(bbox.max_coord.y - bbox.min_coord.y + 1);
		size_t width = static_cast<size_t>(bbox.max_coord.x - bbox.min_coord.x + 1);

		std::vector<T> grid_data(depth * height * width,
														 sparse_data_->get_default_value());

		for (const auto& [coord, voxel] : *sparse_data_) {
			size_t grid_x = static_cast<size_t>(coord.x - bbox.min_coord.x);
			size_t grid_y = static_cast<size_t>(coord.y - bbox.min_coord.y);
			size_t grid_z = static_cast<size_t>(coord.z - bbox.min_coord.z);
			size_t grid_idx = grid_z * height * width + grid_y * width + grid_x;
			grid_data[grid_idx] = voxel.value;
		}

		return std::make_shared<Grid3D<T>>(grid_data, depth, height, width);
	}

	/**
	 * Convert Grid3D to sparse representation
	 */
	static std::shared_ptr<Sparse3DContainer<T>> from_dense_grid3d(
		const Grid3D<T>& grid,
		const T& skip_value = T{},
		const Coord3D& offset = Coord3D(0, 0, 0))
	{

		auto sparse_data = std::make_shared<Sparse3DContainer<T>>(skip_value);

		for (size_t d = 0; d < grid.depth(); ++d) {
			for (size_t r = 0; r < grid.rows(); ++r) {
				/**
				 * Convert Grid3D to sparse representation
				 */
				static std::shared_ptr<Sparse3DContainer<T>> from_dense_grid3d(
					const Grid3D<T>& grid,
					const T& skip_value = T{},
					const Coord3D& offset = Coord3D(0, 0, 0))
				{

					auto sparse_data = std::make_shared<Sparse3DContainer<T>>(skip_value);

					for (size_t d = 0; d < grid.depth(); ++d) {
						for (size_t r = 0; r < grid.rows(); ++r) {
							for (size_t c = 0; c < grid.cols(); ++c) {
								T value = grid(d, r, c);
								if (value != skip_value) {
									Coord3D coord(static_cast<int64_t>(c) + offset.x,
																static_cast<int64_t>(r) + offset.y,
																static_cast<int64_t>(d) + offset.z);
									sparse_data->set_value(coord, value);
								}
							}
						}
					}

					return sparse_data;
				}

				/**
				 * Convert Grid2D to sparse representation (as XY plane at Z=0)
				 */
				static std::shared_ptr<Sparse3DContainer<T>> from_grid2d(
					const Grid2D<T>& grid,
					int64_t z_coord = 0,
					const T& skip_value = T{},
					const Coord3D& offset = Coord3D(0, 0, 0))
				{

					auto sparse_data = std::make_shared<Sparse3DContainer<T>>(skip_value);

					for (size_t r = 0; r < grid.rows(); ++r) {
						for (size_t c = 0; c < grid.cols(); ++c) {
							T value = grid(r, c);
							if (value != skip_value) {
								Coord3D coord(static_cast<int64_t>(c) + offset.x,
															static_cast<int64_t>(r) + offset.y,
															z_coord + offset.z);
								sparse_data->set_value(coord, value);
							}
						}
					}

					return sparse_data;
				}

				// ======================
				// PROJECTION OPERATIONS
				// ======================

				/**
				 * Project along Z axis (top-down view)
				 */
				enum class ProjectionMode
				{
					MAX,	 // Maximum value along projection axis
					MIN,	 // Minimum value along projection axis
					MEAN,	 // Mean value along projection axis
					SUM,	 // Sum of values along projection axis
					COUNT, // Count of non-default values
					FIRST, // First non-default value encountered
					LAST	 // Last non-default value encountered
				};

				std::shared_ptr<Grid2D<T>> project_z(ProjectionMode mode =
																							 ProjectionMode::MAX) const
				{
					if (sparse_data_->empty()) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid2D<T>>(data, 1, 1);
					}

					// Get all unique X and Y coordinates
					auto unique_x = sparse_data_->get_unique_x_coords();
					auto unique_y = sparse_data_->get_unique_y_coords();

					if (unique_x.empty() || unique_y.empty()) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid2D<T>>(data, 1, 1);
					}

					std::sort(unique_x.begin(), unique_x.end());
					std::sort(unique_y.begin(), unique_y.end());

					int64_t min_x = unique_x.front(), max_x = unique_x.back();
					int64_t min_y = unique_y.front(), max_y = unique_y.back();

					size_t width = static_cast<size_t>(max_x - min_x + 1);
					size_t height = static_cast<size_t>(max_y - min_y + 1);

					std::vector<T> grid_data(width * height);
					std::vector<bool> has_data(width * height, false);

					// Initialize based on projection mode
					T init_value = sparse_data_->get_default_value();
					if (mode == ProjectionMode::MAX) {
						std::fill(grid_data.begin(),
											grid_data.end(),
											std::numeric_limits<T>::lowest());
					} else if (mode == ProjectionMode::MIN) {
						std::fill(grid_data.begin(),
											grid_data.end(),
											std::numeric_limits<T>::max());
					} else {
						std::fill(grid_data.begin(), grid_data.end(), init_value);
					}

					// Project values
					std::unordered_map<size_t, std::vector<T>> projection_data;
					std::unordered_map<size_t, size_t> count_data;

					for (const auto& [coord, voxel] : *sparse_data_) {
						size_t grid_x = static_cast<size_t>(coord.x - min_x);
						size_t grid_y = static_cast<size_t>(coord.y - min_y);
						size_t grid_idx = grid_y * width + grid_x;

						T value = voxel.value;

						switch (mode) {
							case ProjectionMode::MAX:
								if (!has_data[grid_idx] || value > grid_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::MIN:
								if (!has_data[grid_idx] || value < grid_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::MEAN:
							case ProjectionMode::SUM:
								projection_data[grid_idx].push_back(value);
								break;

							case ProjectionMode::COUNT:
								count_data[grid_idx]++;
								break;

							case ProjectionMode::FIRST:
								if (!has_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::LAST:
								grid_data[grid_idx] = value;
								break;
						}

						has_data[grid_idx] = true;
					}

					// Finalize projection calculations
					if (mode == ProjectionMode::MEAN) {
						for (const auto& [idx, values] : projection_data) {
							if (!values.empty()) {
								T sum = std::accumulate(values.begin(), values.end(), T{});
								grid_data[idx] = sum / static_cast<T>(values.size());
							}
						}
					} else if (mode == ProjectionMode::SUM) {
						for (const auto& [idx, values] : projection_data) {
							if (!values.empty()) {
								grid_data[idx] =
									std::accumulate(values.begin(), values.end(), T{});
							}
						}
					} else if (mode == ProjectionMode::COUNT) {
						for (const auto& [idx, count] : count_data) {
							grid_data[idx] = static_cast<T>(count);
						}
					}

					// Set default values for cells without data
					for (size_t i = 0; i < grid_data.size(); ++i) {
						if (!has_data[i]) {
							grid_data[i] = init_value;
						}
					}

					return std::make_shared<Grid2D<T>>(grid_data, height, width);
				}

				/**
				 * Project along Y axis (front view)
				 */
				std::shared_ptr<Grid2D<T>> project_y(ProjectionMode mode =
																							 ProjectionMode::MAX) const
				{
					if (sparse_data_->empty()) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid2D<T>>(data, 1, 1);
					}

					auto unique_x = sparse_data_->get_unique_x_coords();
					auto unique_z = sparse_data_->get_unique_z_coords();

					if (unique_x.empty() || unique_z.empty()) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid2D<T>>(data, 1, 1);
					}

					std::sort(unique_x.begin(), unique_x.end());
					std::sort(unique_z.begin(), unique_z.end());

					int64_t min_x = unique_x.front(), max_x = unique_x.back();
					int64_t min_z = unique_z.front(), max_z = unique_z.back();

					size_t width = static_cast<size_t>(max_x - min_x + 1);
					size_t depth = static_cast<size_t>(max_z - min_z + 1);

					std::vector<T> grid_data(width * depth);
					std::vector<bool> has_data(width * depth, false);

					// Similar projection logic as project_z but for Y axis
					T init_value = sparse_data_->get_default_value();
					if (mode == ProjectionMode::MAX) {
						std::fill(grid_data.begin(),
											grid_data.end(),
											std::numeric_limits<T>::lowest());
					} else if (mode == ProjectionMode::MIN) {
						std::fill(grid_data.begin(),
											grid_data.end(),
											std::numeric_limits<T>::max());
					} else {
						std::fill(grid_data.begin(), grid_data.end(), init_value);
					}

					std::unordered_map<size_t, std::vector<T>> projection_data;
					std::unordered_map<size_t, size_t> count_data;

					for (const auto& [coord, voxel] : *sparse_data_) {
						size_t grid_x = static_cast<size_t>(coord.x - min_x);
						size_t grid_z = static_cast<size_t>(coord.z - min_z);
						size_t grid_idx = grid_z * width + grid_x;

						T value = voxel.value;

						switch (mode) {
							case ProjectionMode::MAX:
								if (!has_data[grid_idx] || value > grid_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::MIN:
								if (!has_data[grid_idx] || value < grid_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::MEAN:
							case ProjectionMode::SUM:
								projection_data[grid_idx].push_back(value);
								break;

							case ProjectionMode::COUNT:
								count_data[grid_idx]++;
								break;

							case ProjectionMode::FIRST:
								if (!has_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::LAST:
								grid_data[grid_idx] = value;
								break;
						}

						has_data[grid_idx] = true;
					}

					// Finalize calculations
					if (mode == ProjectionMode::MEAN) {
						for (const auto& [idx, values] : projection_data) {
							if (!values.empty()) {
								T sum = std::accumulate(values.begin(), values.end(), T{});
								grid_data[idx] = sum / static_cast<T>(values.size());
							}
						}
					} else if (mode == ProjectionMode::SUM) {
						for (const auto& [idx, values] : projection_data) {
							if (!values.empty()) {
								grid_data[idx] =
									std::accumulate(values.begin(), values.end(), T{});
							}
						}
					} else if (mode == ProjectionMode::COUNT) {
						for (const auto& [idx, count] : count_data) {
							grid_data[idx] = static_cast<T>(count);
						}
					}

					for (size_t i = 0; i < grid_data.size(); ++i) {
						if (!has_data[i]) {
							grid_data[i] = init_value;
						}
					}

					return std::make_shared<Grid2D<T>>(grid_data, depth, width);
				}

				/**
				 * Project along X axis (side view)
				 */
				std::shared_ptr<Grid2D<T>> project_x(ProjectionMode mode =
																							 ProjectionMode::MAX) const
				{
					if (sparse_data_->empty()) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid2D<T>>(data, 1, 1);
					}

					auto unique_y = sparse_data_->get_unique_y_coords();
					auto unique_z = sparse_data_->get_unique_z_coords();

					if (unique_y.empty() || unique_z.empty()) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid2D<T>>(data, 1, 1);
					}

					std::sort(unique_y.begin(), unique_y.end());
					std::sort(unique_z.begin(), unique_z.end());

					int64_t min_y = unique_y.front(), max_y = unique_y.back();
					int64_t min_z = unique_z.front(), max_z = unique_z.back();

					size_t height = static_cast<size_t>(max_y - min_y + 1);
					size_t depth = static_cast<size_t>(max_z - min_z + 1);

					std::vector<T> grid_data(height * depth);
					std::vector<bool> has_data(height * depth, false);

					T init_value = sparse_data_->get_default_value();
					if (mode == ProjectionMode::MAX) {
						std::fill(grid_data.begin(),
											grid_data.end(),
											std::numeric_limits<T>::lowest());
					} else if (mode == ProjectionMode::MIN) {
						std::fill(grid_data.begin(),
											grid_data.end(),
											std::numeric_limits<T>::max());
					} else {
						std::fill(grid_data.begin(), grid_data.end(), init_value);
					}

					std::unordered_map<size_t, std::vector<T>> projection_data;
					std::unordered_map<size_t, size_t> count_data;

					for (const auto& [coord, voxel] : *sparse_data_) {
						size_t grid_y = static_cast<size_t>(coord.y - min_y);
						size_t grid_z = static_cast<size_t>(coord.z - min_z);
						size_t grid_idx = grid_z * height + grid_y;

						T value = voxel.value;

						switch (mode) {
							case ProjectionMode::MAX:
								if (!has_data[grid_idx] || value > grid_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::MIN:
								if (!has_data[grid_idx] || value < grid_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::MEAN:
							case ProjectionMode::SUM:
								projection_data[grid_idx].push_back(value);
								break;

							case ProjectionMode::COUNT:
								count_data[grid_idx]++;
								break;

							case ProjectionMode::FIRST:
								if (!has_data[grid_idx]) {
									grid_data[grid_idx] = value;
								}
								break;

							case ProjectionMode::LAST:
								grid_data[grid_idx] = value;
								break;
						}

						has_data[grid_idx] = true;
					}

					// Finalize calculations
					if (mode == ProjectionMode::MEAN) {
						for (const auto& [idx, values] : projection_data) {
							if (!values.empty()) {
								T sum = std::accumulate(values.begin(), values.end(), T{});
								grid_data[idx] = sum / static_cast<T>(values.size());
							}
						}
					} else if (mode == ProjectionMode::SUM) {
						for (const auto& [idx, values] : projection_data) {
							if (!values.empty()) {
								grid_data[idx] =
									std::accumulate(values.begin(), values.end(), T{});
							}
						}
					} else if (mode == ProjectionMode::COUNT) {
						for (const auto& [idx, count] : count_data) {
							grid_data[idx] = static_cast<T>(count);
						}
					}

					for (size_t i = 0; i < grid_data.size(); ++i) {
						if (!has_data[i]) {
							grid_data[i] = init_value;
						}
					}

					return std::make_shared<Grid2D<T>>(grid_data, depth, height);
				}

				// ======================
				// INTERPOLATION AND SAMPLING
				// ======================

				/**
				 * Sample data at regular grid points with interpolation
				 */
				std::shared_ptr<Grid3D<T>> sample_regular_grid(
					const BoundingBox3D& bbox,
					size_t resolution_x,
					size_t resolution_y,
					size_t resolution_z,
					bool use_interpolation = true) const
				{

					if (!bbox.is_valid) {
						std::vector<T> data = { sparse_data_->get_default_value() };
						return std::make_shared<Grid3D<T>>(data, 1, 1, 1);
					}

					std::vector<T> grid_data(resolution_x * resolution_y * resolution_z,
																	 sparse_data_->get_default_value());

					double step_x =
						(resolution_x > 1)
							? static_cast<double>(bbox.max_coord.x - bbox.min_coord.x) /
									(resolution_x - 1)
							: 0.0;
					double step_y =
						(resolution_y > 1)
							? static_cast<double>(bbox.max_coord.y - bbox.min_coord.y) /
									(resolution_y - 1)
							: 0.0;
					double step_z =
						(resolution_z > 1)
							? static_cast<double>(bbox.max_coord.z - bbox.min_coord.z) /
									(resolution_z - 1)
							: 0.0;

					for (size_t z = 0; z < resolution_z; ++z) {
						for (size_t y = 0; y < resolution_y; ++y) {
							for (size_t x = 0; x < resolution_x; ++x) {
								double world_x = bbox.min_coord.x + x * step_x;
								double world_y = bbox.min_coord.y + y * step_y;
								double world_z = bbox.min_coord.z + z * step_z;

								T value;
								if (use_interpolation) {
									value = interpolate_at_point(world_x, world_y, world_z);
								} else {
									Coord3D nearest_coord(
										static_cast<int64_t>(std::round(world_x)),
										static_cast<int64_t>(std::round(world_y)),
										static_cast<int64_t>(std::round(world_z)));
									value = sparse_data_->get_value(nearest_coord);
								}

								size_t grid_idx =
									z * resolution_y * resolution_x + y * resolution_x + x;
								grid_data[grid_idx] = value;
							}
						}
					}

					return std::make_shared<Grid3D<T>>(
						grid_data, resolution_z, resolution_y, resolution_x);
				}

				/**
				 * Trilinear interpolation at a point
				 */
				T interpolate_at_point(double x, double y, double z) const
				{
					// Get the 8 surrounding integer coordinates
					int64_t x0 = static_cast<int64_t>(std::floor(x));
					int64_t y0 = static_cast<int64_t>(std::floor(y));
					int64_t z0 = static_cast<int64_t>(std::floor(z));
					int64_t x1 = x0 + 1;
					int64_t y1 = y0 + 1;
					int64_t z1 = z0 + 1;

					// Get fractional parts
					double fx = x - x0;
					double fy = y - y0;
					double fz = z - z0;

					// Get values at 8 corners
					T v000 = sparse_data_->get_value(Coord3D(x0, y0, z0));
					T v001 = sparse_data_->get_value(Coord3D(x0, y0, z1));
					T v010 = sparse_data_->get_value(Coord3D(x0, y1, z0));
					T v011 = sparse_data_->get_value(Coord3D(x0, y1, z1));
					T v100 = sparse_data_->get_value(Coord3D(x1, y0, z0));
					T v101 = sparse_data_->get_value(Coord3D(x1, y0, z1));
					T v110 = sparse_data_->get_value(Coord3D(x1, y1, z0));
					T v111 = sparse_data_->get_value(Coord3D(x1, y1, z1));

					// Trilinear interpolation
					T c00 = v000 * (1 - fx) + v100 * fx;
					T c01 = v001 * (1 - fx) + v101 * fx;
					T c10 = v010 * (1 - fx) + v110 * fx;
					T c11 = v011 * (1 - fx) + v111 * fx;

					T c0 = c00 * (1 - fy) + c10 * fy;
					T c1 = c01 * (1 - fy) + c11 * fy;

					return c0 * (1 - fz) + c1 * fz;
				}

				// ======================
				// UTILITIES
				// ======================

				/**
				 * Get the underlying sparse container
				 */
				std::shared_ptr<Sparse3DContainer<T>> get_sparse_container() const
				{
					return sparse_data_;
				}

				/**
				 * Create a new interface with different sparse data
				 */
				void set_sparse_container(
					std::shared_ptr<Sparse3DContainer<T>> new_data)
				{
					if (!new_data) {
						throw std::invalid_argument("Sparse data container cannot be null");
					}
					sparse_data_ = new_data;
				}
			};

			// Type aliases for common use cases
			using Sparse3DGridInterfaceFloat = Sparse3DGridInterface<float>;
			using Sparse3DGridInterfaceDouble = Sparse3DGridInterface<double>;
			using Sparse3DGridInterfaceInt = Sparse3DGridInterface<int32_t>;

		} // namespace dagger2
