#pragma once

#include "dg2_sparse3d_grid_interface.hpp"
#include <algorithm>
#include <chrono>
#include <functional>
#include <memory>
#include <queue>
#include <random>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * Sparse3DAlgorithmsBasic: Basic algorithms for sparse 3D data processing
 *
 * This class provides fundamental algorithms including flood fill, connected
 * components, morphological operations, and distance transforms.
 */
template<typename T>
class Sparse3DAlgorithmsBasic
{
public:
	using SparseContainer = Sparse3DContainer<T>;
	using GridInterface = Sparse3DGridInterface<T>;

	// ======================
	// FLOOD FILL ALGORITHMS
	// ======================

	/**
	 * Flood fill starting from a seed point
	 */
	static std::unordered_set<Coord3D> flood_fill(
		const SparseContainer& container,
		const Coord3D& seed,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		std::unordered_set<Coord3D> filled;
		if (!container.exists(seed) || !condition(container.get_value(seed))) {
			return filled;
		}

		std::queue<Coord3D> queue;
		queue.push(seed);
		filled.insert(seed);

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			// Check all neighbors
			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current + offset;

				if (filled.find(neighbor) == filled.end() &&
						container.exists(neighbor) &&
						condition(container.get_value(neighbor))) {

					filled.insert(neighbor);
					queue.push(neighbor);
				}
			}
		}

		return filled;
	}

	/**
	 * Flood fill with value replacement
	 */
	static void flood_fill_replace(
		SparseContainer& container,
		const Coord3D& seed,
		const T& new_value,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto filled_coords = flood_fill(container, seed, condition, connectivity);

		for (const auto& coord : filled_coords) {
			container.set_value(coord, new_value);
		}
	}

	/**
	 * Flood fill with custom operation
	 */
	static void flood_fill_operation(
		SparseContainer& container,
		const Coord3D& seed,
		std::function<T(const T&)> operation,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto filled_coords = flood_fill(container, seed, condition, connectivity);

		for (const auto& coord : filled_coords) {
			T current_value = container.get_value(coord);
			T new_value = operation(current_value);
			container.set_value(coord, new_value);
		}
	}

	/**
	 * Multi-seed flood fill
	 */
	static std::unordered_set<Coord3D> flood_fill_multi_seed(
		const SparseContainer& container,
		const std::vector<Coord3D>& seeds,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		std::unordered_set<Coord3D> filled;
		std::queue<Coord3D> queue;

		// Initialize with all valid seeds
		for (const auto& seed : seeds) {
			if (container.exists(seed) && condition(container.get_value(seed))) {
				filled.insert(seed);
				queue.push(seed);
			}
		}

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current + offset;

				if (filled.find(neighbor) == filled.end() &&
						container.exists(neighbor) &&
						condition(container.get_value(neighbor))) {

					filled.insert(neighbor);
					queue.push(neighbor);
				}
			}
		}

		return filled;
	}

	// ======================
	// CONNECTED COMPONENTS
	// ======================

	/**
	 * Find all connected components
	 */
	static std::vector<std::unordered_set<Coord3D>> find_connected_components(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		std::vector<std::unordered_set<Coord3D>> components;
		std::unordered_set<Coord3D> visited;

		// Iterate through all voxels
		for (const auto& [coord, voxel] : container) {
			if (visited.find(coord) == visited.end() && condition(voxel.value)) {
				auto component = flood_fill(container, coord, condition, connectivity);

				if (!component.empty()) {
					components.push_back(component);
					visited.insert(component.begin(), component.end());
				}
			}
		}

		return components;
	}

	/**
	 * Find largest connected component
	 */
	static std::unordered_set<Coord3D> find_largest_component(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto components =
			find_connected_components(container, condition, connectivity);

		if (components.empty()) {
			return std::unordered_set<Coord3D>();
		}

		auto largest = std::max_element(
			components.begin(), components.end(), [](const auto& a, const auto& b) {
				return a.size() < b.size();
			});

		return *largest;
	}

	/**
	 * Label connected components
	 */
	static std::shared_ptr<SparseContainer> label_components(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto labeled = std::make_shared<SparseContainer>(static_cast<T>(0));
		auto components =
			find_connected_components(container, condition, connectivity);

		T label = static_cast<T>(1);
		for (const auto& component : components) {
			for (const auto& coord : component) {
				labeled->set_value(coord, label);
			}
			++label;
		}

		return labeled;
	}

	/**
	 * Filter components by size
	 */
	static std::shared_ptr<SparseContainer> filter_components_by_size(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		size_t min_size,
		size_t max_size = std::numeric_limits<size_t>::max(),
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());
		auto components =
			find_connected_components(container, condition, connectivity);

		for (const auto& component : components) {
			if (component.size() >= min_size && component.size() <= max_size) {
				for (const auto& coord : component) {
					result->set_value(coord, container.get_value(coord));
				}
			}
		}

		return result;
	}

	// ======================
	// MORPHOLOGICAL OPERATIONS
	// ======================

	/**
	 * Erosion operation
	 */
	static std::shared_ptr<SparseContainer> erosion(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6,
		size_t iterations = 1)
	{

		auto result = std::make_shared<SparseContainer>(container);
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		for (size_t iter = 0; iter < iterations; ++iter) {
			std::vector<Coord3D> to_remove;

			for (const auto& [coord, voxel] : *result) {
				if (condition(voxel.value)) {
					// Check if all neighbors satisfy condition
					bool all_neighbors_valid = true;

					for (Direction3D dir : directions) {
						Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
						Coord3D neighbor = coord + offset;

						if (!result->exists(neighbor) ||
								!condition(result->get_value(neighbor))) {
							all_neighbors_valid = false;
							break;
						}
					}

					if (!all_neighbors_valid) {
						to_remove.push_back(coord);
					}
				}
			}

			for (const auto& coord : to_remove) {
				result->remove(coord);
			}
		}

		return result;
	}

	/**
	 * Dilation operation
	 */
	static std::shared_ptr<SparseContainer> dilation(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		const T& fill_value,
		Connectivity3D connectivity = Connectivity3D::C6,
		size_t iterations = 1)
	{

		auto result = std::make_shared<SparseContainer>(container);
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		for (size_t iter = 0; iter < iterations; ++iter) {
			std::vector<Coord3D> to_add;

			for (const auto& [coord, voxel] : *result) {
				if (condition(voxel.value)) {
					// Add all valid neighbors
					for (Direction3D dir : directions) {
						Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
						Coord3D neighbor = coord + offset;

						if (!result->exists(neighbor)) {
							to_add.push_back(neighbor);
						}
					}
				}
			}

			for (const auto& coord : to_add) {
				result->set_value(coord, fill_value);
			}
		}

		return result;
	}

	/**
	 * Opening operation (erosion followed by dilation)
	 */
	static std::shared_ptr<SparseContainer> opening(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		const T& fill_value,
		Connectivity3D connectivity = Connectivity3D::C6,
		size_t iterations = 1)
	{

		auto eroded = erosion(container, condition, connectivity, iterations);
		return dilation(*eroded, condition, fill_value, connectivity, iterations);
	}

	/**
	 * Closing operation (dilation followed by erosion)
	 */
	static std::shared_ptr<SparseContainer> closing(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		const T& fill_value,
		Connectivity3D connectivity = Connectivity3D::C6,
		size_t iterations = 1)
	{

		auto dilated =
			dilation(container, condition, fill_value, connectivity, iterations);
		return erosion(*dilated, condition, connectivity, iterations);
	}

	/**
	 * Custom structuring element morphology
	 */
	static std::shared_ptr<SparseContainer> morphology_custom_kernel(
		const SparseContainer& container,
		const std::vector<Coord3D>& kernel_offsets,
		std::function<bool(const T&)> condition,
		bool is_erosion = true)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());

		for (const auto& [coord, voxel] : container) {
			if (!condition(voxel.value))
				continue;

			bool operation_result = true;

			for (const auto& offset : kernel_offsets) {
				Coord3D test_coord = coord + offset;
				bool neighbor_satisfies = container.exists(test_coord) &&
																	condition(container.get_value(test_coord));

				if (is_erosion) {
					// For erosion, all kernel positions must be satisfied
					if (!neighbor_satisfies) {
						operation_result = false;
						break;
					}
				} else {
					// For dilation, any kernel position satisfying condition is enough
					if (neighbor_satisfies) {
						operation_result = true;
						break;
					}
				}
			}

			if (operation_result) {
				result->set_value(coord, voxel.value);
			}
		}

		return result;
	}

	// ======================
	// DISTANCE TRANSFORMS
	// ======================

	/**
	 * Euclidean distance transform
	 */
	static std::shared_ptr<Sparse3DContainer<double>>
	distance_transform_euclidean(
		const SparseContainer& container,
		std::function<bool(const T&)> foreground_condition)
	{

		auto result = std::make_shared<Sparse3DContainer<double>>(
			std::numeric_limits<double>::infinity());

		// Initialize: set foreground voxels to 0, collect them
		std::queue<Coord3D> queue;
		for (const auto& [coord, voxel] : container) {
			if (foreground_condition(voxel.value)) {
				result->set_value(coord, 0.0);
				queue.push(coord);
			}
		}

		// Propagate distances using Dijkstra-like approach
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C26);

		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			double current_dist = result->get_value(current);

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current + offset;
				double edge_length = Direction3DUtils::get_distance(dir);
				double new_dist = current_dist + edge_length;

				// Only process if this gives a shorter distance
				if (new_dist < result->get_value(neighbor)) {
					result->set_value(neighbor, new_dist);
					queue.push(neighbor);
				}
			}
		}

		return result;
	}

	/**
	 * Manhattan distance transform
	 */
	static std::shared_ptr<Sparse3DContainer<int64_t>>
	distance_transform_manhattan(
		const SparseContainer& container,
		std::function<bool(const T&)> foreground_condition)
	{

		auto result = std::make_shared<Sparse3DContainer<int64_t>>(
			std::numeric_limits<int64_t>::max());

		// Initialize
		std::queue<Coord3D> queue;
		for (const auto& [coord, voxel] : container) {
			if (foreground_condition(voxel.value)) {
				result->set_value(coord, 0);
				queue.push(coord);
			}
		}

		// Propagate using 6-connectivity only
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C6);

		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			int64_t current_dist = result->get_value(current);

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current + offset;
				int64_t new_dist = current_dist + 1;

				if (new_dist < result->get_value(neighbor)) {
					result->set_value(neighbor, new_dist);
					queue.push(neighbor);
				}
			}
		}

		return result;
	}

	/**
	 * Chamfer distance transform (3-4-5 distances)
	 */
	static std::shared_ptr<Sparse3DContainer<int32_t>> distance_transform_chamfer(
		const SparseContainer& container,
		std::function<bool(const T&)> foreground_condition)
	{

		auto result = std::make_shared<Sparse3DContainer<int32_t>>(
			std::numeric_limits<int32_t>::max());

		// Chamfer distances: face=3, edge=4, corner=5
		std::vector<std::pair<Coord3D, int32_t>> chamfer_offsets = {
			// Face neighbors (distance 3)
			{ { 1, 0, 0 }, 3 },
			{ { -1, 0, 0 }, 3 },
			{ { 0, 1, 0 }, 3 },
			{ { 0, -1, 0 }, 3 },
			{ { 0, 0, 1 }, 3 },
			{ { 0, 0, -1 }, 3 },

			// Edge neighbors (distance 4)
			{ { 1, 1, 0 }, 4 },
			{ { 1, -1, 0 }, 4 },
			{ { -1, 1, 0 }, 4 },
			{ { -1, -1, 0 }, 4 },
			{ { 1, 0, 1 }, 4 },
			{ { 1, 0, -1 }, 4 },
			{ { -1, 0, 1 }, 4 },
			{ { -1, 0, -1 }, 4 },
			{ { 0, 1, 1 }, 4 },
			{ { 0, 1, -1 }, 4 },
			{ { 0, -1, 1 }, 4 },
			{ { 0, -1, -1 }, 4 },

			// Corner neighbors (distance 5)
			{ { 1, 1, 1 }, 5 },
			{ { 1, 1, -1 }, 5 },
			{ { 1, -1, 1 }, 5 },
			{ { 1, -1, -1 }, 5 },
			{ { -1, 1, 1 }, 5 },
			{ { -1, 1, -1 }, 5 },
			{ { -1, -1, 1 }, 5 },
			{ { -1, -1, -1 }, 5 }
		};

		// Initialize
		std::queue<Coord3D> queue;
		for (const auto& [coord, voxel] : container) {
			if (foreground_condition(voxel.value)) {
				result->set_value(coord, 0);
				queue.push(coord);
			}
		}

		// Propagate
		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			int32_t current_dist = result->get_value(current);

			for (const auto& [offset, dist] : chamfer_offsets) {
				Coord3D neighbor = current + offset;
				int32_t new_dist = current_dist + dist;

				if (new_dist < result->get_value(neighbor)) {
					result->set_value(neighbor, new_dist);
					queue.push(neighbor);
				}
			}
		}

		return result;
	}

	// ======================
	// SURFACE OPERATIONS
	// ======================

	/**
	 * Extract surface voxels (voxels with at least one background neighbor)
	 */
	static std::shared_ptr<SparseContainer> extract_surface(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		for (const auto& [coord, voxel] : container) {
			if (!condition(voxel.value))
				continue;

			bool is_surface = false;

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = coord + offset;

				if (!container.exists(neighbor) ||
						!condition(container.get_value(neighbor))) {
					is_surface = true;
					break;
				}
			}

			if (is_surface) {
				result->set_value(coord, voxel.value);
			}
		}

		return result;
	}

	/**
	 * Extract interior voxels (voxels with all foreground neighbors)
	 */
	static std::shared_ptr<SparseContainer> extract_interior(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		for (const auto& [coord, voxel] : container) {
			if (!condition(voxel.value))
				continue;

			bool is_interior = true;

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = coord + offset;

				if (!container.exists(neighbor) ||
						!condition(container.get_value(neighbor))) {
					is_interior = false;
					break;
				}
			}

			if (is_interior) {
				result->set_value(coord, voxel.value);
			}
		}

		return result;
	}

	/**
	 * Extract boundary voxels (voxels on the boundary of the bounding box)
	 */
	static std::shared_ptr<SparseContainer> extract_boundary(
		const SparseContainer& container,
		std::function<bool(const T&)> condition)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());
		auto bbox = container.get_bounding_box();

		if (!bbox.is_valid)
			return result;

		for (const auto& [coord, voxel] : container) {
			if (!condition(voxel.value))
				continue;

			bool is_boundary =
				(coord.x == bbox.min_coord.x || coord.x == bbox.max_coord.x ||
				 coord.y == bbox.min_coord.y || coord.y == bbox.max_coord.y ||
				 coord.z == bbox.min_coord.z || coord.z == bbox.max_coord.z);

			if (is_boundary) {
				result->set_value(coord, voxel.value);
			}
		}

		return result;
	}
};

// Type aliases for common use cases
using Sparse3DAlgorithmsBasicFloat = Sparse3DAlgorithmsBasic<float>;
using Sparse3DAlgorithmsBasicDouble = Sparse3DAlgorithmsBasic<double>;
using Sparse3DAlgorithmsBasicInt = Sparse3DAlgorithmsBasic<int32_t>;

/**
 * Sparse3DAlgorithmsAdvanced: Advanced algorithms for sparse 3D data processing
 *
 * This class provides sophisticated algorithms including filtering,
 * skeletonization, pathfinding, geometric analysis, and specialized processing
 * operations.
 */
template<typename T>
class Sparse3DAlgorithmsAdvanced
{
public:
	using SparseContainer = Sparse3DContainer<T>;
	using GridInterface = Sparse3DGridInterface<T>;
	using BasicAlgorithms = Sparse3DAlgorithmsBasic<T>;

	// ======================
	// FILTERING OPERATIONS
	// ======================

	/**
	 * 3D Gaussian filter (approximated using separable kernels)
	 */
	static std::shared_ptr<SparseContainer> gaussian_filter(
		const SparseContainer& container,
		double sigma,
		int kernel_size = -1)
	{

		if (kernel_size <= 0) {
			kernel_size = static_cast<int>(std::ceil(6 * sigma));
			if (kernel_size % 2 == 0)
				kernel_size++;
		}

		// Generate 1D Gaussian kernel
		std::vector<double> kernel(kernel_size);
		int center = kernel_size / 2;
		double sum = 0.0;

		for (int i = 0; i < kernel_size; ++i) {
			int offset = i - center;
			kernel[i] = std::exp(-(offset * offset) / (2 * sigma * sigma));
			sum += kernel[i];
		}

		// Normalize kernel
		for (auto& k : kernel) {
			k /= sum;
		}

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());

		// Apply separable filtering
		for (const auto& [coord, voxel] : container) {
			double filtered_value = 0.0;
			double weight_sum = 0.0;

			// Apply 3D kernel (simplified, not truly separable)
			for (int dx = -center; dx <= center; ++dx) {
				for (int dy = -center; dy <= center; ++dy) {
					for (int dz = -center; dz <= center; ++dz) {
						Coord3D sample_coord = coord + Coord3D(dx, dy, dz);
						double weight =
							kernel[dx + center] * kernel[dy + center] * kernel[dz + center];

						T sample_value = container.get_value(sample_coord);
						filtered_value += static_cast<double>(sample_value) * weight;
						weight_sum += weight;
					}
				}
			}

			if (weight_sum > 0) {
				result->set_value(coord, static_cast<T>(filtered_value / weight_sum));
			}
		}

		return result;
	}

	/**
	 * Median filter
	 */
	static std::shared_ptr<SparseContainer> median_filter(
		const SparseContainer& container,
		int radius = 1,
		Connectivity3D connectivity = Connectivity3D::C26)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());

		for (const auto& [coord, voxel] : container) {
			std::vector<T> values;

			// Collect values in neighborhood
			for (int dx = -radius; dx <= radius; ++dx) {
				for (int dy = -radius; dy <= radius; ++dy) {
					for (int dz = -radius; dz <= radius; ++dz) {
						if (connectivity == Connectivity3D::C6 &&
								(std::abs(dx) + std::abs(dy) + std::abs(dz)) > 1) {
							continue;
						}
						if (connectivity == Connectivity3D::C18 &&
								std::max({ std::abs(dx), std::abs(dy), std::abs(dz) }) > 1) {
							continue;
						}

						Coord3D sample_coord = coord + Coord3D(dx, dy, dz);
						values.push_back(container.get_value(sample_coord));
					}
				}
			}

			if (!values.empty()) {
				std::sort(values.begin(), values.end());
				T median_value = values[values.size() / 2];
				result->set_value(coord, median_value);
			}
		}

		return result;
	}

	/**
	 * Bilateral filter (edge-preserving)
	 */
	static std::shared_ptr<SparseContainer> bilateral_filter(
		const SparseContainer& container,
		double spatial_sigma,
		double intensity_sigma,
		int radius = 2)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());

		for (const auto& [coord, voxel] : container) {
			double filtered_value = 0.0;
			double weight_sum = 0.0;
			T center_value = voxel.value;

			// Apply bilateral kernel
			for (int dx = -radius; dx <= radius; ++dx) {
				for (int dy = -radius; dy <= radius; ++dy) {
					for (int dz = -radius; dz <= radius; ++dz) {
						Coord3D sample_coord = coord + Coord3D(dx, dy, dz);
						T sample_value = container.get_value(sample_coord);

						// Spatial weight
						double spatial_dist = std::sqrt(dx * dx + dy * dy + dz * dz);
						double spatial_weight =
							std::exp(-(spatial_dist * spatial_dist) /
											 (2 * spatial_sigma * spatial_sigma));

						// Intensity weight
						double intensity_diff = static_cast<double>(sample_value) -
																		static_cast<double>(center_value);
						double intensity_weight =
							std::exp(-(intensity_diff * intensity_diff) /
											 (2 * intensity_sigma * intensity_sigma));

						double total_weight = spatial_weight * intensity_weight;
						filtered_value += static_cast<double>(sample_value) * total_weight;
						weight_sum += total_weight;
					}
				}
			}

			if (weight_sum > 0) {
				result->set_value(coord, static_cast<T>(filtered_value / weight_sum));
			}
		}

		return result;
	}

	/**
	 * Anisotropic diffusion filter
	 */
	static std::shared_ptr<SparseContainer> anisotropic_diffusion(
		const SparseContainer& container,
		double kappa,
		double lambda = 0.25,
		int iterations = 10)
	{

		auto result = std::make_shared<SparseContainer>(container);

		// Face neighbor offsets for gradient computation
		std::vector<Coord3D> face_offsets = { { 1, 0, 0 }, { -1, 0, 0 },
																					{ 0, 1, 0 }, { 0, -1, 0 },
																					{ 0, 0, 1 }, { 0, 0, -1 } };

		for (int iter = 0; iter < iterations; ++iter) {
			auto temp = std::make_shared<SparseContainer>(*result);

			for (const auto& [coord, voxel] : *temp) {
				double center_value = static_cast<double>(voxel.value);
				double diffusion_sum = 0.0;

				// Compute diffusion for each face neighbor
				for (const auto& offset : face_offsets) {
					Coord3D neighbor = coord + offset;
					double neighbor_value =
						static_cast<double>(temp->get_value(neighbor));

					double gradient = neighbor_value - center_value;
					double gradient_mag = std::abs(gradient);

					// Diffusion coefficient (edge-stopping function)
					double c = std::exp(-(gradient_mag / kappa) * (gradient_mag / kappa));

					diffusion_sum += c * gradient;
				}

				double new_value = center_value + lambda * diffusion_sum;
				result->set_value(coord, static_cast<T>(new_value));
			}
		}

		return result;
	}

	// ======================
	// SKELETONIZATION & THINNING
	// ======================

	/**
	 * Simple 3D thinning algorithm
	 */
	static std::shared_ptr<SparseContainer> skeletonize(
		const SparseContainer& container,
		std::function<bool(const T&)> condition)
	{

		auto result = std::make_shared<SparseContainer>(container);
		bool changed = true;

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C26);

		while (changed) {
			changed = false;
			std::vector<Coord3D> to_remove;

			for (const auto& [coord, voxel] : *result) {
				if (!condition(voxel.value))
					continue;

				// Count foreground neighbors
				int fg_neighbors = 0;
				std::vector<Coord3D> neighbor_coords;

				for (Direction3D dir : directions) {
					Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
					Coord3D neighbor = coord + offset;
					neighbor_coords.push_back(neighbor);

					if (result->exists(neighbor) &&
							condition(result->get_value(neighbor))) {
						fg_neighbors++;
					}
				}

				// Simple thinning condition
				if (fg_neighbors >= 2 && fg_neighbors <= 6) {
					if (is_simple_point(*result, coord, condition, neighbor_coords)) {
						to_remove.push_back(coord);
						changed = true;
					}
				}
			}

			for (const auto& coord : to_remove) {
				result->remove(coord);
			}
		}

		return result;
	}

	/**
	 * Distance-based skeleton
	 */
	static std::shared_ptr<SparseContainer> distance_skeleton(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		double threshold_ratio = 0.8)
	{

		// Compute distance transform
		auto distance_map =
			BasicAlgorithms::distance_transform_euclidean(container, condition);

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());

		// Find local maxima in distance map
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C26);

		for (const auto& [coord, dist_data] : *distance_map) {
			double center_dist = dist_data.value;
			if (center_dist <= 0)
				continue;

			bool is_local_maximum = true;
			double max_neighbor_dist = 0;

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = coord + offset;
				double neighbor_dist = distance_map->get_value(neighbor);

				max_neighbor_dist = std::max(max_neighbor_dist, neighbor_dist);

				if (neighbor_dist > center_dist) {
					is_local_maximum = false;
					break;
				}
			}

			if (is_local_maximum ||
					center_dist >= threshold_ratio * max_neighbor_dist) {
				result->set_value(coord, container.get_value(coord));
			}
		}

		return result;
	}

	// ======================
	// PATHFINDING & NAVIGATION
	// ======================

	/**
	 * Find shortest path between two points using A*
	 */
	static std::vector<Coord3D> shortest_path(
		const SparseContainer& container,
		const Coord3D& start,
		const Coord3D& goal,
		std::function<bool(const T&)> walkable_condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		struct Node
		{
			Coord3D coord;
			double g_cost; // Distance from start
			double h_cost; // Heuristic distance to goal
			double f_cost() const { return g_cost + h_cost; }
			Coord3D parent;

			Node(const Coord3D& c, double g, double h, const Coord3D& p)
				: coord(c)
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
		std::unordered_set<Coord3D> closed_set;
		std::unordered_map<Coord3D, Node> nodes;

		if (!container.exists(start) ||
				!walkable_condition(container.get_value(start)) ||
				!container.exists(goal) ||
				!walkable_condition(container.get_value(goal))) {
			return {};
		}

		double h_start = start.distance_to(goal);
		open_set.emplace(start, 0.0, h_start, Coord3D(-1, -1, -1));
		nodes.emplace(start, Node(start, 0.0, h_start, Coord3D(-1, -1, -1)));

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		while (!open_set.empty()) {
			Node current = open_set.top();
			open_set.pop();

			if (current.coord == goal) {
				// Reconstruct path
				std::vector<Coord3D> path;
				Coord3D trace = goal;
				while (trace != Coord3D(-1, -1, -1)) {
					path.push_back(trace);
					trace = nodes[trace].parent;
				}
				std::reverse(path.begin(), path.end());
				return path;
			}

			closed_set.insert(current.coord);

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current.coord + offset;

				if (closed_set.count(neighbor) || !container.exists(neighbor) ||
						!walkable_condition(container.get_value(neighbor))) {
					continue;
				}

				double edge_cost = Direction3DUtils::get_distance(dir);
				double tentative_g = current.g_cost + edge_cost;
				double h = neighbor.distance_to(goal);

				auto it = nodes.find(neighbor);
				if (it == nodes.end() || tentative_g < it->second.g_cost) {
					nodes[neighbor] = Node(neighbor, tentative_g, h, current.coord);
					open_set.emplace(neighbor, tentative_g, h, current.coord);
				}
			}
		}

		return {}; // No path found
	}

	/**
	 * Generate random walk from a starting point
	 */
	static std::vector<Coord3D> random_walk(
		const SparseContainer& container,
		const Coord3D& start,
		size_t max_steps,
		std::function<bool(const T&)> valid_condition,
		Connectivity3D connectivity = Connectivity3D::C6,
		uint32_t seed = 42)
	{

		std::vector<Coord3D> path;
		if (!container.exists(start) ||
				!valid_condition(container.get_value(start))) {
			return path;
		}

		std::mt19937 rng(seed);
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		Coord3D current = start;
		path.push_back(current);

		for (size_t step = 0; step < max_steps; ++step) {
			std::vector<Coord3D> valid_neighbors;

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current + offset;

				if (container.exists(neighbor) &&
						valid_condition(container.get_value(neighbor))) {
					valid_neighbors.push_back(neighbor);
				}
			}

			if (valid_neighbors.empty())
				break;

			std::uniform_int_distribution<size_t> neighbor_dist(
				0, valid_neighbors.size() - 1);
			current = valid_neighbors[neighbor_dist(rng)];
			path.push_back(current);
		}

		return path;
	}

	/**
	 * Find multiple paths with different constraints
	 */
	static std::vector<std::vector<Coord3D>> find_multiple_paths(
		const SparseContainer& container,
		const Coord3D& start,
		const Coord3D& goal,
		std::function<bool(const T&)> walkable_condition,
		size_t num_paths = 3,
		double path_separation = 2.0,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		std::vector<std::vector<Coord3D>> paths;
		std::unordered_set<Coord3D> blocked_coords;

		for (size_t i = 0; i < num_paths; ++i) {
			// Create temporary container with blocked coordinates
			auto temp_container = std::make_shared<SparseContainer>(container);
			for (const auto& blocked : blocked_coords) {
				temp_container->remove(blocked);
			}

			auto path = shortest_path(
				*temp_container, start, goal, walkable_condition, connectivity);

			if (path.empty())
				break;

			paths.push_back(path);

			// Block coordinates around this path for next iteration
			for (const auto& coord : path) {
				for (int dx = -static_cast<int>(path_separation);
						 dx <= static_cast<int>(path_separation);
						 ++dx) {
					for (int dy = -static_cast<int>(path_separation);
							 dy <= static_cast<int>(path_separation);
							 ++dy) {
						for (int dz = -static_cast<int>(path_separation);
								 dz <= static_cast<int>(path_separation);
								 ++dz) {
							Coord3D block_coord = coord + Coord3D(dx, dy, dz);
							if (block_coord != start && block_coord != goal) {
								blocked_coords.insert(block_coord);
							}
						}
					}
				}
			}
		}

		return paths;
	}

	// ======================
	// GEOMETRIC ANALYSIS
	// ======================

	/**
	 * Compute convex hull (simplified 3D version)
	 */
	static std::vector<Coord3D> convex_hull(
		const SparseContainer& container,
		std::function<bool(const T&)> condition)
	{

		std::vector<Coord3D> points;
		for (const auto& [coord, voxel] : container) {
			if (condition(voxel.value)) {
				points.push_back(coord);
			}
		}

		if (points.size() < 4)
			return points;

		// Simplified extremal points approach
		std::vector<Coord3D> hull;

		auto min_x = std::min_element(
			points.begin(), points.end(), [](const Coord3D& a, const Coord3D& b) {
				return a.x < b.x;
			});
		auto max_x = std::max_element(
			points.begin(), points.end(), [](const Coord3D& a, const Coord3D& b) {
				return a.x < b.x;
			});
		auto min_y = std::min_element(
			points.begin(), points.end(), [](const Coord3D& a, const Coord3D& b) {
				return a.y < b.y;
			});
		auto max_y = std::max_element(
			points.begin(), points.end(), [](const Coord3D& a, const Coord3D& b) {
				return a.y < b.y;
			});
		auto min_z = std::min_element(
			points.begin(), points.end(), [](const Coord3D& a, const Coord3D& b) {
				return a.z < b.z;
			});
		auto max_z = std::max_element(
			points.begin(), points.end(), [](const Coord3D& a, const Coord3D& b) {
				return a.z < b.z;
			});

		hull = { *min_x, *max_x, *min_y, *max_y, *min_z, *max_z };

		// Remove duplicates
		std::sort(hull.begin(), hull.end());
		hull.erase(std::unique(hull.begin(), hull.end()), hull.end());

		return hull;
	}

	/**
	 * Compute centroid
	 */
	static Coord3D compute_centroid(const SparseContainer& container,
																	std::function<bool(const T&)> condition)
	{

		int64_t sum_x = 0, sum_y = 0, sum_z = 0;
		size_t count = 0;

		for (const auto& [coord, voxel] : container) {
			if (condition(voxel.value)) {
				sum_x += coord.x;
				sum_y += coord.y;
				sum_z += coord.z;
				count++;
			}
		}

		if (count == 0)
			return Coord3D(0, 0, 0);

		return Coord3D(sum_x / static_cast<int64_t>(count),
									 sum_y / static_cast<int64_t>(count),
									 sum_z / static_cast<int64_t>(count));
	}

	/**
	 * Compute weighted centroid
	 */
	static Coord3D compute_weighted_centroid(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		std::function<double(const T&)> weight_function)
	{

		double sum_x = 0, sum_y = 0, sum_z = 0;
		double total_weight = 0;

		for (const auto& [coord, voxel] : container) {
			if (condition(voxel.value)) {
				double weight = weight_function(voxel.value);
				sum_x += coord.x * weight;
				sum_y += coord.y * weight;
				sum_z += coord.z * weight;
				total_weight += weight;
			}
		}

		if (total_weight == 0)
			return Coord3D(0, 0, 0);

		return Coord3D(static_cast<int64_t>(sum_x / total_weight),
									 static_cast<int64_t>(sum_y / total_weight),
									 static_cast<int64_t>(sum_z / total_weight));
	}

	/**
	 * Compute principal component analysis
	 */
	struct PCAResult
	{
		std::array<std::array<double, 3>, 3> eigenvectors;
		std::array<double, 3> eigenvalues;
		std::array<double, 3> centroid;
	};

	static PCAResult compute_pca(const SparseContainer& container,
															 std::function<bool(const T&)> condition)
	{

		PCAResult result = {};

		// Compute centroid
		std::vector<std::array<double, 3>> points;
		double cx = 0, cy = 0, cz = 0;
		size_t count = 0;

		for (const auto& [coord, voxel] : container) {
			if (condition(voxel.value)) {
				points.push_back({ static_cast<double>(coord.x),
													 static_cast<double>(coord.y),
													 static_cast<double>(coord.z) });
				cx += coord.x;
				cy += coord.y;
				cz += coord.z;
				count++;
			}
		}

		if (count == 0)
			return result;

		cx /= count;
		cy /= count;
		cz /= count;

		result.centroid = { cx, cy, cz };

		// Compute covariance matrix
		std::array<std::array<double, 3>, 3> covariance = {
			{ { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }
		};

		for (const auto& point : points) {
			double dx = point[0] - cx;
			double dy = point[1] - cy;
			double dz = point[2] - cz;

			covariance[0][0] += dx * dx;
			covariance[0][1] += dx * dy;
			covariance[0][2] += dx * dz;
			covariance[1][1] += dy * dy;
			covariance[1][2] += dy * dz;
			covariance[2][2] += dz * dz;
		}

		// Symmetric matrix
		covariance[1][0] = covariance[0][1];
		covariance[2][0] = covariance[0][2];
		covariance[2][1] = covariance[1][2];

		// Normalize
		for (auto& row : covariance) {
			for (auto& val : row) {
				val /= count;
			}
		}

		// Simplified eigenvalue computation (diagonal approximation)
		result.eigenvalues = { covariance[0][0],
													 covariance[1][1],
													 covariance[2][2] };
		result.eigenvectors = { { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };

		return result;
	}

	/**
	 * Compute surface area approximation
	 */
	static double compute_surface_area(const SparseContainer& container,
																		 std::function<bool(const T&)> condition)
	{

		double surface_area = 0.0;
		auto face_offsets =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C6);

		for (const auto& [coord, voxel] : container) {
			if (!condition(voxel.value))
				continue;

			// Count exposed faces
			for (Direction3D dir : face_offsets) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = coord + offset;

				if (!container.exists(neighbor) ||
						!condition(container.get_value(neighbor))) {
					surface_area += 1.0; // Each face has area 1
				}
			}
		}

		return surface_area;
	}

	/**
	 * Compute volume
	 */
	static size_t compute_volume(const SparseContainer& container,
															 std::function<bool(const T&)> condition)
	{

		size_t volume = 0;
		for (const auto& [coord, voxel] : container) {
			if (condition(voxel.value)) {
				volume++;
			}
		}
		return volume;
	}

	// ======================
	// SPECIALIZED ALGORITHMS
	// ======================

	/**
	 * Watershed segmentation
	 */
	static std::shared_ptr<SparseContainer> watershed_segmentation(
		const SparseContainer& container,
		const std::vector<Coord3D>& seeds,
		std::function<double(const T&)> gradient_function)
	{

		auto result = std::make_shared<SparseContainer>(static_cast<T>(0));

		// Priority queue for watershed expansion
		struct WatershedPixel
		{
			Coord3D coord;
			double priority;
			T label;

			bool operator>(const WatershedPixel& other) const
			{
				return priority > other.priority;
			}
		};

		std::priority_queue<WatershedPixel,
												std::vector<WatershedPixel>,
												std::greater<WatershedPixel>>
			pq;
		std::unordered_set<Coord3D> processed;

		// Initialize seeds
		T current_label = static_cast<T>(1);
		for (const auto& seed : seeds) {
			if (container.exists(seed)) {
				result->set_value(seed, current_label);
				double gradient = gradient_function(container.get_value(seed));
				pq.push({ seed, gradient, current_label });
				processed.insert(seed);
				++current_label;
			}
		}

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C6);

		// Expand watersheds
		while (!pq.empty()) {
			WatershedPixel current = pq.top();
			pq.pop();

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current.coord + offset;

				if (processed.count(neighbor) || !container.exists(neighbor)) {
					continue;
				}

				double neighbor_gradient =
					gradient_function(container.get_value(neighbor));
				result->set_value(neighbor, current.label);
				pq.push({ neighbor, neighbor_gradient, current.label });
				processed.insert(neighbor);
			}
		}

		return result;
	}

	/**
	 * Region growing
	 */
	static std::shared_ptr<SparseContainer> region_growing(
		const SparseContainer& container,
		const Coord3D& seed,
		std::function<bool(const T&, const T&)> similarity_function,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());

		if (!container.exists(seed)) {
			return result;
		}

		std::queue<Coord3D> queue;
		std::unordered_set<Coord3D> processed;

		T seed_value = container.get_value(seed);
		queue.push(seed);
		processed.insert(seed);
		result->set_value(seed, seed_value);

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			T current_value = container.get_value(current);

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D neighbor = current + offset;

				if (processed.count(neighbor) || !container.exists(neighbor)) {
					continue;
				}

				T neighbor_value = container.get_value(neighbor);

				if (similarity_function(current_value, neighbor_value)) {
					result->set_value(neighbor, neighbor_value);
					queue.push(neighbor);
					processed.insert(neighbor);
				}
			}
		}

		return result;
	}

	/**
	 * Level set evolution (simplified)
	 */
	static std::shared_ptr<SparseContainer> level_set_evolution(
		const SparseContainer& container,
		std::function<bool(const T&)> initial_condition,
		std::function<double(const Coord3D&)> speed_function,
		double dt = 0.1,
		int iterations = 50)
	{

		// Create signed distance function
		auto distance_map = std::make_shared<Sparse3DContainer<double>>(1000.0);

		// Initialize level set
		for (const auto& [coord, voxel] : container) {
			if (initial_condition(voxel.value)) {
				distance_map->set_value(coord, -1.0); // Inside
			} else {
				distance_map->set_value(coord, 1.0); // Outside
			}
		}

		// Evolve level set
		for (int iter = 0; iter < iterations; ++iter) {
			auto new_distance_map =
				std::make_shared<Sparse3DContainer<double>>(*distance_map);

			for (const auto& [coord, dist_data] : *distance_map) {
				double current_dist = dist_data.value;
				double speed = speed_function(coord);

				// Simple upwind scheme for gradient computation
				double grad_x = 0, grad_y = 0, grad_z = 0;

				if (speed > 0) {
					// Forward differences
					grad_x =
						distance_map->get_value(coord + Coord3D(1, 0, 0)) - current_dist;
					grad_y =
						distance_map->get_value(coord + Coord3D(0, 1, 0)) - current_dist;
					grad_z =
						distance_map->get_value(coord + Coord3D(0, 0, 1)) - current_dist;
				} else {
					// Backward differences
					grad_x =
						current_dist - distance_map->get_value(coord + Coord3D(-1, 0, 0));
					grad_y =
						current_dist - distance_map->get_value(coord + Coord3D(0, -1, 0));
					grad_z =
						current_dist - distance_map->get_value(coord + Coord3D(0, 0, -1));
				}

				double grad_mag =
					std::sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);
				double new_dist = current_dist - dt * speed * grad_mag;

				new_distance_map->set_value(coord, new_dist);
			}

			distance_map = new_distance_map;
		}

		// Extract final segmentation
		auto result =
			std::make_shared<SparseContainer>(container.get_default_value());
		for (const auto& [coord, dist_data] : *distance_map) {
			if (dist_data.value < 0 && container.exists(coord)) {
				result->set_value(coord, container.get_value(coord));
			}
		}

		return result;
	}

private:
	// ======================
	// PRIVATE HELPER METHODS
	// ======================

	/**
	 * Check if a point is simple (can be removed without changing topology)
	 */
	static bool is_simple_point(const SparseContainer& container,
															const Coord3D& coord,
															std::function<bool(const T&)> condition,
															const std::vector<Coord3D>& neighbors)
	{

		std::unordered_set<Coord3D> fg_neighbors;
		for (const auto& neighbor : neighbors) {
			if (container.exists(neighbor) &&
					condition(container.get_value(neighbor))) {
				fg_neighbors.insert(neighbor);
			}
		}

		if (fg_neighbors.size() <= 1)
			return true;

		// Check connectivity among foreground neighbors
		std::unordered_set<Coord3D> visited;
		std::queue<Coord3D> queue;

		auto first_neighbor = *fg_neighbors.begin();
		queue.push(first_neighbor);
		visited.insert(first_neighbor);

		auto directions =
			Direction3DUtils::get_directions_for_connectivity(Connectivity3D::C26);

		while (!queue.empty()) {
			Coord3D current = queue.front();
			queue.pop();

			for (Direction3D dir : directions) {
				Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
				Coord3D next = current + offset;

				if (fg_neighbors.count(next) && !visited.count(next)) {
					visited.insert(next);
					queue.push(next);
				}
			}
		}

		return visited.size() == fg_neighbors.size();
	}
};

// ======================
// COMBINED ALGORITHM INTERFACE
// ======================

/**
 * Sparse3DAlgorithms: Complete algorithm collection combining basic and
 * advanced operations
 */
template<typename T>
class Sparse3DAlgorithms
	: public Sparse3DAlgorithmsBasic<T>
	, public Sparse3DAlgorithmsAdvanced<T>
{
public:
	using SparseContainer = Sparse3DContainer<T>;
	using GridInterface = Sparse3DGridInterface<T>;
	using BasicAlgorithms = Sparse3DAlgorithmsBasic<T>;
	using AdvancedAlgorithms = Sparse3DAlgorithmsAdvanced<T>;

	// ======================
	// COMBINED WORKFLOWS
	// ======================

	/**
	 * Complete segmentation pipeline
	 */
	static std::shared_ptr<SparseContainer> segment_objects(
		const SparseContainer& container,
		std::function<bool(const T&)> object_condition,
		size_t min_object_size = 10,
		size_t max_object_size = std::numeric_limits<size_t>::max(),
		bool apply_smoothing = true,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		auto result = std::make_shared<SparseContainer>(container);

		if (apply_smoothing) {
			result = AdvancedAlgorithms::gaussian_filter(*result, 1.0, 3);
		}

		auto objects = BasicAlgorithms::filter_components_by_size(*result,
																															object_condition,
																															min_object_size,
																															max_object_size,
																															connectivity);

		objects = BasicAlgorithms::closing(
			*objects, object_condition, static_cast<T>(1), connectivity, 1);

		return objects;
	}

	/**
	 * Surface analysis pipeline
	 */
	struct SurfaceAnalysis
	{
		std::shared_ptr<SparseContainer> surface_voxels;
		std::shared_ptr<SparseContainer> interior_voxels;
		std::shared_ptr<SparseContainer> skeleton;
		double surface_area;
		size_t volume;
		double surface_to_volume_ratio;
	};

	static SurfaceAnalysis analyze_surface(
		const SparseContainer& container,
		std::function<bool(const T&)> condition,
		Connectivity3D connectivity = Connectivity3D::C6)
	{

		SurfaceAnalysis analysis;

		analysis.surface_voxels =
			BasicAlgorithms::extract_surface(container, condition, connectivity);
		analysis.interior_voxels =
			BasicAlgorithms::extract_interior(container, condition, connectivity);
		analysis.skeleton = AdvancedAlgorithms::skeletonize(container, condition);
		analysis.surface_area =
			AdvancedAlgorithms::compute_surface_area(container, condition);
		analysis.volume = AdvancedAlgorithms::compute_volume(container, condition);
		analysis.surface_to_volume_ratio =
			(analysis.volume > 0) ? analysis.surface_area / analysis.volume : 0.0;

		return analysis;
	}

	/**
	 * Comprehensive shape analysis
	 */
	struct ShapeMetrics
	{
		Coord3D centroid;
		typename AdvancedAlgorithms::PCAResult pca;
		std::vector<Coord3D> convex_hull;
		double surface_area;
		size_t volume;
		double compactness;
		double sphericity;
		std::array<double, 3> principal_lengths;
	};

	static ShapeMetrics analyze_shape(const SparseContainer& container,
																		std::function<bool(const T&)> condition)
	{

		ShapeMetrics metrics;

		metrics.centroid =
			AdvancedAlgorithms::compute_centroid(container, condition);
		metrics.pca = AdvancedAlgorithms::compute_pca(container, condition);
		metrics.convex_hull = AdvancedAlgorithms::convex_hull(container, condition);
		metrics.surface_area =
			AdvancedAlgorithms::compute_surface_area(container, condition);
		metrics.volume = AdvancedAlgorithms::compute_volume(container, condition);

		if (metrics.volume > 0) {
			metrics.compactness =
				std::pow(metrics.surface_area, 1.5) / metrics.volume;
			double pi = 3.14159265358979323846;
			metrics.sphericity =
				(std::pow(pi, 1.0 / 3.0) * std::pow(6 * metrics.volume, 2.0 / 3.0)) /
				metrics.surface_area;
		} else {
			metrics.compactness = 0.0;
			metrics.sphericity = 0.0;
		}

		for (size_t i = 0; i < 3; ++i) {
			metrics.principal_lengths[i] = std::sqrt(metrics.pca.eigenvalues[i]);
		}

		return metrics;
	}

	// ======================
	// UTILITY FUNCTIONS
	// ======================

	/**
	 * Benchmark algorithm performance
	 */
	struct BenchmarkResult
	{
		std::string algorithm_name;
		double execution_time_ms;
		size_t memory_usage_bytes;
		bool success;
		std::string error_message;
	};

	template<typename Func, typename... Args>
	static BenchmarkResult benchmark_algorithm(const std::string& name,
																						 Func&& func,
																						 Args&&... args)
	{

		BenchmarkResult result;
		result.algorithm_name = name;
		result.error_message = "";

		auto start_time = std::chrono::high_resolution_clock::now();

		try {
			auto algorithm_result = func(std::forward<Args>(args)...);
			result.success = true;

			if constexpr (std::is_same_v<decltype(algorithm_result),
																	 std::shared_ptr<SparseContainer>>) {
				result.memory_usage_bytes = algorithm_result->get_memory_usage_bytes();
			} else {
				result.memory_usage_bytes = 0;
			}
		} catch (const std::exception& e) {
			result.success = false;
			result.memory_usage_bytes = 0;
			result.error_message = e.what();
		} catch (...) {
			result.success = false;
			result.memory_usage_bytes = 0;
			result.error_message = "Unknown error occurred";
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			end_time - start_time);
		result.execution_time_ms = duration.count() / 1000.0;

		return result;
	}

	/**
	 * Algorithm recommendation based on data characteristics
	 */
	struct DataCharacteristics
	{
		size_t total_voxels;
		double sparsity_ratio;
		size_t num_components;
		double avg_component_size;
		BoundingBox3D bounding_box;
		double memory_usage_mb;
	};

	static DataCharacteristics analyze_data_characteristics(
		const SparseContainer& container,
		std::function<bool(const T&)> condition)
	{

		DataCharacteristics chars;
		chars.total_voxels = container.size();
		chars.bounding_box = container.get_bounding_box();
		chars.memory_usage_mb =
			container.get_memory_usage_bytes() / (1024.0 * 1024.0);

		size_t valid_voxels = 0;
		for (const auto& [coord, voxel] : container) {
			if (condition(voxel.value)) {
				valid_voxels++;
			}
		}

		if (chars.bounding_box.is_valid) {
			int64_t total_possible = chars.bounding_box.volume();
			chars.sparsity_ratio = static_cast<double>(valid_voxels) / total_possible;
		} else {
			chars.sparsity_ratio = 0.0;
		}

		auto components =
			BasicAlgorithms::find_connected_components(container, condition);
		chars.num_components = components.size();

		if (!components.empty()) {
			size_t total_component_size = 0;
			for (const auto& component : components) {
				total_component_size += component.size();
			}
			chars.avg_component_size =
				static_cast<double>(total_component_size) / components.size();
		} else {
			chars.avg_component_size = 0.0;
		}

		return chars;
	}

	/**
	 * Get algorithm recommendations
	 */
	static std::vector<std::string> recommend_algorithms(
		const DataCharacteristics& chars)
	{

		std::vector<std::string> recommendations;

		if (chars.sparsity_ratio < 0.1) {
			recommendations.push_back(
				"Use morphological operations for noise removal");
			recommendations.push_back("Consider distance transforms for analysis");
		} else if (chars.sparsity_ratio > 0.8) {
			recommendations.push_back(
				"Data is dense - consider converting to Grid3D");
			recommendations.push_back("Use efficient filtering algorithms");
		}

		if (chars.num_components > 100) {
			recommendations.push_back(
				"Many components detected - use component filtering");
			recommendations.push_back("Consider multi-scale analysis");
		} else if (chars.num_components == 1) {
			recommendations.push_back("Single component - focus on shape analysis");
			recommendations.push_back("Skeletonization may be useful");
		}

		if (chars.avg_component_size < 10) {
			recommendations.push_back("Small components - use noise filtering");
			recommendations.push_back("Apply morphological closing");
		} else if (chars.avg_component_size > 1000) {
			recommendations.push_back("Large components - consider subdivision");
			recommendations.push_back("Surface analysis recommended");
		}

		if (chars.memory_usage_mb > 500) {
			recommendations.push_back(
				"High memory usage - consider processing in chunks");
			recommendations.push_back("Use streaming algorithms for large datasets");
		}

		return recommendations;
	}

	/**
	 * Generate comprehensive report
	 */
	static std::string generate_analysis_report(
		const SparseContainer& container,
		std::function<bool(const T&)> condition)
	{

		std::ostringstream report;

		auto chars = analyze_data_characteristics(container, condition);
		auto recommendations = recommend_algorithms(chars);

		report << "=== SPARSE 3D DATA ANALYSIS REPORT ===\n\n";

		report << "Data Characteristics:\n";
		report << "  Total Voxels: " << chars.total_voxels << "\n";
		report << "  Sparsity Ratio: " << (chars.sparsity_ratio * 100.0) << "%\n";
		report << "  Connected Components: " << chars.num_components << "\n";
		report << "  Average Component Size: " << chars.avg_component_size << "\n";
		report << "  Memory Usage: " << chars.memory_usage_mb << " MB\n";

		if (chars.bounding_box.is_valid) {
			report << "  Bounding Box: " << chars.bounding_box.min_coord.to_string()
						 << " to " << chars.bounding_box.max_coord.to_string() << "\n";
			report << "  Bounding Volume: " << chars.bounding_box.volume() << "\n";
		}

		report << "\nRecommendations:\n";
		for (const auto& rec : recommendations) {
			report << "  - " << rec << "\n";
		}

		return report.str();
	}
};

// Type aliases for common use cases
using Sparse3DAlgorithmsFloat = Sparse3DAlgorithms<float>;
using Sparse3DAlgorithmsDouble = Sparse3DAlgorithms<double>;
using Sparse3DAlgorithmsInt = Sparse3DAlgorithms<int32_t>;

using Sparse3DAlgorithmsAdvancedFloat = Sparse3DAlgorithmsAdvanced<float>;
using Sparse3DAlgorithmsAdvancedDouble = Sparse3DAlgorithmsAdvanced<double>;
using Sparse3DAlgorithmsAdvancedInt = Sparse3DAlgorithmsAdvanced<int32_t>;

} // namespace dagger2
