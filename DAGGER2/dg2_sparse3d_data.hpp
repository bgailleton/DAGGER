#pragma once

#include "dg2_sparse3d_core.hpp"
#include <algorithm>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * Voxel data container - stores actual data values with type flexibility
 */
template<typename T>
struct VoxelData
{
	T value;
	uint32_t flags; // For future use (material types, states, etc.)

	VoxelData()
		: value(T{})
		, flags(0)
	{
	}
	VoxelData(const T& val, uint32_t f = 0)
		: value(val)
		, flags(f)
	{
	}

	bool operator==(const VoxelData& other) const
	{
		return value == other.value && flags == other.flags;
	}

	bool operator!=(const VoxelData& other) const { return !(*this == other); }
};

/**
 * Sparse3DContainer: Core sparse data structure for 3D voxel data
 *
 * This class provides efficient storage and retrieval of sparse 3D data
 * using hash maps for O(1) access while maintaining spatial coherence.
 */
template<typename T>
class Sparse3DContainer
{
public:
	using VoxelType = VoxelData<T>;
	using DataMap = std::unordered_map<Coord3D, VoxelType>;
	using Iterator = typename DataMap::iterator;
	using ConstIterator = typename DataMap::const_iterator;

private:
	DataMap data_;
	BoundingBox3D bounding_box_;
	T default_value_;
	uint32_t default_flags_;
	mutable bool bbox_dirty_;

	// Spatial indexing for faster queries
	mutable std::map<int64_t, std::set<Coord3D>> x_index_;
	mutable std::map<int64_t, std::set<Coord3D>> y_index_;
	mutable std::map<int64_t, std::set<Coord3D>> z_index_;
	mutable bool indices_dirty_;

public:
	/**
	 * Constructor
	 */
	Sparse3DContainer(const T& default_val = T{}, uint32_t default_flags = 0)
		: default_value_(default_val)
		, default_flags_(default_flags)
		, bbox_dirty_(true)
		, indices_dirty_(true)
	{
	}

	/**
	 * Copy constructor
	 */
	Sparse3DContainer(const Sparse3DContainer& other)
		: data_(other.data_)
		, bounding_box_(other.bounding_box_)
		, default_value_(other.default_value_)
		, default_flags_(other.default_flags_)
		, bbox_dirty_(other.bbox_dirty_)
		, indices_dirty_(true)
	{
	}

	/**
	 * Move constructor
	 */
	Sparse3DContainer(Sparse3DContainer&& other) noexcept
		: data_(std::move(other.data_))
		, bounding_box_(other.bounding_box_)
		, default_value_(other.default_value_)
		, default_flags_(other.default_flags_)
		, bbox_dirty_(other.bbox_dirty_)
		, indices_dirty_(true)
	{
	}

	/**
	 * Assignment operators
	 */
	Sparse3DContainer& operator=(const Sparse3DContainer& other)
	{
		if (this != &other) {
			data_ = other.data_;
			bounding_box_ = other.bounding_box_;
			default_value_ = other.default_value_;
			default_flags_ = other.default_flags_;
			bbox_dirty_ = other.bbox_dirty_;
			indices_dirty_ = true;
		}
		return *this;
	}

	Sparse3DContainer& operator=(Sparse3DContainer&& other) noexcept
	{
		if (this != &other) {
			data_ = std::move(other.data_);
			bounding_box_ = other.bounding_box_;
			default_value_ = other.default_value_;
			default_flags_ = other.default_flags_;
			bbox_dirty_ = other.bbox_dirty_;
			indices_dirty_ = true;
		}
		return *this;
	}

	// ======================
	// BASIC ACCESS
	// ======================

	/**
	 * Get voxel data (returns default if not present)
	 */
	const VoxelType& get(const Coord3D& coord) const
	{
		auto it = data_.find(coord);
		if (it != data_.end()) {
			return it->second;
		}

		static VoxelType default_voxel;
		default_voxel.value = default_value_;
		default_voxel.flags = default_flags_;
		return default_voxel;
	}

	/**
	 * Get voxel value directly
	 */
	T get_value(const Coord3D& coord) const
	{
		auto it = data_.find(coord);
		return (it != data_.end()) ? it->second.value : default_value_;
	}

	/**
	 * Check if voxel exists (is not default)
	 */
	bool exists(const Coord3D& coord) const
	{
		return data_.find(coord) != data_.end();
	}

	/**
	 * Set voxel data
	 */
	void set(const Coord3D& coord, const VoxelType& voxel)
	{
		data_[coord] = voxel;
		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Set voxel value (uses default flags)
	 */
	void set_value(const Coord3D& coord, const T& value)
	{
		data_[coord] = VoxelType(value, default_flags_);
		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Remove voxel (returns to default state)
	 */
	bool remove(const Coord3D& coord)
	{
		auto it = data_.find(coord);
		if (it != data_.end()) {
			data_.erase(it);
			bbox_dirty_ = true;
			indices_dirty_ = true;
			return true;
		}
		return false;
	}

	/**
	 * Clear all data
	 */
	void clear()
	{
		data_.clear();
		bounding_box_ = BoundingBox3D();
		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	// ======================
	// PROPERTIES
	// ======================

	size_t size() const { return data_.size(); }
	bool empty() const { return data_.empty(); }

	T get_default_value() const { return default_value_; }
	uint32_t get_default_flags() const { return default_flags_; }

	void set_default_value(const T& value) { default_value_ = value; }
	void set_default_flags(uint32_t flags) { default_flags_ = flags; }

	// ======================
	// ITERATORS
	// ======================

	Iterator begin() { return data_.begin(); }
	Iterator end() { return data_.end(); }
	ConstIterator begin() const { return data_.cbegin(); }
	ConstIterator end() const { return data_.cend(); }
	ConstIterator cbegin() const { return data_.cbegin(); }
	ConstIterator cend() const { return data_.cend(); }

	// ======================
	// SPATIAL QUERIES
	// ======================

	/**
	 * Get bounding box of all stored voxels
	 */
	const BoundingBox3D& get_bounding_box() const
	{
		if (bbox_dirty_) {
			update_bounding_box();
		}
		return bounding_box_;
	}

	/**
	 * Get all coordinates within a bounding box
	 */
	std::vector<Coord3D> get_coords_in_bbox(const BoundingBox3D& bbox) const
	{
		std::vector<Coord3D> result;
		result.reserve(std::min(data_.size(), static_cast<size_t>(bbox.volume())));

		for (const auto& [coord, voxel] : data_) {
			if (bbox.contains(coord)) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Get all coordinates within a sphere
	 */
	std::vector<Coord3D> get_coords_in_sphere(const Coord3D& center,
																						double radius) const
	{
		std::vector<Coord3D> result;
		double radius_sq = radius * radius;

		for (const auto& [coord, voxel] : data_) {
			if (coord.distance_to(center) <= radius) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Get all coordinates at specific X plane
	 */
	std::vector<Coord3D> get_coords_at_x(int64_t x) const
	{
		update_indices_if_needed();

		std::vector<Coord3D> result;
		auto it = x_index_.find(x);
		if (it != x_index_.end()) {
			result.reserve(it->second.size());
			for (const auto& coord : it->second) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Get all coordinates at specific Y plane
	 */
	std::vector<Coord3D> get_coords_at_y(int64_t y) const
	{
		update_indices_if_needed();

		std::vector<Coord3D> result;
		auto it = y_index_.find(y);
		if (it != y_index_.end()) {
			result.reserve(it->second.size());
			for (const auto& coord : it->second) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Get all coordinates at specific Z plane
	 */
	std::vector<Coord3D> get_coords_at_z(int64_t z) const
	{
		update_indices_if_needed();

		std::vector<Coord3D> result;
		auto it = z_index_.find(z);
		if (it != z_index_.end()) {
			result.reserve(it->second.size());
			for (const auto& coord : it->second) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Get all unique X coordinates
	 */
	std::vector<int64_t> get_unique_x_coords() const
	{
		update_indices_if_needed();

		std::vector<int64_t> result;
		result.reserve(x_index_.size());
		for (const auto& [x, coords] : x_index_) {
			result.push_back(x);
		}

		return result;
	}

	/**
	 * Get all unique Y coordinates
	 */
	std::vector<int64_t> get_unique_y_coords() const
	{
		update_indices_if_needed();

		std::vector<int64_t> result;
		result.reserve(y_index_.size());
		for (const auto& [y, coords] : y_index_) {
			result.push_back(y);
		}

		return result;
	}

	/**
	 * Get all unique Z coordinates
	 */
	std::vector<int64_t> get_unique_z_coords() const
	{
		update_indices_if_needed();

		std::vector<int64_t> result;
		result.reserve(z_index_.size());
		for (const auto& [z, coords] : z_index_) {
			result.push_back(z);
		}

		return result;
	}

	// ======================
	// NEIGHBOR OPERATIONS
	// ======================

	/**
	 * Get neighbor coordinates for a given connectivity
	 */
	std::vector<Neighbor3D> get_neighbors(const Coord3D& coord,
																				Connectivity3D connectivity) const
	{
		std::vector<Neighbor3D> neighbors;
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);
		neighbors.reserve(directions.size());

		for (Direction3D dir : directions) {
			Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
			Coord3D neighbor_coord = coord + offset;
			double distance = Direction3DUtils::get_distance(dir);
			bool neighbor_exists = exists(neighbor_coord);

			neighbors.emplace_back(neighbor_coord, dir, distance, neighbor_exists);
		}

		return neighbors;
	}

	/**
	 * Get only existing neighbors
	 */
	std::vector<Neighbor3D> get_existing_neighbors(
		const Coord3D& coord,
		Connectivity3D connectivity) const
	{
		std::vector<Neighbor3D> neighbors;
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		for (Direction3D dir : directions) {
			Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
			Coord3D neighbor_coord = coord + offset;

			if (exists(neighbor_coord)) {
				double distance = Direction3DUtils::get_distance(dir);
				neighbors.emplace_back(neighbor_coord, dir, distance, true);
			}
		}

		return neighbors;
	}

	/**
	 * Count existing neighbors
	 */
	size_t count_neighbors(const Coord3D& coord,
												 Connectivity3D connectivity) const
	{
		size_t count = 0;
		auto directions =
			Direction3DUtils::get_directions_for_connectivity(connectivity);

		for (Direction3D dir : directions) {
			Coord3D offset = Direction3DUtils::get_neighbor_offset(dir);
			Coord3D neighbor_coord = coord + offset;
			if (exists(neighbor_coord)) {
				++count;
			}
		}

		return count;
	}

	// ======================
	// BULK OPERATIONS
	// ======================

	/**
	 * Set multiple voxels at once
	 */
	void set_bulk(const std::vector<std::pair<Coord3D, VoxelType>>& voxels)
	{
		for (const auto& [coord, voxel] : voxels) {
			data_[coord] = voxel;
		}
		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Set multiple voxel values at once
	 */
	void set_values_bulk(const std::vector<std::pair<Coord3D, T>>& values)
	{
		for (const auto& [coord, value] : values) {
			data_[coord] = VoxelType(value, default_flags_);
		}
		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Remove multiple voxels at once
	 */
	size_t remove_bulk(const std::vector<Coord3D>& coords)
	{
		size_t removed_count = 0;
		for (const auto& coord : coords) {
			auto it = data_.find(coord);
			if (it != data_.end()) {
				data_.erase(it);
				++removed_count;
			}
		}

		if (removed_count > 0) {
			bbox_dirty_ = true;
			indices_dirty_ = true;
		}

		return removed_count;
	}

	/**
	 * Get multiple voxel values at once
	 */
	std::vector<T> get_values_bulk(const std::vector<Coord3D>& coords) const
	{
		std::vector<T> values;
		values.reserve(coords.size());

		for (const auto& coord : coords) {
			values.push_back(get_value(coord));
		}

		return values;
	}

	/**
	 * Fill a 3D region with a value
	 */
	void fill_region(const BoundingBox3D& bbox,
									 const T& value,
									 uint32_t flags = 0)
	{
		if (!bbox.is_valid)
			return;

		VoxelType voxel(value, flags);

		for (int64_t x = bbox.min_coord.x; x <= bbox.max_coord.x; ++x) {
			for (int64_t y = bbox.min_coord.y; y <= bbox.max_coord.y; ++y) {
				for (int64_t z = bbox.min_coord.z; z <= bbox.max_coord.z; ++z) {
					data_[Coord3D(x, y, z)] = voxel;
				}
			}
		}

		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Fill a spherical region with a value
	 */
	void fill_sphere(const Coord3D& center,
									 double radius,
									 const T& value,
									 uint32_t flags = 0)
	{
		VoxelType voxel(value, flags);
		double radius_sq = radius * radius;

		// Determine bounding box for the sphere
		int64_t r_int = static_cast<int64_t>(std::ceil(radius));
		BoundingBox3D sphere_bbox(
			Coord3D(center.x - r_int, center.y - r_int, center.z - r_int),
			Coord3D(center.x + r_int, center.y + r_int, center.z + r_int));

		for (int64_t x = sphere_bbox.min_coord.x; x <= sphere_bbox.max_coord.x;
				 ++x) {
			for (int64_t y = sphere_bbox.min_coord.y; y <= sphere_bbox.max_coord.y;
					 ++y) {
				for (int64_t z = sphere_bbox.min_coord.z; z <= sphere_bbox.max_coord.z;
						 ++z) {
					Coord3D coord(x, y, z);
					if (coord.distance_to(center) <= radius) {
						data_[coord] = voxel;
					}
				}
			}
		}

		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Remove all voxels in a region
	 */
	size_t clear_region(const BoundingBox3D& bbox)
	{
		if (!bbox.is_valid)
			return 0;

		size_t removed_count = 0;
		std::vector<Coord3D> to_remove;

		for (const auto& [coord, voxel] : data_) {
			if (bbox.contains(coord)) {
				to_remove.push_back(coord);
			}
		}

		for (const auto& coord : to_remove) {
			data_.erase(coord);
			++removed_count;
		}

		if (removed_count > 0) {
			bbox_dirty_ = true;
			indices_dirty_ = true;
		}

		return removed_count;
	}

	// ======================
	// TRANSFORMATION OPERATIONS
	// ======================

	/**
	 * Transform all coordinates using a function
	 */
	void transform_coordinates(
		std::function<Coord3D(const Coord3D&)> transform_func)
	{
		DataMap new_data;
		new_data.reserve(data_.size());

		for (const auto& [coord, voxel] : data_) {
			Coord3D new_coord = transform_func(coord);
			new_data[new_coord] = voxel;
		}

		data_ = std::move(new_data);
		bbox_dirty_ = true;
		indices_dirty_ = true;
	}

	/**
	 * Translate all coordinates by an offset
	 */
	void translate(const Coord3D& offset)
	{
		transform_coordinates(
			[offset](const Coord3D& coord) { return coord + offset; });
	}

	/**
	 * Scale all coordinates by factors
	 */
	void scale(int64_t scale_x, int64_t scale_y, int64_t scale_z)
	{
		transform_coordinates([scale_x, scale_y, scale_z](const Coord3D& coord) {
			return Coord3D(coord.x * scale_x, coord.y * scale_y, coord.z * scale_z);
		});
	}

	/**
	 * Mirror coordinates along specified axes
	 */
	void mirror(bool mirror_x, bool mirror_y, bool mirror_z)
	{
		transform_coordinates([mirror_x, mirror_y, mirror_z](const Coord3D& coord) {
			return Coord3D(mirror_x ? -coord.x : coord.x,
										 mirror_y ? -coord.y : coord.y,
										 mirror_z ? -coord.z : coord.z);
		});
	}

	// ======================
	// STATISTICAL OPERATIONS
	// ======================

	/**
	 * Apply function to all values
	 */
	void transform_values(std::function<T(const T&)> transform_func)
	{
		for (auto& [coord, voxel] : data_) {
			voxel.value = transform_func(voxel.value);
		}
	}

	/**
	 * Get statistics of all stored values
	 */
	struct Statistics
	{
		T min_value;
		T max_value;
		double mean_value;
		double variance;
		size_t count;

		Statistics()
			: min_value(T{})
			, max_value(T{})
			, mean_value(0.0)
			, variance(0.0)
			, count(0)
		{
		}
	};

	Statistics compute_statistics() const
	{
		Statistics stats;
		if (data_.empty())
			return stats;

		stats.count = data_.size();
		auto first_value = data_.begin()->second.value;
		stats.min_value = first_value;
		stats.max_value = first_value;

		double sum = 0.0;
		for (const auto& [coord, voxel] : data_) {
			T value = voxel.value;
			stats.min_value = std::min(stats.min_value, value);
			stats.max_value = std::max(stats.max_value, value);
			sum += static_cast<double>(value);
		}

		stats.mean_value = sum / stats.count;

		// Compute variance
		double variance_sum = 0.0;
		for (const auto& [coord, voxel] : data_) {
			double diff = static_cast<double>(voxel.value) - stats.mean_value;
			variance_sum += diff * diff;
		}
		stats.variance = variance_sum / stats.count;

		return stats;
	}

	/**
	 * Find coordinates with specific value
	 */
	std::vector<Coord3D> find_coords_with_value(const T& target_value) const
	{
		std::vector<Coord3D> result;

		for (const auto& [coord, voxel] : data_) {
			if (voxel.value == target_value) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Find coordinates with values in range
	 */
	std::vector<Coord3D> find_coords_in_range(const T& min_value,
																						const T& max_value) const
	{
		std::vector<Coord3D> result;

		for (const auto& [coord, voxel] : data_) {
			if (voxel.value >= min_value && voxel.value <= max_value) {
				result.push_back(coord);
			}
		}

		return result;
	}

	/**
	 * Count voxels with specific value
	 */
	size_t count_value(const T& target_value) const
	{
		size_t count = 0;
		for (const auto& [coord, voxel] : data_) {
			if (voxel.value == target_value) {
				++count;
			}
		}
		return count;
	}

	// ======================
	// DEBUGGING AND UTILITIES
	// ======================

	/**
	 * Get memory usage estimation
	 */
	size_t get_memory_usage_bytes() const
	{
		size_t base_size = sizeof(*this);
		size_t data_size = data_.size() * (sizeof(Coord3D) + sizeof(VoxelType));
		size_t index_size = 0;

		if (!indices_dirty_) {
			index_size += x_index_.size() * sizeof(int64_t);
			index_size += y_index_.size() * sizeof(int64_t);
			index_size += z_index_.size() * sizeof(int64_t);

			for (const auto& [x, coords] : x_index_) {
				index_size += coords.size() * sizeof(Coord3D);
			}
			for (const auto& [y, coords] : y_index_) {
				index_size += coords.size() * sizeof(Coord3D);
			}
			for (const auto& [z, coords] : z_index_) {
				index_size += coords.size() * sizeof(Coord3D);
			}
		}

		return base_size + data_size + index_size;
	}

	/**
	 * Get summary information
	 */
	std::string get_summary() const
	{
		std::ostringstream ss;
		ss << "Sparse3DContainer Summary:\n";
		ss << "  Size: " << size() << " voxels\n";
		ss << "  Memory: " << (get_memory_usage_bytes() / 1024.0 / 1024.0)
			 << " MB\n";
		ss << "  Default value: " << default_value_ << "\n";

		if (!empty()) {
			auto bbox = get_bounding_box();
			ss << "  Bounding box: " << bbox.min_coord.to_string() << " to "
				 << bbox.max_coord.to_string() << "\n";
			ss << "  Volume: " << bbox.volume() << " voxels\n";
			ss << "  Density: " << (100.0 * size() / bbox.volume()) << "%\n";

			auto stats = compute_statistics();
			ss << "  Value range: [" << stats.min_value << ", " << stats.max_value
				 << "]\n";
			ss << "  Mean: " << stats.mean_value << "\n";
			ss << "  Variance: " << stats.variance << "\n";
		}

		return ss.str();
	}

	/**
	 * Validate internal consistency
	 */
	bool validate() const
	{
		// Check if all stored coordinates are within computed bounding box
		if (!empty()) {
			auto bbox = get_bounding_box();
			for (const auto& [coord, voxel] : data_) {
				if (!bbox.contains(coord)) {
					return false;
				}
			}
		}

		return true;
	}

private:
	// ======================
	// PRIVATE METHODS
	// ======================

	void update_bounding_box() const
	{
		bounding_box_ = BoundingBox3D();

		for (const auto& [coord, voxel] : data_) {
			bounding_box_.expand_to_include(coord);
		}

		bbox_dirty_ = false;
	}

	void update_indices_if_needed() const
	{
		if (!indices_dirty_)
			return;

		x_index_.clear();
		y_index_.clear();
		z_index_.clear();

		for (const auto& [coord, voxel] : data_) {
			x_index_[coord.x].insert(coord);
			y_index_[coord.y].insert(coord);
			z_index_[coord.z].insert(coord);
		}

		indices_dirty_ = false;
	}
};

// Type aliases for common use cases
using Sparse3DFloat = Sparse3DContainer<float>;
using Sparse3DDouble = Sparse3DContainer<double>;
using Sparse3DInt = Sparse3DContainer<int32_t>;
using Sparse3DLong = Sparse3DContainer<int64_t>;
using Sparse3DBool = Sparse3DContainer<bool>;

} // namespace dagger2
