#pragma once

#include "dg2_array.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dagger2 {

/**
 * 3D coordinate structure with efficient hashing and comparison
 */
struct Coord3D
{
	int64_t x, y, z;

	Coord3D()
		: x(0)
		, y(0)
		, z(0)
	{
	}
	Coord3D(int64_t x_, int64_t y_, int64_t z_)
		: x(x_)
		, y(y_)
		, z(z_)
	{
	}

	bool operator==(const Coord3D& other) const
	{
		return x == other.x && y == other.y && z == other.z;
	}

	bool operator!=(const Coord3D& other) const { return !(*this == other); }

	bool operator<(const Coord3D& other) const
	{
		if (x != other.x)
			return x < other.x;
		if (y != other.y)
			return y < other.y;
		return z < other.z;
	}

	Coord3D operator+(const Coord3D& other) const
	{
		return Coord3D(x + other.x, y + other.y, z + other.z);
	}

	Coord3D operator-(const Coord3D& other) const
	{
		return Coord3D(x - other.x, y - other.y, z - other.z);
	}

	Coord3D operator*(int64_t scalar) const
	{
		return Coord3D(x * scalar, y * scalar, z * scalar);
	}

	double distance_to(const Coord3D& other) const
	{
		double dx = static_cast<double>(x - other.x);
		double dy = static_cast<double>(y - other.y);
		double dz = static_cast<double>(z - other.z);
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	int64_t manhattan_distance_to(const Coord3D& other) const
	{
		return std::abs(x - other.x) + std::abs(y - other.y) +
					 std::abs(z - other.z);
	}

	std::string to_string() const
	{
		return "(" + std::to_string(x) + "," + std::to_string(y) + "," +
					 std::to_string(z) + ")";
	}
};

/**
 * Bounding box for 3D coordinates
 */
struct BoundingBox3D
{
	Coord3D min_coord;
	Coord3D max_coord;
	bool is_valid;

	BoundingBox3D()
		: is_valid(false)
	{
	}

	BoundingBox3D(const Coord3D& min_c, const Coord3D& max_c)
		: min_coord(min_c)
		, max_coord(max_c)
		, is_valid(true)
	{
	}

	void expand_to_include(const Coord3D& coord)
	{
		if (!is_valid) {
			min_coord = max_coord = coord;
			is_valid = true;
		} else {
			min_coord.x = std::min(min_coord.x, coord.x);
			min_coord.y = std::min(min_coord.y, coord.y);
			min_coord.z = std::min(min_coord.z, coord.z);
			max_coord.x = std::max(max_coord.x, coord.x);
			max_coord.y = std::max(max_coord.y, coord.y);
			max_coord.z = std::max(max_coord.z, coord.z);
		}
	}

	bool contains(const Coord3D& coord) const
	{
		if (!is_valid)
			return false;
		return coord.x >= min_coord.x && coord.x <= max_coord.x &&
					 coord.y >= min_coord.y && coord.y <= max_coord.y &&
					 coord.z >= min_coord.z && coord.z <= max_coord.z;
	}

	bool intersects(const BoundingBox3D& other) const
	{
		if (!is_valid || !other.is_valid)
			return false;
		return !(
			max_coord.x < other.min_coord.x || min_coord.x > other.max_coord.x ||
			max_coord.y < other.min_coord.y || min_coord.y > other.max_coord.y ||
			max_coord.z < other.min_coord.z || min_coord.z > other.max_coord.z);
	}

	int64_t volume() const
	{
		if (!is_valid)
			return 0;
		return (max_coord.x - min_coord.x + 1) * (max_coord.y - min_coord.y + 1) *
					 (max_coord.z - min_coord.z + 1);
	}

	std::vector<Coord3D> get_corners() const
	{
		if (!is_valid)
			return {};

		return { min_coord,
						 { max_coord.x, min_coord.y, min_coord.z },
						 { min_coord.x, max_coord.y, min_coord.z },
						 { min_coord.x, min_coord.y, max_coord.z },
						 { max_coord.x, max_coord.y, min_coord.z },
						 { max_coord.x, min_coord.y, max_coord.z },
						 { min_coord.x, max_coord.y, max_coord.z },
						 max_coord };
	}
};

} // namespace dagger2

// Hash function for Coord3D
namespace std {
template<>
struct hash<dagger2::Coord3D>
{
	size_t operator()(const dagger2::Coord3D& coord) const
	{
		// Use a good hash combining function
		size_t h1 = std::hash<int64_t>{}(coord.x);
		size_t h2 = std::hash<int64_t>{}(coord.y);
		size_t h3 = std::hash<int64_t>{}(coord.z);

		// Combine hashes using boost-style hash_combine
		size_t seed = h1;
		seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};
}

namespace dagger2 {

/**
 * 3D direction enumeration for neighbor connectivity
 */
enum class Direction3D : uint8_t
{
	// 6-connected (face neighbors)
	EAST = 0,	 // +X
	WEST = 1,	 // -X
	NORTH = 2, // +Y
	SOUTH = 3, // -Y
	UP = 4,		 // +Z
	DOWN = 5,	 // -Z

	// 18-connected (face + edge neighbors)
	NORTHEAST = 6,	 // +X+Y
	NORTHWEST = 7,	 // -X+Y
	SOUTHEAST = 8,	 // +X-Y
	SOUTHWEST = 9,	 // -X-Y
	EAST_UP = 10,		 // +X+Z
	EAST_DOWN = 11,	 // +X-Z
	WEST_UP = 12,		 // -X+Z
	WEST_DOWN = 13,	 // -X-Z
	NORTH_UP = 14,	 // +Y+Z
	NORTH_DOWN = 15, // +Y-Z
	SOUTH_UP = 16,	 // +Y+Z
	SOUTH_DOWN = 17, // +Y-Z

	// 26-connected (face + edge + corner neighbors)
	NORTHEAST_UP = 18,	 // +X+Y+Z
	NORTHEAST_DOWN = 19, // +X+Y-Z
	NORTHWEST_UP = 20,	 // -X+Y+Z
	NORTHWEST_DOWN = 21, // -X+Y-Z
	SOUTHEAST_UP = 22,	 // +X-Y+Z
	SOUTHEAST_DOWN = 23, // +X-Y-Z
	SOUTHWEST_UP = 24,	 // -X-Y+Z
	SOUTHWEST_DOWN = 25, // -X-Y-Z

	CENTER = 26,	// Current voxel (self)
	INVALID = 255 // Invalid direction
};

/**
 * Connectivity type for 3D grids
 */
enum class Connectivity3D : uint8_t
{
	C6 = 6,		// 6-connected (face neighbors only)
	C18 = 18, // 18-connected (face + edge neighbors)
	C26 = 26	// 26-connected (face + edge + corner neighbors)
};

/**
 * Direction utilities for 3D operations
 */
class Direction3DUtils
{
private:
	// Direction offset vectors [dx, dy, dz]
	static constexpr std::array<std::array<int, 3>, 26> direction_offsets_ = { {
		// 6-connected (faces)
		{ 1, 0, 0 },	// EAST
		{ -1, 0, 0 }, // WEST
		{ 0, 1, 0 },	// NORTH
		{ 0, -1, 0 }, // SOUTH
		{ 0, 0, 1 },	// UP
		{ 0, 0, -1 }, // DOWN

		// 12 additional edge neighbors
		{ 1, 1, 0 },	 // NORTHEAST
		{ -1, 1, 0 },	 // NORTHWEST
		{ 1, -1, 0 },	 // SOUTHEAST
		{ -1, -1, 0 }, // SOUTHWEST
		{ 1, 0, 1 },	 // EAST_UP
		{ 1, 0, -1 },	 // EAST_DOWN
		{ -1, 0, 1 },	 // WEST_UP
		{ -1, 0, -1 }, // WEST_DOWN
		{ 0, 1, 1 },	 // NORTH_UP
		{ 0, 1, -1 },	 // NORTH_DOWN
		{ 0, -1, 1 },	 // SOUTH_UP
		{ 0, -1, -1 }, // SOUTH_DOWN

		// 8 additional corner neighbors
		{ 1, 1, 1 },	 // NORTHEAST_UP
		{ 1, 1, -1 },	 // NORTHEAST_DOWN
		{ -1, 1, 1 },	 // NORTHWEST_UP
		{ -1, 1, -1 }, // NORTHWEST_DOWN
		{ 1, -1, 1 },	 // SOUTHEAST_UP
		{ 1, -1, -1 }, // SOUTHEAST_DOWN
		{ -1, -1, 1 }, // SOUTHWEST_UP
		{ -1, -1, -1 } // SOUTHWEST_DOWN
	} };

	// Distances for each direction
	static constexpr std::array<double, 26> distances_ = {
		// 6-connected (distance = 1)
		1.0,
		1.0,
		1.0,
		1.0,
		1.0,
		1.0,

		// 12 edge neighbors (distance = sqrt(2))
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,
		1.4142135623730951,

		// 8 corner neighbors (distance = sqrt(3))
		1.7320508075688772,
		1.7320508075688772,
		1.7320508075688772,
		1.7320508075688772,
		1.7320508075688772,
		1.7320508075688772,
		1.7320508075688772,
		1.7320508075688772
	};

public:
	static Coord3D get_neighbor_offset(Direction3D dir)
	{
		uint8_t idx = static_cast<uint8_t>(dir);
		if (idx >= 26)
			return Coord3D(0, 0, 0);

		const auto& offset = direction_offsets_[idx];
		return Coord3D(offset[0], offset[1], offset[2]);
	}

	static double get_distance(Direction3D dir)
	{
		uint8_t idx = static_cast<uint8_t>(dir);
		if (idx >= 26)
			return 0.0;
		return distances_[idx];
	}

	static std::vector<Direction3D> get_directions_for_connectivity(
		Connectivity3D connectivity)
	{
		std::vector<Direction3D> directions;
		uint8_t max_dir = static_cast<uint8_t>(connectivity);

		directions.reserve(max_dir);
		for (uint8_t i = 0; i < max_dir; ++i) {
			directions.push_back(static_cast<Direction3D>(i));
		}

		return directions;
	}

	static Direction3D get_opposite_direction(Direction3D dir)
	{
		switch (dir) {
			case Direction3D::EAST:
				return Direction3D::WEST;
			case Direction3D::WEST:
				return Direction3D::EAST;
			case Direction3D::NORTH:
				return Direction3D::SOUTH;
			case Direction3D::SOUTH:
				return Direction3D::NORTH;
			case Direction3D::UP:
				return Direction3D::DOWN;
			case Direction3D::DOWN:
				return Direction3D::UP;

			case Direction3D::NORTHEAST:
				return Direction3D::SOUTHWEST;
			case Direction3D::NORTHWEST:
				return Direction3D::SOUTHEAST;
			case Direction3D::SOUTHEAST:
				return Direction3D::NORTHWEST;
			case Direction3D::SOUTHWEST:
				return Direction3D::NORTHEAST;
			case Direction3D::EAST_UP:
				return Direction3D::WEST_DOWN;
			case Direction3D::EAST_DOWN:
				return Direction3D::WEST_UP;
			case Direction3D::WEST_UP:
				return Direction3D::EAST_DOWN;
			case Direction3D::WEST_DOWN:
				return Direction3D::EAST_UP;
			case Direction3D::NORTH_UP:
				return Direction3D::SOUTH_DOWN;
			case Direction3D::NORTH_DOWN:
				return Direction3D::SOUTH_UP;
			case Direction3D::SOUTH_UP:
				return Direction3D::NORTH_DOWN;
			case Direction3D::SOUTH_DOWN:
				return Direction3D::NORTH_UP;

			case Direction3D::NORTHEAST_UP:
				return Direction3D::SOUTHWEST_DOWN;
			case Direction3D::NORTHEAST_DOWN:
				return Direction3D::SOUTHWEST_UP;
			case Direction3D::NORTHWEST_UP:
				return Direction3D::SOUTHEAST_DOWN;
			case Direction3D::NORTHWEST_DOWN:
				return Direction3D::SOUTHEAST_UP;
			case Direction3D::SOUTHEAST_UP:
				return Direction3D::NORTHWEST_DOWN;
			case Direction3D::SOUTHEAST_DOWN:
				return Direction3D::NORTHWEST_UP;
			case Direction3D::SOUTHWEST_UP:
				return Direction3D::NORTHEAST_DOWN;
			case Direction3D::SOUTHWEST_DOWN:
				return Direction3D::NORTHEAST_UP;

			default:
				return Direction3D::INVALID;
		}
	}

	static std::string to_string(Direction3D dir)
	{
		switch (dir) {
			case Direction3D::EAST:
				return "EAST";
			case Direction3D::WEST:
				return "WEST";
			case Direction3D::NORTH:
				return "NORTH";
			case Direction3D::SOUTH:
				return "SOUTH";
			case Direction3D::UP:
				return "UP";
			case Direction3D::DOWN:
				return "DOWN";
			case Direction3D::NORTHEAST:
				return "NORTHEAST";
			case Direction3D::NORTHWEST:
				return "NORTHWEST";
			case Direction3D::SOUTHEAST:
				return "SOUTHEAST";
			case Direction3D::SOUTHWEST:
				return "SOUTHWEST";
			case Direction3D::EAST_UP:
				return "EAST_UP";
			case Direction3D::EAST_DOWN:
				return "EAST_DOWN";
			case Direction3D::WEST_UP:
				return "WEST_UP";
			case Direction3D::WEST_DOWN:
				return "WEST_DOWN";
			case Direction3D::NORTH_UP:
				return "NORTH_UP";
			case Direction3D::NORTH_DOWN:
				return "NORTH_DOWN";
			case Direction3D::SOUTH_UP:
				return "SOUTH_UP";
			case Direction3D::SOUTH_DOWN:
				return "SOUTH_DOWN";
			case Direction3D::NORTHEAST_UP:
				return "NORTHEAST_UP";
			case Direction3D::NORTHEAST_DOWN:
				return "NORTHEAST_DOWN";
			case Direction3D::NORTHWEST_UP:
				return "NORTHWEST_UP";
			case Direction3D::NORTHWEST_DOWN:
				return "NORTHWEST_DOWN";
			case Direction3D::SOUTHEAST_UP:
				return "SOUTHEAST_UP";
			case Direction3D::SOUTHEAST_DOWN:
				return "SOUTHEAST_DOWN";
			case Direction3D::SOUTHWEST_UP:
				return "SOUTHWEST_UP";
			case Direction3D::SOUTHWEST_DOWN:
				return "SOUTHWEST_DOWN";
			case Direction3D::CENTER:
				return "CENTER";
			case Direction3D::INVALID:
				return "INVALID";
			default:
				return "UNKNOWN";
		}
	}
};

/**
 * Neighbor information for 3D sparse grids
 */
struct Neighbor3D
{
	Coord3D coord;
	Direction3D direction;
	double distance;
	bool exists;

	Neighbor3D()
		: direction(Direction3D::INVALID)
		, distance(0.0)
		, exists(false)
	{
	}

	Neighbor3D(const Coord3D& c, Direction3D dir, double dist, bool ex = true)
		: coord(c)
		, direction(dir)
		, distance(dist)
		, exists(ex)
	{
	}
};

} // namespace dagger2
