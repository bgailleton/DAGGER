#pragma once

#include <numeric>
#include <unordered_map>
#include <vector>

/**
 * Unified Union-Find (Disjoint Set) data structure
 *
 * This class provides a unified interface that maintains compatibility with
 * multiple existing UnionFind implementations. It supports both method naming
 * conventions and provides optional component tracking.
 *
 * Features:
 * - Path compression for efficient find operations
 * - Union by rank for balanced tree structure
 * - Component counting (optional, enabled by default)
 * - Multiple method aliases for backward compatibility
 * - Const-correct find operation
 */
class UnionFind
{
private:
	mutable std::vector<size_t> parent_;
	std::vector<size_t> rank_;
	size_t num_components_;
	bool track_components_;

public:
	/**
	 * Constructor
	 * @param n Number of elements
	 * @param track_components Whether to track component count (default: true)
	 */
	explicit UnionFind(size_t n, bool track_components = true)
		: parent_(n)
		, rank_(n, 0)
		, num_components_(n)
		, track_components_(track_components)
	{
		std::iota(parent_.begin(), parent_.end(), 0);
	}

	/**
	 * Find the root of the set containing element x
	 * Uses path compression for efficiency
	 * @param x Element to find
	 * @return Root of the set containing x
	 */
	size_t find(size_t x) const
	{
		if (parent_[x] != x) {
			parent_[x] = find(parent_[x]); // Path compression
		}
		return parent_[x];
	}

	/**
	 * Union two sets containing elements x and y
	 * Uses union by rank for balanced trees
	 * @param x First element
	 * @param y Second element
	 * @return true if union was performed, false if already in same set
	 */
	bool unite(size_t x, size_t y)
	{
		size_t root_x = find(x);
		size_t root_y = find(y);

		if (root_x == root_y)
			return false;

		// Union by rank
		if (rank_[root_x] < rank_[root_y]) {
			parent_[root_x] = root_y;
		} else if (rank_[root_x] > rank_[root_y]) {
			parent_[root_y] = root_x;
		} else {
			parent_[root_y] = root_x;
			rank_[root_x]++;
		}

		if (track_components_) {
			num_components_--;
		}
		return true;
	}

	/**
	 * Alias for unite() - for backward compatibility with Kruskal's algorithm
	 * version
	 * @param x First element
	 * @param y Second element
	 * @return true if union was performed, false if already in same set
	 */
	bool union_sets(size_t x, size_t y) { return unite(x, y); }

	/**
	 * Check if two elements are in the same connected component
	 * @param x First element
	 * @param y Second element
	 * @return true if x and y are connected
	 */
	bool connected(size_t x, size_t y) const { return find(x) == find(y); }

	/**
	 * Get the number of connected components
	 * @return Number of components (only valid if track_components is true)
	 */
	size_t num_components() const
	{
		return track_components_ ? num_components_ : count_components();
	}

	/**
	 * Get the size of the original set
	 * @return Total number of elements
	 */
	size_t size() const { return parent_.size(); }

	/**
	 * Get all connected components as vectors
	 * @return Vector of components, where each component is a vector of element
	 * indices
	 */
	std::vector<std::vector<size_t>> get_components() const
	{
		std::unordered_map<size_t, std::vector<size_t>> component_map;

		for (size_t i = 0; i < parent_.size(); ++i) {
			component_map[find(i)].push_back(i);
		}

		std::vector<std::vector<size_t>> components;
		components.reserve(component_map.size());
		for (const auto& [root, members] : component_map) {
			components.push_back(members);
		}

		return components;
	}

	/**
	 * Get the size of the component containing element x
	 * @param x Element to query
	 * @return Size of component containing x
	 */
	size_t component_size(size_t x) const
	{
		size_t root = find(x);
		size_t count = 0;
		for (size_t i = 0; i < parent_.size(); ++i) {
			if (find(i) == root) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Reset the data structure to initial state (all elements in separate sets)
	 */
	void reset()
	{
		std::iota(parent_.begin(), parent_.end(), 0);
		std::fill(rank_.begin(), rank_.end(), 0);
		if (track_components_) {
			num_components_ = parent_.size();
		}
	}

	/**
	 * Check if the data structure is in a valid state (for debugging)
	 * @return true if all invariants hold
	 */
	bool is_valid() const
	{
		// Check that all parents are valid indices
		for (size_t i = 0; i < parent_.size(); ++i) {
			if (parent_[i] >= parent_.size()) {
				return false;
			}
		}

		// Check component count if tracking is enabled
		if (track_components_) {
			size_t actual_components = count_components();
			if (actual_components != num_components_) {
				return false;
			}
		}

		return true;
	}

private:
	/**
	 * Count components by traversing all elements (expensive)
	 * Used when component tracking is disabled
	 */
	size_t count_components() const
	{
		std::unordered_map<size_t, bool> roots;
		for (size_t i = 0; i < parent_.size(); ++i) {
			roots[find(i)] = true;
		}
		return roots.size();
	}
};

// ======================
// CONVENIENCE ALIASES
// ======================

// For backward compatibility, you can use these typedefs if needed:
using UnionFindWithComponents = UnionFind; // Default behavior
using UnionFindBasic =
	UnionFind; // Can be constructed with track_components=false
