#pragma once

#include "dg2_array.hpp" // For Grid2D

#ifdef _MSC_VER
#include <intrin.h>
#include <malloc.h> // For _aligned_malloc
#else
#include <cstdlib> // For aligned_alloc
#endif

#if defined(__x86_64__) || defined(_M_X64)
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <immintrin.h>
#endif
#define HAS_AVX2
#endif
#include <array>
#include <cassert>
#include <cstring>
#include <memory>
#include <vector>

// Forward declarations from your existing codebase
enum class ConnectivityType : uint8_t
{
	D4 = 4,
	D8 = 8
};
enum class NodeType : uint8_t
{
	NO_DATA = 0,
	NORMAL = 1,
	HAS_TO_OUT = 2,
	CAN_OUT = 3,
	IN = 4,
	PERIODIC = 5,
	REFLECT = 6
};
enum class Direction : uint8_t
{
	NORTH = 0,
	EAST = 1,
	SOUTH = 2,
	WEST = 3,
	NORTHEAST = 4,
	SOUTHEAST = 5,
	SOUTHWEST = 6,
	NORTHWEST = 7,
	INVALID = 255
};

namespace dagger2 {

/**
 * Ultra-optimized FastConnector for high-performance PDE solving
 * Template specializations for D4/D8 connectivity with SIMD optimizations
 */
template<typename T = double, ConnectivityType CONN = ConnectivityType::D8>
class FastConnector
{
private:
	static constexpr size_t NUM_DIRS = static_cast<size_t>(CONN);
	static constexpr size_t CACHE_LINE = 64;
	static constexpr size_t SIMD_WIDTH = sizeof(__m256) / sizeof(T);

	// SoA Memory Layout - cache-aligned
	struct alignas(CACHE_LINE) ConnectorData
	{
		// Grid properties
		size_t rows, cols, size;

		// Aligned coordinate arrays
		T* __restrict__ coord_x;
		T* __restrict__ coord_y;

		// Boundary condition bitfield (4 bits per node)
		uint32_t* __restrict__ boundary_bits;

		// Pre-computed neighbor indices (lazy-allocated)
		mutable size_t* __restrict__ neighbor_indices;
		mutable bool* __restrict__ neighbor_computed;

		// Memory pool for temporary operations
		T* __restrict__ temp_buffer;
		size_t temp_capacity;
	} data_;

	// Lookup Tables (compile-time constants)
	static constexpr std::array<int, 8> DR = { -1, 0, 1, 0, -1, 1, 1, -1 };
	static constexpr std::array<int, 8> DC = { 0, 1, 0, -1, 1, 1, -1, -1 };
	static constexpr std::array<T, 8> DISTANCES = { 1.0,
																									1.0,
																									1.0,
																									1.0,
																									1.4142135623730951,
																									1.4142135623730951,
																									1.4142135623730951,
																									1.4142135623730951 };
	static constexpr std::array<Direction, 8> OPPOSITES = {
		Direction::SOUTH,			Direction::WEST,			Direction::NORTH,
		Direction::EAST,			Direction::SOUTHWEST, Direction::NORTHWEST,
		Direction::NORTHEAST, Direction::SOUTHEAST
	};

	// Cross-platform bit operations
	inline static size_t count_trailing_zeros(size_t x) noexcept
	{
#ifdef _MSC_VER
		unsigned long index;
		return _BitScanForward64(&index, x) ? index : 64;
#else
		return __builtin_ctzl(x);
#endif
	}

	inline static size_t popcount(size_t x) noexcept
	{
#ifdef _MSC_VER
		return __popcnt64(x);
#else
		return __builtin_popcountl(x);
#endif
	}

	// Cross-platform aligned allocation
	static void* aligned_alloc_compat(size_t alignment, size_t size)
	{
#ifdef _MSC_VER
		return _aligned_malloc(size, alignment);
#else
		return std::aligned_alloc(alignment, size);
#endif
	}

	static void aligned_free_compat(void* ptr)
	{
#ifdef _MSC_VER
		_aligned_free(ptr);
#else
		std::free(ptr);
#endif
	}

	// Memory management
	void allocate_aligned_memory()
	{
		size_t total_coords = data_.size * 2;
		size_t total_boundary = (data_.size + 7) / 8; // 4 bits per node, packed
		size_t total_neighbors = data_.size * NUM_DIRS;

		data_.coord_x = static_cast<T*>(
			aligned_alloc_compat(CACHE_LINE, total_coords * sizeof(T)));
		data_.coord_y = data_.coord_x + data_.size;

		data_.boundary_bits = static_cast<uint32_t*>(
			aligned_alloc_compat(CACHE_LINE, total_boundary * sizeof(uint32_t)));

		// Lazy allocation for neighbors
		data_.neighbor_indices = static_cast<size_t*>(
			aligned_alloc_compat(CACHE_LINE, total_neighbors * sizeof(size_t)));
		data_.neighbor_computed = static_cast<bool*>(
			aligned_alloc_compat(CACHE_LINE, data_.size * sizeof(bool)));

		data_.temp_capacity = data_.size * 4; // 4x size for safety
		data_.temp_buffer = static_cast<T*>(
			aligned_alloc_compat(CACHE_LINE, data_.temp_capacity * sizeof(T)));

		// Initialize coordinates
		for (size_t i = 0; i < data_.size; ++i) {
			auto [row, col] = to_2d_fast(i);
			data_.coord_x[i] = static_cast<T>(col);
			data_.coord_y[i] = static_cast<T>(row);
		}

		// Initialize boundary bits (all NORMAL by default)
		std::memset(data_.boundary_bits, 0x11, total_boundary * sizeof(uint32_t));
		std::memset(data_.neighbor_computed, false, data_.size);
	}

	void deallocate_memory()
	{
		aligned_free_compat(data_.coord_x);
		aligned_free_compat(data_.boundary_bits);
		aligned_free_compat(data_.neighbor_indices);
		aligned_free_compat(data_.neighbor_computed);
		aligned_free_compat(data_.temp_buffer);
	}

public:
	// Neighbor structure for return values - matching original interface
	struct Neighbor
	{
		size_t index;
		T distance;
		Direction direction;
		bool is_valid;
		NodeType boundary_type; // Added for PriorityFlood compatibility

		Neighbor()
			: index(SIZE_MAX)
			, distance(0)
			, direction(Direction::INVALID)
			, is_valid(false)
			, boundary_type(NodeType::NO_DATA)
		{
		}
		Neighbor(size_t idx, T dist, Direction dir, NodeType bc = NodeType::NORMAL)
			: index(idx)
			, distance(dist)
			, direction(dir)
			, is_valid(true)
			, boundary_type(bc)
		{
		}
	};

	// FlowTransfer struct for compatibility
	struct FlowTransfer
	{
		size_t destination_index;
		T flux_multiplier;
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

	/**
	 * Constructor
	 */
	FastConnector(size_t rows, size_t cols)
		: data_{ rows,		cols,		 rows * cols, nullptr, nullptr,
						 nullptr, nullptr, nullptr,			nullptr, 0 }
	{
		if (rows == 0 || cols == 0) {
			throw std::invalid_argument("Grid dimensions must be positive");
		}
		allocate_aligned_memory();
	}

	/**
	 * Destructor
	 */
	~FastConnector() { deallocate_memory(); }

	// Copy/move constructors and operators
	FastConnector(const FastConnector&) = delete;
	FastConnector& operator=(const FastConnector&) = delete;
	FastConnector(FastConnector&&) = delete;
	FastConnector& operator=(FastConnector&&) = delete;

	// Basic properties
	inline size_t rows() const noexcept { return data_.rows; }
	inline size_t cols() const noexcept { return data_.cols; }
	inline size_t size() const noexcept { return data_.size; }
	inline ConnectivityType connectivity_type() const noexcept { return CONN; }
	inline size_t num_directions() const noexcept { return NUM_DIRS; }

	// Ultra-fast coordinate conversion with bit operations
	inline size_t to_1d_fast(size_t row, size_t col) const noexcept
	{
		return (row << count_trailing_zeros(data_.cols)) +
					 col; // Assumes cols is power of 2
	}

	inline size_t to_1d_safe(size_t row, size_t col) const noexcept
	{
		return row * data_.cols + col;
	}

	inline std::pair<size_t, size_t> to_2d_fast(size_t index) const noexcept
	{
		// Fast division using bit operations (assumes cols is power of 2)
		size_t col_mask = data_.cols - 1;
		return { index >> count_trailing_zeros(data_.cols), index & col_mask };
	}

	inline std::pair<size_t, size_t> to_2d_safe(size_t index) const noexcept
	{
		return { index / data_.cols, index % data_.cols };
	}

	// Choose optimal conversion based on cols being power of 2
	inline size_t to_1d(size_t row, size_t col) const noexcept
	{
		return (popcount(data_.cols) == 1) ? to_1d_fast(row, col)
																			 : to_1d_safe(row, col);
	}

	inline std::pair<size_t, size_t> to_2d(size_t index) const noexcept
	{
		return (popcount(data_.cols) == 1) ? to_2d_fast(index) : to_2d_safe(index);
	}

	// SIMD-optimized coordinate validation
	inline bool is_valid_coord(int row, int col) const noexcept
	{
		// Branchless validation using bit operations
		uint32_t row_valid = static_cast<uint32_t>(row) < data_.rows;
		uint32_t col_valid = static_cast<uint32_t>(col) < data_.cols;
		return (row_valid & col_valid) != 0;
	}

	inline bool is_valid_coord(size_t row, size_t col) const noexcept
	{
		return (row < data_.rows) & (col < data_.cols);
	}

	inline bool is_valid_index(size_t index) const noexcept
	{
		return index < data_.size;
	}

	// Vectorized batch coordinate validation
	void validate_coordinates_batch(const size_t* rows,
																	const size_t* cols,
																	bool* results,
																	size_t count) const noexcept
	{
#ifdef __AVX2__
		const __m256i max_row = _mm256_set1_epi64x(data_.rows);
		const __m256i max_col = _mm256_set1_epi64x(data_.cols);

		size_t simd_count = count & ~3; // Process 4 at a time
		for (size_t i = 0; i < simd_count; i += 4) {
			__m256i vrows =
				_mm256_loadu_si256(reinterpret_cast<const __m256i*>(&rows[i]));
			__m256i vcols =
				_mm256_loadu_si256(reinterpret_cast<const __m256i*>(&cols[i]));

			__m256i row_valid = _mm256_cmpgt_epi64(max_row, vrows);
			__m256i col_valid = _mm256_cmpgt_epi64(max_col, vcols);
			__m256i valid = _mm256_and_si256(row_valid, col_valid);

			int mask = _mm256_movemask_pd(_mm256_castsi256_pd(valid));
			results[i] = (mask & 1) != 0;
			results[i + 1] = (mask & 2) != 0;
			results[i + 2] = (mask & 4) != 0;
			results[i + 3] = (mask & 8) != 0;
		}

		// Handle remainder
		for (size_t i = simd_count; i < count; ++i) {
			results[i] = is_valid_coord(rows[i], cols[i]);
		}
#else
		for (size_t i = 0; i < count; ++i) {
			results[i] = is_valid_coord(rows[i], cols[i]);
		}
#endif
	}

	// Bit-manipulated boundary conditions
	inline NodeType get_boundary_type(size_t index) const noexcept
	{
		size_t byte_idx = index >> 3;					// Divide by 8
		size_t bit_offset = (index & 7) << 2; // (index % 8) * 4
		uint32_t bits = (data_.boundary_bits[byte_idx] >> bit_offset) & 0xF;
		return static_cast<NodeType>(bits);
	}

	inline NodeType get_boundary_type(size_t row, size_t col) const noexcept
	{
		return get_boundary_type(to_1d(row, col));
	}

	inline void set_boundary_type(size_t index, NodeType type) noexcept
	{
		size_t byte_idx = index >> 3;
		size_t bit_offset = (index & 7) << 2;
		uint32_t mask = ~(0xFU << bit_offset);
		data_.boundary_bits[byte_idx] = (data_.boundary_bits[byte_idx] & mask) |
																		(static_cast<uint32_t>(type) << bit_offset);
	}

	inline void set_boundary_type(size_t row, size_t col, NodeType type) noexcept
	{
		set_boundary_type(to_1d(row, col), type);
	}

	// Fast boundary queries using bit operations
	inline bool is_active_node(size_t index) const noexcept
	{
		NodeType type = get_boundary_type(index);
		return static_cast<uint8_t>(type) != 0; // NO_DATA = 0
	}

	inline bool is_active_node(size_t row, size_t col) const noexcept
	{
		return is_active_node(to_1d(row, col));
	}

	inline bool is_boundary_node(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return (row == 0) | (row == data_.rows - 1) | (col == 0) |
					 (col == data_.cols - 1);
	}

	// Template specialization for connectivity-specific neighbor calculation
	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D4, void>::type
	compute_neighbors_unrolled(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		size_t base_idx = index * NUM_DIRS;

		// Unrolled D4 neighbor computation
		// North
		if (row > 0) {
			data_.neighbor_indices[base_idx] = to_1d(row - 1, col);
		} else {
			data_.neighbor_indices[base_idx] = SIZE_MAX;
		}

		// East
		if (col < data_.cols - 1) {
			data_.neighbor_indices[base_idx + 1] = to_1d(row, col + 1);
		} else {
			data_.neighbor_indices[base_idx + 1] = SIZE_MAX;
		}

		// South
		if (row < data_.rows - 1) {
			data_.neighbor_indices[base_idx + 2] = to_1d(row + 1, col);
		} else {
			data_.neighbor_indices[base_idx + 2] = SIZE_MAX;
		}

		// West
		if (col > 0) {
			data_.neighbor_indices[base_idx + 3] = to_1d(row, col - 1);
		} else {
			data_.neighbor_indices[base_idx + 3] = SIZE_MAX;
		}

		data_.neighbor_computed[index] = true;
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, void>::type
	compute_neighbors_unrolled(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		size_t base_idx = index * NUM_DIRS;

		// Unrolled D8 neighbor computation using lookup tables
		for (size_t d = 0; d < NUM_DIRS; ++d) {
			int new_row = static_cast<int>(row) + DR[d];
			int new_col = static_cast<int>(col) + DC[d];

			if (is_valid_coord(new_row, new_col)) {
				data_.neighbor_indices[base_idx + d] =
					to_1d(static_cast<size_t>(new_row), static_cast<size_t>(new_col));
			} else {
				data_.neighbor_indices[base_idx + d] = SIZE_MAX;
			}
		}

		data_.neighbor_computed[index] = true;
	}

	// Lazy neighbor computation
	inline void ensure_neighbors_computed(size_t index) const noexcept
	{
		if (!data_.neighbor_computed[index]) {
			compute_neighbors_unrolled(index);
		}
	}

	// Vectorized direction calculation
	std::vector<Neighbor> get_neighbors(size_t index) const
	{
		ensure_neighbors_computed(index);

		std::vector<Neighbor> neighbors;
		neighbors.reserve(NUM_DIRS);

		size_t base_idx = index * NUM_DIRS;
		for (size_t d = 0; d < NUM_DIRS; ++d) {
			size_t neighbor_idx = data_.neighbor_indices[base_idx + d];
			if (neighbor_idx != SIZE_MAX && is_active_node(neighbor_idx)) {
				neighbors.emplace_back(
					neighbor_idx, DISTANCES[d], static_cast<Direction>(d));
			}
		}

		return neighbors;
	}

	// Batch neighbor computation for multiple indices
	void get_neighbors_batch(const size_t* indices,
													 size_t count,
													 std::vector<std::vector<Neighbor>>& results) const
	{
		results.resize(count);

		for (size_t i = 0; i < count; ++i) {
			results[i] = get_neighbors(indices[i]);
		}
	}

	// High-performance coordinate transformation for PDEs
	void transform_coordinates_batch(const size_t* input_indices,
																	 size_t count,
																	 T* output_x,
																	 T* output_y) const noexcept
	{
#ifdef __AVX2__
		constexpr size_t simd_width = 4; // 4 doubles per AVX2 register
		size_t simd_count = count & ~(simd_width - 1);

		for (size_t i = 0; i < simd_count; i += simd_width) {
			// Load indices
			__m256i indices =
				_mm256_loadu_si256(reinterpret_cast<const __m256i*>(&input_indices[i]));

			// Gather coordinates
			__m256d x_vals = _mm256_i64gather_pd(data_.coord_x, indices, 8);
			__m256d y_vals = _mm256_i64gather_pd(data_.coord_y, indices, 8);

			// Store results
			_mm256_storeu_pd(&output_x[i], x_vals);
			_mm256_storeu_pd(&output_y[i], y_vals);
		}

		// Handle remainder
		for (size_t i = simd_count; i < count; ++i) {
			size_t idx = input_indices[i];
			output_x[i] = data_.coord_x[idx];
			output_y[i] = data_.coord_y[idx];
		}
#else
		for (size_t i = 0; i < count; ++i) {
			size_t idx = input_indices[i];
			output_x[i] = data_.coord_x[idx];
			output_y[i] = data_.coord_y[idx];
		}
#endif
	}

	// Memory pool access for temporary computations
	T* get_temp_buffer(size_t required_size) const noexcept
	{
		return (required_size <= data_.temp_capacity) ? data_.temp_buffer : nullptr;
	}

	// ===================================================================
	// RAW NEIGHBOR ACCESS - High Performance Index-Only Methods
	// ===================================================================

	// Raw neighbor indices (1D) - valid only, boundary-aware
	void get_neighbors_raw(size_t index,
												 size_t* output,
												 size_t& count) const noexcept
	{
		ensure_neighbors_computed(index);
		count = 0;

		size_t base_idx = index * NUM_DIRS;
		for (size_t d = 0; d < NUM_DIRS; ++d) {
			size_t neighbor_idx = data_.neighbor_indices[base_idx + d];
			if (neighbor_idx != SIZE_MAX && is_active_node(neighbor_idx)) {
				output[count++] = neighbor_idx;
			}
		}
	}

	void get_neighbors_raw(size_t row,
												 size_t col,
												 size_t* output,
												 size_t& count) const noexcept
	{
		get_neighbors_raw(to_1d(row, col), output, count);
	}

	// Raw neighbor coordinates (2D) - valid only, boundary-aware
	void get_neighbors_raw_2d(size_t index,
														size_t* rows_out,
														size_t* cols_out,
														size_t& count) const noexcept
	{
		size_t temp_indices[NUM_DIRS];
		get_neighbors_raw(index, temp_indices, count);

		for (size_t i = 0; i < count; ++i) {
			auto [r, c] = to_2d(temp_indices[i]);
			rows_out[i] = r;
			cols_out[i] = c;
		}
	}

	void get_neighbors_raw_2d(size_t row,
														size_t col,
														size_t* rows_out,
														size_t* cols_out,
														size_t& count) const noexcept
	{
		get_neighbors_raw_2d(to_1d(row, col), rows_out, cols_out, count);
	}

	// Raw effective neighbors (1D) - handles boundary conditions
	void get_effective_neighbors_raw(size_t index,
																	 size_t* output,
																	 size_t& count) const noexcept
	{
		auto [row, col] = to_2d(index);
		get_effective_neighbors_raw(row, col, output, count);
	}

	void get_effective_neighbors_raw(size_t row,
																	 size_t col,
																	 size_t* output,
																	 size_t& count) const noexcept
	{
		count = 0;

		for (size_t d = 0; d < NUM_DIRS; ++d) {
			int target_row = static_cast<int>(row) + DR[d];
			int target_col = static_cast<int>(col) + DC[d];

			if (is_valid_coord(target_row, target_col)) {
				size_t target_idx = to_1d(static_cast<size_t>(target_row),
																	static_cast<size_t>(target_col));
				if (is_active_node(target_idx)) {
					output[count++] = target_idx;
				}
			} else {
				// Handle boundary conditions
				NodeType from_bc = get_boundary_type(row, col);

				if (from_bc == NodeType::PERIODIC) {
					auto [wrapped_row, wrapped_col] =
						apply_periodic_bc(target_row, target_col);
					size_t wrapped_idx = to_1d(wrapped_row, wrapped_col);
					if (is_active_node(wrapped_idx)) {
						output[count++] = wrapped_idx;
					}
				} else if (from_bc == NodeType::REFLECT) {
					size_t from_idx = to_1d(row, col);
					output[count++] = from_idx; // Reflect back to self
				}
				// HAS_TO_OUT, CAN_OUT: no neighbor added (flow exits)
			}
		}
	}

	// Raw effective neighbors (2D)
	void get_effective_neighbors_raw_2d(size_t index,
																			size_t* rows_out,
																			size_t* cols_out,
																			size_t& count) const noexcept
	{
		size_t temp_indices[NUM_DIRS];
		get_effective_neighbors_raw(index, temp_indices, count);

		for (size_t i = 0; i < count; ++i) {
			auto [r, c] = to_2d(temp_indices[i]);
			rows_out[i] = r;
			cols_out[i] = c;
		}
	}

	void get_effective_neighbors_raw_2d(size_t row,
																			size_t col,
																			size_t* rows_out,
																			size_t* cols_out,
																			size_t& count) const noexcept
	{
		get_effective_neighbors_raw_2d(to_1d(row, col), rows_out, cols_out, count);
	}

	// Raw effective valid neighbors (1D) - excludes NO_DATA
	void get_effective_valid_neighbors_raw(size_t index,
																				 size_t* output,
																				 size_t& count) const noexcept
	{
		auto [row, col] = to_2d(index);
		get_effective_valid_neighbors_raw(row, col, output, count);
	}

	void get_effective_valid_neighbors_raw(size_t row,
																				 size_t col,
																				 size_t* output,
																				 size_t& count) const noexcept
	{
		count = 0;

		for (size_t d = 0; d < NUM_DIRS; ++d) {
			int target_row = static_cast<int>(row) + DR[d];
			int target_col = static_cast<int>(col) + DC[d];

			if (is_valid_coord(target_row, target_col)) {
				size_t target_idx = to_1d(static_cast<size_t>(target_row),
																	static_cast<size_t>(target_col));
				NodeType target_bc = get_boundary_type(target_idx);
				if (target_bc != NodeType::NO_DATA) {
					output[count++] = target_idx;
				}
			} else {
				// Handle boundary conditions
				NodeType from_bc = get_boundary_type(row, col);

				if (from_bc == NodeType::PERIODIC) {
					auto [wrapped_row, wrapped_col] =
						apply_periodic_bc(target_row, target_col);
					size_t wrapped_idx = to_1d(wrapped_row, wrapped_col);
					NodeType wrapped_bc = get_boundary_type(wrapped_idx);
					if (wrapped_bc != NodeType::NO_DATA) {
						output[count++] = wrapped_idx;
					}
				} else if (from_bc == NodeType::REFLECT) {
					size_t from_idx = to_1d(row, col);
					output[count++] = from_idx; // Reflect back to self
				}
				// HAS_TO_OUT, CAN_OUT: no neighbor added (flow exits)
			}
		}
	}

	// Raw effective valid neighbors (2D)
	void get_effective_valid_neighbors_raw_2d(size_t index,
																						size_t* rows_out,
																						size_t* cols_out,
																						size_t& count) const noexcept
	{
		size_t temp_indices[NUM_DIRS];
		get_effective_valid_neighbors_raw(index, temp_indices, count);

		for (size_t i = 0; i < count; ++i) {
			auto [r, c] = to_2d(temp_indices[i]);
			rows_out[i] = r;
			cols_out[i] = c;
		}
	}

	void get_effective_valid_neighbors_raw_2d(size_t row,
																						size_t col,
																						size_t* rows_out,
																						size_t* cols_out,
																						size_t& count) const noexcept
	{
		get_effective_valid_neighbors_raw_2d(
			to_1d(row, col), rows_out, cols_out, count);
	}

	// Batch raw neighbor processing (for vectorized algorithms)
	void get_neighbors_raw_batch(const size_t* indices,
															 size_t batch_size,
															 size_t* output,
															 size_t* counts,
															 size_t* offsets) const noexcept
	{
		size_t current_offset = 0;

		for (size_t i = 0; i < batch_size; ++i) {
			offsets[i] = current_offset;
			get_neighbors_raw(indices[i], &output[current_offset], counts[i]);
			current_offset += counts[i];
		}
	}

	// Ultra-fast single direction neighbor lookup
	inline size_t get_neighbor_raw_direction(size_t index,
																					 Direction dir) const noexcept
	{
		auto [row, col] = to_2d(index);
		return get_neighbor_raw_direction(row, col, dir);
	}

	inline size_t get_neighbor_raw_direction(size_t row,
																					 size_t col,
																					 Direction dir) const noexcept
	{
		uint8_t d = static_cast<uint8_t>(dir);
		int target_row = static_cast<int>(row) + DR[d];
		int target_col = static_cast<int>(col) + DC[d];

		if (is_valid_coord(target_row, target_col)) {
			size_t target_idx =
				to_1d(static_cast<size_t>(target_row), static_cast<size_t>(target_col));
			return is_active_node(target_idx) ? target_idx : SIZE_MAX;
		}

		return SIZE_MAX;
	}

	// Cross-platform vectorized coordinate validation (remove duplicate)
	// This method is already defined above, removing duplicate

	// Efficient distance calculation using lookup tables
	inline T get_distance(Direction dir) const noexcept
	{
		return DISTANCES[static_cast<size_t>(dir)];
	}

	inline Direction get_opposite_direction(Direction dir) const noexcept
	{
		return OPPOSITES[static_cast<size_t>(dir)];
	}

	// Performance statistics
	struct PerfStats
	{
		size_t cache_hits = 0;
		size_t cache_misses = 0;
		size_t simd_operations = 0;
		double cache_hit_ratio() const
		{
			return static_cast<double>(cache_hits) / (cache_hits + cache_misses);
		}
	};

	mutable PerfStats perf_stats_;
	const PerfStats& get_performance_stats() const noexcept
	{
		return perf_stats_;
	}
	void reset_performance_stats() const noexcept { perf_stats_ = {}; }

	// ===================================================================
	// MISSING INTERFACE METHODS - Required for FlowGraph and PriorityFlood
	// ===================================================================

	// Constructor with boundary grid (for compatibility)
	FastConnector(size_t rows,
								size_t cols,
								std::shared_ptr<Grid2D<NodeType>> boundary_grid)
		: data_{ rows,		cols,		 rows * cols, nullptr, nullptr,
						 nullptr, nullptr, nullptr,			nullptr, 0 }
	{
		if (rows == 0 || cols == 0) {
			throw std::invalid_argument("Grid dimensions must be positive");
		}
		allocate_aligned_memory();

		// Copy boundary grid data
		if (boundary_grid && boundary_grid->rows() == rows &&
				boundary_grid->cols() == cols) {
			for (size_t i = 0; i < data_.size; ++i) {
				set_boundary_type(i, boundary_grid->data()[i]);
			}
		}
	}

	// Direction-specific neighbor access (D8 connectivity)
	Neighbor north(size_t row, size_t col) const noexcept
	{
		if (row == 0)
			return Neighbor();
		size_t idx = to_1d(row - 1, col);
		return Neighbor(
			idx, DISTANCES[0], Direction::NORTH, get_boundary_type(idx));
	}

	Neighbor north(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return north(row, col);
	}

	Neighbor east(size_t row, size_t col) const noexcept
	{
		if (col >= data_.cols - 1)
			return Neighbor();
		size_t idx = to_1d(row, col + 1);
		return Neighbor(idx, DISTANCES[1], Direction::EAST, get_boundary_type(idx));
	}

	Neighbor east(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return east(row, col);
	}

	Neighbor south(size_t row, size_t col) const noexcept
	{
		if (row >= data_.rows - 1)
			return Neighbor();
		size_t idx = to_1d(row + 1, col);
		return Neighbor(
			idx, DISTANCES[2], Direction::SOUTH, get_boundary_type(idx));
	}

	Neighbor south(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return south(row, col);
	}

	Neighbor west(size_t row, size_t col) const noexcept
	{
		if (col == 0)
			return Neighbor();
		size_t idx = to_1d(row, col - 1);
		return Neighbor(idx, DISTANCES[3], Direction::WEST, get_boundary_type(idx));
	}

	Neighbor west(size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return west(row, col);
	}

	// Diagonal neighbors (D8 only)
	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type northeast(
		size_t row,
		size_t col) const noexcept
	{
		if (row == 0 || col >= data_.cols - 1)
			return Neighbor();
		size_t idx = to_1d(row - 1, col + 1);
		return Neighbor(
			idx, DISTANCES[4], Direction::NORTHEAST, get_boundary_type(idx));
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type northeast(
		size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return northeast(row, col);
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type southeast(
		size_t row,
		size_t col) const noexcept
	{
		if (row >= data_.rows - 1 || col >= data_.cols - 1)
			return Neighbor();
		size_t idx = to_1d(row + 1, col + 1);
		return Neighbor(
			idx, DISTANCES[5], Direction::SOUTHEAST, get_boundary_type(idx));
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type southeast(
		size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return southeast(row, col);
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type southwest(
		size_t row,
		size_t col) const noexcept
	{
		if (row >= data_.rows - 1 || col == 0)
			return Neighbor();
		size_t idx = to_1d(row + 1, col - 1);
		return Neighbor(
			idx, DISTANCES[6], Direction::SOUTHWEST, get_boundary_type(idx));
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type southwest(
		size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return southwest(row, col);
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type northwest(
		size_t row,
		size_t col) const noexcept
	{
		if (row == 0 || col == 0)
			return Neighbor();
		size_t idx = to_1d(row - 1, col - 1);
		return Neighbor(
			idx, DISTANCES[7], Direction::NORTHWEST, get_boundary_type(idx));
	}

	template<ConnectivityType C = CONN>
	typename std::enable_if<C == ConnectivityType::D8, Neighbor>::type northwest(
		size_t index) const noexcept
	{
		auto [row, col] = to_2d(index);
		return northwest(row, col);
	}

	// All neighbors access (includes invalid ones)
	std::vector<Neighbor> get_all_neighbors(size_t row, size_t col) const
	{
		std::vector<Neighbor> neighbors;
		neighbors.reserve(NUM_DIRS);

		if constexpr (NUM_DIRS == 4) {
			neighbors.push_back(north(row, col));
			neighbors.push_back(east(row, col));
			neighbors.push_back(south(row, col));
			neighbors.push_back(west(row, col));
		} else {
			neighbors.push_back(north(row, col));
			neighbors.push_back(east(row, col));
			neighbors.push_back(south(row, col));
			neighbors.push_back(west(row, col));
			neighbors.push_back(northeast(row, col));
			neighbors.push_back(southeast(row, col));
			neighbors.push_back(southwest(row, col));
			neighbors.push_back(northwest(row, col));
		}

		return neighbors;
	}

	std::vector<Neighbor> get_all_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_all_neighbors(row, col);
	}

	// Valid neighbors only (existing implementation already correct)
	std::vector<Neighbor> get_valid_neighbors(size_t row, size_t col) const
	{
		return get_neighbors(to_1d(row, col));
	}

	std::vector<Neighbor> get_valid_neighbors(size_t index) const
	{
		return get_neighbors(index);
	}

	// Vectorized neighbor indices (for performance)
	std::vector<size_t> get_neighbor_indices(size_t index) const
	{
		ensure_neighbors_computed(index);

		std::vector<size_t> indices;
		indices.reserve(NUM_DIRS);

		size_t base_idx = index * NUM_DIRS;
		for (size_t d = 0; d < NUM_DIRS; ++d) {
			indices.push_back(data_.neighbor_indices[base_idx + d]);
		}

		return indices;
	}

	std::vector<size_t> get_valid_neighbor_indices(size_t index) const
	{
		ensure_neighbors_computed(index);

		std::vector<size_t> indices;
		indices.reserve(NUM_DIRS);

		size_t base_idx = index * NUM_DIRS;
		for (size_t d = 0; d < NUM_DIRS; ++d) {
			size_t neighbor_idx = data_.neighbor_indices[base_idx + d];
			if (neighbor_idx != SIZE_MAX && is_active_node(neighbor_idx)) {
				indices.push_back(neighbor_idx);
			}
		}

		return indices;
	}

	// Boundary condition utilities
	std::pair<size_t, size_t> apply_periodic_bc(int target_row,
																							int target_col) const noexcept
	{
		size_t wrapped_row = (target_row < 0) ? data_.rows - 1
												 : (target_row >= static_cast<int>(data_.rows))
													 ? 0
													 : target_row;
		size_t wrapped_col = (target_col < 0) ? data_.cols - 1
												 : (target_col >= static_cast<int>(data_.cols))
													 ? 0
													 : target_col;
		return { wrapped_row, wrapped_col };
	}

	std::vector<size_t> get_boundary_nodes() const
	{
		std::vector<size_t> boundary_nodes;

		for (size_t i = 0; i < data_.size; ++i) {
			if (is_boundary_node(i)) {
				boundary_nodes.push_back(i);
			}
		}

		return boundary_nodes;
	}

	std::vector<size_t> get_all_boundary_nodes() const
	{
		return get_boundary_nodes(); // Same implementation for now
	}

	void set_border_boundary(NodeType type) noexcept
	{
		// Set all border nodes to specified type
		for (size_t r = 0; r < data_.rows; ++r) {
			for (size_t c = 0; c < data_.cols; ++c) {
				if (r == 0 || r == data_.rows - 1 || c == 0 || c == data_.cols - 1) {
					set_boundary_type(r, c, type);
				}
			}
		}
	}

	void set_periodic_boundaries() noexcept
	{
		set_border_boundary(NodeType::PERIODIC);
	}

	void set_reflective_boundaries() noexcept
	{
		set_border_boundary(NodeType::REFLECT);
	}

	// Effective neighbor access (with boundary condition handling)
	Neighbor get_effective_neighbor(size_t row,
																	size_t col,
																	Direction direction) const
	{
		uint8_t d = static_cast<uint8_t>(direction);
		int target_row = static_cast<int>(row) + DR[d];
		int target_col = static_cast<int>(col) + DC[d];

		if (is_valid_coord(target_row, target_col)) {
			size_t target_idx =
				to_1d(static_cast<size_t>(target_row), static_cast<size_t>(target_col));
			NodeType target_bc = get_boundary_type(target_idx);
			return Neighbor(target_idx, DISTANCES[d], direction, target_bc);
		}

		// Handle boundary conditions
		NodeType from_bc = get_boundary_type(row, col);

		switch (from_bc) {
			case NodeType::PERIODIC: {
				auto [wrapped_row, wrapped_col] =
					apply_periodic_bc(target_row, target_col);
				size_t wrapped_idx = to_1d(wrapped_row, wrapped_col);
				NodeType wrapped_bc = get_boundary_type(wrapped_idx);
				return Neighbor(wrapped_idx, DISTANCES[d], direction, wrapped_bc);
			}
			case NodeType::REFLECT: {
				size_t from_idx = to_1d(row, col);
				return Neighbor(from_idx, DISTANCES[d], direction, from_bc);
			}
			default:
				return Neighbor(); // Invalid
		}
	}

	std::vector<Neighbor> get_effective_neighbors(size_t row, size_t col) const
	{
		std::vector<Neighbor> neighbors;
		neighbors.reserve(NUM_DIRS);

		for (size_t d = 0; d < NUM_DIRS; ++d) {
			Neighbor n = get_effective_neighbor(row, col, static_cast<Direction>(d));
			if (n.is_valid) {
				neighbors.push_back(n);
			}
		}

		return neighbors;
	}

	std::vector<Neighbor> get_effective_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_effective_neighbors(row, col);
	}

	std::vector<Neighbor> get_effective_valid_neighbors(size_t row,
																											size_t col) const
	{
		std::vector<Neighbor> neighbors = get_effective_neighbors(row, col);

		// Filter out NO_DATA nodes
		neighbors.erase(std::remove_if(neighbors.begin(),
																	 neighbors.end(),
																	 [](const Neighbor& n) {
																		 return n.boundary_type ==
																						NodeType::NO_DATA;
																	 }),
										neighbors.end());

		return neighbors;
	}

	std::vector<Neighbor> get_effective_valid_neighbors(size_t index) const
	{
		auto [row, col] = to_2d(index);
		return get_effective_valid_neighbors(row, col);
	}

	bool has_boundary_handling() const noexcept
	{
		return true; // FastConnector always handles boundaries
	}

	// Flow transfer method (for advanced flow routing)
	FlowTransfer transfer_flow(size_t from_index,
														 Direction direction,
														 T flux) const
	{
		auto [row, col] = to_2d(from_index);
		NodeType from_bc = get_boundary_type(row, col);

		uint8_t d = static_cast<uint8_t>(direction);
		int target_row = static_cast<int>(row) + DR[d];
		int target_col = static_cast<int>(col) + DC[d];

		if (is_valid_coord(target_row, target_col)) {
			size_t target_idx =
				to_1d(static_cast<size_t>(target_row), static_cast<size_t>(target_col));
			NodeType target_bc = get_boundary_type(target_idx);

			if (target_bc != NodeType::NO_DATA) {
				return FlowTransfer(
					target_idx, static_cast<T>(1.0), false, false, false);
			} else {
				return FlowTransfer(); // Flow lost to NO_DATA
			}
		}

		// Handle out-of-bounds based on source boundary condition
		switch (from_bc) {
			case NodeType::PERIODIC: {
				auto [wrapped_row, wrapped_col] =
					apply_periodic_bc(target_row, target_col);
				size_t wrapped_idx = to_1d(wrapped_row, wrapped_col);
				NodeType wrapped_bc = get_boundary_type(wrapped_idx);

				if (wrapped_bc != NodeType::NO_DATA) {
					return FlowTransfer(
						wrapped_idx, static_cast<T>(1.0), true, false, false);
				} else {
					return FlowTransfer(); // Flow lost
				}
			}
			case NodeType::REFLECT: {
				return FlowTransfer(
					from_index, static_cast<T>(-1.0), false, true, false);
			}
			case NodeType::HAS_TO_OUT: {
				return FlowTransfer(); // Flow exits domain
			}
			case NodeType::CAN_OUT: {
				return FlowTransfer(); // Flow exits domain
			}
			default:
				return FlowTransfer(); // Flow exits domain
		}
	}
};

// Type aliases for common specializations
using FastConnectorD4_F32 = FastConnector<float, ConnectivityType::D4>;
using FastConnectorD4_F64 = FastConnector<double, ConnectivityType::D4>;
using FastConnectorD8_F32 = FastConnector<float, ConnectivityType::D8>;
using FastConnectorD8_F64 = FastConnector<double, ConnectivityType::D8>;

} // namespace dagger2
