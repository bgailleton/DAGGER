#pragma once

#include "dg2_array.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <vector>

namespace dagger2 {

/**
 * Vector2: Simple 2D vector class for gradient calculations
 */
template<typename T>
struct Vector2
{
	T x, y;

	Vector2()
		: x(0)
		, y(0)
	{
	}
	Vector2(T x_val, T y_val)
		: x(x_val)
		, y(y_val)
	{
	}

	Vector2 operator+(const Vector2& other) const
	{
		return Vector2(x + other.x, y + other.y);
	}
	Vector2 operator-(const Vector2& other) const
	{
		return Vector2(x - other.x, y - other.y);
	}
	Vector2 operator*(T scalar) const { return Vector2(x * scalar, y * scalar); }
	T dot(const Vector2& other) const { return x * other.x + y * other.y; }
	T length() const { return std::sqrt(x * x + y * y); }
	Vector2 normalize() const
	{
		T len = length();
		return len > 0 ? Vector2(x / len, y / len) : Vector2();
	}
};

/**
 * PerlinNoiseConfig: Configuration for Perlin noise generation
 */
template<typename T>
struct PerlinNoiseConfig
{
	// Basic parameters
	T frequency = static_cast<T>(1.0);	 // Base frequency
	T amplitude = static_cast<T>(1.0);	 // Base amplitude
	int octaves = 4;										 // Number of octaves for fractal noise
	T persistence = static_cast<T>(0.5); // Amplitude multiplier per octave
	T lacunarity = static_cast<T>(2.0);	 // Frequency multiplier per octave

	// Offset and scaling
	T offset_x = static_cast<T>(0.0); // X offset
	T offset_y = static_cast<T>(0.0); // Y offset
	T scale_x = static_cast<T>(1.0);	// X scale factor
	T scale_y = static_cast<T>(1.0);	// Y scale factor

	// Output range
	T min_value = static_cast<T>(0.0); // Minimum output value
	T max_value = static_cast<T>(1.0); // Maximum output value

	// Quality settings
	uint32_t permutation_size = 256; // Size of permutation table
	uint32_t seed = 12345;					 // Random seed
	bool normalize_gradients = true; // Whether to normalize gradient vectors

	// Advanced features
	bool use_ridged = false;									// Use ridged noise variant
	T ridge_threshold = static_cast<T>(0.0);	// Threshold for ridged noise
	bool use_billow = false;									// Use billow noise variant
	T turbulence_power = static_cast<T>(1.0); // Power for turbulence

	// Boundary conditions
	bool tile_x = false;								 // Tile in X direction
	bool tile_y = false;								 // Tile in Y direction
	T tile_width = static_cast<T>(1.0);	 // Width for tiling
	T tile_height = static_cast<T>(1.0); // Height for tiling
};

/**
 * PerlinNoise: Comprehensive Perlin noise generator
 *
 * This class provides high-quality Perlin noise generation with support for:
 * - Classic Perlin noise
 * - Fractal/octaved noise
 * - Ridged and billow variants
 * - Turbulence effects
 * - Seamless tiling
 * - Custom gradient interpolation
 * - Integration with Grid2D
 */
template<typename T = double>
class PerlinNoise
{
public:
	using ConfigType = PerlinNoiseConfig<T>;
	using Vec2 = Vector2<T>;

private:
	// Configuration
	ConfigType config_;

	// Permutation table for pseudo-random gradients
	std::vector<uint32_t> permutation_;
	std::vector<Vec2> gradients_;

	// Random number generator
	mutable std::mt19937 rng_;

	// Cache for performance
	mutable std::vector<T> fade_cache_;
	mutable bool cache_dirty_;

public:
	/**
	 * Constructor with configuration
	 */
	explicit PerlinNoise(const ConfigType& config = ConfigType())
		: config_(config)
		, rng_(config.seed)
		, cache_dirty_(true)
	{
		initialize_permutation();
		initialize_gradients();
		update_cache();
	}

	/**
	 * Constructor with simple parameters
	 */
	PerlinNoise(T frequency, T amplitude, int octaves = 4, uint32_t seed = 12345)
		: config_()
		, rng_(seed)
		, cache_dirty_(true)
	{
		config_.frequency = frequency;
		config_.amplitude = amplitude;
		config_.octaves = octaves;
		config_.seed = seed;
		initialize_permutation();
		initialize_gradients();
		update_cache();
	}

	/**
	 * Generate noise value at a single point
	 */
	T noise(T x, T y) const
	{
		// Apply transformations
		T transformed_x = (x + config_.offset_x) * config_.scale_x;
		T transformed_y = (y + config_.offset_y) * config_.scale_y;

		if (config_.octaves == 1) {
			return single_octave_noise(transformed_x * config_.frequency,
																 transformed_y * config_.frequency) *
						 config_.amplitude;
		} else {
			return fractal_noise(transformed_x, transformed_y);
		}
	}

	/**
	 * Generate noise for entire Grid2D
	 */
	void generate(Grid2D<T>& grid) const
	{
		generate(
			grid, 0, 0, static_cast<T>(grid.cols()), static_cast<T>(grid.rows()));
	}

	/**
	 * Generate noise for Grid2D with custom world coordinates
	 */
	void generate(Grid2D<T>& grid,
								T world_x,
								T world_y,
								T world_width,
								T world_height) const
	{
		const size_t rows = grid.rows();
		const size_t cols = grid.cols();

		const T dx = world_width / static_cast<T>(cols);
		const T dy = world_height / static_cast<T>(rows);

		std::cout << "DEBUG::5" << std::endl;

		for (size_t r = 0; r < rows; ++r) {
			for (size_t c = 0; c < cols; ++c) {
				T x = world_x + static_cast<T>(c) * dx;
				T y = world_y + static_cast<T>(r) * dy;

				T noise_value = noise(x, y);

				// Normalize to output range
				noise_value =
					config_.min_value +
					(noise_value + 1.0) * 0.5 * (config_.max_value - config_.min_value);

				grid(r, c) = noise_value;
			}
		}
	}

	/**
	 * Generate turbulence (absolute value of noise)
	 */
	T turbulence(T x, T y) const
	{
		T value = 0;
		T amplitude = config_.amplitude;
		T frequency = config_.frequency;

		for (int i = 0; i < config_.octaves; ++i) {
			T transformed_x = (x + config_.offset_x) * config_.scale_x * frequency;
			T transformed_y = (y + config_.offset_y) * config_.scale_y * frequency;

			value +=
				std::abs(single_octave_noise(transformed_x, transformed_y)) * amplitude;

			amplitude *= config_.persistence;
			frequency *= config_.lacunarity;
		}

		return std::pow(value, config_.turbulence_power);
	}

	/**
	 * Generate ridged noise (inverted absolute value)
	 */
	T ridged_noise(T x, T y) const
	{
		T value = 0;
		T amplitude = config_.amplitude;
		T frequency = config_.frequency;
		T weight = 1.0;

		for (int i = 0; i < config_.octaves; ++i) {
			T transformed_x = (x + config_.offset_x) * config_.scale_x * frequency;
			T transformed_y = (y + config_.offset_y) * config_.scale_y * frequency;

			T signal = single_octave_noise(transformed_x, transformed_y);
			signal = std::abs(signal);
			signal = config_.ridge_threshold - signal;
			signal *= signal;
			signal *= weight;

			weight = std::clamp(signal * 2.0, 0.0, 1.0);
			value += signal * amplitude;

			amplitude *= config_.persistence;
			frequency *= config_.lacunarity;
		}

		return value;
	}

	/**
	 * Generate billow noise (absolute value with different scaling)
	 */
	T billow_noise(T x, T y) const
	{
		T value = 0;
		T amplitude = config_.amplitude;
		T frequency = config_.frequency;

		for (int i = 0; i < config_.octaves; ++i) {
			T transformed_x = (x + config_.offset_x) * config_.scale_x * frequency;
			T transformed_y = (y + config_.offset_y) * config_.scale_y * frequency;

			T signal = single_octave_noise(transformed_x, transformed_y);
			signal = 2.0 * std::abs(signal) - 1.0;
			value += signal * amplitude;

			amplitude *= config_.persistence;
			frequency *= config_.lacunarity;
		}

		return value;
	}

	/**
	 * Generate derivative (gradient) information
	 */
	Vec2 gradient(T x, T y) const
	{
		const T h = static_cast<T>(0.001); // Small offset for numerical derivative

		T dx = (noise(x + h, y) - noise(x - h, y)) / (2.0 * h);
		T dy = (noise(x, y + h) - noise(x, y - h)) / (2.0 * h);

		return Vec2(dx, dy);
	}

	/**
	 * Generate noise with custom interpolation function
	 */
	template<typename InterpolationFunc>
	T noise_with_interpolation(T x, T y, InterpolationFunc interpolate) const
	{
		T transformed_x =
			(x + config_.offset_x) * config_.scale_x * config_.frequency;
		T transformed_y =
			(y + config_.offset_y) * config_.scale_y * config_.frequency;

		return single_octave_noise_custom(
						 transformed_x, transformed_y, interpolate) *
					 config_.amplitude;
	}

	/**
	 * Generate seamlessly tiling noise
	 */
	T tiling_noise(T x, T y) const
	{
		if (!config_.tile_x && !config_.tile_y) {
			return noise(x, y);
		}

		T nx = x, ny = y;

		if (config_.tile_x) {
			T s = x / config_.tile_width;
			T t = y / config_.tile_height;

			T dx1 = noise(x, y);
			T dx2 = noise(x - config_.tile_width, y);
			T dx3 = noise(x, y - config_.tile_height);
			T dx4 = noise(x - config_.tile_width, y - config_.tile_height);

			T i1 = interpolate(dx1, dx2, s);
			T i2 = interpolate(dx3, dx4, s);

			return interpolate(i1, i2, t);
		}

		return noise(nx, ny);
	}

	/**
	 * Generate noise using different variants
	 */
	T noise_variant(T x, T y) const
	{
		if (config_.use_ridged) {
			return ridged_noise(x, y);
		} else if (config_.use_billow) {
			return billow_noise(x, y);
		} else {
			return noise(x, y);
		}
	}

	/**
	 * Batch generation for multiple points
	 */
	void generate_batch(const std::vector<Vec2>& points,
											std::vector<T>& output) const
	{
		output.resize(points.size());

		for (size_t i = 0; i < points.size(); ++i) {
			output[i] = noise(points[i].x, points[i].y);
		}
	}

	/**
	 * Generate noise with custom amplitude per octave
	 */
	T noise_custom_amplitudes(T x, T y, const std::vector<T>& amplitudes) const
	{
		T value = 0;
		T frequency = config_.frequency;

		const int octaves =
			std::min(config_.octaves, static_cast<int>(amplitudes.size()));

		for (int i = 0; i < octaves; ++i) {
			T transformed_x = (x + config_.offset_x) * config_.scale_x * frequency;
			T transformed_y = (y + config_.offset_y) * config_.scale_y * frequency;

			value +=
				single_octave_noise(transformed_x, transformed_y) * amplitudes[i];
			frequency *= config_.lacunarity;
		}

		return value;
	}

	// ======================
	// CONFIGURATION METHODS
	// ======================

	/**
	 * Update configuration
	 */
	void set_config(const ConfigType& config)
	{
		bool need_reinit = (config.seed != config_.seed ||
												config.permutation_size != config_.permutation_size);

		config_ = config;
		cache_dirty_ = true;

		if (need_reinit) {
			rng_.seed(config_.seed);
			initialize_permutation();
			initialize_gradients();
		}

		update_cache();
	}

	const ConfigType& get_config() const { return config_; }

	void set_seed(uint32_t seed)
	{
		config_.seed = seed;
		rng_.seed(seed);
		initialize_permutation();
		initialize_gradients();
	}

	void set_octaves(int octaves) { config_.octaves = std::max(1, octaves); }
	void set_frequency(T frequency) { config_.frequency = frequency; }
	void set_amplitude(T amplitude) { config_.amplitude = amplitude; }
	void set_persistence(T persistence) { config_.persistence = persistence; }
	void set_lacunarity(T lacunarity) { config_.lacunarity = lacunarity; }

	void set_offset(T x, T y)
	{
		config_.offset_x = x;
		config_.offset_y = y;
	}

	void set_scale(T x, T y)
	{
		config_.scale_x = x;
		config_.scale_y = y;
	}

	void set_output_range(T min_val, T max_val)
	{
		config_.min_value = min_val;
		config_.max_value = max_val;
	}

	// ======================
	// UTILITY METHODS
	// ======================

	/**
	 * Get noise statistics for a region
	 */
	struct NoiseStats
	{
		T min_value;
		T max_value;
		T mean_value;
		T std_deviation;
		size_t sample_count;
	};

	NoiseStats analyze_region(T x,
														T y,
														T width,
														T height,
														size_t samples_x,
														size_t samples_y) const
	{
		NoiseStats stats;
		std::vector<T> values;
		values.reserve(samples_x * samples_y);

		const T dx = width / static_cast<T>(samples_x - 1);
		const T dy = height / static_cast<T>(samples_y - 1);

		for (size_t j = 0; j < samples_y; ++j) {
			for (size_t i = 0; i < samples_x; ++i) {
				T sample_x = x + static_cast<T>(i) * dx;
				T sample_y = y + static_cast<T>(j) * dy;
				values.push_back(noise(sample_x, sample_y));
			}
		}

		stats.sample_count = values.size();
		auto minmax = std::minmax_element(values.begin(), values.end());
		stats.min_value = *minmax.first;
		stats.max_value = *minmax.second;

		T sum = std::accumulate(values.begin(), values.end(), static_cast<T>(0));
		stats.mean_value = sum / static_cast<T>(values.size());

		T variance = 0;
		for (T value : values) {
			T diff = value - stats.mean_value;
			variance += diff * diff;
		}
		variance /= static_cast<T>(values.size());
		stats.std_deviation = std::sqrt(variance);

		return stats;
	}

	/**
	 * Create Grid2D with noise
	 */
	Grid2D<T> create_noise_grid(size_t rows,
															size_t cols,
															T world_x = 0,
															T world_y = 0,
															T world_width = 1,
															T world_height = 1) const
	{
		std::cout << "DEBUG::3.5" << std::endl;
		std::vector<T> data(rows * cols);
		std::cout << "DEBUG::3.6" << std::endl;
		Grid2D<T> grid(data, rows, cols);
		std::cout << "DEBUG::4" << std::endl;
		generate(grid, world_x, world_y, world_width, world_height);
		return grid;
	}

private:
	// ======================
	// INTERNAL METHODS
	// ======================

	/**
	 * Initialize permutation table
	 */
	void initialize_permutation()
	{
		permutation_.resize(config_.permutation_size);
		std::iota(permutation_.begin(), permutation_.end(), 0);
		std::shuffle(permutation_.begin(), permutation_.end(), rng_);

		// Duplicate for overflow handling
		permutation_.insert(
			permutation_.end(), permutation_.begin(), permutation_.end());
	}

	/**
	 * Initialize gradient vectors
	 */
	void initialize_gradients()
	{
		gradients_.clear();
		gradients_.reserve(config_.permutation_size);

		std::uniform_real_distribution<T> dist(-1.0, 1.0);

		for (uint32_t i = 0; i < config_.permutation_size; ++i) {
			Vec2 grad(dist(rng_), dist(rng_));
			if (config_.normalize_gradients) {
				grad = grad.normalize();
			}
			gradients_.push_back(grad);
		}
	}

	/**
	 * Update performance cache
	 */
	void update_cache()
	{
		if (!cache_dirty_)
			return;

		fade_cache_.clear();
		fade_cache_.reserve(1024);

		for (size_t i = 0; i < 1024; ++i) {
			T t = static_cast<T>(i) / 1023.0;
			fade_cache_.push_back(fade(t));
		}

		cache_dirty_ = false;
	}

	/**
	 * Smooth interpolation function (quintic)
	 */
	T fade(T t) const { return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }

	/**
	 * Linear interpolation
	 */
	T interpolate(T a, T b, T t) const { return a + t * (b - a); }

	/**
	 * Get gradient at permutation index
	 */
	Vec2 get_gradient(int x, int y) const
	{
		int index =
			permutation_[(x + permutation_[y & (config_.permutation_size - 1)]) &
									 (config_.permutation_size - 1)];
		return gradients_[index & (config_.permutation_size - 1)];
	}

	/**
	 * Single octave Perlin noise
	 */
	T single_octave_noise(T x, T y) const
	{
		// Find grid cell containing point
		int x0 = static_cast<int>(std::floor(x));
		int y0 = static_cast<int>(std::floor(y));
		int x1 = x0 + 1;
		int y1 = y0 + 1;

		// Relative coordinates within cell
		T sx = x - static_cast<T>(x0);
		T sy = y - static_cast<T>(y0);

		// Get gradients at corner points
		Vec2 g00 = get_gradient(x0, y0);
		Vec2 g10 = get_gradient(x1, y0);
		Vec2 g01 = get_gradient(x0, y1);
		Vec2 g11 = get_gradient(x1, y1);

		// Distance vectors from corners to point
		Vec2 d00(sx, sy);
		Vec2 d10(sx - 1.0, sy);
		Vec2 d01(sx, sy - 1.0);
		Vec2 d11(sx - 1.0, sy - 1.0);

		// Dot products
		T n00 = g00.dot(d00);
		T n10 = g10.dot(d10);
		T n01 = g01.dot(d01);
		T n11 = g11.dot(d11);

		// Smooth interpolation weights
		T u = fade(sx);
		T v = fade(sy);

		// Interpolate
		T nx0 = interpolate(n00, n10, u);
		T nx1 = interpolate(n01, n11, u);

		return interpolate(nx0, nx1, v);
	}

	/**
	 * Single octave noise with custom interpolation
	 */
	template<typename InterpolationFunc>
	T single_octave_noise_custom(T x,
															 T y,
															 InterpolationFunc interpolate_func) const
	{
		// Find grid cell containing point
		int x0 = static_cast<int>(std::floor(x));
		int y0 = static_cast<int>(std::floor(y));
		int x1 = x0 + 1;
		int y1 = y0 + 1;

		// Relative coordinates within cell
		T sx = x - static_cast<T>(x0);
		T sy = y - static_cast<T>(y0);

		// Get gradients at corner points
		Vec2 g00 = get_gradient(x0, y0);
		Vec2 g10 = get_gradient(x1, y0);
		Vec2 g01 = get_gradient(x0, y1);
		Vec2 g11 = get_gradient(x1, y1);

		// Distance vectors from corners to point
		Vec2 d00(sx, sy);
		Vec2 d10(sx - 1.0, sy);
		Vec2 d01(sx, sy - 1.0);
		Vec2 d11(sx - 1.0, sy - 1.0);

		// Dot products
		T n00 = g00.dot(d00);
		T n10 = g10.dot(d10);
		T n01 = g01.dot(d01);
		T n11 = g11.dot(d11);

		// Custom interpolation
		T nx0 = interpolate_func(n00, n10, sx);
		T nx1 = interpolate_func(n01, n11, sx);

		return interpolate_func(nx0, nx1, sy);
	}

	/**
	 * Fractal (multi-octave) noise
	 */
	T fractal_noise(T x, T y) const
	{
		T value = 0;
		T amplitude = config_.amplitude;
		T frequency = config_.frequency;
		T max_value = 0; // For normalization

		for (int i = 0; i < config_.octaves; ++i) {
			T transformed_x = x * frequency;
			T transformed_y = y * frequency;

			value += single_octave_noise(transformed_x, transformed_y) * amplitude;
			max_value += amplitude;

			amplitude *= config_.persistence;
			frequency *= config_.lacunarity;
		}

		return value / max_value; // Normalize
	}
};

// ======================
// CONVENIENCE FUNCTIONS
// ======================

/**
 * Quick noise generation functions
 */
template<typename T>
Grid2D<T>
generate_perlin_noise(size_t rows,
											size_t cols,
											T frequency = 1.0,
											T amplitude = 1.0,
											int octaves = 4,
											uint32_t seed = 12345)
{
	PerlinNoise<T> noise(frequency, amplitude, octaves, seed);
	return noise.create_noise_grid(rows, cols);
}

template<typename T>
void
add_perlin_noise(Grid2D<T>& grid,
								 T frequency = 1.0,
								 T amplitude = 1.0,
								 int octaves = 4,
								 uint32_t seed = 12345)
{
	PerlinNoise<T> noise(frequency, amplitude, octaves, seed);
	Grid2D<T> noise_grid = noise.create_noise_grid(grid.rows(), grid.cols());

	for (size_t i = 0; i < grid.size(); ++i) {
		grid[i] += noise_grid[i];
	}
}

template<typename T>
Grid2D<T>
generate_turbulence(size_t rows,
										size_t cols,
										T frequency = 1.0,
										T amplitude = 1.0,
										int octaves = 4,
										uint32_t seed = 12345)
{
	PerlinNoise<T> noise(frequency, amplitude, octaves, seed);

	std::vector<T> data(rows * cols);
	Grid2D<T> grid(data, rows, cols);

	const T dx = 1.0 / static_cast<T>(cols);
	const T dy = 1.0 / static_cast<T>(rows);

	for (size_t r = 0; r < rows; ++r) {
		for (size_t c = 0; c < cols; ++c) {
			T x = static_cast<T>(c) * dx;
			T y = static_cast<T>(r) * dy;
			grid(r, c) = noise.turbulence(x, y);
		}
	}

	return grid;
}

// /**
//  * Pybind11 integration helper
//  */
// template<typename T>
// void bind_perlin_noise(py::module& m, const std::string& suffix) {
//     using PerlinType = PerlinNoise<T>;
//     using ConfigType = PerlinNoiseConfig<T>;
//     using Vec2Type = Vector2<T>;
//     using StatsType = typename PerlinType::NoiseStats;

//     // Vector2
//     py::class_<Vec2Type>(m, ("Vector2" + suffix).c_str())
//         .def(py::init<>())
//         .def(py::init<T, T>())
//         .def_readwrite("x", &Vec2Type::x)
//         .def_readwrite("y", &Vec2Type::y)
//         .def("dot", &Vec2Type::dot)
//         .def("length", &Vec2Type::length)
//         .def("normalize", &Vec2Type::normalize);

//     // Configuration
//     py::class_<ConfigType>(m, ("PerlinNoiseConfig" + suffix).c_str())
//         .def(py::init<>())
//         .def_readwrite("frequency", &ConfigType::frequency)
//         .def_readwrite("amplitude", &ConfigType::amplitude)
//         .def_readwrite("octaves", &ConfigType::octaves)
//         .def_readwrite("persistence", &ConfigType::persistence)
//         .def_readwrite("lacunarity", &ConfigType::lacunarity)
//         .def_readwrite("offset_x", &ConfigType::offset_x)
//         .def_readwrite("offset_y", &ConfigType::offset_y)
//         .def_readwrite("scale_x", &ConfigType::scale_x)
//         .def_readwrite("scale_y", &ConfigType::scale_y)
//         .def_readwrite("min_value", &ConfigType::min_value)
//         .def_readwrite("max_value", &ConfigType::max_value)
//         .def_readwrite("seed", &ConfigType::seed)
//         .def_readwrite("use_ridged", &ConfigType::use_ridged)
//         .def_readwrite("use_billow", &ConfigType::use_billow)
//         .def_readwrite("tile_x", &ConfigType::tile_x)
//         .def_readwrite("tile_y", &ConfigType::tile_y);

//     // Noise statistics
//     py::class_<StatsType>(m, ("NoiseStats" + suffix).c_str())
//         .def_readonly("min_value", &StatsType::min_value)
//         .def_readonly("max_value", &StatsType::max_value)
//         .def_readonly("mean_value", &StatsType::mean_value)
//         .def_readonly("std_deviation", &StatsType::std_deviation)
//         .def_readonly("sample_count", &StatsType::sample_count);

//     // Main PerlinNoise class
//     py::class_<PerlinType>(m, ("PerlinNoise" + suffix).c_str())
//         .def(py::init<const ConfigType&>(), py::arg("config") = ConfigType())
//         .def(py::init<T, T, int, uint32_t>(),
//              py::arg("frequency"), py::arg("amplitude"),
//              py::arg("octaves") = 4, py::arg("seed") = 12345)
//         .def("noise", &PerlinType::noise)
//         .def("turbulence", &PerlinType::turbulence)
//         .def("ridged_noise", &PerlinType::ridged_noise)
//         .def("billow_noise", &PerlinType::billow_noise)
//         .def("tiling_noise", &PerlinType::tiling_noise)
//         .def("noise_variant", &PerlinType::noise_variant)
//         .def("gradient", &PerlinType::gradient)
//         .def("generate", py::overload_cast<Grid2D<T>&>(&PerlinType::generate,
//         py::const_)) .def("generate", py::overload_cast<Grid2D<T>&, T, T, T,
//         T>(&PerlinType::generate, py::const_)) .def("create_noise_grid",
//         &PerlinType::create_noise_grid,
//              py::arg("rows"), py::arg("cols"),
//              py::arg("world_x") = 0, py::arg("world_y") = 0,
//              py::arg("world_width") = 1, py::arg("world_height") = 1)
//         .def("analyze_region", &PerlinType::analyze_region)
//         .def("set_config", &PerlinType::set_config)
//         .def("get_config", &PerlinType::get_config,
//         py::return_value_policy::reference_internal) .def("set_seed",
//         &PerlinType::set_seed) .def("set_octaves", &PerlinType::set_octaves)
//         .def("set_frequency", &PerlinType::set_frequency)
//         .def("set_amplitude", &PerlinType::set_amplitude)
//         .def("set_persistence", &PerlinType::set_persistence)
//         .def("set_lacunarity", &PerlinType::set_lacunarity)
//         .def("set_offset", &PerlinType::set_offset)
//         .def("set_scale", &PerlinType::set_scale)
//         .def("set_output_range", &PerlinType::set_output_range);

//     // Convenience functions
//     m.def(("generate_perlin_noise" + suffix).c_str(),
//     &generate_perlin_noise<T>,
//           "Generate Perlin noise grid",
//           py::arg("rows"), py::arg("cols"),
//           py::arg("frequency") = 1.0, py::arg("amplitude") = 1.0,
//           py::arg("octaves") = 4, py::arg("seed") = 12345);

//     m.def(("add_perlin_noise" + suffix).c_str(), &add_perlin_noise<T>,
//           "Add Perlin noise to existing grid",
//           py::arg("grid"), py::arg("frequency") = 1.0, py::arg("amplitude")
//           = 1.0, py::arg("octaves") = 4, py::arg("seed") = 12345);

//     m.def(("generate_turbulence" + suffix).c_str(), &generate_turbulence<T>,
//           "Generate turbulence noise grid",
//           py::arg("rows"), py::arg("cols"),
//           py::arg("frequency") = 1.0, py::arg("amplitude") = 1.0,
//           py::arg("octaves") = 4, py::arg("seed") = 12345);
// }

} // namespace dagger2

/*
Example usage:

#include "dg2_perlin_noise.hpp"

// Basic usage
auto noise = dagger2::PerlinNoise<double>(1.0, 1.0, 4);
auto grid = noise.create_noise_grid(256, 256);

// Advanced configuration
dagger2::PerlinNoiseConfig<double> config;
config.frequency = 0.01;
config.amplitude = 100.0;
config.octaves = 6;
config.persistence = 0.6;
config.lacunarity = 2.1;
config.use_ridged = true;
config.seed = 42;

auto terrain_noise = dagger2::PerlinNoise<double>(config);
auto terrain = terrain_noise.create_noise_grid(512, 512, 0, 0, 1000, 1000);

// Integration with existing Grid2D
std::vector<double> data(128 * 128);
dagger2::Grid2D<double> my_grid(data, 128, 128);
terrain_noise.generate(my_grid, -500, -500, 1000, 1000);

// Turbulence for cloud-like effects
auto clouds = dagger2::generate_turbulence<float>(256, 256, 0.02f, 1.0f, 5);

// Analyze noise characteristics
auto stats = terrain_noise.analyze_region(0, 0, 100, 100, 50, 50);
std::cout << "Noise range: [" << stats.min_value << ", " << stats.max_value <<
"]" << std::endl;

// Python integration example:
PYBIND11_MODULE(landscape, m) {
		bind_array_types<float>(m, "F32");
		bind_array_types<double>(m, "F64");

		bind_perlin_noise<float>(m, "F32");
		bind_perlin_noise<double>(m, "F64");
}
*/
