#pragma once

#include "dg2_visu_helper.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace dagger2 {

// ==============================================
// VISUALIZATION FUNCTIONS BINDINGS
// ==============================================

template<typename T>
void
bind_visualization_functions(py::module& m, const std::string& suffix)
{
	// ======================
	// HILLSHADE RELIEF
	// ======================

	m.def(("create_d8_hillshade_relief" + suffix).c_str(),
				&create_d8_hillshade_relief<T>,
				"Generate D8 hillshade relief using Horn's algorithm with connector "
				"awareness",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("sun_azimuth") = static_cast<T>(315.0),
				py::arg("sun_elevation") = static_cast<T>(45.0),
				py::arg("z_factor") = static_cast<T>(1.0),
				py::arg("cell_size") = static_cast<T>(1.0),
				py::arg("no_data_value") = static_cast<T>(-1.0));

	// ======================
	// TERRAIN ANALYSIS
	// ======================

	m.def(("create_slope_grid" + suffix).c_str(),
				&create_slope_grid<T>,
				"Calculate slope magnitude using D8 connectivity",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("cell_size") = static_cast<T>(1.0),
				py::arg("return_degrees") = true,
				py::arg("no_data_value") = static_cast<T>(-1.0));

	m.def(("create_aspect_grid" + suffix).c_str(),
				&create_aspect_grid<T>,
				"Calculate aspect (slope direction) using D8 connectivity",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("cell_size") = static_cast<T>(1.0),
				py::arg("no_data_value") = static_cast<T>(-1.0));

	m.def(
		("create_drainage_density_grid" + suffix).c_str(),
		&create_drainage_density_grid<T>,
		"Create a simple flow accumulation proxy based on local drainage density",
		py::arg("elevation"),
		py::arg("connector"),
		py::arg("no_data_value") = static_cast<T>(-1.0));

	m.def(("create_roughness_grid" + suffix).c_str(),
				&create_roughness_grid<T>,
				"Calculate terrain roughness index using D8 neighborhood",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("no_data_value") = static_cast<T>(-1.0));

	m.def(("create_plan_curvature_grid" + suffix).c_str(),
				&create_plan_curvature_grid<T>,
				"Calculate plan curvature (contour curvature) using D8 connectivity",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("cell_size") = static_cast<T>(1.0),
				py::arg("no_data_value") = static_cast<T>(-1.0));

	// ======================
	// CREATIVE VISUALIZATION
	// ======================

	m.def(("create_erosion_pattern_grid" + suffix).c_str(),
				&create_erosion_pattern_grid<T>,
				"Generate procedural erosion pattern grid",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("intensity") = static_cast<T>(1.0),
				py::arg("seed") = 12345u,
				py::arg("no_data_value") = static_cast<T>(-1.0));

	m.def(("create_exposure_grid" + suffix).c_str(),
				&create_exposure_grid<T>,
				"Generate a visibility/exposure grid (how 'exposed' each cell is)",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("no_data_value") = static_cast<T>(-1.0));

	m.def(("create_shelter_grid" + suffix).c_str(),
				&create_shelter_grid<T>,
				"Generate a shelter/protection grid (inverse of exposure)",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("no_data_value") = static_cast<T>(-1.0));
}

// ==============================================
// PRESET VISUALIZATION CONFIGURATIONS
// ==============================================

template<typename T>
void
bind_visualization_presets(py::module& m, const std::string& suffix)
{
	// Hillshade presets for different lighting conditions
	m.def(
		("create_morning_hillshade" + suffix).c_str(),
		[](const Grid2D<T>& elevation, const Connector<T>& connector) -> Grid2D<T> {
			return create_d8_hillshade_relief(
				elevation,
				connector,
				static_cast<T>(60.0), // Morning sun from NE
				static_cast<T>(30.0), // Low angle
				static_cast<T>(2.0)); // Enhanced relief
		},
		"Create morning hillshade with low-angle eastern lighting",
		py::arg("elevation"),
		py::arg("connector"));

	m.def(
		("create_noon_hillshade" + suffix).c_str(),
		[](const Grid2D<T>& elevation, const Connector<T>& connector) -> Grid2D<T> {
			return create_d8_hillshade_relief(elevation,
																				connector,
																				static_cast<T>(180.0), // South
																				static_cast<T>(70.0),	 // High angle
																				static_cast<T>(1.0));	 // Normal relief
		},
		"Create noon hillshade with high-angle southern lighting",
		py::arg("elevation"),
		py::arg("connector"));

	m.def(
		("create_dramatic_hillshade" + suffix).c_str(),
		[](const Grid2D<T>& elevation, const Connector<T>& connector) -> Grid2D<T> {
			return create_d8_hillshade_relief(
				elevation,
				connector,
				static_cast<T>(315.0), // NW (classic)
				static_cast<T>(25.0),	 // Very low angle
				static_cast<T>(3.0));	 // High relief exaggeration
		},
		"Create dramatic hillshade with very low-angle lighting and high relief",
		py::arg("elevation"),
		py::arg("connector"));

	// Combined terrain analysis
	m.def(
		("create_terrain_analysis_suite" + suffix).c_str(),
		[](const Grid2D<T>& elevation, const Connector<T>& connector) -> py::dict {
			py::dict result;
			result["hillshade"] = create_d8_hillshade_relief(elevation, connector);
			result["slope"] = create_slope_grid(elevation, connector);
			result["aspect"] = create_aspect_grid(elevation, connector);
			result["roughness"] = create_roughness_grid(elevation, connector);
			result["curvature"] = create_plan_curvature_grid(elevation, connector);
			return result;
		},
		"Create a complete terrain analysis suite with all basic grids",
		py::arg("elevation"),
		py::arg("connector"));

	// Creative visualization suite
	m.def(("create_artistic_terrain_suite" + suffix).c_str(),
				[](const Grid2D<T>& elevation,
					 const Connector<T>& connector,
					 uint32_t seed = 42) -> py::dict {
					py::dict result;
					result["dramatic_hillshade"] =
						create_d8_hillshade_relief(elevation,
																			 connector,
																			 static_cast<T>(315.0),
																			 static_cast<T>(25.0),
																			 static_cast<T>(3.0));
					result["erosion_pattern"] = create_erosion_pattern_grid(
						elevation, connector, static_cast<T>(1.5), seed);
					result["exposure"] = create_exposure_grid(elevation, connector);
					result["shelter"] = create_shelter_grid(elevation, connector);
					result["drainage"] =
						create_drainage_density_grid(elevation, connector);
					return result;
				},
				"Create an artistic terrain visualization suite for creative mapping",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("seed") = 42u);

	// Multi-angle hillshade composite
	m.def(
		("create_multi_angle_hillshade" + suffix).c_str(),
		[](const Grid2D<T>& elevation, const Connector<T>& connector) -> py::dict {
			py::dict result;
			result["northwest"] = create_d8_hillshade_relief(
				elevation, connector, static_cast<T>(315.0), static_cast<T>(45.0));
			result["northeast"] = create_d8_hillshade_relief(
				elevation, connector, static_cast<T>(45.0), static_cast<T>(45.0));
			result["southwest"] = create_d8_hillshade_relief(
				elevation, connector, static_cast<T>(225.0), static_cast<T>(45.0));
			result["southeast"] = create_d8_hillshade_relief(
				elevation, connector, static_cast<T>(135.0), static_cast<T>(45.0));
			return result;
		},
		"Create hillshades from four cardinal directions for composite effects",
		py::arg("elevation"),
		py::arg("connector"));
}

// ==============================================
// ADVANCED VISUALIZATION FUNCTIONS
// ==============================================

template<typename T>
void
bind_advanced_visualization(py::module& m, const std::string& suffix)
{
	// Enhanced hillshade with multiple light sources
	m.def(("create_enhanced_hillshade" + suffix).c_str(),
				[](const Grid2D<T>& elevation,
					 const Connector<T>& connector,
					 const std::vector<T>& azimuths,
					 const std::vector<T>& elevations,
					 const std::vector<T>& weights) -> Grid2D<T> {
					if (azimuths.size() != elevations.size() ||
							azimuths.size() != weights.size()) {
						throw std::invalid_argument(
							"Azimuth, elevation, and weight vectors must have same size");
					}

					if (azimuths.empty()) {
						throw std::invalid_argument(
							"Must provide at least one light source");
					}

					const size_t rows = elevation.rows();
					const size_t cols = elevation.cols();
					std::vector<T> combined_data(rows * cols, static_cast<T>(0.0));

					T total_weight = static_cast<T>(0.0);
					for (const auto& weight : weights) {
						total_weight += weight;
					}

					// Combine multiple hillshades
					for (size_t i = 0; i < azimuths.size(); ++i) {
						auto hillshade = create_d8_hillshade_relief(
							elevation, connector, azimuths[i], elevations[i]);

						T normalized_weight = weights[i] / total_weight;

						for (size_t j = 0; j < combined_data.size(); ++j) {
							combined_data[j] += hillshade[j] * normalized_weight;
						}
					}

					return Grid2D<T>(combined_data, rows, cols);
				},
				"Create enhanced hillshade with multiple weighted light sources",
				py::arg("elevation"),
				py::arg("connector"),
				py::arg("azimuths"),
				py::arg("elevations"),
				py::arg("weights"));

	// Terrain classification grid
	m.def(("create_terrain_classification_grid" + suffix).c_str(),
				[](const Grid2D<T>& elevation,
					 const Connector<T>& connector) -> Grid2D<int> {
					auto slope = create_slope_grid(elevation, connector);
					auto curvature = create_plan_curvature_grid(elevation, connector);
					auto exposure = create_exposure_grid(elevation, connector);

					const size_t rows = elevation.rows();
					const size_t cols = elevation.cols();
					std::vector<int> classification_data(rows * cols);

					for (size_t i = 0; i < rows * cols; ++i) {
						if (connector.get_boundary_type(i) == NodeType::NO_DATA) {
							classification_data[i] = -1;
							continue;
						}

						T s = slope[i];
						T c = curvature[i];
						T e = exposure[i];

						// Simple terrain classification
						if (s < static_cast<T>(5.0)) {
							classification_data[i] = 0; // Flat
						} else if (s < static_cast<T>(15.0)) {
							if (c > static_cast<T>(0.1)) {
								classification_data[i] = 1; // Gentle convex slope
							} else if (c < static_cast<T>(-0.1)) {
								classification_data[i] = 2; // Gentle concave slope
							} else {
								classification_data[i] = 3; // Gentle planar slope
							}
						} else if (s < static_cast<T>(35.0)) {
							if (e > static_cast<T>(10.0)) {
								classification_data[i] = 4; // Exposed moderate slope
							} else {
								classification_data[i] = 5; // Sheltered moderate slope
							}
						} else {
							classification_data[i] = 6; // Steep slope
						}
					}

					return Grid2D<int>(classification_data, rows, cols);
				},
				"Create terrain classification grid based on slope, curvature, and "
				"exposure",
				py::arg("elevation"),
				py::arg("connector"));
}

// ==============================================
// COMPREHENSIVE BINDING FUNCTION
// ==============================================

template<typename T>
void
bind_all_visualization_types(py::module& m, const std::string& suffix)
{
	// Core visualization functions
	bind_visualization_functions<T>(m, suffix);

	// Preset configurations
	bind_visualization_presets<T>(m, suffix);

	// Advanced functions
	bind_advanced_visualization<T>(m, suffix);
}

// ==============================================
// MAIN BINDING FUNCTION
// ==============================================

inline void
bind_visualization_module(py::module& m)
{
	// Bind for different numeric types
	bind_all_visualization_types<float>(m, "F32");
	bind_all_visualization_types<double>(m, "F64");

	// Module documentation
	m.doc() = R"pbdoc(
        🏔️ DAGGER2 TERRAIN VISUALIZATION TOOLKIT 🏔️
        ============================================

        High-performance C++ backend for generating raw Grid2D visualization data!

        🎨 CORE FEATURES:
        • D8 Hillshade Relief using Horn's algorithm with connector awareness
        • Slope and Aspect analysis with proper boundary condition handling
        • Terrain Roughness and Curvature calculations
        • Drainage Density and Flow Accumulation proxies

        🌟 CREATIVE VISUALIZATIONS:
        • Erosion Pattern generation with procedural noise
        • Exposure and Shelter analysis for lighting effects
        • Multi-angle hillshade compositing
        • Enhanced hillshade with multiple light sources
        • Automatic terrain classification

        🚀 USAGE EXAMPLES:
        ```python
        import dagger2

        # Basic hillshade
        hillshade = dagger2.create_d8_hillshade_reliefF64(elevation, connector)

        # Complete terrain analysis
        terrain_suite = dagger2.create_terrain_analysis_suiteF64(elevation, connector)

        # Artistic visualization suite
        artistic = dagger2.create_artistic_terrain_suiteF64(elevation, connector, seed=42)

        # Enhanced multi-light hillshade
        azimuths = [315.0, 45.0, 225.0]
        elevations = [45.0, 30.0, 60.0]
        weights = [0.5, 0.3, 0.2]
        enhanced = dagger2.create_enhanced_hillshadeF64(
            elevation, connector, azimuths, elevations, weights)
        ```

        ✨ Perfect for creating beautiful terrain visualizations and maps! ✨

        All functions return raw Grid2D objects ready for Python visualization!
    )pbdoc";
}

} // namespace dagger2
