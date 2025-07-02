#pragma once

#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include "dg2_hydraulic_erosion.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace dagger2 {

// ==============================================
// EROSION PARAMETERS BINDINGS
// ==============================================

template<typename T>
void
bind_erosion_params(py::module& m, const std::string& suffix)
{
	using ParamsType = typename HydraulicErosionEngine<T>::ErosionParams;

	py::class_<ParamsType>(m, ("ErosionParams" + suffix).c_str())
		.def(py::init<>())
		.def_readwrite("numDroplets", &ParamsType::numDroplets)
		.def_readwrite("inertia", &ParamsType::inertia)
		.def_readwrite("capacity", &ParamsType::capacity)
		.def_readwrite("deposition", &ParamsType::deposition)
		.def_readwrite("erosion", &ParamsType::erosion)
		.def_readwrite("evaporation", &ParamsType::evaporation)
		.def_readwrite("radius", &ParamsType::radius)
		.def_readwrite("minSlope", &ParamsType::minSlope)
		.def_readwrite("maxLifetime", &ParamsType::maxLifetime)
		.def_readwrite("initialWater", &ParamsType::initialWater)
		.def_readwrite("initialSpeed", &ParamsType::initialSpeed)
		.def("__repr__", [](const ParamsType& p) {
			return "ErosionParams(numDroplets=" + std::to_string(p.numDroplets) +
						 ", erosion=" + std::to_string(p.erosion) +
						 ", deposition=" + std::to_string(p.deposition) + ")";
		});
}

// ==============================================
// HYDRAULIC EROSION ENGINE BINDINGS
// ==============================================

template<typename T>
void
bind_hydraulic_erosion_engine(py::module& m, const std::string& suffix)
{
	using EngineType = HydraulicErosionEngine<T>;
	using ParamsType = typename EngineType::ErosionParams;
	using GridType = Grid2D<T>;
	using ConnectorType = Connector<T>;

	py::class_<EngineType>(m, ("HydraulicErosionEngine" + suffix).c_str())
		// Constructors
		.def(py::init<std::shared_ptr<ConnectorType>>(), py::arg("connector"))
		.def(py::init<std::shared_ptr<ConnectorType>, unsigned int>(),
				 py::arg("connector"),
				 py::arg("seed"))

		// Parameter management
		.def("setParams",
				 &EngineType::setParams,
				 "Set erosion parameters",
				 py::arg("params"))
		.def("getParams",
				 &EngineType::getParams,
				 py::return_value_policy::reference_internal,
				 "Get current erosion parameters")

		// Main erosion functions
		.def("erode",
				 &EngineType::erode,
				 "Run hydraulic erosion simulation (in-place)",
				 py::arg("heightmap"))
		.def("erode_copy",
				 &EngineType::erode_copy,
				 "Run erosion and return new heightmap",
				 py::arg("heightmap"))

		// Parameter shortcuts (removed - not in new API)

		// Presets
		.def("preset_light_erosion",
				 &EngineType::preset_light_erosion,
				 "Apply light erosion preset")
		.def("preset_moderate_erosion",
				 &EngineType::preset_moderate_erosion,
				 "Apply moderate erosion preset")
		.def("preset_heavy_erosion",
				 &EngineType::preset_heavy_erosion,
				 "Apply heavy erosion preset")

		// String representation
		.def("__repr__", [](const EngineType& e) {
			const auto& params = e.getParams();
			return "HydraulicErosionEngine(droplets=" +
						 std::to_string(params.numDroplets) +
						 ", erosion=" + std::to_string(params.erosion) +
						 ", deposition=" + std::to_string(params.deposition) + ")";
		});
}

// ==============================================
// CONVENIENCE FUNCTIONS
// ==============================================

template<typename T>
void
bind_erosion_convenience_functions(py::module& m, const std::string& suffix)
{
	using GridType = Grid2D<T>;
	using ConnectorType = Connector<T>;
	using EngineType = HydraulicErosionEngine<T>;

	// Quick erosion function
	m.def(("quick_erosion" + suffix).c_str(),
				[](std::shared_ptr<GridType> heightmap,
					 std::shared_ptr<ConnectorType> connector,
					 int droplets,
					 float intensity) -> std::shared_ptr<GridType> {
					EngineType engine(connector);
					engine.preset_moderate_erosion();
					auto params = engine.getParams();
					params.numDroplets = droplets;
					params.erosion *= intensity;
					params.deposition *= intensity;
					engine.setParams(params);
					return engine.erode_copy(heightmap);
				},
				"Quick erosion with basic parameters",
				py::arg("heightmap"),
				py::arg("connector"),
				py::arg("droplets") = 50000,
				py::arg("intensity") = 1.0f);

	// Batch erosion with different intensities
	m.def(("batch_erosion" + suffix).c_str(),
				[](std::shared_ptr<GridType> heightmap,
					 std::shared_ptr<ConnectorType> connector,
					 const std::vector<float>& intensities)
					-> std::vector<std::shared_ptr<GridType>> {
					std::vector<std::shared_ptr<GridType>> results;
					results.reserve(intensities.size());

					for (float intensity : intensities) {
						EngineType engine(connector);
						engine.preset_moderate_erosion();
						auto params = engine.getParams();
						params.erosion *= intensity;
						params.deposition *= intensity;
						engine.setParams(params);
						results.push_back(engine.erode_copy(heightmap));
					}
					return results;
				},
				"Run erosion with multiple intensity levels",
				py::arg("heightmap"),
				py::arg("connector"),
				py::arg("intensities"));

	// Progressive erosion (multiple passes)
	m.def(("progressive_erosion" + suffix).c_str(),
				[](std::shared_ptr<GridType> heightmap,
					 std::shared_ptr<ConnectorType> connector,
					 int passes,
					 float intensity_per_pass) -> std::shared_ptr<GridType> {
					auto current = heightmap;

					for (int i = 0; i < passes; ++i) {
						EngineType engine(connector);
						engine.preset_light_erosion();
						auto params = engine.getParams();
						params.erosion *= intensity_per_pass;
						params.deposition *= intensity_per_pass;
						engine.setParams(params);
						current = engine.erode_copy(current);
					}
					return current;
				},
				"Apply erosion in multiple progressive passes",
				py::arg("heightmap"),
				py::arg("connector"),
				py::arg("passes") = 3,
				py::arg("intensity_per_pass") = 0.5f);
}

// ==============================================
// PRESET PARAMETER FUNCTIONS
// ==============================================

template<typename T>
void
bind_erosion_presets(py::module& m, const std::string& suffix)
{
	using ParamsType = typename HydraulicErosionEngine<T>::ErosionParams;

	// Preset parameter generators
	m.def(("create_light_erosion_params" + suffix).c_str(),
				[]() -> ParamsType {
					ParamsType params;
					params.numDroplets = 25000;
					params.erosion = 0.3f;
					params.deposition = 0.05f;
					params.capacity = 4.0f;
					return params;
				},
				"Create light erosion parameters");

	m.def(("create_moderate_erosion_params" + suffix).c_str(),
				[]() -> ParamsType {
					ParamsType params;
					params.numDroplets = 50000;
					params.erosion = 0.8f;
					params.deposition = 0.1f;
					params.capacity = 8.0f;
					return params;
				},
				"Create moderate erosion parameters");

	m.def(("create_heavy_erosion_params" + suffix).c_str(),
				[]() -> ParamsType {
					ParamsType params;
					params.numDroplets = 100000;
					params.erosion = 1.2f;
					params.deposition = 0.15f;
					params.capacity = 12.0f;
					return params;
				},
				"Create heavy erosion parameters");

	m.def(
		("create_custom_erosion_params" + suffix).c_str(),
		[](int droplets, float erosion_rate, float deposition_rate) -> ParamsType {
			ParamsType params;
			params.numDroplets = droplets;
			params.erosion = erosion_rate;
			params.deposition = deposition_rate;
			return params;
		},
		"Create custom erosion parameters",
		py::arg("droplets"),
		py::arg("erosion_rate"),
		py::arg("deposition_rate"));
}

// ==============================================
// CONVENIENCE BINDING FUNCTIONS
// ==============================================

template<typename T>
void
bind_all_hydraulic_erosion_types(py::module& m, const std::string& suffix)
{
	bind_erosion_params<T>(m, suffix);
	bind_hydraulic_erosion_engine<T>(m, suffix);
	bind_erosion_convenience_functions<T>(m, suffix);
	bind_erosion_presets<T>(m, suffix);
}

// ==============================================
// MAIN BINDING FUNCTION
// ==============================================

inline void
bind_hydraulic_erosion_module(py::module& m)
{
	// Bind template classes for different types
	bind_all_hydraulic_erosion_types<float>(m, "F32");
	bind_all_hydraulic_erosion_types<double>(m, "F64");

	// Module documentation
	m.doc() = R"pbdoc(
        🌊 DAGGER2 HYDRAULIC EROSION MODULE 🌊
        ====================================

        Advanced particle-based hydraulic erosion engine with:
        - Framework-integrated boundary condition handling
        - Physical sediment transport and deposition
        - Configurable erosion parameters and presets
        - Batch processing and progressive erosion support
        - High-quality terrain generation for rendering

        Key Features:
        - Engine bound to Connector for proper boundary handling
        - Function-level elevation input via shared_ptr<Grid2D>
        - Sub-pixel accuracy with bilinear interpolation
        - Physically-based sediment capacity calculation
        - Multiple erosion presets for different scenarios
    )pbdoc";
}

} // namespace dagger2
