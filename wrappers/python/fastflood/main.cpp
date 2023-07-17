#include "D4connector.hpp"
#include "D8connector.hpp"
#include "fastflood.hpp"
#include "fastflood_recorder.hpp"
#include "graph.hpp"
#include "hillshading.hpp"
#include "wrap_helper_python.hpp"
#include <pybind11/pybind11.h>

#ifdef OPENMP_YOLO
#include <omp.h>
#endif

namespace py = pybind11;

namespace fastflood {

template <typename CONNECTOR_T>
void declare_ff(py::module &m, std::string typestr) {
  py::class_<fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                       DAGGER::numvec<double>>>(m, typestr.c_str())
      .def(py::init<DAGGER::graph<double, CONNECTOR_T> &, CONNECTOR_T &,
                    py::array_t<double, 1> &, py::array_t<double, 1> &>())
      .def("run_SFD", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>,
                                 CONNECTOR_T, DAGGER::numvec<double>>::run_SFD)
      .def("run_SFD_with_erosion",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::run_SFD_with_erosion)
      .def("run_MFD_erosion",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::run_MFD_erosion)
      .def("run_MFD_erosion_B",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::run_MFD_erosion_B)
      .def("run_MFD", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>,
                                 CONNECTOR_T, DAGGER::numvec<double>>::run_MFD)
      .def("run_MFD_dynamic",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::run_MFD_dynamic)
      .def("run_MFD_exp",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::run_MFD_exp)
      .def("get_hw",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template get_hw<py::array>)
      .def("get_spatial_dts",
           &fastflood<
               double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
               DAGGER::numvec<double>>::template get_spatial_dts<py::array>)
      .def("get_Qwin",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template get_Qwin<py::array>)
      .def("get_Qwout",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template get_Qwout<py::array>)
      .def("get_Qs",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template get_Qs<py::array>)
      .def("get_topography",
           &fastflood<
               double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
               DAGGER::numvec<double>>::template get_topography<py::array>)
      .def("add_to_hw",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::add_to_hw)
      .def("set_Qbase", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>,
                                   CONNECTOR_T, DAGGER::numvec<double>>::
                            template set_Qbase<py::array_t<double, 1>>)
      .def("set_Qs_entry_points",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::
               template set_Qs_entry_points<py::array_t<double, 1>,
                                            py::array_t<int, 1>>)
      .def("increment_hw_from_Qbase",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::increment_hw_from_Qbase)
      .def("caesar_lisflood",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::caesar_lisflood)
      .def("set_topological_number",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_topological_number)
      .def("basicFloodos",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::basicFloodos)
      .def("basicFloodos_v2",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::basicFloodos_v2)
      .def("basicFloodos_v3",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::basicFloodos_v3)
      .def("fill_up", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>,
                                 CONNECTOR_T, DAGGER::numvec<double>>::fill_up)
      .def("set_manning",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_mannings)
      .def("testDebugWalk",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::testDebugWalk)
      .def("set_parting_coeff",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_parting_coeff)
      .def("check_SD_val",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::check_SD_val)
      .def("set_out_boundaries_to_permissive",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_out_boundaries_to_permissive)
      .def("set_edges_to_0",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_edges_to_0)
      .def("get_a_eff",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template get_a_eff<py::array>)
      .def("get_w_eff",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template get_w_eff<py::array>)
      .def("get_hydraulic_slope_D8",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::
               template get_hydraulic_slope_D8<py::array>)
      .def("spatial_dt", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>,
                                    CONNECTOR_T, DAGGER::numvec<double>>::
                             template spatial_dt<py::array_t<double, 1>>)
      .def("set_dt", &fastflood<double, DAGGER::graph<double, CONNECTOR_T>,
                                CONNECTOR_T, DAGGER::numvec<double>>::set_dt)
      .def("enable_Afdt",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::enable_Afdt)
      .def("disable_Afdt",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::disable_Afdt)
      .def("config_Afdt",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::config_Afdt)
      .def("enable_hflow",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::enable_hflow)
      .def("disable_hflow",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::disable_hflow)
      .def("set_sensibility_to_flowdepth",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_sensibility_to_flowdepth)
      .def("get_sensibility_to_flowdepth",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::get_sensibility_to_flowdepth)
      .def("fill_topo",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::fill_topo)
      .def("set_stochaslope",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::set_stochaslope)

#ifdef OPENMP_YOLO
      .def("check_devices",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template check_devices,
           py::call_guard<py::gil_scoped_release>())
      .def("caesar_lisflood_OMP",
           &fastflood<double, DAGGER::graph<double, CONNECTOR_T>, CONNECTOR_T,
                      DAGGER::numvec<double>>::template caesar_lisflood_OMP,
           py::call_guard<py::gil_scoped_release>())
#endif

      ;
}

PYBIND11_MODULE(fastflood_dagger, m) {
  m.doc() = R"pbdoc(
      fastflood-dagger example plugin
      -----------------------
      .. currentmodule:: fastflood-dagger
      .. autosummary::
         :toctree: _generate

  )pbdoc";

  declare_ff<DAGGER::D8connector<double>>(m, "FF");
  declare_ff<DAGGER::D4connector<double>>(m, "FFD4");

  py::class_<recorder<double>>(m, "recorder")
      .def(py::init<>())
      .def("enable_edot_recording", &recorder<double>::enable_edot_recording)
      .def("disable_edot_recording", &recorder<double>::disable_edot_recording)
      .def("get_edot", &recorder<double>::get_edot<py::array_t<double, 1>>)
      .def("enable_ddot_recording", &recorder<double>::enable_ddot_recording)
      .def("disable_ddot_recording", &recorder<double>::disable_ddot_recording)
      .def("get_ddot", &recorder<double>::get_ddot<py::array_t<double, 1>>)
      .def("enable_lateral_edot_recording",
           &recorder<double>::enable_lateral_edot_recording)
      .def("disable_lateral_edot_recording",
           &recorder<double>::disable_lateral_edot_recording)
      .def("get_lateral_edot",
           &recorder<double>::get_lateral_edot<py::array_t<double, 1>>)
      .def("enable_lateral_ddot_recording",
           &recorder<double>::enable_lateral_ddot_recording)
      .def("disable_lateral_ddot_recording",
           &recorder<double>::disable_lateral_ddot_recording)
      .def("get_lateral_ddot",
           &recorder<double>::get_lateral_ddot<py::array_t<double, 1>>)
      .def("enable_qs_recording", &recorder<double>::enable_qs_recording)
      .def("disable_qs_recording", &recorder<double>::disable_qs_recording)
      .def("get_qs", &recorder<double>::get_qs<py::array_t<double, 1>>)
      .def("enable_dhw_recording", &recorder<double>::enable_dhw_recording)
      .def("disable_dhw_recording", &recorder<double>::disable_dhw_recording)
      .def("get_dhw", &recorder<double>::get_dhw<py::array_t<double, 1>>);
};

}; // namespace fastflood
