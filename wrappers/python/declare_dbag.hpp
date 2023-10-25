#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_dbag(py::module& m)
{
	py::class_<Hermes<int, double>>(m, "Hermes")
		.def(py::init<>())
		.def("set_surface",
				 &Hermes<int, double>::set_surface<py::array_t<double, 1>>)
		.def("get_surface",
				 &Hermes<int, double>::get_surface<py::array_t<double, 1>>)
		.def("set_hw", &Hermes<int, double>::set_hw<py::array_t<double, 1>>)
		.def("get_hw", &Hermes<int, double>::get_hw<py::array_t<double, 1>>)
		.def("set_Qwin", &Hermes<int, double>::set_Qwin<py::array_t<double, 1>>)
		.def("get_Qwin", &Hermes<int, double>::get_Qwin<py::array_t<double, 1>>)
		.def("set_Qwout", &Hermes<int, double>::set_Qwout<py::array_t<double, 1>>)
		.def("get_Qwout", &Hermes<int, double>::get_Qwout<py::array_t<double, 1>>)

		.def("set_boundaries",
				 &Hermes<int, double>::set_boundaries<py::array_t<std::uint8_t, 1>>)
		.def("get_boundaries",
				 &Hermes<int, double>::get_boundaries<py::array_t<std::uint8_t, 1>>)

		;
}
