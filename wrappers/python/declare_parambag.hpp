#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_param(py::module& m)
{
	py::class_<ParamBag<int, double>>(m, "ParamBag")
		.def(py::init<>())
		.def("set_ke", &ParamBag<int, double>::set_ke)
		.def("enable_gf2_diffuse_Qwin",
				 &ParamBag<int, double>::enable_gf2_diffuse_Qwin)
		.def("disable_gf2_diffuse_Qwin",
				 &ParamBag<int, double>::disable_gf2_diffuse_Qwin)
		.def("enable_gf2_morpho", &ParamBag<int, double>::enable_gf2_morpho)
		.def("disable_gf2_morpho", &ParamBag<int, double>::disable_gf2_morpho)
		.def("set_kd", &ParamBag<int, double>::set_kd)
		.def("get_kd", &ParamBag<int, double>::get_kd)
		.def("enable_TSG_dist", &ParamBag<int, double>::enable_TSG_dist)
		.def("disable_TSG_dist", &ParamBag<int, double>::disable_TSG_dist)
		.def("set_TSG_distmax", &ParamBag<int, double>::set_TSG_distmax)

		;
}
