#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_dbag(py::module& m)
{
	py::class_<ParamBag<int, double>>(m, "ParamBag")
		.def(py::init<>())
		.def("set_ke", &ParamBag<int, double>::set_ke)

		;
}
