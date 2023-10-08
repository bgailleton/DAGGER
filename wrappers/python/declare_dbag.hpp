#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_dbag(py::module& m)
{
	py::class_<Hermes<int, double>>(m, "Hermes").def(py::init<>());
}
