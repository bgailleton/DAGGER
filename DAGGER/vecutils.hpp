#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

namespace VECD {

// Take a cartesian vector and returns its magnitude
template<class f_t>
f_t
vec_magnitude(f_t x, f_t y)
{
	return std::sqrt(std::pow(x, 2) + std::pow(y, 2));
}

// transform a cartesian vector to a unit one
template<class f_t>
void
unit_vector(f_t x, f_t y, f_t& xu, f_t& yu)
{
	f_t mag = vec_magnitude(x, y);
	xu = x / mag;
	yu = y / mag;
}

template<class f_t>
void
vec_addition(f_t x1, f_t y1, f_t x2, f_t y2, f_t& xr, f_t& yr)
{
	xr = x1 + x2;
	yr = y1 + y2;
}

template<class f_t>
void
vec_soustraction(f_t x1, f_t y1, f_t x2, f_t y2, f_t& xr, f_t& yr)
{
	xr = x1 - x2;
	yr = y1 - y2;
}

// Transform a cartesian vector into a polar one
template<class f_t>
void
cart2pol(f_t x, f_t y, f_t& mag, f_t& theta)
{
	mag = vec_magnitude(x, y);
	theta = std::atan2(y, x);
}

// Transform a unit cartesian vector into a polar one
template<class f_t>
void
cart2pol(f_t x, f_t y, f_t& theta)
{
	theta = std::atan2(y, x);
}

};
