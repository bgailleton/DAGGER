#ifndef wrap_helper_HPP
#define wrap_helper_HPP


// STL imports
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <ctime>
#include <fstream>
#include <queue>
#include <stack>
#include <iostream>
#include <numeric>
#include <cmath>
#include <initializer_list>
#include <thread>
#include<stdlib.h>
#include<ctime>

#include "utils.hpp"


// defines all the format_input depnding on the eventual wrapper
#ifdef DAGGER_FT_PYTHON
#include "wrap_helper_python.hpp"
#elif DAGGER_FT_JULIA
#include "wrap_helper_julia.hpp"
#else
#include "wrap_helper_cpp.hpp"
#endif



namespace DAGGER
{



	template<typename in_t> 
	auto format_input(in_t& tin){auto ret = _format_input(tin); return ret;}

	template<class in_t, class out_t> 
	out_t format_output(in_t& tin){out_t ret = _format_output(tin); return ret;}






























}

















































#endif