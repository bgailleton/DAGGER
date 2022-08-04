//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wrap_helper_python_HPP
#define wrap_helper_python_HPP

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

// If compiling for web (using emscripten)
#ifdef __EMSCRIPTEN__
	#include <emscripten.h>
// -> Else I assume you are compiling for python (for c++ only use I'll make a separate file, just need to remvove these lines)
#else
	#include <pybind11/pybind11.h>
	#include <pybind11/stl.h>
	#include <pybind11/numpy.h>
#endif

#ifdef __EMSCRIPTEN__
#else
	namespace py = pybind11;
#endif





template<class T>
class numvec
{
public:

	T* ptr = nullptr;
	int isize = 0;
	size_t usize = 0;


	numvec(){};
	numvec(py::array_t<T,1>& arr)
	{
		auto boeuf_heure = arr.request();
		this->ptr = (T *)boeuf_heure.ptr;
		this->isize = arr.size();
		this->usize = arr.size();
	};

  T & operator [](int i) {return this->ptr[i];}

  void set(int i, T val)
  {(*this)[i] = val;}
  T get(int i)
  {return(*this)[i];}

  std::vector<T> to_vec()
  {
  	std::vector<T> out(this->isize);
  	for(int i=0;i<this->isize; ++i)
  		out[i] = (*this)[i];
  }

  // size_t size(){return this->usize;}
  const size_t size() {return this->usize;}
  


	
};


template<class T>
py::array format_output(std::vector<T>& yolo){return py::array(yolo.size(), yolo.data());}
template<class T>
py::array format_output(pvector<T>& yolo){return py::array(yolo.data->size(), yolo.data->data());}
template<class T>
py::array format_output(numvec<T>& yolo){auto vec = yolo.to_vec();  return py::array(vec.size(), vec.data());}
py::array format_output(py::array& yolo){return yolo;}

template<class T>
std::vector<T> to_vec(numvec<T>& in)
{
	std::vector<T> out(in.size());
	for(size_t i=0;i<in.size(); ++i)
		out[i] = in[i];
	return out;
}

template<class T>
std::vector<T> to_vec(py::array_t<T,1>& tin)
{
	numvec<double> in(tin);
	std::vector<T> out(in.size());
	for(size_t i=0;i<in.size(); ++i)
		out[i] = in[i];
	return out;
}


template<class T>
numvec<T> format_input(py::array_t<T,1>& yolo)
{
	numvec<T> gabul(yolo);
	return gabul;
}

template<class T>
numvec<T> format_input(numvec<T>& yolo)
{
	return yolo;
}















#endif