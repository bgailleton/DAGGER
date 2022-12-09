//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// WIP - At the moment, direct conversion do not really work and this file is not called.
// MATLAB™®© (proprietay, license protected, not free) sometimes (?!) and for no reason encompasses 
// inbuilt std::vector structures from the standard library
// into a very obscure, undocumented, unecessary and unwanted MwCppContainerForVector<std::vector>
// This f*cks my whole system - and any system that would dare relying on calling the operator[] on a std::vector - and explains why the MATLAB™®© bindings for DAGGER are different, uncomplete and less efficient and require a specific data structure
// (yes I wrote that out of frustration after spending days on the problem. 
// For the record I never coded in Julia before and a day was enough to learn the basics about the language and seamlessly bind it to c++ and Julia is FREE OPEN SOURCE AND EVEN BEATS MATLAB performances in many benchmark tests)
// UPDATE: The creation of that MwCppContainerForVector<std::vector> happens in an automatically generated c++ file that appears a fraction of second during the building process meaning it is VIRTUALLY impossible to try and adapt a wrapper for the class
// UPDATE 2: it seems to be MATLAB version dependent too (and ofc switching from a version to another means another license)
// UPDATE 3: unlike pybind11 or Julia CXX.jl, the process is OS-dependent and require more than a few adaptations from one to another
// FINAL UPDATE: I give up, and to the question did we try to get support from mathwork? yes, just asking why the OFFICIAL tuto fails for example remains unanswered: 
// -> https://fr.mathworks.com/matlabcentral/answers/1754720-clibgen-c-on-macos-matlab2019b-cannot-use-standard-library?s_tid=answers_rc1-1_p1_BOTH
// -> https://fr.mathworks.com/matlabcentral/answers/769427-errors-parsing-header-file-when-building-interface-c-library?s_tid=answers_rc1-1_p1_BOTH
// As MATLAB™®© is a paywall and really hurts already depleted academic budgets, I'll keep the bindings to the minimum needed to interact with TopoToolBox and/or other project really needing MATLAB.
// Disclaimer: these comments  are _slightly_ subjective and opinion based. 
// Maybe with the compiler toolbox (not included in the base price ;) ) it works better? 
// but bloody hell, >1000euros a license for something that is outperformed in term of performances and ergonomy by literally all the other free-open-source alternatives...
// Now there might be some solutions or something I did wrong and some of my arguments are probably in bad faith. But again, I spent MULTIPLE FULL DAYS on it, significantly more than for python, Julia and R.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wrap_helper_MATLAB_HPP
#define wrap_helper_MATLAB_HPP


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
#include <stdlib.h>
#include <ctime>

#include "utils.hpp"

// MATLAB™®© (proprietay, license protected, not free) related includes
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

// MATLAB™®© (proprietay, license protected, not free) engine stuff
using namespace matlab::engine;



namespace DAGGER
{



template<class T>
class matlvec
{
public:
	//Create MATLAB ™®© (proprietay, license protected, not free) data array factory
	matlab::data::ArrayFactory factory;

	int isize = 0;
	size_t usize = 0;
	std::vector<T> arr;


	matlvec(){};
	matlvec(matlab::data::TypedArray<T>& arr)
	{
		this->arr = std::vector<T>(arr.begin(),arr.end());
		this->isize = arr.getNumberOfElements();
		this->usize = arr.getNumberOfElements();
	};

	matlvec(matlvec<T>& tin)
	{
		this->arr = tin.arr;
		this->isize = tin.isize;
		this->usize = tin.usize;
	}

	matlvec(const matlvec<T>& tin)
	{
		this->arr = tin.arr;
		this->isize = tin.isize;
		this->usize = tin.usize;
	}

  T & operator [](int i) {return this->arr[i];}

  void set(int i, T val)
  {this->arr[i] = val;}

  T get(int i)
  {return this->arr[i];}

  std::vector<T> to_vec()
  {
  	std::vector<T> out(this->isize);
  	for(int i=0;i<this->isize; ++i)
  		out[i] = this->arr[i];
  }

  // size_t size(){return this->usize;}
  const size_t size() {return this->usize;}
  	
};

template<class A, class T>
matlab::data::TypedArray<T> vec2Marray(A& vec)
{
	matlab::data::ArrayFactory factory;
	auto dim = vec.size();
	matlab::data::TypedArray<T> retr = factory.createArray({1,dim},vec.begin(),vec.end());
	return retr;
}

template<class A, class T>
matlab::data::TypedArray<T> pvec2Marray(A& vec)
{
	matlab::data::ArrayFactory factory;
	auto dim = vec.data->size();
	return factory.createArray({1,dim},vec.data->begin(),vec.data->end());
}


template<class T>
matlab::data::TypedArray<T> _format_output(std::vector<T>& yolo){return vec2Marray<std::vector<T>,T>(yolo);}
template<class T>
matlab::data::TypedArray<T> _format_output(pvector<T>& yolo){return pvec2Marray<pvector<T>,T>(yolo);}
template<class T>
std::vector<T> _format_output(pvector<T>& yolo){return yolo.to_vec();}

// template<class T>
// std::vector<T> _format_output(std::vector<T>& yolo){return yolo;}

template<class T>
matlab::data::TypedArray<T> _format_output(matlvec<T>& yolo){auto vec = yolo.to_vec();  return vec2Marray(vec);}
template<class T>
matlab::data::TypedArray<T> _format_output(matlab::data::TypedArray<T>& yolo){return yolo;}

// template<class T>
// std::vector<T> to_vec(matlvec<T>& in)
// {
// 	std::vector<T> out(in.arr);
// 	return out;
// }

template<class T>
std::vector<T> to_vec(matlvec<T> in)
{
	std::vector<T> out(in.arr);
	return out;
}

template<class T>
std::vector<T> to_vec(matlab::data::TypedArray<T>& tin)
{
	std::vector<T> out(tin.begin(), tin.end());
	return out;
}


template<class A>
std::vector<double> to_vec(A& tin)
{
	matlvec<double> in(tin);
	std::vector<double> out(in.size());
	for(size_t i=0;i<in.size(); ++i)
		out[i] = in[i];
	return out;
}


template<class T>
matlvec<T> _format_input(matlab::data::TypedArray<T>& yolo)
{
	matlvec<T> gabul(yolo);
	return gabul;
}

template<class T>
matlvec<T> _format_input(matlvec<T>& yolo)
{
	return yolo;
}


template<class A>
A _format_input(A& yolo)
{
	return to_vec(yolo);
}




// end of namesapce DAGGER
};








#endif