#include "D8connector.hpp"
#include "MatlabDataArray.hpp"
#include "graph.hpp"
#include "popscape.hpp"
#include "wrap_helper.hpp"
#include "wrap_helper_MATLAB.hpp"
#include <stdio.h>
#include <string>
#include <vector>

// matlab::data::Array add5(std::vector<double> foo)
// {
// 	std::vector<double> bar(foo);
// 	for(auto &v:bar)
// 		v+=5;

// 	matlab::data::ArrayFactory fac;
// 	return fac.createArray({bar.size()}, bar.begin(),bar.end());
// }

// template<class fT>
// matlab::data::TypedArray<double> quick_fluvial_topo(int ncycles, std::string
// boundaries)
matlab::data::Array quick_fluvial_topo(int ncycles, std::string boundaries) {
  auto ttopo = DAGGER::_quick_fluvial_topo<double>(ncycles, boundaries);
  return DAGGER::format_output<decltype(ttopo),
                               matlab::data::TypedArray<double>>(ttopo);
}

// fT is the generic floating point type
// int_t is the generic floating point type
template <class fT, class int_t> class daggerFD {
public:
  int_t nx, ny, nxy;
  fT dx, dy, xmin, ymin;
  std::vector<int_t> ix, ixc;
  std::vector<fT> fraction, distances;
  // DAGGER::D8connector<fT> connector;

  daggerFD() { ; }
  void init(int_t nx, int_t ny, fT dx, fT dy, fT xmin, fT ymin,
            std::string boundary_condition) {
    this->nx = nx;
    this->ny = ny;
    this->nxy = ny * nx;
    this->dx = dx;
    this->dy = dy;
    this->xmin = xmin;
    this->ymin = ymin;
    this->connector = DAGGER::D8connector<fT>(nx, ny, dx, dy, xmin, ymin);
    this->connector.set_default_boundaries(boundary_condition);
    this->graph = DAGGER::graph<fT, DAGGER::D8connector<fT>>(this->connector);
  }

  matlab::data::TypedArray<fT> compute(matlab::data::TypedArray<fT> &ttopo,
                                       bool SFD) {
    // std::vector<fT> topo = DAGGER::to_vec(ttopo);
    matlab::data::TypedArray<fT> ret =
        this->graph.template compute_graph<matlab::data::TypedArray<fT>,
                                           matlab::data::TypedArray<fT>>(
            ttopo, SFD, false);
    if (SFD) {
      this->ix = std::vector<int_t>(this->nxy, 0);
      this->ixc = std::vector<int_t>(this->nxy, 0);
      this->fraction = std::vector<fT>(this->nxy, 1.);
      this->distances = std::vector<fT>(this->connector.Sdistance2receivers);
      for (int i = 0; i < this->nxy; ++i) {
        this->ix[i] = this->graph.Sstack[i];
        this->ixc[i] = this->connector.Sreceivers[this->graph.Sstack[i]];
      }
    }

    return ret;
  }

  matlab::data::TypedArray<fT> get_DA() {
    return this->graph.template accumulate_constant_downstream_SFD<
        matlab::data::TypedArray<fT>>(this->dx * this->dy);
  }

private:
  DAGGER::D8connector<fT> connector;
  DAGGER::graph<fT, DAGGER::D8connector<fT>> graph;
};

// template<class fT>
// matlab::data::TypedArray<fT> quick_fluvial_topo(int ncycles, std::string
// boundaries)
// {
// 	auto ttopo = DAGGER::_quick_fluvial_topo<fT>(ncycles,boundaries);
// 	return DAGGER::format_output<decltype(ttopo),
// matlab::data::TypedArray<fT> >(ttopo);
// }

// // The following functions are required to simulate the instanciation of the
// different function and force the black-box obscure matlab clibgen to consider
// all the needed functions
// // THis is the only I found to do so, yes it seems dirty, but eh, template
// c++ has only be around since the late 90s (see the top of
// wrap_helper_MATLAB.hpp to understand my bitterness)

// template<class fT, class int_t>
// void _instanciator(fT XX, int_t YY)
// {
// 	int_t nx = 50, ny = 50;
// 	fT dx =1., dy = 1., xmin = 0., ymin = 0.;
// 	std::vector<fT> randtopo(nx*ny,0.);
// 	matlab::data::ArrayFactory factory;
// 	DAGGER::add_noise_to_vector(randtopo,0.,1.);
// 	unsigned long long size = nx*ny;
// 	matlab::data::TypedArray<double> ttopo = factory.createArray( {size},
// randtopo.begin(), randtopo.end()); 	daggerFD<fT, int_t> DAGDAG;
// 	DAGDAG.init(nx, ny, dx, dy, xmin, ymin, "4edges");
// 	auto dagba = DAGDAG.compute(ttopo, true);

// }

// void instanciator()
// {
// 	_instanciator<double, int>(30.,50);
// }
