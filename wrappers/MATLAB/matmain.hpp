#include "D8connector.hpp"
#include "D4connector.hpp"
#include "graph.hpp"
#include <vector>
#include <string>
#include <stdio.h>
#include "wrap_helper.hpp"
#include "wrap_helper_MATLAB.hpp"
#include "MatlabDataArray.hpp"


// DAGGER::D8connector<double> createConnectorD8(int nx, int ny, double dx, double dy, double xmin, double ymin){return DAGGER::D8connector<double>(nx,ny,dx,dy,xmin,ymin);}
// DAGGER::graph<double> createGraph(int nnodes, int n_neighbours){return DAGGER::graph<double>(nnodes,n_neighbours);}
// std::string dummytest()
// {

// 	DAGGER::easyRand ttbt;
// 	std::vector<double> yolo(2500,0.);
// 	DAGGER::D8connector<double> d8(50,50,1,1,0.,0.);
// 	DAGGER::graph<double>tg(2500,8);
// 	for(int i=0; i<2500; ++i)
// 		yolo[i]+= ttbt.get();;
// 	auto fill = tg.compute_graph<DAGGER::D8connector<double>, std::vector<double>, std::vector<double>>(yolo, d8, true, true);

// 	return "OK";
// };

matlab::data::Array add5(std::vector<double> foo)
{
	std::vector<double> bar(foo);
	for(auto &v:bar)
		v+=5;

	matlab::data::ArrayFactory fac;
	return fac.createArray({bar.size()}, bar.begin(),bar.end());
}




// float_t is the generic floating point type
// int_t is the generic floating point type
template<class float_t, class int_t>
class daggerFD
{
public:

	int_t nx, ny, nxy;
	float_t dx, dy, xmin, ymin;
	std::vector<int_t> ix,ixc;
	std::vector<float_t> fraction, distances;


	daggerFD(){;}
	void init(int_t nx, int_t ny, float_t dx, float_t dy, float_t xmin, float_t ymin, std::string boundary_condition)
	{
		this->nx = nx;
		this->ny = ny;
		this->nxy = ny*nx;
		this->dx = dx;
		this->dy = dy;
		this->xmin = xmin;
		this->ymin = ymin;
		this->connector = DAGGER::D8connector<float_t>(nx,ny,dx,dy,xmin,ymin);
		this->connector.set_default_boundaries(boundary_condition);
		this->graph = DAGGER::graph<float_t>(nx*ny,8);
		this->graph.init_graph(this->connector);
	}

	std::vector<float_t> compute(matlab::data::TypedArray<float_t>& ttopo, bool SFD)
	{
		std::vector<float_t> topo = DAGGER::to_vec(ttopo);
		std::vector<float_t> ret = this->graph.template compute_graph< DAGGER::D8connector<float_t> , std::vector<float_t>, std::vector<float_t> >(topo, this->connector, SFD, true) ;
		if(SFD)
		{
			this->ix = std::vector<int_t>(this->nxy, 0);
			this->ixc = std::vector<int_t>(this->nxy, 0);
			this->fraction = std::vector<float_t>(this->nxy, 1.);
			this->distances = std::vector<float_t>(this->graph.Sdistance2receivers);
			for(int i=0;i<this->nxy;++i)
			{
				this->ix[i] = this->graph.Sstack[i];
				this->ixc[i] = this->graph.Sreceivers[this->graph.Sstack[i]];
			}
		}

		return DAGGER::to_vec(ret);
	}

	std::vector<float_t> get_DA(){std::vector<float_t> out = this->graph.template accumulate_constant_downstream_SFD< DAGGER::D8connector<float_t>, std::vector<float_t> >(this->connector,this->dx * this->dy);return out;}




private:
	
	DAGGER::D8connector<float_t> connector;
	DAGGER::graph<float_t> graph;




};



// The following functions are required to simulate the instanciation of the different function and force the black-box obscure matlab clibgen to consider all the needed functions
// THis is the only I found to do so, yes it seems dirty, but eh, template c++ has only be around since the late 90s (see the top of wrap_helper_MATLAB.hpp to understand my bitterness)

template<class float_t, class int_t>
void _instanciator(float_t XX, int_t YY)
{
	int_t nx = 50, ny = 50;
	float_t dx =1., dy = 1., xmin = 0., ymin = 0.;
	std::vector<float_t> randtopo(nx*ny,0.);
	matlab::data::ArrayFactory factory;
	DAGGER::add_noise_to_vector(randtopo,0.,1.);
	unsigned long long size = nx*ny;
	matlab::data::TypedArray<double> ttopo = factory.createArray( {size}, randtopo.begin(), randtopo.end());
	daggerFD<float_t, int_t> DAGDAG;
	DAGDAG.init(nx, ny, dx, dy, xmin, ymin, "4edges");
	auto dagba = DAGDAG.compute(ttopo, true);

}

void instanciator()
{
	_instanciator<double, int>(30.,50);
}


