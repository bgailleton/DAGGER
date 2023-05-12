#ifndef POPSCAPE_HPP
#define POPSCAPE_HPP

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
#include <omp.h>

// local includes 
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D8connector.hpp"
#include "D4connector.hpp"
#include "popscape_utils.hpp"
// #include "D8connector.hpp"
#include "graph.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// #include <pybind11/numpy.h>

// namespace py = pybind11;


namespace DAGGER
{

enum class PARAM_MODE
{
	CONSTANT,
	VARIABLE,
};



template<class float_t,class Graph_t, class Connector_t>
class popscape
{

public:

	// stored as direct children: they will be changed anyway
	Graph_t graph;
	Connector_t connector;

	std::vector<float_t> QA;
	std::vector<float_t> topography;

	PARAM_MODE Kbase_mode = PARAM_MODE::CONSTANT;
	std::vector<float_t> _Kbase = {1e-3};
	PARAM_MODE Kmod_mode = PARAM_MODE::CONSTANT;
	std::vector<float_t> _Kmod = {1.};
	PARAM_MODE m_mode = PARAM_MODE::CONSTANT;
	std::vector<float_t> _m = {0.45};
	PARAM_MODE n_mode = PARAM_MODE::CONSTANT;
	std::vector<float_t> _n = {1.};
	PARAM_MODE precip_mode = PARAM_MODE::CONSTANT;
	std::vector<float_t> _precip = {1.};

	PARAM_MODE UE_mode = PARAM_MODE::CONSTANT;
	std::vector<float_t> _UE = {1e-3};

	std::string boundary_string = "periodic_EW";



	popscape(){;}
	popscape(Graph_t& graph, Connector_t& con)
	{
		
		this->connector = _create_connector(con.nx,con.ny,con.dx,con.dy,0.,0.);
		_create_graph(graph.nnodes, this->connector,this->graph);
		
		// this->graph = graph;

	}

	template<class in_t>
	void set_topo(in_t& itopo)
	{
		auto ftopo = DAGGER::format_input(itopo);
		this->topography = to_vec(ftopo);
	}

	template<class out_t>
	out_t get_topo(){return DAGGER::format_output<std::vector<float_t>,out_t>(this->topography);}
	template<class out_t>
	out_t get_QA(){return DAGGER::format_output<std::vector<float_t>,out_t>(this->QA);}

	template<class out_t>
	out_t get_chistar(){std::vector<float_t> chistar = this->_chi_star();return DAGGER::format_output<std::vector<float_t>,out_t>(chistar);}


	// Parameters:
	float_t Kbase(int i)
	{
		if(this->Kbase_mode == PARAM_MODE::CONSTANT)
			return this->_Kbase[0];
		else
			return this->_Kbase[i];
	}

	float_t Kmod(int i)
	{
		if(this->Kmod_mode == PARAM_MODE::CONSTANT)
			return this->_Kmod[0];
		else
			return this->_Kmod[i];
	}

	float_t m(int i)
	{
		if(this->m_mode == PARAM_MODE::CONSTANT)
			return this->_m[0];
		else
			return this->_m[i];
	}

	float_t n(int i)
	{
		if(this->n_mode == PARAM_MODE::CONSTANT)
			return this->_n[0];
		else
			return this->_n[i];
	}

	// Parameters:
	float_t precip(int i)
	{
		if(this->precip_mode == PARAM_MODE::CONSTANT)
			return this->_precip[0];
		else
			return this->_precip[i];
	}
	float_t UE(int i)
	{
		if(this->UE_mode == PARAM_MODE::CONSTANT)
			return this->_UE[0];
		else
			return this->_UE[i];
	}

	void _init_vecs(){this->QA = std::vector<float_t>(this->connector.nnodes,0.);}


	void StSt(int n_iterations)
	{
		// running the code for n_iterations
		for(int nit = 0; nit < n_iterations; ++nit)
		{
			this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
			this->graph._compute_graph(this->topography,true, false);
			this->_init_vecs();
			this->QA = this->graph._accumulate_constant_downstream_SFD(this->connector.get_area_at_node(0));
			for(int i = 0; i<this->graph.nnodes; ++i)
			{
				int node = this->graph.Sstack[i];

				if(this->connector.flow_out_or_pit(node))
					continue;

				int rec = this->connector.Sreceivers[node];
				this->topography[node] = this->topography[rec] + this->connector.Sdistance2receivers[node] * 
					std::pow(this->UE(node)/(this->Kbase(node) * this->Kmod(node)),1./this->n(node)) *
					std::pow(this->QA[node],-this->m(node)/this->n(node));
			}

		}
		
	}

	void restriction(float_t noise_intensity)
	{
		
		int nxy = 4 * this->graph.nnodes;

		int nx = this->connector.nx * 2;
		int ny = this->connector.ny * 2;

		std::vector<float_t> ntopo(nxy,0.);


		for(int i=0;i<this->graph.nnodes;++i)
		{
			int row,col;
			this->connector.rowcol_from_node_id(i, row, col);
			int nid = row*2 * nx + col*2;
			ntopo[nid] = this->topography[i] + this->connector.randu->get() * noise_intensity;
			nid = row*2 * nx + col*2+1;
			ntopo[nid] = this->topography[i] + this->connector.randu->get() * noise_intensity;
			nid = (row*2+1) * nx + col*2+1;
			ntopo[nid] = this->topography[i] + this->connector.randu->get() * noise_intensity;
			nid = (row*2+1) * nx + col*2;
			ntopo[nid] = this->topography[i] + this->connector.randu->get() * noise_intensity;
		}

		
		float_t dx = this->connector.dx/2;
		float_t dy = this->connector.dy/2;
		// init connector
		this->connector = _create_connector(nx,ny,dx,dy,0.,0.);
		this->connector.set_default_boundaries(this->boundary_string);
		
		// init graph
		_create_graph(nxy, this->connector,this->graph);
		// graph.
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();
		this->border20();

	}

	// TODO
	void interpolation()
	{
		// int nxy = floor(this->graph.nnodes/4);

		int nx = floor(this->connector.nx / 2);
		int ny = floor(this->connector.ny / 2);
		int nxy = nx * ny;

		std::vector<float_t> ntopo(nxy,0.);

		D8connector<double> ncon = _create_connector(nx,ny,this->connector.dx*2,this->connector.dy*2,0.,0.);
		ncon.set_default_boundaries(this->boundary_string);



		for(int i=0;i<nxy;++i)
		{
			int row,col;
			ncon.rowcol_from_node_id(i, row, col);
			int nid = row*2 * this->connector.nx + col*2;
			float_t tntopo = 0.;
			tntopo += this->topography[nid];
			nid = row*2 * this->connector.nx + col*2+1;
			tntopo += this->topography[nid];
			nid = (row*2+1) * this->connector.nx + col*2+1;
			tntopo += this->topography[nid];
			nid = (row*2+1) * this->connector.nx + col*2;
			tntopo += this->topography[nid];

			ntopo[i] = tntopo/4;

		}
		// init connector
		this->connector = std::move(ncon);
		// init graph
		_create_graph(nxy, this->connector,this->graph);
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();
		this->border20();

	}

	void smooth(float_t rr)
	{
		// this->topography = On_gaussian_blur<float_t>(rr, this->topography, this->connector.nx, this->connector.ny);
		this->topography = On_gaussian_blur(rr, this->topography, this->connector.nx, this->connector.ny);
	}

	void border20()
	{
		for(int i=0; i<this->connector.nnodes; ++i)
		{
			if(this->connector.flow_out_model(i))
				this->topography[i] = 0;
		}
	}


	std::vector<float_t> _chi_star()
	{
		std::vector<float_t> A = this->graph._accumulate_constant_downstream_SFD(this->connector.get_area_at_node(0));
		std::vector<float_t> chistar(this->graph.nnodes, 0.);
		float_t chimax = 0;
		for(int i=0;i<this->graph.nnodes; ++i)
		{
			int node = this->graph.Sstack[i];
			int rec = this->connector.Sreceivers[node];
			if(node == rec)
				continue;

			// chistar[node] = chistar[rec] + this->connector.Sdistance2receivers[node] * (std::pow(1/A[node], this->m(node)/this->n(node)));
			chistar[node] = chistar[rec] + this->connector.Sdistance2receivers[node] * (std::pow(1/A[node], 0.2) );
			if(chistar[node] > chimax)
				chimax = chistar[node];
		}

		for(auto&v:chistar)
			v/=chimax;

		return chistar;
	}

	std::vector<float_t> _z_star()
	{
		std::vector<float_t> zstar(this->topography);
		float_t zmax = 0;
		for(int i=0;i<this->graph.nnodes; ++i)
		{
			if(zstar[i] > zmax)
				zmax = zstar[i];
		}

		for(auto&v:zstar)
			v/=zmax;

		return zstar;
	}


	void simple_Kfchi(float_t tkmod, float_t chimin, float_t chimax)
	{

		auto chistar = this->_chi_star();

		this->Kmod_mode = PARAM_MODE::VARIABLE;
		this->_Kmod = std::vector<float_t>(this->graph.nnodes, 1.);
		for(int i = 0; i<this->graph.nnodes; ++i)
		{
			if(chistar[i] > chimin && chistar[i] <chimax)
				this->_Kmod[i] = tkmod;
		}

		this->StSt(1);

		this->Kmod_mode = PARAM_MODE::CONSTANT;
		this->_Kmod = {1};

	}


	void simple_Kfz(float_t tkmod, float_t zmin, float_t zmax)
	{

		auto zstar = this->_z_star();

		this->Kmod_mode = PARAM_MODE::VARIABLE;
		this->_Kmod = std::vector<float_t>(this->graph.nnodes, 1.);
		for(int i = 0; i<this->graph.nnodes; ++i)
		{
			if(zstar[i] > zmin && zstar[i] <zmax)
				this->_Kmod[i] = tkmod;
		}

		this->StSt(1);

		this->Kmod_mode = PARAM_MODE::CONSTANT;
		this->_Kmod = {1};

	}


	// std::vector<float_t> chiculations()
	// {
		
	// }


};























// Former Popscape, to be deprecated soon
template<class float_t,class Graph_t, class Connector_t>
class popscape_old
{
public:
	
	// the graph is directly integrated into popscape
	Graph_t graph;
	// And so is the connector as popscape is optimised to pop a final landscape.
	Connector_t connector;
	
	// Topography
	std::vector<float_t> topography, QA;

	// randomiser helper
	std::shared_ptr<DAGGER::easyRand> randu = std::make_shared<DAGGER::easyRand>();


	// empty constructor
	popscape_old(){};

	// constructor 1:
	// -> noise type is from the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	popscape_old(RANDNOISE noisetype, int start_nx,  int nit, float_t dx, float_t dy)
	{

		if(start_nx % 8 != 0)
			throw std::runtime_error("target nx and start nx needs to be a multiple of 8");
		
		// total number of nodes (so far assuming only regular grids)
		int nx = start_nx, ny = start_nx;
		int nxy = nx * ny;
		
		// init the topo to 0
		this->topography = std::vector<float_t>(nxy,0.);
		
		// init connector
		this->connector = _create_connector(nx,ny,dx,dy,0.,0.);
		
		// init graph
		_create_graph(nxy, this->connector,this->graph);
		
		// init random noise
		if(noisetype == RANDNOISE::WHITE)
			DAGGER::add_noise_to_vector(this->topography,0.,1.);

		this->connector.set_default_boundaries("periodic_EW");

		for(int i = 0; i<nit; ++i)
		{	
			// std::cout << "REFINING::" << i << std::endl;
			int nnit = (i==0)?50:5;
			// std::cout << "REFINING::A" << i << std::endl;
			this->solve_generic(0.45, 1.11, nnit);
			// std::cout << "REFINING::B" << i << std::endl;
			if(i<nit-1)
			{
				// std::cout << "REFINING::C" << i << std::endl;

				this->double_res(10, noisetype);
			}
			// std::cout << "done" << std::endl;
		}


	}

	void _init_vecs(){this->QA = std::vector<float_t>(this->graph.nnodes,0.);}


	void solve_generic(float_t m, float_t n, int n_iterations)
	{
		// std::cout << "REFINING::|||||" <<n_iterations << std::endl;
		// running the code for n_iterations
		for(int nit = 0; nit < n_iterations; ++nit)
		{
			// std::cout << nit << "|A" << std::endl;

			this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
			this->graph._compute_graph(this->topography,true, false);
			this->_init_vecs();

			this->QA = this->graph._accumulate_constant_downstream_SFD(this->connector.get_area_at_node(0));

			// std::cout << nit << "|C" << std::endl;

			for(int i = 0; i<this->graph.nnodes; ++i)
			{
				int node = this->graph.Sstack[i];
				if(this->connector.flow_out_or_pit(node))
					continue;
				int rec = this->connector.Sreceivers[node];
				this->topography[node] = this->topography[rec] + this->connector.Sdistance2receivers[node] * std::pow(1e2,1./n)/std::pow(this->QA[node],m/n);
				// this->topography[node] = this->topography[rec] + 1;
			}
			// std::cout << nit << "|D" << std::endl;

		}
			// std::cout << std::endl;
		
	}

	void double_res(float_t noise_intensity, RANDNOISE noisetype)
	{
		
		int nxy = 4 * this->graph.nnodes;

		int nx = this->connector.nx * 2;
		int ny = this->connector.ny * 2;

		std::vector<float_t> ntopo(nxy,0.);


		if(noisetype == RANDNOISE::WHITE)
		for(int i=0;i<this->graph.nnodes;++i)
		{
			int row,col;
			this->connector.rowcol_from_node_id(i, row, col);
			int nid = row*2 * nx + col*2;
			ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
			nid = row*2 * nx + col*2+1;
			ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
			nid = (row*2+1) * nx + col*2+1;
			ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
			nid = (row*2+1) * nx + col*2;
			ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
		}

		
		float_t dx = this->connector.dx/2;
		float_t dy = this->connector.dy/2;
		// init connector
		this->connector = _create_connector(nx,ny,dx,dy,0.,0.);
		
		// init graph
		_create_graph(nxy, this->connector,this->graph);
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();

	}

	void compute_graph(bool SFD_only)
	{
		this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		this->graph._compute_graph(this->topography, SFD_only, false);
	}

	void compute_DA_SFD(){this->QA = this->graph._accumulate_constant_downstream_SFD(this->connector.get_area_at_node(0));}

	void solve_SFD_SPL_imp(float_t m, float_t n, float_t K, float_t dt)
	{

    for (int i=0; i< this->graph.nnodes; ++i)
  	{

	    int node = this->graph.Sstack[i];
	    int rec = this->connector.Sreceivers[node];


	    if (!this->connector.flow_out_or_pit(node) == false)
				continue;

	    float_t factor = K * dt * std::pow(this->QA[node],m) / std::pow(this->connector.Sdistance2receivers[node],n);

	    float_t ielevation = this->topography[node];
	    float_t irec_elevation = this->topography[rec];

	    float_t elevation_k = ielevation;
	    float_t elevation_prev = std::numeric_limits<float_t>::max();
	    float_t tolerance = 1e-4;

	    while (abs(elevation_k - elevation_prev) > tolerance)
	    {
				elevation_prev = elevation_k;
				float_t slope = std::max(elevation_k - irec_elevation,1e-6);
				float_t diff =  (elevation_k - ielevation + factor * std::pow(slope,n)) / ( 1. + factor * n * std::pow(slope, n - 1) )  ;
				elevation_k -= diff;
	    }

	    this->topography[node] = elevation_k;
	  }
	}


	template<class out_t>
	out_t get_topo(){return DAGGER::format_output<std::vector<float_t>,out_t>(this->topography);}
	template<class out_t>
	out_t get_QA(){return DAGGER::format_output<std::vector<float_t>,out_t>(this->QA);}

	void apply_uplift(float_t dt, float_t U)
	{
		for(int i=0; i < this->graph.nnodes;++i)
		{
			if(!this->connector.flow_out_or_pit(i))
			  this->topography[i]+=U*dt;
		}
	}

	template<class out_t>
	void apply_variable_uplift(float_t dt, out_t tU)
	{
		auto U = DAGGER::format_input(tU);
		for(int i=0; i < this->graph.nnodes;++i)
		{
			if(!this->connector.flow_out_or_pit(i))
			  this->topography[i]+=U[i]*dt;
		}
	}


	void hydraulic_erosion_v0(int n_particules, float_t erosor)
	{

		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator
		std::uniform_int_distribution<> distr(0, this->graph.nnodes-1); // define the range

		std::vector<float_t> newtopo(this->topography);

		// for(int tnp=0; tnp<n_particules; ++tnp)
		// {
		// 	// spawning the particle
		// 	Particle tpart(distr(gen));
		// 	// std::cout << tpart.pos << std::endl;
		// 	std::vector<int> nelinks(8,0);
		// 	while(true)
		// 	{
		// 		// Increment the speed
		// 		int NnN = this->connector.get_neighbour_idx_links(tpart.pos, nelinks);
		// 		int nrecs = 0;
		// 		for(int i =0; i<NnN ; ++i)
		// 		{
		// 			// std::cout << "A" << std::endl;
		// 			int li = nelinks[i];
		// 			int n1 = this->graph.linknodes[li*2];
		// 			int n2 = this->graph.linknodes[li*2 + 1];
		// 			// std::cout << "B" << std::endl;
		// 			if((n1 == tpart.pos && this->topography[n1] <= this->topography[n2]) || (n2 == tpart.pos && this->topography[n2] < this->topography[n1]))
		// 				continue;
		// 			++nrecs;
		// 			auto dxdy = this->connector.get_directed_dxdy_from_links_idx( li, tpart.pos, n1, n2);
		// 			auto grad = std::abs(this->topography[n1] + this->topography[n2])/this->connector.get_dx_from_links_idx(li);
		// 			// std::cout << "C" << std::endl;
		// 			tpart.speed_up(dxdy,grad);
		// 		}
		// 		// int col,row; this->connector.rowcol_from_node_id(tpart.pos, row, col);
		// 		auto dirxy = tpart.get_normalised_speed();
		// 		int newpos = this->connector.get_neighbour_idx_from_normalised_dxdy(tpart.pos,dirxy.first, dirxy.second);
		// 		// std::cout << newpos << "|" << dirxy.first << "|"  << dirxy.second << std::endl;
				
		// 		// temp erosion
		// 		newtopo[tpart.pos] -= erosor;
		// 		if(!this->graph.flow_out_or_pit(tpart.pos,this->connector) == false || nrecs == 0 || newpos < 0 || newpos > this->graph.nnodes - 1 )
		// 			break;

		// 		tpart.pos = newpos;
		// 	}
		// }

		this->topography = std::move(newtopo);

	}


	// void precipitSPL(float_t dt, float_t K, float_t m, float_t n, int n_prec)
	// {

	// 	float_t Aprec = std


	// }


	void normalise_topography()
	{
		float_t mini = std::numeric_limits<float_t>::max();
		float_t maxi = std::numeric_limits<float_t>::min();
		for(auto v:this->topography)
		{
			if(v < mini) mini = v;
			if(v > maxi) maxi = v;
		}
		maxi -= mini;
		for(auto& v:this->topography)
			v = (v - mini)/maxi;
	}


	void run_SFD_exp_latmag(float_t K, float_t m, float_t n, float_t Kl, float_t dt)
	{
		// this->compute_graph(true);
		// this->compute_DA_SFD();
		// std::vector<float_t> vmot(this->graph.nnodes,0);

		// for(int i=this->graph.nnodes-1; i>=0; --i)
		// {
		// 	int node = this->graph.Sstack[i];
		// 	if(!this->graph.flow_out_or_pit(node, this->connector) == false) continue;

		// 	int rec = this->graph.Sreceivers[node];

		// 	float_t S = std::max((this->topography[node] - this->topography[rec])/this->graph.Sdistance2receivers[node], 1e-6);

		// 	float_t E = std::pow(S,n) * std::pow(this->QA[node],m) * K;

		// 	auto orthonodes = this->connector.get_orthogonal_nodes(node,rec);

		// 	if(this->topography[orthonodes.first] - this->topography[node] > 0) vmot[orthonodes.first] -= Kl * E * dt;
		// 	if(this->topography[orthonodes.second] - this->topography[node] > 0) vmot[orthonodes.second] -= Kl * E * dt;

		// 	vmot[node] -= E * dt;
		// }

		// for(int i=0; i<this->graph.nnodes; ++i)
		// 	this->topography[i] += vmot[i];
	}

	
	


};

















// End of namespace fastflood
};





























#endif