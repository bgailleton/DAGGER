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

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


namespace DAGGER
{

template<class float_t,class Graph_t, class Connector_t>
class popscape
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
	popscape(){};

	// constructor 1:
	// -> noise type is from the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	popscape(RANDNOISE noisetype, int start_nx,  int nit, float_t dx, float_t dy)
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
		this->graph.init_graph(this->connector);
		
		// init random noise
		if(noisetype == RANDNOISE::WHITE)
			DAGGER::add_noise_to_vector(this->topography,0.,1.);

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
			this->graph._compute_graph(this->topography, this->connector,true, false);
			this->_init_vecs();
			// std::cout << nit << "|B" << std::endl;

			this->QA = this->graph._accumulate_constant_downstream_SFD(this->connector,this->connector.get_area_at_node(0));

			// std::cout << nit << "|C" << std::endl;

			for(int i = 0; i<this->graph.nnodes; ++i)
			{
				int node = this->graph.Sstack[i];
				if(!this->graph.flow_out_or_pit(node,this->connector) == false)
					continue;
				int rec = this->graph.Sreceivers[node];
				this->topography[node] = this->topography[rec] + this->graph.Sdistance2receivers[node] * std::pow(1e2,1./n)/std::pow(this->QA[node],m/n);
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
		this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();

	}

	void compute_graph(bool SFD_only)
	{
		this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		this->graph._compute_graph(this->topography, this->connector, SFD_only, false);
	}

	void compute_DA_SFD(){this->QA = this->graph._accumulate_constant_downstream_SFD(this->connector,this->connector.get_area_at_node(0));}

	void solve_SFD_SPL_imp(float_t m, float_t n, float_t K, float_t dt)
	{

    for (int i=0; i< this->graph.nnodes; ++i)
  	{

	    int node = this->graph.Sstack[i];
	    int rec = this->graph.Sreceivers[node];


	    if (!this->graph.flow_out_or_pit(node,this->connector) == false)
       continue;

	    float_t factor = K * dt * std::pow(this->QA[node],m) / std::pow(this->graph.Sdistance2receivers[node],n);

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
			if(!this->graph.flow_out_or_pit(i,this->connector))
			  this->topography[i]+=U*dt;
		}
	}

	template<class out_t>
	void apply_variable_uplift(float_t dt, out_t tU)
	{
		auto U = DAGGER::format_input(tU);
		for(int i=0; i < this->graph.nnodes;++i)
		{
			if(!this->graph.flow_out_or_pit(i,this->connector))
			  this->topography[i]+=U[i]*dt;
		}
	}


	void hydraulic_erosion_v0(int n_particules, float_t erosor)
	{

		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator
		std::uniform_int_distribution<> distr(0, this->graph.nnodes-1); // define the range

		std::vector<float_t> newtopo(this->topography);

		for(int tnp=0; tnp<n_particules; ++tnp)
		{
			// spawning the particle
			Particle tpart(distr(gen));
			// std::cout << tpart.pos << std::endl;
			std::vector<int> nelinks(8,0);
			while(true)
			{
				// Increment the speed
				int NnN = this->connector.get_neighbour_idx_links(tpart.pos, nelinks);
				int nrecs = 0;
				for(int i =0; i<NnN ; ++i)
				{
					// std::cout << "A" << std::endl;
					int li = nelinks[i];
					int n1 = this->graph.linknodes[li*2];
					int n2 = this->graph.linknodes[li*2 + 1];
					// std::cout << "B" << std::endl;
					if((n1 == tpart.pos && this->topography[n1] <= this->topography[n2]) || (n2 == tpart.pos && this->topography[n2] < this->topography[n1]))
						continue;
					++nrecs;
					auto dxdy = this->connector.get_directed_dxdy_from_links_idx( li, tpart.pos, n1, n2);
					auto grad = std::abs(this->topography[n1] + this->topography[n2])/this->connector.get_dx_from_links_idx(li);
					// std::cout << "C" << std::endl;
					tpart.speed_up(dxdy,grad);
				}
				// int col,row; this->connector.rowcol_from_node_id(tpart.pos, row, col);
				auto dirxy = tpart.get_normalised_speed();
				int newpos = this->connector.get_neighbour_idx_from_normalised_dxdy(tpart.pos,dirxy.first, dirxy.second);
				// std::cout << newpos << "|" << dirxy.first << "|"  << dirxy.second << std::endl;
				
				// temp erosion
				newtopo[tpart.pos] -= erosor;
				if(!this->graph.flow_out_or_pit(tpart.pos,this->connector) == false || nrecs == 0 || newpos < 0 || newpos > this->graph.nnodes - 1 )
					break;

				tpart.pos = newpos;
			}
		}

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
		this->compute_graph(true);
		this->compute_DA_SFD();
		std::vector<float_t> vmot(this->graph.nnodes,0);

		for(int i=this->graph.nnodes-1; i>=0; --i)
		{
			int node = this->graph.Sstack[i];
			if(!this->graph.flow_out_or_pit(node, this->connector) == false) continue;

			int rec = this->graph.Sreceivers[node];

			float_t S = std::max((this->topography[node] - this->topography[rec])/this->graph.Sdistance2receivers[node], 1e-6);

			float_t E = std::pow(S,n) * std::pow(this->QA[node],m) * K;

			auto orthonodes = this->connector.get_orthogonal_nodes(node,rec);

			if(this->topography[orthonodes.first] - this->topography[node] > 0) vmot[orthonodes.first] -= Kl * E * dt;
			if(this->topography[orthonodes.second] - this->topography[node] > 0) vmot[orthonodes.second] -= Kl * E * dt;

			vmot[node] -= E * dt;
		}

		for(int i=0; i<this->graph.nnodes; ++i)
			this->topography[i] += vmot[i];
	}

	
	


};

















// End of namespace fastflood
};





























#endif