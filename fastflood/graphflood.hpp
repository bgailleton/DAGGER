#ifndef GRAPHFLOOD_HPP
#define GRAPHFLOOD_HPP

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

// local includes 
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D8connector.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"
#include "fastflood_recorder.hpp"

namespace DAGGER
{

enum class HYDRO // : std::uint8_t
{
	GRAPH_SFD,
	GRAPH_MFD,
};


enum class MORPHO // : std::uint8_t
{
	NONE,
	TL,
};

enum class PARAM_KE // : std::uint8_t
{
	CONSTANT,
	VARIABLE,
	EROSION,
};

enum class PARAM_DT_HYDRO // : std::uint8_t
{
	CONSTANT,
	VARIABLE,
	COURANT,
};

enum class PARAM_DT_MORPHO // : std::uint8_t
{
	CONSTANT,
	VARIABLE,
	COURANT,
	HYDRO,
};


enum class HYDROGRAPH_LM
{
	IGNORE,
	FILL,
};


enum class MFD_PARTITIONNING
{
	PROPOSLOPE,
	SQRTSLOPE,
	PROPOREC,
};

template<class float_t, class Graph_t, class Connector_t>
class graphflood
{
public:

	// Underlying grid:
	std::shared_ptr<Graph_t> graph;
	std::shared_ptr<Connector_t> connector;

	// Global modes
	HYDRO hydromode = HYDRO::GRAPH_MFD;
	MORPHO morphomode = MORPHO::NONE;
	DEPRES depression_resolver = DEPRES::priority_flood;
	HYDROGRAPH_LM depression_management = HYDROGRAPH_LM::FILL;
	MFD_PARTITIONNING weight_management = MFD_PARTITIONNING::PROPOSLOPE;

	// Global constants:
	const float_t GRAVITY = 9.81, FIVETHIRD = 5./3., minslope = 1e-6;

	bool stochaslope = false;
	float_t stochaslope_coeff = 1.;

	bool hydrostationary = true;

	// vecotr of node size
	// #  Hydraulic Surface elevation (bedrock + sed + water)
	std::vector<float_t> _surface;

	// # Water depth
	std::vector<float_t> _hw;
	// # Water discharge
	std::vector<float_t> _Qw;

	// # Sediment discahrge
	std::vector<float_t> _Qs;
	// # Sediment height
	std::vector<float_t> _hs;


	// ###################################### 
	// Parameters for morpho ################
	// ###################################### 

	// # a exponent for erosion
	bool mode_aexp = false;
	std::vector<float_t> _aexp;

	// # ke (coefficient for erosion)
	PARAM_KE mode_ke = PARAM_KE::CONSTANT;
	std::vector<float_t> _ke = {1e-4};

	// # ke (coefficient for erosion)
	bool mode_ke_lateral = false;
	std::vector<float_t> _ke_lateral = {0.1};

	// # ke (coefficient for erosion)
	bool mode_kd = false;
	std::vector<float_t> _kd = {100};;

	// # kd (coefficient for erosion)
	bool mode_kd_lateral = false;
	std::vector<float_t> _kd_lateral = {0.1};;
	
	// # ke (coefficient for erosion)
	PARAM_DT_MORPHO mode_dt_morho = PARAM_DT_MORPHO::HYDRO;
	std::vector<float_t> _dt_morho = {1e-3};

	// Randomiser helper
	DAGGER::easyRand randu;

	// ###################################### 
	// Parameters for Hydrograph ############
	// ######################################

	// # ke (coefficient for erosion)
	bool mode_mannings = false;
	std::vector<float_t> _mannings = {0.033};

	// float_t mannings = 0.033
	float_t topological_number = 4./8;

	// # ke (coefficient for erosion)
	PARAM_DT_HYDRO mode_dt_hydro = PARAM_DT_HYDRO::CONSTANT;
	std::vector<float_t> _dt_hydro = {1e-3};


	bool hflow = false;



	graphflood(){};
	graphflood(Graph_t& graph, Connector_t& connector)
	{
		// Ingesting graph and connectors
		this->graph = std::make_shared<Graph_t>(graph);
		this->connector = std::make_shared<Connector_t>(connector);
	}


	float_t hw(int i){return this->_hw[i];}
	float_t Qw(int i){return this->_Qw[i];}
	float_t Qs(int i){return this->_Qs[i];}
	float_t surface(int i){return this->_surface[i];}


	float_t aexp(int i)
	{
		if(this->mode_aexp)
			return this->_aexp[i];
		else
			return this->_aexp[0];
	}

	float_t ke(int i)
	{
		if(PARAM_KE::VARIABLE == this->mode_ke)
			return this->_ke[i];
		else
			return this->_ke[0];
	}

	float_t ke_lateral(int i)
	{
		if(this->mode_ke_lateral)
			return this->_ke_lateral[i];
		else
			return this->_ke_lateral[0];
	}

	float_t kd(int i)
	{
		if (this->mode_kd)
			return this->_kd[i];
		else
			return this->_kd[0];
	}

	float_t kd_lateral(int i)
	{
		if (this->mode_kd_lateral)
			return this->_kd_lateral[i];
		else
			return this->_kd_lateral[0];
	}

	float_t dt_hydro(int i)
	{
		if (this->mode_dt_hydro == PARAM_DT_HYDRO::VARIABLE)
			return this->_dt_hydro[i];
		else
			return this->_dt_hydro[0];
	}

	float_t dt_morpho(int i)
	{
		if (this->mode_dt_morpho == PARAM_DT_MORPHO::VARIABLE)
			return this->_dt_morpho[i];
		else if (this->mode_dt_morpho == PARAM_DT_MORPHO::HYDRO)
			return this->_dt_hydro(i);
		else
			return this->_dt_morpho[0];
	}


	float_t mannings(int i)
	{
		if (this->mode_kd)
			return this->_kd[i];
		else
			return this->_kd[0];
	}


	// Main running function:
	void run()
	{
		// Graph Processing
		this->graph_automator();

		std::vector<float_t> vmot, vmot_hw(this->graph->nnodes,0.);
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<float_t>(this->graph->nnodes,0.);

		// main loop
		std::vector<int> receivers = this->connector->get_empty_neighbour();
		std::vector<float_t> weights(receivers.size(),0.), slopes(receivers.size(),0.);

		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			// getting next node in line
			int node = this->get_istack_node(i);
			// Getting the receivers
			int nrecs = (SF) ? 1 : this->connector->get_receivers_idx_links(node,receivers);
			float_t Smax;
			if(SF == false)
			{
				Smax = this->weights_automator(receivers, weights, slopes, node, nrecs);
			}
			else
			{
				Smax = this->get_Sw(node,this->connector->receivers[node],this->connector->Sdistance2receivers[node],this->minslope);
			}


			// Initialising the total Qout
			float_t total_Qout = 0.;
			float_t Qwin = this->_Qw[node];
			// precalculating the power
			float_t pohw = std::pow(this->hw[node], this->FIVETHIRD);
			// Squarerooting Smax
			Smax = std::sqrt(Smax);
			// going through the receiver(s)
			for(int j =0; j<nrecs; ++j)
			{
				// Hydro for the link
				// universal params
				int rec = (SF)?this->connector->Sreceivers[node]:this->connector->get_to_links(receivers[j]);
				float_t dx = (SF)?this->connector->Sdistance2receivers[node]:this->connector->get_dx_from_links_idx(receivers[j]);
				// float_t dl = this->connector->get_travers_dy_from_dx(dx);
				float_t Sw = (SF)?this->get_Sw(node,rec,dx,this->minslope):slopes[j];
				float_t tQout = (SF) ? dx * this->mannings(node) * pohw * Smax: this->topological_number * pohw / dx * Sw/Smax * this->connector->get_area_at_node(node)/this->mannings(node);
				total_Qout += tQout;
				if(this->hydrostationary)
					this->_Qw[rec] += (SF)?Qwin: weights[i] * Qwin; 
				else
					this->_Qw[rec] += tQout;

				// Now the morpho
				// if(this->morphomode != MORPHO::NONE)
				// {
				// 	float_t tau = 
				// }				
			}

			vmot_hw[node] += (this->_Qw[node] - total_Qout)/this->connector->get_area_at_node(node);
		}


		// Applying vmots
		for(int i=0; i<this->graph->nnodes; ++i)
		{
			this->hw[i] += vmot_hw[i] * this->dt_hydro(i);
			this->surface[i] += vmot_hw[i] * this->dt_hydro(i);
			if(this->morphomode != MORPHO::NONE)
				this->surface[i] += vmot[i] * this->dt_morpho(i);
		}

	}


	void graph_automator()
	{

		// is SS?
		bool only_SD = (this->hydromode == HYDRO::GRAPH_MFD) ? false: true;
		// making sure it has the right depression solver (SHOULD BE MOVED TO THE GRAPH MANAGEMENT LATER)
		this->graph->set_LMR_method(this->depression_resolver);
		// preformatting post-topo
		std::vector<float_t> post_topo(this->surface);
		this->graph->_compute_graph(post_topo, only_SD, false);

		// fill water where depressions have been solved
		if(this->depression_management == HYDROGRAPH_LM::FILL)
		{
			for(int i=0; i<this->graph->nnodes; ++i)
			{
				if(this->surface[i] < post_topo[i])
				{
					this->hw[i] += post_topo[i] - this->surface[i];
					this->surface[i] = post_topo[i];
				}
			}
		}
	}

	float_t weights_automator(std::vector<int>& receivers, std::vector<float_t>& weights, std::vector<float_t>& slopes, int& node, int& nrecs)
	{

		float_t sumw = 0., Smax = this->minslope;
		for(int i = 0; i < nrecs; ++i)
		{
			int lix = receivers[i];
			// int rec = this->connector->get_to_links(lix);
			slopes[i] = this->get_Sw(lix, this->minslope);
			if(slopes[i]>Smax)
				Smax = slopes[i];
			// slopes[i] = slope;
			if(this->weight_management == MFD_PARTITIONNING::PROPOSLOPE)
				weights[i] = slopes[i];
			else if(this->weight_management == MFD_PARTITIONNING::SQRTSLOPE)
				weights[i] = std::sqrt(slopes[i]);
			else if(this->weight_management == MFD_PARTITIONNING::PROPOREC)
				weights[i] = 1.;

			if(this->stochaslope)
				weights[i] *= this->stochaslope_coeff * this->randu.get();

			sumw += weights[i];
		}

		for(auto&v:weights)
			v = v/sumw;

		return Smax;
	}

	// small helper function returning node index from the right stack
	int get_istack_node(int i){return (this->hydromode == HYDRO::GRAPH_MFD)?this->graph->stack[i]:this->graph->Sstack[i];}

	float_t get_Sw(int lix, float_t minslope)
	{
		auto no = this->connector->get_from_to_links(lix);
		return std::min((this->surface[no.first] - this->surface[no.second])/this->connector->get_dx_from_links_idx(),minslope); 
	}

	float_t get_Sw(int node, int rec, float_t dx, float_t minslope)
	{
		return std::min((this->surface[node] - this->surface[rec])/dx,minslope); 
	}



	// GET DATA OUT
	template<class out_t>
	out_t get_hw(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_hw) ;}
	template<class out_t>
	out_t get_surface_topo(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->surface) ;}






};
// end of graphflood class






} // End of DAGGER namespace

































#endif