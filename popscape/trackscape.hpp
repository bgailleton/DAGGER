//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Trackscape main code
// 
// Trackscape is a static CHONK implementation using DAGGER as backend
// I tuned it for tracking provenance and sediment timing  at relatively long term.
// Fluvial and hillslope processes are using laws adapted from CIDRE (Carretier et al. 2016)
//
// B. Gailleton - 2022
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


#ifndef TRACKSCAPE_HPP
#define TRACKSCAPE_HPP

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
#include "trackscape_utils.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

// #include "sedcells.hpp"



namespace DAGGER
{


// templates are: floating point type, graph type (DAGGER graph) and connector type (DAGGER connector)
template<class float_t,class Graph_t, class Connector_t>
class trackscape
{
public:

	// Surface topography (top of bedrock + sediment)
	std::vector<float_t> z_surf;
	// Sediment height
	std::vector<float_t> h_sed;
	// Water Discharge (as P * drainage area) 
	std::vector<float_t> Qw;
	// Hillslope sediment flux
	std::vector<float_t> Qs_hs;
	// Fluvial sediment flux
	std::vector<float_t> Qs_sed;
	
	// Spatial erodibility (sediments)
	std::vector<float_t> _Ks;
	// Spatial erodibility (bedrock)
	std::vector<float_t> _Kr;
	// Spatial diffusivity (sediments)
	std::vector<float_t> _kappa_s;
	// Spatial diffusivity (bedrock)
	std::vector<float_t> _kappa_r;
	// Spatial precipitations
	std::vector<float_t> _precipitations;
	// Spatial deposition coefficient
	std::vector<float_t> _depcoeff;
	// Spatial critical slope
	std::vector<float_t> _Sc;
	//Spatial marine _Sc
	std::vector<float_t> _Sc_M;
	//Spatial marine Ke
	std::vector<float_t> _Ke;
	//Spatial marine lamdba
	std::vector<float_t> _lambda;

	// SPL exponents
	float_t mexp = 0.45,nexp = 1;

	// Switching on and/or off the hillslope and/or fluvial 
	bool hillslopes = true, fluvial = true, marine = false;

	// Switching on and/or off the spatiality of the different parameters
	bool variable_Ks = false;
	bool variable_Kr = false;
	bool variable_kappa_s = false;
	bool variable_kappa_r = false;
	bool variable_precipitations = false;
	bool variable_depcoeff = false;
	bool variable_Sc = false;
	bool variable_Sc_M = false;
	bool variable_Ke = false;
	bool variable_lambda = false;


	// Module TSP: 
	// TRACKING SINGLE SOURCE
	// # IS it activated
	bool TSP_module = false;
	// # source concentration (bedrock)
	std::vector<float_t> TSP_concentrations;
	// # tracking the suspended concentrations
	std::vector<BasePropStorer<float_t> > TSP_Qsf, TSP_Qsh;
	// # Stratigraphic info
	VerticalStorer<float_t, BasePropStorer<float_t> > TSP_store;


	// Module Ch_MTSI
	// CHRONO MEAN TIME SINCE INCISION
	// # Activated?? 
	bool Ch_MTSI = false;
	// # tracking the mobile timer (hs vs fl) 
	std::vector<BasePropStorer<float_t> > Ch_MTSI_Qsf, Ch_MTSI_Qsh;
	// # tracking the strati timer (hs vs fl) 
	VerticalStorer<float_t, BasePropStorer<float_t> > Ch_MTSI_store;

	// The graph (DAGGER)
	Graph_t graph;

	// The connector (DAGGER)
	Connector_t connector;

	// empty constructor with default values for the different param (default = non spatial)
	trackscape()
	{
		this->_Ks = {2e-5};
		this->_Kr = {1e-5};
		this->_kappa_s = {2e-2};
		this->_kappa_r = {1e-2};
		this->_precipitations = {1.};
		this->_depcoeff = {10.};
		this->_Sc = {0.6};
		this->_Sc_M = {0.6};
		this->_Ke = {2e-2};
		this->_lambda = {1000};

	};

	// constructor 1:
	// -> noise type is from the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	void init_random(RANDNOISE noisetype, int nx, int ny, float_t dx, float_t dy, std::string boundary_conditions)
	{
		
		// total number of nodes (so far assuming only regular grids)
		int nxy = nx * ny;
		
		// init the topo to 0
		this->z_surf = std::vector<float_t>(nxy,0.);
		
		// init connector
		this->connector = _create_connector(nx,ny,dx,dy,0.,0.);

		// boundary conditions:
		this->connector.set_default_boundaries(boundary_conditions);
		
		// init graph
		_create_graph(nxy, this->connector,this->graph);
		this->graph.init_graph(this->connector);
		
		// init random noise
		if(noisetype == RANDNOISE::WHITE)
			DAGGER::add_noise_to_vector(this->z_surf,0.,1.);

		this->h_sed = std::vector<float_t>(this->graph.nnodes, 0.);
	}




	// switch on and off the hillslopes and fluvial processes
	void hillslopes_on(){this->hillslopes = true;}
	void hillslopes_off(){this->hillslopes = false;}
	void fluvial_on(){this->fluvial = true;}
	void fluvial_off(){this->fluvial = false;}
	void marine_on(){this->marine = true;}
	void marine_off(){this->marine = false;}

	// Set the single parameters
	void set_single_Ks(float_t tks){this->_Ks = {tks}; this->variable_Ks = false;}
	void set_single_Kr(float_t tkr){this->_Kr = {tkr}; this->variable_Kr = false;}
	void set_single_depcoeff(float_t tdepcoeff){this->_depcoeff = {tdepcoeff}; this->variable_depcoeff = false;}
	void set_single_precipitations(float_t tprecipitations){this->_precipitations = {tprecipitations}; this->variable_precipitations = false;}
	void set_single_kappa_s(float_t tkappa_s){this->_kappa_s = {tkappa_s}; this->variable_kappa_s = false;}
	void set_single_kappa_r(float_t tkappa_r){this->_kappa_r = {tkappa_r}; this->variable_kappa_r = false;}
	void set_single_Sc_M(float_t tSc_M){this->_Sc_M = {tSc_M}; this->variable_Sc_M = false;}
	void set_single_Sc(float_t tSc){this->_Sc = {tSc}; this->variable_Sc = false;}
	void set_single_Ke(float_t TKe){this->_Ke = {TKe}; this->variable_Ke = false;}
	void set_single_lambda(float_t tlambda){this->_lambda = {tlambda}; this->variable_lambda = false;}

	void set_m(float_t m){this->mexp = m;}
	void set_n(float_t n){this->nexp = n;}

	template<class in_t>
	void set_variable_Kr(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Kr.clear(); this->_Kr.reserve(this->connector.nnodes);
		this->variable_Kr = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_Kr[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Ks(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Ks.clear(); this->_Ks.reserve(this->connector.nnodes);
		this->variable_Ks = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_Ks[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_precipitations(in_t& prec)
	{
		auto tprec = DAGGER::format_input(prec);
		this->_precipitations.clear(); this->_precipitations.reserve(this->connector.nnodes);
		this->variable_precipitations = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_precipitations[i] = tprec[i];
		}
	}

	template<class in_t>
	void set_variable_depcoeff(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_depcoeff.clear(); this->_depcoeff.reserve(this->connector.nnodes);
		this->variable_depcoeff = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_depcoeff[i] = tarr[i];
		}
	}


	template<class in_t>
	void set_variable_kappa_s(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_kappa_s.clear(); this->_kappa_s.reserve(this->connector.nnodes);
		this->variable_kappa_s = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_kappa_s[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_kappa_r(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_kappa_r.clear(); this->_kappa_r.reserve(this->connector.nnodes);
		this->variable_kappa_r = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_kappa_r[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Sc(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Sc.clear(); this->_Sc.reserve(this->connector.nnodes);
		this->variable_Sc = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_Sc[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Sc_M(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Sc_M.clear(); this->_Sc_M.reserve(this->connector.nnodes);
		this->variable_Sc_M = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_Sc_M[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Ke(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Ke.clear(); this->_Ke.reserve(this->connector.nnodes);
		this->variable_Ke = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_Ke[i] = tarr[i];
		}
	}


	template<class in_t>
	void set_variable_lambda(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_lambda.clear(); this->_lambda.reserve(this->connector.nnodes);
		this->variable_lambda = true;
		for(int i=0; i < this->graph.nnodes; ++i)
		{
				this->_lambda[i] = tarr[i];
		}
	}






	// ------------------------------------------------
	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   Access to param value           /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to accessing parameter values in a generic way
	// No matter which parametrisation method is used, these functions are called to get the parameter from its node indice
	// If spatial is not activated it returns the single value;
	// if spatial it returns the value at teh location
	// other option would be also called from here (e.g.K f(Qs))
	// ------------------------------------------------

	float_t Ks(int i)
	{
		if(variable_Ks == false)
			return this->_Ks[0];
		else
			return this->_Ks[i];
	}

	float_t Kr(int i)
	{
		if(variable_Kr == false)
			return this->_Kr[0];
		else
			return this->_Kr[i];
	}

	float_t kappa_s(int i)
	{
		if(variable_kappa_s == false)
			return this->_kappa_s[0];
		else
			return this->_kappa_s[i];
	}

	float_t kappa_r(int i)
	{
		if(variable_kappa_r == false)
			return this->_kappa_r[0];
		else
			return this->_kappa_r[i];
	}

	float_t precipitations(int i)
	{
		if(this->variable_precipitations == false)
			return this->_precipitations[0];
		else
			return this->_precipitations[i];
	}

	float_t depcoeff(int i)
	{
		if(variable_depcoeff == false)
			return this->_depcoeff[0];
		else
			return this->_depcoeff[i];
	}

	float_t Sc(int i)
	{
		if(variable_Sc == false)
			return this->_Sc[0];
		else
			return this->_Sc[i];
	}


	float_t Sc_M(int i)
	{
		if(variable_Sc_M == false)
			return this->_Sc_M[0];
		else
			return this->_Sc_M[i];
	}


	float_t Ke(int i)
	{
		if(variable_Ke == false)
			return this->_Ke[0];
		else
			return this->_Ke[i];
	}

	float_t lambda(int i)
	{
		if(variable_lambda == false)
			return this->_lambda[0];
		else
			return this->_lambda[i];
	}



	// ------------------------------------------------
	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   Initialisation METHODS          /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|
	//
	// All the methods related to initialising things
	// ------------------------------------------------

	void fill_up()
	{
		this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_fill;
		// std::vector<bool> testlinkchange(this->graph.links), nnodes2change(this->graph.nnodes,false);
		this->graph._compute_graph(this->z_surf, this->connector, true, false);
	}


	void init_vectors()
	{
		this->Qw = std::vector<float_t>(this->graph.nnodes,0.);

		if(this->fluvial)
		{
			this->Qs_sed = std::vector<float_t>(this->graph.nnodes,0.);
		}

		if(this->hillslopes || this->marine)
		{
			this->Qs_hs = std::vector<float_t>(this->graph.nnodes,0.);
		}


		if(this->TSP_module)
			this->init_Qs_TSP();
		
		if(this->Ch_MTSI)
			this->init_Qs_Ch_MTSI();
	}



	void run_SFD(float_t dt)
	{

		if(this->Ch_MTSI)
			this->Ch_MTSI_age(dt);

		this->init_vectors();
		this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		// std::vector<bool> testlinkchange(this->graph.links), nnodes2change(this->graph.nnodes,false);
		this->graph._compute_graph(this->z_surf, this->connector, true, false);

		// int n_links_changed = 0, n_links = int(this->graph.links.size()), nnodeschanged = 0;;
		// for(size_t i=0; i <this->graph.links.size();++i)
		// {
		// 	if(this->graph.links[i] != testlinkchange[i])
		// 	{
		// 		nnodes2change[this->graph.linknodes[i*2]] = true;
		// 		nnodes2change[this->graph.linknodes[i*2 + 1]] = true;
		// 		++n_links_changed;
		// 	}
		// }
		// for(auto v:nnodes2change)
		// {
		// 	if(v)++nnodeschanged;
		// }



		// std::cout << n_links_changed << "/" << n_links << " links modified - " << nnodeschanged << "/" << this->graph.nnodes << " nodes to recompute" << std::endl;

		for(int i= this->graph.nnodes -1; i>=0; --i)
		{
			// Getting geometrical info
			// # location
			int node = this->graph.Sstack[i];
			// # Aborting if outnode
			if(!this->graph.flow_out_or_pit(node,this->connector) == false)
				continue;
			// # receiver
			int rec = this->graph.Sreceivers[node];
			// # distance to receiver
			float_t dx = this->graph.Sdistance2receivers[node];
			// # cell area
			float_t cellarea = this->connector.get_area_at_node(node);
			// # local gradient
			float_t S = std::max((this->z_surf[node] - this->z_surf[rec])/dx, 1e-6);

			// Hydrology
			// # local addition
			this->Qw[node] += cellarea * precipitations(node);
			

			// preparing processes
			float_t fEs = 0, fEr = 0, fDs = 0, hEs = 0., hEr = 0., hDs = 0.;

			// Running the hillslope processes
			if(this->hillslopes)
			{
				// Storing the proportion of bedrock-power used
				float_t propused = 0;

				// If the slope is bellow the critical values
				if(S <= Sc(node) - 1e-6)
				{
					// If I have sediments
					if(this->h_sed[node] > 0)
					{
						
						// erosion power
						hEs = Ks(node) * S;

						// Checking I am not strapping more than the sediment pile
						if(hEs * dt > this->h_sed[node])
						{
							// Correcting to the max sed height removal possible and calculating the proportion of sediment used
							propused = this->h_sed[node]/dt /hEs;
							hEs = this->h_sed[node]/dt;
						}
						else
						{
							// I don't have any sed? -> no propused
							propused = 1.;
						}

					// end of the check for excessing sed height
					}

					// Remaining applied to bedrock
					hEr = (1. - propused) * Kr(node) * S;
					float_t L = cellarea/(1 - std::pow(S/Sc(node),2));
					hDs = this->Qs_hs[node]/L;

				}
				else
				{
					float_t tothE = (z_surf[node] - (dx * Sc(node) + z_surf[rec]) )/dt;
					if(tothE > this->h_sed[node])
					{
						hEs = this->h_sed[node]/dt;
						hEr = tothE - hEs;
					}
					else
						hEs = tothE;
				}

				if(this->TSP_module)
					this->apply_TSP(node, rec, hEs, hEr, hDs, dt, true);

				if(this->Ch_MTSI)
					this->apply_Ch_MTSI_SFD(node, rec, hEs, hEr, hDs, dt, true);

				this->Qs_hs[node] += (hEs + hEr - hDs) * cellarea;

				// Applying the fluxes
				this->h_sed[node] += (hDs - hEs) * dt;
				this->z_surf[node] += (hDs - hEs - hEr) * dt;
				if(this->Qs_hs[node] < 0) this->Qs_hs[node] = 0;
				this->Qs_hs[rec] += this->Qs_hs[node];

			}


			// Fluvial processes
			if(this->fluvial)
			{
				// Erosion
				float_t stream_power = std::pow(this->Qw[node],this->mexp) * std::pow(S,this->nexp);
				float_t propused = 0;
				if(this->h_sed[node] > 0)
				{
					fEs = Ks(node) * stream_power;
					if(this->h_sed[node] < fEs * dt)
					{
						propused = this->h_sed[node] / dt / fEs;
						fEs = this->h_sed[node] / dt;
					}
					else
					{
						propused = 1.;
					}
				}
				fEr = (1. - propused) * Kr(node) * stream_power;

				// Deposition
				float_t L = std::max(this->depcoeff(node) * this->Qw[node],cellarea);

				fDs = this->Qs_sed[node]/L;
				
				if(this->TSP_module)
					this->apply_TSP(node, rec, fEs, fEr, fDs, dt, false);

				if(this->Ch_MTSI)
					this->apply_Ch_MTSI_SFD(node, rec, fEs, fEr, fDs, dt, false);
			
				// Applying the fluxes
				this->h_sed[node] += (fDs - fEs) * dt;
				this->z_surf[node] += (fDs - fEs - fEr) * dt;
				this->Qs_sed[node] += (fEs + fEr - fDs) * cellarea;

				if(this->Qs_sed[node] < 0) this->Qs_sed[node] = 0;

				// Transferring to receivers
				this->Qw[rec] += this->Qw[node];
				this->Qs_sed[rec] += this->Qs_sed[node];

			}
		}

		
	}

	template<class out_t>
	out_t get_topo(){return DAGGER::format_output<std::vector<float_t>,out_t>(this->z_surf);}

	template<class out_t>
	out_t get_h_sed(){return DAGGER::format_output<std::vector<float_t>,out_t>(this->h_sed);}

	template<class in_t>
	void external_uplift(in_t& tU, float_t dt, bool apply_to_edges = false)
	{
		auto U = DAGGER::format_input(tU);
		for(int i=0; i < this->graph.nnodes; ++i)
		{
			if(!this->connector.boundaries.can_out(i) || apply_to_edges)
				this->z_surf[i] += U[i] * dt;
		}

	}

	void block_uplift(float_t rate, float_t dt)
	{
		for(int i=0; i < this->graph.nnodes; ++i)
		{
			if(this->connector.boundaries.can_out(i) == false)
				this->z_surf[i] += rate * dt;
		}
	}


	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// TSP module: Tracking Single Provenance
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

	template<class in_t>
	void init_TSP_module(float_t dz, in_t& tconcent)
	{
		this->TSP_module = true;
		auto concent = DAGGER::format_input(tconcent);
		this->TSP_concentrations = DAGGER::to_vec(concent);
		this->TSP_store = VerticalStorer<float_t, BasePropStorer<float_t> >(dz, this->graph.nnodes);
		this->init_Qs_TSP();

	}


	// Function reinitialising the trackers for the TSP module
	void init_Qs_TSP()
	{
		if(this->fluvial)
			this->TSP_Qsf = std::vector<BasePropStorer<float_t> >(this->graph.nnodes,BasePropStorer<float_t>()); 
		if(this->hillslopes)
			this->TSP_Qsh = std::vector<BasePropStorer<float_t> >(this->graph.nnodes,BasePropStorer<float_t>());
	}

	void apply_TSP(int i, int ir, float_t Es, float_t Er, float_t Ds, float_t dt, bool thillslopes)
	{

		Es *= dt;
		Er *= dt;
		Ds *= dt;

		// Removing from the sediment pile
		auto rem1 = this->TSP_store.remove(i, Es);
		
		// Depositing on the sediment pile
		BasePropStorer<float_t>& tpropr = (thillslopes) ? this->TSP_Qsh[i] : this->TSP_Qsf[i];
		this->TSP_store.pile_up(i,Ds,tpropr);

		// Modifying the fluxes
		float_t zfQ = (thillslopes) ? this->Qs_hs[i]/this->connector.get_area_at_node(i) : this->Qs_sed[i]/this->connector.get_area_at_node(i);
		zfQ *= dt;
		BasePropStorer<float_t>::mix(zfQ,tpropr, Es, rem1);
		zfQ += Es;
		BasePropStorer<float_t> bdr(this->TSP_concentrations[i]);
		BasePropStorer<float_t>::mix(zfQ,tpropr, Er, bdr);
		zfQ += Er;
		zfQ -= Ds;
		if(zfQ < 0)
			zfQ = 0;

		if(!this->graph.flow_out_or_pit(ir,this->connector))
		{
			if(thillslopes)
				BasePropStorer<float_t>::mix(this->Qs_hs[ir]/this->connector.get_area_at_node(i) * dt, this->TSP_Qsh[ir], zfQ, tpropr);
			else
				BasePropStorer<float_t>::mix(this->Qs_sed[ir]/this->connector.get_area_at_node(i) * dt, this->TSP_Qsf[ir], zfQ, tpropr);

		}
	
		



	}


	template<class out_t>
	out_t get_TSP_surface_concentrations()
	{
		if(this->TSP_module == false)
			throw std::runtime_error("Cannot return surface TSP if there is no TSP module activated (yo!)");

		std::vector<float_t> props(this->graph.nnodes, 0.);
		for(int i=0; i<this->graph.nnodes; ++i)
		{
			if(!this->graph.flow_out_or_pit(i,this->connector) == false)
				continue;

			if(this->TSP_store.pile[i].size() > 0)
				props[i] = this->TSP_store.pile[i].back().prop;

		}

		return DAGGER::format_output<std::vector<float_t>, out_t >(props);
	}

	template<class out_t>
	out_t sample_carrot_TSP(int row, int col)
	{
		int i = this->connector.nodeid_from_row_col(row,col);

		std::vector<float_t> carrot(this->TSP_store.pile[i].size());
		for(size_t j=0; j < carrot.size(); ++j)
		{
			carrot[j] = this->TSP_store.pile[i][j].prop;
		}

		return DAGGER::format_output<std::vector<float_t>, out_t >(carrot);

	}



	// Returns a cross section in the X (rwo = true) or Y (row = False) drection for nz number of depth cells
	template<class out_t>
	out_t get_transect_TSP(int rowcol, int nz, bool row = true)
	{
		// number of nodes in the desired dimentsion
		int nxy = (row)?this->connector.nx:this->connector.ny;
		// total number of nodes
		int nn = nxy * nz;

		std::vector<float_t> transect(nn,0.);
		int trow = rowcol;
		int tcol = rowcol;

		float_t maxz = std::numeric_limits<float_t>::min();
		for(int r=0;r<nxy;++r)
		{
			trow = (row)?rowcol:r;
			tcol = (row)?r:rowcol;
			int i = this->connector.nodeid_from_row_col(trow,tcol);
			if(this->z_surf[i] > maxz)
				maxz = this->z_surf[i];
		}

		
		for(int r=0;r<nxy;++r)
		{
			trow = (row)?rowcol:r;
			tcol = (row)?r:rowcol;
			int i = this->connector.nodeid_from_row_col(trow,tcol);

			float_t tz = maxz;
			int j2 = int(this->TSP_store.pile[i].size() - 1);;
			for(int j=0; j < nz; ++j)
			{
				int idxj = r * nz + j;
				if(tz>this->z_surf[i])
				{
					transect[idxj] = -1.;
				}
				else
				{
					transect[idxj] = (j2>0) ? this->TSP_store.pile[i][j2].prop : -1.;
					--j2;
				}

				tz-=this->TSP_store.dz;
			}
		}

		return DAGGER::format_output<std::vector<float_t>, out_t >(transect);

	}




	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// Ch_MTSI module: timing the time since incision
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~



	void init_Ch_MTSI(float_t dz)
	{
		this->Ch_MTSI = true;
		this->Ch_MTSI_store = VerticalStorer<float_t, BasePropStorer<float_t> >(dz, this->graph.nnodes);
		this->init_Qs_TSP();
	}


	// Function reinitialising the trackers for the Ch_MTSI module
	void init_Qs_Ch_MTSI()
	{
		if(this->fluvial)
			this->Ch_MTSI_Qsf = std::vector<BasePropStorer<float_t> >(this->graph.nnodes,BasePropStorer<float_t>()); 
		if(this->hillslopes)
			this->Ch_MTSI_Qsh = std::vector<BasePropStorer<float_t> >(this->graph.nnodes,BasePropStorer<float_t>());
	}


	void apply_Ch_MTSI_SFD(int i, int ir, float_t Es, float_t Er, float_t Ds, float_t dt, bool thillslopes)
	{

		Es *= dt;
		Er *= dt;
		Ds *= dt;

		// Removing from the sediment pile
		auto rem1 = this->Ch_MTSI_store.remove(i, Es);
		
		// Depositing on the sediment pile
		BasePropStorer<float_t>& tpropr = (thillslopes) ? this->Ch_MTSI_Qsh[i] : this->Ch_MTSI_Qsf[i];
		this->Ch_MTSI_store.pile_up(i,Ds,tpropr);

		// Modifying the fluxes
		float_t zfQ = (thillslopes) ? this->Qs_hs[i]/this->connector.get_area_at_node(i) : this->Qs_sed[i]/this->connector.get_area_at_node(i);
		zfQ *= dt;
		BasePropStorer<float_t>::mix(zfQ,tpropr, Es, rem1);
		zfQ += Es;
		BasePropStorer<float_t> bdr(0.);
		BasePropStorer<float_t>::mix(zfQ,tpropr, Er, bdr);
		zfQ += Er;
		zfQ -= Ds;
		if(zfQ < 0)
			zfQ = 0;

		if(!this->graph.flow_out_or_pit(ir,this->connector))
		{
			if(thillslopes)
				BasePropStorer<float_t>::mix(this->Qs_hs[ir]/this->connector.get_area_at_node(i) * dt, this->Ch_MTSI_Qsh[ir], zfQ, tpropr);
			else
				BasePropStorer<float_t>::mix(this->Qs_sed[ir]/this->connector.get_area_at_node(i) * dt, this->Ch_MTSI_Qsf[ir], zfQ, tpropr);
		}

	}

	void Ch_MTSI_age(float_t dt)
	{
		for(int i=0; i< this->graph.nnodes;++i)
		{
			for(size_t j=0; j<this->Ch_MTSI_store.pile[i].size(); ++j)
				this->Ch_MTSI_store.pile[i][j].prop += dt;
		}
	}



	template<class out_t>
	out_t get_Ch_MTIS_surface_age()
	{
		if(this->Ch_MTSI == false)
			throw std::runtime_error("Cannot return surface Ch_MTSI if there is no Ch_MTSI module activated (yo!)");

		std::vector<float_t> timeryolo(this->graph.nnodes, 0.);
		for(int i=0; i<this->graph.nnodes; ++i)
		{
			if(!this->graph.flow_out_or_pit(i,this->connector) == false)
				continue;

			if(this->Ch_MTSI_store.pile[i].size() > 0)
				timeryolo[i] = this->Ch_MTSI_store.pile[i].back().prop;

		}

		return DAGGER::format_output<std::vector<float_t>, out_t >(timeryolo);
	}

	template<class out_t>
	out_t sample_carrot_Ch_MTSI(int row, int col)
	{
		int i = this->connector.nodeid_from_row_col(row,col);

		std::vector<float_t> carrot(this->Ch_MTSI_store.pile[i].size());
		for(size_t j=0; j < carrot.size(); ++j)
		{
			carrot[j] = this->Ch_MTSI_store.pile[i][j].prop;
		}

		return DAGGER::format_output<std::vector<float_t>, out_t >(carrot);

	}


	template<class out_t>
	out_t get_transect_Ch_MTSI(int rowcol, int nz, bool row = true)
	{
		// number of nodes in the desired dimentsion
		int nxy = (row)?this->connector.nx:this->connector.ny;
		// total number of nodes
		int nn = nxy * nz;

		std::vector<float_t> transect(nn,0.);
		int trow = rowcol;
		int tcol = rowcol;

		float_t maxz = std::numeric_limits<float_t>::min();
		for(int r=0;r<nxy;++r)
		{
			trow = (row)?rowcol:r;
			tcol = (row)?r:rowcol;
			int i = this->connector.nodeid_from_row_col(trow,tcol);
			if(this->z_surf[i] > maxz)
				maxz = this->z_surf[i];
		}

		
		for(int r=0;r<nxy;++r)
		{
			trow = (row)?rowcol:r;
			tcol = (row)?r:rowcol;
			int i = this->connector.nodeid_from_row_col(trow,tcol);

			float_t tz = maxz;
			int j2 = int(this->Ch_MTSI_store.pile[i].size() - 1);;
			for(int j=0; j < nz; ++j)
			{
				int idxj = r * nz + j;
				if(tz>this->z_surf[i])
				{
					transect[idxj] = -1.;
				}
				else
				{
					transect[idxj] = (j2>0) ? this->Ch_MTSI_store.pile[i][j2].prop : -1.;
					--j2;
				}

				tz-=this->Ch_MTSI_store.dz;
			}
		}

		return DAGGER::format_output<std::vector<float_t>, out_t >(transect);

	}


	










	void rise_boundary_by(std::string wish, float_t woosh)
	{
		for(int i=0; i<this->graph.nnodes; ++i)
		{
			if(wish == "N" && this->connector.is_on_top_row(i))
				this->z_surf[i] += woosh;
			else if(wish == "E" && this->connector.is_on_rightest_col(i))
				this->z_surf[i] += woosh;
			else if(wish == "W" && this->connector.is_on_leftest_col(i))
				this->z_surf[i] += woosh;
			else if(wish == "S" && this->connector.is_on_bottom_row(i))
				this->z_surf[i] += woosh;
		}
	}


};



}
























#endif