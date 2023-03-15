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
	REROUTE,
	FILL,
};


enum class MFD_PARTITIONNING
{
	PROPOSLOPE,
	SQRTSLOPE,
	PROPOREC,
};


enum class WATER_INPUT
{
	PRECIPITATIONS_CONSTANT,
	PRECIPITATIONS_VARIABLE,
	ENTRY_POINTS_H,
};

enum class SED_INPUT
{
	NONE,
	ENTRY_POINTS_Q,
};


enum class BOUNDARY_HW
{
	FIXED_HW,
	FIXED_SLOPE,
};

template<class float_t, class Graph_t, class Connector_t>
class graphflood
{
public:

	// Underlying grid:
	// std::shared_ptr<Graph_t> graph;
	// std::shared_ptr<Connector_t> connector;
	Graph_t* graph;
	Connector_t* connector;

	// Global modes
	HYDRO hydromode = HYDRO::GRAPH_MFD;
	MORPHO morphomode = MORPHO::NONE;
	// DEPRES depression_resolver = DEPRES::priority_flood; // MANAGED BY THE GRAPH!
	HYDROGRAPH_LM depression_management = HYDROGRAPH_LM::FILL;
	MFD_PARTITIONNING weight_management = MFD_PARTITIONNING::PROPOSLOPE;




	// Global constants:
	const float_t GRAVITY = 9.81, FIVETHIRD = 5./3., minslope = 1e-6;

	bool stochaslope = false;
	float_t stochaslope_coeff = 1.;
	void set_stochaslope(float_t val)
	{
		this->stochaslope = true;
		this->stochaslope_coeff = val;
		this->connector->set_stochaticiy_for_SFD(val);
	}
	void disable_stochaslope(){this->stochaslope = false;}

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

	BOUNDARY_HW boundhw = BOUNDARY_HW::FIXED_HW;
	float_t bou_fixed_val = 0.;
	void set_fixed_hw_at_boundaries(float_t val){this->boundhw = BOUNDARY_HW::FIXED_HW; this->bou_fixed_val = val;}
	void set_fixed_slope_at_boundaries(float_t val){this->boundhw = BOUNDARY_HW::FIXED_SLOPE; this->bou_fixed_val = val;}


	




	// ###################################### 
	// ###### Hydro monitorers ##############
	// ###################################### 

		// Low cost monitoring parameters
	float_t tot_Qw_input = 0;
	float_t get_tot_Qw_input()const{return this->tot_Qw_input;}

	float_t tot_Qwin_output = 0;
	float_t get_tot_Qwin_output()const{return this->tot_Qwin_output;}

	float_t tot_Qw_output = 0;
	float_t get_tot_Qw_output()const{return this->tot_Qw_output;}

	float_t tot_Qs_output = 0;
	float_t get_tot_Qs_output()const{return this->tot_Qs_output;}

	bool record_Qw_out = false;
	std::vector<float_t> _rec_Qwout;
	void enable_Qwout_recording(){this->record_Qw_out = true;};
	void disable_Qwout_recording(){this->record_Qw_out = false; this->_rec_Qwout.clear();};
	template<class out_t>
	out_t get_Qwout_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_Qwout) ;}

	bool record_Sw = false;
	std::vector<float_t> _rec_Sw;
	void enable_Sw_recording(){this->record_Sw = true;};
	void disable_Sw_recording(){this->record_Sw = false; this->_rec_Sw.clear();};
	template<class out_t>
	out_t get_Sw_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_Sw) ;}

	bool record_dhw = false;
	std::vector<float_t> _rec_dhw;
	void enable_dhw_recording(){this->record_dhw = true;};
	void disable_dhw_recording(){this->record_dhw = false; this->_rec_dhw.clear();};
	template<class out_t>
	out_t get_dhw_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_dhw) ;}

	bool record_filling = false;
	std::vector<float_t> _rec_filling;
	void enable_filling_recording(){this->record_filling = true;};
	void disable_filling_recording(){this->record_filling = false; this->_rec_filling.clear();};
	template<class out_t>
	out_t get_filling_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_filling) ;}


	// ###################################### 
	// Parameters for morpho ################
	// ###################################### 

	// # a exponent for erosion
	bool mode_aexp = false;
	std::vector<float_t> _aexp = {1.5};

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
	bool mode_rho = false;
	std::vector<float_t> _rho = {1000};;

	// # ke (coefficient for erosion)
	bool mode_tau_c = false;
	std::vector<float_t> _tau_c = {6};;
	
	// # ke (coefficient for erosion)
	PARAM_DT_MORPHO mode_dt_morpho = PARAM_DT_MORPHO::HYDRO;
	std::vector<float_t> _dt_morpho = {1e-3};

	SED_INPUT sed_input_mode = SED_INPUT::NONE;
	std::vector<int> _sed_entry_nodes;
	std::vector<float_t> _sed_entries;





	// Randomiser helper
	DAGGER::easyRand randu;

	// ###################################### 
	// Parameters for Hydrograph ############
	// ######################################

	// # ke (coefficient for erosion)
	bool mode_mannings = false;
	std::vector<float_t> _mannings = {0.033};

	// # ke (coefficient for erosion)
	WATER_INPUT water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
	std::vector<float_t> _precipitations = {1e-4};
	std::vector<int> _water_entry_nodes;
	std::vector<float_t> _water_entries;

	// float_t mannings = 0.033 
	float_t topological_number = 4./8;

	// # ke (coefficient for erosion)
	PARAM_DT_HYDRO mode_dt_hydro = PARAM_DT_HYDRO::CONSTANT;
	std::vector<float_t> _dt_hydro = {1e-3};
	float_t get_dt_hydro() {return this->_dt_hydro[0];}


	bool hflow = false;



	graphflood(){};
	graphflood(Graph_t& graph, Connector_t& connector)
	{
		// Ingesting graph and connectors
		this->graph = &graph;
		this->connector = &connector;
	}

	template<class topo_t>
	void set_topo(topo_t& topo)
	{
		auto tin = DAGGER::format_input(topo);
		std::vector<double> temp = DAGGER::to_vec(tin);
		if(this->_hw.size() == 0) this->_hw = std::vector<float_t>(this->graph->nnodes,0);
		this->_surface = std::vector<float_t>(temp);
		for(int i=0; i < this->graph->nnodes; ++i)
			this->_surface[i] -= this->_hw[i];
	}

	template<class topo_t>
	void set_hw(topo_t& thw)
	{
		auto tin = DAGGER::format_input(thw);
		std::vector<double> temp = DAGGER::to_vec(tin);
		
		if(this->_hw.size() == 0)
			this->_hw = std::move(temp);
		else
		{
			for(int i=0; i < this->graph->nnodes; ++i)
				this->_surface[i] += temp[i] - this->_hw[i];
			this->_hw = std::move(temp);
		}
		
	}


	void set_dt_hydro(float_t tdt)
	{
		this->mode_dt_hydro = PARAM_DT_HYDRO::CONSTANT;
		this->_dt_hydro = {tdt};
	}

	void set_dt_morpho(float_t tdt)
	{
		this->mode_dt_morpho = PARAM_DT_MORPHO::CONSTANT;
		this->_dt_morpho = {tdt};
	}

	template<class out_t, class in_t>
	void set_water_input_by_entry_points(out_t& hw_entry, in_t& hw_indices)
	{
		// preformatting the inputs
		auto tin = DAGGER::format_input(hw_entry);
		auto tin_idx = DAGGER::format_input(hw_indices);

		// Setting the general mode
		this->water_input_mode = WATER_INPUT::ENTRY_POINTS_H;

		this->_water_entry_nodes = DAGGER::to_vec(tin_idx);
		this->_water_entries = DAGGER::to_vec(tin);

		for(size_t i =0; i< this->_water_entry_nodes.size(); ++i)
		{
			this->connector->boundaries.codes[this->_water_entry_nodes[i]] = BC::FORCE_IN;
		}
	}

	// template<class out_t>
	void set_water_input_by_constant_precipitation_rate(float_t precipitations)
	{
		this->water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
		this->_precipitations = {precipitations};
	}

	template<class out_t>
	void set_water_input_by_variable_precipitation_rate(out_t& precipitations)
	{
		this->water_input_mode = WATER_INPUT::PRECIPITATIONS_VARIABLE;
		auto tin = format_input(precipitations);
		this->_precipitations = DAGGER::to_vec(tin);
	}

	template<class out_t, class in_t>
	void set_sed_input_by_entry_points(out_t& sed_entry, in_t& sed_indices)
	{
		// preformatting the inputs
		auto tin = DAGGER::format_input(sed_entry);
		auto tin_idx = DAGGER::format_input(sed_indices);

		// Setting the general mode
		this->sed_input_mode = SED_INPUT::ENTRY_POINTS_Q;

		this->_sed_entry_nodes = DAGGER::to_vec(tin_idx);
		this->_sed_entries = DAGGER::to_vec(tin);

	}




	float_t hw(int i){return this->_hw[i];}
	float_t Qw(int i){return this->_Qw[i];}
	float_t Qs(int i){return this->_Qs[i];}
	float_t surface(int i){return this->_surface[i];}


	void enable_MFD(){this->hydromode = HYDRO::GRAPH_MFD;}
	void enable_SFD(){this->hydromode = HYDRO::GRAPH_SFD;}

	void fill_minima(){this->depression_management = HYDROGRAPH_LM::FILL;}
	void reroute_minima(){this->depression_management = HYDROGRAPH_LM::REROUTE;}
	void ignore_minima(){this->depression_management = HYDROGRAPH_LM::IGNORE;}

	void enable_morpho(){this->morphomode = MORPHO::TL;}
	void disable_morpho(){this->morphomode = MORPHO::NONE;}


	void set_single_aexp(float_t ta){this->mode_aexp = false; this->_aexp = {ta};}
	void set_single_ke(float_t ta){this->mode_ke = PARAM_KE::CONSTANT; this->_ke = {ta};}
	void set_single_ke_lateral(float_t ta){this->mode_ke_lateral = false; this->_ke_lateral = {ta};}
	void set_single_kd(float_t ta){this->mode_kd = false; this->_kd = {ta};}
	void set_single_kd_lateral(float_t ta){this->mode_kd_lateral = false; this->_kd_lateral = {ta};}
	void set_single_tau_c(float_t ta){this->mode_tau_c = false; this->_tau_c = {ta};}


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
			return this->dt_hydro(i);
		else
			return this->_dt_morpho[0];
	}


	float_t tau_c(int i)
	{
		if (this->mode_tau_c)
			return this->_tau_c[i];
		else
			return this->_tau_c[0];
	}

	float_t rho(int i)
	{
		if (this->mode_rho)
			return this->_rho[i];
		else
			return this->_rho[0];
	}

	float_t mannings(int i)
	{
		if (this->mode_mannings)
			return this->_mannings[i];
		else
			return this->_mannings[0];
	}

	float_t precipitations(int i)
	{
		if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
			return this->_precipitations[i];
		else
			return this->_precipitations[0];
	}


	// Main running function:
	void run()
	{

		// Initialise the water discharge fields according to water input condition and other monitoring features:
		this->init_Qw();
		//
		if(this->morphomode != MORPHO::NONE)
			this->init_Qs();
	
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// Graph Processing
		this->graph_automator();


		std::vector<std::uint8_t> topological_number_helper;
		if(SF == false)
			topological_number_helper = std::vector<std::uint8_t>(this->connector->nnodes, 0);

		std::vector<float_t> vmot, vmot_hw(this->graph->nnodes,0.);

		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<float_t>(this->graph->nnodes,0.);

		// main loop
		auto receivers = this->connector->get_empty_neighbour();
		// auto neighbours = this->connector->get_empty_neighbour();
		std::vector<float_t> weights(receivers.size(),0.), slopes(receivers.size(),0.);

		// std::cout << "main loop" << std::endl;

		// std::cout << "SF::" << SF << std::endl;
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			// std::cout << "i::" << i << std::endl;
			// getting next node in line
			int node = this->get_istack_node(i);
			// std::cout << "node::" << node << std::endl;

			if(this->connector->boundaries.no_data(node) || this->connector->flow_out_or_pit(node)) 
			{
				this->tot_Qwin_output += this->_Qw[node];
				if(this->connector->boundaries.force_giving(node))
					throw std::runtime_error("Force_giving_is_pit?");
				continue;
			}

			int nvalidneighb = this->connector->get_n_potential_receivers(node);

			// Getting the receivers
			int nrecs,nrecs_eff = 0; 
			if(SF)
			{ 
				nrecs = 1;
			}
			else
			{
				nrecs = this->connector->get_receivers_idx_links(node,receivers);
			}


			float_t Smax;
			if(SF == false)
			{
				Smax = this->weights_automator(receivers, weights, slopes, node, nrecs);
				for(int j=0; j<nrecs; ++j)
				{
					if(weights[j] > 0 && this->_hw[node] > 0)
						++nrecs_eff;
				}
			}
			else
			{
				Smax = this->get_Sw(node,this->connector->Sreceivers[node],this->connector->Sdistance2receivers[node],this->minslope);
				// this->catch_nan(Smax, "hw node " + std::to_string(this->_hw[node]) + " and rec " + std::to_string(this->_hw[this->connector->Sreceivers[node]]) + " node: " + std::to_string(node) + " rec: " + std::to_string(this->connector->Sreceivers[node]));
			}

			if(this->_hw[node] == 0)
				Smax = 0;

			// std::cout << nrecs << "|" << std::endl;;

			// Initialising the total Qout
			float_t total_Qout = 0.;
			float_t Qwin = this->_Qw[node];

			// precalculating the power
			float_t pohw = std::pow(this->_hw[node], this->FIVETHIRD);
			
			float_t ttoponum = 1.;
			if(SF == false)
			{
				// int neff = nrecs_eff + int(topological_number_helper[node]);
				if(nvalidneighb > 4)
				{
					ttoponum = 4./nvalidneighb;
					// std::cout << "node: " << node << " toponum = " << ttoponum << std::endl;
				}
			}

			// Squarerooting Smax
			// float_t debug_S = Smax;
			Smax = std::sqrt(Smax);

			// going through the receiver(s)
			for(int j =0; j<nrecs; ++j)
			{

				// Hydro for the link
				// universal params
				int rec;
				if(SF)
				{
					rec=this->connector->Sreceivers[node]; 
				}
				else
				{
					rec = this->connector->get_to_links(receivers[j]);
					if(this->connector->boundaries.force_giving(rec))
						throw std::runtime_error("forcer is receiver?");
				}
				
				if(rec < 0) continue;

				// this->_Qw[rec] += this->_Qw[node];
				// continue;

				float_t dx; 
				if(SF)
				{
					dx=this->connector->Sdistance2receivers[node];
				}
				else
				{ 
					dx=this->connector->get_dx_from_links_idx(receivers[j]);
				}

				float_t dw = this->connector->get_travers_dy_from_dx(dx);


				// TEMPORARY OVERWRITTING
				// dx = this->connector->dx;
				// dw = this->connector->dy;


				float_t Sw; 
				if(this->connector->flow_out_or_pit(rec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
				{
					Sw = this->bou_fixed_val;
					Smax = std::sqrt(Sw);
					// std::cout << "rec is " << rec << " and Sw is " << Sw << " hw is " << this->_hw[node] <<  std::endl;
					if(this->record_Sw)
					{
						dx = this->connector->dx;
						dw = this->connector->dy;
						if(SF)
							this->_rec_Sw[node] += Sw;
						else
							this->_rec_Sw[node] += Sw * weights[j];
					}
				}
				else if(SF) 
				{
					Sw=this->get_Sw(node,rec,dx,this->minslope);
					if(this->record_Sw)
						this->_rec_Sw[node] += Sw;
				}
				else
				{
					Sw = slopes[j];
					if(this->record_Sw)
						this->_rec_Sw[node] += Sw * weights[j];
				}
				
				float_t tQout = 0.; 
				if(SF)
				{ 
					tQout = dw / this->mannings(node) * pohw * std::sqrt(Sw);
					// this->catch_nan(this->mannings(node), "this->mannings(node)");
					// this->catch_nan(pohw, "pohw");
					// this->catch_nan(Smax, "Smax");
					// this->catch_nan(tQout, "tQwout");

				}
				else
				{

					// tQout = (Smax > 0) ? this->topological_number * pohw * dw/this->mannings(node) * Sw/Smax:0;
					tQout = (Smax > 0) ? ttoponum * pohw * dw/this->mannings(node) * Sw/Smax:0;
					// this->catch_nan(this->mannings(node), "this->mannings(node)");
					// this->catch_nan(pohw, "pohw");
					// this->catch_nan(Smax, "Smax");
					// this->catch_nan(Sw, "Sw");
					// this->catch_nan(tQout, "tQwout " + std::to_string(Smax));
				}

				total_Qout += tQout;

				if(this->hydrostationary)
				{
					if(SF)
					{
						this->_Qw[rec] += Qwin; 
					}
					else if(weights[j] > 0 && Qwin > 0)
					{						

						// if(this->connector->boundaries.force_giving(node) && Qwin > 0)
						// 	std::cout << weights[j] * Qwin << " from " << node << " to " << rec << std::endl;

						this->_Qw[rec] += weights[j] * Qwin;
						++topological_number_helper[rec]; 

					}

					if(this->connector->flow_out_model(rec))
						this->tot_Qw_output += tQout;
				}
				else
				{
					this->_Qw[rec] += tQout;
				}

				// Now the morpho
				if(this->morphomode != MORPHO::NONE)
				{
					// initialising all the variables. e = erosion, d = deposition, _l is for the lateral ones and AB are the 2 lateral nodes
					float_t edot = 0., ddot = 0., eldot_A = 0., eldot_B = 0., dldot_A = 0., dldot_B = 0.;
					// Getting the transversal (perpendicular) dx to apply lateral erosion/deposition
					float_t tdl = this->connector->get_traverse_dx_from_links_idx(receivers[j]);
					// And gathering the orthogonal nodes
					std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,rec);

					// Calculating the shear stress for the given link
					float_t tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Sw;

					// Calculating the sediment flux going through this link
					float_t tQs = (SF) ? this->_Qs[node] : this->_Qs[node] * weights[j];


					// Double checking the orthogonal nodes and if needs be to process them
					int oA = orthonodes.first;
					if(this->connector->is_in_bound(oA))
					{
						if(this->connector->boundaries.forcing_io(oA)) oA = -1;
						if(this->connector->boundaries.can_receive(oA) == false || this->connector->boundaries.can_give(oA) == false) oA = -1;
					}
					else oA = -1;

					int oB = orthonodes.second;
					if(this->connector->is_in_bound(oB))
					{
						if(this->connector->boundaries.forcing_io(oB)) oB = -1;
						if(this->connector->boundaries.can_receive(oB) == false || this->connector->boundaries.can_give(oB) == false) oB = -1;
					}
					else oB = -1;
									

					// Calculating local erosion rates, which depends on whether the shear stress exceeds the critical one
					if( tau > this->tau_c(node))
						edot += this->ke(node) * std::pow(tau - this->tau_c(node),this->aexp(node));
					
					// And the local deposition, which depends on the transport distance
					ddot = tQs/this->kd(node);

					// Dealing with lateral deposition if lS > 0 and erosion if lS <0
					if(oA >= 0 )
					{
						float_t tSwl = this->get_Sw(node,oA, tdl);
						if(tSwl > 0)
						{
							dldot_A = tSwl * this->kd_lateral(node) * ddot;
						}
						else
						{
							eldot_A = std::abs(tSwl) * this->ke_lateral(node) * edot;
						}

					}

					if(oB >= 0 )
					{
						float_t tSwl = this->get_Sw(node,oB, tdl);
						if(tSwl > 0)
						{
							dldot_B = tSwl * this->kd_lateral(node) * ddot;
						}
						else
						{
							eldot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
						}
					}

					//
					float_t fbatch = (ddot + dldot_B + dldot_A - edot - eldot_A - eldot_B) * dx;
					tQs -= fbatch;

					float_t corrector = 1.;
					float_t totd = ddot + dldot_B + dldot_A;
					if(tQs < 0) 
					{
						// std::cout << "happend?????????" << std::endl;
						corrector = totd/(totd + corrector); 
						tQs = 0;
					}

					vmot[node] += ddot*corrector - edot;
					vmot[oA] += dldot_A*corrector;
					vmot[oA] -= eldot_A;
					vmot[oB] += dldot_B*corrector;
					vmot[oB] -= eldot_B;

					this->_Qs[rec] += tQs;

					if(this->connector->flow_out_or_pit(rec))
						this->tot_Qs_output += tQs;
				}				
			
			}

			vmot_hw[node] += (this->_Qw[node] - total_Qout)/this->connector->get_area_at_node(node);

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += total_Qout;

		}

		// std::cout << "Apply motions" << std::endl;


		// Applying vmots
		for(int i=0; i<this->graph->nnodes; ++i)
		{

			if(this->connector->flow_out_or_pit(i) && this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if(this->connector->boundaries.forcing_io(i)) continue;
			
			float_t tvh = vmot_hw[i] * this->dt_hydro(i);
			if(tvh < - this->_hw[i])
			{
				// std::cout << "HAPPENS::" << tvh << " vs " << this->_hw[i] << std::endl;
				tvh = - this->_hw[i];
			}

			this->_hw[i] += tvh;
			if(this->record_dhw)
				this->_rec_dhw[i] = tvh;
			
			this->_surface[i] += tvh;

			if(this->morphomode != MORPHO::NONE)
				this->_surface[i] += vmot[i] * this->dt_morpho(i);
		}

		// std::cout << "LOL" << std::endl;

	}

	void init_Qw()
	{	
		// resetting Qwin (needed to be stored at all time)
		this->_Qw = std::vector<float_t>(this->graph->nnodes,0.);
			
		// resetting global monitors (low cost)
		this->tot_Qw_input = 0;
		this->tot_Qw_output = 0;
		this->tot_Qwin_output = 0;

		if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT || this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
		{
			for(int i=0; i < this->graph->nnodes; ++i)
			{
				if(this->connector->boundaries.can_give(i))
				{
					this->_Qw[i] += this->precipitations(i) * this->connector->get_area_at_node(i);
					this->tot_Qw_input += this->_Qw[i];
				}
			}
		}
		else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H)
		{
			for(size_t i=0; i<this->_water_entries.size(); ++i)
			{
				int node = this->_water_entry_nodes[i];				
				this->_Qw[node] += this->_water_entries[i] * this->connector->get_area_at_node(node);
				this->tot_Qw_input += this->_Qw[node];
				// std::cout << node << " is given " << this->_water_entries[i] * this->connector->get_area_at_node(node);
			}
		}

		if(this->record_Qw_out)
			this->_rec_Qwout = std::vector<float_t>(this->graph->nnodes,0.);

		if(this->record_Sw)
			this->_rec_Sw = std::vector<float_t>(this->graph->nnodes,0.);
		if(this->record_dhw)
				this->_rec_dhw = std::vector<float_t>(this->graph->nnodes,0.);;
	}

	void init_Qs()
	{
		this->_Qs = std::vector<float_t>(this->graph->nnodes,0.);
		if(this->sed_input_mode == SED_INPUT::ENTRY_POINTS_Q)
		{
			for(size_t i=0; i<this->_sed_entries.size(); ++i)
			{
				int node = this->_sed_entry_nodes[i];				
				this->_Qs[node] += this->_sed_entries[i];
				// this->tot_Qw_input += this->_Qs[node];

			}
		}

		this->tot_Qs_output = 0.;
	}


	void graph_automator()
	{

		// is SS?
		bool only_SD = (this->hydromode == HYDRO::GRAPH_SFD);

		// making sure it has the right depression solver (SHOULD BE MOVED TO THE GRAPH MANAGEMENT LATER)
		if(this->depression_management == HYDROGRAPH_LM::IGNORE)
			this->graph->set_LMR_method(DEPRES::none);

		if(this->record_filling)
			this->_rec_filling = std::vector<float_t>(this->graph->nnodes,0.);

		// std::cout << "COMPUTING GRAPH " << only_SD << std::endl;

		// preformatting post-topo
		std::vector<float_t> post_topo(this->_surface.size(),0);
		for(int i=0; i<this->graph->nnodes; ++i)
		{
			post_topo[i] = this->_surface[i];
		}

		this->graph->_compute_graph(post_topo, only_SD, false);

		// fill water where depressions have been solved
		if(this->depression_management == HYDROGRAPH_LM::FILL)
		{
			for(int i=0; i<this->graph->nnodes; ++i)
			{
				if(this->connector->boundaries.no_data(i) || this->connector->boundaries.forcing_io(i))
					continue;

				if(this->_surface[i] < post_topo[i])
				{
					
					float_t ddhhee = post_topo[i] - this->_surface[i];
					
					if(this->record_filling)
						this->_rec_filling[i] = ddhhee;

					this->_hw[i] += ddhhee;

					this->_surface[i] = post_topo[i];
				}
			}
		}
	}

	float_t weights_automator(std::array<int,8>& receivers, std::vector<float_t>& weights, std::vector<float_t>& slopes, int& node, int& nrecs)
	{

		float_t sumw = 0., Smax = this->minslope;
		// int yolo1 = node;
		// bool isit = yolo1 == 300;
		for(int i = 0; i < nrecs; ++i)
		{
			int lix = receivers[i];

			
			if(this->connector->is_link_valid(lix) == false)
			{

				continue;
			}

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

		float_t sumf = 0.;
		for(int i = 0; i<nrecs;++i)
		{
			int lix = receivers[i];

			if(this->connector->is_link_valid(lix) == false)
				continue;

			weights[i] = weights[i]/sumw;
			sumf += weights[i];
		}


		return Smax;
	}

	// small helper function returning node index from the right stack
	int get_istack_node(int i)
	{
		if(this->hydromode == HYDRO::GRAPH_MFD)
			return this->graph->stack[i];
		else 
			return this->graph->Sstack[i];
	}

	float_t get_Sw(int lix, float_t minslope)
	{
		int from,to; this->connector->from_to_from_link_index(lix, from, to);
		return std::max((this->_surface[from] - this->_surface[to])/this->connector->get_dx_from_links_idx(lix),minslope); 
	}

	float_t get_Sw(int node, int rec, float_t dx, float_t minslope)
	{
		return std::max((this->_surface[node] - this->_surface[rec])/dx,minslope); 
	}

	float_t get_Sw(int node, int rec, float_t dx)
	{
		return (this->_surface[node] - this->_surface[rec])/dx; 
	}



	// GET DATA OUT
	template<class out_t>
	out_t get_hw(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_hw) ;}
	
	template<class out_t>
	out_t get_surface_topo(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_surface) ;}

	template<class out_t>
	out_t get_bedrock_topo()
	{
		std::vector<float_t> diff(this->_surface);
		
		for(int i=0; i< this->graph->nnodes; ++i)
			diff[i] -= this->_hw[i];

		return DAGGER::format_output<std::vector<float_t>, out_t >(diff) ;

	}

	template<class out_t>
	out_t get_Qwin(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_Qw) ;}

	template<class out_t>
	out_t get_SSTACKDEBUG(){ return DAGGER::format_output<std::vector<size_t>, out_t >(this->graph->Sstack) ;}



	void catch_nan(float_t testval, std::string error_message)
	{
		if(std::isfinite(testval) == false)
			throw std::runtime_error(error_message);
	}


};
// end of graphflood class






} // End of DAGGER namespace

































#endif