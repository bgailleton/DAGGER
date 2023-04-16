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
	GRAPH_HYBRID,
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

	float_t courant_number_max  = 5e-4;
	float_t max_courant_dt_hydro = 5e-2;
	void set_courant_numer(float_t val){this->courant_number_max = val;};
	void set_max_courant_dt_hydro(float_t val){this->max_courant_dt_hydro = val;};
	float_t courant_dt_hydro = 1e-4;
	float_t get_courant_dt_hydro(){return this->courant_dt_hydro;}
	void enable_courant_dt_hydro(){this->mode_dt_hydro = PARAM_DT_HYDRO::COURANT;}


	MFD_PARTITIONNING weight_management = MFD_PARTITIONNING::PROPOSLOPE;
	void set_partition_method(MFD_PARTITIONNING& tmffmeth){this->weight_management = tmffmeth;}



	// Global constants:
	const float_t GRAVITY = 9.81, FIVETHIRD = 5./3., TWOTHIRD = 2./3., minslope = 0.;

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
	float_t Qwin_crit = 0.;
	void set_Qwin_crit(float_t val) {this->Qwin_crit = val; this->hydromode = HYDRO::GRAPH_HYBRID;}

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

	bool debugntopo = true;
	std::vector<float_t> DEBUGNTOPO;
	std::vector<float_t> get_nT(){return this->DEBUGNTOPO;}

	float_t debug_CFL = 0.;
	




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


	bool record_flowvec = false;
	std::vector<float_t> _rec_flowvec;
	void enable_flowvec_recording(){this->record_flowvec = true;};
	void disable_flowvec_recording(){this->record_flowvec = false; this->_rec_flowvec.clear();};
	template<class out_t>
	out_t get_flowvec_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_flowvec) ;}



	// EXPERIMENTAL
	std::vector<float_t> last_dt_prec;
	std::vector<float_t> last_dt_prec_e;
	std::vector<float_t> last_sw_prec;
	std::vector<float_t> last_dx_prec;
	float_t current_dt_prec = 0.;
	float_t current_dt_prec_e = 0.;
	float_t Vp = 0.;
	float_t Vps = 0.;
	// Create random device and Mersenne Twister engine
	std::random_device rd;
	std::mt19937 gen;
	// Create uniform integer distribution
	std::uniform_int_distribution<> dis;

	// ###################################### 
	// ###### Morpho monitorers #############
	// ###################################### 

	bool record_edot = false;
	std::vector<float_t> _rec_edot;
	void enable_edot_recording(){this->record_edot = true;};
	void disable_edot_recording(){this->record_edot = false; this->_rec_edot.clear();};
	template<class out_t>
	out_t get_edot_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_edot) ;}

	bool record_ddot = false;
	std::vector<float_t> _rec_ddot;
	void enable_ddot_recording(){this->record_ddot = true;};
	void disable_ddot_recording(){this->record_ddot = false; this->_rec_ddot.clear();};
	template<class out_t>
	out_t get_ddot_recording(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->_rec_ddot) ;}




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
	float_t dt_morpho_multiplier = 1.;
	void set_dt_morpho_multiplier(float_t val){this->dt_morpho_multiplier = val;}
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
	void set_topological_number(float_t val){this->topological_number = val;};
	float_t get_topological_number(){return this->topological_number;};

	// # ke (coefficient for erosion)
	PARAM_DT_HYDRO mode_dt_hydro = PARAM_DT_HYDRO::COURANT;
	std::vector<float_t> _dt_hydro = {1e-3};
	float_t get_dt_hydro() {return this->_dt_hydro[0];}


	bool hflow = false;



	graphflood(){this->gen = std::mt19937(this->rd()); };
	graphflood(Graph_t& graph, Connector_t& connector)
	{
		// Ingesting graph and connectors
		this->graph = &graph;
		this->connector = &connector;

		this->gen = std::mt19937(this->rd()); 
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

	template<class topo_t>
	void set_variable_ke(topo_t& variable_ke)
	{
		auto tke = format_input(variable_ke);
		this->_ke = to_vec(tke);
		this->mode_ke = PARAM_KE::VARIABLE;
	}


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
		else if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			return this->courant_dt_hydro;
		}
		else
			return this->_dt_hydro[0];
	}

	float_t dt_morpho(int i)
	{
		return this->dt_hydro(i) * this->dt_morpho_multiplier;
		// if (this->mode_dt_morpho == PARAM_DT_MORPHO::VARIABLE)
		// 	return this->_dt_morpho[i];
		// else if (this->mode_dt_morpho == PARAM_DT_MORPHO::HYDRO)
		// 	return this->dt_hydro(i);
		// else
		// 	return this->_dt_morpho[0];
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


	void block_uplift(float_t val)
	{
		for (int i=0;i<this->connector->nnodes;++i)
		{
			if(this->connector->boundaries.can_out(i) == false)
				this->_surface[i] += this->dt_morpho(i) * val;
		}
	}


	template<class topo_t>
	void variable_uplift(topo_t& ival)
	{	
		auto val = format_input(ival);
		for (int i=0;i<this->connector->nnodes;++i)
		{
			// if(this->connector->boundaries.can_out(i) == false)
			this->_surface[i] += this->dt_morpho(i) * val[i];
		}
	}


	// Main running function:
	void _run()
	{

		if(this->debugntopo)
			this->DEBUGNTOPO = std::vector<float_t>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition and other monitoring features:
		this->init_Qw();

		// reinitialising the sediments if needed
		if(this->morphomode != MORPHO::NONE)
			this->init_Qs();
		

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);


		// To be used if courant dt hydro is selected
		float_t tcourant_dt_hydro = std::numeric_limits<float_t>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<float_t> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<float_t>(this->graph->nnodes,0.);

		
		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::vector<float_t> weights(receivers.size(),0.), slopes(receivers.size(),0.);
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			if(this->connector->boundaries.no_data(node) || this->connector->flow_out_or_pit(node)) 
			{
				this->tot_Qwin_output += this->_Qw[node];
				if(this->connector->flow_out_model(node) == false && this->connector->boundaries.no_data(node) ==false )
				{
					int nn = this->connector->get_neighbour_idx(node,receivers);
					float_t hzh = this->_surface[node];
					for(int j=0; j < nn; ++j)
					{
						if(this->_surface[receivers[j]] > hzh)
							hzh = this->_surface[receivers[j]];
					}
					// vmot_hw[node] += this->_Qw[node]/this->connector->get_area_at_node(node);
					vmot_hw[node] += hzh - this->_surface[node];
				}
				continue;
			}


			// CFL calculator
			float_t sum_ui_over_dxi = 0.;

			// Deprecated test to switch dynamically between SFD and MFD
			if(this->hydromode == HYDRO::GRAPH_HYBRID)
				SF = this->_Qw[node] < this->Qwin_crit;


			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Slope max
			float_t Smax;

			//NOTE:
			// No need to calculate the topological number anymore
			float_t topological_number_v2 = 0.;

			// Calculating the slopes/weights and topological numbers
			if(SF == false)
			{
				Smax = this->weights_automator(receivers, weights, slopes, node, nrecs, topological_number_v2);
				if(topological_number_v2 == 0)
				{
					topological_number_v2 = 1.;
				}
				if(this->debugntopo)
				{
					this->DEBUGNTOPO[node] = topological_number_v2;
				}
			}
			else
			{
				Smax = this->get_Sw(node,this->connector->Sreceivers[node],this->connector->Sdistance2receivers[node],this->minslope);
			}

			// I can only have an hydraulic slope
			if(this->_hw[node] == 0)
				Smax = 0;

			// std::cout << nrecs << "|" << std::endl;;

			// Initialising the total Qout
			float_t total_Qout = 0.;
			float_t Qwin = this->_Qw[node];

			// precalculating the power
			float_t pohw = std::pow(this->_hw[node], this->FIVETHIRD);

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

				// if(topological_number_helper[rec] == 1)
				// 	throw std::runtime_error("ALREADY");

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
				
				if(SF) 
				{
					Sw=this->get_Sw(node,rec,dx,this->minslope);
					if(this->connector->flow_out_model(rec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
					{
						dx = this->connector->dx;
						dw = this->connector->dy;
						Sw = this->bou_fixed_val;
					}

					if(this->record_Sw)
						this->_rec_Sw[node] += Sw;
				}
				else
				{
					Sw = slopes[j];
					if(this->record_Sw)
						this->_rec_Sw[node] += Sw * weights[j];
				


					if(this->connector->flow_out_or_pit(rec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
					{
						// Sw = this->bou_fixed_val;
						// std::cout << "rec is " << rec << " and Sw is " << Sw << " hw is " << this->_hw[node] <<  std::endl;
						dx = this->connector->dx;
						dw = this->connector->dy;
					}
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
					tQout = (Smax > 0) ? topological_number_v2 * pohw * dw/this->mannings(node) * Sw/Smax:0;
					

					// std::cout << topological_number_v2 << "|";

					// tQout = (Smax > 0) ? ttoponum * pohw * dw/this->mannings(node) * Sw/Smax:0;
					// this->catch_nan(this->mannings(node), "this->mannings(node)");
					// this->catch_nan(pohw, "pohw");
					// this->catch_nan(Smax, "Smax");
					// this->catch_nan(Sw, "Sw");
					// this->catch_nan(tQout, "tQwout " + std::to_string(Smax));
				}	



				if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
					// sum_ui_over_dxi += tQout/(this->_hw[node] * dw)/dx;
					sum_ui_over_dxi += std::pow(this->_hw[node],2./3.)/this->mannings(node) * Sw/Smax /dx;


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
						// ++topological_number_helper[rec]; 

					}

				}
				else
				{
					this->_Qw[rec] += tQout;
				}

				if(this->connector->flow_out_model(rec))
					this->tot_Qw_output += tQout;


				// Now the morpho
				if(this->morphomode != MORPHO::NONE)
				{

					// Calculating the sediment flux going through this link
					float_t tQs = (SF) ? this->_Qs[node] : this->_Qs[node] * weights[j];

					if (this->connector->boundaries.forcing_io(node))
					{
						this->_Qs[rec] += tQs;
						if(this->connector->flow_out_or_pit(rec))
							this->tot_Qs_output += tQs;
						continue;
					}

					// initialising all the variables. e = erosion, d = deposition, _l is for the lateral ones and AB are the 2 lateral nodes
					float_t edot = 0., ddot = 0., eldot_A = 0., eldot_B = 0., dldot_A = 0., dldot_B = 0.;
					// Getting the transversal (perpendicular) dx to apply lateral erosion/deposition
					float_t tdl = this->connector->get_traverse_dx_from_links_idx(receivers[j]);
					// And gathering the orthogonal nodes
					std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,rec);

					// Calculating the shear stress for the given link
					float_t tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Sw;

					


					// Double checking the orthogonal nodes and if needs be to process them
					int oA = orthonodes.first;
					if(this->connector->boundaries.forcing_io(oA) || this->connector->is_in_bound(oA) == false ||  this->connector->boundaries.no_data(oA))
						oA = -1;
					int oB = orthonodes.second;
					if(this->connector->boundaries.forcing_io(oB) || this->connector->is_in_bound(oB) == false ||  this->connector->boundaries.no_data(oB))
						oB = -1;

					// Calculating local erosion rates, which depends on whether the shear stress exceeds the critical one
					if( tau > this->tau_c(node) && this->_hw[node] < 10)
					{
						edot += this->ke(node) * std::pow(tau - this->tau_c(node),this->aexp(node));
						
						if(this->record_edot)
							this->_rec_edot[node] += edot;
					}
					
					// And the local deposition, which depends on the transport distance
					ddot = tQs/this->kd(node);

					if(this->record_ddot)
						this->_rec_ddot[node] += ddot;

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

			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && sum_ui_over_dxi > 0)
			{
				float_t provisional_dt = this->courant_number_max/(sum_ui_over_dxi);
				// if(provisional_dt < tcourant_dt_hydro)
				tcourant_dt_hydro = std::min(provisional_dt, this->max_courant_dt_hydro);
			}

			vmot_hw[node] += (this->_Qw[node] - total_Qout)/this->connector->get_area_at_node(node);

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += total_Qout;
			
			// this->debug_CFL = std::max(this->debug_CFL, sum_ui_over_dxi * this->dt_hydro(node));

		}

		// std::cout << "Apply motions" << std::endl;

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<float_t>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}




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

		// std::cout << "CFL::" << this->debug_CFL << std::endl;

	}

	// Main running function (version 2)
	void run()
	{

		// Saving the topological number if needed
		if(this->debugntopo)
			this->DEBUGNTOPO = std::vector<float_t>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition and other monitoring features:
		this->init_Qw();

		// reinitialising the sediments if needed
		if(this->morphomode != MORPHO::NONE)
			this->init_Qs();
		
		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// To be used if courant dt hydro is selected
		float_t tcourant_dt_hydro = std::numeric_limits<float_t>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<float_t> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<float_t>(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<float_t,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers, vmot_hw)) continue;

			// CFL calculator
			float_t sum_ui_over_dxi = 0.;

			// Deprecated test to switch dynamically between SFD and MFD (does not really add anything and is buggy in rivers)
			if(this->hydromode == HYDRO::GRAPH_HYBRID)
				SF = this->_Qw[node] < this->Qwin_crit;

			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max
			float_t Smax;
			float_t dw0max;
			float_t dx;
			int recmax = node;

			//NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			float_t topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			float_t Qwin = this->_Qw[node];

			// precalculating the power
			float_t pohw = std::pow(this->_hw[node], this->TWOTHIRD);

			// Squarerooting Smax
			// float_t debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			float_t u_flow = pohw * sqrtSmax/this->mannings(node);
			// Volumetric discahrge
			float_t Qwout = dw0max * this->_hw[node] * u_flow;			

			// Eventually recording Smax
			if(this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// temp calc for courant
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] > 0)
				sum_ui_over_dxi = u_flow/dx;


			if(this->morphomode != MORPHO::NONE && this->connector->boundaries.forcing_io(node) == false)
			{
				float_t edot = 0., ddot = 0., eldot_A = 0., eldot_B = 0., dldot_A = 0., dldot_B = 0.;
				// float_t this->_Qs[node] = this->_Qs[node];
				float_t tdl = this->connector->get_travers_dy_from_dx(dx);
				// And gathering the orthogonal nodes
				std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,recmax);
				float_t tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Smax;
				// Double checking the orthogonal nodes and if needs be to process them
				int oA = orthonodes.first;
				if(this->connector->boundaries.forcing_io(oA) || this->connector->is_in_bound(oA) == false ||  this->connector->boundaries.no_data(oA))
					oA = -1;
				int oB = orthonodes.second;
				if(this->connector->boundaries.forcing_io(oB) || this->connector->is_in_bound(oB) == false ||  this->connector->boundaries.no_data(oB))
					oB = -1;

				if( tau > this->tau_c(node) && this->_hw[node] < 10)
				{
					edot += this->ke(node) * std::pow(tau - this->tau_c(node),this->aexp(node));
					
					if(this->record_edot)
						this->_rec_edot[node] += edot;
				}

				ddot = this->_Qs[node]/this->kd(node);


				// Dealing with lateral deposition if lS > 0 and erosion if lS <0
				if(oA >= 0 )
				{
					float_t tSwl = this->get_Stopo(node,oA, tdl);
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
					float_t tSwl = this->get_Stopo(node,oB, tdl);
					if(tSwl > 0)
					{
						dldot_B = tSwl * this->kd_lateral(node) * ddot;
					}
					else
					{
						eldot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
					}
				}
				float_t fbatch = (ddot + dldot_B + dldot_A - edot - eldot_A - eldot_B) * dx;
				this->_Qs[node] -= fbatch;

				// float_t corrector = 1.;
				float_t totd = ddot + dldot_B + dldot_A;
				if(this->_Qs[node] < 0) 
				{
					// std::cout << "happend?????????" << std::endl;
					// corrector = totd/(totd + corrector); 
					this->_Qs[node] = 0;
				}

				vmot[node] += ddot - edot;
				vmot[oA] += dldot_A;
				vmot[oA] -= eldot_A;
				vmot[oB] += dldot_B;
				vmot[oB] -= eldot_B;

			}



			// going through the receiver(s)
			for(int j =0; j<nrecs; ++j)
			{

				// Hydro for the link
				// universal params
				int rec;
				if(SF)
					rec = recmax; 
				else
					rec = this->connector->get_to_links(receivers[j]);
				
				if(rec < 0) continue;

				if(this->connector->flow_out_model(rec))
				{
					if(SF)
					{
						this->tot_Qw_output += Qwout; 
					}
					else if(weights[j] > 0 && Qwin > 0)
					{						
						this->tot_Qw_output += weights[j] * Qwout;
					}
				}

				

				if(this->hydrostationary)
				{
					if(SF)
					{
						this->_Qw[rec] += Qwin; 
					}
					else if(weights[j] > 0 && Qwin > 0)
					{						
						this->_Qw[rec] += weights[j] * Qwin;
					}

				}
				else
				{
					this->_Qw[rec] += Qwout * weights[j];
				}

				if(this->morphomode != MORPHO::NONE)
				{
					this->_Qs[rec] += (SF == false)?weights[j] * this->_Qs[node]:this->_Qs[node] ;
				}

			}

			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && sum_ui_over_dxi > 0)
			{
				float_t provisional_dt = this->courant_number_max/(sum_ui_over_dxi);
				// if(provisional_dt < tcourant_dt_hydro)
				tcourant_dt_hydro = std::min(provisional_dt, this->max_courant_dt_hydro);
			}

			vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;


			
		}



		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<float_t>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}


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
				if(this->connector->boundaries.can_give(i) && this->connector->flow_out_or_pit(i) == false)
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
				this->_rec_dhw = std::vector<float_t>(this->graph->nnodes,0.);

		if(record_flowvec)
			this->_rec_flowvec = std::vector<float_t>(this->graph->nnodes * 2,0.);
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


		// Initialising the Qs recorders

		if(this->record_edot)
			this->_rec_edot = std::vector<float_t>(this->connector->nnodes, 0.);

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

	// initial check for boundary conditions and eventually applying relevant changes
	bool _initial_check_boundary_pit(int& node, std::array<int,8>& receivers, std::vector<float_t>& vmot_hw)
	{
		if(this->connector->boundaries.no_data(node) || this->connector->flow_out_or_pit(node)) 
		{
			// Checking mass conservations
			this->tot_Qwin_output += this->_Qw[node];

			// I encountered a pit, that happened when LM are not preprocessed
			// Then I fill it slightly, but it does not really work
			if(this->connector->flow_out_model(node) == false && this->connector->boundaries.no_data(node) ==false )
			{
				int nn = this->connector->get_neighbour_idx(node,receivers);
				float_t hzh = this->_surface[node];
				for(int j=0; j < nn; ++j)
				{
					if(this->_surface[receivers[j]] > hzh)
						hzh = this->_surface[receivers[j]];
				}
				// vmot_hw[node] += this->_Qw[node]/this->connector->get_area_at_node(node);
				vmot_hw[node] += hzh - this->_surface[node];
			}
			return true;
		}
		return false;
	}


	// offsetting some of the calculations from the run to make it more undertandable adn reusable
	void _compute_slopes_weights_et_al(int& node, bool& SF, float_t& Smax, std::array<float_t,8>& slopes, 
	std::array<float_t,8>& weights,int& nrecs, std::array<int,8>& receivers, int& recmax, float_t& dx, float_t& dw0max, float_t& topological_number_v2)
	{
		// Calculating the slopes/weights and topological numbers
		if(SF == false)
		{

			// ROUTINES FOR Multiple FLow Directions
			// fill the arrays of slopes/weights/.. and calculate Smax, dmax, recmas and all
			Smax = this->weights_automator_v2(receivers, weights, slopes, node, nrecs, topological_number_v2, dw0max, recmax,dx);
			
			// Debugging, I guess deprecated TODO::CHECK
			if(topological_number_v2 == 0)
			{
				topological_number_v2 = 1.;
			}

			// Check the recording of topological number
			if(this->debugntopo)
			{
				this->DEBUGNTOPO[node] = topological_number_v2;
			}

		}
		else
		{

			// ROUTINES FOR Single FLow Directions
			// -> Slope
			Smax = this->get_Sw(node,this->connector->Sreceivers[node],this->connector->Sdistance2receivers[node],this->minslope);
			// -> rec
			recmax = this->connector->Sreceivers[node];
			// -> dy (integrated width)
			dw0max = this->connector->get_travers_dy_from_dx(this->connector->Sdistance2receivers[node]);
			// -> dx (flow distance)
			dx = this->connector->Sdistance2receivers[node];
			
			// boundary case
			if(this->connector->flow_out_model(recmax) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
			{
				dw0max = this->connector->dy;
				Smax = this->bou_fixed_val;
				dx = this->connector->dx;
			}
		}

		// I can only have an hydraulic slope if I have water
		if(this->_hw[node] == 0)
			Smax = 0;

		// Done
	}

	float_t weights_automator(std::array<int,8>& receivers, 
		std::vector<float_t>& weights, std::vector<float_t>& slopes, 
		int& node, 
		int& nrecs,
		float_t& topological_number_v2)
	{

		float_t sumw = 0., Smax = this->minslope, dw0max= 0 , sumSdw = 0.;;

		// int_dw = 0;
		// int yolo1 = node;
		// bool isit = yolo1 == 300;
		for(int i = 0; i < nrecs; ++i)
		{
			int lix = receivers[i];

			
			if(this->connector->is_link_valid(lix) == false)
			{
				continue;
			}

			// int_dw += this->connector->get_dx_from_links_idx(lix);
			float_t tdw = this->connector->get_traverse_dx_from_links_idx(lix);
			int rec = this->connector->get_to_links(lix);
			if(this->connector->flow_out_or_pit(rec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
			{
				slopes[i] = this->bou_fixed_val;		
				tdw	= this->connector->dy;
			}
			else
				slopes[i] = this->get_Sw(lix, this->minslope);

			sumSdw += slopes[i] * tdw;
			
			if(slopes[i]>Smax)
			{
				Smax = slopes[i];
				dw0max = tdw;
			}
			
			// slopes[i] = slope;
			if(this->weight_management == MFD_PARTITIONNING::PROPOSLOPE)
				weights[i] = slopes[i] * tdw;
			
			else if(this->weight_management == MFD_PARTITIONNING::SQRTSLOPE)
				weights[i] = std::sqrt(slopes[i] * tdw);
			
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

		topological_number_v2 = (Smax*dw0max)/sumSdw;

		return Smax;
	}

	float_t weights_automator_v2(std::array<int,8>& receivers, 
		std::array<float_t,8>& weights, std::array<float_t,8>& slopes, 
		int& node, 
		int& nrecs,
		float_t& topological_number_v2,
		float_t& dw0max,
		int& recmax,
		float_t& dx
		)
	{

		// Placeholders
		float_t sumw = 0., Smax = this->minslope, sumSdw = 0.;;
		std::pair<float_t,float_t> fvech;

		for(int i = 0; i < nrecs; ++i)
		{
			int lix = receivers[i];

			
			if(this->connector->is_link_valid(lix) == false)
			{
				continue;
			}

			// int_dw += this->connector->get_dx_from_links_idx(lix);
			float_t tdw = this->connector->get_traverse_dx_from_links_idx(lix);
			int rec = this->connector->get_to_links(lix);
			bool isboud = false;
			if(this->connector->flow_out_or_pit(rec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
			{
				slopes[i] = this->bou_fixed_val;		
				tdw	= this->connector->dy;
				isboud = true;
			}
			else
				slopes[i] = this->get_Sw(lix, this->minslope);

			sumSdw += slopes[i] * tdw;

			if(this->record_flowvec)
			{
				this->connector->get_dxdy_from_links_idx( lix, node, fvech, false);
				this->_rec_flowvec[node * 2] += slopes[i] * tdw * fvech.first;
				this->_rec_flowvec[node * 2 + 1] += slopes[i] * tdw * fvech.second;
			}

			
			if(slopes[i]>Smax)
			{
				Smax = slopes[i];
				dw0max = tdw;
				dx = (isboud)?this->connector->dx : this->connector->get_dx_from_links_idx(lix);
				recmax = rec;
			}
			
			// // slopes[i] = slope;
			// if(this->stochaslope)
			// 	weights[i] = this->randu.get();

			if(this->weight_management == MFD_PARTITIONNING::PROPOSLOPE)
				weights[i] = slopes[i] * tdw;
			
			else if(this->weight_management == MFD_PARTITIONNING::SQRTSLOPE)
				weights[i] = std::sqrt(slopes[i] * tdw);
			
			else if(this->weight_management == MFD_PARTITIONNING::PROPOREC)
				weights[i] = 1.;

			if(this->stochaslope)
				weights[i] += this->stochaslope_coeff * this->randu.get();


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

		if(this->record_flowvec)
		{

			float_t length = std::sqrt(std::pow(this->_rec_flowvec[node * 2],2) + std::pow(this->_rec_flowvec[node * 2 + 1],2));
			if(length!=0)
			{
				this->_rec_flowvec[node * 2] /= length;
				this->_rec_flowvec[node * 2 + 1] /= length;
			}
		}


		topological_number_v2 = (Smax*dw0max)/sumSdw;

		return Smax;
	}

	// small helper function returning node index from the right stack
	int get_istack_node(int i)
	{
		if(this->hydromode != HYDRO::GRAPH_SFD)
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

	float_t get_Stopo(int node, int rec, float_t dx)
	{
		return (this->_surface[node] - this->_hw[node] - (this->_surface[rec]  - this->_hw[rec]))/dx; 
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


	void run_precipitions(int N, float_t precdt)
	{
		if(this->last_dt_prec.size() == 0)
		{
			this->last_dt_prec = std::vector<float_t>(this->connector->nnodes,0.);
			this->last_dt_prec_e = std::vector<float_t>(this->connector->nnodes,0.);
			this->last_sw_prec = std::vector<float_t>(this->connector->nnodes,0.);
			this->last_dx_prec = std::vector<float_t>(this->connector->nnodes,0.);
		}

		this->compute_Vp(precdt);


		// current_dt_prec
		auto neighbours = this->connector->get_empty_neighbour();
		auto slopes = this->connector->get_empty_neighbour();


		for(int _ = 0; _<N; ++_)
		{

			int node = this->spawn_precipition();
			// std::cout << "init at " << node << std::endl;
			
			this->current_dt_prec += precdt;

			if(this->morphomode != MORPHO::NONE)
				this->current_dt_prec_e += precdt;

			float_t tQs = this->Vps * this->connector->dy;

			int NMAX = this->connector->nnodes * 2;
			while(true)
			{
				--NMAX;
				if(this->connector->boundaries.has_to_out(node) || NMAX == 0)
					break;

				// update step from Qwout
				float_t tdt = this->current_dt_prec - this->last_dt_prec[node];
				// this->last_dt_prec[node] = this->current_dt_prec;

				int rec = node;
				float_t Sw = 0.;
				float_t Sw_sel = 0.;
				float_t dx = 0.;
				float_t dy = 0.;
				int nn = this->connector->get_neighbour_idx_links(node,neighbours);

				for(int j=0; j<nn; ++j)
				{
					int trec = this->connector->get_other_node_from_links(neighbours[j], node);
					if(this->connector->boundaries.no_data(trec))
						continue;

					// if(this->last_dt_prec[trec] != this->current_dt_prec)
					this->update_hw_prec_and_dt(trec);
				}
				this->update_hw_prec_and_dt(node);

				float_t deltah = 0.;
				if(true)
				{
					deltah = this->Vp;
					this->_hw[node] += deltah;
					if(this->_hw[node] < 0.) 
					{
						deltah += this->_hw[node];
						this->_hw[node] = 0.;
					}
				}
				this->_surface[node] += deltah;
				

				for(int j=0; j<nn; ++j)
				{
					int trec = this->connector->get_other_node_from_links(neighbours[j], node);

					if(this->connector->boundaries.no_data(trec))
						continue;

					if(this->_surface[trec] <= this->_surface[node])
					{
						float_t tsw = (this->_surface[node] - this->_surface[trec])/this->connector->get_dx_from_links_idx(neighbours[j]);
						float_t tsw_sel = tsw * this->randu.get() * this->stochaslope;
						if(tsw_sel > Sw_sel)
						{
							Sw_sel = tsw_sel;
							Sw = tsw;
							rec = trec;
							dx = this->connector->get_dx_from_links_idx(neighbours[j]);
							dy = this->connector->get_travers_dy_from_dx(dx);
						}

						if(this->connector->boundaries.can_out(trec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
						{
							Sw_sel = this->bou_fixed_val;
							Sw = this->bou_fixed_val;
							rec = trec;
							dx = this->connector->dx;
							dy = this->connector->dy;
						}
					}
				}



				this->last_sw_prec[node] = Sw;
				this->last_dx_prec[node] = dx;


				
				if( this->morphomode != MORPHO::NONE && this->connector->boundaries.forcing_io(node) == false)
				{
					if(Sw>0 && rec != node)
					{
						tdt = (this->current_dt_prec_e - this->last_dt_prec_e[node]) * this->dt_morpho_multiplier;
						this->last_dt_prec_e[node] = this->current_dt_prec_e;
						float_t tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Sw;
						float_t edot = 0., ddot = 0., edot_A = 0., edot_B = 0, ddot_A = 0, ddot_B = 0.;
						if(tau > this->tau_c(node))
							edot += this->ke(node) * std::pow(tau - this->tau_c(node),this->aexp(node));
						
						ddot = tQs/this->kd(node);
						if(ddot * dx > tQs)
							ddot = tQs/dx;

						tQs += (edot - ddot) * dx;

						// if(edot > 0)
						// 	std::cout << "|e|" << edot;
						// if(ddot > 0)
						// 	std::cout << "|d|" << ddot;

						this->_surface[node] += (ddot - edot) * tdt;

						std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,rec);
						int oA = orthonodes.first;
						if(this->connector->boundaries.forcing_io(oA) || this->connector->is_in_bound(oA) == false ||  this->connector->boundaries.no_data(oA))
							oA = -1;
						int oB = orthonodes.second;
						if(this->connector->boundaries.forcing_io(oB) || this->connector->is_in_bound(oB) == false ||  this->connector->boundaries.no_data(oB))
							oB = -1;

						// Dealing with lateral deposition if lS > 0 and erosion if lS <0
						if(oA >= 0 )
						{
							float_t tSwl = this->get_Stopo(node,oA, dy);
							if(tSwl > 0)
							{
								ddot_A = tSwl * this->kd_lateral(node) * ddot;
							}
							else
							{
								edot_A = std::abs(tSwl) * this->ke_lateral(node) * edot;
							}

							this->_surface[oA] += (ddot_A -  edot_A) * tdt;
							tQs += (edot_A - ddot_A) * dy;
							// std::cout << edot_A << "|" << ddot_A  << std::endl;

						}

						if(oB >= 0 )
						{
							float_t tSwl = this->get_Stopo(node,oB, dy);
							if(tSwl > 0)
							{
								ddot_B = tSwl * this->kd_lateral(node) * ddot;
							}
							else
							{
								edot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
							}

							this->_surface[oB] += (ddot_B -  edot_B) * tdt;
							tQs += (edot_B - ddot_B) * dy;
						}


						if(tQs <0) tQs = 0;


					}
					// else
					// {
					// 	this->_surface[node] += this->Vp;
					// }

				}

				


				if(node == rec && this->connector->boundaries.can_out(node))
					break;

				node = rec;
			}

		}


	}

	void update_hw_prec_and_dt(int node)
	{

		float_t tdt = this->current_dt_prec - this->last_dt_prec[node];
		this->last_dt_prec[node] = this->current_dt_prec;

		if(tdt == 0 || this->last_sw_prec[node] == 0 || this->last_dx_prec[node] == 0)
		{
			return;
		}
		
		auto delta = tdt * std::pow(this->_hw[node],5./3.) * std::sqrt(this->last_sw_prec[node])/this->mannings(node)/this->last_dx_prec[node];
		this->_hw[node] -= delta;
		this->_surface[node] -= delta; 
	}

	void compute_Vp(float_t precdt)
	{
		this->Vp = 0. ; 
		if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT || this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
		{
		
			for(int i=0;i<this->connector->nnodes; ++i)
			{
				this->Vp+=this->precipitations(i) * precdt;
			}
			this->dis = std::uniform_int_distribution<> (0,this->connector->nnodes-1);
		}
		else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H)
		{
			for(size_t i=0; i<this->_water_entries.size(); ++i)
			{
				this->Vp+=this->_water_entries[i] * precdt;
			}
			this->dis = std::uniform_int_distribution<> (0,this->_water_entries.size()-1);
		
		}

		this->Vps = 0;
		if (this->sed_input_mode == SED_INPUT::ENTRY_POINTS_Q)
		{
			for(size_t i=0; i<this->_sed_entries.size(); ++i)
			{
				this->Vps+=this->_sed_entries[i] * precdt * this->dt_morpho_multiplier;
			}
			// this->dis = std::uniform_int_distribution<> (0,this->_water_entries.size()-1);
		
		}


	}


	int spawn_precipition()
	{

		int node;

		do
		{
			
			if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT || this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
			{
				node =  this->dis(this->gen);
			}
			else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H)
			{
				node =  this->_water_entry_nodes[this->dis(this->gen)];
			}
			else
				throw std::runtime_error("NOT DEFINED");

		}while(this->connector->boundaries.no_data(node));

		return node;

	}


};
// end of graphflood class






} // End of DAGGER namespace

































#endif