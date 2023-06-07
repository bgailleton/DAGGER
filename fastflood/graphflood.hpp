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

template<class fT, class Graph_t, class Connector_t>
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

	fT courant_number  = 5e-4;
	fT max_courant_dt_hydro = 1e3;
	fT min_courant_dt_hydro = 1e-4;
	void set_courant_numer(fT val){this->courant_number = val;};
	void set_max_courant_dt_hydro(fT val){this->max_courant_dt_hydro = val;};
	void set_min_courant_dt_hydro(fT val){this->min_courant_dt_hydro = val;};

	fT courant_dt_hydro = -1;
	fT get_courant_dt_hydro(){return this->courant_dt_hydro;}
	void enable_courant_dt_hydro(){this->mode_dt_hydro = PARAM_DT_HYDRO::COURANT;}


	MFD_PARTITIONNING weight_management = MFD_PARTITIONNING::PROPOSLOPE;
	void set_partition_method(MFD_PARTITIONNING& tmffmeth){this->weight_management = tmffmeth;}



	// Global constants:
	const fT GRAVITY = 9.81, FIVETHIRD = 5./3., TWOTHIRD = 2./3., minslope = 0.;

	bool stochaslope = false;
	fT stochaslope_coeff = 1.;
	void set_stochaslope(fT val)
	{
		this->stochaslope = true;
		this->stochaslope_coeff = val;
		this->connector->set_stochaticiy_for_SFD(val);
	}
	void disable_stochaslope(){this->stochaslope = false;}

	bool hydrostationary = true;
	void enable_hydrostationnary(){this->hydrostationary = true;};
	void disable_hydrostationnary(){this->hydrostationary = false;};
	fT Qwin_crit = 0.;
	void set_Qwin_crit(fT val) {this->Qwin_crit = val; this->hydromode = HYDRO::GRAPH_HYBRID;}

	// vecotr of node size
	// #  Hydraulic Surface elevation (bedrock + sed + water)
	std::vector<fT> _surface;

	// # Water depth
	std::vector<fT> _hw;
	// # Water discharge
	std::vector<fT> _Qw;

	// # Sediment discahrge
	std::vector<fT> _Qs;
	// # Sediment height
	std::vector<fT> _hs;

	BOUNDARY_HW boundhw = BOUNDARY_HW::FIXED_HW;
	fT bou_fixed_val = 0.;
	void set_fixed_hw_at_boundaries(fT val){this->boundhw = BOUNDARY_HW::FIXED_HW; this->bou_fixed_val = val;}
	void set_fixed_slope_at_boundaries(fT val){this->boundhw = BOUNDARY_HW::FIXED_SLOPE; this->bou_fixed_val = val;}

	bool debugntopo = true;
	std::vector<fT> DEBUGNTOPO;
	std::vector<fT> get_nT(){return this->DEBUGNTOPO;}

	fT debug_CFL = 0.;
	

	std::vector<std::uint8_t> converged;




	// ###################################### 
	// ###### Hydro monitorers ##############
	// ###################################### 

		// Low cost monitoring parameters
	fT tot_Qw_input = 0;
	fT get_tot_Qw_input()const{return this->tot_Qw_input;}

	fT tot_Qwin_output = 0;
	fT get_tot_Qwin_output()const{return this->tot_Qwin_output;}

	fT tot_Qw_output = 0;
	fT get_tot_Qw_output()const{return this->tot_Qw_output;}

	fT tot_Qs_output = 0;
	fT get_tot_Qs_output()const{return this->tot_Qs_output;}

	bool record_Qw_out = false;
	std::vector<fT> _rec_Qwout;
	void enable_Qwout_recording(){this->record_Qw_out = true;};
	void disable_Qwout_recording(){this->record_Qw_out = false; this->_rec_Qwout.clear();};
	template<class out_t>
	out_t get_Qwout_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_Qwout) ;}

	bool record_Sw = false;
	std::vector<fT> _rec_Sw;
	void enable_Sw_recording(){this->record_Sw = true;};
	void disable_Sw_recording(){this->record_Sw = false; this->_rec_Sw.clear();};
	template<class out_t>
	out_t get_Sw_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_Sw) ;}

	bool record_dhw = false;
	std::vector<fT> _rec_dhw;
	void enable_dhw_recording(){this->record_dhw = true;};
	void disable_dhw_recording(){this->record_dhw = false; this->_rec_dhw.clear();};
	template<class out_t>
	out_t get_dhw_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_dhw) ;}

	bool record_filling = false;
	std::vector<fT> _rec_filling;
	void enable_filling_recording(){this->record_filling = true;};
	void disable_filling_recording(){this->record_filling = false; this->_rec_filling.clear();};
	template<class out_t>
	out_t get_filling_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_filling) ;}


	bool record_flowvec = false;
	std::vector<fT> _rec_flowvec;
	void enable_flowvec_recording(){this->record_flowvec = true;};
	void disable_flowvec_recording(){this->record_flowvec = false; this->_rec_flowvec.clear();};
	template<class out_t>
	out_t get_flowvec_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_flowvec) ;}

	fT tau_max = 0.;

	// EXPERIMENTAL STUFF
	// EXPE OTHER STUFF
	std::vector<fT> last_Smax;

	// # EXPE PRECIPITIONS
	std::vector<fT> last_dt_prec;
	std::vector<fT> last_dt_prec_e;
	std::vector<fT> last_sw_prec;
	std::vector<fT> last_dx_prec;
	fT current_dt_prec = 0.;
	fT current_dt_prec_e = 0.;
	fT Vp = 0.;
	fT Vps = 0.;
	// Create random device and Mersenne Twister engine
	std::random_device rd;
	std::mt19937 gen;
	// Create uniform integer distribution
	std::uniform_int_distribution<> dis;

	// ###################################### 
	// ###### Morpho monitorers #############
	// ###################################### 

	bool record_edot = false;
	std::vector<fT> _rec_edot;
	void enable_edot_recording(){this->record_edot = true;};
	void disable_edot_recording(){this->record_edot = false; this->_rec_edot.clear();};
	template<class out_t>
	out_t get_edot_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_edot) ;}

	bool record_ddot = false;
	std::vector<fT> _rec_ddot;
	void enable_ddot_recording(){this->record_ddot = true;};
	void disable_ddot_recording(){this->record_ddot = false; this->_rec_ddot.clear();};
	template<class out_t>
	out_t get_ddot_recording(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_rec_ddot) ;}




	// ###################################### 
	// Parameters for morpho ################
	// ###################################### 

	// # a exponent for erosion
	bool mode_aexp = false;
	std::vector<fT> _aexp = {1.5};

	// # ke (coefficient for erosion)
	PARAM_KE mode_ke = PARAM_KE::CONSTANT;
	std::vector<fT> _ke = {1e-4};

	// # ke (coefficient for erosion)
	bool mode_ke_lateral = false;
	std::vector<fT> _ke_lateral = {0.1};

	// # ke (coefficient for erosion)
	bool mode_kd = false;
	std::vector<fT> _kd = {100};;

	// # kd (coefficient for erosion)
	bool mode_kd_lateral = false;
	std::vector<fT> _kd_lateral = {0.1};;

	// # ke (coefficient for erosion)
	bool mode_rho = false;
	std::vector<fT> _rho = {1000};;

	// # ke (coefficient for erosion)
	bool mode_tau_c = false;
	std::vector<fT> _tau_c = {6};;
	
	// # ke (coefficient for erosion)
	PARAM_DT_MORPHO mode_dt_morpho = PARAM_DT_MORPHO::HYDRO;
	fT dt_morpho_multiplier = 1.;
	void set_dt_morpho_multiplier(fT val){this->dt_morpho_multiplier = val;}
	std::vector<fT> _dt_morpho = {1e-3};

	SED_INPUT sed_input_mode = SED_INPUT::NONE;
	std::vector<int> _sed_entry_nodes;
	std::vector<fT> _sed_entries;





	// Randomiser helper
	DAGGER::easyRand randu;

	// ###################################### 
	// Parameters for Hydrograph ############
	// ######################################

	// # ke (coefficient for erosion)
	bool mode_mannings = false;
	std::vector<fT> _mannings = {0.033};

	// # ke (coefficient for erosion)
	WATER_INPUT water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
	std::vector<fT> _precipitations = {1e-4};
	std::vector<int> _water_entry_nodes;
	std::vector<fT> _water_entries;

	// fT mannings = 0.033 
	fT topological_number = 4./8;
	void set_topological_number(fT val){this->topological_number = val;};
	fT get_topological_number(){return this->topological_number;};

	// # ke (coefficient for erosion)
	PARAM_DT_HYDRO mode_dt_hydro = PARAM_DT_HYDRO::COURANT;
	std::vector<fT> _dt_hydro = {1e-3};
	fT get_dt_hydro() {return this->_dt_hydro[0];}


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
		if(this->_hw.size() == 0) this->_hw = std::vector<fT>(this->graph->nnodes,0);
		this->_surface = std::vector<fT>(temp);
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


	void set_dt_hydro(fT tdt)
	{
		this->mode_dt_hydro = PARAM_DT_HYDRO::CONSTANT;
		this->_dt_hydro = {tdt};
	}

	void set_dt_morpho(fT tdt)
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
	void set_water_input_by_constant_precipitation_rate(fT precipitations)
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




	fT hw(int i){return this->_hw[i];}
	fT Qw(int i){return this->_Qw[i];}
	fT Qs(int i){return this->_Qs[i];}
	fT surface(int i){return this->_surface[i];}


	void enable_MFD(){this->hydromode = HYDRO::GRAPH_MFD;}
	void enable_SFD(){this->hydromode = HYDRO::GRAPH_SFD;}

	void fill_minima(){this->depression_management = HYDROGRAPH_LM::FILL;}
	void reroute_minima(){this->depression_management = HYDROGRAPH_LM::REROUTE;}
	void ignore_minima(){this->depression_management = HYDROGRAPH_LM::IGNORE;}

	void enable_morpho(){this->morphomode = MORPHO::TL;}
	void disable_morpho(){this->morphomode = MORPHO::NONE;}


	void set_single_aexp(fT ta){this->mode_aexp = false; this->_aexp = {ta};}
	void set_single_ke(fT ta){this->mode_ke = PARAM_KE::CONSTANT; this->_ke = {ta};}
	void set_single_ke_lateral(fT ta){this->mode_ke_lateral = false; this->_ke_lateral = {ta};}
	void set_single_kd(fT ta){this->mode_kd = false; this->_kd = {ta};}
	void set_single_kd_lateral(fT ta){this->mode_kd_lateral = false; this->_kd_lateral = {ta};}
	void set_single_tau_c(fT ta){this->mode_tau_c = false; this->_tau_c = {ta};}

	template<class topo_t>
	void set_variable_ke(topo_t& variable_ke)
	{
		auto tke = format_input(variable_ke);
		this->_ke = to_vec(tke);
		this->mode_ke = PARAM_KE::VARIABLE;
	}


	fT aexp(int i)
	{
		if(this->mode_aexp)
			return this->_aexp[i];
		else
			return this->_aexp[0];
	}

	fT ke(int i)
	{
		if(PARAM_KE::VARIABLE == this->mode_ke)
			return this->_ke[i];
		else
			return this->_ke[0];
	}

	fT ke_lateral(int i)
	{
		if(this->mode_ke_lateral)
			return this->_ke_lateral[i];
		else
			return this->_ke_lateral[0];
	}

	fT kd(int i)
	{
		if (this->mode_kd)
			return this->_kd[i];
		else
			return this->_kd[0];
	}

	fT kd_lateral(int i)
	{
		if (this->mode_kd_lateral)
			return this->_kd_lateral[i];
		else
			return this->_kd_lateral[0];
	}

	fT dt_hydro(int i)
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

	fT dt_morpho(int i)
	{
		return this->dt_hydro(i) * this->dt_morpho_multiplier;
		// if (this->mode_dt_morpho == PARAM_DT_MORPHO::VARIABLE)
		// 	return this->_dt_morpho[i];
		// else if (this->mode_dt_morpho == PARAM_DT_MORPHO::HYDRO)
		// 	return this->dt_hydro(i);
		// else
		// 	return this->_dt_morpho[0];
	}


	fT tau_c(int i)
	{
		if (this->mode_tau_c)
			return this->_tau_c[i];
		else
			return this->_tau_c[0];
	}

	fT rho(int i)
	{
		if (this->mode_rho)
			return this->_rho[i];
		else
			return this->_rho[0];
	}

	fT mannings(int i)
	{
		if (this->mode_mannings)
			return this->_mannings[i];
		else
			return this->_mannings[0];
	}

	fT precipitations(int i)
	{
		if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
			return this->_precipitations[i];
		else
			return this->_precipitations[0];
	}


	void block_uplift(fT val)
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
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

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
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<fT>(this->graph->nnodes,0.);

		
		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::vector<fT> weights(receivers.size(),0.), slopes(receivers.size(),0.);
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			if(this->connector->boundaries.no_data(node) || this->connector->flow_out_or_pit(node) || (this->hydrostationary == false && this->_hw[node] == 0 && this->_Qw[node] == 0)) 
			{
				this->tot_Qwin_output += this->_Qw[node];
				if(this->connector->flow_out_model(node) == false && this->connector->boundaries.no_data(node) ==false )
				{
					// int nn = this->connector->get_neighbour_idx(node,receivers);
					// fT hzh = this->_surface[node];
					// for(int j=0; j < nn; ++j)
					// {
					// 	if(this->_surface[receivers[j]] > hzh)
					// 		hzh = this->_surface[receivers[j]];
					// }
					// vmot_hw[node] += this->_Qw[node]/this->connector->get_area_at_node(node);
					vmot_hw[node] = this->_Qw[node]/this->connector->get_area_at_node(node);
					if(vmot_hw[node]>0) std::cout << vmot_hw[node] * this->dt_hydro(node) << "|" << this->_Qw[node] << std::endl;; 
				}
				continue;
			}


			// CFL calculator
			fT sum_ui_over_dxi = 0.;

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
			fT Smax;

			//NOTE:
			// No need to calculate the topological number anymore
			fT topological_number_v2 = 0.;

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
			fT total_Qout = 0.;
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->FIVETHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
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

				fT dx; 
				if(SF)
				{
					dx=this->connector->Sdistance2receivers[node];
				}
				else
				{ 
					dx=this->connector->get_dx_from_links_idx(receivers[j]);
				}

				fT dw = this->connector->get_travers_dy_from_dx(dx);


				// TEMPORARY OVERWRITTING
				// dx = this->connector->dx;
				// dw = this->connector->dy;


				fT Sw; 
				
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
				
				fT tQout = 0.; 
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
					fT tQs = (SF) ? this->_Qs[node] : this->_Qs[node] * weights[j];

					if (this->connector->boundaries.forcing_io(node))
					{
						this->_Qs[rec] += tQs;
						if(this->connector->flow_out_or_pit(rec))
							this->tot_Qs_output += tQs;
						continue;
					}

					// initialising all the variables. e = erosion, d = deposition, _l is for the lateral ones and AB are the 2 lateral nodes
					fT edot = 0., ddot = 0., eldot_A = 0., eldot_B = 0., dldot_A = 0., dldot_B = 0.;
					// Getting the transversal (perpendicular) dx to apply lateral erosion/deposition
					fT tdl = this->connector->get_traverse_dx_from_links_idx(receivers[j]);
					// And gathering the orthogonal nodes
					std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,rec);

					// Calculating the shear stress for the given link
					fT tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Sw;

					


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
						fT tSwl = this->get_Sw(node,oA, tdl);
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
						fT tSwl = this->get_Sw(node,oB, tdl);
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
					fT fbatch = (ddot + dldot_B + dldot_A - edot - eldot_A - eldot_B) * dx;
					tQs -= fbatch;

					fT corrector = 1.;
					fT totd = ddot + dldot_B + dldot_A;
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
				fT provisional_dt = this->courant_number/(sum_ui_over_dxi);
				// if(provisional_dt < tcourant_dt_hydro)
				tcourant_dt_hydro = std::min(provisional_dt, this->max_courant_dt_hydro);
				tcourant_dt_hydro = std::max(provisional_dt, this->min_courant_dt_hydro);
			}

			vmot_hw[node] += (this->_Qw[node] - total_Qout)/this->connector->get_area_at_node(node);

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += total_Qout;
			
			// this->debug_CFL = std::max(this->debug_CFL, sum_ui_over_dxi * this->dt_hydro(node));

		}

		// std::cout << "Apply motions" << std::endl;

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<fT>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}




		// Applying vmots
		for(int i=0; i<this->graph->nnodes; ++i)
		{

			if(this->connector->flow_out_or_pit(i) && this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if(this->connector->boundaries.forcing_io(i)) continue;
			
			fT tvh = vmot_hw[i] * this->dt_hydro(i);
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

		// if(this->hydrostationary)
		if(true)
		{
			// std::cout << "_run_hydrostationnary" << std::endl;
			this->_run_hydrostationnary();
			// this->_run_hydrostationnary_lock();
			// this->_run_hydrostationnary_subgraph_test();
		}
		else
			this->_run_dynamics();

	}

	void _run_hydrostationnary()
	{
		// Saving the topological number if needed
		if(this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		this->tau_max = 0.;

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
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<fT>(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers, vmot_hw)) continue;

			// if(this->hydrostationary == false)
			// {
			// 	if(this->dt_hydro(node) > 0)
			// 		this->_hw[node] += this->dt_hydro(node) * this->precipitations(node);
			// 	else
			// 		this->_hw[node] += 1e-3 * this->precipitations(node);

			// }

			// CFL calculator
			// fT sum_ui_over_dxi = 0.;

			// Deprecated test to switch dynamically between SFD and MFD (does not really add anything and is buggy in rivers)
			// if(this->hydromode == HYDRO::GRAPH_HYBRID)
			// 	SF = this->_Qw[node] < this->Qwin_crit;

			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			//NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);
			// std:: cout << "pohw:" << pohw << "|" << this->_hw[node] << std::endl;

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax/this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;			
			// std::cout << Qwout << "|";

			// Eventually recording Smax
			if(this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// temp calc for courant
			// if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] > 0 && dx >0)
			// 	sum_ui_over_dxi = u_flow/dx;

			// Automates the computation of morpho LEM if needed
			this->_compute_morpho(node, recmax, dx, Smax, vmot);


			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(nrecs, recmax, node,  SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0)
			{
				fT provisional_dt =  (this->courant_number * dx)/u_flow;
				provisional_dt = std::min(provisional_dt, this->max_courant_dt_hydro);
				provisional_dt = std::max(provisional_dt, this->min_courant_dt_hydro);
				tcourant_dt_hydro = std::min(provisional_dt,tcourant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			if(this->hydrostationary)
				vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);
			else
				vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;
		}

		if(this->tau_max > 500)
			std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			// std::cout << "tcourant_dt_hydro::" << tcourant_dt_hydro << std::endl;
			if(this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<fT>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}


		// Applying vmots with the right dt
		this->_compute_vertical_motions(vmot_hw, vmot);
	}







	void _run_hydrostationnary_lock()
	{
		// Saving the topological number if needed
		if(this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

		if(this->converged.size() == 0)
			this->converged = std::vector<std::uint8_t>(connector->nnodes, 10);

		this->debug_CFL = 0.;

		this->tau_max = 0.;

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
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<fT>(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers, vmot_hw)) continue;

			// if(this->hydrostationary == false)
			// {
			// 	if(this->dt_hydro(node) > 0)
			// 		this->_hw[node] += this->dt_hydro(node) * this->precipitations(node);
			// 	else
			// 		this->_hw[node] += 1e-3 * this->precipitations(node);

			// }

			// CFL calculator
			// fT sum_ui_over_dxi = 0.;

			// Deprecated test to switch dynamically between SFD and MFD (does not really add anything and is buggy in rivers)
			// if(this->hydromode == HYDRO::GRAPH_HYBRID)
			// 	SF = this->_Qw[node] < this->Qwin_crit;

			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			//NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);
			// std:: cout << "pohw:" << pohw << "|" << this->_hw[node] << std::endl;

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax/this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;

			if(std::abs(1 - Qwout/this->_Qw[node]) < 0.05 && this->converged[node] > 0)
			{
				--this->converged[node];
				
			}
			// std::cout << Qwout << "|";

			// Eventually recording Smax
			if(this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// temp calc for courant
			// if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] > 0 && dx >0)
			// 	sum_ui_over_dxi = u_flow/dx;

			// Automates the computation of morpho LEM if needed
			this->_compute_morpho(node, recmax, dx, Smax, vmot);


			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(nrecs, recmax, node,  SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0)
			{
				fT provisional_dt =  (this->courant_number * dx)/u_flow;
				provisional_dt = std::min(provisional_dt, this->max_courant_dt_hydro);
				provisional_dt = std::max(provisional_dt, this->min_courant_dt_hydro);
				tcourant_dt_hydro = std::min(provisional_dt,tcourant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			if(this->hydrostationary)
			{
				if(this->converged[node] >0)
					vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);
			}
			else
			{
				if(this->converged[node] >0)
					vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);
			}

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;
		}

		if(this->tau_max > 500)
			std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			// std::cout << "tcourant_dt_hydro::" << tcourant_dt_hydro << std::endl;
			if(this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<fT>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}


		// Applying vmots with the right dt
		this->_compute_vertical_motions(vmot_hw, vmot);
	}




	void _run_hydrostationnary_averaged_test()
	{
		// Saving the topological number if needed
		if(this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		this->tau_max = 0.;

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
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<fT>(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers, vmot_hw)) continue;

			// if(this->hydrostationary == false)
			// {
			// 	if(this->dt_hydro(node) > 0)
			// 		this->_hw[node] += this->dt_hydro(node) * this->precipitations(node);
			// 	else
			// 		this->_hw[node] += 1e-3 * this->precipitations(node);

			// }

			// CFL calculator
			// fT sum_ui_over_dxi = 0.;

			// Deprecated test to switch dynamically between SFD and MFD (does not really add anything and is buggy in rivers)
			// if(this->hydromode == HYDRO::GRAPH_HYBRID)
			// 	SF = this->_Qw[node] < this->Qwin_crit;

			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			//NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);
			// std:: cout << "pohw:" << pohw << "|" << this->_hw[node] << std::endl;

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax/this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;			
			// std::cout << Qwout << "|";

			// Eventually recording Smax
			if(this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// temp calc for courant
			// if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] > 0 && dx >0)
			// 	sum_ui_over_dxi = u_flow/dx;

			// Automates the computation of morpho LEM if needed
			this->_compute_morpho(node, recmax, dx, Smax, vmot);


			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(nrecs, recmax, node,  SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0)
			{
				fT provisional_dt =  (this->courant_number * dx)/u_flow;
				tcourant_dt_hydro = std::min(provisional_dt, this->max_courant_dt_hydro);
				tcourant_dt_hydro = std::max(provisional_dt, this->min_courant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			if(this->hydrostationary)
				vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);
			else
				vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;
		}

		if(this->tau_max > 500)
			std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			// std::cout << "tcourant_dt_hydro::" << tcourant_dt_hydro << std::endl;
			if(this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<fT>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}


		// Applying vmots with the right dt
		this->_compute_vertical_motions_averaged_test(vmot_hw, vmot);
	}

	void _run_hydrostationnary_subgraph_test()
	{/*
		// 

		std::stack<int, std::vector<int> > nQ;
		std::vector<std::uint8_t> donsdone(this->connector->nnodes, 0);
		std::vector<std::uint8_t> recsdone(this->connector->nnodes, 0);
		std::vector<std::uint8_t> isproc(this->connector->nnodes, false);
		std::priority_queue< PQ_helper<int,double>, std::vector<PQ_helper<int,double> >, std::greater<PQ_helper<int,double> > > 
		std::vector<int> stack(this->connector->nnodes, false);

		auto neighbours = this->connector.get_empty_neighbour();
		for (int i=0; i<this->connector->nnodes; ++i)
		{
			if(this->connector->boundaries.no_data(i)) {isproc[i] = true; continue;};
			int nn = this->connector.get_neighbour_idx(i, neighbours);
			for(int j=0; j<nn;++j)
			{
				if(this->_surface[neighbours[j]] == this->_surface[i]) this->_surface[i] += this->connector->randu->get() * 1e-8;

				if(this->_surface[neighbours[j]] > this->_surface[i]) ++donsdone[i]; else ++recsdone[i];
			}

			if(donsdone[i] == 0)	nQ.emplace(i);
		}

		// STOPPED HERE
		std::array<fT,8> weights, slopes;
		
		// main loop
		while(nQ.empty() == false)
		{
			
			// Getting next node in line
			int node = nQ.yop(i);
			nQ.pop();

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax/this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;			
			// std::cout << Qwout << "|";

			// Eventually recording Smax
			if(this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// temp calc for courant
			// if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] > 0 && dx >0)
			// 	sum_ui_over_dxi = u_flow/dx;

			// Automates the computation of morpho LEM if needed
			this->_compute_morpho(node, recmax, dx, Smax, vmot);


			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(nrecs, recmax, node,  SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0)
			{
				fT provisional_dt =  (this->courant_number * dx)/u_flow;
				tcourant_dt_hydro = std::min(provisional_dt, this->max_courant_dt_hydro);
				tcourant_dt_hydro = std::max(provisional_dt, this->min_courant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			if(this->hydrostationary)
				vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node) * tcourant_dt_hydro;
			else
				vmot_hw[node] += (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node) * tcourant_dt_hydro;

			if(this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;
		}

		if(this->tau_max > 500)
			std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			// std::cout << "tcourant_dt_hydro::" << tcourant_dt_hydro << std::endl;
			if(this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<fT>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}


		// Applying vmots with the right dt
		this->_compute_vertical_motions(vmot_hw, vmot, false);
	
	*/}


	void _run_dynamics()
	{
		// Saving the topological number if needed
		if(this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);
		// std::cout << "adgfskjlfakmj" << std::endl;

		this->debug_CFL = 0.;

		this->tau_max = 0.;

		// Graph Processing
		// this->graph_automator();


		if(this->_Qw.size() == 0)
		{
			// std::cout << "adgfskjlfakmj" << std::endl;
			// throw std::runtime_error("BAGULLE");
			this->graph_automator();
			this->init_Qw();
			// this->_Qw = std::vector<fT>(this->connector->nnodes, 0.);
			for(int i=0; i<this->connector->nnodes; ++i)
				this->_Qw[i] = this->precipitations(i)  * this->connector->get_area_at_node(i);
			
		}
		

		// this->connector->_reallocate_vectors();
		// this->connector->update_links(this->_surface);
		this->graph_automator();
		for(int i=0; i<this->connector->nnodes; ++i)
				this->_Qw[i] = this->precipitations(i)  * this->connector->get_area_at_node(i);



		// Initialise the water discharge fields according to water input condition and other monitoring features:
		// this->init_Qw();
		std::vector<fT> nQ(this->connector->nnodes, 0.);
		// std::vector<fT> nQ(this->_Qw);

		// // reinitialising the sediments if needed
		// if(this->morphomode != MORPHO::NONE)
		// 	this->init_Qs();
		
		// Am I in SFD or MFD
		// bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// To be used if courant dt hydro is selected
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> vmot, vmot_hw(this->graph->nnodes,0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if(this->morphomode != MORPHO::NONE)
			vmot = std::vector<fT>(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		// auto receivers = this->connector->get_empty_neighbour();
		// std::array<fT,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			// if (this->_initial_check_boundary_pit(node, receivers, vmot_hw)) continue;
			if(this->connector->boundaries.no_data(node) || this->connector->flow_out_model(node)) continue;

			// // CFL calculator
			fT sum_ui_over_dxi = 0.;

			// // Deprecated test to switch dynamically between SFD and MFD (does not really add anything and is buggy in rivers)
			// // if(this->hydromode == HYDRO::GRAPH_HYBRID)
			// // 	SF = this->_Qw[node] < this->Qwin_crit;

			// // Getting the receivers
			int nrecs; 
			// if(SF)
				nrecs = 1;
			// else
			// 	nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max

			int recmax = this->connector->Sreceivers[node];
			if(recmax == node || this->connector->flow_out_model(recmax))
			{
				if(recmax == node)
					vmot_hw[node] += this->dt_hydro(node) * this->_Qw[node]/this->connector->get_area_at_node(node);
					// vmot_hw[node] = 0;
				
				continue;
			}

			fT dx = this->connector->Sdistance2receivers[node];
			fT Smax = this->get_Sw(node,recmax,dx);

			if(Smax<0)
			{
				std::cout <<  this->_surface[node] <<"|" << this->_surface[recmax] << std::endl;;
				continue;
			}
			fT dw0max = this->connector->get_travers_dy_from_dx(dx);


			// NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			// fT topological_number_v2 = 0.;

			// this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			// fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax/this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;

			// // Eventually recording Smax
			// if(this->record_Sw)
			// 	this->_rec_Sw[node] = Smax;

			// temp calc for courant
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] > 0 && dx >0)
				sum_ui_over_dxi = u_flow/dx;

			// Automates the computation of morpho LEM if needed
			// this->_compute_morpho(node, recmax, dx, Smax, vmot);


			// transfer fluxes
			// Computing the transfer of water and sed
			// this->_compute_transfers(nrecs, recmax, node,  SF, receivers, weights, Qwin, Qwout);
			nQ[recmax] += Qwout;
			// nQ[node] += Qwout;

			// Computing courant based dt
			fT provisional_dt = this->dt_hydro(node);
			if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && sum_ui_over_dxi > 0)
			{
				provisional_dt =  this->courant_number/(sum_ui_over_dxi * nrecs);
				tcourant_dt_hydro = std::min(provisional_dt, this->max_courant_dt_hydro);
				tcourant_dt_hydro = std::max(provisional_dt, this->min_courant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			vmot_hw[node] += provisional_dt * (this->_Qw[node] - Qwout)/this->connector->get_area_at_node(node);
			// this->_Qw[node] = this->precipitations(node)  * this->connector->get_area_at_node(node);
			// vmot_hw[node] = 0;
			// if(vmot_hw[node] < 0) 
			// {
			// 	std::cout << "RECASTING" << std::endl;
			// 	vmot_hw[node] = 0;
			// }

			if(std::isfinite(vmot_hw[node]) == false)
			{
				std::cout << recmax << "|" << node << "|" << Qwout << "|" << this->_Qw[node] << "|" << Smax << "|" << dx << std::endl;
				throw std::runtime_error("BIBIBIBITE");
			}

			// if(this->record_Qw_out)
			// 	this->_rec_Qwout[node] += Qwout;
		}
		// throw std::runtime_error("break");

		// END OF MAIN LOOP

		// Computing final courant based dt

		if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT)
		{
			if(this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if(tcourant_dt_hydro > 0 && tcourant_dt_hydro != std::numeric_limits<fT>::max() )
				this->courant_dt_hydro = tcourant_dt_hydro;
		}


		// Applying vmots with the right dt
		this->_compute_vertical_motions(vmot_hw, vmot,  false);
		// this->_Qw = std::move(nQ);
		for(int i=0; i<this->connector->nnodes; ++i)
		{
			// if(this->vmot_hw[i] > 0)
				this->_Qw[i] = nQ[i] + this->precipitations(i) * this->connector->get_area_at_node(i) ; 
			// std::max(this->_Qw[i],nQ[i]);
		}
	}




	void run_exp()
	{
		if(this->hydrostationary)
			this->_run_exp_stationnary();
		else
			this->_run_dynamics();
	}



	// Main running function (experimental 2)
	void _run_exp_stationnary()
	{
		

	}


	// ##########################
	// ##########################
	// ##########################
	// HELPER FUNCTIONS #########
	// ##########################
	// ##########################
	// ##########################


	void init_Qw()
	{	
		// resetting Qwin (needed to be stored at all time)
		this->_Qw = std::vector<fT>(this->graph->nnodes,0.);
			
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
			this->_rec_Qwout = std::vector<fT>(this->graph->nnodes,0.);

		if(this->record_Sw)
			this->_rec_Sw = std::vector<fT>(this->graph->nnodes,0.);
		if(this->record_dhw)
				this->_rec_dhw = std::vector<fT>(this->graph->nnodes,0.);

		if(record_flowvec)
			this->_rec_flowvec = std::vector<fT>(this->graph->nnodes * 2,0.);
	}

	void init_Qs()
	{
		this->_Qs = std::vector<fT>(this->graph->nnodes,0.);
		if(this->sed_input_mode == SED_INPUT::ENTRY_POINTS_Q)
		{
			for(size_t i=0; i<this->_sed_entries.size(); ++i)
			{
				int node = this->_sed_entry_nodes[i];				
				this->_Qs[node] += this->_sed_entries[i] * this->connector->dy;
				// this->tot_Qw_input += this->_Qs[node];

			}
		}


		// Initialising the Qs recorders

		if(this->record_edot)
			this->_rec_edot = std::vector<fT>(this->connector->nnodes, 0.);

		this->tot_Qs_output = 0.;
	}



	/// Automates the processing of the graph
	/// Function of global parameters for flow topology, local minima management, ...
	void graph_automator()
	{

		// is SS?
		bool only_SD = (this->hydromode == HYDRO::GRAPH_SFD);

		// making sure it has the right depression solver (SHOULD BE MOVED TO THE GRAPH MANAGEMENT LATER)
		if(this->depression_management == HYDROGRAPH_LM::IGNORE)
			this->graph->set_LMR_method(DEPRES::none);

		if(this->record_filling)
			this->_rec_filling = std::vector<fT>(this->graph->nnodes,0.);

		// preformatting post-topo
		std::vector<fT> post_topo(this->_surface.size(),0);
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
					
					fT ddhhee = post_topo[i] - this->_surface[i];
					
					if(this->record_filling)
						this->_rec_filling[i] = ddhhee;

					this->_hw[i] += ddhhee;

					this->_surface[i] = post_topo[i];
				}
			}
		}
	}

	// initial check for boundary conditions and eventually applying relevant changes
	bool _initial_check_boundary_pit(int& node, std::array<int,8>& receivers, std::vector<fT>& vmot_hw)
	{
		if(this->connector->boundaries.no_data(node) || this->connector->flow_out_or_pit(node) || (this->hydrostationary == false && this->_hw[node] == 0 && this->_Qw[node] == 0) ) 
		{
			// Checking mass conservations
			this->tot_Qwin_output += this->_Qw[node];

			// I encountered a pit, that happened when LM are not preprocessed
			// Then I fill it slightly, but it does not really work
			if(this->connector->flow_out_model(node) == false && this->connector->boundaries.no_data(node) ==false )
			{
				// int nn = this->connector->get_neighbour_idx(node,receivers);
				// fT hzh = this->_surface[node];
				// for(int j=0; j < nn; ++j)
				// {
				// 	if(this->_surface[receivers[j]] > hzh)
				// 		hzh = this->_surface[receivers[j]];
				// }
				vmot_hw[node] = this->_Qw[node]/this->connector->get_area_at_node(node);
				// if(vmot_hw[node]>0)
				// {
				// 	std::cout << this->_Qw[node] << "??::" << vmot_hw[node] * this->dt_hydro(node) << " ||| " << this->dt_hydro(node) << std::endl;;
				// 	if(vmot_hw[node]> 1e6)
				// 		throw std::runtime_error("flkflsdfkjdsf");
				// }
				// vmot_hw[node] += hzh - this->_surface[node];
			}
			return true;
		}
		return false;
	}

	// initial check for boundary conditions WITHOUT applying relevant changes
	bool _initial_check_boundary_pit(int& node, std::array<int,8>& receivers)
	{
		if(this->connector->boundaries.no_data(node) || this->connector->flow_out_or_pit(node)) 
		{
			// Checking mass conservations
			this->tot_Qwin_output += this->_Qw[node];
			return true;
		}
		return false;
	}


	// offsetting some of the calculations from the run to make it more undertandable adn reusable
	void _compute_slopes_weights_et_al(int& node, bool& SF, fT& Smax, std::array<fT,8>& slopes, 
	std::array<fT,8>& weights,int& nrecs, std::array<int,8>& receivers, int& recmax, fT& dx, fT& dw0max, fT& topological_number_v2)
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

	fT weights_automator(std::array<int,8>& receivers, 
		std::vector<fT>& weights, std::vector<fT>& slopes, 
		int& node, 
		int& nrecs,
		fT& topological_number_v2)
	{

		fT sumw = 0., Smax = this->minslope, dw0max= 0 , sumSdw = 0.;;

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
			fT tdw = this->connector->get_traverse_dx_from_links_idx(lix);
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

		fT sumf = 0.;
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

	fT weights_automator_v2(std::array<int,8>& receivers, 
		std::array<fT,8>& weights, std::array<fT,8>& slopes, 
		int& node, 
		int& nrecs,
		fT& topological_number_v2,
		fT& dw0max,
		int& recmax,
		fT& dx
		)
	{

		// Placeholders
		fT sumw = 0., Smax = this->minslope, sumSdw = 0.;;
		std::pair<fT,fT> fvech;

		for(int i = 0; i < nrecs; ++i)
		{
			int lix = receivers[i];

			
			if(this->connector->is_link_valid(lix) == false)
			{
				continue;
			}

			// int_dw += this->connector->get_dx_from_links_idx(lix);
			fT tdw = this->connector->get_traverse_dx_from_links_idx(lix);
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

		fT sumf = 0.;
		for(int i = 0; i<nrecs;++i)
		{
			int lix = receivers[i];

			if(this->connector->is_link_valid(lix) == false)
				continue;

			if(sumw > 0)
				weights[i] = weights[i]/sumw;
			else
				weights[i] = 1/nrecs;

			sumf += weights[i];
		}

		if(this->record_flowvec)
		{

			fT length = std::sqrt(std::pow(this->_rec_flowvec[node * 2],2) + std::pow(this->_rec_flowvec[node * 2 + 1],2));
			if(length!=0)
			{
				this->_rec_flowvec[node * 2] /= length;
				this->_rec_flowvec[node * 2 + 1] /= length;
			}
		}


		topological_number_v2 = (Smax*dw0max)/sumSdw;

		return Smax;
	}


	void _compute_morpho(int& node, int& recmax, fT& dx, fT& Smax, std::vector<fT>& vmot)
	{
		if(this->morphomode != MORPHO::NONE && this->connector->boundaries.forcing_io(node) == false)
		{

			// precaching rates
			fT edot = 0., ddot = 0., eldot_A = 0., eldot_B = 0., dldot_A = 0., dldot_B = 0.;
			
			// Lateral dx for lat e/d
			fT dy = this->connector->get_travers_dy_from_dx(dx);
			
			// And gathering the orthogonal nodes
			std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,recmax);

			// Calculating sheer stress
			fT tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Smax;

			if(tau>tau_max)
				this->tau_max = tau;

			// if(tau>150)
			// 	tau = 150;
			
			// Double checking the orthogonal nodes and if needs be to process them (basically if flow outs model or has boundary exceptions)
			int oA = orthonodes.first;
			if(this->connector->boundaries.forcing_io(oA) || this->connector->is_in_bound(oA) == false ||
			  this->connector->boundaries.no_data(oA) || this->connector->flow_out_model(oA))
				oA = -1;
			int oB = orthonodes.second;
			if(this->connector->boundaries.forcing_io(oB) || this->connector->is_in_bound(oB) == false ||
			  this->connector->boundaries.no_data(oB) || this->connector->flow_out_model(oB))
				oB = -1;

			// Calculating erosion if sheer stress is above critical
			if( tau > this->tau_c(node))
			{
				// Basal erosion rates
				edot = this->ke(node) * std::pow(tau - this->tau_c(node),this->aexp(node));
				// if recording stuff
				if(this->record_edot)
					this->_rec_edot[node] += edot;
			}

			// Claculating the deposition rates
			ddot = this->_Qs[node]/this->kd(node);

			// Dealing with lateral deposition if lS > 0 and erosion if lS <0
			if(oA >= 0 )
			{
				fT tSwl = this->get_Stopo(node,oA, dy);
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
				fT tSwl = this->get_Stopo(node,oB, dy);
				if(tSwl > 0)
				{
					dldot_B = tSwl * this->kd_lateral(node) * ddot;
				}
				else
				{
					eldot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
				}
			}

			// Am I depositing more than I can chew?
			fT totdqs = dx * (dldot_B + dldot_A + ddot);
			if(totdqs > this->_Qs[node])
			{
				// std::cout << " happens??? " << totdqs;
				fT corrqs = this->_Qs[node]/totdqs;
				dldot_B *= corrqs;
				dldot_A *= corrqs;
				ddot *= corrqs;
			}
			fT sqs = this->_Qs[node];
			fT fbatch = (ddot + dldot_B + dldot_A - edot - eldot_A - eldot_B) * dx;
			this->_Qs[node] -= fbatch;
			if(std::isfinite(this->_Qs[node]) == false)
			{
				std::cout << "QS NAN:" << this->_Qs[node] << " vs " << sqs << std::endl;
				throw std::runtime_error("BITE");
			}

			if(this->_Qs[node] < 0) 
			{
				this->_Qs[node] = 0;
			}

			vmot[node] += ddot - edot;
			if(oA>=0)
			{
				vmot[oA] += dldot_A;
				vmot[oA] -= eldot_A;
			}
			if(oB>=0)
			{
				vmot[oB] += dldot_B;
				vmot[oB] -= eldot_B;	
			}

			if(std::isfinite(vmot[node]) == false)
			{
				std::cout << "edot:" << edot << " ddot" << ddot << std::endl;
				std::cout << "qs:" << sqs << " tau" << tau << std::endl;
				throw std::runtime_error("Non finite vmot");
			}
		}

	}

	void _compute_transfers(int& nrecs, int& recmax, int& node, bool& SF, std::array<int,8>& receivers, std::array<fT,8>& weights, fT& Qwin, fT& Qwout)
	{
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
				// this->_Qw[rec] += std::min(Qwout * weights[j], Qwin * weights[j]);
				this->_Qw[rec] += Qwout * weights[j];
				// std::cout << "transferring:" << Qwout << std::endl;
			}

			if(this->morphomode != MORPHO::NONE)
			{
				if(std::isfinite(this->_Qs[node]) == false)
					throw std::runtime_error("QS NAN");
				if(std::isfinite(this->_Qs[rec]) == false)
					throw std::runtime_error("QSREC NAN");
				this->_Qs[rec] += (SF == false)?weights[j] * this->_Qs[node]:this->_Qs[node] ;
				if(std::isfinite(this->_Qs[rec]) == false)
				{
					std::cout << weights[j] << std::endl;;
					throw std::runtime_error("QSREC NAN AFTER");
				}
			}

		}
	}

	void _compute_transfers_exp(int& nrecs, int& recmax, int& node, bool& SF, std::array<int,8>& receivers, std::array<fT,8>& weights, fT& Qwin)
	{
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


			// if(this->hydrostationary)
			// {
			if(SF)
			{
				this->_Qw[rec] += Qwin; 
			}
			else if (weights[j] > 0 && Qwin > 0)
			{						
				this->_Qw[rec] += weights[j] * Qwin;
			}
			
			// if(this->morphomode != MORPHO::NONE)
			// {
			// 	if(std::isfinite(this->_Qs[node]) == false)
			// 		throw std::runtime_error("QS NAN");
			// 	if(std::isfinite(this->_Qs[rec]) == false)
			// 		throw std::runtime_error("QSREC NAN");
			// 	this->_Qs[rec] += (SF == false)?weights[j] * this->_Qs[node]:this->_Qs[node] ;
			// 	if(std::isfinite(this->_Qs[rec]) == false)
			// 	{
			// 		std::cout << weights[j] << std::endl;;
			// 		throw std::runtime_error("QSREC NAN AFTER");
			// 	}
			// }

		}
	}

	void _compute_vertical_motions(std::vector<fT>& vmot_hw, std::vector<fT>& vmot, bool use_dt = true)
	{

		for(int i=0; i<this->graph->nnodes; ++i)
		{

			if(this->connector->flow_out_model(i) && this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if(this->connector->boundaries.forcing_io(i)) continue;
			
			fT tvh = vmot_hw[i];
			if(use_dt)
			{
				tvh *= this->dt_hydro(i);
				// std::cout << this->dt_hydro(i) << std::endl;
			}
			if(tvh < - this->_hw[i])
			{
				tvh = - this->_hw[i];
			}

			this->_hw[i] += tvh;
			if(this->record_dhw)
				this->_rec_dhw[i] = tvh;
			
			this->_surface[i] += tvh;

			if(this->morphomode != MORPHO::NONE)
				this->_surface[i] += vmot[i] * this->dt_morpho(i);

			if(this->_hw[i] < 0)
				throw std::runtime_error("hw < 0???");
		}

	}

	void _compute_vertical_motions_averaged_test(std::vector<fT>& vmot_hw, std::vector<fT>& vmot, bool use_dt = true)
	{

		for(int i=0; i<this->graph->nnodes; ++i)
		{

			if(this->connector->flow_out_model(i) && this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if(this->connector->boundaries.forcing_io(i)) continue;
			
			fT tvh = vmot_hw[i];
			if(use_dt)
			{
				tvh *= this->dt_hydro(i);
				// std::cout << this->dt_hydro(i) << std::endl;
			}
			if(tvh < - this->_hw[i])
			{
				tvh = - this->_hw[i];
			}

			fT newhw = this->_hw[i] + tvh;
			fT prevhw = this->_hw[i];
			this->_hw[i] = (this->_hw[i] + newhw)/2;

			if(this->record_dhw)
				this->_rec_dhw[i] = tvh;
			
			this->_surface[i] += newhw - prevhw;

			if(this->morphomode != MORPHO::NONE)
				this->_surface[i] += vmot[i] * this->dt_morpho(i);

			if(this->_hw[i] < 0)
				throw std::runtime_error("hw < 0???");
		}

	}



	// small helper function returning node index from the right stack
	int get_istack_node(int i)
	{
		if(this->hydromode != HYDRO::GRAPH_SFD)
			return this->graph->stack[i];
		else 
			return this->graph->Sstack[i];
	}

	fT get_Sw(int lix, fT minslope)
	{
		int from,to; this->connector->from_to_from_link_index(lix, from, to);
		return std::max((this->_surface[from] - this->_surface[to])/this->connector->get_dx_from_links_idx(lix),minslope); 
	}

	fT get_Sw(int node, int rec, fT dx, fT minslope)
	{
		return std::max((this->_surface[node] - this->_surface[rec])/dx,minslope); 
	}

	fT get_Sw(int node, int rec, fT dx)
	{
		return (this->_surface[node] - this->_surface[rec])/dx; 
	}

	fT get_Stopo(int node, int rec, fT dx)
	{
		return (this->_surface[node] - this->_hw[node] - (this->_surface[rec]  - this->_hw[rec]))/dx; 
	}



	// GET DATA OUT
	template<class out_t>
	out_t get_hw(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_hw) ;}
	
	template<class out_t>
	out_t get_surface_topo(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_surface) ;}

	template<class out_t>
	out_t get_bedrock_topo()
	{
		std::vector<fT> diff(this->_surface);
		
		for(int i=0; i< this->graph->nnodes; ++i)
			diff[i] -= this->_hw[i];

		return DAGGER::format_output<std::vector<fT>, out_t >(diff) ;

	}

	template<class out_t>
	out_t get_Qwin(){ return DAGGER::format_output<std::vector<fT>, out_t >(this->_Qw) ;}

	template<class out_t>
	out_t get_SSTACKDEBUG(){ return DAGGER::format_output<std::vector<size_t>, out_t >(this->graph->Sstack) ;}



	void catch_nan(fT testval, std::string error_message)
	{
		if(std::isfinite(testval) == false)
			throw std::runtime_error(error_message);
	}



	// Run a basic version of floodos, for testing purposes
	// - N is the number of precipitons to launch
	// - precdt is the time step between 2 successive precipitons
	void run_precipitions(int N, fT precdt)
	{
		// check if this is the first lauch and allocate the vectors if needed	
		if(this->last_dt_prec.size() == 0)
		{
			this->last_dt_prec = std::vector<fT>(this->connector->nnodes,0.);
			this->last_dt_prec_e = std::vector<fT>(this->connector->nnodes,0.);
			this->last_sw_prec = std::vector<fT>(this->connector->nnodes,0.);
			this->last_dx_prec = std::vector<fT>(this->connector->nnodes,0.);
		}


		// Computes the total volume of Qin.dt for all entry points
		// -> if precipitation: sums all the precipitations * cellarea * dt
		// -> if entry Qwin: sum all the entryQwin*dt
		// this also sets the stochastic spawner
		this->compute_Vp(precdt);

		// current_dt_prec
		// caching neighbours and slopes
		auto neighbours = this->connector->get_empty_neighbour();

		// launching n precipitons
		for(int _ = 0; _<N; ++_)
		{

			// Spawn location
			int node = this->spawn_precipition();
			// caching the starting point for debugging
			// int onode = node;

			// increment the global chronometer
			this->current_dt_prec += precdt;

			// and for erosion if activated
			if(this->morphomode != MORPHO::NONE)
				this->current_dt_prec_e += precdt;

			// sediment flux
			fT tQs = this->Vps * this->connector->dy;

			// adding a safeguard stopping the precipiton if it runs for too long (mostly used to avoid infinite loop when divergence occurs)
			int NMAX = this->connector->nnodes * 2;
			while(true)
			{
				--NMAX;

				// stops the precipitons if it has reached model edge
				if(this->connector->boundaries.has_to_out(node) || NMAX == 0)
					break;


				// local time step
				fT tdt = this->current_dt_prec - this->last_dt_prec[node];
				// this->last_dt_prec[node] = this->current_dt_prec;

				// caching receivers/slope/...
				int rec = node;// will be rec
				fT Sw = 0.; // will be hydraulic slope
				fT Sw_sel = 0.; // stochastic slope to select next receiver
				// spatial steps
				fT dx = 0.;
				fT dy = 0.;

				// getting the current neighbours
				int nn = this->connector->get_neighbour_idx_links(node,neighbours);
				for(int j=0; j<nn; ++j)
				{
					// local neighbour
					int trec = this->connector->get_other_node_from_links(neighbours[j], node);
					// ignore no data
					if(this->connector->boundaries.no_data(trec))
						continue;

					// if(this->last_dt_prec[trec] != this->current_dt_prec)
					// decrease the water height of neighbour based on previous sw
					this->update_hw_prec_and_dt(trec);
				}

				// decrease node's hw based on previous sw
				this->update_hw_prec_and_dt(node);


				// Takes care of adding the precipiton's hw (making sure water is not negative after the decrease too)
				fT deltah = 0.;
				if(this->connector->boundaries.can_out(node) == false)
				{
					deltah = this->Vp;
					this->_hw[node] += deltah;
					if(this->_hw[node] < 0.) 
					{
						deltah += this->_hw[node];
						this->_hw[node] = 0.;
					}
				}

				// propagating to the surface
				this->_surface[node] += deltah;
				

				// Determining the next path and updating the slopw/time/...
				for(int j=0; j<nn; ++j)
				{

					// testing next receiver
					int trec = this->connector->get_other_node_from_links(neighbours[j], node);

					if(this->connector->boundaries.no_data(trec))
					{
						continue;
					}

					// if potential receiver
					if(this->_surface[trec] <= this->_surface[node])
					{
						
						// calculating slope
						fT tsw = (this->_surface[node] - this->_surface[trec])/this->connector->get_dx_from_links_idx(neighbours[j]);
						// stochasticity f(sqrt(tsw))
						fT tsw_sel = std::sqrt(tsw) * this->randu.get() *this->connector->get_traverse_dx_from_links_idx(neighbours[j]);
						
						// Selecting these info if selected
						if(tsw_sel > Sw_sel)
						{
							Sw_sel = tsw_sel;
							rec = trec;
							dx = this->connector->get_dx_from_links_idx(neighbours[j]);
							dy = this->connector->get_travers_dy_from_dx(dx);
						}

						// finding Smax
						if(tsw > Sw)
							Sw = tsw;

						// managing preboundary conditions
						if(this->connector->boundaries.can_out(trec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
						{
							Sw_sel = this->bou_fixed_val;
							Sw = this->bou_fixed_val;
							rec = trec;
							dx = this->connector->dx;
							dy = this->connector->dy;
							break;
						}
					}

				}


				// I have the Smax and other thingies, saving them
				this->last_sw_prec[node] = Sw;
				this->last_dx_prec[node] = dx;


				// Dealing with morpho if enabled
				if( this->morphomode != MORPHO::NONE && this->connector->boundaries.forcing_io(node) == false)
				{
					if(Sw>0 && rec != node && dx > 0)
					{
						tdt = (this->current_dt_prec_e - this->last_dt_prec_e[node]) * this->dt_morpho_multiplier;
						this->last_dt_prec_e[node] = this->current_dt_prec_e;
						fT tau = this->rho(node) * this->_hw[node] * this->GRAVITY * Sw;
						fT edot = 0., ddot = 0., edot_A = 0., edot_B = 0, ddot_A = 0, ddot_B = 0.;
						if(tau > this->tau_c(node))
							edot += this->ke(node) * std::pow(tau - this->tau_c(node),this->aexp(node));
						
						ddot = tQs/this->kd(node);

						if(tQs <0)
						{
							// std::cout << "happens"; 
							tQs = 0;
						}


						std::pair<int,int> orthonodes = this->connector->get_orthogonal_nodes(node,rec);
						int oA = orthonodes.first;
						if(this->connector->boundaries.forcing_io(oA) || this->connector->is_in_bound(oA) == false ||  this->connector->boundaries.no_data(oA) || this->connector->flow_out_model(oA))
							oA = -1;
						int oB = orthonodes.second;
						if(this->connector->boundaries.forcing_io(oB) || this->connector->is_in_bound(oB) == false ||  this->connector->boundaries.no_data(oB) || this->connector->flow_out_model(oB))
							oB = -1;

						// Dealing with lateral deposition if lS > 0 and erosion if lS <0
						if(oA >= 0 )
						{
							fT tSwl = this->get_Stopo(node,oA, dy);
							if(tSwl > 0)
							{
								ddot_A = tSwl * this->kd_lateral(node) * ddot;
							}
							else
							{
								edot_A = std::abs(tSwl) * this->ke_lateral(node) * edot;
							}

							// std::cout << edot_A << "|" << ddot_A  << std::endl;

						}

						if(oB >= 0 )
						{
							fT tSwl = this->get_Stopo(node,oB, dy);
							if(tSwl > 0)
							{
								ddot_B = tSwl * this->kd_lateral(node) * ddot;
							}
							else
							{
								edot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
							}
						}

						fT totd = ddot + ddot_B + ddot_A;
						totd *= dx;
						if(totd>tQs)
						{
							auto cor = tQs/totd;
							ddot*= cor;
							ddot_B*= cor;
							ddot_A*= cor;
						}

						tQs += (edot - ddot) * dx;
						tQs += (edot_B - ddot_B) * dy;
						tQs += (edot_A - ddot_A) * dy;
						if(oA>=0)
							this->_surface[oA] += (ddot_A -  edot_A) * tdt;
						if(oB >=0)
							this->_surface[oB] += (ddot_B -  edot_B) * tdt;
	
						this->_surface[node] += (ddot - edot) * tdt;

					}
					// else
					// {
					// 	this->_surface[node] += this->Vp;
					// }

				}

				// stopping the loop if node is leaving the model				
				if(node == rec && this->connector->boundaries.can_out(node))
					break;

				// updating positions
				node = rec;

			} // end of while loop for a single precipiton

		} // end of the for loop for all the launches


	} // end of the precipiton function

	std::vector<fT> get_precipitations_vector()
	{
		std::vector<fT> out(this->connector->nnodes,0.);
		for(int i=0; i<this->connector->nnodes; ++i)
		{
			out[i] = this->precipitations(i);
		}
		return out;
	}



	void define_precipitations_Ath(fT Ath)
	{
		std::vector<BC> nBCs(this->connector->boundaries.codes);
		std::vector<fT> fake(this->_surface);
		this->graph->_compute_graph(fake,false,false);
		auto DA = this->graph->_accumulate_constant_downstream_SFD(this->connector->get_area_at_node(0));
		std::vector<std::uint8_t> vis(this->connector->nnodes,false);
		std::vector<fT> tprecs = this->get_precipitations_vector();

		for(int i = this->connector->nnodes-1; i>= 0; --i)
		{

			int node = this->graph->Sstack[i];
			
			if(this->connector->boundaries.no_data(node))
				continue;

			// if(vis[node]) continue;
			
			int rec = this->connector->Sreceivers[node];
			if(Ath<DA[node])
			{
				vis[node] = true;
			}
			else
			{
				tprecs[rec] += tprecs[node];
			}

			vis[rec] = vis[node];
		}

		
		std::cout << "LKDFJDLKF" << std::endl;

		auto receivers = this->connector->get_empty_neighbour();
		for(int i = this->connector->nnodes-1; i>= 0; --i)
		{
			int node = this->graph->stack[i];

			if(vis[node])
			{
				int nn = this->connector->get_receivers_idx(node,receivers);
				for(int j=0; j< nn; ++j)
					vis[receivers[j]] = true;
			}

		}

		for(int i = this->connector->nnodes-1; i>= 0; --i)
		{
			if(vis[i] == false)
				nBCs[i] = BC::NO_FLOW;
		}


		std::cout << "??" << std::endl;
		this->connector->set_custom_boundaries(nBCs);
		this->set_water_input_by_variable_precipitation_rate(tprecs);
		std::cout << "!!" << std::endl;

	}


	void run_graphipiton(int N, fT precdt, fT Ath)
	{
		if(this->last_dt_prec.size() == 0)
		{
			this->last_dt_prec = std::vector<fT>(this->connector->nnodes,0.);
			this->last_dt_prec_e = std::vector<fT>(this->connector->nnodes,0.);
			this->last_sw_prec = std::vector<fT>(this->connector->nnodes,0.);
			this->last_dx_prec = std::vector<fT>(this->connector->nnodes,0.);
		}


		std::vector<fT> fake(this->_surface);
		this->graph->_compute_graph(fake,true,false);
		auto DA = this->graph->_accumulate_constant_downstream_SFD(this->connector->get_area_at_node(0));
		std::vector<std::uint8_t> vis(this->connector->nnodes,false);
		std::vector<int> starters;starters.reserve(1000);
		for(int i = this->connector->nnodes-1; i>= 0; --i)
		{
			int node = this->graph->Sstack[i];
			if(vis[node]) continue;
			int rec = this->connector->Sreceivers[node];
			vis[node] = true;
			if(Ath<DA[node])
			{
				vis[rec] = true;
				starters.emplace_back(node);
			}

		}

		this->compute_Vp(precdt);

		this->dis = std::uniform_int_distribution<> (0, starters.size()-1);



		// current_dt_prec
		auto neighbours = this->connector->get_empty_neighbour();

		for(int _ = 0; _<N; ++_)
		{

			// int node = this->spawn_precipition();
			int node = starters[this->dis(this->gen)];
			// int onode = node;
			// std::cout << "init at " << node << std::endl;
			
			this->current_dt_prec += precdt;

			if(this->morphomode != MORPHO::NONE)
				this->current_dt_prec_e += precdt;

			// fT tQs = this->Vps * this->connector->dy;

			int NMAX = this->connector->nnodes * 2;
			while(true)
			{
				--NMAX;
				if(this->connector->boundaries.has_to_out(node) || NMAX == 0)
					break;

				// // update step from Qwout
				// fT tdt = this->current_dt_prec - this->last_dt_prec[node];
				// // this->last_dt_prec[node] = this->current_dt_prec;

				int rec = node;
				fT Sw = 0.;
				fT Sw_sel = 0.;
				fT dx = 0.;
				// fT dy = 0.;
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

				fT deltah = 0.;
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
					{
						continue;
					}

					if(this->_surface[trec] <= this->_surface[node])
					{
						
						fT tsw = (this->_surface[node] - this->_surface[trec])/this->connector->get_dx_from_links_idx(neighbours[j]);
						fT tsw_sel = tsw * this->randu.get();
						// std::cout << tsw_sel << std::endl;
						if(tsw_sel > Sw_sel)
						{
							Sw_sel = tsw_sel;
							rec = trec;
							dx = this->connector->get_dx_from_links_idx(neighbours[j]);
							// dy = this->connector->get_travers_dy_from_dx(dx);
						}

						if(tsw > Sw)
							Sw = tsw;

						if(this->connector->boundaries.can_out(trec) && this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
						{
							Sw_sel = this->bou_fixed_val;
							Sw = this->bou_fixed_val;
							rec = trec;
							dx = this->connector->dx;
							// dy = this->connector->dy;
						}
					}

				}



				this->last_sw_prec[node] = Sw;
				this->last_dx_prec[node] = dx;

				
				// std::cout << onode << "||" << node << " vs " << rec << "||";

				if(node == rec && this->connector->boundaries.can_out(node))
					break;

				// if(node != rec)
				// 	throw std::runtime_error("DONE!");
				node = rec;
			}
			// throw std::runtime_error("BITE2");

		}


	}

	void update_hw_prec_and_dt(int node)
	{

		fT tdt = this->current_dt_prec - this->last_dt_prec[node];
		this->last_dt_prec[node] = this->current_dt_prec;

		if(tdt == 0 || this->last_sw_prec[node] == 0 || this->last_dx_prec[node] == 0)
		{
			return;
		}
		
		auto delta = tdt * std::pow(this->_hw[node],5./3.) * std::sqrt(this->last_sw_prec[node])/this->mannings(node)/this->last_dx_prec[node];
		this->_hw[node] -= delta;
		this->_surface[node] -= delta; 
	}

	void update_hw_prec_and_dt_exp_2(int node, std::vector<fT>& gHwin, fT tHwin, fT precdt)
	{

		fT tdt = this->current_dt_prec - this->last_dt_prec[node];
		this->last_dt_prec[node] = this->current_dt_prec;

		// SHAKING THINGS HERE
		tdt = precdt;

		if(tdt == 0 || this->last_sw_prec[node] == 0 || this->last_dx_prec[node] == 0)
		{
			// fT dHwin =  std::max(tHwin, gHwin[node]);

			// this->_hw[node] += dtfill * dHwin ;
			// this->_surface[node] += dtfill * dHwin ;
			return;
		}

		fT dHwin =  std::max(tHwin, gHwin[node]);
		// std::cout << "+:" << dHwin * tdt;
		
		auto delta = tdt * std::pow(this->_hw[node],5./3.) * std::sqrt(this->last_sw_prec[node])/this->mannings(node)/this->last_dx_prec[node];
		// std::cout << " vs - :" << delta << std::endl;;
		delta -= dHwin * tdt;
		if(delta < - this->_hw[node])
			delta = - this->_hw[node];

		this->_hw[node] -= delta;
		this->_surface[node] -= delta;
		// if(delta > 0)
		// 	std::cout <<"delta::" << delta << std::endl;
	}

	void update_hw_prec_and_dt_exp(int node, std::vector<int>& baslab, std::vector<fT>& basdt)
	{

		fT tdt = basdt[baslab[node]] - this->last_dt_prec[node];
		this->last_dt_prec[node] =  basdt[baslab[node]];

		if(tdt == 0 || this->last_sw_prec[node] == 0 || this->last_dx_prec[node] == 0)
		{
			return;
		}
		
		auto delta = tdt * std::pow(this->_hw[node],5./3.) * std::sqrt(this->last_sw_prec[node])/this->mannings(node)/this->last_dx_prec[node];
		this->_hw[node] -= delta;
		this->_surface[node] -= delta; 
	}

	void compute_Vp(fT precdt)
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
				this->Vps+=this->_sed_entries[i] * this->connector->dy * precdt * this->dt_morpho_multiplier;
			}
			// this->dis = std::uniform_int_distribution<> (0,this->_water_entries.size()-1);
		
		}


	}


	void run_precipitions_exp(int N, fT precdt, fT dtfill)
	{
		// place holder for debugging
		return;
	}


	void run_precipitions_exp_1(int N, fT precdt)
	{
		// place holder for debugging
		return;
	}




















	// TOPO ANALYSIS


	// computes u, q or Q
	// dim: 1,2 or 3
	template<class out_t>
	out_t compute_tuqQ(int dimensions)
	{
		if(dimensions < 0 || dimensions > 3)
			throw std::runtime_error("Invalid number of dimensions. Needs to be 0 for the topological number (dfimensionless coefficient of paritionning), 1 (velocity in m/s), 2 (discharge per unit width in m^2/s) or 3 (volumetric discarghe in m^3/s)");

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition and other monitoring features:
		this->init_Qw();

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> out(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);
			
			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers)) continue;

			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			//NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);


			// Metric to calculate:
			fT metric = topological_number_v2;
			if(dimensions > 0)
				metric = pohw * sqrtSmax/this->mannings(node);
			if(dimensions > 1)
				metric *= this->_hw[node];
			if(dimensions > 2)
				metric *= dw0max;

			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(nrecs, recmax, node, SF, receivers, weights, Qwin, metric);


			out[node] = metric;


			// Volumetric discahrge
			// fT Qwout = dw0max * this->_hw[node] * u_flow;

		}

	
		return format_output<decltype(out),out_t>(out);
	}

	
	template<class in_t, class out_t>
	out_t compute_elemental_transfer(in_t& _production, fT total_time, fT concentration_max)
	{


		auto production = format_input(_production);

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition and other monitoring features:
		this->init_Qw();

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> out(this->graph->nnodes,0.), remaining_t(this->graph->nnodes,0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT,8> weights, slopes;
		
		// main loop
		for(int i = this->graph->nnodes-1; i>=0; --i)
		{
			
			// Getting next node in line
			int node = this->get_istack_node(i);

			production[node]  = std::min(production[node], this->_Qw[node] * concentration_max);

			if(production[node] > 0)
			{
				remaining_t[node] = (remaining_t[node] * out[node] + production[node] *total_time)/(out[node] + production[node]);
				out[node] +=  production[node]/this->_hw[node];
			}

			out[node] = std::min(out[node], this->_Qw[node] * concentration_max);

			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need to be processed
			// THis is where all the boundary treatment happens, if you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers)) continue;

			// Getting the receivers
			int nrecs; 
			if(SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node,receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			//NOTE:
			// No need to calculate the topological number anymore, but keeping it for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al( node, SF, Smax, slopes, weights, nrecs, receivers, recmax, dx, dw0max, topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], this->TWOTHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);


			// Metric to calculate:
			fT transfer_time = dx/(pohw * sqrtSmax/this->mannings(node));
			remaining_t[node] -= transfer_time;

			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(nrecs, recmax, node,  SF, receivers, weights, Qwin, transfer_time);

			if(out[node] > 0 && remaining_t[node] > 0)
			{
				for (int j=0; j<nrecs; ++j)
				{
					int trec = this->connector->get_to_links(receivers[j]);
					fT tadd = out[node] * weights[j];

					remaining_t[trec] = (remaining_t[trec] * out[trec] + remaining_t[node] * tadd )/(tadd + out[trec]);
					
					out[trec] = (tadd + out[trec]);

					out[trec] = std::min(out[trec], this->_Qw[trec] * concentration_max);
				}
			}

			// out[node] = metric;


			// Volumetric discahrge
			// fT Qwout = dw0max * this->_hw[node] * u_flow;

		}


		for(int i=0; i<this->connector->nnodes; ++i)
			out[i] *= remaining_t[i];
	
		return format_output<decltype(out),out_t>(out);
	}



	// Experimental drainage area stuff
	// template<class out_t








	int spawn_precipition()
	{

		int node;

		do
		{
			bool OK = false;
			
			if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT || this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
			{
				node =  this->dis(this->gen);
				if(this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT)
					OK = true;
				else if (this->precipitations(node) > 0)
					OK = true;
			}
			else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H)
			{
				node =  this->_water_entry_nodes[this->dis(this->gen)];
				OK = true;
			}
			else
				throw std::runtime_error("NOT DEFINED");

			if(this->connector->boundaries.no_data(node) == false && OK) break;

		} while(true);
		// std::cout << node << std::endl;

		return node;

	}








};
// end of graphflood class






} // End of DAGGER namespace

































#endif