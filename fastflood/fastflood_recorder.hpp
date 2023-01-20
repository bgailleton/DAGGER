#ifndef FASTFLOOD_RECORDER_HPP
#define FASTFLOOD_RECORDER_HPP

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

#ifdef OPENMP_YOLO
#include <omp.h>
#endif

// local includes 
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D8connector.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

namespace fastflood
{

	template<class float_t>
	class recorder
	{
	public:

		recorder(){;}
		recorder(int nnodes, int nlinks){this->nnodes = nnodes; this->nlinks = nlinks;}

		// nnodes for the recording
		int nnodes = 0;
		int nlinks = 0;

		// stuff to record, eventually
		bool edot2record = false;
		std::vector<float_t> edot;
		void enable_edot_recording(){this->edot2record = true; this->edot = std::vector<float_t>(this->nnodes,0.);};
		void disable_edot_recording(){this->edot2record = false; this->edot = std::vector<float_t>();}
		template<class out_t>
		out_t get_edot(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->edot) ;}

		bool ddot2record = false;
		std::vector<float_t> ddot;
		void enable_ddot_recording(){this->ddot2record = true; this->ddot = std::vector<float_t>(this->nnodes,0.);};
		void disable_ddot_recording(){this->ddot2record = false; this->ddot = std::vector<float_t>();}
		template<class out_t>
		out_t get_ddot(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->ddot) ;}

		bool lateral_edot2record = false;
		std::vector<float_t> lateral_edot;
		void enable_lateral_edot_recording(){this->lateral_edot2record = true; this->lateral_edot = std::vector<float_t>(this->nnodes,0.);};
		void disable_lateral_edot_recording(){this->lateral_edot2record = false; this->lateral_edot = std::vector<float_t>();}
		template<class out_t>
		out_t get_lateral_edot(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->lateral_edot) ;}

		bool lateral_ddot2record = false;
		std::vector<float_t> lateral_ddot;
		void enable_lateral_ddot_recording(){this->lateral_ddot2record = true; this->lateral_ddot = std::vector<float_t>(this->nnodes,0.);};
		void disable_lateral_ddot_recording(){this->lateral_ddot2record = false; this->lateral_ddot = std::vector<float_t>();}
		template<class out_t>
		out_t get_lateral_ddot(){ return DAGGER::format_output<std::vector<float_t>, out_t >(this->lateral_ddot) ;}




	};


}


#endif