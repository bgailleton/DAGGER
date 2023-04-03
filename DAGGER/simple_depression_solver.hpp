/*
This header file extends the graph to provide routines to process depressions using Cordonnier et al, 2019
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef SIMPLE_DEP_SOLVER_HPP
#define SIMPLE_DEP_SOLVER_HPP

// STL imports
#include <limits>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <fstream>
#include <queue>
#include <iostream>
#include <numeric>
#include <cmath>
#include <initializer_list>
#include <chrono>
#include <unordered_map>

#include "utils.hpp"
// #include "graph.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


namespace DAGGER
{


	template<class Connector_t, class float_t>
	bool simple_depression_solver(Connector_t* connector, std::vector<float_t>& topography, std::vector<size_t>& Sstack)
	{


		// STEP 1 LABEL ALL BASINS and fetch the one giving out
		std::vector<int> basin_key(this->connector->nnodes, -1);
		int lab = 0;
		for(auto node : Sstack)
		{
			int rec = connector->Sreceivers[node];
			if(connector->flow_out_model(node))
			{
				basin_key[node] = 0;
			}
			else if (node == rec)
			{
				++lab;
				basin_key[node] = lab;
			}
			else
			{
				basin_key[node] = basin_key[rec];
			}
		}

		// if no outside basins, returning false (i.e. does not need reprocessing)
		if(lab == 0)
			return false;

		int n_basins = lab+1;
		std::vector<int> connode(n_basins, -1);
		std::vector<std::uint8_t> isopened(n_basins, 0);
		isopened[0] = 1;
		std::vector<float_t> connode_z(n_basins, std::numeric_limits<float_t>::max());

		while( lab > 0 )
		{
			for(int i = 0; i < int(connector->links.size()) ++i)
			{
				if(connector->is_link_valid(i) == false)
					continue;
				int a,b; node_idx_from_link_idx_nocheck(i,a,b);
				int ba = basin_key[a], bb = basin_key[b];
				if(isopened[ba] != isopened[bb])
				{
					int tclosed = (isopened[ba] == 1) ? bb : ba;
					float_t tmitopo = std::min(topography[a],topography[b]);
					if(connode_z[tclosed] > tmitopo)
					{
						connode//...
					}
				}
			}
		}


	}




} // end of namespace










































#endif