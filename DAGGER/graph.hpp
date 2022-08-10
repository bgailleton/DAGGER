//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef graph_HPP
#define graph_HPP


/*
This file contains the graph class.
The graph manages anything linked to the DAG and is the main interface to the code.
B.G. 2022
*/



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

// defines all the format_input depnding on the eventual wrapper
#ifdef DAGGER_FT_PYTHON
#include "wrap_helper_python.hpp"
#else
#include "wrap_helper_cpp.hpp"
#endif


namespace DAGGER
{


template<class float_t>
class graph
{

// Everything goes public, more straighforward
public:

	// Number of nodes in the graph
	int nnodes;
	
	// Number of neighbours by nodes
	int n_neighbours;
	
	// bool vector for each link: true is receiver direction and false is donor
	// The meaning of the index depends on the connector
	std::vector<bool> links;
	
	// integer vector of 2*links size with the node indices of each link
	// for example, the nodes of link #42 would be indices 84 and 85
	std::vector<int> linknodes;
	
	// Single graph receivers
	// -> Sreceivers: steepest recervers (nnodes size), 
	// -> number of donors (nnodes size),
	// -> Steepest donors (nnodes * n_neighbours size)
	// --> Sdonors of node i are located from index i*n_neighbours to index i*n_neighbours + nSdonors[i] not included
	std::vector<int> Sreceivers,nSdonors,Sdonors;

	// Single graph distance to receivers
	std::vector<float_t> Sdistance2receivers;

	// Steepest slope
	std::vector<float_t> SS;
	
	// Topological order and Single graph topological order
	// While the MD stack can be used for single flow topology, 
	// the SD graph sensu Braun and Willett 2013 is built in a very comprehensive way 
	// making operations such as watershed labelling or connected component gathering particularly efficient
	std::vector<size_t> stack, Sstack;

	
	// default empty constructor
	graph(){};

	// Classic constructor, simply giving the number of nodes and the number of neighbours per node
	graph(int nnodes, int n_neighbours){this->nnodes = nnodes; this->n_neighbours = n_neighbours;}


	// Initialisation of the graph structure
	// Uses the nnodes, n_neighbours and connector to alloate the single vectors and the irec/linknodes attribute
	template<class Connector_t>
	void init_graph(Connector_t& connector)
	{
		// Allocate vectors to node size
		this->_allocate_vectors();
		// Calling the cionnector to allocate the linknodes to the right size
		connector.fill_linknodes(this->linknodes);
	}

	// used to reinitialise the the vectors
	template<class Connector_t>
	void reinit_graph(Connector_t& connector)
	{
		// reinitialise the vector without reallocating the full memory
		this->_reallocate_vectors();
	}

	// Compute the graph using the cordonnier metod to solve the depressions
	// template arguments are the connector type, the wrapper input type for topography and the wrapper output type for topography
	template<class Connector_t,class topo_t, class out_t>
	out_t compute_graph(
		std::string depression_solver, // String switching the type of depression solver: "cordonnier_carve", "cordonnier_fill", "cordonnier_simple" or "priority_flood"
	  topo_t& ttopography, // the input topography
	  Connector_t& connector, // the input connector to use (e.g. D8connector)
	  bool only_SD, // only computes the single flow graph if true
	  bool quicksort // computes the MF toposort with a quicksort algo if true, else uses a BFS-based algorithm (which one is better depends on many things)
	  )
	{
		// Formatting the input to match all the wrappers
		auto topography = format_input(ttopography);

		// Checking if the depression method is cordonnier or node
		bool isCordonnier = this->is_method_cordonnier(depression_solver);

		// Formatting the output topo
		std::vector<float_t> faketopo(to_vec(topography));

		// if the method is not Cordonnier -> apply the other first
		if(isCordonnier == false)
		{
			// filling the topography with a minimal slope using Wei et al., 2018
			faketopo = connector.PriorityFlood_Wei2018(topography);
		}

		// Making sure the graph is not inheriting previous values
		this->reinit_graph(connector);

		// Updates the links vector and the Srecs vector by checking each link new elevation
		this->update_recs(faketopo, connector);

		// Compute the topological sorting for single stack
		// Braun and willett 2014 (modified)
		this->topological_sorting_SF();


		// manages the Cordonnier method if needed
		if(isCordonnier)
		{
			
			// LMRerouter is the class managing the different cordonnier's mthod
			LMRerouter<float_t> depsolver;
			// Execute the local minima solving, return true if rerouting was necessary, meaning that some element needs to be recomputed
			// Note that faketopo are modified in place.
			bool need_recompute = depsolver.run(depression_solver, faketopo, connector, this->Sreceivers, this->Sdistance2receivers, this->Sstack, this->linknodes);		

			// Right, if reomputed needs to be
			if(need_recompute)
			{
				
				// Re-inverting the Sreceivers into Sdonors
				this->recompute_SF_donors_from_receivers();
				
				// Recomputing Braun and willett 2014 (modified)
				this->topological_sorting_SF();

				// This is a bit confusing and needs to be changed but filling in done in one go while carving needs a second step here
				if(depression_solver == "carve")
					this->carve_topo_v2(1e-5, connector, faketopo);

				// My work here is done if only SD is needed
				if(only_SD)
					return format_output(faketopo);
				
				// Otherwise, conducting the topological sorting
				if(quicksort)
					this->topological_sorting_quicksort(faketopo);
				else
					this->topological_sorting_dag(connector);

				// And updating the multiple flow receivers (! careful not to touch the Sreceivers which are conditionned to Cordonnier solver)
				this->update_Mrecs(faketopo,connector);
			}
		}

		// if there is no need to recompute neighbours, then I can only calculate the topological sorting
		// for multiple as the toposort for SF is already done
		else if (only_SD == false)
		{
			if(quicksort)
					this->topological_sorting_quicksort(faketopo);
				else
					this->topological_sorting_dag(connector);
		}


		// Finally I format the output topography to the right wrapper
		return format_output(faketopo);
	}



	// Function updating ONLY the MFD receivers
	// This is useful in the cases where SFD recs are conditionned by an other mean
	// and cannot be touched (e.g. Cordonnier)
	template<class Connector_t,class topo_t>
	void update_Mrecs(topo_t& topography, Connector_t& connector)
	{
		// iterating though every links
		for(size_t i = 0; i<this->links.size(); ++i)
		{
			// Getting hte 2 nodes of the current link
			int from = this->linknodes[i*2];
			int to = this->linknodes[i*2 + 1];
			
			// Checking the validity of the link
			if(connector.is_in_bound(from) == false || connector.is_in_bound(to) == false)
				continue;

			// by convention true -> topo1 > topo2
			if(topography[from] > topography[to])
				this->links[i] = true;
			else
				this->links[i] = false;
		}
		// done
	}

	// Updates all the link and the SFD info
	template<class Connector_t,class topo_t>
	void update_recs(topo_t& topography, Connector_t& connector)
	{
		// iterating through all the nodes
		for(size_t i = 0; i<this->links.size(); ++i)
		{

			// Getting ht etwo nodes of the links
			int from = this->linknodes[i*2];
			int to = this->linknodes[i*2 + 1];
			
			// Checking the validity of the link
			if(connector.is_in_bound(from) == false || connector.is_in_bound(to) == false)
			{
				continue;
			}

			// getting the link infos
			// -> dx
			float_t dx = connector.get_dx_from_links_idx(i);
			// -> slope
			float_t slope = (topography[from] - topography[to])/dx;

			// if slope is positive, to is the receiver by convention
			if(slope>0)
			{
				// Conventional direction
				this->links[i] = true;
				// if Steepest Slope is higher than the current recorded one
				if(this->SS[from]<slope)
				{
					// saving the Sreceivers info as temporary best choice
					this->Sreceivers[from] = to;
					this->Sdistance2receivers[from] = dx;
					this->SS[from] = slope;
				}
			}
			else
			{
				// Otherwise the convention is inverted:
				// isrec is falese and to is giving to from
				this->links[i] = false;
				// NOte that slope is absolute values
				slope = std::abs(slope);
				if(this->SS[to]<slope)
				{
					this->Sreceivers[to] = from;
					this->Sdistance2receivers[to] = dx;
					this->SS[to] = slope;
				}
			}

		}

		// Finally inverting the Sreceivers into the Sdonors info
		// Required for several routines
		this->compute_SF_donors_from_receivers();
	}



	// This is a debugging function checking the stack
	// You an ignore
	template<class out_t>
	out_t test_Srecs()
	{
		std::vector<int> OUT(this->nnodes,0);
		for(int i=0; i< this->nnodes;++i)
		{
			if(i !=  this->Sreceivers[i])
				++OUT[i];
		}

		return format_output(OUT);
	}








	// Helper functions to allocate and reallocate vectors when computing/recomputing the graph
	void _allocate_vectors()
	{
		this->links = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		this->linknodes = std::vector<int>(int(this->nnodes * this->n_neighbours), 0);
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		this->Sstack = std::vector<size_t>(this->nnodes,0);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		this->Sdistance2receivers = std::vector<float_t >(this->nnodes,-1);
		this->SS = std::vector<float_t>(this->nnodes,0.);
	}
	void _reallocate_vectors()
	{
		for(int i=0;i<this->nnodes; ++i)
		{
			this->Sreceivers[i] = i;
			this->Sdistance2receivers[i] = 0;
			this->SS[i] = 0;
		}
	}


	// Fucntion inverting the SFD receivers into donors
	void compute_SF_donors_from_receivers()
	{
		// Initialising the graph dimesions for the donors
		// All of thenm have the graph dimension
		this->Sdonors = std::vector<int>(this->nnodes * this->n_neighbours,-1);
		this->nSdonors = std::vector<int>(this->nnodes,0);

		// iterating through all the nodes
		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;

			// feeding the Sdonors array at rec position with current node and...
			this->Sdonors[trec * this->n_neighbours  + this->nSdonors[trec]] = i;
			// ... incrementing hte number of Sdonors
			this->nSdonors[trec] += 1;
		}
		// done
	}


	// Sma function than above but without reallocating hte memory (can save a bit of time depending on the context)
	void recompute_SF_donors_from_receivers()
	{

		for(int i=0; i < this->nnodes; ++i)
		{
			for(int j=0; j < this->n_neighbours; ++j)
				this->Sdonors[i * this->n_neighbours + j] = -1;
			this->nSdonors[i] = 0;
		}

		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;
			this->Sdonors[trec * this->n_neighbours  + this->nSdonors[trec]] = i;
			this->nSdonors[trec] += 1;
		}

	}



	// This is my implementation of Braun and willett 2014
	// Slightly modified:
	// - First it is based on the fastscapelib version, which bypass the need of the delta vectors and all by simply applying the recursion directly feeding the stack
	// - Secondly, recursion is not the best practice for most languages so I am using a stack data structure instead
	// Results are identical and performances similar (marginally better, but that is linked to c++ not being a heavy recursion friendly language)
	void topological_sorting_SF()
	{
		// The stack container helper
		std::stack<size_t, std::vector<size_t> > stackhelper;
		// std::vector<bool> isdone(this->nnodes,false);
		// going through all the nodes
		int istack = 0;
		for(int i=0; i<this->nnodes; ++i)
		{
			// if they are base level I include them in the stack
			if(this->Sreceivers[i] == i)
			{
				stackhelper.emplace(i);
				// ++istack;
			}

			// While I still have stuff in the stack helper
			while(stackhelper.empty() == false)
			{
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();stackhelper.pop();
				this->Sstack[istack] = nextnode;
				++istack;

				// as well as all its donors which will be processed next
				for( int j = 0; j < this->nSdonors[nextnode]; ++j)
				{
					stackhelper.emplace(this->Sdonors[nextnode*this->n_neighbours + j]);
				}

			}

		}
	}

	template< class Connector_t>
	void topological_sorting_dag(Connector_t& connector)
	{
		std::vector<int> nrecs(this->nnodes,0);
		std::queue<int> toproc;
		this->stack.clear();
		this->stack.reserve(this->nnodes);
		for(int i = 0; i < int( this->links.size() ) ; ++i)
		{
			
			if(connector.can_flow_even_go_there(i) == false)
			{
				this->stack.emplace_back(i);
				continue;
			}

			if(this->links[i])
				++nrecs[this->linknodes[i*2]];
			else
				++nrecs[this->linknodes[i*2+1]];
		}

		for(int i=0; i<this->nnodes; ++i)
		{
			if(nrecs[i] == 0 && connector.can_flow_even_go_there(i))
				toproc.emplace(i);
		}

		while(toproc.empty() == false)
		{
			int next = toproc.front();
			toproc.pop();
			this->stack.emplace_back(next);
			auto dons = this->get_donors_idx(next, connector);
			for(auto d : dons)
			{
				--nrecs[d];
				if(nrecs[d] == 0)
					toproc.emplace(d);
			}
		}

	}

	// Multiple flow topological sorting using quicksort algorithm
	// Topography is simply sorted by absolute elevation keeping track of the indices
	template<class topo_t>
	void topological_sorting_quicksort(topo_t& ttopography)
	{
		// Formatting hte topographic input from the wrapper
		auto topography = format_input(ttopography);
		// Dortng by index
		auto yolo = sort_indexes(topography);
		// the sorted index is now the stack
		this->stack = std::move(yolo);
	}

	

	/// this function enforces minimal slope 
	template<class Connector_t, class topo_t>
	std::vector<int> carve_topo_v2(float_t slope, Connector_t& connector, topo_t& topography)
	{

		std::cout << std::setprecision(8);
		std::vector<int> to_recompute;
		to_recompute.reserve(1000);
		for(int i=this->nnodes-1; i >= 0; --i)
		{
			int node  = this->Sstack[i];
			if(connector.can_flow_out_there(node) || connector.can_flow_even_go_there(node) == false)
				continue;
			int rec = this->Sreceivers[node];
			float_t dz = topography[node] - topography[rec];
			if(dz <= 0)
			{
				topography[rec] = topography[node] - slope + connector.randu.get() * 1e-7;// * d2rec;
				to_recompute.emplace_back(rec);
			}

		}
		return to_recompute;
	}
	/// this function enforces minimal slope 
	template<class Connector_t, class topo_t>
	std::vector<int> fill_topo_v2(float_t slope, Connector_t& connector, topo_t& topography)
	{
		std::vector<int> to_recompute;
		for(int i=0; i < this->nnodes; ++i)
		{
			int node  = this->Sstack[i];
			if(connector.can_flow_out_there(node) || connector.can_flow_even_go_there(node) == false)
				continue;

			int rec = this->Sreceivers[node];
			float_t dz = topography[node] - topography[rec];

			if(dz <= 0)
			{
				topography[node] = topography[rec] + slope + connector.randu.get() * 1e-6;// * d2rec;
				to_recompute.emplace_back(node);
			}
		}
		return to_recompute;
	}


	
	template<class Connector_t>
	std::vector<int> get_receivers_idx(int i, Connector_t& connector)
	{
		std::vector<int> recs; recs.reserve(this->n_neighbours);
		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		for(auto li:linkidx)
		{
			if(i == this->linknodes[li*2] && this->links[li])
				recs.emplace_back(this->linknodes[li*2 + 1]);
			else if(this->links[li] == false)
				recs.emplace_back(this->linknodes[li*2]);
		}
		return recs;
	}

	template<class Connector_t>
	std::vector<int> get_receivers_idx_links(int i, Connector_t& connector)
	{
		std::vector<int> recs; recs.reserve(this->n_neighbours);
		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		for(auto li:linkidx)
		{
			if(i == this->linknodes[li*2] && this->links[li])
				recs.emplace_back(li);
			else if(this->links[li] == false)
				recs.emplace_back(li);
		}
		return recs;
	}

	template<class Connector_t>
	std::vector<std::pair<int,int> > get_node_receivers_idx_pair(int i, Connector_t& connector)
	{
		std::vector<std::pair<int,int>> recs; recs.reserve(this->n_neighbours);

		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		for(auto li:linkidx)
		{
			if(i == this->linknodes[li*2] && this->links[li])
				recs.emplace_back(std::make_pair<int,int>(i, this->linknodes[li*2 + 1]));
			else if(this->links[li] == false)
				recs.emplace_back(std::make_pair<int,int>(i, this->linknodes[li*2]));
		}

		return recs;
	}

	template<class Connector_t>
	std::vector<int> get_donors_idx(int i, Connector_t& connector)
	{
		std::vector<int> dons; dons.reserve(this->n_neighbours);
		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		for(auto li:linkidx)
		{
			if(i == this->linknodes[li*2] && this->links[li] == false)
				dons.emplace_back(this->linknodes[li*2 + 1]);
			else if(this->links[li])
				dons.emplace_back(this->linknodes[li*2]);
		}
		return dons;
	}

	template<class Connector_t>
	std::vector<int> get_donors_idx_links(int i, Connector_t& connector)
	{
		std::vector<int> dons; dons.reserve(this->n_neighbours);
		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		for(auto li:linkidx)
		{
			if(i == this->linknodes[li*2] && this->links[li] == false)
				dons.emplace_back(li);
			else if(this->links[li])
				dons.emplace_back(li);
		}
		return dons;
	}

	template<class Connector_t>
	std::vector<std::pair<int,int> > get_node_donors_idx_pair(int i, Connector_t& connector)
	{
		std::vector<std::pair<int,int>> dons; dons.reserve(this->n_neighbours);

		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		for(auto li:linkidx)
		{
			if(i == this->linknodes[li*2] && this->links[li] == false)
				dons.emplace_back(std::make_pair<int,int>(i, this->linknodes[li*2 + 1]));
			else if(this->links[li])
				dons.emplace_back(std::make_pair<int,int>(i, this->linknodes[li*2]));
		}

		return dons;
	}







	template<class Connector_t, class topo_t, class out_t>
	out_t get_DA_proposlope(Connector_t& connector, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		std::vector<float_t> DA(connector.nnodes,0.);
		for(int i = connector.nnodes - 1; i>=0; --i)
		{
			int node = this->stack[i];
			DA[node] += connector.get_area_at_node(node);

			if(connector.is_active(node))
			{
				auto recs = this->get_receivers_idx(node, connector);

				std::vector<float_t> slopes(recs.size());
				float_t sumslopes = 0;
				for(size_t j = 0;j < recs.size(); ++j)
				{
					int rec = recs[j];
					slopes[j] = (topography[node] - topography[rec])/connector.dx;
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(size_t j = 0;j < recs.size(); ++j)
				{
					int rec = recs[j];
					DA[rec] += DA[node] * slopes[j]/sumslopes;
				}

			}

			// std::cout << std::endl;;

		}

		return format_output(DA);
	}


	template<class Connector_t, class topo_t, class out_t>
	out_t get_DA_SS(Connector_t& connector, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		std::vector<float_t> DA(connector.nnodes,0.);
		for(int i = connector.nnodes - 1; i>=0; --i)
		{
			int node = this->Sstack[i];
			DA[node] += connector.get_area_at_node(node);

			if(connector.can_flow_even_go_there(node) && node != this->Sreceivers[node])
			{
				int Srec = this->Sreceivers[node];
				DA[Srec] += DA[node];
			}
		}

		return format_output(DA);
	}


	template<class Connector_t>
	std::vector<int> get_rowcol_Sreceivers(int row, int col,  Connector_t& connector)
	{
		int node = connector.nodeid_from_row_col(row,col);
		std::vector<int> out_receivers;
		int trow,tcol;
		connector.rowcol_from_node_id(this->Sreceivers[node],trow,tcol);
		out_receivers = std::vector<int>{trow,tcol};
		
		std::cout << "Srec is " << this->Sreceivers[node] << " node was " << node << std::endl;
		return out_receivers;
	}


	template<class Connector_t, class topo_t>
	void print_receivers(int i,Connector_t& connector, topo_t& ttopography)
	{
		std::cout << std::setprecision(12);
		auto topography = format_input(ttopography);
		auto receivers = this->get_receivers_idx(i, connector);

		std::cout << "Topography is " << topography[i] << "# receivers: " << receivers.size() << std::endl;
		for(auto r: receivers)
		{
			int row,col;
			connector.rowcol_from_node_id(r,row,col);
			std::cout << "Rec " << r << " row " << row << " col " << col << " topo " << topography[r] << std::endl;

		}


		auto neighbours = connector.get_neighbours_idx(i);
		std::cout << "Neighbours are :" << std::endl;

		for(auto r: neighbours)
		{
			int row,col;
			connector.rowcol_from_node_id(r,row,col);
			std::cout << "Neighbour " << r << " row " << row << " col " << col << " topo " << topography[r] << std::endl;
		}
	}


	int get_rec_array_size(){return int(this->links.size());}


	/// Takes an array of nnodes size and sum the values at the outlets
	/// This can be useful for checking mass balances for example
	/// if true, include_internal_pits allow the code to add internal unprocessed pits, wether they are on purpose or not
	template<class Connector_t,class array_t, class T>
	T sum_at_outlets(Connector_t& connector, array_t& tarray, bool include_internal_pits = true)
	{
		auto array = format_input(tarray);
		T out = 0;
		for(int i=0; i<this->nnodes; ++i)
		{
			if (this->Sreceivers[i] == i)
			{
				if(include_internal_pits)
				{
					out += array[i];
				}
				else if(connector.can_flow_out_there(i) )
				{
					out += array[i];
				}
			}
		}
		return out;

	}

	/// Takes an array of nnodes size and sum the values at the outlets
	/// This can be useful for checking mass balances for example
	/// if true, include_internal_pits allow the code to add internal unprocessed pits, wether they are on purpose or not
	template<class Connector_t,class array_t, class out_t>
	out_t keep_only_at_outlets(Connector_t& connector, array_t& tarray, bool include_internal_pits = true)
	{
		auto array = format_input(tarray);
		std::vector<float_t> out = std::vector<float_t> (this->nnodes,0);
		for(int i=0; i<this->nnodes; ++i)
		{
			if (this->Sreceivers[i] == i)
			{
				if(include_internal_pits)
					out[i] = array[i];
				else if(connector.can_flow_out_there(i) )
					out[i] = array[i];
			}
		}
		return format_output(out);

	}

	bool is_Sstack_full()
	{
		if(int(this->Sstack.size()) != this->nnodes)
		{
			std::cout << "stack size (" << this->Sstack.size() << ") is invalid." << std::endl;
			return false;
		}
		std::vector<int> ntimenodes(this->nnodes,0);

		for(auto v:this->Sstack)
		{
			++ntimenodes[v];
		}

		int n_0 = 0,n_p1 = 0;
		for(int i=0; i<this->nnodes; ++i)
		{
			if(ntimenodes[i] == 0)
				++n_0;
			else if(ntimenodes[i]>1)
				++n_p1;
		}

		if(n_0 > 0 || n_p1 > 0)
		{
			std::cout << "Stack issue: " << n_p1 << " nodes appearing more than once and " << n_0 << " nodes not appearing" << std::endl;
			return false;
		}

		std::vector<bool> isdone(this->nnodes,false);
		for(int i = this->nnodes - 1; i>=0; --i)
		{
			auto v = this->Sstack[i];
			isdone[v] = true;
			if(int(v) !=  this->Sreceivers[v])
			{
				if(isdone[this->Sreceivers[v]])
				{
					std::cout << "Receiver processed before node stack is fucked" << std::endl;
					return false;
				}
			}
		}

		return true;

	}


	std::vector<bool> has_Srecs()
	{
		std::vector<bool> haSrecs(this->nnodes,true);
		for(int i=0;i<this->nnodes;++i)
		{
			if(this->Sreceivers[i] == i)
				haSrecs[i] = false;
		}
		return haSrecs;
	}


	bool is_method_cordonnier(std::string method)
	{
		if(method == "cordonnier_fill" || method == "cordonnier_carve")
			return true;
		else
			return false;
	}





/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	                      . - ~ ~ ~ - .
      ..     _      .-~               ~-.
     //|     \ `..~                      `.
    || |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
 \`.-~  o      /       }       |        /    \
 (__          |       /        |       /      `.
  `- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
              |     /          |     /     ~-.     ~- _
              |_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Functions to access sets of nodes draining to/from a single point
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	template<class Connector_t,class out_t>
	out_t get_all_nodes_upstream_of(Connector_t& connector, int node, bool use_Sgraph = true, bool only_SD = false)
	{
		std::vector<int> out;
		if(use_Sgraph)
		{
			out = this->_get_all_nodes_upstream_of_using_graph(connector ,node, only_SD);
		}
		else
			throw std::runtime_error("graph::get_all_nodes_upstream_of::error not implemented yet without graph");

		return format_output(out);
	}

	template< class Connector_t>
	std::vector<int> _get_all_nodes_upstream_of_using_graph(Connector_t& connector,int node, bool only_SD)
	{

		// Formatting the output vector
		std::vector<int> out;
		out.reserve(round(this->nnodes/4));
		// Creating a visited vector tracking which node has been visited
		std::vector<bool> vis(this->nnodes,false);
		// marking the initial node as true
		vis[node] = true;

		for(auto tnode:this->Sstack)
		{
			// ignoring the not ode and outlets
			if(connector.is_active(tnode))
			{
				// Getting the receiver
				int rec = this->Sreceivers[tnode];
				// checkng if receiver is visited but not node
				if(vis[rec] && rec != node)
				{
					// current noer is visited
					vis[tnode] = true;
					// and is draining to this node
					out.emplace_back(tnode);
				}
			}
		}

		// if only steepest descent is needed, we stop there
		if(only_SD)
			return out;

		// else, we have to use a queue to add all the donors
		std::queue<int> toproc;

		// first checking if all the steepest descent nodes I already have there have a not-SD donor
		for(auto v:out)
		{
			// gettign the donors
			auto donors = this->get_donors_idx(v,connector);
			// for all donors of dat nod
			for(auto d:donors)
			{
				// if not visited
				if(vis[d] == false)
				{
					// becomes visited
					vis[d] = true;
					// and I emplace it in the queues
					toproc.emplace(d);
				}
			}
		}		 

		// once this is done I work until the queue is empty
		while(toproc.empty() == false)
		{
			// getting the next node in line
			int next = toproc.front();
			toproc.pop();
			// recording it as draining to the original node
			out.emplace_back(next);
			// getting all its donors
			auto donors = this->get_donors_idx(next,connector);
			for(auto d:donors)
			{
				// if not visited including it (see above)
				if(vis[d] == false)
				{
					vis[d] = true;
					toproc.emplace(d);
				}
			}

		}

		// LOK queue is empty and I have everything I need

		return out;
	}

	template<class Connector_t,class out_t>
	out_t get_all_nodes_downstream_of(Connector_t& connector, int node, bool use_Sgraph = true, bool only_SD = false)
	{
		std::vector<int> out;
		if(use_Sgraph)
		{
			out = this->_get_all_nodes_downstream_of_using_graph(connector ,node, only_SD);
		}
		else
			throw std::runtime_error("graph::get_all_nodes_downstream_of::error not implemented yet without graph");

		return format_output(out);
	}

	template< class Connector_t>
	std::vector<int> _get_all_nodes_downstream_of_using_graph(Connector_t& connector,int node, bool only_SD)
	{

		// Formatting the output vector
		std::vector<int> out;
		out.reserve(round(this->nnodes/4));
		// Creating a visited vector tracking which node has been visited
		std::vector<bool> vis(this->nnodes,false);
		// marking the initial node as true
		vis[node] = true;

		for(int i = this->nnodes - 1; i >=0; --i )
		{
			int tnode = this->Sstack[i];
			// ignoring the not ode and outlets
			if(connector.is_active(tnode))
			{
				// Getting the receiver
				int rec = this->Sreceivers[tnode];
				// checkng if receiver is visited but not node
				if(vis[node] && rec != node && tnode != node)
				{
					// current noer is visited
					vis[rec] = true;
					// and is draining to this node
					out.emplace_back(tnode);
				}
			}
		}

		// if only steepest descent is needed, we stop there
		if(only_SD)
			return out;

		// else, we have to use a queue to add all the receivers
		std::queue<int> toproc;

		// first checking if all the steepest descent nodes I already have there have a not-SD rec
		for(auto v:out)
		{
			// gettign the receivers
			auto receivers = this->get_receivers_idx(v,connector);
			// for all receivers of dat nod
			for(auto r:receivers)
			{
				// if not visited
				if(vis[r] == false)
				{
					// becomes visited
					vis[r] = true;
					// and I emplace it in the queues
					toproc.emplace(r);
				}
			}
		}		 

		// once this is done I work until the queue is empty
		while(toproc.empty() == false)
		{
			// getting the next node in line
			int next = toproc.front();
			toproc.pop();
			// recording it as draining to the original node
			out.emplace_back(next);
			// getting all its receivers
			auto receivers = this->get_receivers_idx(next,connector);
			for(auto r:receivers)
			{
				// if not visited including it (see above)
				if(vis[r] == false)
				{
					vis[r] = true;
					toproc.emplace(r);
				}
			}

		}

		// LOK queue is empty and I have everything I need

		return out;
	}



	template< class Connector_t>
	std::vector<int> _get_flow_acc(Connector_t& connector)
	{
		std::vector<int> flowacc(this->nnodes,0);
		for(int i = this->nnodes-1; i>=0; --i)
		{
			int node = this->Sstack[i];
			if(connector.can_flow_even_go_there(node))
				continue;
			int rec = this->Sreceivers[node];
			if(connector.can_flow_out_there(node) == false)
			{
				flowacc[rec] += flowacc[node] + 1;
			}
		}
		return flowacc;
	}



/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	                      . - ~ ~ ~ - .
      ..     _      .-~               ~-.
     //|     \ `..~                      `.
    || |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
 \`.-~  o      /       }       |        /    \
 (__          |       /        |       /      `.
  `- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
              |     /          |     /     ~-.     ~- _
              |_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Functions to access to bulk receivers/links/donors/...
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/


	template<class out_t>
	out_t get_SFD_receivers()
	{return format_output(this->Sreceivers);}

	template<class out_t>
	out_t get_SFD_dx()
	{return format_output(this->Sdistance2receivers);}

	template<class out_t>
	out_t get_SFD_ndonors()
	{return format_output(this->nSdonors);}

	template<class out_t>
	out_t get_SFD_donors_flat()
	{return format_output(this->Sdonors);}

	template<class out_t>
	out_t get_SFD_donors_list()
	{
		std::vector<std::vector<int> > out(this->nnodes);
		for(int i=0; i < this->nnodes; ++i)
		{
			std::vector<int> tvec;
			for (int j=0; j<this->nSdonors[i]; ++j)
				tvec.emplace_back(this->Sdonors[i * this->n_neighbours +j]);
			out[i] = tvec;
		}

		return out;
	}

	template<class out_t>
	out_t get_SFD_stack()
	{return format_output(this->Sstack);}



	template<class out_t>
	out_t get_MFD_stack()
	{return format_output(this->stack);}

	template<class out_t>
	out_t get_links()
	{return this->links;}

	template<class out_t>
	out_t get_linknodes_flat()
	{return format_output(this->linknodes);}

	template<class out_t>
	out_t get_linknodes_list()
	{
		std::vector<std::vector<int> > out(this->links.size());
		for(size_t i=0; i<this->links.size();++i)
		{
			out[i] = std::vector<int>{this->linknodes[i*2], this->linknodes[i*2+1]};
		}
		return out;
	}

	template<class out_t>
	out_t get_linknodes_list_oriented()
	{
		std::vector<std::vector<int> > out(this->links.size());
		for(size_t i=0; i<this->links.size();++i)
		{
			out[i] = (this->links[i])? std::vector<int>{this->linknodes[i*2], this->linknodes[i*2+1]} : std::vector<int>{this->linknodes[i*2 + 1], this->linknodes[i*2]};
		}
		return out;
	}






	int get_SFD_receivers_at_node(int i)
	{return this->Sreceivers[i];}

	int get_SFD_dx_at_node(int i)
	{return this->Sdistance2receivers[i];}

	int get_SFD_ndonors_at_node(int i)
	{return this->nSdonors[i];}

	template<class out_t>
	out_t get_SFD_donors_at_node(int i)
	{
		std::vector<int> out(this->n_neighbours);
		for (int j=0; j<this->nSdonors[i]; ++j)
			out.emplace_back(this->Sdonors[i * this->n_neighbours +j]);
		return out;
	}


/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	                      . - ~ ~ ~ - .
      ..     _      .-~               ~-.
     //|     \ `..~                      `.
    || |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
 \`.-~  o      /       }       |        /    \
 (__          |       /        |       /      `.
  `- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
              |     /          |     /     ~-.     ~- _
              |_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Functions computing gradients
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	template<class out_t, class topo_t>
	out_t get_SFD_gradient(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		auto gradient = this->_get_SFD_gradient(topography);
		return format_output(gradient);
	}

	template<class topo_t>
	std::vector<float_t> _get_SFD_gradient(topo_t& topography)
	{
		std::vector<float_t> gradient(this->nnodes,0.);
		for(int i=0; i<this->nnodes;++i)
		{
			if(this->Sreceivers[i] != i)
				gradient[i] = (topography[i] - topography[this->Sreceivers[i]])/this->Sdistance2receivers[i];
		}
		return gradient;
	}

	template<class Connector_t,class out_t, class topo_t>
	out_t get_links_gradient(Connector_t& connector, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		std::vector<float_t> gradient = this->_get_links_gradient(connector, topography);
		return format_output(gradient);
	}


	template<class Connector_t, class topo_t>
	std::vector<float_t> _get_links_gradient(Connector_t& connector, topo_t& topography)
	{

		std::vector<float_t> gradient = std::vector<float_t>(this->links.size(), 0);

		for(size_t i = 0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				gradient[i] = std::abs(topography[this->linknodes[i*2] - this->linknodes[i*2 + 1]])/connector.get_dx_from_links_idx(i);
			}
		}

		return gradient;
	}

	template<class Connector_t,class out_t, class topo_t>
	out_t get_MFD_mean_gradient(Connector_t& connector,topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		auto gradient = this->_get_MFD_mean_gradient(connector,topography);
		return format_output(gradient);
	}

	template<class Connector_t, class topo_t>
	std::vector<float_t> _get_MFD_mean_gradient(Connector_t& connector, topo_t& topography)
	{

		std::vector<float_t> gradient = std::vector<float_t>(this->nnodes,0);
		std::vector<int> ngradient = std::vector<int>(this->nnodes,0);

		for(size_t i = 0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				float_t this_gradient = std::abs(topography[this->linknodes[i*2] - this->linknodes[i*2 + 1]])/connector.get_dx_from_links_idx(i);
				auto frto = this->get_from_to_links(i);
				gradient[frto.first] += this_gradient;
				++ngradient[frto.first];
			}
		}

		for(int i=0; i< this->nnodes; ++i)
		{
			if(ngradient[i] > 0)
				gradient[i] = gradient[i]/ngradient[i];
		}

		return gradient;
	}


	template<class Connector_t,class out_t, class topo_t>
	out_t get_MFD_weighted_gradient(Connector_t& connector,topo_t& ttopography, topo_t& tweights)
	{
		auto topography = format_input(ttopography);
		auto weights = format_input(tweights);
		auto gradient = this->_get_MFD_weighted_gradient(connector,topography, weights);
		return format_output(gradient);
	}

	template<class Connector_t, class topo_t>
	std::vector<float_t> _get_MFD_weighted_gradient(Connector_t& connector, topo_t& topography, topo_t& weights)
	{

		std::vector<float_t> gradient = std::vector<float_t>(this->nnodes,0);
		std::vector<float_t> wgradient = std::vector<float_t>(this->nnodes,0.);

		for(size_t i = 0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				float_t this_gradient = std::abs(topography[this->linknodes[i*2] - this->linknodes[i*2 + 1]])/connector.get_dx_from_links_idx(i);
				auto frto = this->get_from_to_links(i);
				gradient[frto.first] += this_gradient * weights[i];
				wgradient[frto.first] += weights[i];
			}
		}

		for(int i=0; i< this->nnodes; ++i)
		{
			if(wgradient[i] > 0)
				gradient[i] = gradient[i]/wgradient[i];
		}

		return gradient;
	}



	/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	                      . - ~ ~ ~ - .
      ..     _      .-~               ~-.
     //|     \ `..~                      `.
    || |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
 \`.-~  o      /       }       |        /    \
 (__          |       /        |       /      `.
  `- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
              |     /          |     /     ~-.     ~- _
              |_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Links utility functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/



	template<class ti_t>
	bool is_link_valid(ti_t i){return (this->linknodes[i*2]>=0)?true:false; }

	template<class ti_t>
	std::pair<ti_t,ti_t> get_from_to_links(ti_t i)
	{
		if(this->links[i])
			return std::make_pair(this->linknodes[i*2], this->linknodes[i*2 + 1]);
		else
			return std::make_pair(this->linknodes[i*2 + 1], this->linknodes[i*2]);

	}



// end of the graph class
};









// end of Dagger namespace
};




















#endif
