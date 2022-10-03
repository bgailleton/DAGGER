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
#include <iomanip>

// local includes 
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D8connector.hpp"

#include "wrap_helper.hpp"


namespace DAGGER
{


template<class float_t, class dummy_t = int> // the class type dummy_t is there to bypass pybind11 issue of not being able to bind the same object twice
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


	// What depression resolver to use when computing the graph
	DEPRES depression_resolver = DEPRES::cordonnier_carve;

	// hte minimum slope to impose on the fake topography
	float_t minimum_slope = 1e-4;
	float_t slope_randomness = 1e-6;

	
	// default empty constructor
	graph(){};

	// Classic constructor, simply giving the number of nodes and the number of neighbours per node
	graph(int nnodes, int n_neighbours){this->nnodes = nnodes; this->n_neighbours = n_neighbours;}


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
	Graph functions	- everything about computing, updating and managing the graph, links, nodes and local minima
	The core of the code in other words.
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/


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
	  topo_t& ttopography, // the input topography
	  Connector_t& connector, // the input connector to use (e.g. D8connector)
	  bool only_SD, // only computes the single flow graph if true
	  bool quicksort // computes the MF toposort with a quicksort algo if true, else uses a BFS-based algorithm (which one is better depends on many things)
	  )
	{
		// Formatting the input to match all the wrappers
		auto topography = format_input<topo_t>(ttopography);
		// Formatting the output topo
		std::vector<float_t> faketopo(this->nnodes,0);

		for(int i=0; i<this->nnodes; ++i)
		{
			faketopo[i] = topography[i];
		}

		this->_compute_graph(faketopo,connector,only_SD,quicksort);


		// Finally I format the output topography to the right wrapper
		return format_output<decltype(faketopo), out_t >(faketopo);
	}


	template<class Connector_t>
	void _compute_graph(
		std::vector<float_t>& faketopo,
	  Connector_t& connector, // the input connector to use (e.g. D8connector)
	  bool only_SD, // only computes the single flow graph if true
	  bool quicksort // computes the MF toposort with a quicksort algo if true, else uses a BFS-based algorithm (which one is better depends on many things)

		)
	{

		// std::cout << "DEBUG::STOPOP::1" << std::endl;
		// Checking if the depression method is cordonnier or node
		bool isCordonnier = this->is_method_cordonnier();

		// std::cout << "DEBUG::STOPOP::2" << std::endl;
		// if the method is not Cordonnier -> apply the other first
		if(isCordonnier == false && this->depression_resolver != DEPRES::none)
		{
			// filling the topography with a minimal slope using Wei et al., 2018
			if(this->depression_resolver == DEPRES::priority_flood)
				faketopo = connector.PriorityFlood_Wei2018(faketopo);
			else
				faketopo = connector.PriorityFlood(faketopo);
		}

		// std::cout << "DEBUG::STOPOP::3" << std::endl;

		// Making sure the graph is not inheriting previous values
		this->reinit_graph(connector);

		// std::cout << "DEBUG::STOPOP::4" << std::endl;

		// Updates the links vector and the Srecs vector by checking each link new elevation
		this->update_recs(faketopo, connector);

		// std::cout << "DEBUG::STOPOP::5" << std::endl;
		// Compute the topological sorting for single stack
		// Braun and willett 2014 (modified)
		this->topological_sorting_SF();

		// std::cout << "DEBUG::STOPOP::6" << std::endl;

		// manages the Cordonnier method if needed
		if(isCordonnier)
		{
			
			// LMRerouter is the class managing the different cordonnier's mthod
			LMRerouter<float_t> depsolver;
			depsolver.minimum_slope = this->minimum_slope;
			depsolver.slope_randomness = this->slope_randomness;
			// Execute the local minima solving, return true if rerouting was necessary, meaning that some element needs to be recomputed
			// Note that faketopo are modified in place.
			// std::cout << "wulf" << std::endl;
		// std::cout << "DEBUG::STOPOP::7" << std::endl;
			bool need_recompute = depsolver.run(this->depression_resolver, faketopo, connector, this->Sreceivers, this->Sdistance2receivers, this->Sstack, this->linknodes);		
		// std::cout << "DEBUG::STOPOP::8" << std::endl;

			// Right, if reomputed needs to be
			if(need_recompute)
			{
				
				// Re-inverting the Sreceivers into Sdonors
				this->recompute_SF_donors_from_receivers();
				
				// Recomputing Braun and willett 2014 (modified)
				this->topological_sorting_SF();

				// This is a bit confusing and needs to be changed but filling in done in one go while carving needs a second step here
				if(this->depression_resolver == DEPRES::cordonnier_carve)
					this->carve_topo_v2(connector, faketopo);

				// My work here is done if only SD is needed
				if(only_SD)
					return;
				
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
		if (only_SD == false)
		{
			if(quicksort)
					this->topological_sorting_quicksort(faketopo);
				else
					this->topological_sorting_dag(connector);
		}

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
			// std::cout << from << "|" << to << "||";
			
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

		return format_output<decltype(OUT), out_t>(OUT);
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
	Admin functions	managing attribute and variable initialisations
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/


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



	void set_LMR_method(DEPRES method){this->depression_resolver = method;}
	void set_minimum_slope_for_LMR(float_t slope){this->minimum_slope = slope;}
	void set_slope_randomness_for_LMR(float_t slope)
	{
		if(slope >= this->minimum_slope)
			throw std::runtime_error("slope randomness cannot be >= to the minimum slope and is even reccomended to be at least one or two order of magnitude lower");
		this->slope_randomness = slope;
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
	Functions performing topological sorting
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/



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



	// performs the topological sorting using a BFS algorithm
	// It starts from the receivers-less nodes, then uses a queue to add all the donors of these nodes.
	// Each time a donor pops our of the queue, it increment a visited array recording the number of time a node is visited.
	// if the number of visits equals the number of receiver of the node, it is saved in the stack/
	// The result is a stack of node from the most downstream to the most upstream one
	template< class Connector_t>
	void topological_sorting_dag(Connector_t& connector)
	{
		// nrecs tracks the number of receivers for each nodes
		std::vector<int> nrecs(this->nnodes,0);
		// The queueß
		std::queue<int> toproc;

		// preparing the stack
		this->stack.clear();
		this->stack.reserve(this->nnodes);

		// Iterating through the links
		for(int i = 0; i < int( this->links.size() ) ; ++i)
		{
			// checking the validity of the links
			if(connector.can_flow_even_go_there(this->linknodes[i*2]) == false)
			{
				// if the flow cannot go there, we emplace it in the stack (ultimately it does not matter where they are in the stack)
				this->stack.emplace_back(this->linknodes[i*2]);
				continue;
			}

			// Otherwise incrementing the number of receivers for the right link
			if(this->links[i])
				++nrecs[this->linknodes[i*2]];
			else
				++nrecs[this->linknodes[i*2+1]];
		}


		// Now checking the receiverless nodes and emplacing them in the queueß
		for(int i=0; i<this->nnodes; ++i)
		{
			if(nrecs[i] == 0 && connector.can_flow_even_go_there(i))
				toproc.emplace(i);
		}

		// then as lon g as there are nodes in the queue:
		auto donors = connector.get_empty_neighbour();
		while(toproc.empty() == false)
		{
			// getting the next node
			int next = toproc.front();
			// (and popping the node from the queue)
			toproc.pop();
			// if the node is poped out of the queue -> then it is ready to be in the stack
			// Because we are using a FIFO queue, they are sorted correctly in the queue
			this->stack.emplace_back(next);
			// getting the idx of the donors
			int nn = this->get_donors_idx(next, connector, donors);
			for(int td=0;td<nn;++td)
			{
				int d = donors[td];
				// Decrementing the number of receivers (I use it as a visited vector)
				--nrecs[d];
				// if it has reached 0, the rec is ready to be put in the FIFO
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
		auto topography = format_input<topo_t>(ttopography);
		// Dortng by index
		auto yolo = sort_indexes(topography);
		// the sorted index is now the stack
		this->stack = std::move(yolo);
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
	Functions affecting topography from corrected receivers
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	/// this function enforces minimal slope 
	/// It starts from the most upstream part of the landscapes and goes down following the Sreceiver route
	/// it carve on the go, making sure the topography of a receiver is lower
	template<class Connector_t, class topo_t>
	void carve_topo_v2(Connector_t& connector, topo_t& topography)
	{

		// Traversing the (SFD) stack on the reverse direction
		for(int i=this->nnodes-1; i >= 0; --i)
		{
			// Getting the node
			int node  = this->Sstack[i];
			// Checking its validiyt AND if it is not a base level
			if(connector.can_flow_out_there(node) || connector.can_flow_even_go_there(node) == false)
				continue;
			// Getting the single receiver info
			int rec = this->Sreceivers[node];
			// Checking the difference in elevation
			float_t dz = topography[node] - topography[rec];
			// if the difference in elevation is bellow 0 I need to carve
			if(dz <= 0)
			{
				// And I do ! Note that I add some very low-grade randomness to avoid flat links
				topography[rec] = topography[node] - this->minimum_slope + connector.randu->get() * this->slope_randomness;// * d2rec;
			}
		}
	}

	/// Opposite of the above function
	/// It starts from the most dowstream nodes and climb its way up.
	/// when a node is bellow its receiver, we correct the slope
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
				topography[node] = topography[rec] + slope + connector.randu->get() * 1e-6;// * d2rec;
				to_recompute.emplace_back(node);
			}
		}
		return to_recompute;
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
	Accessing receivers/donors/... for a single node
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/


	// get receivers of node i and put them in the recs vector fed in
	// It returns the number of receivers in the recs vector
	// THis whole process optimises repeated receiver fetching, by never reallocating/initialising the vector recs
	template<class Connector_t>
	int get_receivers_idx(int i, Connector_t& connector, std::vector<int>& recs)
	{
		// getting the related links stroing them temporarily in the rec vec
		int nli = connector.get_neighbour_idx_links(i,recs);

		// going through the linksß
		// The idx is the idx of insertion in the recs vectors
		// the idx of insertion is always <= of the index of reading (both receivers and links-to-assess are stored in the recs vector)	
		int idx = 0;
		// counter used to keep track of the number of receivers: starts at the number of links related to the given node and get decremented at each donor link 
		int newli = nli;

		// Iterating through the links
		for(int ti=0;ti<nli;++ti)
		{
			// current link index
			int li = recs[ti];

			// this link is a rec if the node is the first of the linknode and links is true
			if(i == this->linknodes[li*2] && this->links[li])
			{
				// in which case the receiver of the current node is the +  1
				recs[idx] = this->linknodes[li*2 + 1];
				++idx;
			}
			// OR if the current node is the +1 and the links false
			else if(this->links[li] == false && i == this->linknodes[li*2 + 1])
			{
				// in which case the receivers is the 0 node
				recs[idx] = this->linknodes[li*2];
				++idx;
			}
			else
			{
				// otherwise, it's not a rec and we decrease the newli
				--newli;
			}
		}
		// recs is changed in place, and we return the number of recs newli
		return newli;
	}


	// Getting the id of the receivers in the links array
	// see get_receivers_idx for full comments about the section
	template<class Connector_t>
	int get_receivers_idx_links(int i, Connector_t& connector, std::vector<int>& recs)
	{
		// getting the related links
		int nli = connector.get_neighbour_idx_links(i,recs);

		// going through the linksß
		int idx = 0;
		int newli = nli;

		for(int ti=0;ti<nli;++ti)
		{
			// checking the orientation
			int li = recs[ti];
			if(i == this->linknodes[li*2] && this->links[li])
			{
				recs[idx] = li;
				++idx;
			}
			else if(this->links[li] == false && i == this->linknodes[li*2 + 1])
			{
				recs[idx] = li;
				++idx;
			}
			else
			{
				--newli;
			}
		}
		return newli;
	}


	// // returns array of pair {node, receivers} ofr a single node i
	// // It can be useful if you wanna calculate  gradient or something similar
	// template<class Connector_t>
	// int get_node_receivers_idx_pair(int i, Connector_t& connector, std::vector<std::pair<int,int> >& recs)
	// {
	// 	// getting the related links
	// 	int nli = connector.get_neighbour_idx_links(i,recs);

	// 	// going through the linksß
	// 	int idx = 0;
	// 	int newli = nli;

	// 	for(int ti=0;ti<nli;++ti)
	// 	{
	// 		// checking the orientation
	// 		int li = recs[ti];

	// 		if(i == this->linknodes[li*2] && this->links[li])
	// 		{
	// 			recs[idx] = std::make_pair(i, this->linknodes[li*2 + 1]);
	// 			++idx;
	// 		}
	// 		else if(this->links[li] == false && i == this->linknodes[li*2 + 1])
	// 		{
	// 			recs[idx] = std::make_pair(i, this->linknodes[li*2]);
	// 			++idx;
	// 		}
	// 		else
	// 		{
	// 			--newli;
	// 		}
	// 	}

	// 	return newli;
	// }

	// Getting donor indicies
	// see get_receivers_idx for full comments about the section
	template<class Connector_t>
	int get_donors_idx(int i, Connector_t& connector, std::vector<int>& dons)
	{
		// getting the related links
		int nli = connector.get_neighbour_idx_links(i,dons);

		// going through the linksß
		int idx = 0;
		int newli = nli;

		for(int ti=0;ti<nli;++ti)
		{
			// checking the orientation
			int li = dons[ti];
			if(i == this->linknodes[li*2] && this->links[li] == false)
			{
				dons[idx] = this->linknodes[li*2 + 1];
				++idx;
			}
			else if(this->links[li] == true && i == this->linknodes[li*2 + 1])
			{
				dons[idx] = this->linknodes[li*2];
				++idx;
			}
			else
			{
				--newli;
			}
		}
		return newli;
	}

	// getting links indices of hte donors (in the links array)
	// see get_receivers_idx for full comments about the section
	template<class Connector_t>
	int get_donors_idx_links(int i, Connector_t& connector, std::vector<int> & dons)
	{
		// getting the related links
		int nli = connector.get_neighbour_idx_links(i,dons);

		// going through the linksß
		int idx = 0;
		int newli = nli;

		for(int ti=0;ti<nli;++ti)
		{
			// checking the orientation
			int li = dons[ti];
			if(i == this->linknodes[li*2] && this->links[li] == false)
			{
				dons[idx] = li;
				++idx;
			}
			else if(this->links[li] == true && i == this->linknodes[li*2 + 1])
			{
				dons[idx] = li;
				++idx;
			}
			else
			{
				--newli;
			}
		}
		return newli;
	}


	// Deprecated???
	template<class Connector_t, class topo_t, class out_t>
	out_t get_DA_proposlope(Connector_t& connector, topo_t& ttopography)
	{
		auto topography = format_input<topo_t>(ttopography);

		std::vector<float_t> DA(connector.nnodes,0.);
		auto reclinks = connector.get_empty_neighbour();
		std::vector<float_t> slopes(reclinks.size(),0);

		for(int i = connector.nnodes - 1; i>=0; --i)
		{
			int node = this->stack[i];
			DA[node] += connector.get_area_at_node(node);

			if(connector.is_active(node))
			{
				int nn = this->get_receivers_idx_links(node, connector,reclinks);

				float_t sumslopes = 0;

				for(int j = 0; j < nn; ++j)
				{
					int li = reclinks[j];
					int rec = this->get_to_links(li);
					slopes[j] = (topography[node] - topography[rec])/connector.get_dx_from_links_idx(li);
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(int j = 0;j < nn; ++j)
				{
					int li = reclinks[j];
					int rec = this->get_to_links(li);
					DA[rec] += DA[node] * slopes[j]/sumslopes;
				}

			}

		}

		return format_output<decltype(DA), out_t>(DA);
	}


	// Deprecated???
	template<class Connector_t, class topo_t, class out_t>
	out_t get_DA_SS(Connector_t& connector, topo_t& ttopography)
	{
		auto topography = format_input<topo_t>(ttopography);

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

		return format_output<decltype(DA), out_t>(DA);
	}

	// Debug function printing to the promt the single receiver of a node
	// WIll probably get deprecated
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
		auto topography = format_input<topo_t>(ttopography);
		
		auto receivers = connector.get_empty_neighbour();
		int nn = this->get_receivers_idx(i, connector, receivers);

		std::cout << "Topography is " << topography[i] << "# receivers: " << nn << std::endl;
		for(int tr = 0; tr<nn; ++tr)
		{
			int r = receivers[tr];
			int row,col;
			connector.rowcol_from_node_id(r,row,col);
			std::cout << "Rec " << r << " row " << row << " col " << col << " topo " << topography[r] << std::endl;

		}


		auto neighbours = connector.get_empty_neighbour();
		nn = connector.get_neighbour_idx(i, neighbours);
		std::cout << "Neighbours are :" << std::endl;

		for(int tr = 0; tr<nn; ++tr)
		{
			int r = neighbours[tr];
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
		auto array = format_input<array_t>(tarray);
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
		auto array = format_input<array_t>(tarray);
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
		return format_output<decltype(out), out_t>(out);

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


	bool is_method_cordonnier()
	{
		if(this->depression_resolver == DEPRES::cordonnier_carve || this->depression_resolver == DEPRES::cordonnier_fill)
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

		return format_output<decltype(out), out_t>(out);
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
		auto donors = connector.get_empty_neighbour();
		for(auto v:out)
		{
			// gettign the donors
			int nn = this->get_donors_idx(v,connector, donors);
			// for all donors of dat nod
			for(int td=0; td<nn;++td)
			{
				int d = donors[td];
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
			int nn = this->get_donors_idx(next,connector, donors);
			// for all donors of dat nod
			for(int td=0; td<nn;++td)
			{
				int d = donors[td];
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

		return format_output<decltype(out), out_t>(out);
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

		auto receivers = connector.get_empty_neighbour();

		// first checking if all the steepest descent nodes I already have there have a not-SD rec
		for(auto v:out)
		{
			// gettign the receivers
			int nn = this->get_receivers_idx(v,connector, receivers);
			// for all receivers of dat nod
			for(int tr = 0; tr < nn; ++tr)
			{
				int r = receivers[tr];
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
			int nn = this->get_receivers_idx(next,connector, receivers);
			for(int tr = 0; tr < nn; ++tr)
			{
				int r = receivers[tr];
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
	{return format_output<std::vector<int>, out_t>(this->Sreceivers);}

	template<class out_t>
	out_t get_SFD_dx()
	{return format_output<std::vector<float_t>, out_t>(this->Sdistance2receivers);}

	template<class out_t>
	out_t get_SFD_ndonors()
	{return format_output<std::vector<int>, out_t>(this->nSdonors);}

	template<class out_t>
	out_t get_SFD_donors_flat()
	{return format_output<std::vector<int>, out_t>(this->Sdonors);}

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
	{return format_output<std::vector<size_t>, out_t>(this->Sstack);}



	template<class out_t>
	out_t get_MFD_stack()
	{return format_output<std::vector<size_t>, out_t>(this->stack);}

	template<class out_t>
	out_t get_links()
	{return this->links;}

	template<class out_t>
	out_t get_linknodes_flat()
	{return format_output<std::vector<int>, out_t>(this->linknodes);}

	template<class out_t>
	out_t get_linknodes_flat_D4()
	{
		std::vector<int> temp(int(this->linknodes.size()/2),-1);
		int j = 0;
		int counter = -1;
		for(int i=0; i< int(this->linknodes.size()); i += 2)
		{
			++counter;
			if(counter == 0 || counter == 2)
			{
				temp[j] = this->linknodes[i];
				++j;
				temp[j] = this->linknodes[i+1];
				++j;
			}

			if(counter == 3)
				counter = -1;
		}

		return format_output<std::vector<int>, out_t>(temp);

	}

	template<class out_t, class Connector_t>
	out_t get_linkdx_flat_D4(Connector_t& connector)
	{
		std::vector<float_t> temp(int(this->links.size()/2),-1);
		int j = 0;
		int counter = -1;
		for(int i=0; i< int(this->links.size()); ++i)
		{
			++counter;

			if(counter == 0 || counter == 2)
			{
				float_t dx = connector.get_dx_from_links_idx(i);
				temp[j] = dx;
				++j;
			}

			if(counter == 3)
				counter = -1;
		}

		return format_output<std::vector<float_t>, out_t>(temp);

	}

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
		auto topography = format_input<topo_t>(ttopography);
		auto gradient = this->_get_SFD_gradient(topography);
		return format_output<decltype(gradient), out_t>(gradient);
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
		auto topography = format_input<topo_t>(ttopography);
		std::vector<float_t> gradient = this->_get_links_gradient(connector, topography);
		return format_output<decltype(gradient), out_t>(gradient);
	}


	template<class Connector_t, class topo_t>
	std::vector<float_t> _get_links_gradient(Connector_t& connector, topo_t& topography)
	{

		std::vector<float_t> gradient = std::vector<float_t>(this->links.size(), 0);

		for(size_t i = 0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				gradient[i] = std::abs(topography[this->linknodes[i*2]] - topography[this->linknodes[i*2 + 1]])/connector.get_dx_from_links_idx(i);
			}
		}

		return gradient;
	}

	template<class Connector_t,class out_t, class topo_t>
	out_t get_MFD_mean_gradient(Connector_t& connector,topo_t& ttopography)
	{
		auto topography = format_input<topo_t>(ttopography);
		auto gradient = this->_get_MFD_mean_gradient(connector,topography);
		return format_output<decltype(gradient), out_t>(gradient);
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
		auto topography = format_input<topo_t>(ttopography);
		auto weights = format_input<topo_t>(tweights);
		auto gradient = this->_get_MFD_weighted_gradient(connector,topography, weights);
		return format_output<decltype(gradient), out_t>(gradient);
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


	template<class topo_t>
	std::vector<float_t> _get_max_val_link_array(topo_t& array)
	{
		std::vector<float_t> tmax(this->nnodes,0);
		for(size_t i=0; i < this->links.size(); ++i)
		{
			if(this->is_link_valid(i) ==false)
				continue;
			int go = this->get_from_links(i);
			if(tmax[go]<array[i]) tmax[go] = array[i];
		}
		return tmax;
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
	Functions to calculate weights
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/	

	template<class out_t, class topo_t>
	out_t get_link_weights(topo_t& tgradient, float_t exp)
	{
		std::vector<float_t> weights(this->links.size(),0.);
		auto gradient = format_input<topo_t>(tgradient);

		if(exp <= 0)
		{
			this->_get_link_weights_f_nrecs(weights);
		}
		else if(exp == 1)
		{
			this->_get_link_weights_proposlope(weights, gradient);
		}
		else
		{
			this->_get_link_weights_exp(weights, gradient, exp);
		}

		return format_output<decltype(weights), out_t>(weights);
	}

	void _get_link_weights_f_nrecs(std::vector<float_t>& weights)
	{
		auto nrecs = this->get_n_receivers();
		for(size_t i=0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				int nr = nrecs[this->get_from_links(i)];
				if(nr > 0)
				{
					weights[i] = 1./nr;
				}
				else
					weights[i] = 1.;
			}
		}

	}

	template<class topo_t>
	void _get_link_weights_proposlope(std::vector<float_t>& weights, topo_t& gradient)
	{
		std::vector<float_t> sumgrad(this->nnodes,0.);
		for(size_t i=0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				sumgrad[this->get_from_links(i)] += gradient[i];
			}
		}

		for(size_t i=0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
				weights[i] = gradient[i]/sumgrad[this->get_from_links(i)];
		}


	}

	template<class topo_t>
	void _get_link_weights_exp(std::vector<float_t>& weights, topo_t& gradient, float_t exp)
	{
		std::vector<float_t> sumgrad(this->nnodes,0.);
		for(size_t i=0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				sumgrad[this->get_from_links(i)] += std::pow(gradient[i],exp);
			}
		}

		for(size_t i=0; i< this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
				weights[i] = std::pow(gradient[i],exp)/sumgrad[this->get_from_links(i)];
		}
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
	Functions to propagate signal upstream and downstream
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	template<class Connector_t,class out_t>
	out_t accumulate_constant_downstream_SFD(Connector_t& connector,float_t var)
	{
		std::vector<float_t> out = this->_accumulate_constant_downstream_SFD(connector,var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class Connector_t>
	std::vector<float_t> _accumulate_constant_downstream_SFD(Connector_t& connector, float_t var)
	{
		std::vector<float_t> out(this->nnodes, 0);
		for(int i = this->nnodes - 1; i>=0; --i)
		{
			int node = this->Sstack[i];
			if(connector.can_flow_even_go_there(node) == false)
				continue;

			out[node] += var;

			if(connector.can_flow_out_there(node))
				continue;

			out[this->Sreceivers[node]] += out[node];
			
		}

		return out;
	}


	template<class Connector_t,class out_t, class topo_t>
	out_t accumulate_variable_downstream_SFD(Connector_t& connector,topo_t& tvar)
	{
		auto var = format_input<topo_t>(tvar);
		std::vector<float_t> out = this->_accumulate_variable_downstream_SFD(connector,var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class Connector_t, class topo_t>
	std::vector<float_t> _accumulate_variable_downstream_SFD(Connector_t& connector, topo_t& var)
	{
		std::vector<float_t> out(this->nnodes, 0);
		for(int i = this->nnodes - 1; i>=0; --i)
		{
			int node = this->Sstack[i];
			if(connector.can_flow_even_go_there(node) == false)
				continue;

			out[node] += var[node];

			if(connector.can_flow_out_there(node))
				continue;

			out[this->Sreceivers[node]] += out[node];
			
		}

		return out;
	}


	template<class Connector_t, class topo_t, class out_t>
	out_t accumulate_constant_downstream_MFD(Connector_t& connector, topo_t& tweights,float_t var)
	{
		auto weights = format_input<topo_t>(tweights);
		std::vector<float_t> out = this->_accumulate_constant_downstream_MFD(connector, weights ,var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class Connector_t, class topo_t>
	std::vector<float_t> _accumulate_constant_downstream_MFD(Connector_t& connector, topo_t& weights, float_t var)
	{
		std::vector<float_t> out(this->nnodes, 0);
		auto reclinks = connector.get_empty_neighbour();
 		for(int i = this->nnodes - 1; i>=0; --i)
		{

			int node = this->stack[i];
			if(connector.can_flow_even_go_there(node) == false)
				continue;

			out[node] += var;

			if(connector.can_flow_out_there(node))
				continue;

			int nn = this->get_receivers_idx_links(node,connector, reclinks);
			for (int ttl = 0; ttl< nn; ++ttl)
			{
				int ti = reclinks[ttl];
				int rec = this->get_to_links(ti);
				if(connector.is_in_bound(rec))
					out[rec] += out[node] * weights[ti];
			}
			
		}

		return out;
	}

	template<class Connector_t, class topo_t, class out_t>
	out_t accumulate_variable_downstream_MFD(Connector_t& connector, topo_t& tweights, topo_t& tvar)
	{
		auto weights = format_input<topo_t>(tweights);
		auto var = format_input<topo_t>(tvar);
		std::vector<float_t> out = this->_accumulate_variable_downstream_MFD(connector, weights ,var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class Connector_t, class topo_t>
	std::vector<float_t> _accumulate_variable_downstream_MFD(Connector_t& connector, topo_t& weights, topo_t& var)
	{
		std::vector<float_t> out(this->nnodes, 0);
		auto reclinks = connector.get_empty_neighbour();
		for(int i = this->nnodes - 1; i>=0; --i)
		{
			int node = this->stack[i];
			if(connector.can_flow_even_go_there(node) == false)
				continue;

			out[node] += var[node];

			if(connector.can_flow_out_there(node))
				continue;

			int nn = this->get_receivers_idx_links(node,connector,reclinks);
			for (int tr=0; tr<nn;++tr)
			{
				int ti = reclinks[tr];
				int rec = this->get_to_links(ti);
				if(connector.is_in_bound(rec))
					out[rec] += out[node] * weights[ti];
			}
			
		}

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

	template<class ti_t>
	ti_t get_from_links(ti_t i)
	{
		if(this->links[i])
			return this->linknodes[i*2];
		else
			return this->linknodes[i*2 + 1];
	}

	template<class ti_t>
	ti_t get_to_links(ti_t i)
	{
		if(this->links[i] == false)
			return this->linknodes[i*2];
		else
			return this->linknodes[i*2 + 1];
	}

	template<class ti_t>
	ti_t get_other_node_from_links(ti_t li, ti_t ni)
	{
		if(this->linknodes[li*2 + 1] == ni)
			return this->linknodes[li*2];
		else if (this->linknodes[li*2] == ni)
			return this->linknodes[li*2 + 1];
		else
			throw std::runtime_error("DAGGER::graph::get_other_node_from_links::trying to access other nodes from a link that does not contains the reference node");
	}

	std::vector<int> get_n_receivers()
	{
		std::vector<int> nrecs(this->nnodes,0);
		for(size_t i = 0; i<this->links.size(); ++i)
		{
			if(this->is_link_valid(i))
			{
				auto frto = this->get_from_to_links(i);
				++nrecs[frto.first];
			}
		}
		return nrecs;
	}

	template<class Connector_t>
	void speed_test_links(Connector_t& connector)
	{
		ocarina epona;
		epona.tik();
		int nrecs = 0;
		// for(int i =0; i < this->nnodes; ++i)
		// {
		// 	auto alllinks = this->get_receivers_idx_links(i, connector);
		// 	nrecs += alllinks.size();
		// }
		// epona.tok("Getting recs");
		// std::cout << "I have " << nrecs << std::endl;
		
		// nrecs = 0;

		// epona.tik();
		// for(int i =0; i < this->nnodes; ++i)
		// {
		// 	auto alllinks = connector.get_ilinknodes_from_node(i);
		// 	nrecs += alllinks.size();
		// }
		// epona.tok("Getting links");
		// std::cout << "I have " << nrecs << std::endl;

		nrecs = 0;

		epona.tik();
		for(int i =0; i < this->nnodes; ++i)
		{
			auto alllinks = connector.get_ilinknodes_from_nodev2(i);
			nrecs += alllinks.size();
		}

		epona.tok("Getting linksv2");
		std::cout << "I have " << nrecs << std::endl;

		nrecs = 0;

		epona.tik();
		std::vector<std::pair<int,bool> > these = {std::make_pair(0,false),std::make_pair(0,false),std::make_pair(0,false),std::make_pair(0,false),std::make_pair(0,false),std::make_pair(0,false),std::make_pair(0,false),std::make_pair(0,false)};
		for(int i =0; i < this->nnodes; ++i)
		{
			connector.get_ilinknodes_from_nodev3(i,these);
			nrecs += these.size();
		}
		epona.tok("Getting linksv3");
		std::cout << "I have " << nrecs << std::endl;

		// epona.tik();
		// std::vector<int> these2 = {0,0,0,0,0,0,0,0};
		// for(int i =0; i < this->nnodes; ++i)
		// {
		// 	auto nn = connector.get_ilinknodes_from_nodev3_light(i,these2);
		// 	nrecs += nn;
		// }
		// epona.tok("Getting linksv3.2");
		// std::cout << "I have " << nrecs << std::endl;
		
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
	Distance utility functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/


	/// this function computes the flow distance from model outlets using the siungle direction graph
	template<class Connector_t, class out_t>
	out_t get_SFD_distance_from_outlets(Connector_t& connector)
	{
		std::vector<float_t> distfromoutlet(this->nnodes,0.);
		this->_get_SFD_distance_from_outlets(connector,distfromoutlet);
		return format_output<decltype(distfromoutlet), out_t >(distfromoutlet);
	}

	template<class Connector_t>
	void _get_SFD_distance_from_outlets(Connector_t& connector, std::vector<float_t>& distfromoutlet)
	{
		// just iterating through the Sstack in the upstream direction adding dx to the receiver
		for(int i=0; i<this->nnodes; ++i)
		{
			// next node in the stack
			int node = this->Sstack[i];
			// checking if active
			if(connector.is_active(node) == false)
				continue;
			// Getting the receiver
			int rec = this->Sreceivers[node];
			// And integrating the distance from outlets
			distfromoutlet[node] = distfromoutlet[rec] + this->Sdistance2receivers[node];
		}

	}

	/// this function computes the flow distance from model outlets using the siungle direction graph
	template<class Connector_t, class out_t>
	out_t get_SFD_min_distance_from_sources(Connector_t& connector)
	{
		std::vector<float_t> distfromsources(this->nnodes,0.);
		this->_get_SFD_min_distance_from_sources(connector,distfromsources);
		return format_output<decltype(distfromsources), out_t >(distfromsources);
	}

	template<class Connector_t>
	void _get_SFD_min_distance_from_sources(Connector_t& connector, std::vector<float_t>& distfromsources)
	{
		// just iterating through the Sstack in the upstream direction adding dx to the receiver
		for(int i=this->nnodes - 1; i>=0; --i)
		{
			// next node in the stack
			int node = this->Sstack[i];
			// checking if active
			if(connector.is_active(node) == false)
				continue;
			int rec = this->Sreceivers[node];
			if(distfromsources[rec] == 0 || distfromsources[rec] > distfromsources[node] + this->Sdistance2receivers[node])
				distfromsources[rec] = distfromsources[node] + this->Sdistance2receivers[node];
		}

	}


	/// this function computes the flow distance from model outlets using the siungle direction graph
	template<class Connector_t, class out_t>
	out_t get_SFD_max_distance_from_sources(Connector_t& connector)
	{
		std::vector<float_t> distfromsources(this->nnodes,0.);
		this->_get_SFD_max_distance_from_sources(connector,distfromsources);
		return format_output<decltype(distfromsources), out_t >(distfromsources);
	}

	template<class Connector_t>
	void _get_SFD_max_distance_from_sources(Connector_t& connector, std::vector<float_t>& distfromsources)
	{
		// just iterating through the Sstack in the upstream direction adding dx to the receiver
		for(int i=this->nnodes - 1; i>=0; --i)
		{
			// next node in the stack
			int node = this->Sstack[i];
			// checking if active
			if(connector.is_active(node) == false)
				continue;
			int rec = this->Sreceivers[node];
			if(distfromsources[rec] == 0 || distfromsources[rec] < distfromsources[node] + this->Sdistance2receivers[node])
				distfromsources[rec] = distfromsources[node] + this->Sdistance2receivers[node];
		}

	}


// end of the graph class
};









// end of Dagger namespace
};




















#endif
