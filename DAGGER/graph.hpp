//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef graph_HPP
#define graph_HPP


/*

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
	// -> Steepest donors (nnodes * 8 size)
	// --> Sdonors of node i are located from index i*8 to index i*8 + nSdonors[i] not included
	std::vector<int> Sreceivers,nSdonors,Sdonors;

	// Single graph distance to receivers
	std::vector<float_t> Sdistance2receivers;

	// Steepest slope
	std::vector<float_t> SS;
	
	// Topological order and Single graph topological order
	std::vector<size_t> stack, Sstack;

	
	// default constructor
	graph(){};

	// Classic constructor, simply giving the number of nodes and the number of neighbours per node
	graph(int nnodes, int n_neighbours){this->nnodes = nnodes; this->n_neighbours = n_neighbours;}


	// Initialisation of the graph structure
	// Uses the nnodes, n_neighbours and connector to alloate the single vectors and the irec/linknodes attribute
	template<class Connector_t>
	void init_graph(Connector_t& connector)
	{
		// Allocate vectors
		this->_allocate_vectors();
		connector.fill_linknodes(this->linknodes);
	}

	// used to reinitialise the the vectors
	template<class Connector_t>
	void reinit_graph(Connector_t& connector)
	{
		// Allocate vectors
		this->_reallocate_vectors();
	}

	// Compute the graph using the cordonnier metod to solve the depressions
	// depression_solver is a string of "carve", "fill" or "simple" setting the rerouting method
	// topography is the vetor-like topography
	// connector is the connector
	template<class Connector_t,class topo_t, class out_t>
	out_t compute_graph(
		std::string depression_solver,
	  topo_t& ttopography, 
	  Connector_t& connector, 
	  bool only_SS
	  )
	{
		// Formatting the input to match all the wrappers
		auto topography = format_input(ttopography);

		// Formatting hte output topo
		std::vector<float_t> faketopo(to_vec(topography));


		if(depression_solver == "priority_flood")
			faketopo = connector.PriorityFlood_Wei2018(topography);

		// Making sure the graph is not inheriting previous values
		this->reinit_graph(connector);

		// Updates the links vector and the Srecs vector by checking each link new elevation
		this->update_recs(faketopo, connector);
		// std::cout << "DEBUGGRAPH6::3" << std::endl;
		

		// std::cout << "DEBUGGRAPH6::4" << std::endl;
		
		this->compute_TO_SF_stack_version();
		// std::cout << "DEBUGGRAPH6::5" << std::endl;


		bool need_recompute = false;

		if(depression_solver == "carve" || depression_solver == "fill")
		{
	
			LMRerouter depsolver;
			// std::cout << "DEBUGGRAPH6::prerun" << std::endl;
			need_recompute = depsolver.run(depression_solver, faketopo, connector, this->Sreceivers, this->Sdistance2receivers, this->Sstack, this->linknodes);
			// std::cout << "DEBUGGRAPH6::postrun" << std::endl;
		}

		if(need_recompute)
		{
		
			this->recompute_SF_donors_from_receivers();
		
			this->compute_TO_SF_stack_version();

			if(depression_solver == "carve")
				this->carve_topo_v2(1e-5, connector, faketopo);

			if(only_SS)
				return format_output(faketopo);
					
			this->compute_MF_topological_order_insort(faketopo);
			this->update_Mrecs(faketopo,connector);


			return format_output(faketopo);
			
		}
		else if (only_SS == false)
		{
			// std::cout << "nodep" << std::endl;
			this->compute_MF_topological_order_insort(faketopo);
			return format_output(faketopo);	
		}
		else
		{
			return format_output(faketopo);
		}

	}

	template<class Connector_t,class topo_t, class out_t>
	out_t compute_graph_SS(std::string depression_solver, topo_t& ttopography, Connector_t& connector)
	{
		// std::cout << "DEBUGGRAPH6::1" << std::endl;
		auto topography = format_input(ttopography);
		this->reinit_graph(connector);
		// std::cout << "DEBUGGRAPH6::2" << std::endl;
		this->update_recs(topography,connector);
		// std::cout << "DEBUGGRAPH6::3" << std::endl;
		
		this->compute_SF_donors_from_receivers();
		// std::cout << "DEBUGGRAPH6::4" << std::endl;
		
		this->compute_TO_SF_stack_version();
		// std::cout << "DEBUGGRAPH6::5" << std::endl;

		std::vector<float_t> faketopo(to_vec(topography));
		
		LMRerouter depsolver;
		// std::cout << "DEBUGGRAPH6::prerun" << std::endl;
		bool need_recompute = depsolver.run(depression_solver, faketopo, connector, this->Sreceivers, this->Sdistance2receivers, this->Sstack, this->linknodes);
		// std::cout << "DEBUGGRAPH6::postrun" << std::endl;

		if(need_recompute)
		{
		
			this->recompute_SF_donors_from_receivers();
		
			this->compute_TO_SF_stack_version();

			if(depression_solver == "carve")
				this->carve_topo_v2(1e-5, connector, faketopo);

			return format_output(faketopo);
		}
		else
		{
			// std::cout << "nodep" << std::endl;
			// this->compute_MF_topological_order_insort(faketopo);
			return format_output(faketopo);	
		}

	}

	template<class Connector_t,class topo_t, class out_t>
	out_t compute_graph_nodep(topo_t& ttopography, Connector_t& connector)
	{
		auto topography = format_input(ttopography);
		this->reinit_graph(connector);
		this->update_recs(topography,connector);
		
		this->compute_SF_donors_from_receivers();
		
		this->compute_TO_SF_stack_version();

		std::vector<float_t> faketopo(to_vec(topography));
		this->compute_MF_topological_order_insort(faketopo);

	
		return format_output(faketopo);	

	}

	template<class Connector_t,class topo_t, class out_t>
	out_t compute_graph_PQ( topo_t& ttopography, Connector_t& connector)
	{
		// std::cout << "DEBUGGRAPH6::1" << std::endl;
		auto topography = format_input(ttopography);

		std::vector<float_t> faketopo = connector.PriorityFlood_Wei2018(topography);

		this->reinit_graph(connector);
		// std::cout << "DEBUGGRAPH6::2" << std::endl;
		this->update_recs(faketopo,connector);
		// std::cout << "DEBUGGRAPH6::3" << std::endl;
		
		this->compute_SF_donors_from_receivers();
		// std::cout << "DEBUGGRAPH6::4" << std::endl;
		
		this->compute_TO_SF_stack_version();

		this->compute_MF_topological_order_insort(faketopo);
		// std::cout << "DEBUGGRAPH6::5" << std::endl;

		return format_output(faketopo);	


	}


	template<class Connector_t,class topo_t>
	void update_Mrecs(topo_t& topography, Connector_t& connector)
	{
		for(size_t i = 0; i<this->links.size(); ++i)
		{
			int from = this->linknodes[i*2];
			int to = this->linknodes[i*2 + 1];
			
			if(connector.is_in_bound(from) == false || connector.is_in_bound(to) == false)
				continue;

			if(topography[from] > topography[to])
				this->links[i] = true;
			else
				this->links[i] = false;
		}
	}
	template<class Connector_t,class topo_t>
	void update_some_Mrecs(topo_t& topography, Connector_t& connector, std::vector<int>& some)
	{
		for(size_t j = 0; j<some.size(); ++j)
		{
			int node = some[j];
			auto ilinknodes = connector.get_ilinknodes_from_node(node);
			for(auto i: ilinknodes)
			{
				int from = this->linknodes[i*2];
				int to = this->linknodes[i*2 + 1];
				
				if(connector.is_in_bound(from) == false || connector.is_in_bound(to) == false)
					continue;

				if(topography[from] > topography[to])
					this->links[i] = true;
				else
					this->links[i] = false;
			}
		}
	}

	template<class Connector_t,class topo_t>
	void update_recs(topo_t& topography, Connector_t& connector)
	{
		for(size_t i = 0; i<this->links.size(); ++i)
		{
			int from = this->linknodes[i*2];
			int to = this->linknodes[i*2 + 1];
			
			if(connector.is_in_bound(from) == false || connector.is_in_bound(to) == false)
			{
				continue;
			}

			// if(connector.is_active(from) == false || connector.is_active(to) == false )
			// 	continue;

			float_t dx = connector.get_dx_from_links_idx(i);
			float_t slope = (topography[from] - topography[to])/dx;
			// std::cout << from << "|" << to << "|";

			if(slope>0)
			{
				this->links[i] = true;
				if(this->SS[from]<slope)
				{
					this->Sreceivers[from] = to;
					this->Sdistance2receivers[from] = dx;
					this->SS[from] = slope;
				}
			}
			else
			{
				this->links[i] = false;
				slope = std::abs(slope);
				if(this->SS[to]<slope)
				{
					this->Sreceivers[to] = from;
					this->Sdistance2receivers[to] = dx;
					this->SS[to] = slope;
				}
			}

		}

		this->compute_SF_donors_from_receivers();



	}

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






#ifdef DAGGER_FT_PYTHON
	py::array_t<int,1> get_Sreceivers(){return py::array_t<int,1>(this->Sreceivers.size(),this->Sreceivers.data() ) ;} 
	py::array_t<float_t,1> get_dx_array(){return py::array_t<float_t,1>(this->Sdistance2receivers.size(),this->Sdistance2receivers.data() ) ;} 
#endif




	template<class topo_t>
	void compute_MF_topological_order_insort(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);


		auto yolo = sort_indexes(topography);
		this->stack = std::move(yolo);

	}

	void compute_SF_donors_from_receivers()
	{
		// Initialising the graph dimesions for the donors
		// All of thenm have the graph dimension
		this->Sdonors = std::vector<int>(this->nnodes * 8,-1);
		this->nSdonors = std::vector<int>(this->nnodes,0);

		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;

			this->Sdonors[trec * 8  + this->nSdonors[trec]] = i;
			this->nSdonors[trec] += 1;
		}

	}

	void recompute_SF_donors_from_receivers()
	{

		for(int i=0; i < this->nnodes; ++i)
		{
			for(int j=0; j<8; ++j)
				this->Sdonors[i * 8 + j] = -1;
			this->nSdonors[i] = 0;
		}

		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;
			this->Sdonors[trec * 8  + this->nSdonors[trec]] = i;
			this->nSdonors[trec] += 1;
		}

	}


	void compute_TO_SF_stack_version()
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
					stackhelper.emplace(this->Sdonors[nextnode*8 + j]);
				}

			}

		}
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
			
			// if(node == 148880)
				// std::cout << "ASSESSED" << std::endl;

			if(connector.can_flow_out_there(node) || connector.can_flow_even_go_there(node) == false)
				continue;
			// if(node == 148880)
				// std::cout << "PASSED" << std::endl;
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
	std::vector<int> get_receiver_indices(int i, Connector_t& connector)
	{
		std::vector<int> recs; recs.reserve(8);
		std::vector<int> linkidx = connector.get_neighbour_idx_links(i); 
		int i_r = connector.get_id_right_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r])
				recs.emplace_back(connector.get_right_index(i));
		}
		i_r = connector.get_id_bottomright_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r])
				recs.emplace_back(connector.get_bottomright_index(i));
		}
		i_r = connector.get_id_bottom_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r])
				recs.emplace_back(connector.get_bottom_index(i));
		}
		i_r = connector.get_id_bottomleft_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r])
				recs.emplace_back(connector.get_bottomleft_index(i));
		}
		i_r = connector.get_id_left_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r] == false)
				recs.emplace_back(connector.get_left_index(i));
		}
		i_r = connector.get_id_topleft_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r] == false)
				recs.emplace_back(connector.get_topleft_index(i));
		}
		i_r = connector.get_id_top_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r] == false)
				recs.emplace_back(connector.get_top_index(i));
		}
		i_r = connector.get_id_topright_SMG(i);
		if(i_r >=0 || i_r < connector.nnodes * 4)
		{
			if(this->links[i_r] == false)
				recs.emplace_back(connector.get_topright_index(i));
		}
		return recs;

	}

	template<class Connector_t>
	std::vector<std::pair<int,int> > get_receiver_link_indices(int i, Connector_t& connector)
	{
		std::vector<std::pair<int,int>> recs; recs.reserve(8);

		int i_r = connector.get_id_right_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_right_index(i);
			if(this->links[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_bottomright_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_bottomright_index(i);
			if(this->links[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_bottom_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_bottom_index(i);
			if(this->links[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_bottomleft_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_bottomleft_index(i);
			if(this->links[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_left_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_left_index(i);
			if(this->links[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_topleft_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_topleft_index(i);
			if(this->links[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_top_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_top_index(i);
			if(this->links[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = connector.get_id_topright_SMG(i);
		if(i_r >=0 && i_r < connector.nnodes * 4)
		{
			int ti = connector.get_topright_index(i);
			if(this->links[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		return recs;

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
				auto receivers = this->get_receiver_indices(node, connector);

				std::vector<float_t> slopes(receivers.size());
				float_t sumslopes = 0;
				for(size_t j = 0;j < receivers.size(); ++j)
				{
					int rec = receivers[j];
					slopes[j] = (topography[node] - topography[rec])/connector.dx;
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(size_t j = 0;j < receivers.size(); ++j)
				{
					int rec = receivers[j];
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
		auto receivers = this->get_receiver_indices(i, connector);

		std::cout << "Topography is " << topography[i] << "# receivers: " << receivers.size() << std::endl;
		for(auto r: receivers)
		{
			int row,col;
			connector.rowcol_from_node_id(r,row,col);
			std::cout << "Rec " << r << " row " << row << " col " << col << " topo " << topography[r] << std::endl;

		}


		auto neighbours = connector.get_neighbours_only_id(i);
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







};






// end of Dagger namespace
};




















#endif
