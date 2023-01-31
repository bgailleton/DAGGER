//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef priority_Q_HPP
#define priority_Q_HPP

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


// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"





namespace DAGGER
{

	// Useful namespace for my priority queue
	using PQ_i_d =  std::priority_queue< PQ_helper<int,double>, std::vector<PQ_helper<int,double> >, std::greater<PQ_helper<int,double> > >;


	// Largely adapted from RichDEM
	template <class topo_t, class Connector_t, class Graph_t>
	std::vector<double> PriorityFlood(topo_t& ttopography, Connector_t& connector, Graph_t& graph)
	{

	  PQ_i_d open;
	  std::queue<PQ_helper<int,double> > pit;

	  std::random_device rd; // obtain a random number from hardware
	  std::mt19937 gen(rd()); // seed the generator
	  std::uniform_real_distribution<> distr(1e-7, 1e-6); // define the range

	  auto topography = format_input(ttopography);

	  std::vector<int8_t> closed(connector.nnodes,false);

	  for(int i=0; i<connector.nnodes;++i)
	  {
	  	if(connector.boundaries.can_out(i))
	  	{
	  		open.emplace(PQ_helper<int,double>(i,topography[i]));
	  		closed[i] = true;
	  	}
	    
	  }
	  
		auto neighbours = connector.get_empty_neighbour();
	  while(open.size()>0 || pit.size()>0)
	  {
	    PQ_helper<int,double> c;
	    if(pit.size()>0)
	    {
	      c=pit.front();
	      pit.pop();
	    } 
	    else 
	    {
	      c=open.top();
	      open.pop();
	    }

	    int nn = connector.get_neighbour_idx(c.node,neighbours);

	    for(int j=0;j<nn;++j)
	    {
	      int n = neighbours[j];

	      if(connector.boundaries.no_data(n) || connector.boundaries.force_giving(n)) continue;
	      
	      if(closed[n])
	        continue;

	      closed[n]=true;
	      
	      if(topography[n]<=c.score)
	      {
          topography[n] = c.score + 1e-6 + distr(gen);
	        pit.emplace(PQ_helper<int,double>(n,topography[n]));
	      } 
	      else
	        open.emplace(PQ_helper<int,double>(n,topography[n]));
	    }
	  }
	  return to_vec(topography);

	}



	// Largely adapted from RichDEM
	template <class topo_t, class Connector_t, class Graph_t>
	std::vector<double> PriorityFlood_broken(topo_t& ttopography, Connector_t& connector, Graph_t& graph)
	{
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator
		std::uniform_real_distribution<> distr(1e-7, 1e-6); // define the range

		auto tttopography = format_input(ttopography);
		std::vector<double> topography = to_vec(tttopography);

	  PQ_i_d open;
	  std::queue<PQ_helper<int,double> > pit;

	  uint64_t processed_cells = 0;
	  uint64_t pitc            = 0;
	  int      false_pit_cells = 0;
	  double PitTop = -9999;

	  std::vector<bool> closed(connector.nnodes,false);

	  // RDLOG_PROGRESS<<"Adding cells to the priority queue...";
	  for(int i = 0; i<connector.nnodes; ++i)
	  {
	  	// if(connector.can_flow_even_go_there(i) && i == graph.Sreceivers[i])
	  	if(connector.boundaries.can_out(i))
	  	{
	  		closed[i] = true; 
	  		open.emplace(i,topography[i]);
	  	}
	  }

	  // RDLOG_PROGRESS<<"Performing Priority-Flood+Epsilon...";
	  auto neighbours = connector.get_empty_neighbour();


	  // throw std::runtime_error("???A");
	  // int ii = 0;
	  while(open.size()>0 || pit.size()>0)
	  {
	  	// ++ii;
	  	// if(ii % 100000 == 0) std::cout << open.size() << "||" << pit.size() << std::flush << std::endl;

	    PQ_helper<int,double> c;

	    if(pit.size()>0 && open.size()>0 && open.top().score == pit.front().score)
	    {
	      c=open.top(); open.pop();
	      PitTop=-9999;
	    } 
	    else if(pit.size()>0)
	    {
	      c=pit.front(); pit.pop();
	      if(PitTop==-9999)
	        PitTop=topography[c.node];
	    } 
	    else 
	    {
	      c=open.top(); open.pop();
	      PitTop=-9999;
	    }
	    processed_cells++;

	    int nn = connector.get_neighbour_idx(c.node, neighbours);
	    for(int tnn = 0; tnn<nn; ++nn)
	    {
	    	int n = neighbours[tnn];
	      if(connector.can_flow_even_go_there(n) == false) 
	      	continue;

	      if(closed[n])
	        continue;

	      closed[n]=true;

	      if(topography[n] <= std::nextafter(c.score,std::numeric_limits<double>::infinity()))
	      {
	        if(PitTop!=-9999 && PitTop<topography[n] && std::nextafter(c.score,std::numeric_limits<double>::infinity())>=topography[n])
	          ++false_pit_cells;
	      
	        ++pitc;
	        // topography[n]=std::nextafter(c.score,std::numeric_limits<double>::infinity());
	        topography[n] = c.score + distr(gen);
	        pit.emplace(n,topography[n]);
	      } else
	        open.emplace(n,topography[n]);
	    }
	  }
	
	  return topography;
	}


}



#endif