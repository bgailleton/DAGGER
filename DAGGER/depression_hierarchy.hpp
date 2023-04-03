//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef DEPHIERCH_HPP
#define DEPHIERCH_HPP

// STL imports
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
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
	using PQH = PQ_helper<int,double>;

	template<class Connector_t>
	int reroute(
		Connector_t* connector, 
		// std::vector<int>& baslab, 
		std::vector<int>& Sreceivers, 
		int from, 
		int pass 
		)
	{
		// std::cout << "rerouting" << std::endl;
		int A = from;
		int B = Sreceivers[A];
		bool goon = true;
		while(goon)
		{
			int C = Sreceivers[B];
			if(C == B)
				goon = false;
			Sreceivers[B] = A;
			A = B;
			B = C;
		}

		// int i=0;
		// do
		// {
		// 	// ++i;
		// 	int C = Sreceivers[B];
		// 	Sreceivers[B] = A;
		// 	connector->debug_print_row_col(B);
		// 	// std::cout << " now gives to "; connector->debug_print_row_col(A); std::cout << "||";
		// 	A = B;
		// 	B = C;
		// 	// if(i >1000)
		// 		// std::cout << A  << "|" << B << std::endl;
		// }while(Sreceivers[B] != B);

		if(pass != -1) Sreceivers[from] = pass;
		else Sreceivers[pass] = pass;

		return B;
		
		// if(passout) Sreceivers[pass] = pass;
	}




	template<class float_t, class topo_t, class Connector_t>
	bool simple_depression_hierarchy(topo_t& topography, Connector_t* connector, std::vector<size_t>& Sstack, std::vector<int>& Sreceivers)
	{

		// PROBLEM IDENTIFIED:: WHEN REROUTING THE SF, THE SFD DOES NOT NECESSARILY DRAINS TO THE BASIN, ITSELF CALCULATED IN MFD LIKE
	
		// std::cout << "running?" << std::endl;

		// record the basin labels
		std::vector<int> baslab(connector->nnodes, -1);

		// the main priority queue
		PQ_i_d PQ;


		// std::unordered_map<int,int> pit2outlet;
		
		// Labelling all the nodes draining to the sea, pushing the local pits to the priority queue
		int lab = 1;
		for(int i =0; i<connector->nnodes;++i)
		{
			int node = int(Sstack[i]);

			if(connector->boundaries.no_data(node))
				continue;

			int rec = Sreceivers[node];

			if(node == rec)
			{
				
				if(connector->flow_out_model(node))
				{
					// std::cout << "OUT!!! :: ";
					// connector->debug_print_row_col(node); 
					baslab[node] = 0;
				}
				else
				{
					baslab[node] = lab;
					PQ.emplace(PQH(node,topography[node]));
					std::cout << "EMPLACING !!! :: ";
					connector->debug_print_row_col(node); 
					++lab;
				}
			}

			if(baslab[rec] == 0)
				baslab[node] = 0;
		}

		if(PQ.empty())
			return false;

		
		std::vector<std::uint8_t> basopen(lab, 0);
		basopen[0] = true;
		std::cout << "BASOPEN::" << std::to_string(basopen.size()) << std::endl;

		std::vector<std::vector<std::array<int,2> > > bas2passes(lab, std::vector<std::array<int,2> >());


		// std::cout << "PQ size = " << PQ.size() << std::endl;;

		auto neighbours = connector->get_empty_neighbour();
		// int yolo = 0;
		while(PQ.empty() == false)
		{
			// std::cout << PQ.size() << std::endl;
			// getting the next node
			auto next = PQ.top(); PQ.pop();

			int tbas = baslab[next.node];

			// std::cout << std::endl << "processing ";
			// connector->debug_print_row_col(next.node);
			// std::cout << ":: basin " << tbas << " (out:" << std::to_string(basopen[tbas]) << ")";
			// if(tbas == 0)
			// 	throw std::runtime_error("should not happen");

			if(basopen[tbas] == true)
			{
				// std::cout << " pass!";
				continue;
			}

			if(connector->boundaries.can_out(next.node))
			{
				// std::cout << "A" << std::endl;
				// std::cout << " reroute from edges";
				reroute(connector, Sreceivers, next.node, -1);
				// std::cout << " done";
				// std::cout << "B" << std::endl;
				basopen[tbas] = true;
				continue;
			}

			int nn = connector->get_neighbour_idx(next.node, neighbours);

			bool does_it_outlets = false;
			PQH node_outlet(-1, topography[next.node] );
			int restor_Srec = -1;
			// if(tbas == 2)
			// 	std::cout << "Processing node "; connector->debug_print_row_col(next.node);

			for(int j = 0; j<nn; ++j)
			{
				int nei = neighbours[j];
				int bnei = baslab[nei];

				// std::cout << tbas << "|" << bnei << std::endl;

				// if(basopen[bnei] == true || (topography[nei] < node_outlet.score && bnei != tbas))
				if(bnei >= 0)
				{
					if(bnei != tbas)
					{

						does_it_outlets = true;
						node_outlet.node = nei;
						node_outlet.score = topography[nei];
						restor_Srec = Sreceivers[nei];
					}
				}

				if(bnei == -1)
				{
					PQ.emplace(PQH(nei,topography[nei]));
					baslab[nei] = tbas;
					bnei = tbas;
					Sreceivers[nei] = next.node;

				}
				

			}

			// if(does_it_outlets && tbas ==2)
			// {
			// 	std::cout << "Processed " << next.node << " | outlet: " << node_outlet.node << std::endl;  
			// }

			// found an outlet, rereouting
			if(does_it_outlets)
			{
				// std::cout << " potential outlet. ";

				int bnei = baslab[node_outlet.node];
				Sreceivers[node_outlet.node] = restor_Srec;
				// if(bnei == 0)
				// 	throw

				if(basopen[bnei] == true)
				{
					// std::cout << "C::" << "|" <<  node_outlet.node << "|" << next.node  << "|" << basopen[bnei] << "|" << std::to_string(bnei == tbas)<< std::endl;
					// std::cout << " rerouting: link is " ; connector->debug_print_row_col(node_outlet.node); std::cout << " to ";
					// connector->debug_print_row_col(next.node);
					int lnode = reroute(connector, Sreceivers, next.node, node_outlet.node);
					// std::cout << "D" << std::endl;
					// basopen[tbas] = true;
					std::queue<int> basb;

					basb.emplace(tbas);

					while(basb.empty() == false)
					{
						int nextbas = basb.front(); basb.pop();
						basopen[nextbas] = true;

						for(auto& link:bas2passes[nextbas])
						{
							// if(baslab[link[0]] == baslab[link[1]])
							// 	std::cout << "HAPPENS::" << baslab[link[0]]  << "|" << baslab[link[1]]  << std::endl;

							if(basopen[baslab[link[1]]] == true)
								continue;

							// std::cout << "C1::" << "|" <<  baslab[link[0]] << "|" << next.node << std::endl;
							lnode = reroute(connector, Sreceivers, link[1],link[0]);
							// std::cout << "D" << std::endl;
							basb.emplace(baslab[lnode]);
						}
					}

				}
				else
				{
					// if(tbas == bnei)
					// 	throw std::runtime_error("???");

					// if(baslab[next.node] == baslab[node_outlet.node])
					// 	throw std::runtime_error("???");

					std::cout << "Link::" << tbas << "|" << bnei << std::endl;

					bas2passes[tbas].emplace_back(std::array<int,2>{next.node,node_outlet.node});
					bas2passes[bnei].emplace_back(std::array<int,2>{node_outlet.node,next.node});
				}



				
			}

			// for(int j = 0; j<nn; ++j)
			// {
			// 	int nei = neighbours[j];
				
			// }	

			// double continue;		
			// cnt:;
		}
		
		int nopen = 0, nclose = 0;
		for(int i=0; i<basopen.size(); ++i)
		{
			if(basopen[i] == false)
			{
				++nclose;
				std::cout << "UNREROUTED BASIN ::" <<i << std::endl;
			}
			else
			{
				++nopen;
			}
		}

		std::cout << "yolo::" << nclose << " vs " << nopen << std::endl;
		return true;



	};



















}
























































#endif