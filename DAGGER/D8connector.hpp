//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef D8connector_HPP
#define D8connector_HPP

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
#include<stdlib.h>
#include<ctime>


// local includes 
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"


// defines all the format_input depnding on the eventual wrapper
#ifdef DAGGER_FT_PYTHON
#include "wrap_helper_python.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#else
#include "wrap_helper_cpp.hpp"
#endif




namespace DAGGER
{

// Useful namespace for my priority queue
using PQ_i_d =  std::priority_queue< PQ_helper<int,double>, std::vector<PQ_helper<int,double> >, std::greater<PQ_helper<int,double> > >;



template<class T>
class D8connector
{
public:

		// General informations about the graph
	// #-> number of nodes (integer and unsized int for loops or list innit).
	int nnodes = 0;
	size_t nnodes_t = 0;
	// #-> number of nodes in x direction.
	int nx = 0;
	// #-> number of nodes in x direction.
	int ny = 0;
	// #-> length in the x direction.
	T dx;
	// #-> length in the y direction.
	T dy;
	T dxy;
	// #-> cell area
	T cellarea;
	T Xmin;
	T Ymin;
	T Xmax;
	T Ymax;
	int not_a_node;

	// Bearer of values
	// #->boundary: index: node index, value: boundary value
	// #-->The possible values are:
	// #---> -1 = no in, no out, not process, internal or external
	// #--->  0 = in, no out, all fluxes leave the system, internal or external
	// #--->  1 = "normal" cell, in, out, internal
	// #--->  2 = external periodic boundary
	std::vector<int> boundary;

	// Helpers for neighbouring operations
	// -> neighbourer holdsthe indices to loop through for each boundary condition
	std::vector<std::vector<int> > neighbourer;
	// -> lengthener is the dx on each directions
	std::vector<T > lengthener;

	// Coordinate stuff
	// Xs and Ys are vectors of nx and ny size converting row to y and col to X
	// Extents holds the cxmin,xmax,ymin,ymax (extent is an option in matplotlib imshow plots)
	std::vector<T> Xs,Ys, extents;

	easyRand randu;

	D8connector(){};

	D8connector(int nx, int ny, T dx, T dy, T xmin, T ymin)
	{
		int nnodes = nx * ny;
		this->set_dimensions(nx,ny,nnodes,dx,dy,xmin,ymin);
		this->set_default_boundaries("4edges");
	}
	
	/// @Description: Initialising the "neighbourer", a data structure managing the iterations though the neighbours of a given node
	/// Function of boundary conditions, one can feed the neighbourer with an indice which returns the index to ADD to the node to get its neighbour.
	/// The neighbour order is (in table referential, be careful if Y axis is inverted) top-left, top, top-right, left, right, bottom-left, bottom, bottom-right
	/// @Authors: B.G.
	/// @Date: 2021
	void initialise_neighbourer()
	{
		T diag = std::sqrt(std::pow(dx,2) + std::pow(dy,2) );
		this->dxy = diag;

		this->lengthener = std::initializer_list<T>{diag,dy,diag,dx,dx,diag,dy,diag};
		this->neighbourer.clear();

		// these vectors are additioned to the node indice to test the neighbors
		this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - this->nx + 1, -1,1,this->nx - 1, this->nx, this->nx + 1 }); // internal node 0
		this->neighbourer.emplace_back(std::initializer_list<int>{(this->ny - 1) * this->nx - 1, (this->ny - 1) * this->nx, (this->ny - 1) * this->nx + 1, -1,1,this->nx - 1, this->nx, this->nx + 1 });// periodic_first_row 1
		this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - this->nx + 1, -1,1,- (this->ny - 1) * this->nx - 1, - (this->ny - 1) * this->nx, - (this->ny - 1) * this->nx + 1 });// periodic_last_row 2
		this->neighbourer.emplace_back(std::initializer_list<int>{- 1, - this->nx, - this->nx + 1, (this->nx - 1),1, 2 * this->nx - 1, this->nx, this->nx + 1 }); // periodic_first_col 3
		this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - 2 * this->nx + 1, -1,-this->nx + 1, this->nx - 1, this->nx, 1 }); // periodic last_col 4
		this->neighbourer.emplace_back(std::initializer_list<int>{this->not_a_node, this->not_a_node, this->not_a_node, -1, 1, this->nx - 1, this->nx, this->nx + 1 }); // normal_first_row 5
		this->neighbourer.emplace_back(std::initializer_list<int>{- this->nx - 1, - this->nx, - this->nx + 1, -1,1, this->not_a_node, this->not_a_node, this->not_a_node}); // normal_last_row 6
		this->neighbourer.emplace_back(std::initializer_list<int>{this->not_a_node, - this->nx, - this->nx + 1, this->not_a_node, 1, this->not_a_node,  this->nx, this->nx + 1 }); // normal_first_col 7
		this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, this->not_a_node, -1, this->not_a_node, this->nx - 1, this->nx, this->not_a_node }); // normal_last_col 8
		this->neighbourer.emplace_back(std::initializer_list<int>{this->not_a_node, this->not_a_node, this->not_a_node, this->not_a_node, 1, this->not_a_node, this->nx, this->nx + 1 }); // normal_top_left 9
		this->neighbourer.emplace_back(std::initializer_list<int>{ this->not_a_node, this->not_a_node, this->not_a_node, -1, this->not_a_node, this->nx - 1, this->nx, this->not_a_node}); // normal_top_right 10
		this->neighbourer.emplace_back(std::initializer_list<int>{ this->not_a_node,- this->nx, - this->nx + 1,this->not_a_node, 1, this->not_a_node, this->not_a_node, this->not_a_node}); // normal_bottom_left 11
		this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, this->not_a_node, -1, this->not_a_node, this->not_a_node, this->not_a_node, this->not_a_node}); // normal_bottom_right 12
		this->neighbourer.emplace_back(std::initializer_list<int>{this->ny * this->nx -1, (this->ny - 1) * this->nx, (this->ny - 1) * this->nx + 1, this->nx - 1, 1, 2 * this->nx - 1, this->nx, this->nx+1}); // top_left_periodic 13
		this->neighbourer.emplace_back(std::initializer_list<int>{(this->ny - 1) * this->nx - 1, (this->ny - 1) * this->nx, (this->ny - 2) * this->nx + 1, -1, - this->nx + 1,this->nx - 1, this->nx, 1 }); // top_right_periodic 14
		this->neighbourer.emplace_back(std::initializer_list<int>{- 1, - this->nx, - this->nx + 1, this->nx - 1,1,- (this->ny - 2) * this->nx - 1, - (this->ny - 1) * this->nx, - (this->ny - 1) * this->nx + 1 });// periodic_bottom_left 15
		this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - this->nx + 1, -1, 1 - this->nx + 1, - (this->ny - 1) * this->nx - 1, - (this->ny - 1) * this->nx, - (this->ny) * this->nx + 1 });// periodic_bottom_right 16
	
	}




	/// @description: Sets the dimension of the graph and initialise the neighbourer
	void set_dimensions(int nx, int ny, int nnodes, T dx, T dy, T xmin, T ymin)
	{
		this->nx = nx;
		this->ny = ny;
		this->nnodes = nnodes;
		this->nnodes_t = size_t(nnodes);
		this->dx = dx;
		this->dy = dy;
		this->cellarea = this->dx * this->dy;
		this->Xmin = xmin;
		this->Ymin = ymin;

		// Not a node is utilised to detect when a neighbouring operation returns not a node
		this->not_a_node = - nx * ny - 10;

		this->initialise_neighbourer();

		// Initialise coordinate stuff
		this->Xs = std::vector<T>(this->nx);
		this->Ys = std::vector<T>(this->ny);
		// this->extents = std::vector<T>(4,0);

		for(int i=0; i<this->nx; ++i)
			this->Xs[i] = this->Xmin + this->dx/2 + i * this->dx;
		for(int i=this->ny -1; i>=0; --i)
		{
			this->Ys[i] = this->Ymin + this->dy/2 + i * this->dy;
		}

		this->Xmax = this->Xs.back() + this->dx/2;
		this->Ymax = this->Ys.back() + this->dy/2;

		std::reverse(this->Ys.begin(), this->Ys.end());

		this->extents = {this->Xmin, this->Xmin +  (this->nx + 1) * this->dx, this->Ymin, this->Ymin + (this->ny + 1) * this->dy};

	}

	void set_default_boundaries(std::string bountype)
	{

		this->boundary = std::vector<int>(this->nnodes_t,1);

		if(bountype == "4edges")
		{
			for(size_t i = 0; i < this->nnodes_t; ++i)
			{
				if(this->is_on_dem_edge(i))
					this->boundary[i] = 0;
			}
		}
		else if(bountype == "periodic_EW")
		{
			for(int i = 0; i < this->nnodes; ++i)
			{
				if(this->is_on_top_row(i) || this->is_on_bottom_row(i))
					this->boundary[i] = 0;
				else if(this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
					this->boundary[i] = 2;
			}
		}
		else if(bountype == "periodic_NS")
		{
			for(int i = 0; i < this->nnodes; ++i)
			{
				if(this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
					this->boundary[i] = 0;
				else if(this->is_on_top_row(i) || this->is_on_bottom_row(i))
					this->boundary[i] = 2;
			}
		}
		else
		{
			throw std::runtime_error("invalid periodic boundaries");
		}
	}

	template<class bou_t>
	void set_custom_boundaries(bou_t& tbound)
	{
		auto bound = format_input(tbound);
		this->boundary = to_vec(bound);
	}

	int get_boundary_at_node(int i){return this->boundary[i];}


	

	bool is_in_bound(int i){return (i>=0 && i<this->nnodes)? true:false;}

	void fill_linknodes_SMG(std::vector<int>& linknodes)
	{
		for(int i=0; i<this->nnodes; ++i)
		{
			int o = this->get_right_index(i); 
			linknodes[i*8] = (this->is_in_bound(o) ? i:-1);
			linknodes[i*8 + 1] = (this->is_in_bound(o) ? o:-1);
			o = this->get_bottomright_index(i); 
			linknodes[i*8 + 2] = (this->is_in_bound(o) ? i:-1);
			linknodes[i*8 + 3] = (this->is_in_bound(o) ? o:-1);
			o = this->get_bottom_index(i); 
			linknodes[i*8 + 4] = (this->is_in_bound(o) ? i:-1);
			linknodes[i*8 + 5] = (this->is_in_bound(o) ? o:-1);
			o = this->get_bottomleft_index(i); 
			linknodes[i*8 + 6] = (this->is_in_bound(o) ? i:-1);
			linknodes[i*8 + 7] = (this->is_in_bound(o) ? o:-1);
		}
	}


	// This function feed is links bool vector of a SMGraph
	template<class topo_t>
	void build_smgraph_only_MF(topo_t& topography, std::vector<bool>& links)
	{		
		// looping thorugh row col to being able to deal with BC
		for(int row = 0; row < this->ny; ++row)
		{
			for(int col = 0; col < this->nx; ++col)
			{

				// Getting node ID
				int i = row * this->nx + col;

				// cannot be a neighbour anyway, abort
				if(this->can_flow_even_go_there(i) == false)
				{
					continue;
				}

				double this_topo = topography[i];

				int n = this->get_id_right_SMG(i);
				int tn = this->get_right_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}
				n = this->get_id_bottomright_SMG(i);
				tn = this->get_bottomright_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}
				n = this->get_id_bottom_SMG(i);
				tn = this->get_bottom_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}
				n = this->get_id_bottomleft_SMG(i);
				tn = this->get_bottomleft_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}

			}
		}

	}

	
	// This function feed is links bool vector of a SMGraph
	template<class topo_t>
	void buildSrecfromlinks(topo_t& topography, std::vector<bool>& links, std::vector<int>& Sreceivers, std::vector<double>& SS)
	{		

		// looping thorugh row col to being able to deal with BC
		for(int row = 0; row < this->ny; ++row)
		{
			for(int col = 0; col < this->nx; ++col)
			{
				// Getting node ID
				int i = row * this->nx + col;
				// cannot be a neighbour anyway, abort
				if(this->can_flow_even_go_there(i) == false)
				{
					continue;
				}

				double this_topo = topography[i];

				int n = this->get_id_right_SMG(i);
				int tn = this->get_right_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dx;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_bottomright_SMG(i);
				tn = this->get_bottomright_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dxy;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_bottom_SMG(i);
				tn = this->get_bottom_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dy;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_bottomleft_SMG(i);
				tn = this->get_bottomleft_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dxy;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
			}
		}
	}

	// This function feed is links bool vector of a SMGraph
	template<class topo_t>
	void build_smgraph_only_MF_mask(topo_t& topography, std::vector<bool>& links, std::vector<bool>& mask)
	{		

		

		// looping thorugh row col to being able to deal with BC
		for(int row = 0; row < this->ny; ++row)
		{
			for(int col = 0; col < this->nx; ++col)
			{

				// Getting node ID
				int i = row * this->nx + col;

				if(mask[i] == false)
					continue;

				// cannot be a neighbour anyway, abort
				if(this->can_flow_even_go_there(i) == false)
				{
					continue;
				}

				double this_topo = topography[i];

				int n = this->get_id_right_SMG(i);
				int tn = this->get_right_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}
				n = this->get_id_bottomright_SMG(i);
				tn = this->get_bottomright_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}
				n = this->get_id_bottom_SMG(i);
				tn = this->get_bottom_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}
				n = this->get_id_bottomleft_SMG(i);
				tn = this->get_bottomleft_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = true;
					else
						links[n] = false;
				}

				n = this->get_id_left_SMG(i);
				tn = this->get_left_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = false;
					else
						links[n] = true;
				}
				n = this->get_id_topleft_SMG(i);
				tn = this->get_topleft_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = false;
					else
						links[n] = true;
				}
				n = this->get_id_top_SMG(i);
				tn = this->get_top_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = false;
					else
						links[n] = true;
				}
				n = this->get_id_topright_SMG(i);
				tn = this->get_topright_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					if(topography[tn] < this_topo)
						links[n] = false;
					else
						links[n] = true;
				}

			}
		}

	}

	// This function feed is links bool vector of a SMGraph
	template<class topo_t>
	void buildSrecfromlinks_mask(topo_t& topography, std::vector<bool>& links, std::vector<int>& Sreceivers, std::vector<double>& SS, std::vector<bool>& mask)
	{		

		// looping thorugh row col to being able to deal with BC
		for(int row = 0; row < this->ny; ++row)
		{
			for(int col = 0; col < this->nx; ++col)
			{
				// Getting node ID
				int i = row * this->nx + col;
				// cannot be a neighbour anyway, abort
				if(this->can_flow_even_go_there(i) == false || mask[i] == false)
				{
					continue;
				}

				double this_topo = topography[i];

				int n = this->get_id_right_SMG(i);
				int tn = this->get_right_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dx;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_bottomright_SMG(i);
				tn = this->get_bottomright_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dxy;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_bottom_SMG(i);
				tn = this->get_bottom_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dy;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_bottomleft_SMG(i);
				tn = this->get_bottomleft_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dxy;
					if(links[n])
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}

				n = this->get_id_left_SMG(i);
				tn = this->get_left_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dx;
					if(links[n] == false)
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_topleft_SMG(i);
				tn = this->get_topleft_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dxy;
					if(links[n] == false)
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_top_SMG(i);
				tn = this->get_top_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dy;
					if(links[n] == false)
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
				n = this->get_id_topright_SMG(i);
				tn = this->get_topright_index(i);
				if(n>=0 && n< this->nnodes * 4 && tn >=0 && tn < this->nnodes)
				{
					double slope = std::abs(this_topo - topography[tn])/this->dxy;
					if(links[n] == false)
					{
						if(SS[i] < slope)
						{
							SS[i] = slope;
							Sreceivers[i] = tn;
						}
					}
					else
					{
						if(SS[tn] < slope)
						{
							SS[tn] = slope;
							Sreceivers[tn] = i;
						}
					}
				}
			}
		}
	}


	// std::vector<int> get_

	int get_id_right_SMG(int i){return i*4;}
	int get_id_bottomright_SMG(int i){return i*4 + 1;}
	int get_id_bottom_SMG(int i){return i*4 + 2;}
	int get_id_bottomleft_SMG(int i){return i*4 + 3;}
	int get_id_left_SMG(int i){int yolo = this->get_left_index(i); if (yolo>=0){ return this->get_id_right_SMG(yolo);} else{ return this->not_a_node;}}
	int get_id_topleft_SMG(int i){int yolo = this->get_topleft_index(i); if (yolo>=0){ return this->get_id_bottomright_SMG(yolo);} else{ return this->not_a_node;}}
	int get_id_top_SMG(int i){int yolo = this->get_top_index(i); if (yolo>=0){ return this->get_id_bottom_SMG(yolo);} else{ return this->not_a_node;}}
	int get_id_topright_SMG(int i){int yolo = this->get_topright_index(i); if (yolo>=0){ return this->get_id_bottomleft_SMG(yolo);} else{ return this->not_a_node;}}

	std::vector<int> get_neighourer_indices_SMG(int i)
	{
		return std::vector<int>{this->get_id_right_SMG(i),this->get_id_bottomright_SMG(i),this->get_id_bottom_SMG(i),this->get_id_bottomleft_SMG(i),this->get_id_left_SMG(i),this->get_id_topleft_SMG(i),this->get_id_top_SMG(i),this->get_id_topright_SMG};
	}

		// ------------------------------------------------

	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   Conversion METHODS        /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to geometrical conversions
	// ------------------------------------------------
	// WILL NEED SOME WORK HERE DEPENDING ON THE GRAPH ORIENTATION AND ALL

	inline void rowcol_from_node_id(int node_index, int& row, int& col)
	{
		col = node_index % this->nx;
		row = (int)std::floor(node_index/this->nx);
	}

	inline int nodeid_from_row_col(int row, int col)
	{
		return row * this->nx + col;
	}

	inline int nodeid_from_XY(T X, T Y)
	{
		int col = floor((X - this->Xmin)/this->dx);
		int row = this->ny - floor((Y - this->Ymin)/this->dy);
		return this->nodeid_from_row_col(row,col);
	}
	// void rowcol_from_XY()


	inline void XY_from_nodeid(int node_index, T& tX, T& tY)
	{
		int row,col;
		this->rowcol_from_node_id(node_index, row, col);
		this->XY_from_rowcol(row, col, tX, tY);
	}

	inline void XY_from_rowcol(int row, int col, T& tX, T& tY)
	{
		tX = this->Xs[col];
		tY = this->Ys[row];
	}

	inline T get_X_from_i(int i)
	{
		int col = i % this->nx;
		return this->Xs[col];
	}

	inline T get_Y_from_i(int i)
	{
		int row = (int)std::floor(i/this->nx);
		return this->Ys[row];
	}


	// ------------------------------------------------

	//	                             	            __
	//                                             / _)
	//                                    _.----._/ /
	//   Neighbouring METHODS            /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to accessing and calculating neighbours
	// ------------------------------------------------

	std::vector<Neighbour<int,T> > get_neighbours(int i, bool ignore_nodata = false)
	{

		// preformatting the output
		std::vector<Neighbour<int,T> > neighbours;

		// Reserving size depending on the flow topology
		neighbours.reserve(8);

		size_t id_neighbourer = this->_get_neighbourer_id(i);

		// Now I need to determine the index of the neighbourer vector
		// Which provides adder to determine the neighbours f(boundary_conditions)
		// internal node 0
		// periodic_first_row 1
		// periodic_last_row 2
		// periodic_first_col 3
		// periodic last_col 4
		// normal_first_row 5
		// normal_last_row 6
		// normal_first_col 7
		// normal_last_col 8
		// normal_top_left 9
		// normal_top_right 10
		// normal_bottom_left 11
		// normal_bottom_right 12
		// top_left_periodic 13
		// top_right_periodic 14
		// periodic_bottom_left 15
		// periodic_bottom_right 16

	
			this->set_neighbours_in_vector(neighbours, i, id_neighbourer, ignore_nodata);

		// and return it
		return neighbours;
	}


	// This function is the helper of the neighbouring function: it fills the neighbours vector with the actul values.
	void set_neighbours_in_vector(std::vector<Neighbour<int,T> >& neighbours, int& i, size_t& id_neighbourer, bool ignore_nodata)
	{
		// node index of the current neighbour
		int tn;

		// If you wonder why I am not iterating with a vector and everything here, it is for the small perf gain of doing it this way
		tn = this->neighbourer[id_neighbourer][1];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[1]));
		tn = this->neighbourer[id_neighbourer][3];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[3]));
		tn = this->neighbourer[id_neighbourer][4];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[4]));
		tn = this->neighbourer[id_neighbourer][6];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[6]));
		tn = this->neighbourer[id_neighbourer][0];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[0]));
		tn = this->neighbourer[id_neighbourer][2];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[2]));
		tn = this->neighbourer[id_neighbourer][5];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[5]));
		tn = this->neighbourer[id_neighbourer][7];
		if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
			neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[7]));				
	
	}

	inline size_t _get_neighbourer_id(int i)
	{

		size_t id_neighbourer = -1;
		// internal node, so neighbourer is 0
		if(this->boundary[i] == 1)
			id_neighbourer = 0;
		else
		{
			// Case top left corner
			if(i == 0)
			{
				// Case Periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 13;
				// Case Noraml
				else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
					id_neighbourer = 9;
			}
			// Case top right corner 
			else if (i == this->nx - 1)
			{
				// Case Periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 14;
				// Case Noraml
				else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
					id_neighbourer = 10;
			}
			// Case bottom left corner 
			else if (i == (this->nnodes - this->nx) )
			{
				// Case Periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 15;
				// Case Noraml
				else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
					id_neighbourer = 11;
			}
			// Case bottom right corner 
			else if (i == (this->nnodes - 1) )
			{
				// Case Periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 16;
				// Case Noraml
				else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
					id_neighbourer = 12;
			}
			// Cases first row (no corner)
			else if(i < this->nx - 1)
			{
				// Case periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 1;
				// Case normal
				else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
					id_neighbourer = 5;
			}
			// Cases last row (no corner)
			else if(i > (this->ny - 1) * this->nx)
			{
				// Case periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 2;
				// Case normal
				else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
					id_neighbourer = 6;
			}
			// Cases first col (no corner)
			else if(i % this->nx == 0)
			{
				// Case periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 3;
				// Case normal
				else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
					id_neighbourer = 7;
			}
			// Cases last col (no corner)
			else if(i % this->nx == this->nx - 1)
			{
				// Case periodic
				if(this->boundary[i] == 2)
					id_neighbourer = 4;
				// Case normal
				else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
					id_neighbourer = 8;
			}
			else
				id_neighbourer = 0;
		}

		// if(id_neighbourer == -1)
		// {
		// 	int row,col;
		// 	this->rowcol_from_node_id(i,row,col);
		// 	std::cout << "boundary was " << this->boundary[i] << " i: " << i << "/" << this->nnodes << " row: " << row << " col: " << col << std::endl;
		// 	throw std::runtime_error("neighbouring issue");
		// }

		return id_neighbourer;
	}

	// D4 single neighbour routines
	Neighbour<int,T> get_left_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][3] + i, this->dx);
	}

	Neighbour<int,T> get_top_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][1] + i, this->dy);
	}

	Neighbour<int,T> get_right_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][4] + i, this->dx);
	}

	Neighbour<int,T> get_bottom_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][6] + i, this->dy);
	}

	// D8 single neighbour extra routines
	Neighbour<int,T> get_topleft_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][0] + i, this->dxy);
	}

	Neighbour<int,T> get_topright_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][2] + i, this->dxy);
	}

	Neighbour<int,T> get_bottomleft_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][5] + i, this->dxy);
	}

	Neighbour<int,T> get_bottomright_neighbour(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return Neighbour<int,T>(this->neighbourer[id_neighbourer][7] + i, this->dxy);
	}

	// D4 single neighbour routines
	int get_left_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][3] + i;
	}

	int get_top_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][1] + i;
	}

	int get_right_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][4] + i;
	}

	int get_bottom_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][6] + i;
	}

	// D8 single neighbour extra routines
	int get_topleft_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][0] + i;
	}

	int get_topright_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][2] + i;
	}

	int get_bottomleft_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][5] + i;
	}

	int get_bottomright_index(int i)
	{
		size_t id_neighbourer = this->_get_neighbourer_id(i);	
		return this->neighbourer[id_neighbourer][7] + i;
	}

	std::vector<int> get_D4_neighbours_only_id(int i)
	{
		std::vector<int> neighs;neighs.reserve(4);
		int tn = this->get_left_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_top_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_right_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_bottom_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		return neighs;
	}

	std::vector<int> get_neighbours_only_id(int i)
	{
		std::vector<int> neighs;neighs.reserve(8);
		int tn = this->get_left_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_top_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_right_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_bottom_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_topleft_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_topright_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_bottomright_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);
		tn = this->get_bottomleft_index(i);
		if(tn >= 0 )
			neighs.emplace_back(tn);


		return neighs;
	}


	// Some of my algorithm are adapted from richdem and require iterating through neighbours the way they do
	// starting from left and going clockwise around the node
	std::vector<int> get_neighbours_only_id_richdemlike(int i)
	{
		std::vector<int> neighs;neighs.reserve(8);
		neighs.emplace_back(this->get_left_index(i));
		neighs.emplace_back(this->get_topleft_index(i));
		neighs.emplace_back(this->get_top_index(i));
		neighs.emplace_back(this->get_topright_index(i));
		neighs.emplace_back(this->get_right_index(i));
		neighs.emplace_back(this->get_bottomright_index(i));
		neighs.emplace_back(this->get_bottom_index(i));
		neighs.emplace_back(this->get_bottomleft_index(i));
		return neighs;
	}





	// Method to test whether a node can outlet flow OUT of the model
	inline bool can_flow_out_there(int i)
	{
		if(this->boundary[i] == 0)
			return true;
		else
		{
			return false;
		}
	}

	// method to check if a node can even accept flow
	inline bool can_flow_even_go_there(int i)
	{
		if(this->boundary[i] < 0)
			return false;
		else
		{
			return true;
		}
	}

	inline bool is_active(int i)
	{
		if(this->can_flow_out_there(i))
			return false;
		if(this->can_flow_even_go_there(i))
			return true;
		return false;
	}

	// Returns true if the node is on the dem edge
	// (i.e. one of the 4 lines)
	bool is_on_dem_edge(int i)
	{
		if( this->is_on_top_row(i) || this->is_on_bottom_row(i) || this->is_on_leftest_col(i) || this->is_on_rightest_col(i) )
			return true;
		return false;
	}

	bool is_on_top_row(int i)
	{
		if(i < this->nx)
			return true;
		return false;
	}

	bool is_on_bottom_row(int i)
	{
		if(i >= this->nnodes - this->nx)
			return true;
		else
			return false;
	}

	bool is_on_leftest_col(int i)
	{
		if(i % this->nx == 0)
			return true;
		else
			return false;
	}

	bool is_on_rightest_col(int i)
	{
		if(i % this->nx == this->nx - 1)
			return true;
		else
			return false;
	}

	void print_dim()
	{
		std::cout << "nx:" << this->nx << std::endl;
		std::cout << "ny:" << this->ny << std::endl;
		std::cout << "nnodes:" << this->nnodes << std::endl;
		std::cout << "dx:" << this->dx << std::endl;
		std::cout << "dy:" << this->dy << std::endl;
	}


	template<class topo_t, class out_t>
	out_t get_HS(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		double altitude = 45;
		double azimuth = 315;
		double z_factor = 1;
		double pi = 3.1415926;

		// std::vector<double> hillshade(ptr, ptr + this->nnodes);
		std::vector<double> hillshade(this->nnodes,0.);


		//convert zenith and azimuth into radians for calculation
		double zenith_rad = (90 - altitude) * pi / 180.0;
		double azimuth_math = 360-azimuth + 90;
		if (azimuth_math >= 360.0) azimuth_math = azimuth_math - 360;
		double azimuth_rad = azimuth_math * pi /180.0;

		for(int i = 0; i<this->nnodes; ++i)
		{
			// Ignoring no data
			if(this->boundary[i] < 0)
				continue;

			double slope_rad = 0;
			double aspect_rad = 0;
			double dzdx = 0;
			double dzdy = 0;

			double ij = topography[i];
			double ijp1 = topography[this->get_right_neighbour(i).node];
			double ip1j = topography[this->get_bottom_neighbour(i).node];
			double ip1jp1 = topography[this->get_bottomright_neighbour(i).node];
			double im1jm1 = topography[this->get_topleft_neighbour(i).node];
			double im1j = topography[this->get_top_neighbour(i).node];
			double im1jp1 = topography[this->get_topright_neighbour(i).node];
			double ijm1 = topography[this->get_left_neighbour(i).node];
			double ip1jm1 = topography[this->get_bottomleft_neighbour(i).node];


			if (ij > 0 )
			{
				dzdx = ((ijp1 + 2*ip1j + ip1jp1) - (im1jm1 + 2*im1j + im1jp1)) / (8 * this->dx);
				dzdy = ((im1jp1 + 2*ijp1 + ip1jp1) - (im1jm1 + 2*ijm1 + ip1jm1)) / (8 * this->dy);
				slope_rad = atan(z_factor * sqrt((dzdx*dzdx) + (dzdy*dzdy)));
				if (dzdx != 0)
				{
					aspect_rad = std::atan2(dzdy, (dzdx*-1));
					if (aspect_rad < 0) aspect_rad = 2*pi + aspect_rad;
				}
				else
				{
					if (dzdy > 0) aspect_rad = pi/2;
					else if (dzdy < 0) aspect_rad = 2 * pi - pi/2;
				}

				hillshade[i] = ((std::cos(zenith_rad) * std::cos(slope_rad)) + (std::sin(zenith_rad) * std::sin(slope_rad) * std::cos(azimuth_rad - aspect_rad)));
				// std::cout << hillshade[i] << "|";
				if (hillshade[i] < 0) hillshade[i] = 0;
			}

		}

		return format_output(hillshade);

	}

	// Largely adapted from RichDEM
	template <class topo_t>
	std::vector<double> fill_barne_2014(topo_t& ttopography)
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

	  std::vector<bool> closed(this->nnodes,false);

	  // RDLOG_PROGRESS<<"Adding cells to the priority queue...";
	  for(int i = 0; i<this->nnodes; ++i)
	  {
	  	if(this->can_flow_out_there(i))
	  	{
	  		closed[i] = true; 
	  		open.emplace(i,topography[i]);
	  	}
	  }

	  // RDLOG_PROGRESS<<"Performing Priority-Flood+Epsilon...";

	  while(open.size()>0 || pit.size()>0)
	  {
	    PQ_helper<int,double> c;

	    if(pit.size()>0 && open.size()>0 && open.top().score==pit.front().score)
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

	    auto neighbours = this->get_neighbours_only_id(c.node);
	    for(auto n:neighbours)
	    {
	      if(this->is_active(n) == false) 
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
	        topography[n] = c.score +1e-3 + distr(gen);
	        pit.emplace(n,topography[n]);
	      } else
	        open.emplace(n,topography[n]);
	    }
	  }
	
	  return topography;
	}


	template<class topo_t>
	void InitPriorityQue(
	  topo_t& topography,
	  std::vector<bool>& flag,
	  PQ_i_d& priorityQueue
	)
	{

	  std::queue<PQ_helper<int,double> > depressionQue;

	  // push border cells into the PQ
	  for(int i=0; i < this->nnodes; ++i)
	  {
	    if (flag[i]) continue;

	    if (this->can_flow_even_go_there(i) == false) 
	    {
	      flag[i]=true;
	      auto neighs = this->get_neighbours_only_id(i);
	      for (auto n:neighs)
	      {
	        if (flag[n])
	          continue;
	        if (this->can_flow_even_go_there(n))
	        {
	          priorityQueue.emplace(n,topography[n]);
	          flag[n]=true;
	        }
	      }
	    } 
	    else if(this->can_flow_out_there(i))
	    {
	      //on the DEM border
	      priorityQueue.emplace(i,topography[i]);
	      flag[i]=true;
	    }
	  }
	}



	template<class topo_t>
	void ProcessTraceQue(
	  topo_t& topography,
	  std::vector<bool>& flag,
	  std::queue<PQ_helper<int,double> >& traceQueue,
	  PQ_i_d& priorityQueue
	)
	{
	  std::queue<PQ_helper<int,double>  > potentialQueue;
	  int indexThreshold=2;  //index threshold, default to 2
	  while (!traceQueue.empty())
	  {
	    const auto node = traceQueue.front();
	    traceQueue.pop();

	    bool Mask[5][5]={{false},{false},{false},{false},{false}};
	    
	    auto neighs = this->get_neighbours_only_id(node.node);
	    for (auto n: neighs)
	    {

	      if(flag[n])
	        continue;

	      if (topography[n] > node.score)
	      {
	        traceQueue.emplace(n, topography[n]);
	        flag[n]=true;
	      } 
	      else 
	      {
	        //initialize all masks as false
	        bool have_spill_path_or_lower_spill_outlet=false; //whether cell n has a spill path or a lower spill outlet than node if n is a depression cell
	        auto nneighs = this->get_neighbours_only_id_richdemlike(n);
	        int row,col; this->rowcol_from_node_id(node.node,row,col);
	        int nrow,ncol; this->rowcol_from_node_id(n,nrow,ncol);
	        int incr=0;
	        for(auto nn:nneighs)
	        {
	        	++incr;
	        	// Checking node validity
	        	if(nn<0 || nn>= this->nnodes)
	        		continue;

	        	int nnrow,nncol;
	        	this->rowcol_from_node_id(nn,nnrow,nncol);
	          if( (Mask[nnrow-row+2][nncol-col+2]) || (flag[nn] && topography[nn] < node.score) )
	          {
	            Mask[nrow-row+2][ncol-col+2]=true;
	            have_spill_path_or_lower_spill_outlet=true;
	            break;
	          }
	        }

	        if(!have_spill_path_or_lower_spill_outlet)
	        {
	          if (incr<indexThreshold) potentialQueue.push(node);
	          else
	            priorityQueue.push(node);
	          break; // make sure node is not pushed twice into PQ
	        }
	      }
	    }//end of for loop
	  }

	  while (!potentialQueue.empty())
	  {
	    const auto node = potentialQueue.front();
	    potentialQueue.pop();

	    //first case
	    auto neigh = this->get_neighbours_only_id_richdemlike(node.node);
	    int incr = 0;
	    for (auto n:neigh)
	    {
	    	++incr;

	      if(flag[n])
	        continue;

	      priorityQueue.push(node);
	      break;
	    }
	  }
	}



	template<class topo_t>
	void ProcessPit(
	  topo_t& topography,
	  std::vector<bool>& flag,
	  std::queue<PQ_helper<int,double> >& depressionQue,
	  std::queue<PQ_helper<int,double> >& traceQueue
	)
	{
	  while (!depressionQue.empty())
	  {
	    auto node = depressionQue.front();
	    depressionQue.pop();
	    auto neigh = this->get_neighbours_only_id(node.node);

	    for (auto n: neigh)
	    {

	      if (flag[n])
	        continue;

	      const auto iSpill = topography[n];
	      if (iSpill > node.score)
	      { // Dlope cell
	        flag[n]=true;
	        traceQueue.emplace(n,iSpill);
	        continue;
	      }

	      // Depression cell
	      flag[n] = true;
	      topography[n] = node.score + 1e-3 + 1e-6 * this->randu.get();
	      depressionQue.emplace(n,topography[n]);
	    }
	  }
	}



	template<class topo_t>
	std::vector<double> PriorityFlood_Wei2018(topo_t &ttopography)
	{
	  std::queue<PQ_helper<int,double> > traceQueue;
	  std::queue<PQ_helper<int,double> > depressionQue;

	  std::vector<bool> flag(this->nnodes,false);
	  std::vector<double> topography = to_vec(ttopography);

	  PQ_i_d priorityQueue;

	  this->InitPriorityQue(topography,flag,priorityQueue);

	  while (!priorityQueue.empty())
	  {
	    const auto tmpNode = priorityQueue.top();
	    priorityQueue.pop();

	    auto neigh = this->get_neighbours_only_id(tmpNode.node);

	    for(auto n:neigh)
	    {

	      if(flag[n])
	        continue;

	      auto iSpill = topography[n];

	      if (iSpill <= tmpNode.score)
	      {
	        //depression cell
	        topography[n] = tmpNode.score + 1e-3 + 1e-6 * this->randu.get();
	        flag[n] = true;
	        depressionQue.emplace(n,topography[n]);
	        this->ProcessPit(topography,flag,depressionQue,traceQueue);
	      } 
	      else 
	      {
	        //slope cell
	        flag[n] = true;
	        traceQueue.emplace(n,iSpill);
	      }
	      ProcessTraceQue(topography,flag,traceQueue,priorityQueue);
	    }
	  }

	  return topography;
	}


	// this function makes the connector generic and allow other connector to have varying cellarea
	template<class i_t>
	T get_area_at_node(i_t i)
	{return this->cellarea;}

	// template<>
	// void fill_neighbour_matrices()

	template<class i_t>
	T get_dx_from_links_idx( i_t i)
	{
		if(i%4 == 0)
			return this->dx;
		else if(i%4 == 1)
			return this->dxy;
		else if(i%4 == 2)
			return this->dy;
		else if(i%4 == 3)
			return this->dxy;
		else
			return this->dx;
	}

	template<class i_t>
	T get_dx_from_linknodes_idx( i_t i)
	{
		int j = std::floor(i/2);
		return this->get_dx_from_links_idx(j);
	}

	std::vector<int> get_ilinknodes_from_node(int i)
	{
		std::vector<int> out;out.reserve(8);
		int tn = this->get_id_right_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		tn = this->get_id_bottomright_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		tn = this->get_id_bottom_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		tn = this->get_id_bottomleft_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		tn = this->get_id_left_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		tn = this->get_id_topleft_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		tn = this->get_id_top_SMG(i);
		if(tn >0 && tn < this->nnodes * 4)
			out.emplace_back(tn);
		return out;
	}

	std::vector<bool> get_mask_array()
	{
		std::vector<bool> mask(this->nnodes,true);
		for (int i=0; i< this->nnodes; ++i)
		{
			if(this->is_active(i) == false)
				mask[i] = false;
		}
		return mask;
	}


};











// end of namespace
};











































#endif