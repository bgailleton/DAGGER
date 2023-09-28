#pragma once

#include "boundary_conditions.hpp"
#include "databag.hpp"
#include "enumutils.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class Connector8
{

public:
	// Grid dimensions
	// # total number of nodes
	int _nxy = 0;
	int nxy() const { return this->_nxy; }
	// # Number of columns
	int _nx = 0;
	// # Number of rows
	int _ny = 0;
	// # Spacing in the x dimension
	f_t _dx = 1.;
	// # Spacing in the y dimension
	f_t _dy = 1.;
	// # Spacing in the diagonal dimension
	f_t _dxy = std::sqrt(2.);
	// # Total Length in the x dimension
	f_t _lx = 0.;
	// # Total Length in the y dimension
	f_t _ly = 0.;

	// Data storer
	Hermes<i_t, f_t>* data;

	// Param sheet
	// # Flow topology to compute
	CONFLOWTOPO flowtopo = CONFLOWTOPO::ALL;
	// # Boundary conditions
	CONBOU boutype = CONBOU::EDGES;

	// Default Constructor
	Connector8() { ; }
	
	// Constructor
	Connector8(int nx, int ny, f_t dx, f_t dy, Hermes<i_t, f_t>& data)
	{
		this->_nx = nx;
		this->_ny = ny;
		this->_nxy = nx * ny;
		this->_dx = dx;
		this->_dy = dy;
		this->_dxy = std::sqrt(dx * dx + dy * dy);
		this->_lx = (nx + 1) * dx;
		this->_ly = (ny + 1) * dy;
		this->data = &data;

		// Initialisatio of the lookup table
		this->data->LK8 =
			lookup8<i_t, f_t>(this->_nx, this->_ny, this->_dx, this->_dy);

	}


	// Initialisation of the model
	// The initialisation has to be run after the construction as it depends on the different params of the parameter sheet
	void init()
	{

		// Initialising the boundary conditions
		// # If boundaries are EDGES, automatically sets them to 
		if (this->boutype == CONBOU::EDGES) {
			this->data->_boundaries = std::vector<BC>(this->nxy(), BC::FLOW);
			for (int i = 0; i < this->_nx; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = 0; i < this->_ny; ++i) {
				this->data->_boundaries[i * this->_nx] = BC::OUT;
				this->data->_boundaries[i * this->_nx + this->_nx - 1] = BC::OUT;
			}

		} else
			throw std::runtime_error("boutype NOT IMPLOEMENTED YET");

		// Computing the neighbour code
		this->_compute_neighbours();

		// Initialising the data concerned by the different topologies
		if (this->flowtopo == CONFLOWTOPO::ALL) {
			this->data->_Sreceivers = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_Sdonors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_donors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_receivers = std::vector<std::uint8_t>(this->_nxy, 0);
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			this->data->_donors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_receivers = std::vector<std::uint8_t>(this->_nxy, 0);
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			this->data->_Sreceivers = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_Sdonors = std::vector<std::uint8_t>(this->_nxy, 0);
		}

	}

	// Operation reinitialising the data. 
	// More efficient than recalling init at every updates as it skips the memory allocation bits
	void reinit()
	{	
		if (this->flowtopo == CONFLOWTOPO::ALL) {
			fillvec(this->data->_Sreceivers, 0);
			fillvec(this->data->_Sdonors, 0);
			fillvec(this->data->_donors, 0);
			fillvec(this->data->_receivers, 0);
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			fillvec(this->data->_donors, 0);
			fillvec(this->data->_receivers, 0);
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			fillvec(this->data->_Sreceivers, 0);
			fillvec(this->data->_Sdonors, 0);
		}
	}


	// one-off computing operation: it computes the neighbours code to loop through 
	// the different neighbours of each node
	void _compute_neighbours()
	{

		// Initialising the neighbour array
		this->data->_neighbours = std::vector<std::uint8_t>(this->_nxy, 0);

		// Setting up some refs for code clarity
		std::vector<std::uint8_t>& neighbours = this->data->_neighbours;
		std::vector<BC>& boundaries = this->data->_boundaries;
		auto& LK8 = this->data->LK8;

		// Starting with internal nodes
		// | | | |
		// | |X| |
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			for (int c = 1; c < this->_nx - 1; ++c) {
				// Getting the node index
				int node = r * this->_nx + c;

				// Default to nonei
				neighbours[node] = 0;

				// if node is no flow - skip
				if (boundaries[node] == BC::NO_FLOW)
					continue;

				// Going through theoretical neighbours
				for (auto TTN : NeighbourerMask8) {
					// getting the index
					int TN = node + LK8.NeighbourerD8[TTN];
					// if neighbour is no flow ignore it
					if (boundaries[TN] == BC::NO_FLOW)
						continue;
					// Else it's basically a neighbour and that's it
					neighbours[node] |= TTN;
				}
			}
		}

		// Maninging topleft corner
		// |X| | |
		// | | | |
		// | | | |
		if (boundaries[0] != BC::NO_FLOW){
			if(boundaries[0] != BC::PERIODIC_BORDER){
				neighbours[0] = LK8.TopLeft_normal_boundary();
			} 
		}

		// Maninging topright corner
		// | | |X|
		// | | | |
		// | | | |
		if (boundaries[this->_nx-1] != BC::NO_FLOW){
			if(boundaries[this->_nx-1] != BC::PERIODIC_BORDER){
				neighbours[this->_nx-1] = LK8.TopRight_normal_boundary();
			} 
		}

		// Mananging bottomleft corner
		// | | | |
		// | | | |
		// |X| | |
		if (boundaries[this->nxy()-this->_nx] != BC::NO_FLOW){
			if(boundaries[this->nxy()-this->_nx] != BC::PERIODIC_BORDER){
				neighbours[this->nxy()-this->_nx] = LK8.BottomLeft_normal_boundary();
			} 
		}

		// Mananging bottomright corner
		// | | | |
		// | | | |
		// | | |X|
		if (boundaries[this->nxy()-1] != BC::NO_FLOW){
			if(boundaries[this->nxy()-1] != BC::PERIODIC_BORDER){
				neighbours[this->nxy()-1] = LK8.BottomRight_normal_boundary();
			} 
		}

		// Managing first row
		// | |X| |
		// | | | |
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			// Getting the node index
			int node = r;
			if(boundaries[node] != BC::NO_FLOW && boundaries[node] != BC::PERIODIC_BORDER)
				neighbours[node] = LK8.Top_normal_boundary();
		}

		

		// Managing last row
		// | | | |
		// | | | |
		// | |X| |
		for (int r = this->nxy() - this->_nx + 1; r < this->nxy() - 1; ++r) {
			// Getting the node index
			int node = r;
			if(boundaries[node] != BC::NO_FLOW && boundaries[node] != BC::PERIODIC_BORDER)
				neighbours[node] = LK8.Bottom_normal_boundary();
		}
		
		// Managing Left and right columns
		// | | | |
		// |X| |X|
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			// Getting the node index
			int n1 = r*this->_nx+0;
			if(boundaries[n1] != BC::NO_FLOW && boundaries[n1] != BC::PERIODIC_BORDER)
				neighbours[n1] = LK8.Left_normal_boundary();
			n1 = r*this->_nx+this->_nx-1;
			if(boundaries[n1] != BC::NO_FLOW && boundaries[n1] != BC::PERIODIC_BORDER)
				neighbours[n1] = LK8.Right_normal_boundary();
		}


	}

	void compute(){
		if (this->flowtopo == CONFLOWTOPO::NONE) continue;

		if (this->flowtopo == CONFLOWTOPO::ALL) {
			// this->_compute_all();
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			// this->_compute_mfd_only();
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			// this->_compute_sfd_only();
		}


	}

};

}
