#pragma once

#include "boundary_conditions.hpp"
#include "databag.hpp"
#include "enumutils.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, f_t>
class Connector8
{

public:
	// dimensions
	int _nxy = 0;
	int nxy() const { return this->_nxy; }
	int _nx = 0;
	int _ny = 0;
	f_t _dx = 1.;
	f_t _dy = 1.;
	f_t _dxy = std::sqrt(2.);
	f_t _lx = 0.;
	f_t _ly = 0.;

	// Dimensions
	DataBag<i_t, f_t>* data;

	// Param sheet
	CONFLOWTOPO flowtopo = CONFLOWTOPO::ALL;
	CONBOU boutype = CONBOU::EDGES;

	Connector8() { ; }
	Connector8(int nx, int ny, f_t dx, f_t dy, DataBag<i_t, f_t>& data)
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
	}

	void init()
	{

		this->data->LK8 =
			lookup8<class i_t, f_t>(this->_nx, this->_ny, this->_dx, this->_dy);

		if (this->boutype == CONBOU::EDGES) {
			this->data->_boundaries = std::vector<uint8_t>(this->nxy(), BC::FLOW);
			for (int i = 0; i < this->_nx; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = 0; i < this->_ny; ++i) {
				this->data->_boundaries[i * this->_nx] = BC::OUT;
				this->data->_boundaries[i * this->_nx + this->_nx - 1] = BC::OUT;
			}

		} else
			throw std::runtime("boutype NOT IMPLOEMENTED YET");

		this->_compute_neighbours();

		if (this->flowtopo == CONFLOWTOPO::ALL) {
			this->data->_Sreceivers = std::vector<std::uint8>(this->_nxy, 0);
			this->data->_Sdonors = std::vector<std::uint8>(this->_nxy, 0);
			this->data->_donors = std::vector<std::uint8>(this->_nxy, 0);
			this->data->_receivers = std::vector<std::uint8>(this->_nxy, 0);
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			this->data->_donors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_receivers = std::vector<std::uint8_t>(this->_nxy, 0);
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			this->data->_Sreceivers = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_Sdonors = std::vector<std::uint8_t>(this->_nxy, 0);
		}

		// this->reinit();
	}

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

	void _compute_neighbours()
	{

		// Initialising the neighbour array
		this->data->_neighbours = std::vector<std::uint8_t>(this->_nxy, 0);

		// Setting up some refs for code clarity
		std::vector<std::uint8_t>& neighbours = this->data->_neighbours;
		std::vector<std::uint8_t>& boundaries = this->data->_boundaries;
		auto& LK8 = this->data->LK8;

		// Starting with internal nodes
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
		if (boundaries[0] ==)

			// Managing first row
			for (int r = 1; r < this->_ny - 1; ++r) {
			}
	}
};

}
