#pragma once

#include "boundary_conditions.hpp"
#include "databag.hpp"
#include "dodcontexts.hpp"
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
	f_t _area = 1.;
	// # Total Length in the x dimension
	f_t _lx = 0.;
	// # Total Length in the y dimension
	f_t _ly = 0.;

	f_t area(i_t i) const { return this->_area; }

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
		this->_area = this->_dx * this->_dy;

		// Initialisatio of the lookup table
		this->data->LK8 =
			lookup8<i_t, f_t>(this->_nx, this->_ny, this->_dx, this->_dy);
	}

	// Initialisation of the model
	// The initialisation has to be run after the construction as it depends on
	// the different params of the parameter sheet
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

		} else if (this->boutype == CONBOU::PEW) {
			this->data->_boundaries = std::vector<BC>(this->nxy(), BC::FLOW);
			for (int i = 0; i < this->_ny; ++i) {
				this->data->_boundaries[i * this->_nx] = BC::PERIODIC_BORDER;
				this->data->_boundaries[i * this->_nx + this->_nx - 1] =
					BC::PERIODIC_BORDER;
			}
			for (int i = 0; i < this->_nx; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
				this->data->_boundaries[i] = BC::OUT;

		}

		else
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
	// More efficient than recalling init at every updates as it skips the memory
	// allocation bits
	void reinit()
	{
		if (this->flowtopo == CONFLOWTOPO::ALL) {
			fillvec(this->data->_Sreceivers, static_cast<std::uint8_t>(0));
			fillvec(this->data->_Sdonors, static_cast<std::uint8_t>(0));
			fillvec(this->data->_donors, static_cast<std::uint8_t>(0));
			fillvec(this->data->_receivers, static_cast<std::uint8_t>(0));
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			fillvec(this->data->_donors, static_cast<std::uint8_t>(0));
			fillvec(this->data->_receivers, static_cast<std::uint8_t>(0));
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			fillvec(this->data->_Sreceivers, static_cast<std::uint8_t>(0));
			fillvec(this->data->_Sdonors, static_cast<std::uint8_t>(0));
		}
	}

	i_t Sreceivers(i_t i) const
	{
		return i + this->data->LK8.NeighbourerD8[this->data->LK8.BC2idAdder(
								 i, this->data->_boundaries[i])][this->data->_Sreceivers[i]];
	}

	i_t Sdonors(i_t i, std::array<i_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourer[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_Sdonors[i]];
		i_t nn = this->nSdonors(i);
		for (int j = 0; j < nn; ++j)
			arr[j] += i;
		return nn;
	}

	i_t nSdonors(i_t i) const
	{
		return this->data->LK8.NeighbourerNN[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_Sdonors[i]];
	}

	// one-off computing operation: it computes the neighbours code to loop
	// through the different neighbours of each node
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
					int TN = node + LK8.NeighbourerD8[this->data->LK8.BC2idAdder(
														node, this->data->_boundaries[node])][TTN];
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
		if (boundaries[0] != BC::NO_FLOW) {
			if (boundaries[0] != BC::PERIODIC_BORDER) {
				neighbours[0] = LK8.TopLeft_normal_boundary();
			} else {
				neighbours[0] = AllMask8;
			}
		}

		// Maninging topright corner
		// | | |X|
		// | | | |
		// | | | |
		if (boundaries[this->_nx - 1] != BC::NO_FLOW) {
			if (boundaries[this->_nx - 1] != BC::PERIODIC_BORDER) {
				neighbours[this->_nx - 1] = LK8.TopRight_normal_boundary();
			} else {
				neighbours[this->_nx - 1] = AllMask8;
			}
		}

		// Mananging bottomleft corner
		// | | | |
		// | | | |
		// |X| | |
		if (boundaries[this->nxy() - this->_nx] != BC::NO_FLOW) {
			if (boundaries[this->nxy() - this->_nx] != BC::PERIODIC_BORDER) {
				neighbours[this->nxy() - this->_nx] = LK8.BottomLeft_normal_boundary();
			} else {
				neighbours[this->nxy() - this->_nx] = AllMask8;
			}
		}

		// Mananging bottomright corner
		// | | | |
		// | | | |
		// | | |X|
		if (boundaries[this->nxy() - 1] != BC::NO_FLOW) {
			if (boundaries[this->nxy() - 1] != BC::PERIODIC_BORDER) {
				neighbours[this->nxy() - 1] = LK8.BottomRight_normal_boundary();
			} else {
				neighbours[this->nxy() - 1] = AllMask8;
			}
		}

		// Managing first row
		// | |X| |
		// | | | |
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			// Getting the node index
			int node = r;
			if (boundaries[node] != BC::NO_FLOW)
				if (boundaries[node] != BC::PERIODIC_BORDER)
					neighbours[node] = LK8.Top_normal_boundary();
				else
					neighbours[node] = AllMask8;
		}

		// Managing last row
		// | | | |
		// | | | |
		// | |X| |
		for (int r = this->nxy() - this->_nx + 1; r < this->nxy() - 1; ++r) {
			// Getting the node index
			int node = r;
			if (boundaries[node] != BC::NO_FLOW)
				if (boundaries[node] != BC::PERIODIC_BORDER)
					neighbours[node] = LK8.Bottom_normal_boundary();
				else
					neighbours[node] = AllMask8;
		}

		// Managing Left and right columns
		// | | | |
		// |X| |X|
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			// Getting the node index
			int n1 = r * this->_nx + 0;
			if (boundaries[n1] != BC::NO_FLOW) {
				if (boundaries[n1] != BC::PERIODIC_BORDER)
					neighbours[n1] = LK8.Left_normal_boundary();
				else
					neighbours[n1] = AllMask8;
			}
			n1 = r * this->_nx + this->_nx - 1;
			if (boundaries[n1] != BC::NO_FLOW) {
				if (boundaries[n1] != BC::PERIODIC_BORDER)
					neighbours[n1] = LK8.Right_normal_boundary();
				else
					neighbours[n1] = AllMask8;
			}
		}
	}

	void compute()
	{
		if (this->flowtopo == CONFLOWTOPO::NONE)
			return;

		this->reinit();

		if (this->flowtopo == CONFLOWTOPO::ALL) {
			this->_compute_all();
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			this->_compute_mfd_only();
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			// this->_compute_sfd_only();
		}
	}

	void _compute_all()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		// Setting up context
		CT_neighbourer_1<i_t, f_t> ctx;

		// Setting up prefetchers

		for (int i = 0; i < this->_nxy; ++i) {
			ctx.update(i, *this);
			f_t SS = 0;
			std::uint8_t rcode = 0;
			std::uint8_t srecode = 0;
			std::uint8_t dcode = 0;

			for (int j = 0; j < ctx.nn; ++j) {
				f_t dz = ctx.topo - ctx.neighboursTopo[j];

				if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false)
					continue;

				// First asserting the connectivity
				if (dz < 0)
					dcode |= ctx.neighboursBits[j];

				else if (dz > 0) {
					rcode |= ctx.neighboursBits[j];
					f_t tSS = dz / ctx.neighboursDx[j];
					if (tSS > SS) {
						SS = tSS;
						srecode = ctx.neighboursBits[j];
					}
				}
			}

			this->data->_Sreceivers[ctx.node] = srecode;
			this->data->_receivers[ctx.node] = rcode;
			this->data->_donors[ctx.node] = dcode;
			this->data->_Sdonors[this->Sreceivers(ctx.node)] |= invBits(srecode);
		}
	}

	void _compute_mfd_only()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		// Setting up context
		CT_neighbourer_1<i_t, f_t> ctx;

		// Setting up prefetchers

		for (int i = 0; i < this->_nxy; ++i) {
			ctx.update(i, *this);
			std::uint8_t rcode = 0;
			std::uint8_t dcode = 0;

			for (int j = 0; j < ctx.nn; ++j) {
				f_t dz = ctx.topo - ctx.neighboursTopo[j];

				if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false)
					continue;

				// First asserting the connectivity
				if (dz < 0)
					dcode |= ctx.neighboursBits[j];
				else if (dz > 0) {
					rcode |= ctx.neighboursBits[j];
				}
			}
			this->data->_receivers[ctx.node] = rcode;
			this->data->_donors[ctx.node] = dcode;
		}
	}

	void _quickSstack()
	{
		// The stack container helper
		this->data->_Sstack = std::vector<int>(this->_nxy, 0);

		std::stack<int, std::vector<int>> stackhelper;
		std::array<i_t, 8> tsdon;
		// std::vector<bool> isdone(this->nnodes,false);
		// going through all the nodes
		int istack = 0;
		for (int i = 0; i < this->_nxy; ++i) {
			// if they are base level I include them in the stack
			if (this->Sreceivers(i) == i) {
				stackhelper.emplace(i);
				// ++istack;
			}

			// While I still have stuff in the stack helper
			while (stackhelper.empty() == false) {
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();
				stackhelper.pop();
				// std::cout << istack << "/" << this->_nxy << std::endl;
				this->data->_Sstack[istack] = nextnode;
				++istack;

				// as well as all its donors which will be processed next
				i_t nn = this->Sdonors(nextnode, tsdon);
				for (int j = 0; j < nn; ++j) {
					stackhelper.emplace(tsdon[j]);
				}
			}
		}
	}
};

} // End of Namespace
