#pragma once
#include "boundary_conditions.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class CT_neighbourer_1
{

public:
	i_t node = 0;
	f_t topo;
	BC boundary;

	i_t nn = 0;
	std::uint8_t idAdder = 0;
	std::array<i_t, 8> neighbours;
	std::array<std::uint8_t, 8> neighboursBits;
	std::array<BC, 8> neighboursCode;
	std::array<f_t, 8> neighboursDx;
	std::array<f_t, 8> neighboursTopo;

	template<class CONNECTOR_T>
	void update(i_t node, CONNECTOR_T& con)
	{

		this->node = node;

		this->topo = con.data->_surface[node];
		this->boundary = con.data->_boundaries[node];

		std::uint8_t nid = con.data->_neighbours[node];

		this->idAdder = con.data->LK8.BC2idAdder(node, con.data->_boundaries[node]);

		this->nn = con.data->LK8.NeighbourerNN[this->idAdder][nid];

		this->neighbours = con.data->LK8.Neighbourer[this->idAdder][nid];

		this->neighboursDx = con.data->LK8.Neighbourerdx[this->idAdder][nid];
		this->neighboursBits = con.data->LK8.NeighbourerBits[this->idAdder][nid];

		for (size_t i = 0; i < this->nn; ++i)
			this->neighbours[i] += node;

		for (size_t i = 0; i < nn; ++i)
			this->neighboursTopo[i] = con.data->_surface[this->neighbours[i]];

		for (size_t i = 0; i < nn; ++i)
			neighboursCode[i] = con.data->_boundaries[this->neighbours[i]];
	}
};

} // end of namespace
