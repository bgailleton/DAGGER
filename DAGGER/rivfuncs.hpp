#pragma once

#include "RivNet1D.hpp"
#include "databag.hpp"
#include "dodconnector.hpp"

namespace DAGGER {

template<class i_t, class f_t>
void
add_river_network_from_threshold(f_t Ath, Connector8<i_t, f_t>& con)
{

	con.PFcompute_all(false);

	std::vector<f_t> DA(con.nxy(), 0);
	std::vector<std::uint8_t> isDone(con.nxy(), false);

	for (int i = con.nxy() - 1; i >= 0; --i) {
		int node = con.data->_Sstack[i];

		if (nodata(con.data->_boundaries[node]))
			continue;

		DA[node] += con.area(node);

		int rec = con.Sreceivers(node);

		if (rec == node)
			continue;

		DA[rec] += con.data->_DA[node];
	}

	std::vector<i_t> outi;

	for (int i = con.nxy() - 1; i >= 0; --i) {
		int node = con.data->_Sstack[i];

		if (nodata(con.data->_boundaries[node]))
			continue;
		int rec = con.Sreceivers(node);

		if (rec == node)
			continue;

		if (DA[node] > Ath && isDone[rec] == false) {
			int tnode = node;
			int trec = con.Sreceivers(node);
			outi.emplace_back(node);

			while (tnode != trec) {
				isDone[tnode] = true;
				tnode = trec;
				trec = con.Sreceivers(tnode);
			}
		}
	}

	RivNet1D<i_t, f_t> triv();
	triv.build_from_sources(outi);

	this->data->river_networks.emplace_back(triv);
}

}
