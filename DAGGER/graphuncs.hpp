#pragma once

#include "databag.hpp"
#include "dodconnector.hpp"
#include "parambag.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<f_t>
_compute_min_distance_from_outlets(CONNECTOR_T& con)
{

	// preformatting the output
	std::vector<f_t> out(con.nxy(), 0.);

	std::array<f_t, 8> dxs;
	std::array<i_t, 8> recs;

	for (auto node : con.data->_stack) {

		int nd = con.Receivers(node, recs);

		con.ReceiversDx(node, dxs);

		for (int j = 0; j < nd; ++j) {
			out[node] = (out[node] == 0) ? out[recs[j]] + dxs[j]
																	 : std::min(out[node], out[recs[j]] + dxs[j]);
		}
	}

	return out;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
compute_min_distance_from_outlets(bool compute_PQ,
																	std::string name,
																	CONNECTOR_T& con)
{

	if (compute_PQ) {
		con.PFcompute_all(false);
	}

	std::vector<f_t> dx2out =
		_compute_min_distance_from_outlets<i_t, f_t, CONNECTOR_T>(con);

	con.data->fbag[name] = dx2out;
}

}; // end of namespace
