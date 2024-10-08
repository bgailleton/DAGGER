/*
This file is experimental && contains equivalent neighbouring operation from
Riverdale GPU model

*/

#pragma once

#include <array>
#include <functional>
#include <stack>
#include <vector>

#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyvectorize.hpp"

namespace DAGGER {

template<class GRID_T, class F_T>
void
compute_full_SF_graph(GRID_T& grid,
											xt::pytensor<F_T, 2>& topography,
											xt::pytensor<int, 2>& Sreceivers,
											xt::pytensor<int, 2>& Ndonors,
											xt::pytensor<int, 3>& donors,
											xt::pytensor<int, 1>& Stack,
											xt::pytensor<std::uint8_t, 2>& BCs)
{

	// First computing the steepest receivers and filling with 0s
	Ndonors.fill(0);
	Stack.fill(-1);
	donors.fill(-1);

	// Calculating the steepest receivers
	for (int r = 0; r < grid.ny; ++r) {
		for (int c = 0; c < grid.nx; ++c) {

			// By default, no receiver = its own receiver
			Sreceivers(r, c) = grid.rowcol2index(r, c);

			// if cannot receive: skip
			if (grid.can_receive(r, c, BCs, &grid) == false)
				continue;

			// elev of steepest rec
			double tZ = topography(r, c);
			// ID of steepest rec
			int tR = -1;

			// Going through 4 or 8 neighbours
			for (int k = 0; k < grid.nk; ++k) {

				// getting the neighbour
				auto nei = grid.neighbours(r, c, k, BCs, &grid);

				// chekcing the neighbour
				if (nei[0] == -1)
					continue;

				// Can the neighbour eceive flux
				if (grid.can_receive(nei[0], nei[1], BCs, &grid) == false)
					continue;

				// finally is the neighbour lower than the steepest neighbour
				if (topography(nei[0], nei[1]) < tZ) {
					tZ = topography(nei[0], nei[1]);
					tR = grid.rowcol2index(nei[0], nei[1]);
				}
			}

			// if there IS a steepest receivers
			if (tR != -1) {

				// Saving it
				Sreceivers(r, c) = tR;

				// getting the 2D indices
				int rr, cc;
				grid.id2rowcol(tR, rr, cc);

				// Registering the Nth donor
				donors(rr, cc, Ndonors(rr, cc)) = grid.rowcol2index(r, c);

				// incrementing the number of donors
				++Ndonors(rr, cc);
			}
		}
	}

	// The stack container helper
	std::stack<size_t, std::vector<size_t>> stackhelper;

	// going through all the nodes
	int istack = 0;
	for (int i = 0; i < grid.nxy; ++i) {
		// Getting the row/col indices
		int r, c;
		grid.id2rowcol(i, r, c);
		// if they are base level I include them in the stack
		if (Sreceivers(r, c) == i) {
			stackhelper.emplace(i);
		}

		// While I still have stuff in the stack helper
		while (stackhelper.empty() == false) {
			// I get the next node and pop it from the stack helper
			int nextnode = stackhelper.top();
			stackhelper.pop();

			// registering it to the stack
			Stack[istack] = nextnode;
			++istack;

			// 2D indices
			int rr, cc;
			grid.id2rowcol(nextnode, rr, cc);

			// as well as all its donors which will be processed next
			for (int j = 0; j < Ndonors(rr, cc); ++j) {
				stackhelper.emplace(donors(rr, cc, j));
			}
		}
	}
}

};
