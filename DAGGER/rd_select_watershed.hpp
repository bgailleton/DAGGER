/*
This file is experimental && contains equivalent neighbouring operation from
Riverdale GPU model

*/

#pragma once

#include <array>
#include <functional>
#include <stack>
#include <tuple>
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

/*
Create a mask of No Data
*/
template<class GRID_T, class F_T>
void
BCs_to_mask(GRID_T& grid,
						xt::pytensor<std::uint8_t, 2>& mask,
						xt::pytensor<std::uint8_t, 2>& BCs)
{

	for (int r = 0; r < grid.ny; ++r) {
		for (int c = 0; c < grid.nx; ++c) {
			mask(r, c) = grid.is_active(r, c, BCs, &grid);
		}
	}
}

/*
Take a mask, adapts the boundary condition array to redraw the boundaries
Warning, invalidates the graph related to the old boundary conditions
*/
template<class GRID_T, class F_T>
void
mask_to_BCs(GRID_T& grid,
						xt::pytensor<std::uint8_t, 2>& mask,
						xt::pytensor<std::uint8_t, 2>& BCs,
						bool force_out)
{

	// What is the code for BCs out
	std::uint8_t tbc = (force_out) ? 4 : 3;

	// Second pass: design_outlet
	for (int r = 0; r < grid.ny; ++r) {
		for (int c = 0; c < grid.nx; ++c) {

			//
			if (grid.is_active(r, c, BCs, &grid) && mask(r, c)) {

				//
				for (int k = 0; k < grid.nk; ++k) {
					// getting the neighbour
					auto nei = grid.neighbours(r, c, k, BCs, &grid);

					// chekcing the neighbour
					if (nei[0] == -1)
						continue;

					if (mask(nei[0], nei[1]) == false)
						BCs(r, c) = tbc;
				}
			}
		}
	}

	// First pass: 0 where no mask or original BCs false
	for (int r = 0; r < grid.ny; ++r) {
		for (int c = 0; c < grid.nx; ++c) {
			if (mask(r, c) == 0)
				BCs(r, c) = 0;
		}
	}
}

/*
correct_
*/
template<class GRID_T, class F_T>
void
correct_Sreceivers_from_mask(GRID_T& grid,
														 xt::pytensor<int, 2>& Sreceivers,
														 xt::pytensor<std::uint8_t, 2>& mask,
														 xt::pytensor<std::uint8_t, 2>& BCs)
{

	for (int node = 0; node < grid.nxy; ++node) {

		int r, c;
		grid.id2rowcol(node, r, c);
		int rr, cc;
		grid.id2rowcol(Sreceivers(r, c), rr, cc);
		if (grid.is_active(r, c, BCs, &grid) == false || mask(r, c) == false)
			Sreceivers(r, c) = node;
		else {
			if (mask(r, c) && mask(rr, cc) == false) {
				Sreceivers(r, c) = node;
			}
		}
	}
}

/*
Takes a mask and a graph and label eah individual basins
*/
template<class GRID_T, class F_T>
void
label_watersheds_mask(GRID_T& grid,
											xt::pytensor<int, 2>& basin_label,
											xt::pytensor<int, 1>& Stack,
											xt::pytensor<int, 2>& Sreceivers,
											xt::pytensor<std::uint8_t, 2>& mask,
											xt::pytensor<std::uint8_t, 2>& BCs)
{

	int lab = 1;

	for (auto node : Stack) {

		int r, c;
		grid.id2rowcol(node, r, c);
		int rr, cc;
		grid.id2rowcol(Sreceivers(r, c), rr, cc);

		if (mask(r, c) == false)
			continue;

		if (grid.is_active(r, c, BCs, &grid) == false)
			continue;

		if (Sreceivers(r, c) == node || mask(rr, cc) == false) {
			++lab;
			basin_label(r, c) = lab;
		} else {
			basin_label(r, c) = basin_label(rr, cc);
		}
	}
}

template<class GRID_T, class F_T>
void
mask_watersheds_min_area(GRID_T& grid,
												 xt::pytensor<int, 1>& Stack,
												 xt::pytensor<int, 2>& Sreceivers,
												 xt::pytensor<std::uint8_t, 2>& mask,
												 xt::pytensor<std::uint8_t, 2>& BCs,
												 F_T area_min)
{

	std::vector<F_T> area(grid.nxy, 0.);

	// first calculating drainage area
	for (auto it = Stack.rbegin(); it != Stack.rend(); ++it) {

		int node = *it;
		int r, c;
		grid.id2rowcol(node, r, c);

		if (grid.is_active(r, c, BCs, &grid) == false)
			continue;

		mask(r, c) = false;
		int rec = Sreceivers(r, c);
		area[node] += grid.dx * grid.dy;

		if (rec != node) {
			area[rec] += area[node];
		}
	}

	for (auto node : Stack) {
		int r, c;
		grid.id2rowcol(node, r, c);

		if (grid.is_active(r, c, BCs, &grid) == false)
			continue;

		int rr, cc;
		grid.id2rowcol(Sreceivers[node], rr, cc);
		if (area[node] >= area_min)
			mask(r, c) = true;
		else
			mask(r, c) = mask(rr, cc);
	}
}

template<class GRID_T, class F_T>
void
mask_watersheds_above_elevations(GRID_T& grid,
																 xt::pytensor<F_T, 2>& topography,
																 xt::pytensor<int, 1>& Stack,
																 xt::pytensor<int, 2>& Sreceivers,
																 xt::pytensor<std::uint8_t, 2>& mask,
																 xt::pytensor<std::uint8_t, 2>& BCs,
																 F_T elevation_min,
																 bool exclude_if_does_not_reach)
{

	for (auto node : Stack) {
		int r, c;
		grid.id2rowcol(node, r, c);
		mask(r, c) == false;

		if (grid.is_active(r, c, BCs, &grid) == false) {
			continue;
		}

		if (topography(r, c) <= elevation_min) {
			continue;
		}

		int rr, cc;
		grid.id2rowcol(Sreceivers(r, c), rr, cc);

		if (Sreceivers(r, c) == node) {
			if (exclude_if_does_not_reach == false) {
				mask(r, c) = true;
			}

		} else if (topography(rr, cc) <= elevation_min) {
			mask(r, c) = true;
		} else {
			mask(r, c) = mask(rr, cc);
		}
	}
}

template<class GRID_T, class F_T>
std::tuple<int, int, int, int, int, int>
bounding_box_from_label(GRID_T& grid,
												xt::pytensor<std::int32_t, 2>& labels,
												std::int32_t lab,
												xt::pytensor<std::int32_t, 2>& Sreceivers)
{

	int row_max = 0;
	int row_min = grid.ny - 1;
	int col_max = 0;
	int col_min = grid.nx - 1;
	int row_out = 0;
	int col_out = 0;

	for (int i = 0; i < grid.ny; ++i) {
		for (int j = 0; j < grid.nx; ++j) {

			if (labels(i, j) == lab) {
				int node = grid.rowcol2index(i, j);

				row_max = std::max(i, row_max);
				row_min = std::min(i, row_min);
				col_max = std::max(j, col_max);
				col_min = std::min(j, col_min);

				if (Sreceivers(i, j) == node) {
					row_out = i;
					col_out = j;
				}
			}
		}
	}

	return std::make_tuple(row_min, row_max, col_min, col_max, row_out, col_out);
}

template<class GRID_T, class F_T>
void
mask_upstream_MFD(GRID_T& grid,
									xt::pytensor<std::uint8_t, 2>& mask,
									xt::pytensor<std::int32_t, 2>& topography,
									xt::pytensor<std::uint8_t, 2>& BCs,
									int row,
									int col)
{

	mask.fill(0);

	std::queue<int> mQ;
	mask(row, col) = 1;

	mQ.emplace(grid.rowcol2index(row, col));

	while (mQ.empty() == false) {

		int next = mQ.front();

		mQ.pop();

		int r, c;
		grid.id2rowcol(next, r, c);

		if (grid.is_active(r, c, BCs, &grid)) {

			//
			for (int k = 0; k < grid.nk; ++k) {
				// getting the neighbour
				auto nei = grid.neighbours(r, c, k, BCs, &grid);

				// chekcing the neighbour
				if (nei[0] == -1)
					continue;

				if (mask(nei[0], nei[1]) == false &&
						topography(r, c) <= topography(nei[0], nei[1])) {
					mask(nei[0], nei[1]) = 1;
					mQ.emplace(grid.rowcol2index(nei[0], nei[1]));
				}
			}
		}
	}
}

};
