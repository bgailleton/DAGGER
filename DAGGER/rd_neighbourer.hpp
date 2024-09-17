/*
This file is experimental && contains equivalent neighbouring operation from
Riverdale GPU model

*/

#pragma once

#include <array>
#include <functional>
#include <vector>

enum class boundaries : std::uint8_t
{
	normal = 0,
	periodicEW = 1,
	periodicNS = 2,
	customs = 3
};

// template<class i_t>

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_top_row_normal(I_T i,
											I_T j,
											I_T k,
											CONTAINER_BCs_T& BCs,
											PAR* par,
											bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the top
	row Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 30/04/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st row

	if (i == 0) {
		// # Checking all the different cases: fretnei[0]s, last cols && the middle
		if (k == 0)
			valid = false;
	}

	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_leftest_col_normal(I_T i,
													I_T j,
													I_T k,
													CONTAINER_BCs_T& BCs,
													PAR* par,
													bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the
	leftest column Caution: this is optimised for neighbouring checks && ignores
	the top && bottom rows Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 30/04/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st col
	if (j == 0)
		if (k == 1)
			valid = false;
	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_rightest_col_normal(I_T i,
													 I_T j,
													 I_T k,
													 CONTAINER_BCs_T& BCs,
													 PAR* par,
													 bool valid)
{

	/*
	Internal function to check if neighbouring is possible for nodes at the
	rightest column Caution: this is optimised for neighbouring checks && ignores
	the top && bottom rows Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 30/04/2024)
	*/
	// # Only checking if it actually is on the fretnei[0]st col
	if (j == par->nx - 1) {
		if (k == 2) {
			valid = false;
		}
	}
	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_bottom_row_normal(I_T i,
												 I_T j,
												 I_T k,
												 CONTAINER_BCs_T& BCs,
												 PAR* par,
												 bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the bottom
	row Caution: this is optimised for neighbouring checks && ignores the top &&
	bottom rows Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 30/04/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st row
	if (i == par->ny - 1) {
		// # Checking all the different cases: fretnei[0]s, last cols && the middle
		if (k == 3)
			valid = false;
	}
	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
std::array<I_T, 2>
_cast_neighbour_normal(I_T i,
											 I_T j,
											 I_T k,
											 bool valid,
											 CONTAINER_BCs_T& BCs,
											 PAR* par)
{
	/*
	Internal function that cast the neighbours to the right values in the case of
	normal boundary conditions Caution: this is optimised for neighbouring checks
	&& should not be used on its own Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
		- valid: a boolean from previous checks
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 30/04/2024)
	*/

	// # Preformat the output
	std::array<I_T, 2> retnei;
	retnei[0] = -1;
	retnei[1] = -1;

	// # if the neighbouring operation is still valid after that:
	if (valid) {
		if (k == 0) {
			retnei[0] = i - 1;
			retnei[1] = j;
		} else if (k == 1) {
			retnei[0] = i;
			retnei[1] = j - 1;
		} else if (k == 2) {
			retnei[0] = i;
			retnei[1] = j + 1;
		} else if (k == 3) {
			retnei[0] = i + 1;
			retnei[1] = j;
		}
	}

	return retnei;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
std::array<I_T, 2>
_neighbours_normal(I_T i, I_T j, I_T k, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
	GPU function returning the neighbours of a given pixel
	Arguments:\
		- i,j are the row && col indices
		- k is the nth neighbour (4 in D4) following riverdale's convention (see top
	of this module)
		- BCs: boundary conditions code. Note that for this function it does not do
	anything && won't be used but it complies to the standard Returns:
		- (-1,-1) if hte neighbour is not valid (e.g. normal boundaries at the left
	border has no left neighbour k=1)
		- the indices of the row/col of the neighbours
	Authors:
		- B.G. (last modification on the 29th of April)
	*/

	// # I fretnei[0]st assume this mneighbour is valid
	bool valid = true;

	// # Breaking down the checks
	valid = _check_top_row_normal(i, j, k, BCs, par, valid);
	valid = _check_leftest_col_normal(i, j, k, BCs, par, valid);
	valid = _check_rightest_col_normal(i, j, k, BCs, par, valid);
	valid = _check_bottom_row_normal(i, j, k, BCs, par, valid);

	// # getting the actual neighbours
	return _cast_neighbour_normal(i, j, k, valid, BCs, par);
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_is_active_normal(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Quick utility function determining if a node is active || not for normal
		boundaries Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True if the node is active
			- False if inactive (i.e in this case outs)
	*/
	bool valid = true;
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_can_receive_normal(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Standard complying function for the normal boundaries
		Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True, all the nodes can receive in the normal boundary conditions
	*/
	bool valid = true;
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_can_give_normal(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Standard complying function for the normal boundaries
		Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True, all the nodes can receive in the normal boundary conditions
	*/
	bool valid = true;
	if (i == 0 || j == 0 || i == par->ny - 1 || j == par->nx - 1)
		valid = false;
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_can_out_normal(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Standard complying function for the normal boundaries
		Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True, all the nodes can receive in the normal boundary conditions
	*/
	bool valid = false;
	if (i == 0 || j == 0 || i == par->ny - 1 || j == par->nx - 1)
		valid = true;
	return valid;
}

/*

#################################################################################################
############################## Customs Boundaries
###############################################
#################################################################################################

'''
Reminder, I am using the DAGGER convention
// Cannot flow at all = nodata
NO_FLOW = 0,

// Internal Node (can flow in every dretnei[0]ections)
FLOW = 1,

// Internal Node (can flow in every dretnei[0]ections) BUT neighbours a special
flow
// condition && may need specific care
FLOW_BUT = 2,

// flow can out there but can also flow to downstream neighbours
CAN_OUT = 3,

// flow can only out from this cell
OUT = 4,

// Not only flow HAS to out there: neighbouring flows will be drained there no
// matter what
FORCE_OUT = 5,

// Flows through the cell is possible, but the cell CANNOT out fluxes from
// this boundary (reserved to model edges, internal boundaries wont give to
// nodata anyway)
CANNOT_OUT = 6,

// Flow can only flow to potential receivers
IN = 7,

// Forced INFLOW: flow will flow to all neighbours (except other FORCE_IN)
FORCE_IN = 8,

// periodic border
PERIODIC_BORDER = 9
'''
*/

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_top_row_customs(I_T i,
											 I_T j,
											 I_T k,
											 CONTAINER_BCs_T& BCs,
											 PAR* par,
											 bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the top
	row Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 02/05/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st row
	if (i == 0) {
		// # Checking all the different cases: fretnei[0]s, last cols && the middle
		if ((j == 0 && k <= 1) || (j == par->nx - 1 && (k == 0 || k == 2)) ||
				(k == 0))
			valid = false;
	}
	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_leftest_col_customs(I_T i,
													 I_T j,
													 I_T k,
													 CONTAINER_BCs_T& BCs,
													 PAR* par,
													 bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the
	leftest column Caution: this is optimised for neighbouring checks && ignores
	the top && bottom rows Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 02/05/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st col
	if (j == 0)
		if (k == 1)
			valid = false;
	// # Done
	return valid;
}
template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_rightest_col_customs(I_T i,
														I_T j,
														I_T k,
														CONTAINER_BCs_T& BCs,
														PAR* par,
														bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the
	rightest column Caution: this is optimised for neighbouring checks && ignores
	the top && bottom rows Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 02/05/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st col
	if (j == par->nx - 1)
		if (k == 2)
			valid = false;
	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_check_bottom_row_customs(I_T i,
													I_T j,
													I_T k,
													CONTAINER_BCs_T& BCs,
													PAR* par,
													bool valid)
{
	/*
	Internal function to check if neighbouring is possible for nodes at the bottom
	row Caution: this is optimised for neighbouring checks && ignores the top &&
	bottom rows Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 02/05/2024)
	*/

	// # Only checking if it actually is on the fretnei[0]st row
	if (i == par->ny - 1)
		// # Checking all the different cases: fretnei[0]s, last cols && the middle
		if ((j == 0 && (k == 1 || k == 3)) ||
				(j == par->nx - 1 && (k == 3 || k == 2)) || (k == 3))
			valid = false;
	// # Done
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
std::array<I_T, 2>
_cast_neighbour_customs(I_T i,
												I_T j,
												I_T k,
												bool valid,
												CONTAINER_BCs_T& BCs,
												PAR* par)
{
	/*
	Internal function that cast the neighbours to the right values in the case of
	normal boundary conditions Caution: this is optimised for neighbouring checks
	&& should not be used on its own Arguments:
		- i: Row index
		- j: column index
		- k: neighbour number (See top of this module for explanations)
		- valid: a boolean from previous checks
	Returns:
		- a boolean: True = neighbour is valid, False: not a neighbour
	Authors:
		- B.G. (last modification 02/05/2024)
	*/

	// # Preformat the output
	std::array<I_T, 2> retnei;
	retnei[0] = -1;
	retnei[1] = -1;

	// # if the neighbouring operation is still valid after that:
	if (valid) {
		if (k == 0) {
			retnei[0] = i - 1;
			retnei[1] = j;
		}
		if (k == 1) {
			retnei[0] = i;
			retnei[1] = j - 1;
		}
		if (k == 2) {
			retnei[0] = i;
			retnei[1] = j + 1;
		}
		if (k == 3) {
			retnei[0] = i + 1;
			retnei[1] = j;
		}
	}
	if (BCs(i, j) == 0 || (retnei[0] != -1 && BCs(retnei[0], retnei[1]) == 0)) {
		retnei[0] = -1;
		retnei[1] = -1;
	}

	return retnei;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
std::array<I_T, 2>
_neighbours_customs(I_T i, I_T j, I_T k, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
	GPU function returning the neighbours of a given pixel
	Arguments:\
		- i,j are the row && col indices
		- k is the nth neighbour (4 in D4) following riverdale's convention (see top
	of this module)
		- BCs: boundary conditions code. Note that for this function it does not do
	anything && won't be used but it complies to the standard Returns:
		- (-1,-1) if hte neighbour is not valid (e.g. normal boundaries at the left
	border has no left neighbour k=1)
		- the indices of the row/col of the neighbours
	Authors:
		- B.G. (last modification 02/05/2024)
	TODO:
		- adding the Periodic boundary management in the checks
	*/
	bool valid = true;

	// # Breaking down the checks
	valid = _check_top_row_customs(i, j, k, BCs, par, valid);
	valid = _check_leftest_col_customs(i, j, k, BCs, par, valid);
	valid = _check_rightest_col_customs(i, j, k, BCs, par, valid);
	valid = _check_bottom_row_customs(i, j, k, BCs, par, valid);

	// # getting the actual neighbours
	return _cast_neighbour_customs(i, j, k, valid, BCs, par);
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_can_receive_customs(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Standard complying function for the normal boundaries
		Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True, all the nodes can receive in the normal boundary conditions
	*/
	bool valid = true;
	if (BCs(i, j) == 6 || BCs(i, j) == 7 || BCs(i, j) == 8 || BCs(i, j) == 0)
		valid = false;
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_can_give_customs(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Standard complying function for the normal boundaries
		Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True, all the nodes can receive in the normal boundary conditions
		Authors:
		- B.G. (last modification 02/05/2024)
	*/
	bool valid = false;
	if (BCs(i, j) == 1 || BCs(i, j) == 6 || BCs(i, j) == 7 || BCs(i, j) == 8 ||
			BCs(i, j) == 9)
		valid = true;
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_can_out_customs(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Standard complying function for the normal boundaries
		Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True, all the nodes can receive in the normal boundary conditions
		Authors:
		- B.G. (last modification 02/05/2024)
	*/
	bool valid = false;
	if (BCs(i, j) == 3 || BCs(i, j) == 4 || BCs(i, j) == 5)
		valid = true;
	return valid;
}

template<class I_T, class CONTAINER_BCs_T, class PAR>
bool
_is_active_customs(I_T i, I_T j, CONTAINER_BCs_T& BCs, PAR* par)
{
	/*
		Quick utility function determining if a node is active || not for normal
		boundaries Arguments:
			- i: the row index
			- j: the column index
			- BCs: a dummy field to keep the standard consistent
		Returns:
			- True if the node is active
			- False if inactive (i.e in this case outs)
		Authors:
		- B.G. (last modification 02/05/2024)
	*/
	bool valid = true;

	if (BCs(i, j) == 0)
		valid = false;
	return valid;
}

template<class I_T, class F_T, class CONTAINER_BCs_T>
class GridCPP
{

public:
	GridCPP() { ; }
	GridCPP(I_T nx, I_T ny, F_T dx, F_T dy, std::uint8_t bouco)
	{
		this->nx = nx;
		this->ny = ny;
		this->dx = dx;
		this->dy = dy;
		this->nxy = nx * ny;
		this->bound = static_cast<boundaries>(bouco);
		this->assign_function();
	}

	I_T nx = 0;
	I_T ny = 0;
	F_T dx = 0.;
	F_T dy = 0.;
	I_T nxy = 0;
	I_T nk = 4;

	boundaries bound = boundaries::normal;

	// Standard API
	std::function<std::array<
		I_T,
		2>(I_T, I_T, I_T, CONTAINER_BCs_T&, GridCPP<I_T, F_T, CONTAINER_BCs_T>*)>
		neighbours;
	std::function<
		bool(I_T, I_T, CONTAINER_BCs_T&, GridCPP<I_T, F_T, CONTAINER_BCs_T>*)>
		can_receive;
	std::function<
		bool(I_T, I_T, CONTAINER_BCs_T&, GridCPP<I_T, F_T, CONTAINER_BCs_T>*)>
		can_give;
	std::function<
		bool(I_T, I_T, CONTAINER_BCs_T&, GridCPP<I_T, F_T, CONTAINER_BCs_T>*)>
		can_out;
	std::function<
		bool(I_T, I_T, CONTAINER_BCs_T&, GridCPP<I_T, F_T, CONTAINER_BCs_T>*)>
		is_active;

	void assign_function()
	{
		if (this->bound == boundaries::normal) {
			this->neighbours = _neighbours_normal<I_T,
																						CONTAINER_BCs_T,
																						GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->can_receive =
				_can_receive_normal<I_T,
														CONTAINER_BCs_T,
														GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->can_give = _can_give_normal<I_T,
																				CONTAINER_BCs_T,
																				GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->can_out = _can_out_normal<I_T,
																			CONTAINER_BCs_T,
																			GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->is_active = _is_active_normal<I_T,
																					CONTAINER_BCs_T,
																					GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
		} else if (this->bound == boundaries::customs) {
			this->neighbours =
				_neighbours_customs<I_T,
														CONTAINER_BCs_T,
														GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->can_receive =
				_can_receive_customs<I_T,
														 CONTAINER_BCs_T,
														 GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->can_give = _can_give_customs<I_T,
																				 CONTAINER_BCs_T,
																				 GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->can_out = _can_out_customs<I_T,
																			 CONTAINER_BCs_T,
																			 GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
			this->is_active = _is_active_customs<I_T,
																					 CONTAINER_BCs_T,
																					 GridCPP<I_T, F_T, CONTAINER_BCs_T>>;
		}
	}

	template<class T>
	void id2rowcol(T index, T& row, T& col)
	{
		col = index % this->nx;
		row = std::floor(double(index) / this->nx);
	}

	template<class T>
	T rowcol2index(T row, T col)
	{
		return row * this->nx + col;
	}
};
