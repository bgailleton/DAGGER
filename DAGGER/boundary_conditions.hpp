//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

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
#include <stdlib.h>
#include <ctime>




// // Difines the different boundary conditions for a single node

// class BoundaryCondition
// {

// public:
// 	BoundaryCondition();

// 	// Can flow leave the model (if true, flow paths can stop there)
// 	template<class i_t>
// 	virtual bool can_out(i_t i) = 0;
// 	// if true, the flow is stopped there
// 	template<class i_t>
// 	virtual bool force_out(i_t i) = 0;
// 	// if true node is nodata and ignored
// 	template<class i_t>
// 	virtual bool nodata(i_t i) = 0;
// 	// if true, the node is and will always be receiver of its neighbours
// 	template<class i_t>
// 	virtual bool forced_receiver(i_t i) = 0;
// 	// if truem the node is and will always be a donor
// 	template<class i_t>
// 	virtual bool forced_donor(i_t i) = 0;


// };

namespace DAGGER
{


enum class BC: std::uint8_t
{
	// Cannot flow at all = nodata
	NO_FLOW,

	// Internal Node (can flow in every directions)
	FLOW,

	// Internal Node (can flow in every directions) BUT neighbours a special flow condition and may need specific care
	FLOW_BUT,

	// flow can out there but can also flow to downstream neighbours
	CAN_OUT,

	// flow can only out from this cell
	OUT,

	// Not only flow HAS to out there: neighbouring flows will be drained there no matter what
	FORCE_OUT,

	// Flow can only flow to potential receivers
	IN,

	// Forced INFLOW: flow will flow to all neighbours (except other FORCE_IN)
	FORCE_IN,

	// periodic border
	PERIODIC_BORDER,

};

std::string BC2str(BC tbc)
{
	if(tbc == NO_FLOW) return "NO_FLOW";
	if(tbc == FLOW) return "FLOW";
	if(tbc == FLOW_BUT) return "FLOW_BUT";
	if(tbc == CAN_OUT) return "CAN_OUT";
	if(tbc == OUT) return "OUT";
	if(tbc == FORCE_OUT) return "FORCE_OUT";
	if(tbc == FORCE_IN) return "FORCE_IN";
	if(tbc == PERIODIC_BORDER) return "PERIODIC_BORDER";

	return "UNREGISTERED BC";
}


// Function testing if a boundary condition object is valid
// template<class TBC>
// void test_TBC(TBC)





// Class complying with the boundary conditions standards
class CodeBC
{
public:

	CodeBC(){;}

	std::vector<BC> boundaries;

	// return true if the node is an internal "classic" node and no specific conditions are required
	// this will be the case of the huge majority of nodes in most cases, that is why it is important
	// to separate and prioritise their detection for the sake of efficiency
	template<class i_t>
	bool is_normal_node(i_t i){BC tbc = this->boundary[i]; if(tbc == BC::FLOW) return true; else return false;}

	// return true if the node is an internal "classic" node but borders a boundary condition requiring specific care]
	// this also is an import test to optimise and prioritise: a receivers of a neigh
	template<class i_t>
	bool is_indirect_bc(i_t i){BC tbc = this->boundary[i]; if(tbc == BC::FLOW_BUT) return true; else return false;}

	template<class i_t>
	bool is_bc(i_t i){BC tbc = this->boundary[i]; if(tbc == BC::FLOW_BUT || tbc == BC::FLOW || tbc == BC::PIT) return false; else return true;}


	template<class i_t>
	bool can_receive(i_t i)
	{
		BC tbc = this->boundary[i];
		if(tbc == BC::NO_FLOW || tbc == BC::FORCE_IN || tbc == BC::IN)
			return false;
		return true;
	}

	template<class i_t>
	bool can_give(i_t i)
	{
		BC tbc = this->boundary[i];
		if(tbc == BC::NO_FLOW || tbc == BC::FORCE_OUT || tbc == BC::OUT)
			return false;
		return true;
	}

	template<class i_t>
	bool can_flow_through(i_t i)
	{
		if(this->can_give(i) && this->can_receive(i))
			return true;
		return false;
	}

};











};


























































#endif