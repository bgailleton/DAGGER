#pragma once

#include "boundary_conditions.hpp"
#include "enumutils.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class Hermes
{
public:
	Hermes(){};

	// Storing lookup tables
	lookup8<i_t, f_t> LK8;

	// Storing Connector-related stuffies
	std::vector<uint8_t> _neighbours;
	std::vector<uint8_t> _Sreceivers;
	std::vector<uint8_t> _Sdonors;
	std::vector<uint8_t> _receivers;
	std::vector<uint8_t> _donors;
	std::vector<BC> _boundaries;

	// Storing Graph-related stuffies
	std::vector<i_t> _stack;
	std::vector<i_t> _Sstack;

	// Universal data
	std::vector<f_t> _surface;

	std::vector<std::vector<f_t>> fbag;
	std::vector<std::vector<i_t>> ibag;
	std::vector<std::vector<std::uint8_t>> u8bag;

	std::shared_ptr<easyRand> randu = std::make_shared<easyRand>();

}; // end of class Hermes

} //  end of namespace
