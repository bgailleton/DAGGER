#pragma once

#include "boundary_conditions.hpp"
#include "enumutils.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"
#include "wrap_helper.hpp"

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
	template<class arrin_t>
	void set_boundaries(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_boundaries = std::vector<BC>(arr.size());
		for (int i = 0; i < arr.size(); ++i)
			this->_boundaries[i] = static_cast<BC>(arr[i]);
	}
	template<class out_t>
	out_t get_boundaries()
	{
		std::vector<std::uint8_t> out(this->_boundaries.size(), 0);
		for (int i = 0; i < this->_boundaries.size(); ++i)
			out[i] = static_cast<uint8_t>(this->_boundaries[i]);
		return format_output<decltype(out), out_t>(out);
	}

	// Storing Graph-related stuffies
	std::vector<i_t> _stack;
	std::vector<i_t> _Sstack;

	// Universal data
	// #Topographic surface
	std::vector<f_t> _surface;
	template<class arrin_t>
	void set_surface(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_surface = to_vec(arr);
	}
	template<class out_t>
	out_t get_surface()
	{
		return format_output<decltype(this->_surface), out_t>(this->_surface);
	}

	// #Water surface
	std::vector<f_t> _hw;
	template<class arrin_t>
	void set_hw(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_hw = to_vec(arr);
	}
	template<class out_t>
	out_t get_hw()
	{
		return format_output<decltype(this->_hw), out_t>(this->_hw);
	}

	// #Time tracker for graphflood
	std::vector<f_t> _timetracker;

	// LEM using water data
	std::vector<f_t> _vmot_hw;
	std::vector<f_t> _Qwin;
	template<class arrin_t>
	void set_Qwin(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_Qwin = to_vec(arr);
	}
	template<class out_t>
	out_t get_Qwin()
	{
		return format_output<decltype(this->_Qwin), out_t>(this->_Qwin);
	}

	std::vector<f_t> _Qwout;
	template<class arrin_t>
	void set_Qwout(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_Qwout = to_vec(arr);
	}
	template<class out_t>
	out_t get_Qwout()
	{
		return format_output<decltype(this->_Qwout), out_t>(this->_Qwout);
	}

	std::vector<f_t> _Qsin;
	template<class arrin_t>
	void set_Qsin(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_Qsin = to_vec(arr);
	}
	template<class out_t>
	out_t get_Qsin()
	{
		return format_output<decltype(this->_Qsin), out_t>(this->_Qsin);
	}

	std::vector<f_t> _Qsout;
	template<class arrin_t>
	void set_Qsout(arrin_t& tarr)
	{
		auto arr = format_input<arrin_t>(tarr);
		this->_Qsout = to_vec(arr);
	}
	template<class out_t>
	out_t get_Qsout()
	{
		return format_output<decltype(this->_Qsout), out_t>(this->_Qsout);
	}

	std::vector<f_t> _debug;

	std::vector<std::vector<f_t>> fbag;
	std::vector<std::vector<i_t>> ibag;
	std::vector<std::vector<std::uint8_t>> u8bag;

	std::shared_ptr<easyRand> randu = std::make_shared<easyRand>();

}; // end of class Hermes

} //  end of namespace
