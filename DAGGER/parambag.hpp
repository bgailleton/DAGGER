#pragma once

#include "boundary_conditions.hpp"
#include "enumutils.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"
#include "wrap_helper.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class ParamBag
{
public:
	ParamBag(){};

	f_t ke = 1e-3;
	void set_ke(f_t val) { this->ke = val; }
	f_t get_ke() { return this->ke; }

	f_t kd = 10;
	void set_kd(f_t val) { this->ke = val; }
	f_t get_kd() { return this->kd; }
};

};
