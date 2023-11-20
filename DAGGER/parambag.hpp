#pragma once

#include "boundary_conditions.hpp"
#include "enumutils.hpp"
#include "graphflood_enums.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"
#include "wrap_helper.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class ParamBag
{
public:
	ParamBag(){};

	// GRAPHFLOOD PARAMETRY

	f_t GRAVITY = 9.8;
	f_t rho_sed = 1000;
	f_t tau_c = 4.;
	f_t alpha = 1.5;

	BOUNDARY_HW gf2Bmode = BOUNDARY_HW::FIXED_SLOPE;
	f_t gf2Bbval = 1e-2;

	f_t TSG_dist = false;
	void enable_TSG_dist() { this->TSG_dist = true; }
	void disable_TSG_dist() { this->TSG_dist = false; }
	f_t TSG_distmax = 1e9;
	void set_TSG_distmax(f_t val) { this->TSG_distmax = val; }

	bool gf2_diffuse_Qwin = false;
	void enable_gf2_diffuse_Qwin() { this->gf2_diffuse_Qwin = true; }
	void disable_gf2_diffuse_Qwin() { this->gf2_diffuse_Qwin = false; }

	bool gf2_morpho = false;
	void enable_gf2_morpho() { this->gf2_morpho = true; }
	void disable_gf2_morpho() { this->gf2_morpho = false; }

	f_t ke = 1e-3;
	void set_ke(f_t val) { this->ke = val; }
	f_t get_ke() { return this->ke; }

	f_t kd = 10;
	void set_kd(f_t val) { this->ke = val; }
	f_t get_kd() { return this->kd; }
};

};
