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
	// Empty constructor
	ParamBag(){};

	// GRAPHFLOOD PARAMETERS

	// CONSTANTS
	f_t GRAVITY = 9.8;
	f_t rho_sed = 1000;

	// Stuff I'll vary later
	f_t tau_c = 4.;
	f_t alpha = 1.5;

	BOUNDARY_HW gf2Bmode = BOUNDARY_HW::FIXED_SLOPE;
	f_t gf2Bbval = 1e-2;
	void set_gf2Bbval(f_t val) { this->gf2Bbval = val; }
	f_t get_gf2Bbval() { return this->gf2Bbval; }

	bool gf2_morpho = false;
	void enable_gf2_morpho() { this->gf2_morpho = true; }
	void disable_gf2_morpho() { this->gf2_morpho = false; }

	f_t time_dilatation_morpho = 1.;
	void set_time_dilatation_morpho(f_t val)
	{
		this->time_dilatation_morpho = val;
	}
	f_t get_time_dilatation_morpho() { return this->time_dilatation_morpho; }

	f_t _ke = 1e-3;
	std::vector<f_t> _ke_v;
	PARAMTYPE p_ke = PARAMTYPE::CTE;
	void set_ke(f_t val)
	{
		this->p_ke = PARAMTYPE::CTE;
		this->_ke = val;
	}
	f_t get_ke() { return this->_ke; }
	f_t ke(i_t node)
	{
		return (this->p_ke == PARAMTYPE::CTE) ? this->_ke : this->_ke_v[node];
	}

	template<class arrin_t>
	void set_variable_ke(arrin_t& tarr)
	{
		this->p_ke = PARAMTYPE::VAR;
		auto arr = format_input<arrin_t>(tarr);
		this->_ke_v = to_vec(arr);
	}

	f_t kd = 10;
	void set_kd(f_t val) { this->kd = val; }
	f_t get_kd() { return this->kd; }

	bool bank_erosion = false;
	void enable_bank_erosion() { this->bank_erosion = true; }
	void disable_bank_erosion() { this->bank_erosion = false; }
	f_t kel = 0.1;
	void set_kel(f_t val) { this->kel = val; }
	f_t get_kel() { return this->kel; }

	f_t kdl = 0.1;
	void set_kdl(f_t val) { this->kdl = val; }
	f_t get_kdl() { return this->kdl; }

	// # Graphflood experimental stuff

	f_t TSG_dist = false;
	void enable_TSG_dist() { this->TSG_dist = true; }
	void disable_TSG_dist() { this->TSG_dist = false; }
	f_t TSG_distmax = 1e9;
	void set_TSG_distmax(f_t val) { this->TSG_distmax = val; }

	bool gf2_diffuse_Qwin = false;
	void enable_gf2_diffuse_Qwin() { this->gf2_diffuse_Qwin = true; }
	void disable_gf2_diffuse_Qwin() { this->gf2_diffuse_Qwin = false; }

	bool transient_flow = false;
	void enable_transient_flow() { this->transient_flow = true; }
	void disable_transient_flow() { this->transient_flow = false; }
};

};
