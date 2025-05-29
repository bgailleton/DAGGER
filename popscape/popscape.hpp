#ifndef POPSCAPE_HPP
#define POPSCAPE_HPP

// STL imports
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <stack>
#include <stdlib.h>
#include <string>
#include <thread>
#include <vector>

// local includes
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D4connector.hpp"
#include "D8connector.hpp"
#include "popscape_utils.hpp"
// #include "D8connector.hpp"
#include "graph.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// #include <pybind11/numpy.h>

// namespace py = pybind11;

namespace DAGGER {

// Enumeration to define parameter modes for the landscape evolution model
// Parameters can either be constant across the entire domain or variable
// (spatially distributed)
enum class PARAM_MODE_POP
{
	CONSTANT, // Single value applied everywhere
	VARIABLE, // Array of values, one per grid node
};

/**
 * Template class for landscape evolution modeling using the "popscape" approach
 * This appears to implement a stream power law erosion model with various
 * parameters
 *
 * Template parameters:
 * - fT: Floating point type (float, double, etc.)
 * - Graph_t: Graph/network type for flow routing
 * - Connector_t: Grid connectivity/topology handler
 */
template<class fT, class Graph_t, class Connector_t>
class popscape
{

public:
	// Core components of the landscape model
	Graph_t graph;				 // Handles flow routing and network topology
	Connector_t connector; // Manages grid connectivity and spatial relationships

	// Primary state variables
	std::vector<fT> QA;					// Discharge/drainage area at each node
	std::vector<fT> topography; // Elevation values at each grid node

	// Erosion rate coefficient (K) parameters
	// Kbase is the base erodibility coefficient, Kmod is a modifier
	PARAM_MODE_POP Kbase_mode = PARAM_MODE_POP::CONSTANT;
	std::vector<fT> _Kbase = { 1e-3 }; // Default: 0.001
	PARAM_MODE_POP Kmod_mode = PARAM_MODE_POP::CONSTANT;
	std::vector<fT> _Kmod = { 1. }; // Default: 1.0 (no modification)

	// Stream power law exponents
	// m: drainage area exponent, n: slope exponent
	PARAM_MODE_POP m_mode = PARAM_MODE_POP::CONSTANT;
	std::vector<fT> _m = { 0.45 }; // Typical value for area exponent
	PARAM_MODE_POP n_mode = PARAM_MODE_POP::CONSTANT;
	std::vector<fT> _n = { 1. }; // Linear slope dependency

	// Precipitation parameters (affects discharge calculation)
	PARAM_MODE_POP precip_mode = PARAM_MODE_POP::CONSTANT;
	std::vector<fT> _precip = { 1. }; // Uniform precipitation

	// Uplift/Erosion equilibrium parameter
	// UE represents the balance between tectonic uplift and erosion
	PARAM_MODE_POP UE_mode = PARAM_MODE_POP::CONSTANT;
	std::vector<fT> _UE = { 1e-3 }; // Default uplift rate

	// Boundary condition specification (as string for flexibility)
	std::string boundary_string =
		"periodic_EW"; // Periodic in East-West direction

	// Grid resolution tracking for multi-resolution operations
	int OG_nx;				 // Original number of columns
	int OG_ny;				 // Original number of rows
	fT fac_resize = 1; // Resize factor for multi-grid operations

	// Default constructor
	popscape() { ; }

	/**
	 * Constructor that initializes the popscape with existing graph and connector
	 * Creates new connector and graph objects based on input dimensions
	 */
	popscape(Graph_t& graph, Connector_t& con)
	{
		// Create a new connector with same dimensions and spacing as input
		this->connector =
			_create_connector(con.nx, con.ny, con.dx, con.dy, fT(0.), fT(0.));
		this->connector.set_default_boundaries(con.boundary_string);

		// Initialize graph with the number of nodes and new connector
		_create_graph(graph.nnodes, this->connector, this->graph);

		this->boundary_string = this->connector.boundary_string;

		// Store original grid dimensions for reference
		this->OG_nx = con.nx;
		this->OG_ny = con.ny;
	}

	/**
	 * Set topography from input data of any compatible type
	 * Uses DAGGER library formatting utilities for type conversion
	 */
	template<class in_t>
	void set_topo(in_t& itopo)
	{
		auto ftopo = DAGGER::format_input(itopo); // Convert to internal format
		this->topography = to_vec(ftopo);					// Convert to vector
	}

	/**
	 * Get topography in specified output format
	 * Template allows flexible output type conversion
	 */
	template<class out_t>
	out_t get_topo()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->topography);
	}

	/**
	 * Get discharge/drainage area in specified output format
	 */
	template<class out_t>
	out_t get_QA()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->QA);
	}

	/**
	 * Calculate and return chi-star values in specified output format
	 * Chi is a transformed coordinate used in landscape evolution analysis
	 */
	template<class out_t>
	out_t get_chistar()
	{
		std::vector<fT> chistar = this->_chi_star();
		return DAGGER::format_output<std::vector<fT>, out_t>(chistar);
	}

	/**
	 * Convert current node index to original grid index
	 * Used when working with resized grids to map back to original coordinates
	 */
	int rindex(int i)
	{
		int row, col;
		// Get row and column from current node ID
		this->connector.rowcol_from_node_id(i, row, col);
		// Calculate corresponding index in original grid
		int nid = floor(row * fac_resize) * this->OG_nx + floor(col * fac_resize);
		return nid;
	}

	// Parameter accessor methods - handle both constant and variable modes

	/**
	 * Get base erodibility coefficient for node i
	 * Returns single value if constant mode, or node-specific value if variable
	 */
	fT Kbase(int i)
	{
		if (this->Kbase_mode == PARAM_MODE_POP::CONSTANT)
			return this->_Kbase[0];
		else
			return this->_Kbase[rindex(i)];
	}

	/**
	 * Get erodibility modifier for node i
	 */
	fT Kmod(int i)
	{
		if (this->Kmod_mode == PARAM_MODE_POP::CONSTANT)
			return this->_Kmod[0];
		else
			return this->_Kmod[rindex(i)];
	}

	/**
	 * Get drainage area exponent for node i
	 */
	fT m(int i)
	{
		if (this->m_mode == PARAM_MODE_POP::CONSTANT)
			return this->_m[0];
		else
			return this->_m[rindex(i)];
	}

	/**
	 * Get slope exponent for node i
	 */
	fT n(int i)
	{
		if (this->n_mode == PARAM_MODE_POP::CONSTANT)
			return this->_n[0];
		else
			return this->_n[rindex(i)];
	}

	// Parameter setter methods for constant values

	/**
	 * Set drainage area exponent to constant value
	 */
	void set_m(fT val) { this->_m[0] = val; }

	/**
	 * Set slope exponent to constant value
	 */
	void set_n(fT val) { this->_n[0] = val; }

	/**
	 * Set base erodibility to constant value
	 */
	void set_Kbase(fT val)
	{
		this->Kbase_mode = PARAM_MODE_POP::CONSTANT;
		this->_Kbase[0] = val;
	}

	/**
	 * Set erodibility modifier to constant value
	 */
	void set_Kmod(fT val)
	{
		this->Kmod_mode = PARAM_MODE_POP::CONSTANT;
		this->_Kmod[0] = val;
	}

	/**
	 * Set spatially variable erodibility modifier
	 * Converts input to internal vector format
	 */
	template<class in_t>
	void set_Kmod_variable(in_t& ival)
	{
		this->Kmod_mode = PARAM_MODE_POP::VARIABLE;
		auto fval = DAGGER::format_input(ival);
		this->_Kmod = to_vec(fval);
	}

	/**
	 * Get precipitation rate for node i
	 */
	fT precip(int i)
	{
		if (this->precip_mode == PARAM_MODE_POP::CONSTANT)
			return this->_precip[0];
		else
			return this->_precip[rindex(i)];
	}

	/**
	 * Set precipitation to constant value
	 */
	void set_precip(fT val)
	{
		this->precip_mode = PARAM_MODE_POP::CONSTANT;
		this->_precip[0] = val;
	}

	/**
	 * Set spatially variable precipitation
	 */
	template<class in_t>
	void set_precip_variable(in_t& ival)
	{
		this->precip_mode = PARAM_MODE_POP::VARIABLE;
		auto fval = DAGGER::format_input(ival);
		this->_precip = to_vec(fval);
	}

	/**
	 * Get uplift/erosion rate for node i
	 * UE stands for Uplift Erosion - represents the equilibrium rate
	 * between tectonic uplift and erosional lowering
	 */
	fT UE(int i)
	{
		if (this->UE_mode == PARAM_MODE_POP::CONSTANT)
			return this->_UE[0];
		else
			return this->_UE[rindex(i)];
	}

	/**
	 * Set uplift/erosion rate to constant value
	 */
	void set_UE(fT val)
	{
		this->UE_mode = PARAM_MODE_POP::CONSTANT;
		this->_UE[0] = val;
	}

	/**
	 * Set spatially variable uplift/erosion rates
	 */
	template<class in_t>
	void set_UE_variable(in_t& ival)
	{
		this->UE_mode = PARAM_MODE_POP::VARIABLE;
		auto fval = DAGGER::format_input(ival);
		this->_UE = to_vec(fval);
	}

	/**
	 * Initialize discharge/drainage area vector with zeros
	 * Called before each computation to reset state
	 */
	void _init_vecs() { this->QA = std::vector<fT>(this->connector.nnodes, 0.); }

	/**
	 * Run steady-state calculation for specified number of iterations
	 * This is the main computational engine of the landscape evolution model
	 */
	void StSt(int n_iterations)
	{
		// Main iteration loop for landscape evolution
		for (int nit = 0; nit < n_iterations; ++nit) {
			// Set depression filling algorithm
			this->graph.depression_resolver = DAGGER::DEPRES::priority_flood;

			// Compute flow routing graph from current topography
			// Parameters: topography, compute_Sreceivers=true, compute_donors=false
			this->graph._compute_graph(this->topography, true, false);

			// Reset discharge/drainage area vector
			this->_init_vecs();

			// Calculate drainage area based on precipitation mode
			if (this->precip_mode == PARAM_MODE_POP::CONSTANT) {
				// Uniform precipitation: accumulate using node area
				this->QA = this->graph._accumulate_constant_downstream_SFD(
					this->connector.get_area_at_node(0));
			} else {
				// Variable precipitation: use spatially distributed rates
				this->QA =
					this->graph._accumulate_variable_downstream_area_SFD(this->_precip);
			}

			// Update elevation for each node using stream power law
			// Process nodes in stack order (from downstream to upstream)
			for (int i = 0; i < this->graph.nnodes; ++i) {
				int node = this->graph.Sstack[i]; // Current node in processing order

				// Skip outlet nodes and pits
				if (this->connector.flow_out_or_pit(node))
					continue;

				// Get downstream receiver node
				int rec = this->connector._Sreceivers[node];

				// Apply stream power law equation:
				// E = K * A^m * S^n, where S = (z_up - z_down) / distance
				// Rearranged to solve for upstream elevation:
				// z_up = z_down + distance * (U/K)^(1/n) * A^(-m/n)
				this->topography[node] =
					this->topography[rec] +											// Downstream elevation
					this->connector.Sdistance2receivers[node] * // Distance to receiver
						std::pow(this->UE(node) / (this->Kbase(node) * this->Kmod(node)),
										 1. / this->n(node)) * // (Uplift/Erodibility)^(1/n)
						std::pow(this->QA[node],
										 -this->m(node) / this->n(node)); // Area^(-m/n)
			}
		}
	}

	/**
	 * Helper function for grid restriction (refinement)
	 * Copies values from coarse grid to fine grid with optional noise
	 * Each coarse cell becomes 4 fine cells (2x2 subdivision)
	 */
	void _restr_helper(std::vector<fT>& oldarr,
										 std::vector<fT>& newarr,
										 fT noise_intensity)
	{

		int nxy = 4 * this->graph.nnodes; // New grid has 4x more nodes

		int nx = this->connector.nx * 2; // Double resolution in x
		int ny = this->connector.ny * 2; // Double resolution in y

		// For each node in the current (coarse) grid
		for (int i = 0; i < this->graph.nnodes; ++i) {
			int row, col;
			this->connector.rowcol_from_node_id(i, row, col);

			// Map to 4 fine grid cells with added noise
			// Bottom-left cell
			int nid = row * 2 * nx + col * 2;
			newarr[nid] = oldarr[i] + this->connector.randu->get() * noise_intensity;

			// Bottom-right cell
			nid = row * 2 * nx + col * 2 + 1;
			newarr[nid] = oldarr[i] + this->connector.randu->get() * noise_intensity;

			// Top-right cell
			nid = (row * 2 + 1) * nx + col * 2 + 1;
			newarr[nid] = oldarr[i] + this->connector.randu->get() * noise_intensity;

			// Top-left cell
			nid = (row * 2 + 1) * nx + col * 2;
			newarr[nid] = oldarr[i] + this->connector.randu->get() * noise_intensity;
		}
	}

	/**
	 * Perform grid restriction (mesh refinement)
	 * Doubles the resolution in both x and y directions
	 * Adds random noise to break symmetry
	 */
	void restriction(fT noise_intensity)
	{
		int nxy = 4 * this->graph.nnodes; // New total number of nodes

		int nx = this->connector.nx * 2;
		int ny = this->connector.ny * 2;

		// Create refined topography grid
		std::vector<fT> ntopo(nxy, 0.);
		this->_restr_helper(this->topography, ntopo, noise_intensity);

		// Commented code would handle variable parameter refinement
		// Currently only topography is refined

		// Calculate new grid spacing (half the original)
		fT dx = this->connector.dx / 2;
		fT dy = this->connector.dy / 2;

		// Create new connector for refined grid
		this->connector = _create_connector(nx, ny, dx, dy, fT(0.), fT(0.));
		this->connector.set_default_boundaries(this->boundary_string);

		// Create new graph for refined grid
		_create_graph(nxy, this->connector, this->graph);

		// Update state with refined topography
		this->topography = std::move(ntopo);
		this->_init_vecs();
		this->border20();			 // Apply boundary conditions
		this->fac_resize /= 2; // Update resize factor
	}

	/**
	 * Helper function for grid interpolation (coarsening)
	 * Averages 4 fine cells to create 1 coarse cell
	 */
	void _interp_helper(std::vector<fT>& oldarr,
											std::vector<fT>& newarr,
											D8connector<fT>& ncon)
	{

		int nx = floor(this->connector.nx / 2); // Half resolution in x
		int ny = floor(this->connector.ny / 2); // Half resolution in y
		int nxy = nx * ny;

		// For each node in the new (coarse) grid
		for (int i = 0; i < nxy; ++i) {
			int row, col;
			ncon.rowcol_from_node_id(i, row, col);

			// Average the 4 corresponding fine grid cells
			int nid = row * 2 * this->connector.nx + col * 2;
			fT tntopo = 0.;
			tntopo += oldarr[nid]; // Bottom-left

			nid = row * 2 * this->connector.nx + col * 2 + 1;
			tntopo += oldarr[nid]; // Bottom-right

			nid = (row * 2 + 1) * this->connector.nx + col * 2 + 1;
			tntopo += oldarr[nid]; // Top-right

			nid = (row * 2 + 1) * this->connector.nx + col * 2;
			tntopo += oldarr[nid]; // Top-left

			newarr[i] = tntopo / 4; // Average of 4 cells
		}
	}

	/**
	 * Perform grid interpolation (mesh coarsening)
	 * Reduces resolution by factor of 2 in both directions
	 * Uses averaging to preserve overall landscape character
	 */
	void interpolation()
	{
		int nx = floor(this->connector.nx / 2);
		int ny = floor(this->connector.ny / 2);
		int nxy = nx * ny;

		// Create new connector for coarsened grid
		D8connector<fT> ncon = _create_connector(
			nx, ny, this->connector.dx * 2, this->connector.dy * 2, fT(0.), fT(0.));
		ncon.set_default_boundaries(this->boundary_string);

		// Create coarsened topography
		std::vector<fT> ntopo(nxy, 0.);
		this->_interp_helper(this->topography, ntopo, ncon);

		// Commented code would handle variable parameter coarsening

		// Update connector and graph for coarsened grid
		this->connector = std::move(ncon);
		_create_graph(nxy, this->connector, this->graph);

		// Update state with coarsened topography
		this->topography = std::move(ntopo);
		this->_init_vecs();
		this->border20();			 // Apply boundary conditions
		this->fac_resize *= 2; // Update resize factor
	}

	/**
	 * Apply Gaussian smoothing to topography
	 * Reduces small-scale roughness while preserving large-scale features
	 * rr: smoothing radius parameter
	 */
	void smooth(fT rr)
	{
		this->topography = On_gaussian_blur(
			rr, this->topography, this->connector.nx, this->connector.ny);
	}

	/**
	 * Set boundary elevations to zero
	 * Creates a "sea level" boundary condition around the model domain
	 */
	void border20()
	{
		for (int i = 0; i < this->connector.nnodes; ++i) {
			if (this->connector.flow_out_model(i))
				this->topography[i] = 0;
		}
	}

	/**
	 * Calculate chi-star values for landscape analysis
	 * Chi is a transformed coordinate that linearizes river profiles
	 * when the landscape is in steady state with uniform parameters
	 */
	std::vector<fT> _chi_star()
	{
		// Calculate drainage area with constant precipitation
		std::vector<fT> A = this->graph._accumulate_constant_downstream_SFD(
			this->connector.get_area_at_node(0));

		std::vector<fT> chistar(this->graph.nnodes, 0.);
		fT chimax = 0;

		// Calculate chi by integrating upstream from outlets
		for (int i = 0; i < this->graph.nnodes; ++i) {
			int node = this->graph.Sstack[i];						 // Process in stack order
			int rec = this->connector._Sreceivers[node]; // Downstream receiver

			if (node == rec) // Skip outlet nodes
				continue;

			// Chi integration: chi = chi_downstream + distance * A^(-m/n)
			// Using fixed exponent of 0.2 instead of m/n ratio
			chistar[node] = chistar[rec] + this->connector.Sdistance2receivers[node] *
																			 (std::pow(1 / A[node], 0.2));

			// Track maximum chi value for normalization
			if (chistar[node] > chimax)
				chimax = chistar[node];
		}

		// Normalize chi values to range [0,1]
		for (auto& v : chistar)
			v /= chimax;

		return chistar;
	}

	/**
	 * Calculate normalized elevation (z-star) values
	 * Normalizes topography to range [0,1] for analysis
	 */
	std::vector<fT> _z_star()
	{
		std::vector<fT> zstar(this->topography); // Copy topography
		fT zmax = 0;

		// Find maximum elevation
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (zstar[i] > zmax)
				zmax = zstar[i];
		}

		// Normalize by maximum elevation
		for (auto& v : zstar)
			v /= zmax;

		return zstar;
	}

	/**
	 * Apply spatially variable erodibility based on chi values
	 * Modifies erosion rates in areas with chi values between specified bounds
	 * This can simulate lithologic contacts or fault zones
	 */
	void simple_Kfchi(fT tkmod, fT chimin, fT chimax)
	{
		auto chistar = this->_chi_star(); // Calculate normalized chi

		// Switch to variable erodibility mode
		this->Kmod_mode = PARAM_MODE_POP::VARIABLE;
		this->_Kmod = std::vector<fT>(this->graph.nnodes, 1.);

		// Apply modified erodibility where chi is in specified range
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (chistar[i] > chimin && chistar[i] < chimax)
				this->_Kmod[i] = tkmod;
		}

		// Run one iteration to incorporate changes
		this->StSt(1);

		// Reset to constant erodibility mode
		this->Kmod_mode = PARAM_MODE_POP::CONSTANT;
		this->_Kmod = { 1 };
	}

	/**
	 * Apply spatially variable erodibility based on elevation
	 * Similar to simple_Kfchi but uses elevation instead of chi
	 * Useful for simulating elevation-dependent weathering processes
	 */
	void simple_Kfz(fT tkmod, fT zmin, fT zmax)
	{
		auto zstar = this->_z_star(); // Calculate normalized elevation

		// Switch to variable erodibility mode
		this->Kmod_mode = PARAM_MODE_POP::VARIABLE;
		this->_Kmod = std::vector<fT>(this->graph.nnodes, 1.);

		// Apply modified erodibility where elevation is in specified range
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (zstar[i] > zmin && zstar[i] < zmax)
				this->_Kmod[i] = tkmod;
		}

		// Run one iteration to incorporate changes
		this->StSt(1);

		// Reset to constant erodibility mode
		this->Kmod_mode = PARAM_MODE_POP::CONSTANT;
		this->_Kmod = { 1 };
	}

	/**
	 * Normalize topography to range [0,1]
	 * Useful for standardizing elevation data for analysis or visualization
	 */
	void normalise_topography()
	{
		// Find minimum and maximum elevations
		fT mini = std::numeric_limits<fT>::max();
		fT maxi = std::numeric_limits<fT>::min();
		for (auto v : this->topography) {
			if (v < mini)
				mini = v;
			if (v > maxi)
				maxi = v;
		}

		maxi -= mini; // Calculate range

		// Normalize each elevation value
		for (auto& v : this->topography)
			v = (v - mini) / maxi;
	}

	// Accessor methods for current grid dimensions
	int get_active_nx()
	{
		return this->connector.nx;
	} // Current number of columns
	int get_active_ny() { return this->connector.ny; } // Current number of rows

	// Placeholder for potential chi-related calculations
	// std::vector<fT> chiculations()
	// {
	//     // Implementation would go here
	// }
};

// Former Popscape, to be deprecated soon
template<class fT, class Graph_t, class Connector_t>
class popscape_old
{
public:
	// the graph is directly integrated into popscape
	Graph_t graph;
	// And so is the connector as popscape is optimised to pop a final landscape.
	Connector_t connector;

	// Topography
	std::vector<fT> topography, QA;

	// randomiser helper
	std::shared_ptr<DAGGER::easyRand> randu =
		std::make_shared<DAGGER::easyRand>();

	// empty constructor
	popscape_old(){};

	// constructor 1:
	// -> noise type is from the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	popscape_old(RANDNOISE noisetype, int start_nx, int nit, fT dx, fT dy)
	{

		if (start_nx % 8 != 0)
			throw std::runtime_error(
				"target nx and start nx needs to be a multiple of 8");

		// total number of nodes (so far assuming only regular grids)
		int nx = start_nx, ny = start_nx;
		int nxy = nx * ny;

		// init the topo to 0
		this->topography = std::vector<fT>(nxy, 0.);

		// init connector
		this->connector = _create_connector(nx, ny, dx, dy, fT(0.), fT(0.));

		// init graph
		_create_graph(nxy, this->connector, this->graph);

		// init random noise
		if (noisetype == RANDNOISE::WHITE)
			DAGGER::add_noise_to_vector(this->topography, 0., 1.);

		this->connector.set_default_boundaries("periodic_EW");

		for (int i = 0; i < nit; ++i) {
			// std::cout << "REFINING::" << i << std::endl;
			int nnit = (i == 0) ? 50 : 5;
			// std::cout << "REFINING::A" << i << std::endl;
			this->solve_generic(0.45, 1.11, nnit);
			// std::cout << "REFINING::B" << i << std::endl;
			if (i < nit - 1) {
				// std::cout << "REFINING::C" << i << std::endl;

				this->double_res(10, noisetype);
			}
			// std::cout << "done" << std::endl;
		}
	}

	void _init_vecs() { this->QA = std::vector<fT>(this->graph.nnodes, 0.); }

	void solve_generic(fT m, fT n, int n_iterations)
	{
		// std::cout << "REFINING::|||||" <<n_iterations << std::endl;
		// running the code for n_iterations
		for (int nit = 0; nit < n_iterations; ++nit) {
			// std::cout << nit << "|A" << std::endl;

			this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
			this->graph._compute_graph(this->topography, true, false);
			this->_init_vecs();

			this->QA = this->graph._accumulate_constant_downstream_SFD(
				this->connector.get_area_at_node(0));

			// std::cout << nit << "|C" << std::endl;

			for (int i = 0; i < this->graph.nnodes; ++i) {
				int node = this->graph.Sstack[i];
				if (this->connector.flow_out_or_pit(node))
					continue;
				int rec = this->connector._Sreceivers[node];
				this->topography[node] =
					this->topography[rec] + this->connector.Sdistance2receivers[node] *
																		std::pow(1e2, 1. / n) /
																		std::pow(this->QA[node], m / n);
				// this->topography[node] = this->topography[rec] + 1;
			}
			// std::cout << nit << "|D" << std::endl;
		}
		// std::cout << std::endl;
	}

	void double_res(fT noise_intensity, RANDNOISE noisetype)
	{

		int nxy = 4 * this->graph.nnodes;

		int nx = this->connector.nx * 2;
		int ny = this->connector.ny * 2;

		std::vector<fT> ntopo(nxy, 0.);

		if (noisetype == RANDNOISE::WHITE)
			for (int i = 0; i < this->graph.nnodes; ++i) {
				int row, col;
				this->connector.rowcol_from_node_id(i, row, col);
				int nid = row * 2 * nx + col * 2;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
				nid = row * 2 * nx + col * 2 + 1;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
				nid = (row * 2 + 1) * nx + col * 2 + 1;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
				nid = (row * 2 + 1) * nx + col * 2;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
			}

		fT dx = this->connector.dx / 2;
		fT dy = this->connector.dy / 2;
		// init connector
		this->connector = _create_connector(nx, ny, dx, dy, fT(0.), fT(0.));

		// init graph
		_create_graph(nxy, this->connector, this->graph);
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();
	}

	void compute_graph(bool SFD_only)
	{
		this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		this->graph._compute_graph(this->topography, SFD_only, false);
	}

	void compute_DA_SFD()
	{
		this->QA = this->graph._accumulate_constant_downstream_SFD(
			this->connector.get_area_at_node(0));
	}

	void solve_SFD_SPL_imp(fT m, fT n, fT K, fT dt)
	{

		for (int i = 0; i < this->graph.nnodes; ++i) {

			int node = this->graph.Sstack[i];
			int rec = this->connector._Sreceivers[node];

			if (!this->connector.flow_out_or_pit(node) == false)
				continue;

			fT factor = K * dt * std::pow(this->QA[node], m) /
									std::pow(this->connector.Sdistance2receivers[node], n);

			fT ielevation = this->topography[node];
			fT irec_elevation = this->topography[rec];

			fT elevation_k = ielevation;
			fT elevation_prev = std::numeric_limits<fT>::max();
			fT tolerance = 1e-4;
			int nit = 0;
			while (abs(elevation_k - elevation_prev) > tolerance && nit < 10) {
				elevation_prev = elevation_k;
				fT slope =
					std::max(elevation_k - irec_elevation, static_cast<fT>(1e-6));
				fT diff = (elevation_k - ielevation + factor * std::pow(slope, n)) /
									(1. + factor * n * std::pow(slope, n - 1));
				elevation_k -= diff;
				++nit;
			}

			this->topography[node] = elevation_k;
		}
	}

	template<class out_t>
	out_t get_topo()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->topography);
	}
	template<class out_t>
	out_t get_QA()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->QA);
	}

	void apply_uplift(fT dt, fT U)
	{
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (!this->connector.flow_out_or_pit(i))
				this->topography[i] += U * dt;
		}
	}

	template<class out_t>
	void apply_variable_uplift(fT dt, out_t tU)
	{
		auto U = DAGGER::format_input(tU);
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (!this->connector.flow_out_or_pit(i))
				this->topography[i] += U[i] * dt;
		}
	}

	void hydraulic_erosion_v0(int n_particules, fT erosor)
	{
		// old test
	}

	void normalise_topography()
	{
		fT mini = std::numeric_limits<fT>::max();
		fT maxi = std::numeric_limits<fT>::min();
		for (auto v : this->topography) {
			if (v < mini)
				mini = v;
			if (v > maxi)
				maxi = v;
		}
		maxi -= mini;
		for (auto& v : this->topography)
			v = (v - mini) / maxi;
	}
};

// generates a quick and dirty fluvial SFD topo for testing purposes
// uses a multigrid analytical solution to the SPL (see saleve algorithm from
// Steer (2022) on which I applied multigrid methods using random noise and
// cordonnier carving for the projections) final size is 16 * 2^ncycles
// boundaries are "4edges" "periodic_EW" or "periodic_NW"
template<class fT>
std::vector<fT>
_quick_fluvial_topo(int ncycles, std::string boundaries, int nrefine = 10)
{
	// init dx to get a final one = 50
	fT dx = std::pow(2, ncycles) * 50;
	// init connector and boundary conditions
	D8connector<fT> con(16, 16, dx, dx, 0, 0);
	con.set_default_boundaries(boundaries);
	// init graph
	graph<fT, D8connector<fT>> gf(con);
	// init Popscape
	popscape<fT, graph<fT, D8connector<fT>>, D8connector<fT>> psc(gf, con);
	// init randomnoise
	std::vector<fT> topo(16 * 16, 0.);
	add_noise_to_vector(topo, 0, 1);
	psc.set_topo(topo);

	for (int i = 0; i < ncycles + 1; ++i) {
		psc.StSt(nrefine);
		if (i < ncycles)
			psc.restriction(5);
	}

	return psc.topography;
}

template<class fT, class out_t>
out_t
quick_fluvial_topo(int ncycles, std::string boundaries)
{
	std::vector<fT> out = _quick_fluvial_topo<fT>(ncycles, boundaries);
	return DAGGER::format_output<std::vector<fT>, out_t>(out);
}

// End of namespace fastflood
}; // namespace DAGGER

#endif
