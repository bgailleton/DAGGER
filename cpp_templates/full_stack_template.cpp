/*
This template provides the minimum example to use DAGGER to:
1) build a graph from a topo
2) process the node in the stack order (or inverse)
3) it also show how to fecth donors/receivers and relevant informations

Note that it generates a random mature fluvial topography using popscape but you
can replace that part feeding your own vector

To compile with g++/clang: g++ -O3 -std=c++17 -o stacker.exe -I../DAGGER
-I../popscape

DAGGER is header-only, and can be built without any dependencies, so you simply
need to make sure the includes can be found by your compiling toolchain. Here,
for simplicity, I use relative paths.

B.G

*/

#include "D8connector.hpp"
#include "graph.hpp"
#include "popscape.hpp" // only required by _quick_fluvial_topo

int
main(int argc, char const* argv[])
{

	// Step 0: setting up topo/geometry, this can be replaced by your own loading
	// method # DAGGER works with 1D vector which for regular grid are row major
	// vectorisation (node_index = nx * row_index + col_index)

	// # Setting up the geometry
	int nx = 512;
	int ny = 512;
	double dx = 30.;
	double dy = 30.;
	double xmin = 0.;
	double ymin = 0.;

	// # Boundary string: can be periodic_EW, periodic_NS or 4edges (for the 4
	// borders) # Note: more sophisticated boundaries (e.g. with nodata or forcing
	// in/out from specific nodes) is possible via a vector of boundary codes and
	// is explained in another example
	std::string boundary = "periodic_EW";

	std::cout << "Generating topography using Salève" << std::endl;
	// # THIS STEP CAN BE REPLACED BY ANY PROCESS THAT LOAD/CREATE A 1D VECTOR
	// CONTAINING YOUR TOPO # this particular line generates a 512*512 steady
	// state fluvial topography
	std::vector<double> topography =
		DAGGER::_quick_fluvial_topo<double>(5, boundary);

	std::cout << "Initialising data structure" << std::endl;
	// Step 1: initialising the connector and graph objects

	// # Creating the connector, the class managing all the neighbouring and
	// boundaries operations
	DAGGER::D8connector<double> connector(nx, ny, dx, dy, xmin, ymin);
	// # Setting the boundaries
	connector.set_default_boundaries(boundary);

	// # initialising the graph for this connector, graph is an adaptator on grid
	// connectors managing local minimas, topological ordering operations as well
	// as traversal algorithms (e.g. anything accumulating variable upstream or
	// downstream).
	DAGGER::graph<double, decltype(connector)> gf(connector);

	std::cout << "Computing graph" << std::endl;
	// Step 2: Computing the graph

	// # first let's pick a local minima resolver, options are:
	// ## DAGGER::DEPRES::none: internal pit are unresolved
	// ## DAGGER::DEPRES::cordonnier_fill: algorithm 4 from
	// https://esurf.copernicus.org/articles/7/549/2019/
	// ## DAGGER::DEPRES::cordonnier_carve: algorithm 3 from
	// https://esurf.copernicus.org/articles/7/549/2019/
	// ## DAGGER::DEPRES::priority_flood: BArnes 2014 applied on topography and
	// THEN building graph on it
	// ## DAGGER::DEPRES::priority_full_MFD: experimental and not battle tested
	// but fastest: builds the stack and solve the local minima in one go mixing
	// multiple methods
	gf.set_LMR_method(DAGGER::DEPRES::cordonnier_carve);

	// # Computing the graph:
	// ## takes the topo and 2 switches and returns the preprocessed topography
	// ## first switch only computes the Single Flow graph if true
	// ## second determine if the stack is computed by sorting elevation (true) or
	// topological sorting (false)
	// ## IMPORTANT, if the topography changes, you can simply update the graph
	// with this function, no need to rebuild it from scratch
	std::vector<double> PP_topo =
		gf.compute_graph<std::vector<double>, std::vector<double>>(
			topography, false, false);

	std::cout << "Running analysis" << std::endl;
	// Step 3: traverse stack and do stuff - this is a template so I don;t do much
	// here # the neighbouring system works by caching neighbours in a std::array,
	// I benchmarked it and it is significantly faster than re-initialising a
	// vector of neighbour each time # getting a cache of right size:
	auto neighbours = connector.get_empty_neighbour();

	// # Iterating through the stack and through the neighbours ban be done in
	// multiple ways # First one need to select the stack: gf.stack is the
	// multiple flow one and gf.Sstack is the single flow one
	// ## only the SFD one is computed is the compute_graph only_SFD switch is
	// true
	// ## There are few differences: SFD is faster to compute (understandably) but
	// also better organised as watersheds are groupped by chunks (which is not
	// possible in MFD as a node belongs to multiple watersheds) # Iterating
	// through the stacks in normal order goes from downstream to upstream and in
	// inverse order, well the opposite.

	// # normal order
	// for(int i = 0; i<connector->nnodes; ++i)
	// # reverse order
	for (int i = connector.nnodes - 1; i >= 0; --i) {
		// ## Getting the node index
		int node = gf.stack[i]; //  or gf.Sstack[i] for SFD

		// ## Getting SFD receiver is straighforward: int rec =
		// connector.Sreceivers[node]; double dx =
		// connector.Sdistance2receivers[node];

		// ## Getting MFD neighbous:

		// ## Now we want to get the neighbours. This also can be done in many ways
		// depending on the needs

		// ### Node based neighbouring: feeds in place the neighbour array with nn
		// neighbouring nodes
		int nn =
			connector.get_neighbour_idx(node, neighbours); // get ALL the neighbours
		// int nn = connector.get_receivers_idx(node, neighbours); // get the
		// receivers int nn = connector.get_donorss_idx(node, neighbours); // get
		// the donors
		// ### iterating through neighbours
		for (int j = 0; j < nn; ++j) {
			int other_node = neighbours[j];

			// common checks:
			if (connector.boundaries.no_data(other_node)) {
				// node is no data
			}

			if (connector.flow_out_or_pit(other_node)) {
				// flow stops here (unprocessed pit) or escapes the model
			}

			if (connector.flow_out_model(other_node)) {
				// flow escapes the model
			}

			// Ecample of node-based characteristic: surface represented by the cell
			double area = connector.get_area_at_node(other_node);
			// do stuff with the node, for example transferring flux to receivers,
			// checking relationship or whatever Node based neighbouring is efficient
			// for traversal requiring only node location. If more geometry is needed,
			// like distance, perpendiculat distance, ... see link based neighbouring
			// below
		}

		// ### link based neighbouring: feeds in place the neighbour array with nn
		// neighbouring links
		nn = connector.get_neighbour_idx_links(
			node, neighbours); // get ALL the neighbouring links
		// int nn = connector.get_receivers_idx_links(node, neighbours); // get the
		// links to receivers int nn = connector.get_donorss_idx_links(node,
		// neighbours); // get the links to donors
		// ### iterating through links
		for (int j = 0; j < nn; ++j) {
			int lix = neighbours[j]; // link id

			// example of link-based characteristic: dx adn dy
			double tdx = connector.get_dx_from_links_idx(
				lix); // distance represented by the link
			double tdy = connector.get_traverse_dx_from_links_idx(lix);

			// it is of course possible to get the other node id
			int other_node = connector.get_other_node_from_links(lix, node);
			// Other manipulations are possible and may be relevant when iterating
			// through links without node knowledge:
			int receiver = connector.get_to_links(lix);
			int donor = connector.get_from_links(lix);
			int A, B;
			connector.from_to_from_link_index(lix, A, B); // feeds A, B in place.
			// do stuff with it, for example transferring flux to receivers, checking
			// relationship or whatever Node based neighbouring is efficient for
			// traversal requiring only node location. If more geometry is needed,
			// like distance, perpendiculat distance, ... see link based neighbouring
			// below
		}
	}

	std::cout << "Done" << std::endl;

	// BONUS: uncomment to save a vector of node size to a numpy compatible .npy
	// file (using topography as example here)
	DAGGER::save_vec_to_2Dnpy("topo.npy", nx, ny, PP_topo);

	return 0;
}
