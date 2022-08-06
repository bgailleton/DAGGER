#include "../../DAGGER/D8connector.hpp"
#include "../../DAGGER/graph.hpp"
#include "../../DAGGER/utils.hpp"




int main(int argc, char const *argv[])
{
	DAGGER::easyRand er = DAGGER::easyRand();
	DAGGER::graph<double> G(250000,8);
	DAGGER::D8connector<double> con(500, 500, 30, 30, 0, 0);
	std::vector<double> topo(250000);
	for(auto& v:topo) v=er.get();

	G.init_graph(con);
	std::vector<double> PPt = G.compute_graph<DAGGER::D8connector<double>, std::vector<double>, std::vector<double> >("carve",topo,con,false);





	return 0;
}