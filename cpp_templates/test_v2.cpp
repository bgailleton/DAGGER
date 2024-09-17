#include "D8connector.hpp"
#include "databag.hpp"
#include "dodconnector.hpp"
#include "enumutils.hpp"
#include "graph.hpp"
#include "graphflood_v2.hpp"
#include "lookup_neighbourer.hpp"
#include "popscape.hpp"
#include "utils.hpp"

using namespace DAGGER;

int
main(int argc, char const* argv[])
{

	ocarina timer;

	int nx = 1001;
	int ny = 1068;
	double dx = 1.;
	double dy = 1.;
	std::string yolo{ "4edges" };

	std::vector<double> topotest =
		// DAGGER::_quick_fluvial_topo<double>(6, yolo, 10);
		// load_npy<double>("topo1024.npy");
		// load_npy<double>("topo1024_rdnoise.npy");
		load_npy<double>("green_river.npy");

	// for (auto& v : topotest)
	// 	v *= 0.1;

	// save_vec_to_1Dnpy("topo1024.npy",1024,1024, topotest);

	// D8connector<double> clacon(nx, ny, dx, dy, 0., 0.);

	// timer.tik();
	// clacon.update_links_from_topo(topotest);
	// timer.tok("clacon took");

	// graph<double, D8connector<double>> gra(clacon);
	// gra.set_LMR_method(DEPRES::priority_flood);

	// timer.tik();
	// gra._compute_graph(topotest, false, false);
	// timer.tok("clagra took");

	// timer.tik();
	// gra.topological_sorting_SF();
	// timer.tok("clastack");

	Hermes<int, double> dbag;
	dbag._surface = std::move(topotest);
	Connector8<int, double> con(nx, ny, dx, dy, dbag);
	con.boutype = CONBOU::EDGES;
	con.flowtopo = CONFLOWTOPO::ALL;

	con.init();

	// timer.tik();
	// con.PFcompute_all();
	// timer.tok("PFcom took");

	std::cout << "init gf2 " << std::endl;

	Graphflood2<int, double, decltype(con), int, decltype(dbag)> gf(
		con, nx, dbag);
	gf.init();

	std::cout << "init entry points" << std::endl;
	gf.Prate = 1e-4;
	gf._computeEntryPoints_prec(1);

	std::cout << "done, I have " << gf.input_node_Qw.size()
						<< " points, now running gf2" << std::endl;
	gf.dt = 2e-3;

	for (int i = 0; i < 100; ++i) {
		std::cout << "run " << i << std::endl;
		gf.run_subgraphflood();
		// save_vec_to_2Dnpy("hw" + std::to_string(i) +".npy", con._nx, con._ny,
		// dbag._hw);

		// std::cout << "ran " << i << "\r";
	}

	// timer.tik();
	// con.compute();
	// timer.tok("newcon took");

	// timer.tik();
	// con._compute_all_exp1();
	// timer.tok("newcon prefetch took");

	// timer.tik();
	// con._quickLM();
	// timer.tok("labelling took");

	// timer.tik();
	// con._compute_all_exp2();
	// timer.tok("newcon prefetch 2 took");

	save_vec_to_2Dnpy("topt.npy", con._nx, con._ny, dbag._surface);
	save_vec_to_2Dnpy("hwtest.npy", con._nx, con._ny, dbag._hw);
	save_vec_to_2Dnpy("Qwtest.npy", con._nx, con._ny, dbag._Qwin);
	save_vec_to_2Dnpy("debugtest.npy", con._nx, con._ny, dbag._debug);

	// // timer.tik();
	// // con._quickSstack();
	// // timer.tok("newconSstack took");

	// std::vector<double> dA(con.nxy(), 0.);

	// for (int i = con.nxy() - 1; i >= 0; --i) {
	// 	int node = dbag._Sstack[i];
	// 	int rec = con.Sreceivers(node);
	// 	dA[node] += con.area(node);
	// 	if (rec != node)
	// 		dA[rec] += dA[node];
	// }

	// save_vec_to_2Dnpy("DA.npy", con._nx, con._ny, dA);

	// Un comment to check the lookup

	// std::array<int,8> arr;
	// int nn = 0;
	// std::cout << "Top left (11) -> ";
	// nn = con.data->LK8.NeighbourerNN[11];
	// arr = con.data->LK8.Neighbourer[11];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "Top (31) -> ";
	// nn = con.data->LK8.NeighbourerNN[31];
	// arr = con.data->LK8.Neighbourer[31];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "TopRight (22) -> ";
	// nn = con.data->LK8.NeighbourerNN[22];
	// arr = con.data->LK8.Neighbourer[22];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "Left (107) -> ";
	// nn = con.data->LK8.NeighbourerNN[107];
	// arr = con.data->LK8.Neighbourer[107];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "Right (214) -> ";
	// nn = con.data->LK8.NeighbourerNN[214];
	// arr = con.data->LK8.Neighbourer[214];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "Bleft (104) -> ";
	// nn = con.data->LK8.NeighbourerNN[104];
	// arr = con.data->LK8.Neighbourer[104];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "Bottom (248) -> ";
	// nn = con.data->LK8.NeighbourerNN[248];
	// arr = con.data->LK8.Neighbourer[248];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "BottomRight (208) -> ";
	// nn = con.data->LK8.NeighbourerNN[208];
	// arr = con.data->LK8.Neighbourer[208];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	// std::cout << "Internal (255) -> ";
	// nn = con.data->LK8.NeighbourerNN[255];
	// arr = con.data->LK8.Neighbourer[255];
	// for(int j=0;j<nn;++j)
	// 	std::cout << " | " <<  arr[j];
	// std::cout << std::endl;

	std::cout << "done" << std::endl;
	return 0;
}
