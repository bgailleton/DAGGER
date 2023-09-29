#include "D8connector.hpp"
#include "databag.hpp"
#include "dodconnector.hpp"
#include "enumutils.hpp"
#include "graph.hpp"
#include "lookup_neighbourer.hpp"
#include "popscape.hpp"
#include "utils.hpp"

using namespace DAGGER;

int
main(int argc, char const* argv[])
{

	ocarina timer;

	int nx = 1024;
	int ny = 1024;
	double dx = 30.;
	double dy = 30.;
	std::string yolo{ "periodic_EW" };
	std::vector<double> topotest =
		DAGGER::_quick_fluvial_topo<double>(6, yolo, 10);

	D8connector<double> clacon(nx, ny, dx, dy, 0., 0.);

	timer.tik();
	clacon.update_links_from_topo(topotest);
	timer.tok("clacon took");

	graph<double, D8connector<double>> gra(clacon);

	timer.tik();
	gra.topological_sorting_SF();
	timer.tok("clastack");

	Hermes<int, double> dbag;
	dbag._surface = std::move(topotest);
	Connector8<int, double> con(nx, ny, dx, dy, dbag);
	con.boutype = CONBOU::PEW;
	con.flowtopo = CONFLOWTOPO::MFD;

	con.init();

	timer.tik();
	con.compute();
	timer.tok("newcon took");

	// save_vec_to_2Dnpy("neighbours.npy", con._nx, con._ny, dbag._neighbours);

	// timer.tik();
	// con._quickSstack();
	// timer.tok("newconSstack took");

	// std::vector<double> dA(con.nxy(), 0.);

	// for(int i=con.nxy()-1;i >= 0; --i){
	// 	int node = dbag._Sstack[i];
	// 	int rec = con.Sreceivers(node);
	// 	dA[node] += con.area(node);
	// 	if(rec != node)
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
