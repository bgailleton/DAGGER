#include "lookup_neighbourer.hpp"
#include "databag.hpp"
#include "utils.hpp"
#include "dodconnector.hpp"

using namespace DAGGER;

int main(int argc, char const *argv[]){

	int nx = 512;
	int ny = 512;
	float dx = 30.;
	float dy = 30.;
	Hermes<int, float> dbag;
	Connector8<int, float> con(nx,ny,dx,dy, dbag);
	con.init();

	std::cout << "done" << std::endl;
	return 0;
}