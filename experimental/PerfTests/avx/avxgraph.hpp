#ifndef GRAPHAVX_HPP
#define GRAPHAVX_HPP
#include "avxcon.hpp"



template<class Connector_t>
class graphavx
{
public:
	
	Connector_t* connector;
	std::vector<int> stack;
	
	graphavx(){;};

	graphavx(Connector_t& con)
	{
		this->connector = &con;
		this->stack = std::vector<int>(this->connector->nxy);
	}

	// This is my implementation of Braun and willett 2014
	// Slightly modified:
	// - First it is based on the fastscapelib version, which bypass the need of the delta vectors and all by simply applying the recursion directly feeding the stack
	// - Secondly, recursion is not the best practice for most languages so I am using a stack data structure instead
	// Results are identical and performances similar (marginally better, but that is linked to c++ not being a heavy recursion friendly language)
	void toposort_SFD()
	{
		// The stack container helper
		std::stack<int, std::vector<int> > stackhelper;
		// std::vector<bool> isdone(this->nnodes,false);
		// going through all the nodes
		std::array<int,8> arr;
		int istack = 0;
		for(int i=0; i<this->connector->nxy; ++i)
		{
			// if they are base level I include them in the stack
			if(this->connector->Sreceivers(i) == i)
			{
				stackhelper.emplace(i);
				// ++istack;
			}

			// While I still have stuff in the stack helper
			while(stackhelper.empty() == false)
			{
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();stackhelper.pop();
				this->stack[istack] = nextnode;
				++istack;
				// std::cout << istack << "/" << this->connector->nxy << std::endl;
				int nn = this->connector->SDonors(nextnode,arr);
				// as well as all its donors which will be processed next
				for( int j = 0; j < nn; ++j)
				{
					stackhelper.emplace(arr[j]);
				}

			}

		}
	}

	std::vector<float> _compute_drainage_area()
	{
		std::vector<float> tA(this->connector->nxy,0.);
		for(int i=this->connector->nxy-1;i>=0; --i)
		{
			int node = this->stack[i];
			tA[node] += this->connector->get_cell_area(node);
			int rec = this->connector->Sreceivers(node);
			if(rec != node)
			{
				tA[rec] += tA[node];
			}

		}
		return tA;
	}

	void save_drainage_area(std::string filename)
	{
		auto tA = this->_compute_drainage_area();
		const std::vector<std::uint64_t> shape{std::uint64_t(this->connector->ny), std::uint64_t(this->connector->nx)};
		const bool fortran_order{false};
		const std::string path{filename};
		npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), tA);
	}
	
};


















































#endif