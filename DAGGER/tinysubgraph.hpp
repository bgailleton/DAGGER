#pragma once

#include "dodconnector.hpp"
#include "dodcontext.hpp"
#include "utils.hpp"
#include <vector>

namespace DAGGER {

template<class i_t, class f_t, class CON_T, class DATA_T, class PARAM_T>
class TinySubGraph
{
public:
	TinySubGraph(){};
	TinySubGraph(CON_T& con, DATA_T& data, PARAM_T& param)
	{
		this->con = &con;
		this->data = &data;
		this->param = &param;
		this->isDone = std::vector<std::uint8_t>(this->con->nxy(), false);
	};

	~TinySubGraph(){};

	CON_T* con;
	DATA_T* data;
	PARAM_T* param;

	std::vector<i_t> stack;
	std::vector<i_t> nodes;
	std::vector<i_t> baseLevels;
	std::vector<std::uint8_t> isDone;

	template<class CONTAINER_INT>
	void build_simple(CONTAINER& startingNodes)
	{

		this->stack.clear();
		this->nodes.clear();
		this->baseLevels.clear();
		fillvec(this->isDone, false);

		// Initialising a node queue
		std::queue<i_t> tQ;

		// Feeding it witht the starting nodes
		for (auto v : startingNodes) {
			nodes.emplace_back(v);
			isDone[v] = true;
			tQ.emplace(v);
		}

		// Setting up context and helper arrays
		CT_neighbourer_1<i_t, f_t> ctx;
		std::array<i_t, 8> recs;
		std::array<std::uint8_t, 8> recbits;

		while (tQ.empty() == false) {

			// Getting next node and popping it out of the Q
			i_t nextnode = tQ.front();
			tQ.pop();

			// saving the node
			this->nodes.emplace_back(nextnode);

			// reset the node recs and donors
			this->con->reset_node(nextnode);

			// compute receivers only at that specific node
			this->con->__compute_recs_single_node(nextnode, ctx);

			// gathering receivers
			int nr = this->con->Receivers(nextnode, recs);
			this->isDone[nextnode] = nr;

			if (nr == 0) {
				this->baseLevels.emplace_back(nextnode);
				continue;
			}

			// Gathering the receivers into the queue
			for (int i = 0; i < nr; ++i) {
				int trec = recs[i];
				// double checking they are not already in there/processed
				if (this->isDone[trec] == false) {
					tQ.emplace(trec);
					this->isDone[trec] = true;
				}
			}
		}

		// A that point, all the nodes and baselevels for the graph have been
		// gathered And this->isDone[trec] has the number of receivers per node
		//=================

		// Let's invert the receiver codes into donors
		for (auto v : this->nodes)
			this->con->__invert_recs_at_node(v, recbits, recs);

		// Build the stack
		this->stack.reserve(this->nodes.size());

		// the stack will start with the base levels
		for (auto v : this->baseLevels)
			tQ.emplace(v);

		while (tQ.empty() == false) {
			// Getting next node and popping it out of the Q
			int nextnode = tQ.front();
			tQ.pop();
			// if the node is inj Q is ready 4 stack
			this->stack.emplace_back(nextnode);

			// Grabbing its donors subgraphically speaking
			int nd = this->con->Donors(nextnode, recs);
			for (int i = 0; i < nd; ++i) {
				// each donors is visited
				int tdon = recs[i];
				// and I keep track of the number of visits
				--this->isDone[tdon];
				// if the donor has been visited by all its receivers, then it is ready
				if (this->isDone[tdon] == 0) {
					tQ.emplace(tdon);
				}
			}
		}

		// Done
	}
};

}
