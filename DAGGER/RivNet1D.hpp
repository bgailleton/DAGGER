#pragma once

#include "enumutils.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class Connector8;
template<class i_t, class f_t>
class Hermes;

template<class i_t, class f_t>
class RivNet1D
{

public:
	RivNet1D(Connector8<i_t, f_t>& con, Hermes<i_t, f_t>& dbag)
	{
		this->con = &con;
		this->data = &dbag;
	}

	// Connector adn data bag
	Connector8<i_t, f_t>* con;
	Hermes<i_t, f_t>* data;

	// N nodes size data
	std::vector<i_t> nodes;
	std::vector<i_t> source_keys;

	// Specific vectors
	std::vector<i_t> sources;

	void reset()
	{
		this->nodes.clear();
		this->source_keys.clear();
		this->sources.clear();
	}

	void build_from_sources(std::vector<int>& sources)
	{

		if (this->data->Sstack.size() == 0)
			throw std::runtime_error(
				"Single flow stack needed for river extraction yo");

		this->reset();

		this->sources = sources;

		std::vector<std::uint8_t> isDone(this->con->nxy(), false);

		int tsk = -1;
		for (auto source : this->sources) {
			++tsk;
			int node = source;
			while (true) {
				if (isDone[node])
					break;

				this->nodes.emplace_back(node);

				this->source_keys.emplace_back(tsk);

				isDone[node] = true;

				if (can_out(this->data->_boundaries[node]))
					break;

				node = this->con->Sreceivers(node);
			}
		}

		return;
	}
};

};
