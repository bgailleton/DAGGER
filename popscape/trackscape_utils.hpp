#ifndef TRACKSCAPE_UTILS_HPP
#define TRACKSCAPE_UTILS_HPP



namespace DAGGER
{


template<class float_t>
class BasePropStorer
{

public:

	float_t prop = 0;

	BasePropStorer(){;}
	BasePropStorer(float_t prop){this->prop = prop;}

	static BasePropStorer<float_t> create(){return BasePropStorer<float_t>();}
	static BasePropStorer<float_t> create(float_t prop){return BasePropStorer<float_t>(prop);}
	template<class T, class V>
	static void mix(float_t w1, BasePropStorer<T>& prop1, float_t w2, BasePropStorer<V>& prop2)
	{
		if(w1 + w2 == 0) return;
		prop1.prop =( w1 * prop1.prop + w2 * prop2.prop)/(w1 + w2);
	}

};


template<class float_t,class stored_t>
class VerticalStorer
{

public:

	int nnodes = 0;
	float_t dz = 1.;
	std::vector<float_t> topcellz;
	std::vector<std::vector<stored_t> > pile;


	VerticalStorer(){;}
	VerticalStorer(float_t dz, int nnodes)
	{
		this->nnodes = nnodes;
		this->dz = dz;
		this->topcellz = std::vector<float_t>(this->nnodes,0.);
		this->pile = std::vector<std::vector<stored_t> >(this->nnodes, std::vector<stored_t>());
	}


	void pile_up(int i, float_t zadd, stored_t& prop)
	{
		if(zadd == 0)
			return;

		if(this->pile[i].size() == 0)
			this->pile[i].emplace_back(stored_t::create());

		while(zadd > 0)
		{
			float_t nextspace = this->dz - this->topcellz[i];
			float_t tzused = std::min(zadd, nextspace);
			zadd -= tzused;
			stored_t::mix(this->topcellz[i],this->pile[i].back(),tzused, prop);

			this->topcellz[i] += tzused;

			if(zadd > 0)
			{
				this->topcellz[i] = 0;
				this->pile[i].emplace_back(stored_t::create());
			}
		}
	}

	stored_t remove(int i, float_t zrem)
	{
		auto removed = stored_t::create();
		if(zrem == 0)
			return removed;

		float_t cumrem = 0;

		while(zrem>0)
		{
			float_t torem = std::min(zrem,this->topcellz[i]);
			this->topcellz[i] -= torem;
			stored_t::mix(cumrem, removed, torem, this->pile[i].back());
			cumrem += torem;
			if(this->topcellz[i] == 0)
			{
				this->pile[i].pop_back();
				this->topcellz[i] = this->dz;
				if(this->pile[i].size() == 0)
				{
					this->topcellz[i] = 0;
					this->pile[i].emplace_back(stored_t::create());
					return removed;
				}
			}
			zrem -= torem;
		}
		return removed;
	}











};




}


























#endif