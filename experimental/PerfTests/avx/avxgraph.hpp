#ifndef GRAPHAVX_HPP
#define GRAPHAVX_HPP
#include "avxcon.hpp"



template<class Connector_t>
class graphavx
{
public:
	
	Connector_t* connector;
	
	
	graphavx(){;};

	graphavx(Connector_t& con)
	{
		this->connector = &con;
	}


	
};