/*
this file contain the algorithm(s) related to hillshading a dem from a connector object
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef HILLSHADING_HPP
#define HILLSHADING_HPP


#include "utils.hpp"
#include "D8connector.hpp"
#include "graph.hpp"

#include "wrap_helper.hpp"



namespace DAGGER
{

template<class Connector_t,class topo_t, class out_t, class float_t>
out_t hillshade(Connector_t& connector, topo_t& ttopography)
{
	auto topography = format_input(ttopography);
	float_t altitude = 45;
	float_t azimuth = 315;
	float_t z_factor = 1;
	float_t pi = 3.1415926;

	// std::vector<float_t> hillshade(ptr, ptr + connector.nnodes);
	std::vector<float_t> hillshade(connector.nnodes,0.);


	//convert zenith and azimuth into radians for calculation
	float_t zenith_rad = (90 - altitude) * pi / 180.0;
	float_t azimuth_math = 360-azimuth + 90;
	if (azimuth_math >= 360.0) azimuth_math = azimuth_math - 360;
	float_t azimuth_rad = azimuth_math * pi /180.0;

	for(int i = 0; i<connector.nnodes; ++i)
	{
		// Ignoring no data
		if(connector.boundaries.no_data(i))
			continue;

		float_t slope_rad = 0;
		float_t aspect_rad = 0;
		float_t dzdx = 0;
		float_t dzdy = 0;

		float_t ij = topography[i];
		int j = connector.get_right_neighbour(i).node;
		float_t ijp1 = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_bottom_neighbour(i).node;
		float_t ip1j = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_bottomright_neighbour(i).node;
		float_t ip1jp1 = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_topleft_neighbour(i).node;
		float_t im1jm1 = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_top_neighbour(i).node;
		float_t im1j = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_topright_neighbour(i).node;
		float_t im1jp1 = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_left_neighbour(i).node;
		float_t ijm1 = (connector.is_in_bound(j)) ? topography[j]:ij;
		j = connector.get_bottomleft_neighbour(i).node;
		float_t ip1jm1 = (connector.is_in_bound(j)) ? topography[j]:ij;


		if (ij > 0 )
		{
			dzdx = ((ijp1 + 2*ip1j + ip1jp1) - (im1jm1 + 2*im1j + im1jp1)) / (8 * connector.dx);
			dzdy = ((im1jp1 + 2*ijp1 + ip1jp1) - (im1jm1 + 2*ijm1 + ip1jm1)) / (8 * connector.dy);
			slope_rad = atan(z_factor * sqrt((dzdx*dzdx) + (dzdy*dzdy)));
			if (dzdx != 0)
			{
				aspect_rad = std::atan2(dzdy, (dzdx*-1));
				if (aspect_rad < 0) aspect_rad = 2*pi + aspect_rad;
			}
			else
			{
				if (dzdy > 0) aspect_rad = pi/2;
				else if (dzdy < 0) aspect_rad = 2 * pi - pi/2;
			}

			hillshade[i] = ((std::cos(zenith_rad) * std::cos(slope_rad)) + (std::sin(zenith_rad) * std::sin(slope_rad) * std::cos(azimuth_rad - aspect_rad)));
			// std::cout << hillshade[i] << "|";
			if (hillshade[i] < 0) hillshade[i] = 0;
		}

	}

	return format_output<decltype(hillshade), out_t>(hillshade);

}




// End of namespace DAGGEr	
}









































//end of definition if
#endif