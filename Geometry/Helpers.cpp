#include "Helpers.h"

std::vector<vgl_point_3d<double> > GetOPCoords(const std::vector<OrientedPoint> &V)
{
	std::vector<vgl_point_3d<double> > Coords(V.size());
	
	for(unsigned int i = 0; i < V.size(); i++)
	{
		Coords[i]	= V[i].getCoord();
	}
	
	return Coords;
}

std::vector<Color<unsigned char> > GetOPColors(const std::vector<OrientedPoint> &V)
{
	std::vector<Color<unsigned char> > Colors(V.size());
	
	for(unsigned int i = 0; i < V.size(); i++)
	{
		Colors[i]	= V[i].getColor();
	}
	
	return Colors;
}