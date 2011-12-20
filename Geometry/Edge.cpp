#include "Edge.h"

bool Edge::ContainsPoint(const vgl_point_3d<double> &P)
{
	if((P1 == P) || (P2 == P))
		return true;
	else	
		return false;
}

bool operator==(const Edge &E1, const Edge &E2)
{
	//same edge
	if( (E1.P1 == E2.P1) && (E1.P2 == E2.P2) )
		return true;
	else if( (E1.P1 == E2.P2) && (E1.P2 == E2.P1) )
		return true;
	else
		return false;
}


void AddEdge(std::vector<Edge> &Edges, const Edge &E)
{
	if(!ContainsEdge(Edges, E))
	{
		Edges.push_back(E);
	}
}

bool ContainsEdge(const std::vector<Edge> &Edges, const Edge &E)
{
	for(unsigned int i = 0; i < Edges.size(); i++)
		if(Edges[i] == E)
			return true;
	
	return false;
}