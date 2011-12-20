#ifndef EDGE_H
#define EDGE_H

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_distance.h>

#include <vector>

#include <iostream>
#include <assert.h>

class Edge
{
	public:
		Edge() {}
		Edge(const vgl_point_3d<double> &p1, const vgl_point_3d<double> &p2) : P1(p1), P2(p2) {}
		vgl_point_3d<double> P1;
		vgl_point_3d<double> P2;
		double Length(){return vgl_distance(P1, P2);	}
		bool ContainsPoint(const vgl_point_3d<double> &P);
};

bool operator==(const Edge &E1, const Edge &E2);

void AddEdge(std::vector<Edge> &Edges, const Edge &E);
bool ContainsEdge(const std::vector<Edge> &Edges, const Edge &E);		
#endif