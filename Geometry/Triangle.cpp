#include "Triangle.h"

#include <algorithm>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

#include <VXLHelpers/VXLHelpers.h>

//////////////// Operators /////////////////
ostream& operator << (ostream &output, const Triangle &T)
{
  output << "Triangle" << endl
	 << "--------" << endl
	 << "P1: " << T.getVertex(0) << endl
	 << "P2: " << T.getVertex(1) << endl
	 << "P3: " << T.getVertex(2) << endl;

  return output;
}

bool Triangle::PointInside(vgl_point_3d<double> &P) const
{
	 return vgl_triangle_3d_test_inside(P, Vertices_[0], Vertices_[1], Vertices_[2]);
}

double Triangle::DistanceToPoint(vgl_point_3d<double> &P) const
{
	return vgl_triangle_3d_distance(P , Vertices_[0], Vertices_[1], Vertices_[2]);
}
