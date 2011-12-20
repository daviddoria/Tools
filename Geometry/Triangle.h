#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_triangle_3d.h>

#include <vector>

#include <iostream>
#include <assert.h>

class Triangle
{
	private:
		std::vector<vgl_point_3d<double> > Vertices_;
		vgl_plane_3d<double> Plane_;
	
	public:
	
	Triangle() : Vertices_(3){}
		
	Triangle(const std::vector<vgl_point_3d<double> > &points)
	{
		assert(points.size() == 3);
		Vertices_ = points;
		Plane_ = getPlane();
	}
	
	Triangle(const vgl_point_3d<double> &p0, const vgl_point_3d<double> &p1, const vgl_point_3d<double> &p2) : Vertices_(3)
	{
		Vertices_[0] = p0;
		Vertices_[1] = p1;
		Vertices_[2] = p2;
		
		Plane_ = getPlane();
	}
	
	///////////// Accessors /////////
	vgl_point_3d<double> const getVertex(const int ind) const
	{
		assert((ind >= 0) && (ind <= 2));
		return Vertices_[ind];
	}
	
	vgl_plane_3d<double> const getPlane() const {return vgl_plane_3d<double> (Vertices_[0], Vertices_[1], Vertices_[2]);}
	
	vgl_vector_3d<double> getNormal() const {return Plane_.normal();}
	
	//////////// Functions ///////////
	void Draw(bool Normals, bool InMesh) const;
	bool PointInside(vgl_point_3d<double> &P) const;
	double DistanceToPoint(vgl_point_3d<double> &P) const;

};

/////////external operators/////////////
std::ostream & operator << (std::ostream &output, const Triangle &T);

#endif
