#ifndef OCTREE_H
#define OCTREE_H

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkOBBTree.h>
#include <vtkPolyData.h>

#include <Geometry/Ray.h>
#include <Geometry/OrientedPoint.h>

#include <vector>

#include <vgl/vgl_point_3d.h>

class Octree
{
	private:

		vtkSmartPointer<vtkOBBTree> Tree;
		/*
		std::vector<vgl_point_3d<double> > Points_;
		std::vector<std::vector<unsigned int> > VertexLists_;
		*/
	public:
		////////// Constructors //////
		Octree(){Built_ = false;}
		
		///////// Member Variables /////////
		bool Built_;
		
		////////// Functions ///////////
		void Build(const std::vector<vgl_point_3d<double> > &Points, const std::vector<std::vector<unsigned int> > &VertexLists);
		void Build(vtkSmartPointer<vtkPolyData> polydata);
		bool IntersectRay(const Ray &R, OrientedPoint &Intersection) const;
};


#endif
