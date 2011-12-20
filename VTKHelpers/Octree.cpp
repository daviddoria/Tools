#include "Octree.h"

#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_distance.h>

#include <VXLHelpers/VXLHelpers.h>

#include <limits>

void Octree::Build(const std::vector<vgl_point_3d<double> > &Points, const std::vector<std::vector<unsigned int> > &VertexLists)
{
	//store everything
	/*
	Points_ = Points;
	VertexLists_ = VertexLists;
	*/
	
	//create the points
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		vgl_point_3d<double> P = Points[i];
		points->InsertNextPoint(P.x(), P.y(), P.z());
	}
			
	//create the triangles
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

	for(unsigned int i = 0; i < VertexLists.size(); i++)
	{
		std::vector<unsigned int> Verts = VertexLists[i];
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		triangle->GetPointIds()->SetId(0,Verts[0]);
		triangle->GetPointIds()->SetId(1,Verts[1]);
		triangle->GetPointIds()->SetId(2,Verts[2]);
		triangles->InsertNextCell(triangle);
	}
	
	vtkSmartPointer<vtkPolyData> polydata = vtkPolyData::New();

	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
	
	//create the locator
	Tree = vtkSmartPointer<vtkOBBTree>::New();
	Tree->SetDataSet(polydata);
	Tree->BuildLocator();
	
	Built_ = true;
}


void Octree::Build(vtkSmartPointer<vtkPolyData> polydata)
{
	
	//create the locator
	Tree = vtkSmartPointer<vtkOBBTree>::New();
	Tree->SetDataSet(polydata);
	Tree->BuildLocator();
	
	Built_ = true;
}


bool Octree::IntersectRay(const Ray &R, OrientedPoint &ClosestIntersection) const
{
	if(!Built_)
	{
		std::cout << "Error: Octree not built!" << std::endl;
		assert(0);
	}
		
	vtkPoints* IntersectPoints = vtkPoints::New();
				
	vgl_point_3d<double> Origin = R.getOrigin();
	double P0[3] = {Origin.x(), Origin.y(), Origin.z()};
	double len = 1000.0;
	vgl_point_3d<double> p1 = Origin + len*R.getDirection();
	double P1[3] = {p1.x(), p1.y(), p1.z()};
			
	//Tree->IntersectWithLine(P0, P1, IntersectPoints, NULL);
	vtkIdList* cellIds = vtkIdList::New();
	
	Tree->IntersectWithLine(P0, P1, IntersectPoints, cellIds);
	
	//int result = Tree->IntersectWithLine(P0, P1, IntersectPoints, cellIds);
	//std::cout << "Result: " << result << std::endl;
			
	if(IntersectPoints->GetNumberOfPoints() == 0)
	{
		//std::clog << "No intersections." << std::endl;
		return false;
	}
		
	double ClosestDistance = std::numeric_limits<double>::infinity();
	
	bool intersect = false;
	
	//check which intersection is closest and in the same direction as the ray
	for(unsigned int p = 0; p < static_cast<unsigned int>(IntersectPoints->GetNumberOfPoints()); p++)
	{
		double CurrentCoord[3];
		IntersectPoints->GetPoint(p, CurrentCoord);
		
		vgl_point_3d<double> CurrentPoint(CurrentCoord[0], CurrentCoord[1], CurrentCoord[2]);
		
		Ray r(R.getOrigin(), CurrentPoint - R.getOrigin());
		
		double dotprod = VXLHelpers::dot(r.getDirection(), R.getDirection());
		double d = vgl_distance(CurrentPoint, R.getOrigin());
		
		if( (d < ClosestDistance) && (dotprod >= 0.0) )
		{
			ClosestDistance = d;
			
			/*
			//triangle normal
			unsigned int cellId = cellIds->GetId(p);
			
			std::vector<unsigned int> VertexIds = VertexLists_[p];
			vgl_point_3d<double> V1 = Points_[VertexIds[0]];
			vgl_point_3d<double> V2 = Points_[VertexIds[1]];
			vgl_point_3d<double> V3 = Points_[VertexIds[2]];
			double v1[3] = {V1.x(), V1.y(), V1.z()};
			double v2[3] = {V2.x(), V2.y(), V2.z()};
			double v3[3] = {V3.x(), V3.y(), V3.z()};
			double n[3];
			vtkTriangle::ComputeNormal(v1, v2, v3, n);
			vgl_vector_3d<double> Normal(n[0], n[1], n[2]);
			std::cout << "Normal: " << Normal << std::endl;
			ClosestIntersection = OrientedPoint(CurrentPoint, Normal);
			*/
			
			//!!! fake normal
			ClosestIntersection = OrientedPoint(CurrentPoint, vgl_vector_3d<double>(1.0, 0.0, 0.0));
			
			intersect = true;
		}

	}
	
	if(!intersect)
	{
		std::cout << "There were interesections, but no forward intersections." << std::endl;
	}

	return intersect;
}