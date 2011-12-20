#include "Triangulate.h"

#include <algorithm>

#include <Geometry/Angles.h>
#include <Geometry/Geometry.h>

#include <vgl/vgl_point_2d.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDelaunay2D.h>
#include <vtkSmartPointer.h>


void TriangulateLidar(const vgl_point_3d<double> &ScannerLocation, const vector<vgl_point_3d<double> > &Points, vector<vector<unsigned int> > &TriangleVertexList, const double MaxSize)
{
	
	vtkSmartPointer<vtkPoints> points2D = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> points3D = vtkSmartPointer<vtkPoints>::New();
	
	/*
	//convert the 3d points to a grid of theta/phi values
	for ( unsigned int i = 0; i < Points.size(); i++ )
	{
    		//convert to polar grid
		double theta, phi;
		vgl_vector_3d<double> Direction = Points[i] - ScannerLocation;
		//Rect2Sphere(vgl_point_to_vgl_vector(Points[i]), phi, theta);
		Rect2Sphere(Direction, phi, theta);
		
		points2D->InsertNextPoint(theta, phi, 0);
		
		points3D->InsertNextPoint(Points[i].x(), Points[i].y(), Points[i].z());
	}
	*/
	
	vector<vgl_point_2d<double> > Points2d = geom::ScannerImage(Points, ScannerLocation);

	for ( unsigned int i = 0; i < Points.size(); i++ )
	{
    		points2D->InsertNextPoint(Points2d[i].x(), Points2d[i].y(), 0);
		
		points3D->InsertNextPoint(Points[i].x(), Points[i].y(), Points[i].z());
	}
	
	//add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points2D);
	
	//triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInput(polydata);
	delaunay->Update();

	//apparently this sets the 3d points into the triangles instead of the 2d points?
	vtkSmartPointer<vtkPolyData> result = delaunay->GetOutput();

	cout << "NumCells: " << result->GetNumberOfCells() << endl;
	cout << "NumPoints2d: " << result->GetNumberOfPoints() << endl;
	
	result->SetPoints(points3D);

	cout << "NumPoints3d: " << result->GetNumberOfPoints() << endl;
	
	//get the resulting triangles from the triangulation
	vtkSmartPointer<vtkCellArray> cells = result->GetPolys();

	
	cout << "NumTris in Delaunay: " << result->GetNumberOfPolys() << endl; //the result is not stored as Polys, so don't use this!
	cout << "NumTris in Delaunay: " << result->GetNumberOfCells() << endl; 
	
	vtkIdType npts;
	vtkIdType * pts;

	//go through all the triangles of the Delaunay Triangulation
	cells->InitTraversal();
	while (cells->GetNextCell(npts,pts))
	{
		//get the three points of the triangle
		vgl_point_3d<double> x = Points[pts[0]];
		vgl_point_3d<double> y = Points[pts[1]];
		vgl_point_3d<double> z = Points[pts[2]];
		
		//find the length of all the edges
		vector<double> dists(3);
		dists[0] = (x-y).length();
		dists[1] = (x-z).length();
		dists[2] = (y-z).length();

		//throw away triangles that are bigger than a threshold
		//sort the edge lengths
		sort(dists.begin(),dists.end());
		//if the longest edge is larger than the tolerance, don't include this triangle
		if (dists[2] > MaxSize)
			continue;
		
		//newcells->InsertNextCell(npts,pts);

     	///ADD IT AS A VALID TRIAGLE
		vector<unsigned int> V;
		V.push_back(pts[0]);
		V.push_back(pts[1]);
		V.push_back(pts[2]);
		TriangleVertexList.push_back(V);
	}
	
	//result->SetPolys(newcells);

}
