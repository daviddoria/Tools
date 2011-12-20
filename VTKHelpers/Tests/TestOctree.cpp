#include <iostream>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>


#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

#include <VTKHelpers/VTKHelpers.h>
#include <VTKHelpers/Octree.h>

void TestOctree();
void TestOctree(const std::string &Filename);
		
int main()
{
	TestOctree();
	TestOctree("car.vtp");
	
	return 0;
}

void TestOctree()
{
	std::vector<vgl_point_3d<double> > Points;
	Points.push_back(vgl_point_3d<double>(1.0, 0.0, 0.0));
	Points.push_back(vgl_point_3d<double>(0.0, 0.0, 0.0));
	Points.push_back(vgl_point_3d<double>(0.0, 1.0, 0.0));
	
	std::vector<std::vector<unsigned int> > VertexLists;
	std::vector<unsigned int> VertexList;
	VertexList.push_back(0);
	VertexList.push_back(1);
	VertexList.push_back(2);
	VertexLists.push_back(VertexList);
	
	Octree Tree;
	Tree.Build(Points, VertexLists);
	
	OrientedPoint Intersection;
	Ray R(vgl_point_3d<double>(0.5, 0.5, -1.0), vgl_vector_3d<double> (0.0, 0.0, 1.0));
	bool intersect = Tree.IntersectRay(R, Intersection);
	
	std::cout << "Intersect? " << intersect << std::endl;
	if(intersect)
		std::cout << "Intersection: " << Intersection << std::endl;
}

void TestOctree(const std::string &Filename)
{
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(Filename.c_str());
	reader->Update();
	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
	
	Octree Tree;
	Tree.Build(polydata);
		
	OrientedPoint Intersection;
	//Ray R(vgl_point_3d<double>(0.5, 0.5, -1.0), vgl_vector_3d<double> (0.0, 0.0, 1.0));
	Ray R(vgl_point_3d<double>(10.0, 0.0, 0.0), vgl_vector_3d<double> (-1.0, 0.0, 0.0));
	bool intersect = Tree.IntersectRay(R, Intersection);
	
	std::cout << "Intersect? " << intersect << std::endl;
	if(intersect)
		std::cout << "Intersection: " << Intersection << std::endl;
}