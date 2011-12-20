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

void TestAddArray();
		
int main()
{
	TestAddArray();
	
	return 0;
}

void TestAddArray()
{
	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();
	
	for(unsigned int i = 0; i < 3; i++)
	{
		vtkIdType pid[1];
		pid[0] = Points->InsertNextPoint(drand48(), drand48(), drand48());
		Vertices->InsertNextCell(1,pid);
		
	}

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(Points);
	polydata->SetVerts(Vertices);
	
	std::vector<double> V;
	V.push_back(10.0);
	V.push_back(12.0);
	V.push_back(14.0);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInput(polydata);
	writer->SetFileName("WithoutArray.vtp");
	writer->Write();
	
	AddArrayToPolydata<double, vtkDoubleArray>(V, "Test", polydata);
	
	{
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetInput(polydata);
		writer->SetFileName("WithArray.vtp");
		writer->Write();
	}
	
}
