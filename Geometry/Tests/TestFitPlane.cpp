#include <Geometry/Ray.h>
#include <Geometry/Color.h>
#include <Geometry/Transformation.h>
#include <Geometry/Triangle.h>
#include <Geometry/Angles.h>
#include <Geometry/Geometry.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_box_3d.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_inverse.h>

#include <VXLHelpers/VXLHelpers.h>

#include <vtkSmartPointer.h>
#include <vtkPlane.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>

#include <iostream>

int main(int argc, char *argv[])
{
	std::vector<vgl_point_3d<double> > Points;
	Points.push_back(vgl_point_3d<double> (0.0,0.0,0.0));
	Points.push_back(vgl_point_3d<double> (.3,.3,.3));
	Points.push_back(vgl_point_3d<double> (1.0,1.1,1.2));
	Points.push_back(vgl_point_3d<double> (1.0,1.1,1.4));

	vgl_plane_3d<double> Plane = geom::FitPlane(Points);

	std::cout << Plane << std::endl;

	vgl_vector_3d<double> Normal = Plane.normal();

	vgl_point_3d<double> P0 = vgl_point_3d<double>(0.0,0.0,0.0) + Normal * Plane.d();
	
	vtkSmartPointer<vtkPlane> VTKPlane = vtkSmartPointer<vtkPlane>::New();
	VTKPlane->SetNormal(Normal.x(), Normal.y(), Normal.z());
	double origin[3] = {P0.x(), P0.y(), P0.z()};
	VTKPlane->SetOrigin(origin);

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIdList> PointIDs = vtkSmartPointer<vtkIdList>::New();
	

	double p0[3] = {0.0, 0.0, 0.0};
	double p1[3] = {0.3, 0.3, 0.3};
	double p2[3] = {1.0, 1.1, 1.2};
	double p3[3] = {1.0, 1.1, 1.4};
	
	pts->InsertNextPoint(p0);
	PointIDs->InsertNextId(0);
	
	pts->InsertNextPoint(p1);
	PointIDs->InsertNextId(1);
	
	pts->InsertNextPoint(p2);
	PointIDs->InsertNextId(2);

	pts->InsertNextPoint(p3);
	PointIDs->InsertNextId(3);

	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();
	Vertices->InsertNextCell(PointIDs);

	
	vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();

	//add the points to the dataset
	pdata->SetPoints(pts);
	pdata->SetVerts(Vertices);

	//write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInput(pdata);
	writer->SetFileName("TestPlane.vtp");
	writer->Write();
	return 0;
}

