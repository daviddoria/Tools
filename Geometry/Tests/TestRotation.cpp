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
	vnl_double_3x3 A = VXLHelpers::MakeRotation('x', deg2rad(45.0));
	vgl_vector_3d<double> T(0.0,0.0,0.0);

	Transformation Trans(T,A);

	vgl_vector_3d<double> V(0.0,1.0,0.0);
	vgl_vector_3d<double> Vout = Trans.ApplyTransform(V);

	std::cout << "V = " << V << std::endl;
	std::cout << "Vout = " << Vout << std::endl;
	
	return 0;
}

