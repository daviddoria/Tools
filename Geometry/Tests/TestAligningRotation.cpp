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
	vgl_vector_3d<double> a(-1.0, -1.0, 1.0);
	vgl_vector_3d<double> b(1.0, -1.0, 1.0);
	vgl_vector_3d<double> c = cross_product(a, b);
	
	/*
	vnl_double_3x3 R = FindAligningRotation(a,b,c);
	
	std::cout << "a: " << a << std::endl;
	std::cout << "b: " << b << std::endl;
	std::cout << "c: " << c << std::endl;
	
	std::cout << "R: " << R << std::endl;
	
	vgl_vector_3d<double> newa = VXLHelpers::vnl_vector_to_vgl_vector(vnl_inverse(R)*VXLHelpers::vgl_vector_to_vnl_vector(a));
	vgl_vector_3d<double> newb = VXLHelpers::vnl_vector_to_vgl_vector(vnl_inverse(R)*VXLHelpers::vgl_vector_to_vnl_vector(b));
	vgl_vector_3d<double> newc = VXLHelpers::vnl_vector_to_vgl_vector(vnl_inverse(R)*VXLHelpers::vgl_vector_to_vnl_vector(c));
	
	std::cout << "newa: " << newa << std::endl;
	std::cout << "newb: " << newb << std::endl;
	std::cout << "newc: " << newc << std::endl;
	*/
	return 0;
}

