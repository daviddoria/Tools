#include <Geometry/Ray.h>
//#include <Geometry/Color.h>
//#include <Geometry/Transformation.h>
//#include <Geometry/Triangle.h>
//#include <Geometry/Angles.h>
//#include <Geometry/Geometry.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
//#include <vgl/vgl_box_3d.h>

//#include <vnl/vnl_matrix.h>
//#include <vnl/vnl_double_3x3.h>
//#include <vnl/vnl_inverse.h>

//#include <VXLHelpers/VXLHelpers.h>

//#include <GL/glut.h>

//#include <SDLHelpers/SDLHelpers.h>

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
	//green
	{
		Ray R0(vgl_point_3d<double> (0.0, 0.0, 0.0), vgl_vector_3d<double> (1.0, 1.0, 1.0));
		Ray R1(vgl_point_3d<double> (0.0, 0.0, 0.0), vgl_vector_3d<double> (1.0, 1.0, 1.0));
		cout << R0.IsSameDirection(R1) << " (should be 1)" << endl;
	}
	
	//red
	{
		Ray R0(vgl_point_3d<double> (0.0, 0.0, 0.0), vgl_vector_3d<double> (1.0, 1.0, 1.0));
		Ray R1(vgl_point_3d<double> (0.0, 0.0, 0.0), vgl_vector_3d<double> (-1.0, -1.0, -1.0));
		cout << R0.IsSameDirection(R1) << " (should be 0)" << endl;
	}
	return 0;
}

