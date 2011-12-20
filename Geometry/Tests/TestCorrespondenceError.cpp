/*
#include <Geometry/Ray.h>
#include <Geometry/Color.h>
#include <Geometry/Transformation.h>
#include <Geometry/Triangle.h>
#include <Geometry/Angles.h>
#include <Geometry/Geometry.h>
*/

#include <vgl/vgl_point_3d.h>
/*
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_box_3d.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_inverse.h>

#include <VXLHelpers/VXLHelpers.h>
*/

/*
#include <vtkSmartPointer.h>
#include <vtkPlane.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
*/

#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
	//Broken!!!
	//zero distance test (sets of points are identical)
	{
		std::vector<vgl_point_3d<double> > ModelPoints;
		ModelPoints.push_back(vgl_point_3d<double> (1.0, 2.0, 3.0));
		ModelPoints.push_back(vgl_point_3d<double> (3.0, 4.0, 5.0));
		ModelPoints.push_back(vgl_point_3d<double> (6.0, 2.0, 1.0));
	
		std::vector<vgl_point_3d<double> > ScanPoints;
		ScanPoints.push_back(vgl_point_3d<double> (1.0, 2.0, 3.0));
		ScanPoints.push_back(vgl_point_3d<double> (3.0, 4.0, 5.0));
		ScanPoints.push_back(vgl_point_3d<double> (6.0, 2.0, 1.0));
	
		//double Error = CorrespondenceError(ModelPoints, ScanPoints);
		//std::cout << "Error: " << Error << std::endl;
	}

	//nonzero distance test (sets of points are different)
	{
		std::vector<vgl_point_3d<double> > ModelPoints;
		ModelPoints.push_back(vgl_point_3d<double> (1.0, 2.0, 3.0));
		ModelPoints.push_back(vgl_point_3d<double> (3.0, 4.0, 5.0));
		ModelPoints.push_back(vgl_point_3d<double> (6.0, 2.0, 1.0));
	
		std::vector<vgl_point_3d<double> > ScanPoints;
		ScanPoints.push_back(vgl_point_3d<double> (1.0, 2.0, 3.0));
		ScanPoints.push_back(vgl_point_3d<double> (3.0, 4.0, 5.0));
		ScanPoints.push_back(vgl_point_3d<double> (6.0, 3.0, 1.0)); //this point changed
	
		//double Error = CorrespondenceError(ModelPoints, ScanPoints);
		//std::cout << "Error: " << Error << std::endl;
	}
	return 0;
}

