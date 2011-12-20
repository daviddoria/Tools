#include "Ray.h"
#include "Color.h"
#include "Transformation.h"
#include "Triangle.h"
#include "Angles.h"
#include "Geometry.h"

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


	return 0;
}





void TestIntersectRayPlane()
{
	Ray R(vgl_point_3d<double> (0,0,10), vgl_vector_3d<double> (0, 0, -1));	
	vgl_plane_3d<double> P(vgl_vector_3d<double> (0, 0, 1), vgl_point_3d<double> (0,0,0));
	
	cout << "R: " << R << endl;
	cout << "P: " << P << endl;
	vgl_point_3d<double> Intersection;
	bool bIntersect = R.IntersectPlane(P, Intersection);
	
	cout << "Intersect? " << bIntersect << endl;
}


void TestTransformationReadWrite(void)
{
	Transformation Trans;
	Trans.ReadFromFile("TestTrans.txt");
	std::cout << Trans << std::endl;
	Trans.WriteToFile("TestTrans2.txt");
}


void TestWriteHorizontalColorBar()
{
	WriteHorizontalColorBar(5, 100);
}

void TestWriteVerticalColorBar()
{
	//WriteVerticalColorBar(5, 100);
	WriteVerticalColorBar(50, 1000);
}
