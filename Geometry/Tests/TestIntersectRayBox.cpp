#include <Geometry/Ray.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_box_3d.h>

//#include <VXLHelpers/VXLHelpers.h>

#include <iostream>

int main(int argc, char *argv[])
{
	//green
	{
		std::cout << "TestIntersectRayBoxGreen()" << std::endl
				<< "---------------------" << std::endl;
		Ray R(vgl_point_3d<double> (0,0,10), vgl_vector_3d<double> (0, 0, -1));
		vgl_box_3d<double> Box(vgl_point_3d<double>(-1, -1, -1), vgl_point_3d<double> (1,1,1));
	
		bool bIntersect = R.IntersectBox(Box);
	
		std::cout << "R: " << R << std::endl;
		std::cout << "Intersect? " << bIntersect << std::endl;
	
	}

	//red
	{
		std::cout << "TestIntersectRayBoxRed()" << std::endl
				<< "---------------------" << std::endl;
		Ray R(vgl_point_3d<double> (10,0,10), vgl_vector_3d<double> (0, 0, -1));
		vgl_box_3d<double> Box(vgl_point_3d<double>(-1, -1, -1), vgl_point_3d<double> (1,1,1));
	
		bool bIntersect = R.IntersectBox(Box);
	
		std::cout << "R: " << R << std::endl;
		std::cout << "Intersect? " << bIntersect << std::endl;
	
	}
	return 0;
}

