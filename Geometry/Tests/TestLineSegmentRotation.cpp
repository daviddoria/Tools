#include <Geometry/Ray.h>
#include <Geometry/Geometry.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

#include <vnl/vnl_matrix.h>

#include <iostream>

int main(int argc, char *argv[])
{
	{
		Ray R1(vgl_point_3d<double>(0.0,0.0,0.0), vgl_vector_3d<double>(0.0,1.0,0.0));
		Ray R2(vgl_point_3d<double>(1.0,0.0,0.0), vgl_vector_3d<double>(0.0,1.0,1.0));
		vnl_matrix<double> M = geom::CameraTransform(R1, R2);
		std::cout << M << std::endl;
	}

	{
		Ray R1(vgl_point_3d<double>(0.0,0.0,0.0), vgl_vector_3d<double>(0.0,1.0,0.0));
		Ray R2(vgl_point_3d<double>(0.0,0.0,0.0), vgl_vector_3d<double>(0.0,1.0,0.0));
		vnl_matrix<double> M = geom::CameraTransform(R1, R2);
		std::cout << M << std::endl;
	}
	return 0;
}

