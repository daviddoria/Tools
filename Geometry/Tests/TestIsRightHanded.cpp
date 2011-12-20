#include <Geometry/Angles.h>
#include <vgl/vgl_vector_3d.h>

#include <VXLHelpers/VXLHelpers.h>

#include <iostream>

int main(int argc, char *argv[])
{
	//green
	vgl_vector_3d<double> x(1.0,0.0,0.0);
	vgl_vector_3d<double> y(0.0,1.0,0.0);
	vgl_vector_3d<double> z(0.0,0.0,1.0);
	std::cout << IsRightHanded(x,y,z) << std::endl;

	z = vgl_vector_3d<double>(0.0,0.0,-1.0);
	std::cout << IsRightHanded(x,y,z) << std::endl;

	return 0;
}

