
#include <vgl/vgl_vector_3d.h>
#include <iostream>

int main(int argc, char *argv[])
{
	//!!!Broken
	vgl_vector_3d<double> perpendicular(0, 1, 0);
	
	vgl_vector_3d<double> up(0,0,1);
	
	vgl_vector_3d<double> forward(1.0,0.0,0.0);
	
	//vgl_vector_3d<double> v(1.0,0.0,0.0);
	
	//vgl_vector_3d<double> Rotated = AxisAngle(forward, up, M_PI/2.0); //should be = perpendicular
	
	//cout << "Rotated = " << Rotated << endl;
	return 0;
}

