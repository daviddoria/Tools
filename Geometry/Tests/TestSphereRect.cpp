#include <Geometry/Angles.h>

#include <vgl/vgl_vector_3d.h>

#include <iostream>

int main(int argc, char *argv[])
{

	//TestSphere2Rect()
	{
		std::cout << "TestSphereToRect()" << std::endl
				<< "-----------------" << std::endl;
	
		double p = 1;
		double t = .5;
	
		std::cout << "Theta: " << t << std::endl;
		std::cout << "Phi: " << p << std::endl;
	
		vgl_vector_3d<double> Direction = Sphere2Rect(p,t);

		std::cout << "Direction: " << Direction << std::endl;
	}

	//TestRect2Sphere()
	{
		std::cout << "TestRectToSphere()" << std::endl
				<< "-----------------" << std::endl;
	
		vgl_vector_3d<double> v (.259, .4034, .8775);
		double t,p;
	
		std::cout << "V: " << v << std::endl;
		Rect2Sphere(v, p, t);
	
		std::cout << "Phi: " << p << std::endl;
		std::cout << "Theta: " << t << std::endl;
	
		std::cout << "Should be p = 1, t = .5." << std::endl;

	}
	return 0;
}

