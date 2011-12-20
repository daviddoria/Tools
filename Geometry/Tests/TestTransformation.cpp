#include <Geometry/Transformation.h>

#include <vgl/vgl_vector_3d.h>

#include <vnl/vnl_double_3x3.h>

#include <VXLHelpers/VXLHelpers.h>

#include <iostream>

int main(int argc, char *argv[])
{
	vnl_double_3x3 R = VXLHelpers::MakeRotation('z', M_PI/4.0);
	//vnl_vector<double> T;
	vgl_vector_3d<double> T(0.0, 0.0, 0.0);

	Transformation Trans(T, R);

	cout << Trans;
	Trans.WriteToFile("TestTrans.txt");

	vgl_vector_3d<double> orig(0.0, 1.0, 0.0);
	vgl_vector_3d<double> transformed = Trans.ApplyTransform(orig);

	std::cout << "orig: " << orig << std::endl;
	std::cout << "transformed: " << transformed << std::endl;
	return 0;
}

