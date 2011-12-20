#include <iostream>

#include "VXLHelpers.h"
#include <Geometry/Angles.h>

void TestAdd();
void TestIntersectBox();

void TestRotateVector();
void TestRotateVectorWithMatrix();
void TestRotationMatrixCreation();

void TestVectorProjection();
void TestGet4x4();

int main()
{
	//TestAdd();
	
	TestRotationMatrixCreation();
	//TestRotateVector();
	//TestRotateVectorWithMatrix();
		
	//TestVectorProjection();
	//TestGet4x4();
	return 0;
}


void TestAdd()
{
	vgl_point_3d<double> P1(1.0, 2.0, 3.0);
	vgl_point_3d<double> P2(4.0, 5.0, 6.0);
	
	vgl_point_3d<double> Sum = P1 + VXLHelpers::vgl_point_to_vgl_vector(P2);
	
	std::cout << Sum << std::endl;
}

void TestRotationMatrixCreation()
{
	vnl_double_3 axis(0.0, 1.0, 0.0);
	vnl_double_3x3 R = VXLHelpers::RotationMatrix(axis, -M_PI/2.0);
	
	std::cout << "R: " << R << std::endl;
	

}

void TestRotateVectorWithMatrix()
{
	vnl_double_3 axis(0.0, 1.0, 0.0);
	vnl_double_3x3 R = VXLHelpers::RotationMatrix(axis, -M_PI/2.0);
	
	vnl_double_3 V(1.0, 0.0, 0.0);
	
	vnl_double_3 Rotated = R*V;
	
	std::cout << "V: " << V << std::endl;
	std::cout << "Rotated: " << Rotated << std::endl;
	
	vnl_double_3 Rotated2 = VXLHelpers::Rotate(V, axis, -M_PI/2.0);
	
	std::cout << "Rotated2: " << Rotated2 << std::endl;
	
	vgl_vector_3d<double> v(0.0, 1.0, 0.0);
}

void TestRotateVector()
{
	vnl_double_3 axis(0.0, 1.0, 0.0);
	vnl_double_3 V(1.0, 0.0, 0.0);
	vnl_double_3 Rotated = VXLHelpers::Rotate(V, axis, -M_PI/2.0);
	
	std::cout << "Rotated: " << Rotated << std::endl;
	
}

void TestVectorProjection()
{
	vgl_vector_3d<double> A(2.0, 1.0, 0.0);
	vgl_vector_3d<double> B(1.0, 0.0, 0.0);
	vgl_vector_3d<double> proj = VXLHelpers::ProjectAonB(A,B);

	std::cout << "A: " << A << std::endl;
	std::cout << "B: " << B << std::endl;
	std::cout << "projAonB: " << proj << std::endl;
}

void TestGet4x4()
{
	vnl_double_3x3 M(4.5);
	std::cout << M << std::endl;
	/*
	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigne)
	}
	*/
	vnl_matrix<double> M4 = VXLHelpers::Get4x4(M);
	std::cout << M4 << std::endl;
}
