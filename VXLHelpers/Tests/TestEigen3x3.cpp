#include <iostream>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_3.h>

#include <VXLHelpers/VXLHelpers.h>

bool TestEigen3x3();

int main(int argc, char* argv[])
{
	bool pass = TestEigen3x3();
	
	if(pass)
	{
		std::cout << "Success." << std::endl;
		return 0;
	}
	else
	{
		std::cout << "Fail!" << std::endl;
		return -1;
	}
}

bool TestEigen3x3()
{
	//!!! check this - wrong matrix!!!
	vnl_double_3x3 M;
	M(0,0) = 2.0;
	M(0,1) = 2.0;
	M(0,2) = 3.0;
	M(1,0) = 2.0;
	M(1,1) = 0.4;
	M(1,2) = 5.0;
	M(2,0) = 3.0;
	M(2,1) = 5.0;
	M(2,2) = 1.0;
	
	//vectors
	std::vector<vnl_double_3> Vecs = VXLHelpers::EigenVectors(M);
	for(unsigned int i = 0; i < Vecs.size(); i++)
		std::cout << Vecs[i] << std::endl;
		
	//eig vecs (smallest to largest)
	//(0.1433, 0.6850, -0.7143)
	//(.8686, -.4330, -.2410)
	//(.4744, .5859, .6570)
	double eps = 1e-3;
	
	if(!(VXLHelpers::CloseEnough(Vecs[0], vnl_double_3(0.1433, 0.6850, -0.7143), eps) || VXLHelpers::CloseEnough(Vecs[0], vnl_double_3(-0.1433, -0.6850, 0.7143), eps)))
	{
		std::cout << "First evec is " << Vecs[0] << " should be " << vnl_double_3(0.1433, 0.6850, -0.7143) << std::endl;
		return false;
	}
	
	if(!(VXLHelpers::CloseEnough(Vecs[1], vnl_double_3(.8686, -.4330, -.2410), eps) || VXLHelpers::CloseEnough(Vecs[1], vnl_double_3(-.8686, .4330, .2410), eps)))
	{
		std::cout << "Second evec is " << Vecs[1] << " should be " << vnl_double_3(.8686, -.4330, -.2410) << std::endl;
		return false;
	}
	
	if(!(VXLHelpers::CloseEnough(Vecs[2], vnl_double_3(.4744, .5859, .6570), eps) || VXLHelpers::CloseEnough(Vecs[2], vnl_double_3(-.4744, -.5859, -.6570), eps)))
	{
		std::cout << "Third evec is " << Vecs[2] << " should be " << vnl_double_3(.4744, .5859, .6570) << std::endl;
		return false;
	}
	
	
	//values
	std::vector<double> Vals = VXLHelpers::EigenValues(M);
	for(unsigned int i = 0; i < Vals.size(); i++)
		std::cout << Vals[i] << std::endl;
	
	//eig vals (smallest to largest) = -4.3961, -0.8293, 7.6254
	
	if(!VXLHelpers::CloseEnough(Vals[0], -4.3961, eps))
	{
		std::cout << "First eval is " << Vals[0] << " should be " << -4.3961 << std::endl;
		return false;
	}
	
	if(!VXLHelpers::CloseEnough(Vals[1], -0.8293, eps))
	{
		std::cout << "First eval is " << Vals[1] << " should be " << -0.8293 << std::endl;
		return false;
	}
	
	if(!VXLHelpers::CloseEnough(Vals[2], 7.6254, eps))
	{
		std::cout << "First eval is " << Vals[2] << " should be " << 7.6254 << std::endl;
		return false;
	}
	
	//everything went well
	return true;
}