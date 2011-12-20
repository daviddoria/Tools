#include <iostream>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_3.h>

#include <VXLHelpers/VXLHelpers.h>

bool TestEigenBig();

int main(int argc, char* argv[])
{
	bool pass = TestEigenBig();
	
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

bool TestEigenBig()
{
	vnl_vector<double> V;
	vnl_matrix<double> M(10, 10, 0.0);
	
	M(0,1) = 2.0;
	M(0,2) = 3.0;
	M(1,0) = 2.0;
	M(1,1) = 0.4;
	M(1,2) = 5.0;
	M(2,0) = 3.0;
	M(2,1) = 5.0;
	M(2,2) = 1.0;
	
	std::vector<vnl_vector<double> > EVecs = VXLHelpers::EigenVectors(M);

	return true;
}