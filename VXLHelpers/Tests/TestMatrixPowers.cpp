#include <iostream>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <VXLHelpers/VXLHelpers.h>

bool TestMatrixPowers();

int main(int argc, char* argv[])
{
	bool pass = TestMatrixPowers();
	
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

bool TestMatrixPowers()
{
	vnl_vector<double> V;
	vnl_matrix<double> M(3, 3);
	
	M(0,0) = 1.0;
	M(0,1) = 2.0;
	M(0,2) = 3.0;
	M(1,0) = 2.0;
	M(1,1) = 4.0;
	M(1,2) = 3.0;
	M(2,0) = 3.0;
	M(2,1) = 3.0;
	M(2,2) = 1.0;
	
	//should be
	vnl_matrix<double> ShouldBe(3, 3);
	
	ShouldBe(0,0) = 14.0;
	ShouldBe(0,1) = 19.0;
	ShouldBe(0,2) = 12.0;
	ShouldBe(1,0) = 19.0;
	ShouldBe(1,1) = 29.0;
	ShouldBe(1,2) = 21.0;
	ShouldBe(2,0) = 12.0;
	ShouldBe(2,1) = 21.0;
	ShouldBe(2,2) = 19.0;
	
	vnl_matrix<double> Raised = VXLHelpers::MatrixPower(M,2);
	
	std::cout << "M: " << std::endl << M << std::endl;
	std::cout << "ShouldBe: " << std::endl << ShouldBe << std::endl;
	std::cout << "Raised: " << std::endl << Raised << std::endl;
	
	if(!VXLHelpers::CloseEnough(ShouldBe, Raised, 1e-3))
		return false;
	
	return true;
}