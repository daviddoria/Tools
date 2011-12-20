#include <iostream>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <VXLHelpers/VXLHelpers.h>

int main(int argc, char* argv[])
{
	
	vnl_matrix<double> M(3,3);
	M(0,0) = 1.0;
	M(0,1) = 2.0;
	M(0,2) = 3.0;
	M(1,0) = 4.0;
	M(1,1) = 5.0;
	M(1,2) = 6.0;
	M(2,0) = 7.0;
	M(2,1) = 8.0;
	M(2,2) = 9.0;
	
	vnl_vector<double> V(9);
	V(0) = 1.0;
	V(1) = 2.0;
	V(2) = 3.0;
	V(3) = 4.0;
	V(4) = 5.0;
	V(5) = 6.0;
	V(6) = 7.0;
	V(7) = 8.0;
	V(8) = 9.0;
	
	vnl_vector<double> Vectorized = VXLHelpers::Vectorize(M);
	
	if(VXLHelpers::CloseEnough(V, Vectorized, 1e-3))
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
