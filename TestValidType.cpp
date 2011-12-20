#include <iostream>
#include "ValidType.h"

int main()
{
	{
		cout << "Invalid Test:" << endl;
		ValidType<double> Test;
		cout << Test.Valid << endl;
	}
	
	{
		cout << "Valid Test:" << endl;
		ValidType<double> Test(4.5);
		cout << Test.Valid << endl;
		cout << Test.Value << endl;
	}
	
	return 0;
}