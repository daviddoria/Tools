#include <iostream>
#include <string>
#include <stdlib.h>

#include "ModelFile.h"

#include <Tools/Tools.h>

int main(int argc, char *argv[])
{
	//check input arguments
	AssertNumArgs(argc,2);

	//parse input arguments
	std::string InputFilename = argv[1];
	std::string OutputFilename = argv[2];

	ModelFile Model(InputFilename);
	
	ModelFile Points(Model.getPoints());
	
	Points.Write(OutputFilename);

	return 0;
}
