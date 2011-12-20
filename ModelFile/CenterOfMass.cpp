#include <iostream>
#include <string>

#include <vnl/vnl_double_3x3.h>

#include "ModelFile.h"

#include <Geometry/Transformation.h>
#include <Tools/Tools.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map &vm);

int main(int argc, char *argv[])
{
	std::string ModelFilename, OutputFilename;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "produce help message")
			("model", po::value<std::string>(&ModelFilename), "Set model file.")
			;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if(vm.count("help"))
	{
		std::cout << desc << "\n";
		exit(-1);
	}
	CheckRequiredArgs(vm);
	
	if(!Tools::FileExists(ModelFilename))
	{
		std::cout << "Model file does not exist!" << std::endl;
		exit(-1);
	}
	
	ModelFile Model;
	Model.Read(ModelFilename);
	Model.Init();
	std::cout << "Center of mass: " << Model.getCenterOfMass() << std::endl;


	return 0;
}

void CheckRequiredArgs(const po::variables_map &vm)
{
	if(!vm.count("model"))
	{
		std::cout << "Input model was not set." << std::endl;
		exit(-1);
	}
	
}