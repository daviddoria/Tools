#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "ModelFile.h"

#include <Tools/Tools.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string InputFilename, OutputFilename;
	
	double x,y,z;
	bool verify;
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFilename), "Set input file.")
			("output", po::value<std::string>(&OutputFilename), "Set output file.")
			("verify", po::value<bool>(&verify)->default_value(false), "Verify location?")
			("x", po::value<double>(&x), "Set x")
			("y", po::value<double>(&y), "Set y")
			("z", po::value<double>(&z), "Set z")
			;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		exit(-1);
	}
	
	CheckRequiredArgs(vm);
	
	ModelFile Model;
	Model.Read(InputFilename);
	Model.Init();
	Model.setScannerLocation(vgl_point_3d<double> (x,y,z));
	Model.Write(OutputFilename);
	
	if(verify)
	{
		//verify that the location was added correctly
		ModelFile NewModel;
		NewModel.Read(OutputFilename);
		NewModel.Init();
		vgl_point_3d<double> Loc;
		bool IsScan = NewModel.getScannerLocation(Loc);
		std::cout << "The new file has tested " << IsScan << " for having a scanner location and it is"  << Loc << std::endl;
	}
	
	return 0;
}

void CheckRequiredArgs(const po::variables_map vm)
{
	if(!vm.count("input"))
	{
		std::cout << "input is required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("output"))
	{
		std::cout << "output is required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("x") || !vm.count("y") || !vm.count("z") )
	{
		std::cout << "Scanner location (x,y,z) is required!" << std::endl;
		exit(-1);
	}
}
