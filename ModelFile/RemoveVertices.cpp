#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "ModelFile.h"

#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include <Tools/Tools.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string InputFilename, OutputFilename;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFilename), "Set input file")
			("output", po::value<std::string>(&OutputFilename), "Set output file")
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
	assert(Model.IsValid());
	
	//keep everything except vertices
	vtkSmartPointer<vtkPolyData> polydata = Model.CreatePolydata();
	polydata->GetVerts()->Reset();
	
	ModelFile NewModel;
	NewModel.ConstructFromPolydata(polydata);
	NewModel.Write(OutputFilename);
	
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
}