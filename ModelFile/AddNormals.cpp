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
	vtkSmartPointer<vtkPolyData> polydata = Model.CreatePolydata();
	AddNormalsToPolyData(polydata);
	
	ModelFile NewModel(polydata);
		
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

#if 0
//try to compute normals from points - should use PINE instead
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "ModelFile.h"

#include <vgl/vgl_point_3d.h>

#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_det.h>

#include <Geometry/Geometry.h>

int main(int argc, char *argv[])
{
	Tools::AssertNumArgs(argc, 2);
	
	std::string InputFile = argv[1];
	std::string OutputFile = argv[2];
	
	ModelFile Model(InputFile);
	
	Model.ComputeNormals(.01);
	Model.Write(OutputFile);
	
	return 0;
}

#endif