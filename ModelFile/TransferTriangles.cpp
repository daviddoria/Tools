#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "ModelFile.h"

#include <Tools/Tools.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string WithTrianglesFilename, WithoutTrianglesFilename, OutputFilename;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("withtris", po::value<std::string>(&WithTrianglesFilename), "Set file with triangles.")
			("withouttris", po::value<std::string>(&WithoutTrianglesFilename), "Set file without triangles.")
			("output", po::value<std::string>(&OutputFilename), "Set output file.")
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

	ModelFile ModelWithTriangles;
	ModelWithTriangles.Read(WithTrianglesFilename);
	ModelWithTriangles.Init();
	
	ModelFile ModelWithoutTriangles;
	ModelWithoutTriangles.Read(WithoutTrianglesFilename);
	ModelWithoutTriangles.Init();
	
	if(ModelWithTriangles.NumPoints() != ModelWithoutTriangles.NumPoints())
	{
		std::cout << "The two models must have identical point lists!" << std::endl;
		exit(-1);
	}
	
	std::vector<vector<unsigned int> > VertexLists = ModelWithTriangles.getVertexLists();
	ModelWithoutTriangles.setVertexLists(VertexLists);
	
	ModelWithoutTriangles.Write(OutputFilename);
	
	return 0;
}

void CheckRequiredArgs(const po::variables_map vm)
{
	if(!vm.count("withtris"))
	{
		std::cout << "File with triangles is required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("withouttris"))
	{
		std::cout << "File without triangles is required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("output"))
	{
		std::cout << "Output is required!" << std::endl;
		exit(-1);
	}
}