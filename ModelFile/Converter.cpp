#include <iostream>
#include <string>

#include "ModelFile.h"
#include <Tools/Tools.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string InputFile, OutputFile;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFile), "Set input file")
			("output", po::value<std::string>(&OutputFile), "Set output file")
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

	std::cout << "Converting " << InputFile << " to " << OutputFile << "..." << std::endl;
	
	ModelFile Model;
	bool ReadSuccess = Model.Read(InputFile);
	if(!ReadSuccess)
	{
		std::cout << "Reading " << InputFile << " failed!" << std::endl;
		exit(-1);
	}
	
	bool WriteSuccess = Model.Write(OutputFile);
	if(!WriteSuccess)
	{
		std::cout << "Writing " << InputFile << " failed!" << std::endl;
		exit(-1);
	}
	
	std::cout << "Converted " << InputFile << " to " << OutputFile << std::endl;
	
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