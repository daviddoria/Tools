#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "ModelFile.h"

#include <Tools/Tools.h>

#include <vtkAppendPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void Append(ModelFile &Combined, const ModelFile &Single);
void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string OutputFilename;
	std::vector<std::string> InputFiles;
	
	//--output must come BEFORE --inputs because the multitoken list breaks anything after it
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("output", po::value<std::string>(&OutputFilename), "Set output file")
			("inputs", po::value<std::vector<std::string> >(&InputFiles)->multitoken(), "Set input files")
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
	
	if(InputFiles.size() < 2)
	{
		std::cout << "Cannot combine less than 2 files!" << std::endl;
		exit(-1);
	}
	
	vtkSmartPointer<vtkAppendPolyData> CombinedData = vtkSmartPointer<vtkAppendPolyData>::New();;
	
	//add the remaining models to the first one
	for(unsigned int i = 0; i < InputFiles.size(); i++)
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		std::clog << "Reading " << InputFiles[i] << endl;
		reader->SetFileName(InputFiles[i].c_str());
		reader->Update();
		vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
		
		CombinedData->AddInput(polydata);
	}

	//save the result
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInput(CombinedData->GetOutput());
	writer->SetFileName(OutputFilename.c_str());
	writer->Write();

	std::cout << "Finished combining." << std::endl;
	
	return 0;
}


void CheckRequiredArgs(const po::variables_map vm)
{
	if(!vm.count("inputs"))
	{
		std::cout << "input files are required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("output"))
	{
		std::cout << "output is required!" << std::endl;
		exit(-1);
	}
}
