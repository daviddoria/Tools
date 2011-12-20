#include <iostream>
#include <string>

#include <vnl/vnl_double_3x3.h>

#include "ModelFile.h"

#include <Geometry/Transformation.h>
#include <Tools/Tools.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

void TransformFromFile(const std::string &ModelFilename, const std::string &TransformFile, const std::string &OutputFile);

void CheckRequiredArgs(const po::variables_map &vm);

int main(int argc, char *argv[])
{
	std::string ModelFilename, OutputFilename, TransformFilename;
	
	double tx,ty,tz;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "produce help message")
			("model", po::value<std::string>(&ModelFilename), "Set model file.")
			("xform", po::value<std::string>(&TransformFilename), "Set transform file.")
			("tx", po::value<double>(&tx), "Set translation x.")
			("ty", po::value<double>(&ty), "Set translation y.")
			("tz", po::value<double>(&tz), "Set translation z.")
			("output", po::value<std::string>(&OutputFilename), "Set output file.")
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
	

	
	if(vm.count("xform"))
	{
		if(!Tools::FileExists(TransformFilename))
		{
			std::cout << "Transform file does not exist!" << std::endl;
			exit(-1);
		}
		
		TransformFromFile(ModelFilename, TransformFilename, OutputFilename);
	}
	else
	{
		if(vm.count("tx") && vm.count("ty") && vm.count("tz"))
		{
			ModelFile OrigModel;
			OrigModel.Read(ModelFilename);
			OrigModel.Init();
			vgl_vector_3d<double> T(tx, ty, tz);
			Transformation Trans(T);
			std::cout << "Translation: " << Trans.getTranslation() << std::endl;
	
			ModelFile Model = OrigModel;

			Model.Transform(Trans);
	
			Model.Write(OutputFilename);
		}
	}
	
	

	return 0;
}

void TransformFromFile(const std::string &ModelFilename, const std::string &TransformFile, const std::string &OutputFile)
{
	ModelFile OrigModel;
	OrigModel.Read(ModelFilename);
	OrigModel.Init();
	
	//read transformation
	Transformation Trans;
	Trans.ReadFromFile(TransformFile);
	std::cout << "Rotation: " << Trans.getRotation() << std::endl;
	std::cout << "Translation: " << Trans.getTranslation() << std::endl;
	
	ModelFile Model = OrigModel;

	Model.Transform(Trans);
	
	Model.Write(OutputFile);
}


void CheckRequiredArgs(const po::variables_map &vm)
{
	if(!vm.count("model"))
	{
		std::cout << "Input model was not set." << std::endl;
		exit(-1);
	}
	/*
	if(!vm.count("xform"))
	{
		std::cout << "Transform file was not set." << std::endl;
		exit(-1);
	}
	*/
	if(!vm.count("output"))
	{
		std::cout << "Output was not set." << std::endl;
		exit(-1);
	}
	
	
}