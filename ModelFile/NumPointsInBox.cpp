#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "ModelFile.h"
#include "Triangulate.h"

#include <Tools/Tools.h>
#include <Geometry/Color.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_box_3d.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	
	std::string InputFilename;
	double SideLength;

	double CenterX, CenterY, CenterZ;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFilename), "Set input file")
			("sidelength", po::value<double>(&SideLength), "Edge length of box.")
			("CenterX", po::value<double>(&CenterX), "CenterX.")
			("CenterY", po::value<double>(&CenterY), "CenterY.")
			("CenterZ", po::value<double>(&CenterZ), "CenterZ.")
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

	std::cout << "Side length: " << SideLength << std::endl;
	
	ModelFile Model;
	Model.Read(InputFilename);
	Model.Init();
	assert(Model.IsValid());
	
	std::vector<vgl_point_3d<double> > Points = Model.getCoords();
	
	vgl_point_3d<double> Center(CenterX, CenterY, CenterZ);
	vgl_box_3d<double> Box(Center, SideLength, SideLength, SideLength, vgl_box_3d<double>::centre);
	
	unsigned int NumPoints = 0;
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		if(Box.contains(Points[i]))
			NumPoints++;
	}
	
	std::cout << "There are " << NumPoints << " in the box." << std::endl;
	
	std::string Filename = "values.txt";
	std::ofstream fout(Filename.c_str(), std::ios::app);
	fout << NumPoints << std::endl;
	fout.close();
	
	return 0;
}

void CheckRequiredArgs(const po::variables_map vm)
{
	if(!vm.count("input"))
	{
		std::cout << "input is required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("sidelength"))
	{
		std::cout << "Side length is required!" << std::endl;
		exit(-1);
	}
	
	if(!vm.count("CenterX"))
	{
		std::cout << "CenterX is required!" << std::endl;
		exit(-1);
	}
	if(!vm.count("CenterY"))
	{
		std::cout << "CenterY is required!" << std::endl;
		exit(-1);
	}
	if(!vm.count("CenterZ"))
	{
		std::cout << "CenterZ is required!" << std::endl;
		exit(-1);
	}

	
}