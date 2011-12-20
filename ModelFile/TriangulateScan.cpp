#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "ModelFile.h"
#include "Triangulate.h"

#include <Tools/Tools.h>
#include <Geometry/Color.h>


#include <boost/program_options.hpp>
namespace po = boost::program_options;

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	vgl_point_3d<double> ScannerLocation;
	std::string InputFilename, OutputFilename;
	double MaxLength;

	double ScannerX, ScannerY, ScannerZ;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFilename), "Set input file")
			("output", po::value<std::string>(&OutputFilename), "Set output file")
			("maxlength", po::value<double>(&MaxLength)->default_value(.1), "Maximum edge length to keep.")
			("ScannerX", po::value<double>(&ScannerX), "ScannerX.")
			("ScannerY", po::value<double>(&ScannerY), "ScannerY.")
			("ScannerZ", po::value<double>(&ScannerZ), "ScannerZ.")
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

	std::cout << "Max edge length: " << MaxLength << std::endl;
	
	ModelFile Model;
	Model.Read(InputFilename);
	Model.Init();
	assert(Model.IsValid());
	
	bool IsScan = Model.getScannerLocation(ScannerLocation);
	if(!IsScan)
	{
		if(!vm.count("ScannerX") || !vm.count("ScannerY") || !vm.count("ScannerZ"))
		{
			std::cout << "Model does not have scanner info and scanner location was not input!" << std::endl;
		}

		ScannerLocation = vgl_point_3d<double> (ScannerX, ScannerY, ScannerZ);
	}

	std::vector<Color<unsigned char> > Colors = Model.getColors();
	
	std::vector<vector<unsigned int> > VertexList;
	std::vector<vgl_point_3d<double> > Points = Model.getCoords();
	
	TriangulateLidar(ScannerLocation, Points, VertexList, MaxLength);

	ModelFile TriangulatedModel;
	TriangulatedModel.setCoords(Points);
	TriangulatedModel.setVertexLists(VertexList);
	TriangulatedModel.setColors(Colors);	
	
	std::cout << "Num tris: " << TriangulatedModel.NumTriangles() << std::endl;
	
	vtkSmartPointer<vtkPolyData> PD = TriangulatedModel.CreatePolydata();
	AddNormalsToPolyData(PD);
	
	ModelFile ModelWithNormals;
	ModelWithNormals.ConstructFromPolydata(PD);
	ModelWithNormals.Write(OutputFilename);
	
	std::cout << "Finished Triangulating." << std::endl;
	
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