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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//vnl_double_3x3 SampleCovarianceMatrix(const vector<vgl_point_3d<double> > &Points);
ModelFile AxisAlign(const ModelFile &Model);


void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string InputFile, OutputFile;
	double X,Y,Z;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFile), "Set input file")
			("output", po::value<std::string>(&OutputFile), "Set output file")
			("x", po::value<double>(&X)->default_value(0.0), "Center X")
			("y", po::value<double>(&Y)->default_value(0.0), "Center Y")
			("z", po::value<double>(&Z)->default_value(0.0), "Center Z")
			;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	CheckRequiredArgs(vm);
	
	ModelFile Input;
	Input.Read(InputFile);
	
	Input.SetCenterOfMass(vgl_point_3d<double>(X,Y,Z));

	Input.Write(OutputFile);
	//ModelFile AxisAligned = AxisAlign(Input);
	//WriteSuccess = AxisAligned.Write(OutputFile);
	
	//if(WriteSuccess)
	//	cout << endl << "Centered file." << endl;
	
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

/*
vnl_double_3x3 SampleCovarianceMatrix(const vector<vgl_point_3d<double> > &Points)
{
	vgl_point_3d<double> u = CenterOfMass(Points);
	
	vnl_double_3x3 CovMat;
	CovMat.fill(0);
	cout << CovMat << endl;
	
	for (unsigned int i = 0; i < Points.size(); i++)
	{
		vnl_double_3 P = vgl_point_to_vnl_vector(Points[i]);
		CovMat = CovMat + VectorToMatrix(P)*VectorToMatrix(P).transpose();
	}
	
	CovMat = (1.0/(Points.size()-1.0)) * CovMat;

	return CovMat;

}
*/

#if 0
// !!! not working
ModelFile AxisAlign(const ModelFile &Model)
{
	ModelFile AxisAligned = Model;
	
	vnl_double_3x3 CovMat = SampleCovarianceMatrix(Model.getCoords());
	//vector<double> EVals = EigenValues(CovMat);
	vector<vnl_double_3> EVecs = EigenVectors(CovMat);

	vgl_vector_3d<double> Xaxis = vnl_vector_to_vgl_vector(EVecs[0]);
	vgl_vector_3d<double> Yaxis = vnl_vector_to_vgl_vector(EVecs[1]);
	vgl_vector_3d<double> Zaxis = vnl_vector_to_vgl_vector(EVecs[2]);

	normalize(Xaxis);
	normalize(Yaxis);
	normalize(Zaxis);

	vnl_double_3x3 R;
	R.set_column(0, vgl_vector_to_vnl_vector(Xaxis));
	R.set_column(1, vgl_vector_to_vnl_vector(Yaxis));
	R.set_column(2, vgl_vector_to_vnl_vector(Zaxis));

	cout << vnl_det(R) << endl;
	
	bool IsRotation = Tools::VerifyRotationMatrix(R);
	cout << "IsRotation? " << IsRotation << endl;

	WriteAxisFile(Model.getCenterOfMass(), Xaxis, Yaxis, Zaxis, "Axis.vtp");
	
	Transformation Trans(R);
	AxisAligned.Transform(Trans);

	return AxisAligned;
}
#endif