#include <iostream>
#include <string>
#include <vector>

#include "../ModelFile.h"

#include <vgl/vgl_point_3d.h>

#include <Geometry/Color.h>
//#include <Tools/Tools.h>

void TestWriteTriSquare(const std::string &OutputFile);
void TestWritePointSquare(const std::string &OutputFile);

void TestReadWrite(const std::string &InputFile, const string &OutputFile);

void TestColoredPointSquare(const std::string &OutputFile);
void TestColoredTriSquare(const std::string &OutputFile);

void TestCube();
void TestTransformModel(ModelFile &Model, const std::string &OutputFile);

void TestRemoveGround(ModelFile &Model, const std::string &OutputFile);
void TestComputeNormals(ModelFile &Model, const std::string &OutputFile);
void TestProjectOntoXY(ModelFile &Model, const std::string &OutputFile);
void TestSetPointIndices(ModelFile &Model, const std::string &OutputFile);

void TestDownsample(ModelFile &Model, const std::string &OutputFile);
void TestSnapToPoints(ModelFile &Messy, ModelFile &Clean, const std::string &OutputFile);

void TestReadScannerLocation(const ModelFile &Model);

void TestWriteMemoryUsage(const ModelFile &Model);

void TestDeleteTriangles(ModelFile &Model, const std::string &OutputFile);

void TestCopy(const ModelFile &Model, const std::string &OutputFile);

int main(int argc, char* argv[])
{
	std::string InputFile = argv[1];
	string OutputFile = argv[2];

	ModelFile Model;
	Model.Read(InputFile);
	Model.Init();
	
	//Model.Write("AudiA4Medium.txt");
	
	//ModelFile Messy(InputFile);
	//ModelFile Clean(OutputFile);
	//TestSnapToPoints(Messy,Clean,"AudiA4fineSnapped.vtp");
	//Messy.Write("Test.vtp");
	
	//TestWritePointSquare("TestPointSquare.vtp");
	//TestWritePointSquare("TestPointSquare.vtu");
	
	//TestReadWrite("TestPointSquare.vtp", "TestPointSquare2.vtp");
	//TestReadWrite("TestPointSquare.vtu", "TestPointSquare2.vtu");
	
	//TestWriteTriSquare("TestTriSquare.vtp");
	//TestWriteTriSquare("TestTriSquare.vtu");
	
	//TestReadWrite("TestTriSquare.vtp", "TestTriSquare2.vtp");
	//TestReadWrite("TestTriSquare.vtu", "TestTriSquare2.vtu");
	
	//TestReadWrite("Neon.vtp", "NeonOut.vtp");
	//TestReadWrite("parkinglot.vtp", "parkinglotOut.vtp");
	//TestRemoveGround(Model, "parkinglot_noground.vtp");
	//TestProjectOntoXY(Model, "projected.vtp");
	//TestSetPointIndices(Model, "Indexed.vtp");

	//TestDownsample(Model, "Downsampled.vtp");
	
	//TestComputeNormals(Model, OutputFile);
	//ColoredPointSquare(OutputFile);
	//Cube();
	//TestTransformModel(Model, OutputFile);
	
	//TestReadScannerLocation(Model);
	
	//TestWriteMemoryUsage(Model);	
	
	//TestDeleteTriangles(Model, "Points.vtp");
	
	TestCopy(Model, OutputFile);
			
	return 0;
}

void TestCopy(const ModelFile &Model, const std::string &OutputFile)
{
	ModelFile Model2 = Model;
	Model2.Write(OutputFile);
}

void TestDeleteTriangles(ModelFile &Model, const std::string &OutputFile)
{
	Model.DeleteTriangles();
	Model.Write(OutputFile);
}

void TestWriteMemoryUsage(const ModelFile &Model)
{
	for(unsigned int i = 0; i < 10; i++)
	{
		std::cout << i << std::endl;
		Model.Write("test.vtp");
	}
}

void TestReadScannerLocation(const ModelFile &Model)
{
	vgl_point_3d<double> P;
	bool IsScan = Model.getScannerLocation(P);
	std::cout << "IsScan? " << IsScan << std::endl;
	std::cout << "P: " << P << std::endl;
}

void TestWritePointSquare(const std::string &OutputFile)
{		
	std::vector<vgl_point_3d<double> > Points;
	Points.push_back(vgl_point_3d<double> (-.5, -.5, 0));
	Points.push_back(vgl_point_3d<double> (-.5, .5, 0));
	Points.push_back(vgl_point_3d<double> (.5, -.5, 0));
	Points.push_back(vgl_point_3d<double> (.5, .5, 0));
	ModelFile World;
	World.setCoords(Points);
	
	World.Write(OutputFile);

}

void TestReadWrite(const std::string &InputFile, const std::string &OutputFile)
{
	ModelFile Model;
	Model.Read(InputFile);
	Model.Init();
	Model.Write(OutputFile);
}

void TestColoredPointSquare(const std::string &OutputFile)
{		
	std::vector<vgl_point_3d<double> > Points;
	Points.push_back(vgl_point_3d<double> (-.5, -.5, 0));
	Points.push_back(vgl_point_3d<double> (-.5, .5, 0));
	Points.push_back(vgl_point_3d<double> (.5, -.5, 0));
	Points.push_back(vgl_point_3d<double> (.5, .5, 0));
	
	std::vector<Color<unsigned char> > Colors;
	Colors.push_back(Colors::Red());
	Colors.push_back(Colors::Green());
	Colors.push_back(Colors::Blue());
	Colors.push_back(Colors::Black());
	
	ModelFile World;
	World.setCoords(Points);
	World.setColors(Colors);
	
	World.Write(OutputFile);

}

void TestWriteTriSquare(const std::string &OutputFile)
{
	std::vector<vgl_point_3d<double> > Points;
	Points.push_back(vgl_point_3d<double> (-1, 1, 0));
	Points.push_back(vgl_point_3d<double> (1, 1, 0));
	Points.push_back(vgl_point_3d<double> (1, -1, 0));
	Points.push_back(vgl_point_3d<double> (-1, -1, 0));
	
	std::vector<vector<unsigned int> > VertexArray;
	
	std::vector<unsigned int> Tri1;
	std::vector<unsigned int> Tri2;
	
	Tri1.push_back(0);
	Tri1.push_back(1);
	Tri1.push_back(2);
	
	Tri2.push_back(0);
	Tri2.push_back(2);
	Tri2.push_back(3);
	
	VertexArray.push_back(Tri1);
	VertexArray.push_back(Tri2);
	
	
	ModelFile Square;
	Square.setCoords(Points);
	Square.setVertexLists(VertexArray);
	Square.Write(OutputFile);
}


void TestCube()
{
	std::vector<vgl_point_3d<double> > Points;
	Points.push_back(vgl_point_3d<double> (-1, -1, -1));
	Points.push_back(vgl_point_3d<double> (1, -1, -1));
	Points.push_back(vgl_point_3d<double> (1, 1, -1));
	Points.push_back(vgl_point_3d<double> (-1, 1, -1));
	Points.push_back(vgl_point_3d<double> (-1, -1, 1));
	Points.push_back(vgl_point_3d<double> (1, -1, 1));
	Points.push_back(vgl_point_3d<double> (1, 1, 1));
	Points.push_back(vgl_point_3d<double> (-1, 1, 1));
	
	std::vector<std::vector<unsigned int> > VertexArray;
	
	std::vector<unsigned int> Tri0;
	std::vector<unsigned int> Tri1;
	std::vector<unsigned int> Tri2;
	std::vector<unsigned int> Tri3;
	std::vector<unsigned int> Tri4;
	std::vector<unsigned int> Tri5;
	std::vector<unsigned int> Tri6;
	std::vector<unsigned int> Tri7;
	std::vector<unsigned int> Tri8;
	std::vector<unsigned int> Tri9;
	std::vector<unsigned int> Tri10;
	std::vector<unsigned int> Tri11;
	
	Tri0.push_back(0);
	Tri0.push_back(3);
	Tri0.push_back(2);
	
	Tri1.push_back(2);
	Tri1.push_back(1);
	Tri1.push_back(0);
	
	
	Tri2.push_back(2);
	Tri2.push_back(3);
	Tri2.push_back(7);
	
	Tri3.push_back(7);
	Tri3.push_back(6);
	Tri3.push_back(2);
	
	Tri4.push_back(0);
	Tri4.push_back(4);
	Tri4.push_back(7);
	
	Tri5.push_back(7);
	Tri5.push_back(3);
	Tri5.push_back(0);
	
	Tri6.push_back(1);
	Tri6.push_back(2);
	Tri6.push_back(6);
	
	Tri7.push_back(6);
	Tri7.push_back(5);
	Tri7.push_back(1);
	
	Tri8.push_back(4);
	Tri8.push_back(5);
	Tri8.push_back(6);
	
	Tri9.push_back(6);
	Tri9.push_back(7);
	Tri9.push_back(4);
	
	Tri10.push_back(0);
	Tri10.push_back(1);
	Tri10.push_back(5);
	
	Tri11.push_back(5);
	Tri11.push_back(4);
	Tri11.push_back(0);
	
	VertexArray.push_back(Tri0);
	VertexArray.push_back(Tri1);
	VertexArray.push_back(Tri2);
	VertexArray.push_back(Tri3);
	VertexArray.push_back(Tri4);
	VertexArray.push_back(Tri5);
	VertexArray.push_back(Tri6);
	VertexArray.push_back(Tri7);
	VertexArray.push_back(Tri8);
	VertexArray.push_back(Tri9);
	VertexArray.push_back(Tri10);
	VertexArray.push_back(Tri11);
	
	ModelFile Square;
	Square.setCoords(Points);
	Square.setVertexLists(VertexArray);
	Square.Write("Cube.vtp");

	
}


void TestTransformModel(ModelFile &Model, const std::string &OutputFile)
{
	
	//vnl_double_3x3 R = MakeRotation('x', .3);
	//vgl_vector_3d<double> T(0,0,-5);
	//Transformation Trans(T,R);

	vgl_vector_3d<double> T(0,0,-5);
	Transformation Trans(T);
	Model.Transform(Trans);

	Model.Write(OutputFile);

}

void TestRemoveGround(ModelFile &Model, const std::string &OutputFile)
{
	Model.RemoveGround(4.0);
	Model.Write(OutputFile);
}

void TestProjectOntoXY(ModelFile &Model, const std::string &OutputFile)
{
	Model.ProjectOntoXY();
	Model.Write(OutputFile);
}

void TestComputeNormals(ModelFile &Model, const std::string &OutputFile)
{
	Model.ComputeNormals(.1);
	Model.Write(OutputFile);
	
}

void TestSetPointIndices(ModelFile &Model, const std::string &OutputFile)
{
	Model.setPointIndices();
	Model.Write(OutputFile);
}

void TestDownsample(ModelFile &Model, const std::string &OutputFile)
{
	Model.Downsample(.05);
	Model.Write(OutputFile);
}

void TestSnapToPoints(ModelFile &Messy, ModelFile &Clean, const std::string &OutputFile)
{
	Messy.SnapToPoints(Clean.getCoords());
	Messy.Write(OutputFile);
}