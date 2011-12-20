#include <iostream>
#include <string>
#include <stdlib.h>

#include "ModelFile.h"

#include <Tools/Tools.h>

void Append(ModelFile &Combined, const ModelFile &Single);

//string FileExtension(const string &Filename);	

void CombinePoints(vector<OrientedPoint> &P1, const std::vector<OrientedPoint> &P2);
//void CombinePoints(vector<vgl_point_3d<double> > &P1, const vector<vgl_point_3d<double> > &P2);
//void CombineColors(vector<Color<unsigned char> > &C1, const vector<Color<unsigned char> > &C2);
void CombineVertices(std::vector<std::vector<unsigned int> > &V1, const std::vector<std::vector<unsigned int> > &V2);

int main(int argc, char *argv[])
{
	//check input arguments
	if(argc < (2 + 1))
	{
		cout << "Required Arguments:" << endl
				<< "1) Output file" << endl
				<< "2 - n) Input files" << endl;
		exit(0);
	}

	//parse input arguments
	std::string OutputFilename = argv[1];
	
	std::vector<std::string> InputFiles;
	for(int i = 2; i < argc; i++)
	{
		InputFiles.push_back(argv[i]);
	}

	//setup the first model
	ModelFile Combined;
	
	Combined.Read(InputFiles[0]); // load the first model

	//add the remaining models to the first one
	for(unsigned int i = 1; i < InputFiles.size(); i++)
	{
		ModelFile Single;
		Single.Read(InputFiles[i]);
	
		Append(Combined, Single);
	}

	//save the result
	Combined.Write(OutputFilename);

	std::cout << "Finished combining." << std::endl;
	
	return 0;
}


void CombinePoints(std::vector<OrientedPoint> &P1, const std::vector<OrientedPoint> &P2)
{
	//P1 has P2 appended
	for(unsigned int i = 0; i < P2.size(); i++)
	{
		P1.push_back(P2[i]);
	}

}

/*
void CombinePoints(vector<vgl_point_3d<double> > &P1, const vector<vgl_point_3d<double> > &P2)
{
	//P1 has P2 appended
	for(unsigned int i = 0; i < P2.size(); i++)
	{
		//PercentComplete(i, NumPoints, NumPoints/100.0);
		P1.push_back(P2[i]);
	}

}
*/
		
/*
void CombineColors(vector<Color<unsigned char> > &C1, const vector<Color<unsigned char> > &C2)
{
	//if(Single.NumColors() > 0)
	//{
		for(unsigned int i = 0; i < C2.size(); i++)
		{
			//PercentComplete(i, NumPoints, NumPoints/100.0);
			C1.push_back(C2[i]);
		}
	//}
}
*/

void CombineVertices(std::vector<std::vector<unsigned int> > &V1, const std::vector<std::vector<unsigned int> > &V2)
{
	for(unsigned int i = 0; i < V2.size(); i++)
	{
		//PercentComplete(i, NumPoints, NumPoints/100.0);
		vector<unsigned int> TriVerts = V2[i];
		for(unsigned int j = 0; j < 3; j++)
			TriVerts[j] += V1.size();

		V1.push_back(TriVerts);
	}
}

void Append(ModelFile &Combined, const ModelFile &Single)
{
	//combine points
	//vector<vgl_point_3d<double> > CombinedPoints = Combined.getCoords();
	//CombinePoints(CombinedPoints, Single.getPoints());
	std::vector<OrientedPoint> CombinedPoints = Combined.getPoints();
	CombinePoints(CombinedPoints, Single.getPoints());

	//combine colors
	/*
	vector<Color<unsigned char> > CombinedColors = Combined.getColors();
	CombineColors(CombinedColors, Single.getColors());
	*/
	
	//combine triangles
	//vector<vector<int> > CombinedTriVerts = Combined.getVertexList();
	//CombineVertices(CombinedTriVerts, Single.getVertexList());
	
	//Combined = ModelFile(CombinedPoints, CombinedTriVerts, CombinedColors);; //return by reference
	//Combined = ModelFile(CombinedPoints, CombinedColors); //return by reference
	Combined = ModelFile(Combined); //return by reference
	
}