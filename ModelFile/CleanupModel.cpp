#include <iostream>
#include <string>
#include <stdlib.h>

#include "ModelFile.h"

#include <Tools/Tools.h>

ModelFile RemoveInsidePoints(const ModelFile &Model, const ModelFile &Bounding, const double Distance);

int main(int argc, char *argv[])
{
	AssertNumArgs(argc, 3);
	/*
		cout << "Required Arguments:" << endl
				<< "1) Output file" << endl
				<< "2 - n) Input files" << endl;
	*/
	
	std::string InputFilename = argv[1];
	std::string BoundingFilename = argv[2];
	std::string OutputFilename = argv[3];
	
	ModelFile Input(InputFilename);
	ModelFile Bounding(BoundingFilename);

	Input.SetCenterOfMass(vgl_point_3d<double>(0,0,0));
	Bounding.SetCenterOfMass(vgl_point_3d<double>(0,0,0));

	ModelFile Output = RemoveInsidePoints(Input, Bounding, .1);

	Output.Write(OutputFilename);

	std::cout << "Finished cleaning." << std::endl;
	
	return 0;
}


/*
ModelFile RemoveInsidePoints(ModelFile &Model)
{
	vector<bool> PointsToKeep(Model.NumPoints(), true);

	for(int i = 0; i < Model.NumPoints(); i++)
	{
		if(i % 100 == 0)
			cout << i << " out of " << Model.NumPoints() << endl;

		if(PointsToKeep[i] == false)
			continue;

		vgl_point_3d<double> P = Model.getPoint(i);

		for(int j = 0; j < Model.NumPoints(); j++)
		{
			if(PointsToKeep[j] == false)
				continue;

			vgl_point_3d<double> P2;

 			if(i!=j)
				P2 = Model.getPoint(j);
			else
				continue;

			if( (AngleBetween(vgl_point_to_vgl_vector(P), vgl_point_to_vgl_vector(P2)) < deg2rad(1)))
			{
				if(vgl_point_to_vgl_vector(P).length() > vgl_point_to_vgl_vector(P2).length())
					PointsToKeep[j] = false;
				else
					PointsToKeep[i] = false;
			}
		}
	}

	vector<vgl_point_3d<double> > NewPoints;
	for(int i = 0; i < Model.NumPoints(); i++)
	{
		if(PointsToKeep[i] ==  true)
			NewPoints.push_back(Model.getPoint(i));
	}

	ModelFile CleanedModel(NewPoints);

	return CleanedModel;
}
*/

ModelFile RemoveInsidePoints(const ModelFile &Model, const ModelFile &Bounding, const double Distance)
{

	unsigned int NumPoints = Model.NumPoints();
	std::vector<bool> KeepPoints(NumPoints, false);

	for(unsigned int point = 0; point < NumPoints; point++)
	{
		if(point % 100 == 0)
			cout << point << " out of " << NumPoints << endl;

		vgl_point_3d<double> P = Model.getPoint(point);

		for(unsigned int tri = 0; tri < Bounding.NumTriangles(); tri++)
		{

			Triangle T = Bounding.getTriangle(tri);
			if(T.DistanceToPoint(P) < Distance)
			{
				KeepPoints[point] = true;
				continue;
			}
		}
	}

	std::vector<vgl_point_3d<double>* > NewPoints;
	for(unsigned int i = 0; i < NumPoints; i++)
	{
		if(KeepPoints[i])
			NewPoints.push_back(new vgl_point_3d<double>(Model.getPoint(i).x(), Model.getPoint(i).y(), Model.getPoint(i).z()));
	}

	ModelFile CleanedModel(NewPoints);

	return CleanedModel;
}