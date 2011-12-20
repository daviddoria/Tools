#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkLine.h"
#include "vtkXMLPolyDataWriter.h"

#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <fstream>

//test input
/*
0 0 0
1 1 0
2 0 0
3 1 0
4 0 0
5 1 0
6 0 0
*/

struct Point
{
	Point(const double xin, const double yin, const double zin) {x = xin; y = yin; z = zin;}
	double x,y,z;
};

std::vector<Point> ReadPointsFromFile(const std::string &Filename);
vtkSmartPointer<vtkPolyData> ConnectDots(const std::vector<Point> &Points);

int main(int argc, char *argv[])
{
	int NumArgs = argc - 1;
	if(NumArgs != 2)
	{
		std::cout << "Required arguments: PointFilename.txt OutputPrefix" << std::endl;
		exit(-1);
	}
		
	//parse input arguments
	std::string InputFilename = argv[1];
	std::string OutputPrefix = argv[2];
	
	std::vector<Point> Points = ReadPointsFromFile(InputFilename);
	
	//could be much more efficient by using SetLines(lines), writing the polydata, then adding a single point to the line and then calling SetLines again
	std::vector<Point> SubPoints;
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		SubPoints.push_back(Points[i]);
		
		vtkSmartPointer<vtkPolyData> polydata = ConnectDots(SubPoints);
		
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		std::stringstream OutputFilename;
		OutputFilename << OutputPrefix << "_" << i << ".vtp";
		writer->SetFileName(OutputFilename.str().c_str());
		writer->SetInput(polydata);
		writer->Write();	
	}
	
	return 0;
}

std::vector<Point> ReadPointsFromFile(const std::string &Filename)
{
	std::ifstream fin(Filename.c_str());

	if(fin == NULL)
	{
		std::cout << "Cannot open file." << std::endl;
		exit(-1);
	}

	std::vector<std::string> Lines;
	std::string line;

	//read all of the lines in the file
	while(getline(fin, line))
	{
		Lines.push_back(line);
	}

	//create points from the lines
	std::vector<Point> Points;
	for(unsigned int i = 0; i < Lines.size(); i++)
	{
		std::stringstream ssLine;
		ssLine << Lines[i];
		double x,y,z;
		ssLine >> x >> y >> z;
		Point P(x,y,z);
		Points.push_back(P);
	}
	
	return Points;
}

vtkSmartPointer<vtkPolyData> ConnectDots(const std::vector<Point> &Points)
{
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	
	if(Points.size() <= 1)
	{
		return polydata; //empty
	}
	
	vtkSmartPointer<vtkPoints> vtkpoints = vtkSmartPointer<vtkPoints>::New();
	
	//add points to polydata
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		vtkpoints->InsertNextPoint(Points[i].x, Points[i].y, Points[i].z);
	}
	
	polydata->SetPoints(vtkpoints);
	
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	//add line segments
	for(unsigned int i = 1; i < Points.size(); i++) //note i = 1
	{
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		line->GetPointIds()->SetId(0,i - 1); //previous point (start at second point)
		line->GetPointIds()->SetId(1,i); //current point
		lines->InsertNextCell(line);
	}	
			
	polydata->SetLines(lines);
	
	return polydata;
}