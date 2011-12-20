#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

struct Point
{
	Point(const double xin, const double yin, const double zin) {x = xin; y = yin; z = zin;}
	double x,y,z;
};

std::vector<double> ReadValuesFromFile(const std::string &Filename);
std::vector<Point> ReadPointsFromFile(const std::string &Filename);
		
int main(int argc, char *argv[])
{
	int NumArgs = argc - 1;
	if(NumArgs != 3)
	{
		std::cout << "Required arguments: PointFilename.txt ValueFilename.txt OutputFilename.txt" << std::endl;
		exit(-1);
	}
	
	std::string PointFilename = argv[1];
	std::string ValueFilename = argv[2];
	std::string OutputFilename = argv[3];
	
	std::vector<double> Values = ReadValuesFromFile(ValueFilename);
	std::vector<Point> Points = ReadPointsFromFile(PointFilename);
	
	double offset = 3.0;
	double scale = 5.0;
	std::ofstream fout(OutputFilename.c_str());
	
	for(unsigned int i = 0; i < Values.size(); i++)
	{
		Points[i].z += offset + scale * Values[i];
		fout << Points[i].x << " " << Points[i].y << " " << Points[i].z << std::endl;
	}
		
	fout.close();
	
	return 0;
}


std::vector<double> ReadValuesFromFile(const std::string &Filename)
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
	std::vector<double> Values;
	for(unsigned int i = 0; i < Lines.size(); i++)
	{
		std::stringstream ssLine;
		ssLine << Lines[i];
		double v;
		ssLine >> v;
		Values.push_back(v);
	}
	
	return Values;
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