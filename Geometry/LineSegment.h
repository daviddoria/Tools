#ifndef LINESEGMENT_H
#define LINESEGMENT_H

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

#include <vector>

#include <iostream>
#include <assert.h>

class LineSegment
{
private:
	
		
public:

	std::vector<vgl_point_3d<double> > EndPoints;
		
	LineSegment(const vgl_point_3d<double> &P0, const vgl_point_3d<double> &P1)
	{
		EndPoints.push_back(P0);
		EndPoints.push_back(P1);
	}

};

void WriteLineFile(const std::vector<LineSegment> &LineSegments, const std::string &Filename);

#endif
