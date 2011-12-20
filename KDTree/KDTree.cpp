#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkKdTree.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_distance.h>

#include <iostream>
#include <fstream>
#include <string>

#include <Tools/Tools.h>

#include "KDTree.h"

void KDTree::CreateTree(const std::vector<vgl_point_3d<double> > &Points)
{
	vtkSmartPointer<vtkPoints> Points3d = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType pid[1];
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		vgl_point_3d<double> Point = Points[i];
		pid[0] = Points3d->InsertNextPoint(Point.x(), Point.y(), Point.z());
		Vertices->InsertNextCell(1,pid);
	}

	vtkIdType numpoints;
	numpoints = Points3d->GetNumberOfPoints();
	//cout << "NumPoints: " << numpoints << endl;
	NumPoints = numpoints;

	ListOfPoints = Points3d->GetData();
	//cout << "NumTuples: " << ListOfPoints->GetNumberOfTuples() << endl;

	PointTree = vtkSmartPointer<vtkKdTree>::New();
	
	PointTree->BuildLocatorFromPoints(Points3d);
	
	Valid_ = true;
}

vgl_point_3d<double> KDTree::ClosestPoint(const vgl_point_3d<double> &Point) const
{
	if(!Valid_)
	{
		std::cout << "You must build the tree before calling ClosestPoint!" << std::endl;
		exit(-1);
	}
	std::vector<vgl_point_3d<double> > NearestVec = KNearest(Point, 1);
	return NearestVec[0];
}

std::vector<vgl_point_3d<double> > KDTree::KNearest(const vgl_point_3d<double> &Point, const unsigned int k) const
{
	if(!Valid_)
	{
		std::cout << "You must build the tree before calling KNearest!" << std::endl;
		exit(-1);
	}
	
	double TestPoint[3] = {Point.x(), Point.y(), Point.z()};
	vtkSmartPointer<vtkIdList> Result = vtkSmartPointer<vtkIdList>::New();
	
	PointTree->FindClosestNPoints(k, TestPoint, Result);
	
	std::vector<vgl_point_3d<double> > NearestPoints;
	for(unsigned int i = 0; i < k; i++)
	{
		int point_ind = static_cast<int>(Result->GetId(i));
		//double* p = ListOfPoints->GetTuple(point_ind);
		//NearestPoints.push_back(vgl_point_3d<double> (p[0], p[1], p[2]));
		NearestPoints.push_back(Points_[point_ind]);
	}
	
	return NearestPoints;
}

std::vector<unsigned int> KDTree::KNearestIndices(const vgl_point_3d<double> &Point, const unsigned int k) const
{
	if(!Valid_)
	{
		std::cout << "You must build the tree before calling KNearestIndices!" << std::endl;
		exit(-1);
	}
	
	std::vector<unsigned int> Indices(k);
	
	double TestPoint[3] = {Point.x(), Point.y(), Point.z()};
	vtkSmartPointer<vtkIdList> Result = vtkSmartPointer<vtkIdList>::New();
	
	PointTree->FindClosestNPoints(k, TestPoint, Result);
	
	for(unsigned int i = 0; i < k; i++)
	{
		int point_ind = static_cast<int>(Result->GetId(i));
		Indices[i] = point_ind;
	}
	
	return Indices;
}

std::vector<unsigned int> KDTree::IndicesWithinRadius(const double R, const vgl_point_3d<double> &Point) const
{
	if(!Valid_)
	{
		std::cout << "You must build the tree before calling IndicesWithinRadius!" << std::endl;
		exit(-1);
	}
	
	//return a list of the indices of all the points within R of Point
	std::vector<unsigned int> Indices;
	
	double TestPoint[3] = {Point.x(), Point.y(), Point.z()};
	vtkSmartPointer<vtkIdList> Result = vtkSmartPointer<vtkIdList>::New();
		
	PointTree->FindPointsWithinRadius(R, TestPoint, Result);
			
	for(unsigned int i = 0; i < static_cast<unsigned int>(Result->GetNumberOfIds()); i++)
	{
		//int point_ind = static_cast<int>(Result->GetId(i));
		unsigned int point_ind = static_cast<unsigned int>(Result->GetId(i));
		Indices.push_back(point_ind);
	}
	
	return Indices;
}

double KDTree::ClosestPointDistance(const vgl_point_3d<double> &Point)
{
	if(!Valid_)
	{
		std::cout << "You must build the tree before calling ClosestPointDistance!" << std::endl;
		exit(-1);
	}
	
	std::vector<unsigned int> ClosestIndices = this->KNearestIndices(Point, 1);
	unsigned int ClosestIndex = ClosestIndices[0];
			
	double dist = vgl_distance(this->Points_[ClosestIndex], Point);
	return dist;
}
	
double KDTree::ClosestNonZeroDistance(const vgl_point_3d<double> &Point)
{
	//This function should be used if the point you want to find the closest neighbor to is itself in the tree
	if(!Valid_)
	{
		std::cout << "You must build the tree before calling ClosestNonZeroDistance!" << std::endl;
		exit(-1);
	}
	
	std::vector<unsigned int> ClosestIndices = this->KNearestIndices(Point, 2);
	
		//find the closest point
	unsigned int ClosestIndex = ClosestIndices[0];
	
		//find the distance to the closest point
	double dist = vgl_distance(this->Points_[ClosestIndex], Point);
	
	if(dist == 0.0)
	{
		ClosestIndex = ClosestIndices[1];
		dist = vgl_distance(this->Points_[ClosestIndex], Point);
	}
		
	return dist;
}

///////////////////////////////////////////////////

namespace KDFuncs
{
	////////////////////// External Functions ////////////////////////
	std::vector<unsigned int> KClosestPointIndices(const unsigned int k, const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points)
	{
		std::vector<double> Distances(Points.size());
		std::vector<unsigned int> Indices(Points.size());
		Tools::CreateIndexVector(Indices);
			
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			double dist = vgl_distance(Points[i], TestPoint);
			Distances[i] = dist;
		}
			
		Tools::ParallelSort(Distances, Indices);
		
		std::vector<unsigned int> Inds;
		for(unsigned int i = 0; i < k; i++)
		{
			if(i >= Indices.size())
				break;
				
			Inds.push_back(Indices[i]);
		}
			
		return Inds;
	}
	
	double AveragePointDistance(const std::vector<vgl_point_3d<double> > &Points)
	{
		assert(Points.size() > 0);
		if(Points.size() == 1)
			return 0.0;
			
		std::vector<double> Distances(Points.size());
		KDTree ModelTree(Points);
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			std::vector<unsigned int> ClosestIndices = ModelTree.KNearestIndices(Points[i], 2);
			unsigned int ClosestIndex = ClosestIndices[1]; //0 should be the same point!
				
			Distances[i] = vgl_distance(Points[i], Points[ClosestIndex]);
		}
			
		return Tools::VectorAverage(Distances);
	}
	
	double AveragePointDistance(const std::vector<vgl_point_3d<double> > &TestPoints, const std::vector<vgl_point_3d<double> > &Points)
	{
		//Find the distance from every TestPoint to Points, then average them.
		assert(Points.size() > 0);
		assert(TestPoints.size() > 0);
	
		std::cout << "There are " << TestPoints.size() << " test points." << std::endl;
		std::cout << "There are " << Points.size() << " points." << std::endl;
		
		std::vector<double> Distances(TestPoints.size());
		KDTree Tree(Points);
		for(unsigned int i = 0; i < TestPoints.size(); i++)
		{
			double dist = Tree.ClosestPointDistance(TestPoints[i]);
			//std::cout << "dist " << i << ": " << dist << std::endl;
			Distances[i] = dist;
		}
			
		return Tools::VectorAverage(Distances);
	}
	
	double AveragePointDistanceLimited(const std::vector<vgl_point_3d<double> > &TestPoints, const std::vector<vgl_point_3d<double> > &Points, const double Limit)
	{
		//Find the distance from every TestPoint to Points, then average them.
		assert(Points.size() > 0);
		assert(TestPoints.size() > 0);
	
		//std::cout << "There are " << TestPoints.size() << " test points." << std::endl;
		//std::cout << "There are " << Points.size() << " points." << std::endl;
		
		std::vector<double> Distances;
		KDTree Tree(Points);
		for(unsigned int i = 0; i < TestPoints.size(); i++)
		{
			double dist = Tree.ClosestPointDistance(TestPoints[i]);
			if(dist < Limit)
				Distances.push_back(dist);
		}
		
		if(Distances.size() == 0)
		{
			std::cout << "Error: no points were within Limit!" << std::endl;
			return 0.0;
		}
		
		return Tools::VectorAverage(Distances);
	}
	
	double MedianDistance(const std::vector<vgl_point_3d<double> > &Points)
	{
		//this should be used as the "mesh resolution" when an actual mesh is not available
		assert(Points.size() > 0);
		if(Points.size() == 1)
			return 0.0;
			
		std::vector<double> Distances(Points.size());
		KDTree ModelTree(Points);
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			std::vector<unsigned int> ClosestIndices = ModelTree.KNearestIndices(Points[i], 2);
			unsigned int ClosestIndex = ClosestIndices[1]; //0 should be the same point!
				
			Distances[i] = vgl_distance(Points[i], Points[ClosestIndex]);
		}
			
		return Tools::VectorMedian(Distances);
	}
	

} //end KDFuncs namespace