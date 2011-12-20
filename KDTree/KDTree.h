#ifndef KDTREE_H
#define KDTREE_H

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkKdTree.h>

#include <vector>

#include <vgl/vgl_point_3d.h>


class KDTree
{
	bool Valid_;
	unsigned int NumPoints;
	vtkSmartPointer<vtkDataArray> ListOfPoints;
	vtkSmartPointer<vtkKdTree> PointTree;
	
	std::vector<vgl_point_3d<double> > Points_;
	
	public:
		KDTree()
		{
			Valid_ = false;
		}
		
		KDTree(const std::vector<vgl_point_3d<double> > &Points)
		{
			//Model_ = Model;
			Points_ = Points;
			CreateTree(Points);
		}
		
		void CreateTree(const std::vector<vgl_point_3d<double> > &Points);
		std::vector<vgl_point_3d<double> > KNearest(const vgl_point_3d<double> &Point, const unsigned int k) const;
		vgl_point_3d<double> ClosestPoint(const vgl_point_3d<double> &Point) const;
		
		unsigned int ClosestPointIndex(const vgl_point_3d<double> &Point) const
		{
			std::vector<unsigned int> knearest = KNearestIndices(Point, 1);
			if(knearest.size() == 0)
			{
				std::cout << "There are no closest points (probably a bug).\n";
				return 0;
			}
			
			return knearest[0];
		}
		
		std::vector<unsigned int> KNearestIndices(const vgl_point_3d<double> &Point, const unsigned int k) const;
		std::vector<unsigned int> IndicesWithinRadius(const double R, const vgl_point_3d<double> &Point) const;
		//vector<OrientedPoint> KNearest(const OrientedPoint &Point, const unsigned int k) const;
		//vector<OrientedPoint> KNearestOriented(const vgl_point_3d<double> &Point, const unsigned int k) const;
		double ClosestPointDistance(const vgl_point_3d<double> &Point);
		double ClosestNonZeroDistance(const vgl_point_3d<double> &Point);
		
		bool isValid(void) {return Valid_;}
};

namespace KDFuncs
{
	//Use these function if a tree has not already been constructed.
	//They will create a tree then perform the desired operations.
	std::vector<unsigned int> KClosestPointIndices(const unsigned int k, const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points);
	double AveragePointDistance(const std::vector<vgl_point_3d<double> > &Points); //points to themselves distances
	double AveragePointDistance(const std::vector<vgl_point_3d<double> > &TestPoints, const std::vector<vgl_point_3d<double> > &Points); //two point sets
	double AveragePointDistanceLimited(const std::vector<vgl_point_3d<double> > &TestPoints, const std::vector<vgl_point_3d<double> > &Points, const double Limit); //two point sets
	double MedianDistance(const std::vector<vgl_point_3d<double> > &Points);
}


#endif
