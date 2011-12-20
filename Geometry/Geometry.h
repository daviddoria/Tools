#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Ray.h"

#include <vector>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_box_3d.h>

#include <vnl/vnl_double_3x3.h>

namespace geom
{
	vgl_point_3d<double> CenterOfMass(const std::vector<vgl_point_3d<double> *> &Points);
	vgl_point_3d<double> ClosestPoint(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points);
	unsigned int ClosestPointIndex(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points);
	double ClosestPointDist(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points);
	
	vgl_point_3d<double> FarthestPoint(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points);
	vgl_point_3d<double> CenterOfMass(const std::vector<vgl_point_3d<double> > &Points);
	
	vgl_plane_3d<double> FitPlane(const std::vector<vgl_point_3d<double> > &Points);
	vgl_plane_3d<double> BestPlane(const std::vector<vgl_point_3d<double> > &Points);
	
	void WriteAxisFile(const vgl_point_3d<double> &Origin, const vgl_vector_3d<double> &V0, const vgl_vector_3d<double> &V1, const vgl_vector_3d<double> &V2, const std::string &OutputFilename);
	
	std::vector<vgl_point_2d<double> > ScannerImage(const std::vector<vgl_point_3d<double> > &Points, const vgl_point_3d<double> &ScannerLocation);
	
	vgl_vector_3d<double> AverageDirection(const std::vector<vgl_vector_3d<double> > &Vectors);
	
	vnl_double_3x3 FindAligningRotation(const vgl_vector_3d<double> &x, const vgl_vector_3d<double> &y, const vgl_vector_3d<double> &z);
	
	vgl_vector_3d<double> GetOrthogonalVector(const vgl_vector_3d<double> &V);
	vgl_vector_3d<double> GetRandomVec(void);
	
	vgl_vector_3d<double> Project(const vgl_vector_3d<double> &A, const vgl_vector_3d<double> &B);
	vgl_vector_3d<double> Project(const vgl_vector_3d<double> &A, const vgl_plane_3d<double> &P);
	
	std::vector<vgl_point_3d<double> > PointsInsideFrustrum(const vgl_box_3d<double> &BoundingBox, const std::vector<vgl_point_3d<double> > &WorldPoints, const vgl_point_3d<double> &ScannerLocation);
	std::vector<unsigned int> IndicesInsideFrustrum(const vgl_box_3d<double> &BoundingBox, const std::vector<vgl_point_3d<double> > &WorldPoints, const vgl_point_3d<double> &ScannerLocation);
	
	vnl_double_3x3 SampleCovarianceMatrix(const std::vector<vgl_point_3d<double> > &Points);
	double AveragePlaneError(const std::vector<vgl_point_3d<double> > &Points, const vgl_plane_3d<double> &Plane);
	vgl_box_3d<double> BoundingBox(const std::vector<vgl_point_3d<double> > &Points);
	
	std::vector<vgl_point_3d<double> > GetPointsInsideSphere(const std::vector<vgl_point_3d<double> > &Points, const vgl_sphere_3d<double> &Sphere);
	
	void WritePlaneVTP(const vgl_point_3d<double> &Point, const vgl_plane_3d<double> &Plane, const double PlaneSize);
	void WriteEigenVectorVTP(const vgl_point_3d<double> &Point, const std::vector<vgl_vector_3d<double> > &Evecs, const double Length);
	
	double SignedPlaneDistance(const vgl_plane_3d<double> &Plane, const vgl_point_3d<double> &Point);
	
	vgl_point_3d<double> GetFarthestPoint(const std::vector<vgl_point_3d<double> > &Points, const vgl_point_3d<double> &TestPoint);
	
	bool IntersectSphere(const vgl_line_3d_2_points<double> &line, const vgl_sphere_3d<double> &sphere);
	
	double GetNeighborhoodSize(const std::vector<vgl_point_3d<double> > &Points, const double VolumeRatio);
	
	double CorrespondenceError(const std::vector<vgl_point_3d<double> > &ModelPoints, const std::vector<vgl_point_3d<double> > &ScanPoints);
	double CorrespondenceError(const std::vector<std::pair<unsigned int, unsigned int> > &Correspondences, const std::vector<vgl_point_3d<double> > &ModelPoints, const std::vector<vgl_point_3d<double> > &ScanPoints);
			
	vnl_matrix<double> CameraTransform(const Ray &View1, const Ray &View2);

	
	
} // end namespace geom

#endif