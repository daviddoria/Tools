#ifndef VXLHELPERS_H
#define VXLHELPERS_H

#include <iostream>
#include <vector>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_matrix.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_box_3d.h>

#include <vgl/vgl_intersection.h>
#include <vgl/vgl_closest_point.h>

#include <vbl/vbl_array_2d.h>

#include <ValidType.h>

/////////////// Conversions /////////////////
namespace VXLHelpers
{
	vgl_vector_3d<double> SubtractVector(const vgl_vector_3d<double> &V1, const vgl_vector_3d<double> &V2);
	vgl_vector_3d<double> AddVector(const vgl_vector_3d<double> &V1, const vgl_vector_3d<double> &V2);
	
	std::vector<double> vnl_vector_to_vector(const vnl_vector<double> &v);
	void OutputVNLVector(const vnl_vector<double> &V);
	double Length(const vgl_point_3d<double> &p);
	
	vgl_vector_3d<double> vgl_point_to_vgl_vector(const vgl_point_3d<double> &p);
	vgl_point_3d<double> vgl_vector_to_vgl_point(const vgl_vector_3d<double> &v);
	
	//vnl_vector<double> vgl_vector_to_vnl_vector(const vgl_vector_3d<double> &vgl_vec);
	vnl_double_3 vgl_vector_to_vnl_vector(const vgl_vector_3d<double> &vgl_vec);
	//vnl_vector<double> vgl_point_to_vnl_vector(const vgl_point_3d<double> &vgl_point);
	vnl_double_3 vgl_point_to_vnl_vector(const vgl_point_3d<double> &vgl_point);
	
	vgl_vector_3d<double> vnl_vector_to_vgl_vector(const vnl_vector<double> &vnl_vec);
	vgl_point_3d<double> vnl_vector_to_vgl_point(const vnl_vector<double> &vnl_vec);
	
	//////////// Functions //////////////
	vgl_vector_3d<double> cross(const vgl_vector_3d<double> &v1, const vgl_vector_3d<double> &v2);
	double dot(const vgl_vector_3d<double> &v1, const vgl_vector_3d<double> &v2);
	vgl_vector_3d<double> normalize(const vgl_vector_3d<double> &v);
	
	//////////// Transformation Matrices ////////////////////
	vnl_double_3x3 MakeRotation(const char WhichAxis, const double Angle);
	vnl_double_3x3 RotationMatrix(vnl_double_3 axis, const double angle);
	vnl_double_3 Rotate(const vnl_double_3 &V, vnl_double_3 axis, const double angle);
	vgl_vector_3d<double> Rotate(const vgl_vector_3d<double> &V, vgl_vector_3d<double> axis, const double angle);
	vnl_double_3x3 MakeScaleMatrix(const double x, const double y, const double z); 
	vnl_matrix<double> MakeTranslationMatrix(const vgl_vector_3d<double> &T) ;
	vnl_matrix<double> Get4x4(const vnl_double_3x3 &M);
	vnl_double_3 Multiply4x4(const vnl_double_3 &V, const vnl_matrix<double> &M);
	
	vgl_point_3d<double> PointAlong(const vgl_point_3d<double> &P, const vgl_vector_3d<double> &V, const double D);
	
	///////////// Box //////////////
	vgl_plane_3d<double> GetFrontPlane(const vgl_box_3d<double> &B);
	vgl_plane_3d<double> GetBackPlane(const vgl_box_3d<double> &B);
	vgl_plane_3d<double> GetTopPlane(const vgl_box_3d<double> &B);
	vgl_plane_3d<double> GetBottomPlane(const vgl_box_3d<double> &B);
	vgl_plane_3d<double> GetLeftPlane(const vgl_box_3d<double> &B);
	vgl_plane_3d<double> GetRightPlane(const vgl_box_3d<double> &B);
	vector<vgl_plane_3d<double> > GetAllPlanes(const vgl_box_3d<double> &B);
	
	
	bool SameSign(double a, double b);
	
	
	vgl_vector_3d<double> ProjectAonB(const vgl_vector_3d<double> &A, const vgl_vector_3d<double> &B);
	
	vnl_matrix<double> VectorToMatrix(const vnl_double_3 &A);
	
	//3x3
	std::vector<double> EigenValues(const vnl_double_3x3 &M);
	std::vector<vnl_double_3> EigenVectors(const vnl_double_3x3 &M);
	
	// > 3x3
	std::vector<double> EigenValues(const vnl_matrix<double> &M);
	std::vector<vnl_vector<double> > EigenVectors(const vnl_matrix<double> &M);
	
	vnl_double_3x3 OuterProduct(const vgl_vector_3d<double> &A, const vgl_vector_3d<double> &B);
	vnl_double_3x3 OuterProduct(const vnl_double_3 &A, const vnl_double_3 &B);
	vnl_double_3x3 OuterProduct(const vgl_point_3d<double> &A, const vgl_point_3d<double> &B);
	
	vnl_double_3x3 Get3x3SubMatrix(const vnl_matrix<double> &M);
	
	/////////// Convert To Arrays /////////////
	template <typename T>
	void GetArray(const T &VglObject, double* arr)
	{
		arr[0] = VglObject.x();
		arr[1] = VglObject.y();
		arr[2] = VglObject.z();
	}
	
	vnl_matrix<double> Array2DtoMatrix(const vbl_array_2d<double> &Array);
	
	void ImageFromVector(const std::vector<double> &V, const std::string &Filename);
	void WriteMatrixImage(const vnl_matrix<double> &M, const std::string &Filename);
	void WriteMatrixImageScaled(const vnl_matrix<double> &M, const std::string &Filename);
	void Write2DArrayImage(const vbl_array_2d<double> &M, const std::string &Filename);
	//std::vector<double> Vectorize(const vnl_matrix<double> &M);
	vnl_vector<double> Vectorize(const vnl_matrix<double> &M);
	vnl_matrix<double> Reshape(const vnl_vector<double> &V, const unsigned int rows, const unsigned int cols);
	vnl_matrix<double> ReadImageMatrix(const string &Filename);
	
	double DistToSphere(const vgl_point_3d<double> &P, const vgl_sphere_3d<double> &S);
	
	ostream& operator<<(std::ostream& output, const vgl_sphere_3d<double> &P);
	
	//pair<unsigned int, unsigned int> MinLocation(const vnl_matrix<double> &M);
	
	template <typename T>
	std::pair<unsigned int, unsigned int> MinValueLocation(const vbl_array_2d<T> &M, T &MinVal)
	{
		T Smallest;
		
		for(unsigned int r = 0; r < M.rows(); r++)
		{
			for(unsigned int c = 0; c < M.cols(); c++)
			{
				if(M(r,c).Valid == true)
				{
					Smallest = M(r,c);
					Smallest.Valid = true;
					break;
				}
			}
			if(Smallest.Valid == true)
				break;
		}
		
		if(Smallest.Valid == false)
			return pair<unsigned int, unsigned int> (0,0);
		
		std::pair<unsigned int, unsigned int> SmallestIndex;
		
		for(unsigned int r = 0; r < M.rows(); r++)
		{
			for(unsigned int c = 0; c < M.cols(); c++)
			{
				T temp = M(r,c);
				if(temp.Valid == false)
					continue;
				
				if(M(r,c).Value <= Smallest.Value)
				{
					Smallest = M(r,c);
					Smallest.Valid = true;
					SmallestIndex = pair<unsigned int, unsigned int> (r, c);
				}
				
			}
		}
		
		//MinVal = ValidType<T>(Smallest); //return by reference
		MinVal.Value = Smallest.Value;
		MinVal.Valid = Smallest.Valid;
		
		return SmallestIndex;
	}
	
	vbl_array_2d<double> ArraySum(const vbl_array_2d<double> &A1, const vbl_array_2d<double> &A2);
	vbl_array_2d<double> ArrayDifference(const vbl_array_2d<double> &A1, const vbl_array_2d<double> &A2);
	
	bool CloseEnough(const vnl_vector<double> &v1, const vnl_vector<double> &v2, const double eps);
	bool CloseEnough(const vnl_matrix<double> &M1, const vnl_matrix<double> &M2, const double eps);
	bool CloseEnough(const vnl_double_3 &v1, const vnl_double_3 &v2, const double eps);
	bool CloseEnough(const double a, const double b, const double eps);
	
	vnl_matrix<double> MatrixPower(const vnl_matrix<double> &M, const double MatPow);
	
} //end namespace
#endif