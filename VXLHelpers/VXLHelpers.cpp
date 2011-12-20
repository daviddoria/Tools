#include "VXLHelpers.h"

#include <vnl/vnl_cross.h>
#include <vnl/vnl_rotation_matrix.h>
#include <vnl/vnl_double_3.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_sphere_3d.h>
#include <vgl/vgl_distance.h>

#include <vnl/algo/vnl_matrix_inverse.h>

#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_real.h>

#include <vbl/vbl_array_2d.h>

//VIL
#include <vil/vil_rgb.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>

#include <Tools.h>

#include <limits>

namespace VXLHelpers
{
	/////////////// Conversions /////////////////
	vgl_vector_3d<double> SubtractVector(const vgl_vector_3d<double> &V1, const vgl_vector_3d<double> &V2)
	{
		//returns V1 - V2
		vgl_vector_3d<double> V(V1.x() - V2.x(), V1.y() - V2.y(), V1.z() - V2.z());
		return V;
	}
	
	vgl_vector_3d<double> AddVector(const vgl_vector_3d<double> &V1, const vgl_vector_3d<double> &V2)
	{
		//returns V1 + V2
		vgl_vector_3d<double> V(V1.x() + V2.x(), V1.y() + V2.y(), V1.z() + V2.z());
		return V;
	}
	
	vector<double> vnl_vector_to_vector(const vnl_vector<double> &v)
	{
		vector<double> vec;
		for(unsigned int i = 0; i < v.size(); i++)
			vec.push_back(v(i));
	
		return vec;
	}
	
	void OutputVNLVector(const vnl_vector<double> &V)
	{
		for(unsigned int i = 0; i < V.size(); i++)
			std::cout << V[i] << " ";
		
		std::cout << std::endl;
	}
	
	////////////// vgl conversions ////////////////
	double Length(const vgl_point_3d<double> &p)
	{
		return vgl_point_to_vgl_vector(p).length();
	}
	
	vgl_vector_3d<double> vgl_point_to_vgl_vector(const vgl_point_3d<double> &p)
	{
		return vgl_vector_3d<double> (p.x(), p.y(), p.z());
	}
	
	vgl_point_3d<double> vgl_vector_to_vgl_point(const vgl_vector_3d<double> &v)
	{
		return vgl_point_3d<double> (v.x(), v.y(), v.z());
	}
	
	
	////////////// Convert between vgl and vnl ////////////////
	
	//vnl_vector<double> vgl_vector_to_vnl_vector(const vgl_vector_3d<double> &vgl_vec)
	vnl_double_3 vgl_vector_to_vnl_vector(const vgl_vector_3d<double> &vgl_vec)
	{
		vnl_vector<double> vnl_vec(3);
		vnl_vec(0) = vgl_vec.x();
		vnl_vec(1) = vgl_vec.y();
		vnl_vec(2) = vgl_vec.z();
		return vnl_vec;
	}
	
	vgl_vector_3d<double> vnl_vector_to_vgl_vector(const vnl_vector<double> &vnl_vec)
	{
		vgl_vector_3d<double> vgl_vec(vnl_vec(0), vnl_vec(1), vnl_vec(2));
		
		return vgl_vec;
	}
	
	//vnl_vector<double> vgl_point_to_vnl_vector(const vgl_point_3d<double> &vgl_point)
	vnl_double_3 vgl_point_to_vnl_vector(const vgl_point_3d<double> &vgl_point)
	{
		vnl_vector<double> vnl_vec(3);
		vnl_vec(0) = vgl_point.x();
		vnl_vec(1) = vgl_point.y();
		vnl_vec(2) = vgl_point.z();
		return vnl_vec;
	}
	
	vgl_point_3d<double> vnl_vector_to_vgl_point(const vnl_vector<double> &vnl_vec)
	{
		vgl_point_3d<double> vgl_vec(vnl_vec(0), vnl_vec(1), vnl_vec(2));
		
		return vgl_vec;
	}
	
	
	//////////// Functions //////////////
	
	vgl_vector_3d<double> cross(const vgl_vector_3d<double> &v1, const vgl_vector_3d<double> &v2)
	{
		vnl_vector<double> vnl_vec1 = vgl_vector_to_vnl_vector(v1);
		vnl_vector<double> vnl_vec2 = vgl_vector_to_vnl_vector(v2);
		
		vnl_vector<double> Crossed = vnl_cross_3d(vnl_vec1, vnl_vec2);
		return vnl_vector_to_vgl_vector(Crossed);
	}
	
	double dot(const vgl_vector_3d<double> &v1, const vgl_vector_3d<double> &v2)
	{
		vnl_double_3 a = vgl_vector_to_vnl_vector(v1);
		vnl_double_3 b = vgl_vector_to_vnl_vector(v2);
	
		double dp = dot_product(a,b);
		return dp;	
	}
	
	
	vgl_vector_3d<double> normalize(const vgl_vector_3d<double> &v)
	{
		vnl_vector<double> vnl_vec = vgl_vector_to_vnl_vector(v);
		vnl_vec.normalize();
		return vnl_vector_to_vgl_vector(vnl_vec);
	}
	
	
	vnl_double_3x3 MakeRotation(const char WhichAxis, const double Angle)
	{
		//expects Angle in radians
		assert((WhichAxis == 'x') || (WhichAxis == 'y') || (WhichAxis == 'z'));
		
		vnl_matrix<double> R(3,3);
		
		if (WhichAxis == 'x')
		{
			R(0,0) = 1;
			R(0,1) = 0;
			R(0,2) = 0;
			R(1,0) = 0;
			R(1,1) = cos(Angle);
			R(1,2) = sin(Angle);
			R(2,0) = 0;
			R(2,1) = -sin(Angle);
			R(2,2) = cos(Angle);
		}
		else if(WhichAxis == 'y')
		{
			R(0,0) = cos(Angle);
			R(0,1) = 0;
			R(0,2) = -sin(Angle);
			R(1,0) = 0;
			R(1,1) = 1;
			R(1,2) = 0;
			R(2,0) = sin(Angle);
			R(2,1) = 0;
			R(2,2) = cos(Angle);
		}
		else if(WhichAxis == 'z')
		{		
			R(0,0) = cos(Angle);
			R(0,1) = sin(Angle);
			R(0,2) = 0;
			R(1,0) = -sin(Angle);
			R(1,1) = cos(Angle);
			R(1,2) = 0;
			R(2,0) = 0;
			R(2,1) = 0;
			R(2,2) = 1;
		}
	
		return R;
	}
	
	vgl_point_3d<double> PointAlong(const vgl_point_3d<double> &P, const vgl_vector_3d<double> &V, const double D)
	{ 
		return P + V*D; 
	}
		
	
	/////////////////// Box /////////////////////////
	//                                 MaxPosition
	//                       |<--width-->|
	//                       O-----------O  ---
	//                      /           /|   ^
	//                     /     T     / |   |
	//                    O-----------O  | height
	//                    |       o   |  |R  |
	//                  L |  centroid |  |   v
	//                    |    F      |  O  ---
	//     Y              |           | /   /_____depth
	//     |   Z          |           |/   /
	//     |  /           O-----------O  ---
	//     | /         MinPosition
	//     O-----X
	
	vgl_plane_3d<double> GetFrontPlane(const vgl_box_3d<double> &B)
	{
		
		vgl_point_3d<double> P0 = B.min_point();
		vgl_point_3d<double> P1 = P0 + vgl_vector_3d<double> (1,0,0);
		vgl_point_3d<double> P2 = P0 + vgl_vector_3d<double> (0,1,0);
		
		return vgl_plane_3d<double> (P0, P1, P2);
	}
	
	vgl_plane_3d<double> GetBackPlane(const vgl_box_3d<double> &B)
	{
		
		vgl_point_3d<double> P0 = B.max_point();
		vgl_point_3d<double> P1 = P0 + vgl_vector_3d<double> (-1,0,0);
		vgl_point_3d<double> P2 = P0 + vgl_vector_3d<double> (0,-1,0);
	
		return vgl_plane_3d<double> (P0, P1, P2);
	}
	
	vgl_plane_3d<double> GetTopPlane(const vgl_box_3d<double> &B)
	{
		
		vgl_point_3d<double> P0 = B.max_point();
		vgl_point_3d<double> P1 = P0 + vgl_vector_3d<double> (0,0,-1);
		vgl_point_3d<double> P2 = P0 + vgl_vector_3d<double> (-1,0,0);
	
		return vgl_plane_3d<double> (P0, P1, P2);
	}
	
	vgl_plane_3d<double> GetBottomPlane(const vgl_box_3d<double> &B)
	{
		
		vgl_point_3d<double> P0 = B.min_point();
		vgl_point_3d<double> P1 = P0 + vgl_vector_3d<double> (1,0,0);
		vgl_point_3d<double> P2 = P0 + vgl_vector_3d<double> (0,0,1);
	
		return vgl_plane_3d<double> (P0, P1, P2);
	}
	
	vgl_plane_3d<double> GetLeftPlane(const vgl_box_3d<double> &B)
	{
		
		vgl_point_3d<double> P0 = B.min_point();
		vgl_point_3d<double> P1 = P0 + vgl_vector_3d<double> (0,1,0);
		vgl_point_3d<double> P2 = P0 + vgl_vector_3d<double> (0,0,1);
	
		return vgl_plane_3d<double> (P0, P1, P2);
	}
	
	vgl_plane_3d<double> GetRightPlane(const vgl_box_3d<double> &B)
	{
		vgl_point_3d<double> P0 = B.max_point();
		vgl_point_3d<double> P1 = P0 + vgl_vector_3d<double> (0,0,-1);
		vgl_point_3d<double> P2 = P0 + vgl_vector_3d<double> (0,-1,0);
	
		return vgl_plane_3d<double> (P0, P1, P2);
	}
	
	vector<vgl_plane_3d<double> > GetAllPlanes(const vgl_box_3d<double> &B)
	{
		vector<vgl_plane_3d<double> > Planes(6);
		Planes[0] = GetFrontPlane(B);
		Planes[1] = GetBackPlane(B);
		Planes[2] = GetLeftPlane(B);
		Planes[3] = GetRightPlane(B);
		Planes[4] = GetTopPlane(B);
		Planes[5] = GetBottomPlane(B);
		return Planes;
	}
	
	bool SameSign(double a, double b)
	{
		return (a < 0) == (b < 0);
	}
	
	vnl_double_3x3 RotationMatrix(vnl_double_3 axis, const double angle)
	{
		//make the norm of axis = angle
		axis = axis.normalize();
		axis = axis*angle;
		
		vnl_double_3x3 R = vnl_rotation_matrix (axis);
		return R;
	}
	
	vnl_double_3x3 MakeScaleMatrix(const double x, const double y, const double z) 
	{
		vnl_double_3x3 s; 
		s.set_identity();
		s(0,0) = x;
		s(1,1) = y;
		s(2,2) = z;
		return s;
	}
	
	vnl_matrix<double> MakeTranslationMatrix(const vgl_vector_3d<double> &T) 
	{
		vnl_matrix<double> Trans(4,4); 
		Trans.set_identity();
		Trans(0,3) = T.x();
		Trans(1,3) = T.y();
		Trans(2,3) = T.z();
		return Trans;
	}
	
	vnl_matrix<double> Get4x4(const vnl_double_3x3 &M)
	{
		vnl_matrix<double> M4(4,4);
		M4.set_identity();
		//M4 = vnl_matrix<double> (M,1,1);
		M4(0,0) = M(0,0);
		M4(0,1) = M(0,1);
		M4(0,2) = M(0,1);
		M4(1,0) = M(1,0);
		M4(1,1) = M(1,1);
		M4(1,2) = M(1,2);
		M4(2,0) = M(2,0);
		M4(2,1) = M(2,1);
		M4(2,2) = M(2,2);
		return M4;
	}
	
	vnl_double_3 Rotate(const vnl_double_3 &V, vnl_double_3 axis, const double angle)
	{
		vnl_double_3x3 R = RotationMatrix(axis, angle);
		
		vnl_double_3 Rotated = R*V;
		
		return Rotated;
	}
	
	vgl_vector_3d<double> Rotate(const vgl_vector_3d<double> &V, vgl_vector_3d<double> ax, const double angle)
	{
		vnl_double_3 v = vgl_vector_to_vnl_vector(V);
		vnl_double_3 axis = vgl_vector_to_vnl_vector(ax);
		
		return vnl_vector_to_vgl_vector(Rotate(v,axis,angle));
	}
	
	vgl_vector_3d<double> ProjectAonB(const vgl_vector_3d<double> &A, const vgl_vector_3d<double> &B)
	{
		vgl_vector_3d<double> a = A;
		vgl_vector_3d<double> b = normalize(B);
	
		return (dot_product(a,b) * b );
	}
	
	/*
	vgl_vector_3d<double> ProjectAonPlane(const vgl_vector_3d<double> &A, const vgl_plane_3d<double> &P)
	{
		//!!! NOT COMPLETE
		//vgl_line_3d_2_points<double> Line(p1, p2);
	
		//vgl_point_3d<double> intersection = vgl_intersection(Line, P);
	
		//vgl_point_3d<T> vgl_closest_point(vgl_plane_3d<T> const& pl, vgl_point_3d<T> const& p);
		
		return vgl_vector_3d<double> (0,0,0);
	}
	*/
	
	vnl_matrix<double> VectorToMatrix(const vnl_double_3 &A)
	{
		vnl_matrix<double> M(3,1);
		M(0,0) = A(0);
		M(1,0) = A(1);
		M(2,0) = A(2);
	
		//cout << M << endl;
		return M;
	}
	
	std::vector<double> EigenValues(const vnl_double_3x3 &M)
	{
		/*
		vnl_real_eigensystem Eigs(M);
	
		vnl_matrix<vcl_complex<double> > Vals = Eigs.D;
	
		vector<double> EVals;
		
		EVals.push_back(vnl_real(Vals)(0,0));
		EVals.push_back(vnl_real(Vals)(1,1));
		EVals.push_back(vnl_real(Vals)(2,2));
		*/
		
		vnl_symmetric_eigensystem<double> Eigs(M);
	
		std::vector<double> EVals;
		
		for(unsigned int i = 0; i < M.columns(); i++)
			EVals.push_back(Eigs.get_eigenvalue(i));
		
		return EVals;
	
	}
	
	std::vector<double> EigenValues(const vnl_matrix<double> &M)
	{
		vnl_symmetric_eigensystem<double> Eigs(M);
	
		std::vector<double> EVals;
		
		for(unsigned int i = 0; i < M.columns(); i++)
			EVals.push_back(Eigs.get_eigenvalue(i));
		
		return EVals;
	
	}
	
	std::vector<vnl_double_3> EigenVectors(const vnl_double_3x3 &M)
	{
		std::vector<vnl_double_3> EVecs;
		
		//assuming input matrix is symmetric
		vnl_symmetric_eigensystem<double> Eigs(M);
	
		EVecs.push_back(Eigs.get_eigenvector(0));
		EVecs.push_back(Eigs.get_eigenvector(1));
		EVecs.push_back(Eigs.get_eigenvector(2));
		//assuming input matrix is NOT symmetric
		/*
		vnl_real_eigensystem Eigs(M);
		
		vnl_double_3x3 Vecs = Eigs.Vreal;
		
		EVecs.push_back(vnl_double_3 (Vecs(0,0), Vecs(1,0), Vecs(2,0)));
		EVecs.push_back(vnl_double_3 (Vecs(0,1), Vecs(1,1), Vecs(2,1)));
		EVecs.push_back(vnl_double_3 (Vecs(0,2), Vecs(1,2), Vecs(2,2)));
		*/
	
		return EVecs;
	}
	
	std::vector<vnl_vector<double> > EigenVectors(const vnl_matrix<double> &M)
	{
		vnl_symmetric_eigensystem<double> Eigs(M);
		
		std::vector<vnl_vector<double> > EVecs;
		for(unsigned int i = 0; i < M.columns(); i++)
			EVecs.push_back(Eigs.get_eigenvector(i));
		
			return EVecs;
		
	}
	
	vnl_double_3x3 OuterProduct(const vgl_vector_3d<double> &A, const vgl_vector_3d<double> &B)
	{
		vnl_double_3 a = vgl_vector_to_vnl_vector(A);
		vnl_double_3 b = vgl_vector_to_vnl_vector(B);
		vnl_double_3x3 OP = OuterProduct(a, b);
		return OP;
	}
	
	vnl_double_3x3 OuterProduct(const vgl_point_3d<double> &A, const vgl_point_3d<double> &B)
	{
		vnl_double_3 a = vgl_point_to_vnl_vector(A);
		vnl_double_3 b = vgl_point_to_vnl_vector(B);
		vnl_double_3x3 OP = OuterProduct(a, b);
		return OP;
	}
	
	vnl_double_3x3 OuterProduct(const vnl_double_3 &A, const vnl_double_3 &B)
	{
		//vnl_double_3x3 A_mat = VectorToMatrix(A);
		//vnl_double_3x3 B_mat = VectorToMatrix(B);
		vnl_matrix<double> A_mat = VectorToMatrix(A);
		vnl_matrix<double> B_mat = VectorToMatrix(B);
		
		vnl_double_3x3 OP = A_mat*B_mat.transpose();
		
		return OP;
	}
	
	void WriteMatrixImageScaled(const vnl_matrix<double> &M, const string &Filename)
	{
		vil_image_view<vxl_byte> Image(M.rows(), M.columns(), 1, 1); //(ni, nj, n_planes, n_interleaved_planes)
		
		//double Max = Tools::VectorMax(Vectorize(M));
		double Max = M.max_value();
		
		for (unsigned j = 0; j < Image.nj(); j++)
		{
			for (unsigned i = 0; i < Image.ni(); i++)
			{
				Image(i,j) = static_cast<vxl_byte>(255 * M(i,j)/Max);
				//cout << "M: " << M(i,j) << endl;
				//cout << "Image: " << Image(i,j) << endl;
			}
		}
		
		vil_save(Image, Filename.c_str());
	}
	
	void ImageFromVector(const std::vector<double> &V, const string &Filename)
	{
		double MaxValue = 0;
		for(unsigned int i = 0; i < V.size(); i++)
		{
			if(V[i] < 0.0)
			{
				std::cout << "Error: all entries must be positive!" << std::endl;
				exit(-1);
			}
			
			if(V[i] > MaxValue)
			{
				MaxValue = V[i];
			}
		}
		
		if(MaxValue == 0.0)
		{
			std::cout << "Error: no values!" << std::endl;
			exit(-1);
		}
		/*
		vil_image_view<vxl_byte> Image(M.rows(), M.columns(), 1, 1); //(ni, nj, n_planes, n_interleaved_planes)
	
		for (unsigned j = 0; j < Image.nj(); j++)
		{
			for (unsigned i = 0; i < Image.ni(); i++)
			{
				Image(i,j) = static_cast<vxl_byte>(255 * M(i,j));
			}
		}
		
		vil_save(Image, Filename.c_str());
		*/
	}
	
	void WriteMatrixImage(const vnl_matrix<double> &M, const string &Filename)
	{
		vil_image_view<vxl_byte> Image(M.rows(), M.columns(), 1, 1); //(ni, nj, n_planes, n_interleaved_planes)
	
		for (unsigned j = 0; j < Image.nj(); j++)
		{
			for (unsigned i = 0; i < Image.ni(); i++)
			{
				Image(i,j) = static_cast<vxl_byte>(255 * M(i,j));
			}
		}
		
		vil_save(Image, Filename.c_str());
	}
	
	vnl_matrix<double> Array2DtoMatrix(const vbl_array_2d<double> &Array)
	{
		vnl_matrix<double> M(Array.rows(), Array.cols());
		
		for (unsigned r = 0; r < Array.rows(); r++)
		{
			for (unsigned c = 0; c < Array.cols(); c++)
			{
				M(r,c) = Array(r,c);
			}
		}
		
		return M;
	}
	
	void Write2DArrayImage(const vbl_array_2d<double> &Array, const string &Filename)
	{
		vnl_matrix<double> M = Array2DtoMatrix(Array);
		WriteMatrixImage(M.transpose(), Filename);
	}
	
	vnl_matrix<double> ReadImageMatrix(const string &Filename)
	{
		vil_image_view<vxl_byte> Image = vil_load(Filename.c_str());
		
		vnl_matrix<double> M(Image.ni(), Image.nj(), 0);
		
		for (unsigned j = 0; j < Image.nj(); j++)
		{
			for (unsigned i = 0; i < Image.ni(); i++)
			{
				M(i,j) = static_cast<double> (Image(i,j));
			}
		}
		
		return M/M.frobenius_norm();
	}
	
	//std::vector<double> Vectorize(const vnl_matrix<double> &M)
	vnl_vector<double> Vectorize(const vnl_matrix<double> &M)
	{
		
		//std::vector<double> V;
		vnl_vector<double> V(M.rows() * M.columns());
		
		for (unsigned j = 0; j < M.rows(); j++)
		{
			for (unsigned i = 0; i < M.columns(); i++)
			{
				//V.push_back(M(i,j));
				V[M.columns() * j + i] = M(i,j);
			}
		}
		
		return V;
	}
	
	double DistToSphere(const vgl_point_3d<double> &P, const vgl_sphere_3d<double> &S)
	{
		double d = vgl_distance(P, S.centre());
		return d - S.radius(); //will be negative if point is inside sphere
	}
	
	std::ostream& operator<<(std::ostream& output, const vgl_sphere_3d<double> &P)
	{
		output << "Center: " << P.centre() << " radius: " << P.radius() << std::endl;
	
		return output;
	}
	
	vnl_double_3 Multiply4x4(const vnl_double_3 &V, const vnl_matrix<double> &M)
	{
		vnl_vector<double> V4(4);
		V4(0) = V(0);
		V4(1) = V(1);
		V4(2) = V(2);
		V4(3) = 1.0;
		
		vnl_vector<double> Result(4);
		Result = M * V4;
		
		Result(0) /= Result(3);
		Result(1) /= Result(3);
		Result(2) /= Result(3);
		
		vnl_double_3 ReturnVec;
		ReturnVec(0) = Result(0);
		ReturnVec(1) = Result(1);
		ReturnVec(2) = Result(2);
		
		return ReturnVec;
	}
	/*
	pair<unsigned int, unsigned int> MinLocation(const vnl_matrix<double> &M)
	{
		double Smallest = numeric_limits<double>::infinity();
		pair<unsigned int, unsigned int> SmallestIndex;
		
		for(unsigned int r = 0; r < M.rows(); r++)
		{
			for(unsigned int c = 0; c < M.cols(); c++)
			{
				if(M(r,c) < Smallest)
				{
					Smallest = M(r,c);
					SmallestIndex = pair<unsigned int, unsigned int> (r, c);
				}
			}
		}
		
		return SmallestIndex;
	}
	*/
	
	vnl_double_3x3 Get3x3SubMatrix(const vnl_matrix<double> &M)
	{
		vnl_double_3x3 R;
		
		for(unsigned r = 0; r < 3; r++)
		{
			for(unsigned c = 0; c < 3; c++)
			{
				R.put(r,c, M.get(r,c));	
			}
		}
		
		return R;
	}
	
	vbl_array_2d<double> ArraySum(const vbl_array_2d<double> &A1, const vbl_array_2d<double> &A2)
	{
		assert(A1.rows() == A2.rows());
		assert(A1.cols() == A2.cols());
		vbl_array_2d<double> SumArray(A1.rows(), A1.cols());
		
		for(unsigned int r = 0; r < A1.rows(); r++)
		{
			for(unsigned int c = 0; c < A1.cols(); c++)
			{
				SumArray(r,c) = A1(r,c) + A2(r,c);
			}
		}
	
		return SumArray;
	}
	
	vbl_array_2d<double> ArrayDifference(const vbl_array_2d<double> &A1, const vbl_array_2d<double> &A2)
	{
		//return A1 - A2
		
		assert(A1.rows() == A2.rows());
		assert(A1.cols() == A2.cols());
		vbl_array_2d<double> DifferenceArray(A1.rows(), A1.cols());
		
		for(unsigned int r = 0; r < A1.rows(); r++)
		{
			for(unsigned int c = 0; c < A1.cols(); c++)
			{
				DifferenceArray(r,c) = A1(r,c) - A2(r,c);
			}
		}
	
		return DifferenceArray;
	}
	
	bool CloseEnough(const vnl_matrix<double> &M1, const vnl_matrix<double> &M2, const double eps)
	{
		unsigned int NumRows = M1.rows();
		unsigned int NumCols = M1.columns();
		if((M2.rows() != NumRows) || (M2.columns() != NumCols))
		{
			std::cout << "Dimensions do not match!" << std::endl;
			return false;	
		}
		
		for(unsigned int r = 0; r < NumRows; r++)
		{
			for(unsigned int c = 0; c < NumCols; c++)
			{
				if(fabs(M1(r,c) - M2(r,c)) > eps)
				{
					std::cout << "Failed comparison: " << "M1: " << M1(r,c) << " M2: " << M2(r,c) << " diff: " << fabs(M1(r,c) - M2(r,c)) << std::endl;
					return false;
				}
			}
		}
		return true;	
	}
	
	bool CloseEnough(const vnl_vector<double> &v1, const vnl_vector<double> &v2, const double eps)
	{
		for(unsigned int i = 0; i < v1.size(); i++)
		{
			if(fabs(v1[i] - v2[i]) > eps)
				return false;
		}
		return true;
	}
	
	bool CloseEnough(const vnl_double_3 &v1, const vnl_double_3 &v2, const double eps)
	{
		for(unsigned int i = 0; i < 3; i++)
		{
			if(fabs(v1[i] - v2[i]) > eps)
				return false;
		}
		return true;
	}
	
	bool CloseEnough(const double a, const double b, const double eps)
	{
		if(fabs(a-b) > eps)
			return false;
		
		return true;
	}
	
	vnl_matrix<double> Reshape(const vnl_vector<double> &V, const unsigned int rows, const unsigned int cols)
	{
		//This function reshapes a vector into a matrix using row major construction.
		
		//check input sizes
		if(V.size() != rows*cols)
		{
			std::cout << "Data sizes do not match!" << std::endl;
		}
		
		vnl_matrix<double> M(rows, cols);
		
		unsigned int counter = 0;
		for(unsigned int r = 0; r < rows; r++)
		{
			for(unsigned int c = 0; c < cols; c++)
			{
				M(r,c) = V[counter];
				counter++;
			}
		}
		
		return M;
	}
	
	vnl_matrix<double> MatrixPower(const vnl_matrix<double> &M, const double MatPow)
	{
		vnl_symmetric_eigensystem<double> Eigs(M);
		
		vnl_matrix<double> V = Eigs.V;
		vnl_diag_matrix<double> D = Eigs.D;
		
		vnl_matrix<double> DRaised(D.rows(), D.columns(), 0.0);
		for(unsigned int i = 0; i < D.rows(); i++)
		{
			DRaised(i,i) = pow(D(i,i), MatPow);
		}
		
		vnl_matrix<double> Raised = V * DRaised * vnl_matrix_inverse<double>(V);
		
		return Raised;
	}

} //end namespace