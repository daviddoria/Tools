#include "Angles.h"

#include <vnl/vnl_matrix.h>

////// Conversions //////////
double deg2rad(const double deg)
{
  return deg * (M_PI/180);
}

double rad2deg(const double rad)
{
  return rad * (180/M_PI);
}

///// Polar and Rectangular (2D) //////
vgl_vector_2d<double> Polar2Rect(const double Theta)
{
	double x = sin(Theta);
	double y = cos(Theta);
	return vgl_vector_2d<double> (x, y);
}

double Rect2Polar(const vgl_vector_2d<double> &V)
{
	//returns theta in radians
	double theta = atan2(V.y(), V.x());
	return theta;
}


//////// Spherical and Rectangular (3D) ///////

vgl_vector_3d<double> Sphere2Rect(const double Phi, const double Theta)
{
	double x = sin(Theta) * cos(Phi);
	double y = sin(Theta) * sin(Phi);
	double z = cos(Theta);
  	return vgl_vector_3d<double> (x,y,z);
}

void Rect2Sphere(const vgl_vector_3d<double> &V, double &Phi, double &Theta)
{
	Phi = acos(V.x()/sqrt(pow(V.x(),2) + pow(V.y(),2)));
	Theta = acos(V.z()/(sqrt(pow(V.x(),2) + pow(V.y(),2) + pow(V.z(),2))));
}

/*
void AddSpherical(const double Phi1, const double Theta1, const double Phi2, const double Theta2, double &PhiOut, double &ThetaOut)
{
	PhiOut = Phi1 + Phi2;
	
	while(PhiOut > 2*PI)
	{
		PhiOut = PhiOut - 2*PI;
	}
	
	while(PhiOut < -2*PI)
	{
		PhiOut = PhiOut + 2*PI;
	}
	
	ThetaOut = Theta1 + Theta2;
}
*/

vgl_vector_3d<double> AxisAngle(const vgl_vector_3d<double> &V, const vgl_vector_3d<double> &AxisIn, double Angle)
{
	//standing at origin looking along the + axis, positive angles are CCW rotation.
	
	vnl_matrix<double> I(3,3); 
	I.set_identity();
	
	vnl_matrix<double> R(3,3);
	
	vnl_matrix<double> Axis(3,1);
	Axis(0,0) = AxisIn.x();
	Axis(1,0) = AxisIn.y();
	Axis(2,0) = AxisIn.z();
	
	vnl_matrix<double> E(3,3);
	E(1,0) = AxisIn.z();
	E(2,0) = -AxisIn.y();
	E(0,1) = -AxisIn.z();
	E(2,1) = AxisIn.x();
	E(0,2) = AxisIn.y();
	E(1,2) = -AxisIn.x();
		
	R = I * cos(Angle) + (1-cos(Angle))*Axis*Axis.transpose() - E*sin(Angle);
	
	vnl_matrix<double> Ans(3,1);
	vnl_matrix<double> Vec(3,1);
	Vec(0,0) = V.x();
	Vec(1,0) = V.y();
	Vec(2,0) = V.z();
	
	Ans = R*Vec;
	
	vgl_vector_3d<double> Rotated(Ans(0,0), Ans(1,0), Ans(2,0));
	
	return Rotated;
}


double IsRightHanded(const vnl_double_3 &x, const vnl_double_3 &y, const vnl_double_3 &z)
{
	double eps = 1e-3;
	vnl_double_3 crossZ = vnl_cross_3d(x,y);
	vnl_double_3 crossX = vnl_cross_3d(y,z);
	vnl_double_3 crossY = vnl_cross_3d(z,x);

	if( (angle(x,crossX) > eps) || (angle(y,crossY) > eps) || (angle(z,crossZ) > eps) )
		return false;
	else
		return true;

}

double IsRightHanded(const vgl_vector_3d<double> &x, const vgl_vector_3d<double> &y, const vgl_vector_3d<double> &z)
{
	vnl_double_3 vnl_x = VXLHelpers::vgl_vector_to_vnl_vector(x);
	vnl_double_3 vnl_y = VXLHelpers::vgl_vector_to_vnl_vector(y);
	vnl_double_3 vnl_z = VXLHelpers::vgl_vector_to_vnl_vector(z);
	return IsRightHanded(vnl_x, vnl_y, vnl_z);
}

double AngleBetween(const vgl_vector_3d<double> &v1, const vgl_vector_3d<double> &v2)
{
	//return the angle between two vectors in radians.
	return acos(dot_product(v1, v2));
}