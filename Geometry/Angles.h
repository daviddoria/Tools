#ifndef ANGLES_H
#define ANGLES_H

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_vector_2d.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_double_3.h>

#include <VXLHelpers/VXLHelpers.h>

#include <iostream>
#include <cmath>


////// Conversions //////////
double deg2rad(const double deg);
double rad2deg(const double rad);

///// Polar and Rectangular (2D) //////
vgl_vector_2d<double> Polar2Rect(const double Theta);
double Rect2Polar(const vgl_vector_2d<double> &V);

//////// Spherical and Rectangular (3D) ///////
vgl_vector_3d<double> Sphere2Rect(const double Phi, const double Theta);
void Rect2Sphere(const vgl_vector_3d<double> &Rect, double &Phi, double &Theta);
void AddSpherical(const double PhiIn, const double ThetaIn, double &Phi, double &Theta);
vgl_vector_3d<double> AxisAngle(const vgl_vector_3d<double> &v, const vgl_vector_3d<double> &AxisIn, double Angle);

double IsRightHanded(const vnl_double_3 &x, const vnl_double_3 &y, const vnl_double_3 &z);
double IsRightHanded(const vgl_vector_3d<double> &x, const vgl_vector_3d<double> &y, const vgl_vector_3d<double> &z);

double AngleBetween(const vgl_vector_3d<double> &v1, const vgl_vector_3d<double> &v2);

#endif
