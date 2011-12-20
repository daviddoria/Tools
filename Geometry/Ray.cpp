#include "Ray.h"

#include <limits>
#include <cmath>

#include <vgl/vgl_line_3d_2_points.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_triangle_3d.h>
#include <vgl/vgl_intersection.h>
#include <vgl/vgl_distance.h>

#include <VXLHelpers/VXLHelpers.h>

#include <Tools.h>

//////////////// Operators ////////////////
ostream& operator << (ostream &output, const Ray &Ray)
{
  output << "Ray:" << endl
	 << "-----" << endl
	 << "Origin: " << Ray.getOrigin() << endl
	 << "Direction: " << Ray.getDirection() << endl << endl;
  return output;
}

//////////// Constructors //////////////
Ray::Ray(const vgl_point_3d<double> &origin, const vgl_vector_3d<double> &direction)
{
  Origin_ = origin;
  Direction_ = normalized(direction);

  //if for some reason the length of the vector of the ray's direction is not 1, mark the ray as invalid
  if(fabs(Direction_.length() - 1.0) < 1e-6)
  {
    valid_ = true;
  }
  else
  {
    valid_ = false;
  }
}

Ray::Ray(const vgl_point_3d<double> &P1, const vgl_point_3d<double> &P2)
{
  assert(P1 != P2);

  Origin_ = P1;
  vgl_vector_3d<double> d = P2 - P1;
  Direction_ = d/d.length();

  //if for some reason the length of the vector of the ray's direction is not 1, mark the ray as invalid
  if(fabs(Direction_.length() - 1.0) < 1e-6)
  {
    valid_ = true;
  }
  else
  {
    valid_ = false;
  }
}

/*
//!!!when is this used?
double Ray::DistanceAlong(const vgl_point_3d<double> &P) const
{
	
	//if ( fabs((P-Origin_).Dot(Direction_)) > epsilon)
	if ( fabs(dot_product((P-Origin_),Direction_)) > epsilon)
	{
		cout << "Point is not on ray!" << endl;
		return numeric_limits<double>::infinity();
	}
	else
	{
		return (P-Origin_).length();
	}

}
*/

bool Ray::IntersectBox(const vgl_box_3d<double> &Box) const
{
	if(IntersectBoxFrontPlane(Box) || 
		  IntersectBoxBackPlane(Box) || 
		  IntersectBoxTopPlane(Box) || 
		  IntersectBoxBottomPlane(Box) || 
		  IntersectBoxLeftPlane(Box) || 
		  IntersectBoxRightPlane(Box))
	{
		return true;
	}
	else
		return false;
}

bool Ray::IntersectBoxFrontPlane(const vgl_box_3d<double> &B) const
{
	
	vgl_point_3d<double> Intersection; 
	bool bIntersectPlane = IntersectPlane(VXLHelpers::GetFrontPlane(B), Intersection);
	
	if(!bIntersectPlane)
		return false;
	
	if( Tools::IsInf(Intersection.x()) || Tools::IsInf(Intersection.y()) || Tools::IsInf(Intersection.z()))
		return false;
	
	
	//check x and y
	if( (Intersection.x() < B.min_x()) || (Intersection.x() > B.max_x()) || (Intersection.y() < B.min_y()) || (Intersection.y() > B.max_y()) )
		return false;
	else
		return true;
}

bool Ray::IntersectBoxBackPlane(const vgl_box_3d<double> &B) const
{
	
	
	vgl_point_3d<double> Intersection;
	bool bIntersectPlane = IntersectPlane(VXLHelpers::GetBackPlane(B), Intersection);
	
	if(!bIntersectPlane)
		return false;
	
	if( Tools::IsInf(Intersection.x()) || Tools::IsInf(Intersection.y()) || Tools::IsInf(Intersection.z()))
		return false;
	
	//check x and y
	if( (Intersection.x() < B.min_x()) || (Intersection.x() > B.max_x()) || (Intersection.y() < B.min_y()) || (Intersection.y() > B.max_y()) )
		return false;
	else
		return true;
}

bool Ray::IntersectBoxTopPlane(const vgl_box_3d<double> &B) const
{
	vgl_point_3d<double> Intersection;
	bool bIntersectPlane = IntersectPlane(VXLHelpers::GetTopPlane(B), Intersection); 
	
	if(!bIntersectPlane)
		return false;
	
	if( Tools::IsInf(Intersection.x()) || Tools::IsInf(Intersection.y()) || Tools::IsInf(Intersection.z()))
		return false;
	
	//check x and z
	if( (Intersection.x() < B.min_x()) || (Intersection.x() > B.max_x()) || (Intersection.z() < B.min_z()) || (Intersection.z() > B.max_z()) )
		return false;
	else
		return true;
}

bool Ray::IntersectBoxBottomPlane(const vgl_box_3d<double> &B) const
{
	vgl_point_3d<double> Intersection;
	bool bIntersectPlane = IntersectPlane(VXLHelpers::GetBottomPlane(B), Intersection);
			
	if(!bIntersectPlane)
		return false;
	
	if( Tools::IsInf(Intersection.x()) || Tools::IsInf(Intersection.y()) || Tools::IsInf(Intersection.z()))
		return false;
	
	//check x and z
	if( (Intersection.x() < B.min_x()) || (Intersection.x() > B.max_x()) || (Intersection.z() < B.min_z()) || (Intersection.z() > B.max_z()) )
		return false;
	else
		return true;
}

bool Ray::IntersectBoxLeftPlane(const vgl_box_3d<double> &B) const
{
	vgl_point_3d<double> Intersection;
	bool bIntersectPlane = IntersectPlane(VXLHelpers::GetLeftPlane(B), Intersection);

	if(!bIntersectPlane)
		return false;
	
	if( Tools::IsInf(Intersection.x()) || Tools::IsInf(Intersection.y()) || Tools::IsInf(Intersection.z()))
		return false;
	
	//check y and z
	if( (Intersection.y() < B.min_y()) || (Intersection.y() > B.max_y()) || (Intersection.z() < B.min_z()) || (Intersection.z() > B.max_z()) )
		return false;
	else
		return true;
}

bool Ray::IntersectBoxRightPlane(const vgl_box_3d<double> &B) const
{
	vgl_point_3d<double> Intersection;
	bool bIntersectPlane = IntersectPlane(VXLHelpers::GetRightPlane(B), Intersection);
	
	if(!bIntersectPlane)
		return false;
	
	if( Tools::IsInf(Intersection.x()) || Tools::IsInf(Intersection.y()) || Tools::IsInf(Intersection.z()))
		return false;
	
	//check y and z
	if( (Intersection.y() < B.min_y()) || (Intersection.y() > B.max_y()) || (Intersection.z() < B.min_z()) || (Intersection.z() > B.max_z()) )
		return false;
	else
		return true;
}

bool Ray::IntersectTriangle(const Triangle &Tri, double &dist, vgl_point_3d<double> &Intersection) const
{
	//vgl_line_3d_2_points<double> Line(Origin_, Origin_ + Direction_);
	vgl_line_segment_3d<double> Line(Origin_, Origin_ + Direction_);
	
	vgl_triangle_3d_intersection_t Intersect = vgl_triangle_3d_line_intersection(Line, Tri.getVertex(0), Tri.getVertex(1), Tri.getVertex(2), Intersection);

	//set outputs
	dist = (Intersection - Origin_).length();
	return Intersect;
}

bool Ray::IntersectPlane(const vgl_plane_3d<double> &Plane, vgl_point_3d<double> &Intersection) const
{
	//intersection could be "behind" ray direction
	vgl_point_3d<double> P2 = Origin_ + Direction_;
	
	vgl_line_3d_2_points<double> Line(Origin_, P2);

	Intersection = vgl_intersection(Line, Plane);
	return true;
	
	//how to check if intersection failed??
}

vgl_point_3d<double> Ray::PointAlong(const double D) const 
{ 
	return Origin_ + Direction_*D; 
}

double Ray::DistanceToPoint(const vgl_point_3d<double> &P) const
{
	//get the orthogonal distance from a point to the ray
	
	if(!valid_)
		assert(0);
	
	//two points that define a line
	vgl_point_3d<double> LineP0 = Origin_;
	vgl_point_3d<double> LineP1 = PointAlong(1.0);

	vgl_line_3d_2_points<double> Line(LineP0, LineP1);

	double D = vgl_distance(Line, P);
	//cout << "Distance: " << D << endl;
	return D;
}

bool Ray::IsSameDirection(const Ray &R) const
{
	if(dot_product(this->getDirection(), R.getDirection()) < 0.0)
		return false;
	else
		return true;
}
