#ifndef RAY_H
#define RAY_H

#include <iostream>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

#include "Triangle.h"

#ifndef EPSILON
	#define EPSILON
	const double epsilon = 1e-6;
#endif

class Ray
{
	//This class contains a point and a direction.  The direction is always stored normalized.
	
	private:
		vgl_vector_3d<double> Direction_; //always stored normalized
		vgl_point_3d<double> Origin_;
		bool valid_;
		
	public:
	
		//////////// Constructors //////////
		Ray(const vgl_point_3d<double> &origin, const vgl_vector_3d<double> &direction);
		Ray(const vgl_point_3d<double> &P1, const vgl_point_3d<double> &P2);
		
		//////////// Accessors /////////////
		vgl_vector_3d<double> getDirection() const {return Direction_;}
		vgl_point_3d<double> getOrigin() const {return Origin_;}
		
		////////////Mutators //////////////
		void setDirection(const vgl_vector_3d<double> &V) {Direction_ = normalized(V);}
		void setOrigin(const vgl_point_3d<double> &P) {Origin_ = P;}
		
		////////// Functions ////////////
		double DistanceAlong(const vgl_point_3d<double> &P) const;
		
		bool IntersectBoxFrontPlane(const vgl_box_3d<double> &B) const;
		bool IntersectBoxBackPlane(const vgl_box_3d<double> &B) const;
		bool IntersectBoxTopPlane(const vgl_box_3d<double> &B) const;
		bool IntersectBoxBottomPlane(const vgl_box_3d<double> &B) const;
		bool IntersectBoxLeftPlane(const vgl_box_3d<double> &B) const;
		bool IntersectBoxRightPlane(const vgl_box_3d<double> &B) const;
		bool IntersectBox(const vgl_box_3d<double> &Box) const;
		
		//vgl_point_3d<double> IntersectTriangle(const Triangle &Tri, double &dist, bool &Intersect) const;
		//vgl_point_3d<double> IntersectPlane(const vgl_plane_3d<double> &Plane, bool &bIntersect) const;

		bool IntersectTriangle(const Triangle &Tri, double &dist, vgl_point_3d<double> &Intersection) const;
		bool IntersectPlane(const vgl_plane_3d<double> &Plane, vgl_point_3d<double> &bIntersection) const;

		vgl_point_3d<double> PointAlong(const double D) const;
		double DistanceToPoint(const vgl_point_3d<double> &P) const;
					
		void Draw(const double len) const;
		
		bool IsSameDirection(const Ray &R) const;
};
	
	//////////// External Operators ///////////
	std::ostream& operator << (std::ostream &output, const Ray &r);

#endif
