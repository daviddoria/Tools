#ifndef CIRCLE_H
#define CIRCLE_H

#include <iostream>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

class Circle
{
 private:
  vgl_point_3d<double> Center_;
  double Radius_;

 public:
  ////////// Constructors //////////
  //Circle() : Center_(vgl_point_3d<double>), Radius_(0) {};
  Circle(const vgl_point_3d<double> &Center, const double Radius) : Center_(Center), Radius_(Radius) {}

  /////////// Accessors //////////
  double getRadius() const {return Radius_;}
  vgl_point_3d<double> getCenter() const {return Center_;}

  /////////// Mutators ///////////
  void setRadius(const double Radius) {Radius_ = Radius;}
  void setCenter(const vgl_point_3d<double> &Center) {Center_ = Center;}

};

std::ostream & operator << (std::ostream &output, const Circle &Circle);

#endif
