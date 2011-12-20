#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include "OrientedPoint.h"

#include <vgl/vgl_point_3d.h>

std::vector<vgl_point_3d<double> > GetOPCoords(const std::vector<OrientedPoint> &V);
std::vector<Color<unsigned char> > GetOPColors(const std::vector<OrientedPoint> &V);

#endif