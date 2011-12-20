#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include <vector>

#include <vgl/vgl_point_3d.h>

void TriangulateLidar(const vgl_point_3d<double> &ScannerLocation, const std::vector<vgl_point_3d<double> > &Points, std::vector<std::vector<unsigned int> > &TriangleVertexList, const double MaxSize);

#endif