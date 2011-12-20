#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <vnl/vnl_vector.h>

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_3x3.h>

#include <vgl/vgl_vector_3d.h>

#include <VXLHelpers/VXLHelpers.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class Transformation
{
private:
	vnl_double_3 Translation_;
	vnl_double_3x3 Rotation_;
	
public:
	
	Transformation() {}

	Transformation(const vnl_double_3 &T, const vnl_double_3x3 &R) : Translation_(T), Rotation_(R) {}

	Transformation(const vgl_vector_3d<double> &T, const vnl_double_3x3 &R) : Translation_(VXLHelpers::vgl_vector_to_vnl_vector(T)), Rotation_(R) {}
	Transformation(const vgl_vector_3d<double> &T) : Translation_(VXLHelpers::vgl_vector_to_vnl_vector(T)) {Rotation_.set_identity();}
	Transformation(const vnl_double_3x3 &R) : Translation_(vnl_double_3 (0,0,0)), Rotation_(R) {}

	/////// Accessors ////////
	//vnl_vector<double> getTranslation() const {return Translation_;}
	vgl_vector_3d<double> getTranslation() const {return VXLHelpers::vnl_vector_to_vgl_vector(Translation_);}
	vnl_double_3x3 getRotation() const {return Rotation_;}

	///////// Mutators //////////
	void setTranslation(vnl_double_3 &T) {Translation_ = T;}
	void setRotation(vnl_double_3x3 &R) {Rotation_ = R;}
	
	////////// Functions /////////
	vgl_vector_3d<double> ApplyTransform(vgl_vector_3d<double> &V) const;

	///////// Input/Output ///////////
	void WriteToFile(const std::string &Filename);
	void ReadFromFile(const std::string &Filename);
};

	////////// External Functions /////////
int NumTransformations(const std::string &Filename);
std::vector<Transformation> ReadAllTransformations(const std::string &Filename);
vnl_matrix<double> ConstructRowMajor(const vnl_vector<double> &in);

Transformation ReadXFormFile(const std::string &filename);
void WriteXFormFile(const std::string &filename);

std::ostream & operator << (std::ostream &output, const Transformation &T);
#endif