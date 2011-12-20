#ifndef ORIENTEDPOINT_H
#define ORIENTEDPOINT_H

#include <iostream>
#include <map>
#include <string>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>

#include "Color.h"
#include "Ray.h"

class OrientedPoint
{
	private:
		vgl_point_3d<double> Coord_;
		vgl_vector_3d<double> Normal_; //always stored normalized
		Color<unsigned char> Color_;

		std::map<std::string, double> DoubleValues;
		
		unsigned int index_;
		bool Valid_;
		bool ValidColor_;
		
	public:
		/////////// Public Members ////////////
		unsigned int type_;//point type (classification)
		//double ConsistencyDistance_;
		
		//////////// Constructors //////////
		OrientedPoint();
		OrientedPoint(const vgl_point_3d<double> &Coord);
		OrientedPoint(const vgl_point_3d<double> &Coord, const vgl_vector_3d<double> &Normal);
		OrientedPoint(const vgl_point_3d<double> &Coord, const Color<unsigned char> &C);
		OrientedPoint(const vgl_point_3d<double> &Coord, const vgl_vector_3d<double> &Normal, const Color<unsigned char> &C);
		void Init();
				
		//////////// Accessors /////////////
		vgl_point_3d<double> getCoord() const {return Coord_;}
		vgl_vector_3d<double> getNormal() const {return Normal_;}
		Color<unsigned char> getColor() const {return Color_;}
		
		double getDoubleValue(const std::string &ValueName) const;
		
		
		bool getDoubleValue(const std::string &ValueName, double &Value) const;
				
		
		unsigned int getIndex() const {return index_;}
		unsigned int getType() const {return type_;}
		bool IsValid(void) const {return Valid_;}
		bool ValidColor(void) const {return ValidColor_;}
		
		vgl_plane_3d<double> getPlane(void)
		{
			return vgl_plane_3d<double>(Normal_, Coord_);
		}
		Ray getRay()
		{
			return Ray(Coord_, Normal_);
		}
		
		////////////Mutators //////////////
		void setCoord(const vgl_point_3d<double> &C) {Coord_ = C;}
		void setNormal(const vgl_vector_3d<double> &N) {Normal_ = normalized(N);}
		void setColor(const Color<unsigned char> &C) {Color_ = C;}
		//void setTotalInformation(const double Info) {TotalInformation_ = Info; RemainingInformation_ = TotalInformation_;}
		//void setRemainingInformation(const double Info) {RemainingInformation_ = Info;}
		void setIndex(const unsigned int index) {index_ = index;}
		void setType(const unsigned int t) {type_ = t;}
		void setValid(const bool v) {Valid_ = v;}
		void setDoubleValue(const std::string &ValueName, const double Value);
		
		void UseInformationPercent(const double Percent);

		void ResetInformation();
		
};

std::ostream& operator<<(std::ostream& output, const OrientedPoint &P);

std::vector<vgl_point_3d<double> > OPVectorToCoordVector(std::vector<OrientedPoint> &OPs);

#endif