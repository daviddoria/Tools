#include "OrientedPoint.h"

OrientedPoint::OrientedPoint()
{
	Init();
}

OrientedPoint::OrientedPoint(const vgl_point_3d<double> &Coord)
{
	Init();
	Coord_ = Coord;
}

OrientedPoint::OrientedPoint(const vgl_point_3d<double> &Coord, const vgl_vector_3d<double> &Normal)
{
	Init();
	Coord_ = Coord;

	if(Normal.length() != 0.0)
	{
		vgl_vector_3d<double> NormalizedNormal = Normal;
		normalize(NormalizedNormal);
		Normal_ = NormalizedNormal;
	}
	
}

OrientedPoint::OrientedPoint(const vgl_point_3d<double> &Coord, const Color<unsigned char> &C)
{
	Init();
	Coord_ = Coord;
	Color_ = C;
	ValidColor_ = true;
}

OrientedPoint::OrientedPoint(const vgl_point_3d<double> &Coord, const vgl_vector_3d<double> &Normal, const Color<unsigned char> &C)
{
	Init();
	Coord_ = Coord;

	vgl_vector_3d<double> NormalizedNormal = Normal;
	normalize(NormalizedNormal);
	Normal_ = NormalizedNormal;

	Color_ = C;
	ValidColor_ = true;
}

void OrientedPoint::Init()
{
	setDoubleValue("RemainingInformation", 0.0);
	setDoubleValue("TotalInformation", 0.0);
	Color_ = Color<unsigned char> (0,0,0);
	Normal_ = vgl_vector_3d<double> (0.0, 0.0, 0.0);
	Coord_ = vgl_point_3d<double> (0.0, 0.0, 0.0);
	type_ = 0;
	//ConsistencyDistance_ = -100.0;
	Valid_ = false;
	ValidColor_ = false;
}
		
std::ostream& operator<<(std::ostream& output, const OrientedPoint &P)
{
	output << "Point" << std::endl << "----------" << std::endl << P.getCoord() << std::endl;
	output << "Normal" << std::endl << "----------" << std::endl << P.getNormal() << std::endl;
	output << "Color" << std::endl << "----------" << std::endl << P.getColor() << std::endl;

	return output;
}

void OrientedPoint::UseInformationPercent(const double Percent)
{
	assert(Percent >= 0.0);
	assert(Percent <= 1.0);
	
	double total = getDoubleValue("TotalInformation");
	double update = Percent * total;
	
	double remain = getDoubleValue("RemainingInformation");
	
	if(remain < 0.0)
		remain = 0.0;
	
	setDoubleValue("RemainingInformation", remain);
	//RemainingInformation_ -= Percent * TotalInformation_;
	
}

void OrientedPoint::ResetInformation()
{
	double total = getDoubleValue("TotalInformation");
	setDoubleValue("RemainingInformation", total);
	//RemainingInformation_ = TotalInformation_;
}

void OrientedPoint::setDoubleValue(const std::string &ValueName, const double Value)
{
	DoubleValues[ValueName] = Value;
}

double OrientedPoint::getDoubleValue(const std::string &ValueName) const
{
	double value;
	bool valid = getDoubleValue(ValueName, value);	
	if(valid)
	{
		return value;
	}
	else
	{
		std::cout << "There is no value named " << ValueName << "!" << std::endl;
		assert(0);
		exit(-1);
	}
}

bool OrientedPoint::getDoubleValue(const std::string &ValueName, double &Value) const
{
	
	std::map < std::string, double > ::const_iterator MyIter;

	MyIter = DoubleValues.find(ValueName);
	/*
	if(MyIter == DoubleValues.end())
		return false; //value not found!
	else
	{
		Value = MyIter->second;
		return true;
	}
	*/
}


///////// External Functions /////////////
std::vector<vgl_point_3d<double> > OPVectorToCoordVector(std::vector<OrientedPoint> &OPs)
{
	std::vector<vgl_point_3d<double> > Coords;
	for(unsigned int i = 0; i < OPs.size(); i++)
	{
		Coords.push_back(OPs[i].getCoord());
	}
	return Coords;
}

