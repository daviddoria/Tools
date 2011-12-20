#include "Circle.h"

std::ostream & operator << (std::ostream &output, const Circle &Circle)
{
	output << "Circle" << std::endl
			<< "------" << std::endl
			<< "Center: " << Circle.getCenter() << std::endl
			<< "Radius: " << Circle.getRadius() << std::endl;

  return output;

}
