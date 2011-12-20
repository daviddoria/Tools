#ifndef VALIDTYPES_H
#define VALIDTYPES_H

#include <iostream>
using namespace std;

template<class T>
struct ValidType
{
	bool Valid;
	T Value;
	ValidType<T>() : Valid(false){}
	ValidType<T>(const T &val) : Valid(true), Value(val){}
	ValidType<T> operator=(const ValidType<T> &rhs);
};

template<class T>
bool operator<(const ValidType<T> &VT1, const ValidType<T> &VT2)
{
	if(VT1.Value < VT2.Value)
		return true;
	else
		return false;

}

template<class T>
bool operator==(const ValidType<T> &VT1, const ValidType<T> &VT2)
{
	if(VT1.Value == VT2.Value)
		return true;
	else
		return false;

}

template<class T>
ostream& operator<<(std::ostream& output, const ValidType<T> &VT)
{
	output << VT.Value << "(" << VT.Valid << ")";

	return output;
}

template<class T>
ValidType<T> ValidType<T>::operator=(const ValidType<T> &rhs)
{
	Value = rhs.Value;
	Valid = rhs.Valid;
	return *this;
}


#endif