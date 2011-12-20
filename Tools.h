#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>

//STL Containers
#include <vector>
#include <list>
#include <vector>
#include <set>

#include <algorithm>
#include <numeric>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <vul/vul_file.h>
#include <vnl/vnl_double_3x3.h>

#include <vgl/vgl_point_3d.h>

using std::cout; using std::endl;

namespace Tools
{
	class Log
	{
		private:
			std::streambuf* clog_save;
			std::ofstream ofs;
			std::string LogFilename;
			bool Setup;
		public:
			Log(){Setup = false;}
			Log(const std::string &Filename)
			{
				LogFilename = Filename;
				Setup = true;
			}
			
			void Enable()
			{
				if(Setup)
				{
					clog_save = std::clog.rdbuf();
					ofs.open(LogFilename.c_str());
					std::clog.rdbuf(ofs.rdbuf());
				}
				else
				{
					std::cout << "Log filename has not been set." << std::endl;
				}
			}
		
			void setFilename(const std::string &Filename)
			{
				LogFilename = Filename;
				Setup = true;
			}
			
			~Log()
			{
				//reset the buffer
				std::clog.rdbuf(clog_save);
				ofs.close();
			}
	
	};
	
	bool CompareStrings(const std::string &A, const std::string &B);
	void PercentComplete(const unsigned int Current, const unsigned int Total, const unsigned int Every);
	void Pause(void);
	
	template <typename T>
	void RemoveVectorElement(std::vector<T> &V, unsigned int ind)
	{
		//why does this work???
		V.erase(V.begin() + ind);
	}
	
	template <typename T>
	void OutputObject(const T &obj)
	{
		cout << obj << endl;
	}
	
	template <typename T>
	void OutputVector(const std::vector<T> &V)
	{
		for_each (V.begin(), V.end(), OutputObject<T>);

		cout << endl;
	}
	
	template <typename T>
	std::vector<T> UniqueElements(const std::vector<T> &V)
	{
		std::set<T> s;
		s.insert(V.begin(), V.end());
		std::vector<T> Elements;
		
		for(typename std::set<T>::iterator iter = s.begin(); iter != s.end(); iter++)
			Elements.push_back(*iter);
		
		return Elements;
		
	}
	
	unsigned int FirstUnusedNumber(const std::vector<unsigned int> &V);
	
	template <typename T>
	unsigned int NumUniqueElements(const std::vector<T> &V)
	{
		std::vector<T> Unique = UniqueElements(V);
		return Unique.size();
	}

	template <typename T>
	T VectorMedian(const std::vector<T> &V)
	{
		std::vector<T> a = V;
		//decide which element is the "middle" (round to consider even and odd lengths)
		unsigned int element = floor(a.size() / 2.0);

		//the median is the center element of a sorted array
		std::sort(a.begin(), a.end());
		return a[element];
	}
	
	template <typename T>
	T VectorMax(const std::vector<T> &a)
	{
		typename std::vector<T>::const_iterator pos;
		pos = max_element (a.begin(), a.end());
		return *pos;
	}
	
	template <typename T>
	T VectorMin(const std::vector<T> &a)
	{
		typename std::vector<T>::const_iterator pos;
		pos = min_element (a.begin(), a.end());
		return *pos;
	}
	
	template <typename T>
	void NormalizeVector(std::vector<T> &a)
	{
		T MaxElement = VectorMax(a);
		for(unsigned int i = 0; i < a.size(); i++)
			a[i] /= MaxElement;
	}
	
	template <typename T>
	void CombineVectors(std::vector<T> &V1, const std::vector<T> &V2)
	{
		//append V2 to the end of V1
		V1.insert(V1.end(), V2.begin(), V2.end());
	}
	
	template <typename T>
	std::list<T> VectorToList(const std::vector<T> &V)
	{
		std::list<T> List(V.begin(), V.end());
		return List;
	}
	
	template <typename T>
	std::vector<T> ListToVector(const std::list<T> &L)
	{
		std::vector<T> Vector(L.begin(), L.end());
		return Vector;
	}

	template <typename T>
	void SetInfNanToZero(std::vector<T> &V)
	{
		for(unsigned int i = 0; i < V.size(); i++)
		{
			if(IsInf(V[i]) || IsNaN(V[i]))
				V[i] = static_cast<T>(0.0);
		}
	}
	
	void SetZero(double &d, const double epsilon);
	bool IsInf(const double a);
	bool IsMinusInf(const double a);
	bool IsNaN(const double a);

	double VectorProduct(const std::vector<double> &a);

	template <typename T>
	T VectorSum(const std::vector<T> &a)
	{
		T asum = accumulate(a.begin(), a.end(), 0.0);
		return asum;
	}
	
	double VectorAbsoluteSum(const std::vector<double> &a);
	double VectorSumNaturalLogs(const std::vector<double> &a);
	double VectorSumLog10(const std::vector<double> &a);
	std::string RepeatCharacter(const char a, const unsigned int num);
	std::string ZeroPad(const unsigned int num, const unsigned int rep);

	template <typename T>
	T VectorAverage(const std::vector<T> &a)
	{
		T sum = static_cast<T>(0);
		unsigned int n = a.size();
		for(unsigned int i = 0; i < n; i++)
			sum += a[i];
		
		return sum/static_cast<T>(n);
	}
	double VectorAverageLogs(const std::vector<double> &a);
	std::string FileExtension(const std::string &Filename);
	void AssertNumArgs(const unsigned int argc, const unsigned int num);
	unsigned int NumArgs(const unsigned int argc);
	std::string FileNameWithoutExtension(const std::string &A);
	bool VerifyRotationMatrix(const vnl_double_3x3 &M);

	
	template <typename T>
	bool HasRepeatedAddress(const std::vector<T> &V)
	{
		for(unsigned int i = 0; i < V.size(); i++)
		{
			for(unsigned int j = 0; j < V.size(); j++)
			{
				if(i == j) continue; //don't compare an element to itself!
				
				if(V[i] == V[j])
				{
					cout << i << " and " << j << " are the same!" << endl;
					return true;
				}
			}	
		}
		
		return false;
	}

	template <class T>
	std::vector<T> DereferenceVector(const std::vector<T*> &PointerVec) //no reference
	{
		std::vector<T> ObjectVec;
		for(unsigned int i = 0; i < PointerVec.size(); i++)
			ObjectVec.push_back(*PointerVec[i]);

		return ObjectVec;
	}
	
	unsigned int RandomInt(const unsigned int MAX);
	//std::vector<unsigned int> RandomIntVector(const unsigned int MAX);
	double RandomDouble(const double MIN, const double MAX);
	void SeedRandom(void);
	
	template <class T>
	double sign(T Number)
	{
		if(Number >= static_cast<T>(0))
			return static_cast<T>(1);
		else
			return static_cast<T>(-1);
	}
	
	
	template <typename T>
	std::vector<T> ReorderVector(const std::vector<T> &Things, const std::vector<unsigned int> &Indices)
	{
		assert(Things.size() == Indices.size());
		
		std::vector<T> Ordered(Things.size());
		
		for(unsigned int i = 0; i < Things.size(); i++)
		{
			Ordered[i] = Things[Indices[i]];
		}
		
		return Ordered;
	};
	
	////////////////////////////
	template <typename T>
	struct NumberedItem
	{
		unsigned int index;
		T Item;
	};

	template <typename T>
	bool operator<(NumberedItem<T> NI1, NumberedItem<T> NI2)
	{
		//return NI1.index < NI2.index;
		return NI1.Item < NI2.Item;
	}
	
	template <typename T>
	std::vector<unsigned int> ParallelSortIndices(const std::vector<T> &Things)
	{
		//this function returns the new order of the indices after the sort
		std::vector<NumberedItem<T> > Pairs(Things.size());
		//vector<unsigned int> Indices(Things.size());
		for(unsigned int i = 0; i < Things.size(); i++)
		{
			Pairs[i].index = i;
			Pairs[i].Item = Things[i];
		}
		
		std::sort(Pairs.begin(), Pairs.end());
		
		std::vector<unsigned int> SortedIndices(Things.size());
		for(unsigned int i = 0; i < Things.size(); i++)
			SortedIndices[i] = Pairs[i].index;
		
		return SortedIndices;
		
	}
	
	template <class T, class S>
	void ParallelSort(std::vector<T> &Things1, std::vector<S> &Things2)
	{
		//this function sorts Things1 and reorders Things2 so that it is in the same order as the sorted Things1
		assert(Things1.size() == Things2.size());
		
		unsigned int NumThings = Things1.size();
		
		//create the sortable objects
		std::vector<NumberedItem<T> > Pairs(NumThings);
		for(unsigned int i = 0; i < NumThings; i++)
		{
			Pairs[i].index = i;
			Pairs[i].Item = Things1[i];
		}
		
		std::sort(Pairs.begin(), Pairs.end());
		
		std::vector<unsigned int> SortedIndices(NumThings);
		for(unsigned int i = 0; i < NumThings; i++)
			SortedIndices[i] = Pairs[i].index;
		
		std::vector<T> Things1Out(NumThings);
		std::vector<S> Things2Out(NumThings);
		for(unsigned int i = 0; i < NumThings; i++)
		{
			Things1Out[i] = Pairs[i].Item;
			Things2Out[i] = Things2[SortedIndices[i]];
		}
		
		//return by reference
		Things1 = Things1Out;
		Things2 = Things2Out;
	}
	//////////////////////////


	template <class T>
	void OutputObject(const T &Object, const std::string &Filename)
	{
		std::ofstream fout(Filename.c_str());
	
		fout << Object;

		fout.close();
	}

	unsigned int MaxIndex(const std::vector<double> &V);

	template <class T>
	unsigned int CountOccurances(const std::vector<T> &Vec, const T &Object)
	{
		unsigned int occurances = 0;
		for(unsigned int i = 0; i < Vec.size(); i++)
		{
			if(Vec[i] == Object)
				occurances++;
		}

		return occurances;
	}
	
	template <class T>
	void Reorder(std::vector<T> &Vec, const std::vector<unsigned int> &Indices)
	{
		unsigned int NumElements = Vec.size();
		
		if(NumElements != Indices.size())
		{
			std::cout << "Number of elements does not match number of indices!" << std::endl;
			exit(-1);
		}
		std::vector<T> NewOrder(NumElements);
		for(unsigned int i = 0; i < NumElements; i++)
		{
			NewOrder[i] = Vec[Indices[i]];
		}
		Vec = NewOrder; //return by reference
	}
	
	template <class T>
	void AscendingSort(std::vector<T> &Vec)
	{
		std::sort(Vec.begin(), Vec.end());
	}
	
	template <class T>
	void DescendingSort(std::vector<T> &Vec)
	{
		std::sort(Vec.begin(), Vec.end());
		std::reverse(Vec.begin(), Vec.end());
	}
	
	double ZeroMeanGaussian(const double distance, const double variance);

	void CreateIndexVector(std::vector<unsigned int> &V);
	std::vector<unsigned int> UniqueRandomIndices(const unsigned int MAX, const unsigned int Number);

	bool FileExists(const std::string &Filename);
	
	std::string AddExtension(const std::string &Basename, const std::string &Extension);
	
} // end namespace Tools

#endif