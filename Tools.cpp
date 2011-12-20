#include "Tools.h"

#include <cmath>
#include <ctime>
#include <limits>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <set>
#include <vector>
#include <assert.h>

#include <vnl/vnl_double_3.h>

#include <vgl/vgl_distance.h>

#include <vul/vul_file.h>

namespace Tools
{

	void PercentComplete(const unsigned int Current, const unsigned int Total, const unsigned int Every)
	{
		//Every = 10 means "Every 10 percent", so you should see (10%, 20%, 30%, etc)
		if(Current == 0)
			return;
		
		unsigned int OldComplete = 100 * static_cast<double>(Current-1)/static_cast<double>(Total);
		unsigned int Complete = 100 * static_cast<double>(Current)/static_cast<double>(Total);
		
		unsigned int OldGroup = OldComplete/Every;
		unsigned int Group = Complete/Every;
		
		if(Group != OldGroup)
		{
			cout << Complete << "%" << endl;
		}
	}
	
	bool CompareStrings(const std::string &A, const std::string &B)
	{
		//ridiculous convention of == 0 means the strings are equivalent
		if(A.compare(B) == 0)
			return true;
		else
			return false;
	}

	
	void Pause(void)
	{
		std::cin.clear();
		char ch;
		while ( std::cin.get ( ch ) && ch != '\n' );
	}
	
	void SetZero(double &d, const double epsilon)
	{
		if(fabs(d) < epsilon)
			d = 0;
	}
	
	bool IsInf(const double a)
	{
		if(a == std::numeric_limits<double>::infinity())
			return true;
		else
			return false;
	}
	
	bool IsMinusInf(const double a)
	{
		if(-1.0 * a == std::numeric_limits<double>::infinity())
			return true;
		else
			return false;
	}
	
	bool IsNaN(const double a)
	{
		if(a!=a)
			return true;
		return false;
	}

	
	double VectorProduct(const std::vector<double> &a)
	{
		return std::accumulate (a.begin(), a.end(),    // range
				1,                           // initial value
    				std::multiplies<int>());           // operation
	}
	

	
	double VectorAbsoluteSum(const std::vector<double> &a)
	{
		double sum = 0;
		for(unsigned int i = 0; i < a.size(); i++)
			sum += fabs(a[i]);
		
		return sum;
	}
	
	double VectorSumNaturalLogs(const std::vector<double> &a)
	{
		double sum = 0;
		for(unsigned int i = 0; i < a.size(); i++)
			sum += log(a[i]); //natural log
		
		return sum;
	}
	
	
	double VectorSumLog10(const std::vector<double> &a)
	{
		double sum = 0;
		for(unsigned int i = 0; i < a.size(); i++)
			sum += log10(a[i]);
		
		return sum;
	}
	
	double VectorAverageLogs(const std::vector<double> &a)
	{
		double sum = 0;
		unsigned int n = a.size();
		for(unsigned int i = 0; i < n; i++)
			sum += log10(a[i]);
		
		return sum/static_cast<double>(n);
	}
	
	std::string RepeatCharacter(const char a, const unsigned int num)
	{
		std::stringstream Filled;
		Filled << std::setfill(a) << std::setw(num) << a;
			
		return Filled.str();
	}
	
	std::string ZeroPad(const unsigned int num, const unsigned int rep)
	{
		std::stringstream Filled;
		Filled << std::setfill('0') << std::setw(rep) << num;
			
		return Filled.str();
	}
	
	std::string FileExtension(const std::string &Filename)
	{
		std::string ext = Filename.substr(Filename.size() - 3, 3);
		//cout << ext << endl;
		return ext;
	}

/*
	string GetFilename(const string &Filename)
	{
		string ext = Filename.substr(0, Filename.size() - 4);
		cout << ext << endl;
		return ext;
	}
*/

	void Tokenize(const std::string &str, std::vector<std::string>& tokens, const std::string &delimiters = " ")
	{
		// Skip delimiters at beginning.
		std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		// Find first "non-delimiter".
		std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
		
		while (std::string::npos != pos || std::string::npos != lastPos)
		{
			// Found a token, add it to the vector.
			tokens.push_back(str.substr(lastPos, pos - lastPos));
			// Skip delimiters.  Note the "not_of"
			lastPos = str.find_first_not_of(delimiters, pos);
			// Find next "non-delimiter"
			pos = str.find_first_of(delimiters, lastPos);
		}
	}
		
	std::string FileNameWithoutExtension(const std::string &A)
	{
		return vul_file::strip_extension(vul_file::strip_directory(A));
	}
	
	
	
	void AssertNumArgs(const unsigned int argc, const unsigned int num)
	{
		assert(argc == (num+1));
	}
	
	unsigned int NumArgs(const unsigned int argc)
	{
		return argc - 1;
	}
	
	bool VerifyRotationMatrix(const vnl_double_3x3 &M)
	{
		//check row sums
		std::vector<double> RowSums(3);
		for(unsigned int i = 0; i < 3; i++)
		{
			vnl_double_3 r = M.get_row(i);
			double sum = 0;
			for(unsigned int j = 0; j < 3; j++)
				sum += pow(r[j], 2);
			
			RowSums[i] = sum;
			//cout << "Sum of row " << i << " = " << sum << endl;
		}
	
		//check col sums
		std::vector<double> ColSums(3);
		for(unsigned int i = 0; i < 3; i++)
		{
			vnl_double_3 c = M.get_column(i);
			double sum = 0;
			for(unsigned int j = 0; j < 3; j++)
				sum += pow(c[j], 2);
	
			ColSums[i] = sum;
			//cout << "Sum of col " << i << " = " << sum << endl;
		}
	
		double eps = 1e-4;
		for(unsigned int i = 0; i < 3; i++)
		{
			if( fabs(RowSums[i] - 1.0) > eps )
				return false;
			if( fabs(ColSums[i] - 1.0) > eps )
				return false;
		}
			
		return true;
	}

	double RandomDouble(const double MIN, const double MAX)
	{
		assert(MAX > MIN);
		//produce a random double between 0 and 1
		double r = drand48();
		
		return MIN + r*(MAX - MIN);
	}
	
	unsigned int RandomInt(const unsigned int MAX)
	{
		//produce an int from 0 to MAX-1
		return rand() % MAX; 
	}

	std::vector<unsigned int> UniqueRandomIndices(const unsigned int MAX, const unsigned int Number)
	{
		//generate Number unique random indices from 0 to MAX

		SeedRandom();
		assert(Number <= MAX+1); //cannot generate more unique numbers than than the size of the set we are sampling
		
		std::set<unsigned int> S;
		while(S.size() < Number)
			S.insert(RandomInt(MAX));

		std::vector<unsigned int> Indices;
		for(std::set<unsigned int>::iterator iter = S.begin(); iter != S.end(); iter++)
			Indices.push_back(*iter);

		return Indices;
	}
	
	unsigned int FirstUnusedNumber(const std::vector<unsigned int> &Vec)
	{
		std::vector<unsigned int> V = Vec;
		std::sort(V.begin(), V.end());
		for(unsigned int i = 0 ; i < V.size(); i++)
		{
			if(V[i] != i)
				return i;
		}
		
		return V.size();
	}
	
	void SeedRandom(void)
	{
		//if you dont do this, they are the same every time
		srand((unsigned)time(0)); 
		srand48((unsigned)time(0)); 
	}

	double ZeroMeanGaussian(const double distance, const double variance)
	{
		//if the normalizing constant is not used, the Gaussian has value 1 at distance = 0
		//double c = (1.0 / sqrt(variance * 2.0 * M_PI));
		//double e = exp(pow(fabs(distance),2)/(2.0*variance));
		//return c * e;
		
		return exp(-pow(fabs(distance),2)/(2.0*variance));
	}
	
	void CreateIndexVector(std::vector<unsigned int> &V)
	{
		for(unsigned int i = 0; i < V.size(); i++)
			V[i] = i;
	}


	bool FileExists(const std::string &Filename)
	{
		return vul_file::exists(Filename);
	}
	
	std::string AddExtension(const std::string &Basename, const std::string &Extension)
	{
		std::stringstream ss;
		ss << Basename << "." << Extension;
		return ss.str();
	}

	unsigned int MaxIndex(const std::vector<double> &V)
	{
		double m = V[0];
		unsigned int ind = 0;
		for(unsigned int i = 0; i < V.size(); i++)
		{
			if(V[i] > m)
			{
				m = V[i];
				ind = i;
			}
		}

		return ind;
	}
	
} //end namespace Tools	
