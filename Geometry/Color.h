#ifndef GEOM_COLOR_H
#define GEOM_COLOR_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <map>
#include <cmath>
#include <assert.h>

template <typename T>
class Color
{
	private:
		char WriteMode_; //'s' for screen or 'f' for file
		T R_, G_, B_;
	public:
		//////////// Constructors ////////////
		Color() : WriteMode_('s'), R_(0), G_(0), B_(0) {}
		Color(T r, T g, T b) : WriteMode_('s'), R_(r), G_(g), B_(b) {}
		Color(std::vector<T> &C)
		{
			assert(C.size() == 3);
			R_ = C[0];
			G_ = C[1];
			B_ = C[2];
		}
		
		///////// Accessors ////////////
		T getR() const {return R_;}
		T getG() const {return G_;}
		T getB() const {return B_;}
		char getWriteMode() const {return WriteMode_;}
		
		//////////// Operators /////////////
		void operator/= (T D) {R_ /= D; G_ /= D; B_ /= D;}
		void operator*= (T M) {R_ *= M; G_ *= M; B_ *= M;}
		
		// so it can be used in a map
		bool operator< (Color<T> M) const
		{
			return (R_ < M.R_) || (R_ == M.R_ && G_ < M.G_) || (R_ == M.R_ && G_ == M.G_ && B_ < M.B_);
		}
		bool operator== (Color<T> M) const
		{
			return (R_ == M.R_) && (G_ == M.G_) && (B_ == M.B_);
		}
	
		/////// Functions ///////////
		/*
		vector<T> getRGB() const;
		vector<T> getHSV() const;
		T getH() const {return getHSV()[0];}
		T getS() const {return getHSV()[1];}
		T getV() const {return getHSV()[2];}
		*/
		
		
				
		std::vector<T> HSV2RGB(const T H, const T S, const T V);
		std::vector<T> HSV2RGB(const std::vector<T> &hsv);
		std::vector<T> RGB2HSV(const T R, const T G, const T B);
		std::vector<T> RGB2HSV(const std::vector<T> &rgb);

		//unsigned char* CharArray(void) const;
		unsigned char* CharArray(void) const;
		

};

////////// Colors //////////
namespace Colors
{
	Color<unsigned char> Red();
	Color<unsigned char> Green();
	Color<unsigned char> Blue();
	Color<unsigned char> Black();
	Color<unsigned char> White();
	Color<unsigned char> Grey(unsigned char Lightness);
	Color<unsigned char> Grey(double L);
}
////// External Operators //////////
std::ostream & operator << (std::ostream &output, const Color<double> &C);
std::ostream & operator << (std::ostream &output, const Color<unsigned char> &C);
 	 


//unsigned char* CharArray(Color<unsigned char> &C);
void CharArray(const Color<unsigned char> &C, unsigned char *arr);


// C and D are for "c"har and "d"ouble
void ColorC2ColorD(const Color<unsigned char> &C, Color<double> &D);
void ColorC2ColorDVec(const std::vector<Color<unsigned char> > &C, std::vector<Color<double> > &D);
void ColorD2ColorC(const Color<double> &D, Color<unsigned char> &C);
void ColorD2ColorCVec(const std::vector<Color<double> > &D, std::vector<Color<unsigned char> > &C);

void GenerateColors(std::map<Color<unsigned char>, int> &ColorMap, std::vector<Color<unsigned char> > &ColorList, const unsigned int num);
std::vector<Color<unsigned char> > GenerateColors(const unsigned int num);

std::vector<Color<unsigned char> > BlueRedSpectrum(const unsigned int N);
std::vector<Color<unsigned char> > Spectrum(const unsigned int N);
Color<unsigned char> ColorByValue(const double value);
Color<unsigned char> HSV2RGB(const double H, const double S, const double V);

void WriteHorizontalColorBar(const unsigned int width, const unsigned int height);
void WriteVerticalColorBar(const unsigned int width, const unsigned int height);

#endif
