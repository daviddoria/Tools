#include "Color.h"

#include <vil/vil_save.h>
#include <vil/vil_image_view.h>


///////// External Operators ///////////
std::ostream & operator << (std::ostream &output, const Color<double> &C)
{
	if(C.getWriteMode() == 's')
		output << "R: " << C.getR() << " G: " << C.getG() << " B: " << C.getB() << std::endl;
	else if(C.getWriteMode() == 'f')
		output << C.getR() << " " << C.getG() << " " << C.getB() << std::endl;
	
	return output;
}

std::ostream & operator << (std::ostream &output, const Color<unsigned char> &C)
{
	if(C.getWriteMode() == 's')
		output << "R: " << int(C.getR()) << " G: " << int(C.getG()) << " B: " << int(C.getB()) << std::endl;
	else if(C.getWriteMode() == 'f')
		output << int(C.getR()) << " " << int(C.getG()) << " " << int(C.getB()) << std::endl;
	
	return output;
}

///////////// Functions ///////////
/*
template <typename T>
vector<T> Color<T>::getRGB() const
{
	vector<T> rgb(3);
	rgb[0] = R_;
	rgb[1] = G_;
	rgb[2] = B_;

	return rgb;
}
*/

/*
template <typename T>
unsigned char* Color<T>::CharArray(void) const
{
	unsigned char arr[3];
	arr[0] = R_;
	arr[1] = G_;
	arr[2] = B_;
	
	return arr;
	
}
*/
/*
unsigned char* Color<unsigned char>::CharArray(void) const
{
	unsigned char arr[3];
	arr[0] = R_;
	arr[1] = G_;
	arr[2] = B_;
	
	return arr;
	
}
*/

void CharArray(const Color<unsigned char> &C, unsigned char *arr)
{
	arr[0] = C.getR();
	arr[1] = C.getG();
	arr[2] = C.getB();
	
}

/*
template <typename T>
vector<T> Color<T>::HSV2RGB(const T H, const T S, const T V)
{
  vector<T> rgb(3);
  
  T var_h, var_i, var_1, var_2, var_3, var_r, var_g, var_b;

  if (S == 0)
    {
      rgb[0] = V;
      rgb[1] = V;
      rgb[2] = V;
    }
  else
    {
      var_h = H * 6;
      if (var_h == 6)
	{
	  var_h = 0; //  H must be < 1
	}

      var_i = floor( var_h );
      var_1 = V * ( 1 - S );
      var_2 = V * ( 1 - S * ( var_h - var_i ) );
      var_3 = V * ( 1 - S * ( 1 - ( var_h - var_i ) ) );

      if (var_i == 0)
	{
	  var_r = V;
	  var_g = var_3;
	  var_b = var_1;
	}
      else if (var_i == 1)
	{   
	  var_r = var_2 ;
	  var_g = V ;
	  var_b = var_1;
	}
      else if (var_i == 2)
	{
	  var_r = var_1 ;
	  var_g = V;
	  var_b = var_3;
	}
      else if (var_i == 3)
	{
	  var_r = var_1 ;
	  var_g = var_2 ;
	  var_b = V;
	}
      else if (var_i == 4)
	{
	  var_r = var_3 ;
	  var_g = var_1 ;
	  var_b = V;
	}
      else
	{
	  var_r = V;
	  var_g = var_1;
	  var_b = var_2;
	}

      rgb[0] = var_r;
      rgb[1] = var_g;
      rgb[2] = var_b;
    }

  return rgb;
}
*/

/*
template <typename T>
vector<T> Color<T>::HSV2RGB(const vector<T> &HSV)
{
  return HSV2RGB(HSV[0], HSV[1], HSV[2]);
}
*/

/*
template <typename T>
vector<T> Color<T>::RGB2HSV(const T R, const T G, const T B)
{
  T H, S, V;
  T max, min;

  //find min and max
  if (R >= G && R >= B)
    max = R;
  else if (G >= R && G >= B)
    max = G;
  else if (B >= R && B >= G)
    max = B;

  if(R <=G && R <= B)
    min = R;
  else if (G <= R && G <= B)
    min = G;
  else if (B <= R && B <= G)
    min = B;

  //find Value
  V = max;

  //find Saturation
  if ( max == 0)
    S = 0;
  else
    S = (max - min) / max;


  //find Hue
  if (max == min)
    H = 0;
  else if(max == R)
    H = 60 * (G - B) / (max-min);
  else if (max == G)
    H = 60 * (B - R) / (max-min) + 120;
  else if (max == B)
    H = 60 * (R - G) / (max-min) + 240;

  H /= 360; //normalize to between 0 and 1
  vector<T> HSV;
  HSV.push_back(H);
  HSV.push_back(S);
  HSV.push_back(V);

  return HSV;

}
*/

void ColorC2ColorD(const Color<unsigned char> &C, Color<double> &D)
{
	D = Color<double>(C.getR() / 255., C.getG() / 255., C.getB() / 255.);
}

void ColorC2ColorDVec(const std::vector<Color<unsigned char> > &C, std::vector<Color<double> > &D)
{
	for(unsigned int i = 0; i < C.size(); i++)
	{
		Color<double> d;
		ColorC2ColorD(C[i], d);
		D.push_back(d);
	}
}

void ColorD2ColorC(const Color<double> &D, Color<unsigned char> &C)
{
	C = Color<unsigned char>(D.getR() * 255., D.getG() * 255., D.getB() * 255.);
}

void ColorD2ColorCVec(const std::vector<Color<double> > &D, std::vector<Color<unsigned char> > &C)
{
	for(unsigned int i = 0; i < D.size(); i++)
	{
		Color<unsigned char> c;
		ColorD2ColorC(D[i], c);
		C.push_back(c);
	}
}


void GenerateColors(std::map<Color<unsigned char>, int> &ColorMap, std::vector<Color<unsigned char> > &ColorList, const unsigned int num)
{

	Color<unsigned char> junk(-1,-1,-1);
	ColorList = std::vector<Color<unsigned char> > (num + 1, junk);
	unsigned int ColorIndex = 0;
		
	for(unsigned char r = 0; r <= 254; r++)
	{
		//printf("r: %u\n", r);
		for(unsigned char g = 0; g <= 254; g++)
		{
			//printf("g: %u\n", g);
			for(unsigned char b = 0; b <= 254; b++)
			{
				//printf("b: %u\n", b);
				std::vector<unsigned char> color(3,0);
				color[0] = r;
				color[1] = g;
				color[2] = b;
				
				ColorMap[color] = ColorIndex;
				ColorList[ColorIndex] = color;
				//ColorList.push_back(color);
						
				ColorIndex++;
				if(ColorIndex == (num - 1))
					return;
			} //end b
		} //end g
	} // end r
	
}

std::vector<Color<unsigned char> > GenerateColors(const unsigned int num)
{

	std::vector<Color<unsigned char> > ColorList(num);
	unsigned int ColorIndex = 0;
		
	for(unsigned char r = 0; r <= 254; r++)
	{
		for(unsigned char g = 0; g <= 254; g++)
		{
			for(unsigned char b = 0; b <= 254; b++)
			{
				/*
				vector<unsigned char> color(3,0);
				color[0] = r;
				color[1] = g;
				color[2] = b;
				*/
				if(ColorIndex == num)
					return ColorList;
				
				ColorList[ColorIndex] = Color<unsigned char> (r,g,b);
				ColorIndex++;

			} //end b
		} //end g
	} // end r
	
	return ColorList;
}

////////////// Colors /////////////
namespace Colors
{
	Color<unsigned char> Red(){return Color<unsigned char> (255,0,0);}
	Color<unsigned char> Green() { return Color<unsigned char>(0,255,0);}
	Color<unsigned char> Blue() { return Color<unsigned char>(0,0,255);}
	Color<unsigned char> Black() { return Color<unsigned char>(0,0,0);}
	Color<unsigned char> White() { return Color<unsigned char>(255,255,255);}
	Color<unsigned char> Grey(unsigned char L) { return Color<unsigned char>(L, L, L);}
	Color<unsigned char> Grey(double L)
	{
		assert(L>=0 && L<=1);
		unsigned char C = round(255. * L);
		return Color<unsigned char>(C,C,C);
	}
} //end namespace Colors

Color<unsigned char> HSV2RGB(const double H, const double S, const double V)
{
//h is 0 to 360
//s and v are 0 to 1
	Color<unsigned char> C;
  //unsigned char var_h, var_i, var_1, var_2, var_3, var_r, var_g, var_b;
	//double var_h, var_i, var_1, var_2, var_3, var_r, var_g, var_b;

    //int hint = H;
	int hi = (int)floor(H/60) % 6;
	double f = H/60. - floor(H/60.);
	/*
	double p = V*(1-S);
	double q = V * (1-f*S);
	double t = V * (1 - (1-f) * S);
	*/
	unsigned char p = V*(1-S) * 255.;
	unsigned char q = V * (1-f*S) * 255.;
	unsigned char t = V * (1 - (1-f) * S) * 255.;

	unsigned char v = V * 255.;
	//unsigned char s = S * 255.;
	
	switch (hi) 
	{
		case 0: 
			C = Color<unsigned char> (v,t,p);
			break;
		case 1:
			C = Color<unsigned char> (q,v,p);
			break;
		case 2:
			C = Color<unsigned char> (p,v,t);
			break;
		case 3:
			C = Color<unsigned char> (p,q,v);
			break;
		case 4:
			C = Color<unsigned char> (t,p,v);
			break;
		case 5:
			C = Color<unsigned char> (v,p,q);
			break;
	}		
	return C;

}

Color<unsigned char> ColorByValue(const double value)
{
	double v = value;
	/*
	assert(intensity >= 0.0);
	assert(intensity <= 1.0);
	*/
	if(v < 0.0)
		v = 0.0;
	
	if(v > 1.0)
		v = 1.0;
	
	double h = 240. + (360.-240.)*v;
	return HSV2RGB(h, 1.0, 1.0);
}

std::vector<Color<unsigned char> > BlueRedSpectrum(const unsigned int N)
{
	//H=240 is blue //should this be 270?
	//H=360 is red
	//return a spectrum of N colors from blue to red
	double step = (360.-240.)/N;
	std::vector<Color<unsigned char> > Colors(N);
	
	for(unsigned i = 0; i < N; i++)
	{
		Colors[i] = HSV2RGB(240. + i*step, 1.0, 1.0);
	}
	
	return Colors;
	
}

std::vector<Color<unsigned char> > Spectrum(const unsigned int N)
{
	
	double step = (360.-90.)/N;
	std::vector<Color<unsigned char> > Colors(N);
	
	for(unsigned i = 0; i < N; i++)
	{
		Colors[i] = HSV2RGB(0. + i*step, 1.0, 1.0);
	}
	
	return Colors;
	
}

void WriteHorizontalColorBar(const unsigned int size, const unsigned int length)
{
	vil_image_view<vxl_byte> Bar(length,size,1,1);
	
	for (unsigned c = 0; c < Bar.nj(); c++)
	{
		for (unsigned r = 0; r < Bar.ni();r++)
		{
			Bar(r,c) = 255. * r/static_cast<double>(Bar.ni());
		}
	}
	
	vil_save(Bar, "GreyColorBar.jpg");
}

void WriteVerticalColorBar(const unsigned int size, const unsigned int length)
{
	//size is how fat the bar is

	unsigned int NumRows = length;
	unsigned int NumCols = size;

	//ni = width
	//nj = height
	//vil_image_view (ni, nj)
	vil_image_view<vxl_byte> Bar(NumCols, NumRows, 1,1);

	std::cout << "NumRows: " << NumRows << " " << Bar.nj() << std::endl;
	std::cout << "NumCols: " << NumCols << " " << Bar.ni() << std::endl;
	
	for (unsigned c = 0; c < NumCols; c++)
	{
		for (unsigned r = 0; r < NumRows; r++)
		{
			double val;
			
			if(r == 0 || c == 0 || r == NumRows-1 || c == NumCols-1) //draw a black border
				val = 0;
			else
				val = 255. - 255. * r/static_cast<double>(NumRows);
			
			Bar(c,r) = val;
		}
	}
	
	vil_save(Bar, "GreyColorBar.jpg");
}