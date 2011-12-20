#include "Transformation.h"

#include <Tools.h>

///////////////////// Functions ////////////////////////
vgl_vector_3d<double> Transformation::ApplyTransform(vgl_vector_3d<double> &Vin) const
{
  vnl_vector<double> vnl_vec = VXLHelpers::vgl_vector_to_vnl_vector(Vin);
  vnl_vector<double> Vout = Rotation_ * vnl_vec + Translation_;
  return VXLHelpers::vnl_vector_to_vgl_vector(Vout);
}
	

///////////////////////////////////// External Functions ////////////////////////

std::vector<Transformation> ReadAllTransformations(const std::string &Filename)
{
  std::ifstream TransformFile(Filename.c_str());

  std::stringstream NumLinesStream;
  std::string NumLinesString;
  getline(TransformFile, NumLinesString);

  //int NumTrans = NumTransformations(Filename);
  unsigned int NumTrans;
  NumLinesStream << NumLinesString;
  NumLinesStream >> NumTrans;
  assert(NumTrans > 0);

  std::cout << "There are " << NumTrans << " transformations." << std::endl;

  std::vector<Transformation> Transformations(NumTrans);

  for(unsigned int trans = 0; trans < NumTrans; trans++)
  {
          std::stringstream line;
          std::string temp;
          getline(TransformFile, temp);

          line << temp;
          //cout << line.str() << endl;

          vnl_vector<double> Rotation(9);
          vnl_double_3 Translation;

          for(unsigned int t = 0; t < 3; t++)
          {
                  line >> Translation[t];
                  //cout << Translation[t] << endl;
          }

          for(unsigned int r = 0; r < 9; r++)
          {
                  line >> Rotation[r];
                  //cout << Rotation[r] << endl;
          }

          vnl_double_3x3 R = ConstructRowMajor(Rotation);
          //std::cout << "R " << trans << " : " << R << std::endl;
          assert(Tools::VerifyRotationMatrix(R));
          //cout << "Trans " << trans << " is valid." << endl;
          Transformations[trans] = Transformation(Translation, R);
  }

  TransformFile.close();

  return Transformations;

}


int NumTransformations(const std::string &Filename)
{	
  //if there are 9 transformations, this return 9
  std::ifstream TransformFile(Filename.c_str());
  int n = 0;
  std::string line;

  while(getline(TransformFile, line))
  {
    n++;
  }

  TransformFile.close();
  return n;
}

vnl_matrix<double> ConstructRowMajor(const vnl_vector<double> &in)
{
  vnl_matrix<double> out(3, 3);
  out(0,0) = in[0];
  out(0,1) = in[1];
  out(0,2) = in[2];
  out(1,0) = in[3];
  out(1,1) = in[4];
  out(1,2) = in[5];
  out(2,0) = in[6];
  out(2,1) = in[7];
  out(2,2) = in[8];

  return out;
}


std::ostream & operator << (std::ostream &output, const Transformation &T)
{
  output << "T: " << std::endl << T.getTranslation() << std::endl;
  output << "R: " << std::endl << T.getRotation() << std::endl;
  return output;
}

void Transformation::ReadFromFile(const std::string &Filename)
{
  //format is
  //Tx Ty Tz R00 R01 R02 R10 R11 R12 R20 R21 R22
  //where
  //R00 R01 R02
  //R10 R11 R12
  //R20 R21 R22

  std::ifstream TransformStream(Filename.c_str());
  vnl_double_3 T;
  TransformStream >> T;
  vnl_double_3x3 R;
  TransformStream >> R;

  Translation_ = T;
  Rotation_ = R;

  TransformStream.close();

}

void Transformation::WriteToFile(const std::string &Filename)
{
  //format is
  //Tx Ty Tz R00 R01 R02 R10 R11 R12 R20 R21 R22
  //where
  //R00 R01 R02
  //R10 R11 R12
  //R20 R21 R22

  std::ofstream TransformFile(Filename.c_str());
  TransformFile << Translation_ << " " << std::endl;
  TransformFile << Rotation_;

  TransformFile.close();
}

/*
//binary
void ReadXFormFile(string &filename, vector<vgl_point_3d<double> > &points)
{
  FILE *file = fopen(filename.c_str(), "rb");
  while(!feof(file))
  {
          double point[3];
          fread(point, sizeof(double)*3, 1, file);
          points.push_back(vgl_point_3d<double> (point[0], point[1], point[2]));
  }
  fclose(file);
}
*/
