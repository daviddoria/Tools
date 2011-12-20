#include <iostream>
#include <cmath>

#include "Tools.h"
#include <vul/vul_file.h>

using namespace Tools;

void TestPercentComplete();
void TestVectorProduct();
void TestVectorSum();
void TestVectorSumLogs();
void TestVectorSumNaturalLogs();
void TestVectorSumLog10();
void TestRepeatCharacter();
void TestZeroPad();
void TestPause();
void TestCompareStrings();
void Test();
void TestFileNames();
void TestDereferenceVector();
void TestSign();
void TestParallelSortIndices();
void TestParallelSort();
void TestReorderVector();
void TestZeroMeanGaussian();
void TestRandomDouble();
void TestVectorMax();
void TestVectorMin();
void TestOutputVector();
void TestUniqueRandomIndices();
void TestCombineVectors();
void TestRemoveVectorElement();
void TestUniqueElements();
void TestLog();


int main()
{
  //TestPercentComplete();
  //TestVectorProduct();
  //TestVectorSum();
  //TestVectorSumNaturalLogs();
  //TestVectorSumLog10();

  //TestRepeatCharacter();
  //TestZeroPad();
  //TestCompareStrings();
  //TestFileNames();
  //TestDereferenceVector();
  //TestPause();

  //TestSign();
  //Test();
  //TestParallelSortIndices();
  //TestParallelSort();
  //TestReorderVector();

  //TestZeroMeanGaussian();
  //TestRandomDouble();

  //TestVectorMax();
  //TestVectorMin();
  //TestOutputVector();
  //TestUniqueRandomIndices();

  //TestCombineVectors();
  TestRemoveVectorElement();

  //TestUniqueElements();

  //TestLog();

  return 0;
}

void TestLog()
{
  std::string LogFilename = "test.txt";

  Log L(LogFilename);
  cout << "Test cout." << endl;
  std::clog << "Test log." << endl;
}

void TestPercentComplete()
{
  unsigned int total = 1234;

  for(unsigned int i = 0; i < total; i++)
  {
    PercentComplete(i, total, 5);
  }
}


void TestVectorProduct()
{
  std::vector<double> a;
  a.push_back(2);
  a.push_back(3);
  a.push_back(4);

  cout << VectorProduct(a) << endl;
}

void TestVectorSum()
{
  cout << "TestVectorSum()" << endl << "--------" << endl;
  std::vector<double> a;
  a.push_back(2.3);
  a.push_back(3.4);
  a.push_back(4.5);

  cout << "Sum of" << endl;
  OutputVector(a);

  cout << VectorSum(a) << endl;
}


void TestVectorSumNaturalLogs()
{
  std::vector<double> a;
  a.push_back(2);
  a.push_back(3);
  a.push_back(4);

  cout << VectorSumNaturalLogs(a) << endl;
}

void TestVectorSumLog10()
{
  std::vector<double> a;
  a.push_back(2);
  a.push_back(3);
  a.push_back(4);

  cout << VectorSumLog10(a) << endl;
}

void TestRepeatCharacter()
{
  cout << "TestRepeatCharacter()" << endl << "------------" << endl;

  std::string test = RepeatCharacter('0', 5);
  cout << test << endl;
}

void TestZeroPad()
{
  cout << "TestZeroPad()" << endl << "------------" << endl;

  std::string test = ZeroPad(23, 5);
  cout << test << endl;

}
		
void TestPause()
{
  cout << "hello" << endl;
  Pause();
  cout << "goodbye" << endl;
	
}

void Test()
{
  long double P = 3.5/(.1 * sqrt(2.0*3.14159)) * exp(-pow(5.0,2) / (2.0*pow(.1,2)));
  cout << P << endl;
  if (P==0)
  {
    cout << "P is zero!" << endl;
  }
}

void TestCompareStrings()
{
  std::string Hello = "Hello";
  std::string Goodbye = "Goodbye";

  cout << CompareStrings(Hello, Goodbye) << " (Should be false.)" << endl;
  cout << CompareStrings(Hello, Hello) << " (Should be true.)" << endl;
	
}

void TestFileNames()
{
  std::string FullName = "/home/doriad/david.txt";
  cout << vul_file::strip_directory(FullName) << endl; //david.txt
  cout << vul_file::basename(FullName, ".txt") << endl; //david (if you know the extension)
  cout << vul_file::strip_extension(vul_file::strip_directory(FullName)) << endl; //david (if you don't know the extension)
  cout << vul_file::strip_extension(FullName) << endl; // /home/doriad/david

}

void TestDereferenceVector()
{
  std::vector<double*> a(10);
  std::vector<double> b = DereferenceVector(a);

  std::vector<vgl_point_3d<double>*> c(10);
  std::vector<vgl_point_3d<double> > d = DereferenceVector(c);

}

void TestSign()
{
  cout << "TestSign()" << endl << "---------------" << endl;
  double a = -5.2;
  int b = -2;
  double c = 6.7;
  int d = 8;
  cout << a << " : " << sign(a) << endl;
  cout << b << " : " << sign(b) << endl;
  cout << c << " : " << sign(c) << endl;
  cout << d << " : " << sign(d) << endl;
}

void TestParallelSortIndices()
{
  std::vector<double> Numbers;
  Numbers.push_back(3.4);
  Numbers.push_back(4.5);
  Numbers.push_back(1.2);

  std::vector<unsigned int> Indices = Tools::ParallelSortIndices(Numbers);

  Tools::OutputVector(Indices);
}

void TestParallelSort()
{
  std::vector<double> Numbers;
  Numbers.push_back(3.4);
  Numbers.push_back(4.5);
  Numbers.push_back(1.2);

  std::vector<std::string> Names;
  Names.push_back("David");
  Names.push_back("Hayley");
  Names.push_back("Joe");

  Tools::ParallelSort(Numbers, Names);

  /*
  output should be
  1.2 Joe
  3.4 David
  4.5 Hayley
  */

  //Tools::OutputVector(Indices);
  Tools::OutputVector(Numbers);
  Tools::OutputVector(Names);
}


void TestReorderVector()
{
  std::vector<double> Numbers;
  Numbers.push_back(3.4);
  Numbers.push_back(4.5);
  Numbers.push_back(1.2);

  std::vector<unsigned int> Indices = Tools::ParallelSortIndices(Numbers);

  Tools::OutputVector(Indices);

  std::vector<double> OrderedNumbers = Tools::ReorderVector(Numbers, Indices);
  Tools::OutputVector(OrderedNumbers);
	
}

void TestZeroMeanGaussian()
{
  cout << ZeroMeanGaussian(0.0, 1.0) << endl;
  cout << ZeroMeanGaussian(.5, 1.0) << endl;
  cout << ZeroMeanGaussian(-.5, 1.0) << endl;
}

void TestRandomDouble()
{
  for(unsigned int i = 0; i < 10; i++)
  {
    cout << RandomDouble(-5.0, 5.0) << endl;
  }
}

void TestVectorMax()
{
  cout << "TestVectorMaxDouble()" << endl << "-----------" << endl;
  std::vector<double> Numbers;
  Numbers.push_back(3.4);
  Numbers.push_back(4.5);
  Numbers.push_back(1.2);

  OutputVector(Numbers);
  cout << endl;

//	cout << VectorMax(Numbers) << endl;
}

void TestVectorMin()
{
  cout << "TestVectorMinDouble()" << endl << "-----------" << endl;
  std::vector<double> Numbers;
  Numbers.push_back(3.4);
  Numbers.push_back(4.5);
  Numbers.push_back(1.2);

  OutputVector(Numbers);
  cout << endl;

//	cout << VectorMin(Numbers) << endl;
}

void TestOutputVector()
{
  cout << "TestOutputVector()" << endl << "-----------" << endl;
  std::vector<double> Numbers;
  Numbers.push_back(3.4);
  Numbers.push_back(4.5);
  Numbers.push_back(1.2);

  OutputVector(Numbers);

}

void TestUniqueRandomIndices()
{
  std::vector<unsigned int> Indices = UniqueRandomIndices(20, 5);
  OutputVector(Indices);

}

void TestCombineVectors()
{
  cout << "TestCombineVectors()" << endl << "-----------" << endl;
  std::vector<double> A;
  A.push_back(3.4);
  A.push_back(4.5);
  A.push_back(1.2);

  cout << "A: " << endl;
  OutputVector(A);

  std::vector<double> B;
  B.push_back(8.1);
  B.push_back(9.2);
  B.push_back(11.3);

  cout << "B: " << endl;
  OutputVector(B);

  CombineVectors(A,B);
  cout << "Combined: " << endl;
  OutputVector(A);
}

void TestRemoveVectorElement()
{
  cout << "TestRemoveVectorElement()" << endl << "-----------" << endl;
  std::vector<double> A;
  A.push_back(3.4);
  A.push_back(4.5);
  A.push_back(1.2);
  A.push_back(7.3);
  A.push_back(10.4);
  A.push_back(9.5);

  OutputVector(A);
  cout << endl << endl;
  RemoveVectorElement(A,1);
  OutputVector(A);
}

void TestUniqueElements()
{
  cout << "TestRemoveVectorElement()" << endl << "-----------" << endl;
  std::vector<double> A;
  A.push_back(3.4);
  A.push_back(4.5);
  A.push_back(3.4);

  OutputVector(A);
  cout << endl << endl;
  std::vector<double> Unique = UniqueElements(A);
  OutputVector(Unique);
}
