#include "KDTree.h"

#include <ModelFile/ModelFile.h>

#include <Tools.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_distance.h>

#include <string>

#include <vul/vul_timer.h>

void TestNearestPoints(const ModelFile &Model);
void TestTimeDifference(const ModelFile &Model);
vector<vgl_point_3d<double> > TestAll1NearestNaive(const ModelFile &Model);
vector<vgl_point_3d<double> > TestAll1NearestKDTree(const ModelFile &Model, const KDTree &Tree);

vul_timer timer;

int main(int argc, char *argv[])
{
  Tools::AssertNumArgs(argc, 1);
  string InputFilename = argv[1];

  ModelFile Model(InputFilename);

  //TestNearestPoints(Model);
  TestTimeDifference(Model);

  return 0;
}


void TestNearestPoints(const ModelFile &Model)
{
  KDTree Tree(Model.getCoords());

  vgl_point_3d<double> TestPoint(.8, -1.3, -.45);
  vector<vgl_point_3d<double> > NearestPoints = Tree.KNearest(TestPoint, 5);

  Tools::OutputVector(NearestPoints);

}

void TestTimeDifference(const ModelFile &Model)
{
  timer.mark();
  vector<vgl_point_3d<double> > Naive = TestAll1NearestNaive(Model);
  cout << "Naive: " << timer.real() << endl;

  timer.mark();
  KDTree Tree(Model.getCoords());
  cout << "Create tree: " << timer.real() << endl;

  timer.mark();
  vector<vgl_point_3d<double> > KD = TestAll1NearestKDTree(Model, Tree);
  cout << "KDTree: " << timer.real() << endl;

  int NumSame = 0;
  for(unsigned int i = 0; i < Model.NumPoints(); i++)
  {
          if(Naive[i] == KD[i])
                  NumSame++;
  }

  cout << "NumPoints: " << Model.NumPoints() << endl;
  cout << "NumSame: " << NumSame << endl;
}

vector<vgl_point_3d<double> > TestAll1NearestNaive(const ModelFile &Model)
{
  vector<vgl_point_3d<double> > ClosestPoints(Model.NumPoints());

  for(unsigned int i = 0; i < Model.NumPoints(); i++)
  {
          double MinDist = 1e6;
          vgl_point_3d<double> ClosestPoint;

          for(unsigned int j = 0; j < Model.NumPoints(); j++)
          {
                  //if(i == j)
                  //	continue;

                  double dist = vgl_distance(Model.getPoint(i).getCoord(), Model.getPoint(j).getCoord());
                  if(dist < MinDist)
                  {
                          MinDist = dist;
                          ClosestPoint = Model.getPoint(i).getCoord();
                  }
          }
          ClosestPoints[i] = ClosestPoint;
  }

  return ClosestPoints;
}

vector<vgl_point_3d<double> > TestAll1NearestKDTree(const ModelFile &Model, const KDTree &Tree)
{
  vector<vgl_point_3d<double> > ClosestPoints(Model.NumPoints());
  for(unsigned int i = 0; i < Model.NumPoints(); i++)
  {
          vector<unsigned int> NearestPoints = Tree.KNearestIndices(Model.getCoord(i), 1);
          ClosestPoints[i] = Model.getCoord(NearestPoints[0]);
  }

  return ClosestPoints;
}
