#include <iostream>

#include "vtkPoints.h"
#include "vtkPlane.h"

#include <VTKHelpers/VTKHelpers.h>

bool CloseEnough(const double a, const double b);

int main()
{
	
	vtkPoints* Points = vtkPoints::New();	
	
	//create 3 random points
	/*
	for(unsigned int i = 0; i < 3; i++)
	{
		Points->InsertNextPoint(drand48(), drand48(), drand48());
	}
	*/
	
	//create 3 known points
	Points->InsertNextPoint(0.0, 0.0, 0.0);
	Points->InsertNextPoint(1.0, 0.0, 0.0);
	Points->InsertNextPoint(0.0, 1.0, 0.0);
	
	vtkPlane* BestPlane = VTKHelpers::BestFitPlane(Points);
	
	std::cout << "Best Plane:\n" << *BestPlane << "\n";
	
	double PlaneOrigin[3], PlaneNormal[3];
	
	BestPlane->GetNormal(PlaneNormal);
	BestPlane->GetOrigin(PlaneOrigin);
	
	if(CloseEnough(PlaneNormal[0], 0.0) && CloseEnough(PlaneNormal[1], 0.0) && CloseEnough(PlaneNormal[2], 1.0) && CloseEnough(PlaneOrigin[0], 0.33333) && CloseEnough(PlaneOrigin[1], .33333) && CloseEnough(PlaneOrigin[2], 0.0))
		return 0;
	else
	{
		std::cout << PlaneNormal[0] << " (should be 0.0)\n";
		std::cout << PlaneNormal[1] << " (should be 0.0)\n";
		std::cout << PlaneNormal[2] << " (should be 1.0)\n";
		std::cout << PlaneOrigin[0] << " (should be 0.33333)\n";
		std::cout << PlaneOrigin[1] << " (should be 0.33333)\n";
		std::cout << PlaneOrigin[2] << " (should be 0.0)\n";
		return -1;
	}
}

bool CloseEnough(const double a, const double b)
{
	double eps = 1e-3;
	if(fabs(a - b) < eps)
		return true;
	else
		return false;
}