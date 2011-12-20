#ifndef VTKHELPERS_H
#define VTKHELPERS_H

#include <iostream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkPlane.h>
#include <vtkPoints.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_plane_3d.h>

#include <Geometry/Transformation.h>

namespace VTKHelpers
{
    bool FarthestTrianglePoint(double* Point, vtkTriangle* Triangle);
    
	vtkPlane* BestFitPlane(vtkPoints *pts);
	//Transformation VTKMatrixToTransformation(const vtkMatrix4x4* M);
	Transformation VTKMatrixToTransformation(const vtkMatrix4x4* M);
	
	void WritePlaneVTP(const vgl_plane_3d<double> &Plane, const vgl_point_3d<double> &PlaneCenter, const std::string &Filename, const bool Big);
	
	template <typename T, typename vtktype>
	void AddArrayToPolydata(const std::vector<T> &V, const std::string &Name, vtkSmartPointer<vtkPolyData> &polydata)
	{
		vtkIdType idNumPoints = polydata->GetNumberOfPoints();
		unsigned int NumPoints = static_cast<unsigned int>(idNumPoints);
		
		if(V.size() != NumPoints)
		{
			std::cout << "Cannot add vector " << Name << ". Size of V: " << V.size() << ", size of polydata points " << NumPoints << std::endl;
			exit(-1);
		}
			
		vtkSmartPointer<vtktype> Array = vtkSmartPointer<vtktype>::New();
		Array->SetNumberOfComponents(1);
		Array->SetName(Name.c_str());
	
		for(unsigned int i = 0; i < NumPoints; i++)
		{
			Array->InsertNextValue(V[i]);
		}
	
		polydata->GetPointData()->AddArray(Array);
	}
	
	
	template <typename T, typename vtktype>
	void GetArrayFromPolydata(const std::vector<T> &V, const std::string &Name, vtkSmartPointer<vtkPolyData> &polydata)
	{
		vtkIdType NumPoints = polydata->GetNumberOfPoints();
		
		vtktype* Array = vtktype::SafeDownCast(polydata->GetPointData()->GetArray(Name.c_str()));
	
		V.clear();
		if(Array)
		{
			for(unsigned int i = 0; i < static_cast<unsigned int>(NumPoints); i++)
			{
				V.push_back(Array->GetValue(i));
			}
			
		}
		
		//returned by reference
	}
	
	void WritePolyData(vtkSmartPointer<vtkPolyData> &polydata, const std::string &Filename);
}//end namespace
#endif