#include "VTKHelpers.h"

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_3x3.h>

#include <vtkSmartPointer.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>

#include <Geometry/Geometry.h>

namespace VTKHelpers
{
  bool FarthestTrianglePoint(double* Point, vtkTriangle* Triangle)
  {
  
  }
  
	vtkPlane* BestFitPlane(vtkPoints *points)
	{
		vtkIdType NumPoints = points->GetNumberOfPoints();
		
		//find the center of mass of the points
		double Center[3] = {0.0, 0.0, 0.0};
		
		for(vtkIdType i = 0; i < NumPoints; i++)
		{
			double point[3];
			points->GetPoint(i, point);
		
			Center[0] += point[0];
			Center[1] += point[1];
			Center[2] += point[2];
		}
				
		Center[0] = Center[0]/static_cast<double>(NumPoints);
		Center[1] = Center[1]/static_cast<double>(NumPoints);
		Center[2] = Center[2]/static_cast<double>(NumPoints);
		
		//Compute Sample Covariance Matrix (with points centered at the origin)
			
		double *a[3], a0[3], a1[3], a2[3];
		a[0] = a0; a[1] = a1; a[2] = a2; 
		for(unsigned int i = 0; i < 3; i++)
		{
			a0[i] = a1[i] = a2[i] = 0.0;
		}
		
		for(unsigned int pointId = 0; pointId < NumPoints; pointId++ )
		{
			double x[3], xp[3];
			points->GetPoint(pointId, x);
			xp[0] = x[0] - Center[0]; 
			xp[1] = x[1] - Center[1]; 
			xp[2] = x[2] - Center[2];
			for (unsigned int i = 0; i < 3; i++)
			{
				a0[i] += xp[0] * xp[i];
				a1[i] += xp[1] * xp[i];
				a2[i] += xp[2] * xp[i];
			}
		}

		for(unsigned int i = 0; i < 3; i++)
		{
			a0[i] /= static_cast<double>(NumPoints);
			a1[i] /= static_cast<double>(NumPoints);
			a2[i] /= static_cast<double>(NumPoints);
		}
		
		// Extract eigenvectors from covariance matrix
		double *v[3], v0[3], v1[3], v2[3];
		v[0] = v0; v[1] = v1; v[2] = v2; 
		double eigval[3];
		vtkMath::Jacobi(a,eigval,v);
		//Jacobi iteration for the solution of eigenvectors/eigenvalues of a 3x3 real symmetric matrix. Square 3x3 matrix a; output eigenvalues in w; and output eigenvectors in v. Resulting eigenvalues/vectors are sorted in decreasing order; eigenvectors are normalized.
		
		vtkPlane* BestPlane = vtkPlane::New();
		//Set the plane normal to the smallest eigen vector
		BestPlane->SetNormal(v2[0], v2[1], v2[2]);
		
		//Set the plane origin to the center of mass
		BestPlane->SetOrigin(Center[0], Center[1], Center[2]);
		
		return BestPlane;
	}
	
	void WritePlaneVTP(const vgl_plane_3d<double> &Plane, const vgl_point_3d<double> &PlaneCenter, const std::string &Filename, const bool Big)
	{
	
		vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
		if(Big)
		{
			vgl_vector_3d<double> V1 = geom::GetOrthogonalVector(Plane.normal());
			normalize(V1);
			vgl_vector_3d<double> V2 = cross_product(Plane.normal(), V1);
			normalize(V2);
			vgl_point_3d<double> P1 = PlaneCenter + 50. * V1;
			vgl_point_3d<double> P2 = PlaneCenter + 50. * V2;
			plane->SetPoint1(P1.x(), P1.y(), P1.z());
			plane->SetPoint2(P2.x(), P2.y(), P2.z());
		}
		
		plane->SetCenter(PlaneCenter.x(), PlaneCenter.y(), PlaneCenter.z());
		
		plane->SetNormal(Plane.normal().x(), Plane.normal().y(), Plane.normal().z());
			
		vtkSmartPointer<vtkPolyData> pdata = plane->GetOutput();
			
		//write the file
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetInput(pdata);
		writer->SetFileName(Filename.c_str());
		writer->Write();
	}
	
	//Transformation VTKMatrixToTransformation(const vtkMatrix4x4* M)
	Transformation VTKMatrixToTransformation(const vtkMatrix4x4* M)
	{
		
		vnl_double_3 T;
		vnl_double_3x3 R;
	
		//extract R
		for(unsigned int r = 0; r < 3; r++)
		{
			for(unsigned int c = 0; c < 3; c++)
			{
				R(r,c) = M->GetElement(r,c);
				//R(r,c) = icp->GetMatrix()->GetElement(c,r);
			}
		}
	
		//extract T
		T(0) = M->GetElement(0,3);
		T(1) = M->GetElement(1,3);
		T(2) = M->GetElement(2,3);
	
		Transformation Trans(T,R);
		
		return Trans;
		
	}
	
	void WritePolyData(vtkSmartPointer<vtkPolyData> &polydata, const std::string &Filename)
	{
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetFileName(Filename.c_str());
		writer->SetInput(polydata);
		writer->Write();
	}
}//end namespace