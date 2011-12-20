#include "Geometry.h"

#include <vgl/algo/vgl_fit_plane_3d.h>
#include <vgl/vgl_homg_plane_3d.h>
#include <vgl/vgl_line_3d_2_points.h>
#include <vgl/vgl_distance.h>
#include <vgl/vgl_sphere_3d.h>

#include <vnl/vnl_inverse.h>

#include <VXLHelpers/VXLHelpers.h>

#include "Angles.h"
#include "Color.h"

#include <Tools.h>

#include <algorithm>

#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkLandmarkTransform.h>
#include <vtkMatrix4x4.h>

#include <limits>

#include <KDTree/KDTree.h>

namespace geom
{

	vnl_double_3x3 SampleCovarianceMatrix(const std::vector<vgl_point_3d<double> > &Points)
	{
		vgl_point_3d<double> u = CenterOfMass(Points);
		vnl_double_3x3 S;
		S.fill(0);
		
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			vgl_vector_3d<double> V = Points[i] - u;
			vnl_double_3x3 OP = VXLHelpers::OuterProduct(V, V);
			S += OP;
		}
		
		S /= Points.size();
		return S;
	}
	
	double AveragePlaneError(const std::vector<vgl_point_3d<double> > &Points, const vgl_plane_3d<double> &Plane)
	{
		double Total = 0;
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			Total += vgl_distance(Plane, Points[i]);
		}
		
		return Total / Points.size();
	}
	
	vgl_plane_3d<double> FitPlane(const std::vector<vgl_point_3d<double> > &Points)
	{
		//this function assumes the points actually fit a plane (within some tolerance)!
		vgl_fit_plane_3d<double> P;
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			P.add_point(Points[i].x(), Points[i].y(), Points[i].z());
		}
	
		//P.fit(1e-6);
		P.fit(1.0); //allowable error
		
		return vgl_plane_3d<double> (P.get_plane());
	}
	
	vgl_plane_3d<double> BestPlane(const std::vector<vgl_point_3d<double> > &Points)
	{
		//this function finds the eigen vector corresponding to the smallest eigen value and places the plane's P0 at the center of mass of the points.
		
		vgl_point_3d<double> Centroid = CenterOfMass(Points);
		vnl_double_3x3 S = SampleCovarianceMatrix(Points);
		std::vector<vnl_double_3> evecs = VXLHelpers::EigenVectors(S);
		
		std::vector<vgl_vector_3d<double> > EvecsVGL;
		for(unsigned int i = 0; i < 3; i++)
			EvecsVGL.push_back(VXLHelpers::vnl_vector_to_vgl_vector(evecs[i]));
		
		//WriteEigenVectorVTP(Centroid, EvecsVGL, .3);
		
		vgl_plane_3d<double> P(EvecsVGL[0], Centroid); //apparently these are sorted from smallest to largest?
		
		return P;
	}
	
	void WriteAxisFile(const vgl_point_3d<double> &Origin, const vgl_vector_3d<double> &V0, const vgl_vector_3d<double> &V1, const vgl_vector_3d<double> &V2, const std::string &OutputFilename)
	{
		unsigned char r[3], g[3], b[3];
		CharArray(Colors::Red(), r);
		CharArray(Colors::Green(), g);
		CharArray(Colors::Blue(), b);
	
		double XaxisPoint[3], YaxisPoint[3], ZaxisPoint[3];
		VXLHelpers::GetArray(VXLHelpers::PointAlong(Origin, V0, 1.0), XaxisPoint);
		VXLHelpers::GetArray(VXLHelpers::PointAlong(Origin, V1, 1.0), YaxisPoint);
		VXLHelpers::GetArray(VXLHelpers::PointAlong(Origin, V2, 1.0),  ZaxisPoint);
	
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	
		//setup colors
		vtkSmartPointer<vtkUnsignedCharArray> Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		Colors->SetNumberOfComponents(3);
		Colors->SetName("Colors");
		
		Colors->InsertNextTupleValue(r);
		Colors->InsertNextTupleValue(g);
		Colors->InsertNextTupleValue(b);
	
		double o[3] = {0.0, 0.0, 0.0};
		pts->InsertNextPoint(o);
		pts->InsertNextPoint(XaxisPoint);
		pts->InsertNextPoint(YaxisPoint);
		pts->InsertNextPoint(ZaxisPoint);
		
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
		//x axis
		vtkSmartPointer<vtkLine> xaxis = vtkSmartPointer<vtkLine>::New();
		xaxis->GetPointIds()->SetId(0,0);
		xaxis->GetPointIds()->SetId(1,1);
		lines->InsertNextCell(xaxis);
	
		//y axis
		vtkSmartPointer<vtkLine> yaxis = vtkSmartPointer<vtkLine>::New();
		yaxis->GetPointIds()->SetId(0,0);
		yaxis->GetPointIds()->SetId(1,2);
		lines->InsertNextCell(yaxis);
	
		//z axis
		vtkSmartPointer<vtkLine> zaxis = vtkSmartPointer<vtkLine>::New();
		zaxis->GetPointIds()->SetId(0,0);
		zaxis->GetPointIds()->SetId(1,3);
		lines->InsertNextCell(zaxis);
	
		vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
	
		//add the points to the dataset
		pdata->SetPoints(pts);
	
		//add the lines to the dataset
		pdata->SetLines(lines);
	
		//color the lines
		pdata->GetCellData()->AddArray(Colors);
	
		//write the file
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetInput(pdata);
	
		writer->SetFileName(OutputFilename.c_str());
		writer->Write();
	
	}
	
	//vector<vgl_point_2d<double> > PointsOnPlaneTo2D(const vector<vgl_point_3d<double> > &Points)
	
	std::vector<vgl_point_2d<double> > ScannerImage(const vector<vgl_point_3d<double> > &Points, const vgl_point_3d<double> &ScannerLocation)
	{
		std::vector<vgl_point_3d<double> > Intersections(Points.size());
		
		//get all directions
		std::vector<vgl_vector_3d<double> > Directions(Points.size());
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			Directions[i] = Points[i] - ScannerLocation;
		}
		
		vgl_vector_3d<double> ProjectionPlaneNormal = VXLHelpers::normalize(AverageDirection(Directions));
			
		vgl_point_3d<double> ProjectionPlanePoint = ScannerLocation + ProjectionPlaneNormal;
		
		//find the 3d coords of the points projected on the plane
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			vgl_line_3d_2_points<double> Line(ScannerLocation, Points[i]);
	
			//vgl_plane_3d<double> Plane(ScannerForward, ScannerLocation + ScannerForward);
			vgl_plane_3d<double> Plane(ProjectionPlaneNormal, ProjectionPlanePoint);
	
			vgl_point_3d<double> Intersection = vgl_intersection(Line, Plane);
	
			Intersections[i] = Intersection;
		}
		
		vgl_vector_3d<double> b = VXLHelpers::normalize(Intersections[0] - ProjectionPlanePoint);
		
		//WriteAxisFile(vgl_point_3d<double> (0,0,0), cross_product(ProjectionPlaneNormal, b), b, ProjectionPlaneNormal, "BeforeAxis.vtp");
		
		vnl_double_3x3 R = FindAligningRotation(cross_product(ProjectionPlaneNormal, b), b, ProjectionPlaneNormal);
		
		/*
		//write out axes after transformation
		vnl_double_3 a1 = vnl_inverse(R) * vgl_vector_to_vnl_vector(cross_product(ProjectionPlaneNormal, b));
		vnl_double_3 a2 = vnl_inverse(R) * vgl_vector_to_vnl_vector(b);
		vnl_double_3 a3 = vnl_inverse(R) * vgl_vector_to_vnl_vector(ProjectionPlaneNormal);
		WriteAxisFile(vgl_point_3d<double> (0,0,0), vnl_vector_to_vgl_vector(a1), vnl_vector_to_vgl_vector(a2), vnl_vector_to_vgl_vector(a3), "AfterAxis.vtp");
		*/
		
		std::vector<vgl_point_2d<double> > Points2d(Points.size());
		
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			vnl_double_3 temp = vnl_inverse(R) * VXLHelpers::vgl_point_to_vnl_vector(Intersections[i]);
			Points2d[i] = vgl_point_2d<double> (temp(0), temp(1));
		}
		
		return Points2d;
	}
	
	vgl_vector_3d<double> AverageDirection(const std::vector<vgl_vector_3d<double> > &Vectors)
	{
		vgl_vector_3d<double> Average(0.0, 0.0, 0.0);
		for(unsigned int i = 0; i < Vectors.size(); i++)
		{
			Average += Vectors[i];
		}
		
		Average /= Vectors.size();
		
		return Average;
	}
	
	vnl_double_3x3 FindAligningRotation(const vgl_vector_3d<double> &x, const vgl_vector_3d<double> &y, const vgl_vector_3d<double> &z)
	{
		//create a matrix that does what ???
		vnl_double_3x3 R;
		R.set_column(0, VXLHelpers::vgl_vector_to_vnl_vector(x));
		R.set_column(1, VXLHelpers::vgl_vector_to_vnl_vector(y));
		R.set_column(2, VXLHelpers::vgl_vector_to_vnl_vector(z));
		
		return R;
	}
	
	vgl_vector_3d<double> Project(const vgl_vector_3d<double> &A, const vgl_vector_3d<double> &B)
	{
		//the projection of A on B
		//return A.Dot3(B)/pow(B.Length(),2) * B;
		return dot_product(A,B)/pow(B.length(), 2) * B;
	}
	
	vgl_vector_3d<double> Project(const vgl_vector_3d<double> &A, const vgl_plane_3d<double> &P)
	{
		//the projection of A on P
		//Project A onto P's normal. Subtract this from A to get the component of A lying in the plane P.
		vgl_vector_3d<double> N = P.normal();
		vgl_vector_3d<double> An = Project(A, N);
		vgl_vector_3d<double> Ap = A - An;
		normalize(Ap);
		return Ap;
	}
	
	vgl_vector_3d<double> GetRandomVec(void)
	{
		vgl_vector_3d<double> randvec(drand48(), drand48(), drand48());
		normalize(randvec);
		return randvec;
	}
	
	vgl_vector_3d<double> GetOrthogonalVector(const vgl_vector_3d<double> &V)
	{
		//Gram Schmidt Orthogonalization
		vgl_vector_3d<double> RandVec = GetRandomVec();
	
		vgl_vector_3d<double> N = RandVec - Project(RandVec, V);
		normalize(N);
		return N;
	}
	
	
	std::vector<vgl_point_3d<double> > PointsInsideFrustrum(const vgl_box_3d<double> &BoundingBox, const std::vector<vgl_point_3d<double> > &WorldPoints, const vgl_point_3d<double> &ScannerLocation)
	{
		vector<vgl_point_3d<double> > Points;
		
		for(unsigned int i = 0; i < WorldPoints.size(); i++)
		{
			Ray R(ScannerLocation, WorldPoints[i] - ScannerLocation);
			
			if(R.IntersectBox(BoundingBox))
				Points.push_back(WorldPoints[i]);
		}
		
		return Points;
	}
	
	std::vector<unsigned int> IndicesInsideFrustrum(const vgl_box_3d<double> &BoundingBox, const std::vector<vgl_point_3d<double> > &WorldPoints, const vgl_point_3d<double> &ScannerLocation)
	{
		vector<unsigned int> Indices;
		
		for(unsigned int i = 0; i < WorldPoints.size(); i++)
		{
			Ray R(ScannerLocation, WorldPoints[i] - ScannerLocation);
			
			if(R.IntersectBox(BoundingBox))
				Indices.push_back(i);
		}
		
		return Indices;
	}
	
	
	vgl_point_3d<double> ClosestPoint(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points)
	{
		double MinDistance = 1e6;
	
		vgl_point_3d<double> closest;
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			double dist = vgl_distance(Points[i], TestPoint);
			if(dist < MinDistance)
			{
				MinDistance = dist;
				closest = Points[i];
			}
					
		}
	
		return closest;
	}
	
	unsigned int ClosestPointIndex(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points)
	{
		double MinDistance = 1e6;
	
		unsigned int ClosestIndex = 0;
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			double dist = vgl_distance(Points[i], TestPoint);
			if(dist < MinDistance)
			{
				MinDistance = dist;
				ClosestIndex = i;
			}
		}
	
		return ClosestIndex;
	}
	
	double ClosestPointDist(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points)
	{
		double MinDistance = 1e6;
			
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			double dist = vgl_distance(Points[i], TestPoint);
			if(dist < MinDistance)
			{
				MinDistance = dist;
			}
		}
	
		return MinDistance;
	}
	

	
	vgl_point_3d<double> FarthestPoint(const vgl_point_3d<double> &TestPoint, const std::vector<vgl_point_3d<double> > &Points)
	{
		double MaxDistance = 0;
	
		vgl_point_3d<double> farthest;
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			double dist = vgl_distance(Points[i], TestPoint);
			if(dist > MaxDistance)
			{
				MaxDistance = dist;
				farthest = Points[i];
			}
					
		}
	
		return farthest;
	}
	
	
	vgl_point_3d<double> CenterOfMass(const std::vector<vgl_point_3d<double> > &Points)
	{
		double x = 0.0, y = 0.0, z = 0.0;
		unsigned int n = Points.size();
		for(unsigned int i = 0; i < n; i++)
		{
			x += Points[i].x();
			y += Points[i].y();
			z += Points[i].z();
		}
			
		return vgl_point_3d<double> (x/static_cast<double>(n), y/static_cast<double>(n), z/static_cast<double>(n));
			
	}
	
	vgl_box_3d<double> BoundingBox(const std::vector<vgl_point_3d<double> > &Points)
	{
		vgl_box_3d<double> BB;
		std::vector<vgl_point_3d<double>  > P;
		for(unsigned int i = 0; i < Points.size(); i++)
			P.push_back(Points[i]);
			
		vgl_box_3d_bounds(P.begin(), P.end(), BB);
	
		return BB;
	}
	
	std::vector<vgl_point_3d<double> > GetPointsInsideSphere(const std::vector<vgl_point_3d<double> > &Points, const vgl_sphere_3d<double> &Sphere)
	{
		std::vector<vgl_point_3d<double> > PointsInSphere;
		
		for(unsigned int j = 0; j < Points.size(); j++)
		{
			vgl_point_3d<double> CheckPoint = Points[j];
			
			if(Sphere.contains(CheckPoint))
			{
				PointsInSphere.push_back(CheckPoint);
			}
		}
		
		return PointsInSphere;
	}
	
	
	void WritePlaneVTP(const vgl_point_3d<double> &Point, const vgl_plane_3d<double> &Plane, const double PlaneSize)
	{
		vgl_vector_3d<double> V1 = GetOrthogonalVector(Plane.normal());
		vgl_vector_3d<double> V2 = VXLHelpers::cross(V1, Plane.normal());
		normalize(V2);
		
		std::vector<vgl_point_3d<double> > Points;
		//Triangle 1
		Points.push_back(Point + V1 * PlaneSize);
		Points.push_back(Point - V1 * PlaneSize);
		Points.push_back(Point + V2 * PlaneSize);
		
		//Triangle 2
		Points.push_back(Point - V2 * PlaneSize);
		
		//center
		Points.push_back(Point); //4
		
		//other point for normal
		Points.push_back(Point + Plane.normal()); //5
		
		//add the normal line
		std::vector<pair<unsigned int, unsigned int> > Lines;
		Lines.push_back(pair<unsigned int, unsigned int> (4,5));
		
		//add 2 triangles
		std::vector<vector<unsigned int> > VertexList;
		std::vector<unsigned int> Tri1;
		Tri1.push_back(0);
		Tri1.push_back(1);
		Tri1.push_back(2);
		
		vector<unsigned int> Tri2;
		Tri2.push_back(0);
		Tri2.push_back(1);
		Tri2.push_back(3);
		
		VertexList.push_back(Tri1);
		VertexList.push_back(Tri2);
			
		/*
		ModelFile Model;
		Model.setCoords(Points);
		Model.setVertexLists(VertexList);
		Model.setLines(Lines);
		Model.Write("Plane.vtp");
		*/
	}
	
	
	void WriteEigenVectorVTP(const vgl_point_3d<double> &Point, const std::vector<vgl_vector_3d<double> > &Evecs, const double Length)
	{
		
		std::vector<vgl_point_3d<double> > Points;
		
		Points.push_back(Point);
		Points.push_back(Point + Evecs[0]*Length);
		Points.push_back(Point + Evecs[1]*Length);
		Points.push_back(Point + Evecs[2]*Length);
		
		//add the 3 lines
		vector<pair<unsigned int, unsigned int> > Lines;
		Lines.push_back(pair<unsigned int, unsigned int> (0,1));
		Lines.push_back(pair<unsigned int, unsigned int> (0,2));
		Lines.push_back(pair<unsigned int, unsigned int> (0,3));
		/*
		ModelFile Model;
		Model.setCoords(Points);
		Model.setLines(Lines);
		Model.Write("Evecs.vtp");
		*/
	}
	
	double SignedPlaneDistance(const vgl_plane_3d<double> &Plane, const vgl_point_3d<double> &Point)
	{
		double d = vgl_distance(Plane, Point);
		vgl_point_3d<double> ClosestPoint = vgl_closest_point(Plane, Point);
		
		vgl_vector_3d<double> V = Point - ClosestPoint;
		
		if(AngleBetween(Plane.normal(), V) < M_PI/2.0)
			return d;
		else
			return -d;
		
	}
	
	vgl_point_3d<double> GetFarthestPoint(const vector<vgl_point_3d<double> > &Points, const vgl_point_3d<double> &TestPoint)
	{
		double MaxDist = 0.0;
		vgl_point_3d<double> FarthestPoint;
		
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			double d = vgl_distance(Points[i], TestPoint);
			if(d > MaxDist)
			{
				MaxDist = d;
				FarthestPoint = Points[i];
			}
		}
		
		return FarthestPoint;
	}
	
	bool IntersectSphere(const vgl_line_3d_2_points<double> &line, const vgl_sphere_3d<double> &sphere)
	{
		vgl_point_3d<double> p1;
		vgl_point_3d<double> p2;
		return sphere.clip(line, p1, p2);
	}
	
	#if 0
	vgl_point_3d<double> IntersectCloud(const vector<vgl_point_3d<double> > &Points, const vgl_line_3d_2_points<double> &line, const vgl_point_3d<double> &ScannerPosition, bool &Valid)
	{
	
		//!!! Was used in Octree to find points inside a cylinder (fat ray)
		//gets the closest point in Points to ScannerPosition
		
		double MinDist = numeric_limits<double>::infinity();
	
		vgl_point_3d<double> cp;
	
		double DistToScanner, DistToLine;
		
		for(unsigned int point = 0; point < Points.size(); point++)
		{
			DistToLine = vgl_distance(line, Points[point]);
			
			/*
			if(DistToLine < Thresh_)
			{
				DistToScanner = vgl_distance(ScannerPosition, Points[point]);
			}
			else
			{
				continue;
			}
			*/
			
			if (DistToScanner < MinDist)
			{
				MinDist = DistToScanner;
				Valid = true;
				cp = Points[point];
			}
		}
	
		assert(MinDist >= 0);
	
		if(Valid)
		{
			return cp;
		}
	
		Valid = false;
		return ScannerPosition;
	}
	#endif
	
	double GetNeighborhoodSize(const std::vector<vgl_point_3d<double> > &Points, const double VolumeRatio)
	{
		assert(VolumeRatio < 1.0);
		assert(VolumeRatio > 0.0);
		
		vgl_box_3d<double> Box = BoundingBox(Points);
		
		//we want a sphere with X% volume of the box
		double vol = VolumeRatio * Box.volume();
		// V = (4/3) pi r^3 -> r = V ((3/4) (1/pi))^(1/3)
		
		double r = vol * pow((static_cast<double>(3)/static_cast<double>(4)) * (static_cast<double>(1)/M_PI), .33333);
		std::cout << "Neighborhood size: " << r << std::endl;
		return r;
	}
	
	double CorrespondenceError(const std::vector<vgl_point_3d<double> > &ModelPoints, const std::vector<vgl_point_3d<double> > &ScanPoints)
	{
		//this function finds the landmark transform error between two sets of corresponding points
		if(ModelPoints.size() != ScanPoints.size())
		{
			std::cout << "Number of model points does not equal number of scan points!" << std::endl;
			std::cout << "ModelPoints : " << ModelPoints.size() << std::endl;
			std::cout << "ScanPoints : " << ScanPoints.size() << std::endl;
			exit(-1);
		}

		unsigned int NumPoints = ModelPoints.size();
	
		vtkSmartPointer<vtkLandmarkTransform> LandmarkTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
		vtkSmartPointer<vtkPoints> SourcePoints = vtkSmartPointer<vtkPoints>::New(); //model points
		vtkSmartPointer<vtkPoints> TargetPoints = vtkSmartPointer<vtkPoints>::New(); //scene points

		for(unsigned int i = 0; i < NumPoints; i++)
		{
			vgl_point_3d<double> ModelVGLPoint = ModelPoints[i];
			//cout << "Model Point: " << ModelVGLPoint << endl;
			double ModelPoint[3] = {ModelVGLPoint.x(), ModelVGLPoint.y(), ModelVGLPoint.z()};
			SourcePoints->InsertNextPoint(ModelPoint);

			vgl_point_3d<double> SceneVGLPoint = ScanPoints[i];
			//cout << "Scene Point: " << SceneVGLPoint << endl;
			double ScenePoint[3] = {SceneVGLPoint.x(), SceneVGLPoint.y(), SceneVGLPoint.z()};
			TargetPoints->InsertNextPoint(ScenePoint);
		}

		LandmarkTransform->SetSourceLandmarks(SourcePoints);
		LandmarkTransform->SetTargetLandmarks(TargetPoints);
		LandmarkTransform->SetModeToRigidBody();
		LandmarkTransform->Update();
	
		vtkSmartPointer<vtkPoints> TransformedModelPoints = vtkSmartPointer<vtkPoints>::New();
		LandmarkTransform->TransformPoints(SourcePoints, TransformedModelPoints);

		//cout << "There are " << SourcePoints->GetNumberOfPoints() << " source points and " << TransformedModelPoints->GetNumberOfPoints() << " transformed points." << endl;
		vector<vgl_point_3d<double> > NewModelPoints(NumPoints);
		for(unsigned int i = 0; i < NumPoints; i++)
		{
			double origpoint[3];
			TransformedModelPoints->GetPoint(i, origpoint);
			//cout << "original point: " << vgl_point_3d<double>(origpoint) << endl;
		
			double transpoint[3];
			TransformedModelPoints->GetPoint(i, transpoint);
			NewModelPoints[i] = vgl_point_3d<double>(transpoint);
			//cout << "trans point: " << NewModelPoints[i] << endl;
		}
	
		if(!(NewModelPoints.size() == ScanPoints.size()))
		{
			cout << "The sets of points must be the same size!" << endl;
			exit(0);
		}

		double Total = 0.0;
		for(unsigned int i = 0; i < NumPoints; i++)
		{
			Total += vgl_distance(NewModelPoints[i], ScanPoints[i]);
		}

		double Average = Total / static_cast<double>(NumPoints);
		return Average;
	}
	
	
	double CorrespondenceError(const std::vector<std::pair<unsigned int, unsigned int> > &Correspondences, const std::vector<vgl_point_3d<double> > &ModelPoints, const std::vector<vgl_point_3d<double> > &ScanPoints)
	{
		//extracts corresponding points from two point sets, then finds the landmark transform error between the corresponding sets
		
		//<model, scene>
		unsigned int NumPoints = Correspondences.size();
		std::vector<vgl_point_3d<double> > CorrModel(NumPoints);
		std::vector<vgl_point_3d<double> > CorrScan(NumPoints);
		for(unsigned int i = 0; i < NumPoints; i++)
		{
			CorrModel[i] = ModelPoints[Correspondences[i].first];
			CorrScan[i] = ScanPoints[Correspondences[i].second];
		}
		
		//find the landmark transform error between two sets of corresponding points
		return CorrespondenceError(CorrModel, CorrScan);
	}
	
	vnl_matrix<double> CameraTransform(const Ray &View1, const Ray &View2)
	{
		//This function takes two rays (ie 2 sets of view points and view directions), and finds the matrix M between them (from V1 to V2)
		
		vgl_point_3d<double> A1 = View1.getOrigin();
		vgl_point_3d<double> A2 = View1.PointAlong(1.0);
		
		vgl_point_3d<double> B1 = View2.getOrigin();
		vgl_point_3d<double> B2 = View2.PointAlong(1.0);
		
		vtkSmartPointer<vtkLandmarkTransform> LandmarkTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
		vtkSmartPointer<vtkPoints> SourcePoints = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> TargetPoints = vtkSmartPointer<vtkPoints>::New();

		double A1array[3] = {A1.x(), A1.y(), A1.z()};
		SourcePoints->InsertNextPoint(A1array);
		
		double A2array[3] = {A2.x(), A2.y(), A2.z()};
		SourcePoints->InsertNextPoint(A2array);
		
		double B1array[3] = {B1.x(), B1.y(), B1.z()};
		TargetPoints->InsertNextPoint(B1array);
		
		double B2array[3] = {B2.x(), B2.y(), B2.z()};
		TargetPoints->InsertNextPoint(B2array);
		
		LandmarkTransform->SetSourceLandmarks(SourcePoints);
		LandmarkTransform->SetTargetLandmarks(TargetPoints);
		LandmarkTransform->SetModeToRigidBody();
		LandmarkTransform->Update();
	
		vnl_matrix<double> Trans(4,4);
		vtkMatrix4x4* M = LandmarkTransform->GetMatrix();
		
		//cout << M << endl;
		
		for(unsigned int r = 0; r < 4; r++)
		{
			for(unsigned int c = 0; c < 4; c++)
			{
				Trans(r,c) = M->GetElement(r,c);
			}
		}
		
		//cout << Trans << endl;
		
		return Trans;
		
	}
	
}//end namespace geom