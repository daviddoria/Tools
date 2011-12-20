#ifndef MODELFILE_H
#define MODELFILE_H

//STL
#include <vector>
#include <string>
#include <map>

//VGL
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_box_3d.h>

//VNL
#include <vnl/vnl_vector.h>
#include <vnl/vnl_double_3x3.h>

//VTK
#include <vtkPoints.h>
#include <vtkPolyData.h>

//Geometry
#include <Geometry/Triangle.h>
#include <Geometry/Color.h>
#include <Geometry/Angles.h>
#include <Geometry/Transformation.h>
#include <Geometry/Geometry.h>
#include <Geometry/OrientedPoint.h>
#include <Geometry/Helpers.h>

#include <VTKHelpers/Octree.h>
#include <KDTree/KDTree.h>

#include <Tools.h>

double PhiAngle(const vgl_point_3d<double> &P, std::string CoordSystem);
double ThetaAngle(const vgl_point_3d<double> &P, std::string CoordSystem);

class ModelFile
{
	private:
		std::map<std::string, std::vector<double> > DoubleArrays;
		
	public:
		//double epsilon;
		bool Valid_;
		bool ValidColors_;
		std::string Name_; //for sanity checking/tracking purposes
		bool Silent_; // determines if the readers/writers will write information to cout
		double SampleSpacing_;
		
		bool IsScan_;
		vgl_point_3d<double> ScannerLocation_;
		
		//Triangles
		std::vector<std::vector<unsigned int> > VertexLists_;
		std::vector<Triangle> Triangles_;
		
		//Lines
		//vector<pair<vgl_point_3d<double>, vgl_point_3d<double> > > Lines_;
		//std::vector<std::pair<unsigned int, unsigned int> >  Lines_;
		
		//std::vector<double> Scores_;
		
		/////////// Functions //////////////
		
		void UpdateTriangles(void);

	public:
		KDTree KDTree_;
		Octree Octree_;
		
		bool KDTreeBuilt_;
		
		/////////// Functions //////////////
		void BuildOctree(void);
		void BuildKDTree(void);
		
		///////////// Public Variables /////////////
		std::vector<OrientedPoint> Points_;

		///////// Constructors ////////////
		ModelFile() 
		{
			KDTreeBuilt_ = false;
			Valid_ = false;
			ValidColors_ = false;
			IsScan_ = false;
			Name_ = "UNNAMED";
		}

		/*
		ModelFile(const std::string &Filename, bool Silent = false);//read from file
		ModelFile(const std::vector<OrientedPoint> &Points);
		*/
		ModelFile(vtkSmartPointer<vtkPolyData> polydata)
		{
			ConstructFromPolydata(polydata);
		}
		
		void Init();
		
		void ConstructFromPolydata(vtkSmartPointer<vtkPolyData> polydata);
		
		////////// Functions ////////////
		vgl_point_3d<double> getCenterOfMass() const;
		void Translate(const vgl_vector_3d<double>  &T);
		void Translate(const double x, const double y, const double z) {vgl_vector_3d<double>  V(x,y,z); Translate(V);}
		void setPosition(const vgl_point_3d<double>  &Center);
		
		void Rotate(const vnl_double_3x3 &R);
		void RotateAroundCenter(const vnl_double_3x3 &R);
		
		void Transform(const Transformation &Trans);

		vgl_box_3d<double> GetBoundingBox() const {return geom::BoundingBox(GetOPCoords(Points_));}
		void WriteBoundingBox(const std::string &Filename);
		double AverageBBEdge() const;
		void PutOnGround();

		double AverageDistance(const vgl_point_3d<double> &P) const;
				
		void DisplayPoints();
				
		
		
		bool IntersectRay(const Ray &R, OrientedPoint &Intersection) const;

		///////// Mutators ///////////
		void setName(const std::string &Name) {Name_ = Name;}
		void setSilent(const bool silent) {Silent_ = silent;}
		void setCoords(const std::vector<vgl_point_3d<double> > &P);
		void setPoints(const std::vector<OrientedPoint> &P);
		void setVertexLists(const std::vector<vector<unsigned int> > &VertexLists);
		void setColors(std::vector<Color<unsigned char> > &C);
		//void setScores(std::vector<double> &Scores);
		void setNormals(vector<vgl_vector_3d<double> > &N);
		void setPointIndices();
		void setScannerLocation(const vgl_point_3d<double> &loc);
		
		void AddDoubleArray(const std::string &ArrayName, const std::vector<double> &V)
		{
			DoubleArrays[ArrayName] = V;
		}
		
		void DeleteTriangles()
		{
			VertexLists_.clear();
			Triangles_.clear();
		}
		
		//void setLines(const vector<pair<unsigned int, unsigned int > > &Lines) { Lines_ = Lines;}
		
		void SetCenterOfMass(const vgl_point_3d<double> &Center);
		void SetCenterOfMassX(const double X) {SetCenterOfMass(vgl_point_3d<double> (X, 0, 0));}
		void SetCenterOfMassY(const double Y) {SetCenterOfMass(vgl_point_3d<double> (0, Y, 0));}
		void SetCenterOfMassZ(const double Z) {SetCenterOfMass(vgl_point_3d<double> (0, 0, Z));}

		/////////////////////// Accessors /////////////////////
		bool IsValid() const {return Valid_;}
		std::string getName() const {return Name_;}
		double getSampleSpacing() const {return SampleSpacing_;}

		bool getScannerLocation(vgl_point_3d<double> &loc) const
		{
			loc = ScannerLocation_;
			return IsScan_;
		}
		
		//////////// Full Data Accessors //////////
		std::vector<OrientedPoint> getPoints() const {return Points_;}
		std::vector<vgl_point_3d<double> > getCoords() const {return GetOPCoords(Points_);}
		std::vector<std::vector<unsigned int > > getVertexLists() const {return VertexLists_;}
		std::vector<unsigned int> getVertexList(const unsigned int i) const {return VertexLists_[i];}
		std::vector<Color<unsigned char> > getColors() const {	return GetOPColors(Points_);}
		Octree getOctree() const {return Octree_;}
		
		std::vector<vgl_vector_3d<double> > getNormals(void) const 
		{
			std::vector<vgl_vector_3d<double> > Normals(NumPoints());
			for(unsigned int i = 0; i < NumPoints(); i++)
			{
				Normals[i] = Points_[i].getNormal();
			}
			return Normals;
		}
		
		/////////// Individual Data Accessors ////////////////
		OrientedPoint getPoint(const unsigned int i) const { assert(i < NumPoints()); return Points_[i];}
		vgl_point_3d<double> getCoord(const unsigned int i) const
		{
			
			if((i >= NumPoints()))
			{
				std::cout << "Requested point " << i << " of " << NumPoints() << std::endl;
				assert(0);
			}
			
			return Points_[i].getCoord();
		}
		
		vgl_vector_3d<double> getPointNormal(const unsigned int i) const { assert(i < NumPoints());return Points_[i].getNormal();}
		std::vector<unsigned int> getTriangleVertices(const unsigned int n) const {assert(n < VertexLists_.size()); return VertexLists_[n];}
		Color<unsigned char> getColor(const unsigned int i) const {	assert(i < Points_.size());	return Points_[i].getColor();}		
		Triangle getTriangle(const unsigned int n) const  {assert(n < Triangles_.size()); return Triangles_[n]; }
		double getTotalInformation(const unsigned int n) const  
		{
			assert(n < Points_.size()); 
			//return Points_[n].getTotalInformation(); 
			return Points_[n].getDoubleValue("TotalInformation");
		}
		double getRemainingInformation(const unsigned int n) const  
		{
			assert(n < Points_.size()); 
			return Points_[n].getDoubleValue("RemainingInformation");
			//return Points_[n].getRemainingInformation(); 
		}
		
		/////////// Length of Data Accessors /////////////
		unsigned int NumVertices() const {return Points_.size();}
		unsigned int NumPoints() const {return Points_.size();} // alias NumVertices
		unsigned int NumTriangles() const {return VertexLists_.size();}
		unsigned int NumColors() const {return Points_.size(); }//the colors are contained in the OrientedPoints
		
		///////// Text Displays ////////////
		void DisplayVertices();
		void DisplayTriangles();


		////////// Read/Writers ////////////
		bool Read(const std::string &filename);
		bool Write(const std::string &filename) const;
				
		bool ObjRead(const std::string &filename);
		bool ObjWrite(const std::string &filename) const;
				
		bool PlyRead(const std::string &filename);
		bool PlyWrite(const std::string &filename) const;
		
		bool PtxRead(const std::string &filename);
		bool PtxWrite(const std::string &filename) const;

		vtkSmartPointer<vtkPolyData> CreatePolydata(void) const;
		bool VtpRead(const std::string &filename);
		bool VtpWrite(const std::string &filename) const;
		bool SimpleVtpWrite(const std::string &Filename) const;
				
		bool VtuRead(const std::string &filename);
		bool VtuWrite(const std::string &filename) const;
		
		bool TxtCartesianRead(const std::string &Filename); // x y z
		bool TxtWrite(const std::string &filename) const;

		bool TxtSphericalRead(const std::string &Filename); // phi theta radius
		bool TxtSphericalWrite(const std::string &filename) const;

		////////// Functions ////////////
		vtkSmartPointer<vtkPolyData> GetVTKPoints() const;
		void ColorByValues(const std::vector<double> &Values);
		void RemoveGround(double angle);
		void ComputeNormals(double radius);
		void ProjectOntoXY(void);
		void Downsample(const double MinDist);
		void SnapToPoints(const std::vector<vgl_point_3d<double> > &Points);
		
		double MeshResolution();
		double GetAverageEdgeLength();
		void FlipNormal(unsigned int index, std::vector<unsigned int> &Used);
		void ColorRed();
		void ColorWhite();
		
		bool AllNormalsValid();
		
		bool GetDoubleArray(const std::string &ArrayName, std::vector<double> &Array);
};

void AddNormalsToPolyData(vtkSmartPointer<vtkPolyData> &polydata);
std::vector<vgl_vector_3d<double> > GetNormalsFromPolydata(vtkSmartPointer<vtkPolyData> polydata);

std::vector<double> GetDoubleArrayFromPolyData(const std::string &InputFilename, const std::string &ArrayName);
std::vector<float> GetFloatArrayFromPolyData(const std::string &InputFilename, const std::string &ArrayName);

std::ostream & operator<< (std::ostream &output, const ModelFile &Model);

#endif
