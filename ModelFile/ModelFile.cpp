#include "ModelFile.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cstdlib> //for atof()
#include <limits>
#include <algorithm>
#include <cmath>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_sphere_3d.h>
#include <vgl/vgl_distance.h>

//VTK
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkCellData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkCubeSource.h>
#include <vtkOBJReader.h>
#include <vtkOBJExporter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkDelaunay2D.h>
#include <vtkCellLocator.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>

#include <VXLHelpers/VXLHelpers.h>
#include <VTKHelpers/VTKHelpers.h>

#include <Tools/Tools.h>

#include <Geometry/Geometry.h>
#include <Geometry/Edge.h>

#include <boost/progress.hpp>

/////////// Operators //////////////////
std::ostream & operator << (std::ostream &output, const ModelFile &Model)
{
	output << "Num Vertices: " << Model.NumVertices() << std::endl
			<< "Num Triangles: " << Model.NumTriangles() << std::endl;
	return output;
}

//////////////// Constructors //////////////////
/*
//read from file
ModelFile::ModelFile(const string &Filename, bool Silent) : Silent_(Silent)
{
	Read(Filename);
	Init();
}

ModelFile::ModelFile(const vector<OrientedPoint> &Points)
{
	Points_ = Points;
	Init();
}
*/

void ModelFile::BuildOctree(void)
{
	Octree_.Build(OPVectorToCoordVector(Points_), VertexLists_);
}

void ModelFile::BuildKDTree(void)
{
	KDTree_ = KDTree(getCoords());
	KDTreeBuilt_ = true;
}

void ModelFile::Init(void)
{
	if(Points_.size() == 0)
	{
		std::cout << "You must set points before calling Init().";
		exit(-1);
	}
	
	SampleSpacing_ = KDFuncs::AveragePointDistance(GetOPCoords(Points_));
	
	//if this is a giant file, we only want the points
	if(Points_.size() > 1e6)
	{
		return;
	}
	else
	{
		//create the triangles from the points
		UpdateTriangles();
	}
	
	Valid_ = true;
}

void ModelFile::WriteBoundingBox(const std::string &Filename)
{
	vgl_box_3d<double> BB = GetBoundingBox();
	double xlen = BB.width();
	double ylen = BB.height();
	double zlen = BB.depth();
	vgl_point_3d<double> centroid = BB.centroid();
	
	vtkSmartPointer<vtkCubeSource> Cube1 = vtkSmartPointer<vtkCubeSource>::New();
	Cube1->SetCenter(centroid.x(), centroid.y(), centroid.z());
	Cube1->SetXLength(xlen);
	Cube1->SetYLength(ylen);
	Cube1->SetZLength(zlen);
		
	vtkSmartPointer<vtkPolyData> polydata = Cube1->GetOutput();
		
	//write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInput(polydata);

	writer->SetFileName(Filename.c_str());
	writer->Write();
}

//////////////// Accessors //////////////////
double ModelFile::AverageBBEdge() const
{
	vgl_box_3d<double> BB = GetBoundingBox();
	double AvEdge = (BB.width() + BB.height() + BB.depth())/3.0;
	if(AvEdge == 0.0)
	{
		std::cout << "The bounding box is empty!" << std::endl;
		exit(-1);
	}
	return AvEdge;
}




vgl_point_3d<double> ModelFile::getCenterOfMass() const
{
	vector<vgl_point_3d<double> > P = getCoords();
	return geom::CenterOfMass(P);
}


//////////// Mutators /////////////
void ModelFile::setPoints(const vector<OrientedPoint> &P)
{
	Points_ = P;
}

void ModelFile::setCoords(const vector<vgl_point_3d<double> > &P)
{
	/*
	//can't do this because P may not be the same size (ie. we are trying to create a new model)
	assert(P.size() == Points_.size());

	for(unsigned int i = 0; i < P.size(); i++)
		Points_[i].setCoord(P[i]);
	
	UpdateTriangles();	
	*/
	
	
	Points_.clear();
	for(unsigned int i = 0; i < P.size(); i++)
		Points_.push_back(OrientedPoint(P[i]));
	
}

void ModelFile::setVertexLists(const vector<vector<unsigned int> > &VertexLists)
{
	VertexLists_ = VertexLists;
	UpdateTriangles();
}

void ModelFile::setNormals(vector<vgl_vector_3d<double> > &N)
{
	assert(N.size() == NumPoints());
	for(unsigned int i = 0; i < NumPoints(); i++)
		Points_[i].setNormal(N[i]);
}

void ModelFile::setColors(vector<Color<unsigned char> > &C) 
{
	if(C.size() != NumPoints())
	{
		clog << "ModelFile::setColors - C.size()=" << C.size() << ", but NumPoints=" << NumPoints() << endl;
		exit(1);
	}
	
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		Points_[i].setColor(C[i]);
	}
	
	ValidColors_ = true;
}

/*
void ModelFile::setScores(vector<double> &Scores)
{
	assert(Scores.size() == NumPoints());
	Scores_ = Scores;
}
*/

void ModelFile::setPointIndices(void)
{
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		Points_[i].setIndex(i);
	}
}

void ModelFile::setScannerLocation(const vgl_point_3d<double> &Loc)
{
	IsScan_ = true;
	ScannerLocation_ = Loc;
}

void ModelFile::SetCenterOfMass(const vgl_point_3d<double> &Center)
{
	vgl_vector_3d<double>  ReCenter = vgl_point_3d<double>(0,0,0) - getCenterOfMass();
	vgl_vector_3d<double>  Move = Center - vgl_point_3d<double>(0,0,0);
		
	//move back to origin and then to new center
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		Points_[i].setCoord(Points_[i].getCoord() + ReCenter);
		Points_[i].setCoord(Points_[i].getCoord() + Move);
	}
	
	UpdateTriangles();
}


double ModelFile::AverageDistance(const vgl_point_3d<double>  &P) const
{
	double TotalDistance = 0;
	unsigned int NumPoints = Points_.size();
	for(unsigned int i = 0; i < NumPoints; i++)
	{
		vgl_vector_3d<double> V = P - Points_[i].getCoord();
		TotalDistance += V.length(); 
	}
	
	return TotalDistance / static_cast<double>(NumPoints);
	
}

void ModelFile::Translate(const vgl_vector_3d<double> &t)
{
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		Points_[i].setCoord(Points_[i].getCoord() + t);
	}
	
//	PointTree_ = Octree(Points_);
	UpdateTriangles();
}

void ModelFile::Rotate(const vnl_double_3x3 &R)
{
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		vnl_vector<double> oldpoint = VXLHelpers::vgl_point_to_vnl_vector(Points_[i].getCoord());
		vnl_vector<double> newpoint = R * oldpoint;
		 
		//Points_[i] = VXLHelpers::vnl_vector_to_vgl_point(newpoint);
		Points_[i].setCoord(VXLHelpers::vnl_vector_to_vgl_point(newpoint));
	}

//	PointTree_ = Octree(Points_);
	UpdateTriangles();	
}

void ModelFile::RotateAroundCenter(const vnl_double_3x3 &R)
{
	//1) center, 2) rotate, 3) move back

	//1) Center
	vgl_point_3d<double> OriginalCenter = getCenterOfMass();
	SetCenterOfMass(vgl_point_3d<double> (0,0,0));

	//2) rotate all the points
	Rotate(R);

	//3) move back to where it was
	SetCenterOfMass(OriginalCenter);
	
//	PointTree_ = Octree(Points_);
	UpdateTriangles();	
}

void ModelFile::Transform(const Transformation &Trans)
{
	Rotate(Trans.getRotation());
	Translate(Trans.getTranslation());
}


bool ModelFile::IntersectRay(const Ray &R, OrientedPoint &Intersection) const
{
	return Octree_.IntersectRay(R, Intersection);
}


void ModelFile::DisplayPoints()
{
	cout << "Points:" << endl << "-----------" << endl;
	for(unsigned int i = 0; i < Points_.size(); i++)
		cout << Points_[i].getCoord() << endl;
}

void ModelFile::UpdateTriangles()
{
	//clog << "Updating triangles..." << endl;
	
	if(VertexLists_.size() == 0)
		return;
	
	Triangles_.clear();
	
	//if all triangles are known to be valid
	//Triangles_.resize(VertexLists_.size());
	  
	//if we may need to remove bad triangles, do not set the size
	
	for(unsigned int n = 0; n < VertexLists_.size(); n++)
	{
		vector<unsigned int> Tri = VertexLists_[n];
		unsigned int ind0 = Tri[0];
		unsigned int ind1 = Tri[1];
		unsigned int ind2 = Tri[2];
		if(ind0 == ind1 || ind0 == ind2 || ind1 == ind2)
		{
			std::cout << "Two indices in triangle " << n << " are identical!" << std::endl;
			exit(-1);
		}
		
		vgl_point_3d<double> p0 = Points_[ind0].getCoord();
		vgl_point_3d<double> p1 = Points_[ind2].getCoord();
		vgl_point_3d<double> p2 = Points_[ind1].getCoord();
		if(p0 == p1 || p0 == p2 || p1 == p2)
		{
			std::cout << "Two points in triangle " << n << " are identical!" << std::endl;
			exit(-1);
		}
		
		unsigned int BadThings = 0;
		if(p0.x() == p1.x() && p0.x() == p2.x())
		{
			//std::cout << "All x coords for triangle " << n << " are the same" << std::endl;
			BadThings++;
		}
		
		if(p0.y() == p1.y() && p0.y() == p2.y())
		{
			//std::cout << "All y coords for triangle " << n << " are the same" << std::endl;
			BadThings++;
		}
		if(p0.z() == p1.z() && p0.z() == p2.z())
		{
			//std::cout << "All z coords for triangle " << n << " are the same" << std::endl;
			BadThings++;
		}
		
		//std::cerr << "P0: " << p0 << " P1: " << p1 << " P2: " << p2 << std::endl;
		
		//Triangles_[n] = Triangle(p0, p1, p2);
		
		if(BadThings > 1)
		{
			Tools::RemoveVectorElement(VertexLists_, n);
		}
		else
		{
			Triangles_.push_back(Triangle(p0, p1, p2));
		}

	}
}

/*
void ModelFile::ConvertToGraphicsCoordSystem()
{
	// scanner coordinate system is
	// z up, x and y plane is the ground
	
	// graphics coordinate system is
	// -z forward
	// x right
	// y up
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		vgl_point_3d<double> Temp = Points_[i].getCoord();
		Points_[i] = vgl_point_3d<double>(Temp.y(), Temp.z(), Temp.x());
	}
	
	UpdateTriangles();
}

void ModelFile::ConvertToScannerCoordSystem()
{
	// scanner coordinate system is
	// z up, x and y plane is the ground
	
	// graphics coordinate system is
	// -z forward
	// x right
	// y up
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		vgl_point_3d<double> Temp = Points_[i].getCoord();
		Points_[i] = vgl_point_3d<double>(Temp.z(), Temp.x(), Temp.y());
	}
	
	UpdateTriangles();
}
*/

void ModelFile::DisplayTriangles()
{
	std::cout << "Triangles:" << std::endl
			<< "-----------" << std::endl;
	
	for (unsigned int triangle = 0; triangle < VertexLists_.size(); triangle++)
	{
		std::cout << getTriangle(triangle);
	}
}



////////// Read/Writers ////////////
bool ModelFile::Read(const string &filename)
{
	Name_ = Tools::FileNameWithoutExtension(filename);
	//std::clog << "Set name to " << Name_ << std::endl;
	
	bool success;
	std::string ext = Tools::FileExtension(filename);
	if(ext == "vtp")
		success = VtpRead(filename);
	else if(ext == "obj")
		success = ObjRead(filename);
	else if(ext == "ply")
		success = PlyRead(filename);
	else if(ext == "vtu")
		success = VtuRead(filename);
	else if(ext == "ptx")
		success = PtxRead(filename);
	else if(ext == "txt")
		success = TxtCartesianRead(filename);
	else if(ext == "car")
		success = TxtCartesianRead(filename);
	else if(ext == "sph")
		success = TxtSphericalRead(filename);
	else
	{
		cout << "Invalid extension: " << ext << endl;
		exit(-1);
	}
	
	if(!success)
	{
		cout << "Warning: unable to read file " << filename << endl;
		Valid_ = false;
	}
	else
	{
		Valid_ = true;
	}
	
	return success;
}
		
bool ModelFile::Write(const std::string &filename) const
{
	
	std::string ext = Tools::FileExtension(filename);
	if(ext == "vtp")
		VtpWrite(filename);
	else if(ext == "obj")
		ObjWrite(filename);
	else if(ext == "vtu")
		VtuWrite(filename);
	else if(ext == "ptx")
		PtxWrite(filename);
	else if(ext == "txt")
		TxtWrite(filename);
	else if(ext == "ply")
		PlyWrite(filename);
	else
	{
		cout << "Invalid extension: " << ext << endl;
		return false;
	}
	
	return true;
}

bool ModelFile::PtxRead(const std::string &Filename)
{
	//note: ptx uses unsigned chars for color
	int verbosity = 0;

	std::ifstream infile;
	infile.open(Filename.c_str());
	if(!infile)
	{
		std::cout << "Could not open ptx file " << Filename << "!" << std::endl;
		return false;
	}

	std::string line;

	getline(infile, line);
	unsigned int PhiPoints_;
	std::stringstream(line) >> PhiPoints_;

	getline(infile, line);
	unsigned int ThetaPoints_;
	std::stringstream(line) >> ThetaPoints_;

	//skip 8 lines (identity matrices)
	for(int i = 0; i < 8; i++)
	{
		getline(infile, line);
	}

	std::vector<double> Point(3);
	double Intensity;
	std::vector<double> Intensities_;

	std::vector<unsigned char> ColorChar(3);
	std::vector<int> ColorInt(3);

	//while(getline(infile, line))
	for(unsigned int counter = 0; counter < PhiPoints_*ThetaPoints_; counter++)
	{
		getline(infile, line);
				
		vgl_point_3d<double> P;
		std::stringstream ParsedLine(line);
		ParsedLine >> Point[0] >> Point[1] >> Point[2] >> Intensity >> ColorInt[0] >> ColorInt[1] >> ColorInt[2];
	
		for(int i = 0; i < 3; i++)
		{
			ColorChar[i] = (unsigned char)ColorInt[i];
		}
	
		if(verbosity >= 1)
		{
			clog << ColorInt[0] << " " << ColorInt[1] << " " << ColorInt[2] << endl;
			clog << int(ColorChar[0]) << " " << int(ColorChar[1]) << " " << int(ColorChar[2]) << endl;
		}
	
		if( !((Point[0] == 0) && (Point[1] == 0) && (Point[2] == 0)) ) //point is a valid point
		{
			P = vgl_point_3d<double>(Point[0], Point[1], Point[2]);
			//Points_.push_back(P);
			Points_.push_back(OrientedPoint(P, ColorChar));
			Intensities_.push_back(Intensity);
			
		}
	

	}//end while getline


	infile.close();
	return true;

}

bool ModelFile::PtxWrite(const string &Filename) const
{

	ofstream fout(Filename.c_str(), ios::out);

	//fout << ScanParams.PhiPoints_ << endl
	//	<< SacnParams.ThetaPoints_ << endl
	fout << static_cast<unsigned int>(sqrt(NumPoints())) << endl
		<< static_cast<unsigned int>(sqrt(NumPoints())) << endl
		<< "0 0 0" << endl
		<< "1 0 0" << endl
		<< "0 1 0" << endl
		<< "0 0 1" << endl
		<< "1 0 0 0" << endl
		<< "0 1 0 0" << endl
		<< "0 0 1 0" << endl
		<< "0 0 0 1" << endl;
	
	for (unsigned int i = 0; i < NumPoints(); i++)
	{
		//fout << Points_[i] << " " << Intensities_[i] << " " << Colors_[i] << endl;
		if(NumColors() == NumPoints())
		{
			//fout << *Points_[i] << " " << 0 << " " << Colors_[i] << endl;
			fout << Points_[i].getCoord().x() << " " << Points_[i].getCoord().y() << " " << Points_[i].getCoord().z() << " " << 0 << " " << Points_[i].getColor() << endl;
		}
		else
		{
			//fout << *Points_[i] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
			fout << Points_[i].getCoord().x() << " " << Points_[i].getCoord().y() << " " << Points_[i].getCoord().z() << " " << 0 << " " << 0 << " "  << 0 << " "  << 0  << endl;
		}
		
	}
	
	fout.close();
		
	return true;//write successful
}
					   
bool ModelFile::VtuRead(const string &filename)
{	
	//read all the data from the file
	vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkUnstructuredGrid* ug = reader->GetOutput();

	//get points
	vtkIdType NumPoints = ug->GetNumberOfPoints();
		
	if(!(NumPoints > 0) )
		return false;
	
	Points_.resize(NumPoints);
	double point[3];
	
	for(unsigned int i = 0; i < static_cast<unsigned int>(NumPoints); i++) // !!! warning: comparison between signed and unsigned ints
	{
		ug->GetPoint(i, point);
		Points_[i].setCoord(vgl_point_3d<double> (point));
	}
	
	//get triangles
	vtkIdType NumCells = ug->GetNumberOfCells();
	VertexLists_.clear();
	if( (ug->GetCellType(0) == VTK_TRIANGLE) && (NumCells > 0) )//vtkCellType.h
	{
		for(vtkIdType tri = 0; tri < NumCells; tri++)
		{
			vtkCell* cell = ug->GetCell(tri);
			vtkIdList* pts = cell->GetPointIds();
			vector<unsigned int> List(3);
			List[0] = pts->GetId(0);
			List[1] = pts->GetId(1);
			List[2] = pts->GetId(2);
		
			VertexLists_.push_back(List);
		}
		
		
		UpdateTriangles();
	}
	
	
	//get colors
	vtkUnsignedCharArray* ColorsData = vtkUnsignedCharArray::SafeDownCast(ug->GetPointData()->GetArray("Colors"));
	
	if(ColorsData)
	{ 
		unsigned char color[3];
		
		for(unsigned int i = 0; i < static_cast<unsigned int> (NumPoints); i++)
		{
			ColorsData->GetTupleValue(i, color);
			Color<unsigned char> c(color[0], color[1], color[2]);
			Points_[i].setColor(c);
		}
	}
	

	std::clog << "Finished reading vtu file " << filename << "." << endl;
	std::clog << "Points: " << NumPoints << " Colors: " << NumColors() << " Triangles: " << NumTriangles() << endl;

	return true; //file read successfully
}


bool ModelFile::VtuWrite(const string &Filename) const
{

	vtkPoints* points3D = vtkPoints::New();
	vtkCellArray* Vertices = vtkCellArray::New();

	for ( unsigned int i = 0; i < NumPoints(); i++ )
	{	
		vtkIdType pid[1];
		vgl_point_3d<double> Point = Points_[i].getCoord();
		pid[0] = points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		Vertices->InsertNextCell(1,pid);
	}
	
	vtkUnsignedCharArray* Colors = vtkUnsignedCharArray::New();
	Colors->SetNumberOfComponents(3);
	Colors->SetName("Colors");

	for ( unsigned int i = 0; i < Points_.size(); i++ )
	{
		Color<unsigned char> Color = Points_[i].getColor();
		unsigned char ColorArray[3];
		CharArray(Color, ColorArray);
		Colors->InsertNextTupleValue(ColorArray);
	}

	bool HasTriangles;
	if(NumTriangles() > 0)
		HasTriangles = true;
	else
		HasTriangles = false;
	
	vtkCellArray* triangles = vtkCellArray::New();
	if(HasTriangles)
	{
		for(unsigned int i = 0; i < NumTriangles(); i++)
		{
			vector<unsigned int> vlist = getTriangleVertices(i);
			vtkTriangle* triangle = vtkTriangle::New();
			triangle->GetPointIds()->SetId(0,vlist[0]);
			triangle->GetPointIds()->SetId(1,vlist[1]);
			triangle->GetPointIds()->SetId(2,vlist[2]);
			triangles->InsertNextCell(triangle);
		}
	}

	vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();

	ug->SetPoints(points3D);
	
	//add triangles
	if(HasTriangles)
		ug->SetCells(VTK_TRIANGLE, triangles);
	else
		ug->SetCells(VTK_VERTEX, Vertices);

	//add colors
	ug->GetPointData()->AddArray(Colors);

	vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
	writer->SetFileName(Filename.c_str());
	writer->SetInput(ug);
	writer->Write();	

	std::clog << "Finished writing vtu file " << Filename << "." << endl;
	std::clog << "NumPoints: " << NumPoints() << " NumColors: " << NumColors() << endl;
	
	return true;//write finshed ok
}

void ModelFile::ConstructFromPolydata(vtkSmartPointer<vtkPolyData> polydata)
{
	//get points
	vtkIdType idNumPointsInFile = polydata->GetNumberOfPoints();
	unsigned int NumPointsInFile = static_cast<unsigned int> (idNumPointsInFile);
	
	if(!(NumPointsInFile > 0) )
	{
		Valid_ = false;
		return;
	}
	
	Points_.resize(NumPointsInFile);
	
	double point[3];
	for(unsigned int i = 0; i < NumPointsInFile; i++)
	{
		polydata->GetPoint(i, point);
		Points_[i].setCoord(vgl_point_3d<double>(point));
	}
	
	//get triangles
	vtkIdType NumPolys = polydata->GetNumberOfPolys();
	if(NumPolys > 0)
	{
		vtkCellArray* TriangleCells = polydata->GetPolys();
		if(!TriangleCells)
		{
			std::cout << "TriangleCells is null!" << std::endl;
			exit(-1);
		}
		vtkIdType npts;
		vtkIdType *pts;

		VertexLists_.clear();
		
		TriangleCells->InitTraversal();
		while(TriangleCells->GetNextCell(npts, pts))
		{
			vector<unsigned int> List(3);
			List[0] = pts[0];
			List[1] = pts[1];
			List[2] = pts[2];
			
			VertexLists_.push_back(List);
		}
		
		UpdateTriangles();
	}
	
	
	//get colors
	vtkUnsignedCharArray* ColorsData = vtkUnsignedCharArray::SafeDownCast(polydata->GetPointData()->GetArray("Colors"));
	
	if(ColorsData)
	{
		ValidColors_ = true;
		
		unsigned int nColors = polydata->GetPointData()->GetArray("Colors")->GetNumberOfTuples();
		if(nColors != NumPointsInFile)
		{
			cout << "The length of Colors is " << nColors << " and the number of points is " << NumPointsInFile << endl;
			exit(-1);
		}
		unsigned char color[3];
		
		for(unsigned int i = 0; i < NumPoints(); i++)
		{
			ColorsData->GetTupleValue(i, color);
			Color<unsigned char> c(color[0], color[1], color[2]);
			Points_[i].setColor(c);
		}
	}
	
	//get normals
	std::vector<vgl_vector_3d<double> > Normals = GetNormalsFromPolydata(polydata);
	if(Normals.size() == NumPointsInFile)
	{
		clog << "Normals found on first pass." << endl;
		for(unsigned int i = 0; i < Normals.size(); i++)
		{
			Points_[i].setNormal(Normals[i]);
		}
	}
	else
	{
		//if none of the normal types we looked for exist, then create them from the cells
		std::clog << "Normals NOT found on first pass, extracting..." << endl;
		AddNormalsToPolyData(polydata);
		
		std::vector<vgl_vector_3d<double> > NewNormals = GetNormalsFromPolydata(polydata);
		std::clog << "There are now " << NewNormals.size() << " normals." << endl;
		if(NewNormals.size() == NumPointsInFile)
		{
			for(unsigned int i = 0; i < NewNormals.size(); i++)
			{
				Points_[i].setNormal(NewNormals[i]);
			}
		}
		else
		{
			std::cout << "Number of normals extracted does not match number of points!" << std::endl;
			exit(-1);
		}
	}
	
	
	std::clog << "Getting TotalInformation..." << endl;
	////////////// Information /////////////////
	vtkDoubleArray* TotalInformation = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("TotalInformation"));
	if(TotalInformation)
	{
		//unsigned int nInfo = polydata->GetPointData()->GetArray("TotalInformation")->GetNumberOfTuples();
		unsigned int nInfo = TotalInformation->GetNumberOfTuples();
		std::clog << "nInfo: " << nInfo << endl;
		if(nInfo != NumPointsInFile)
		{
			std::cout << "The length of TotalInformation is " << nInfo << " and the number of points is " << NumPointsInFile << endl;
			exit(-1);
		}
		for(unsigned int i = 0; i < NumPointsInFile; i++)
		{	
			double info;
			info = TotalInformation->GetValue(i);
			//Points_[i].setTotalInformation(info);
			Points_[i].setDoubleValue("TotalInformation",info);
		}
	}
	else
	{
		std::clog << "No TotalInformation!" << std::endl;
	}
	
	std::clog << "Getting RemainingInformation..." << endl;
	////////////// Information /////////////////
	vtkDoubleArray* RemainingInformation = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("RemainingInformation"));
	if(RemainingInformation)
	{
		//unsigned int nInfo = polydata->GetPointData()->GetArray("TotalInformation")->GetNumberOfTuples();
		unsigned int nInfo = RemainingInformation->GetNumberOfTuples();
		std::clog << "nInfo: " << nInfo << endl;
		if(nInfo != NumPointsInFile)
		{
			std::cout << "The length of RemainingInformation is " << nInfo << " and the number of points is " << NumPointsInFile << endl;
			exit(-1);
		}
		for(unsigned int i = 0; i < NumPointsInFile; i++)
		{	
			double info;
			info = RemainingInformation->GetValue(i);
			//Points_[i].setRemainingInformation(info);
			Points_[i].setDoubleValue("RemainingInformation", info);
		}
	}
	else
	{
		std::clog << "No TotalInformation!" << std::endl;
	}

	
	/////////////////////// Field Data ///////////////////////////
	//get scanner location
	vtkDataArray* ScannerLocationArray;
	
	vtkFieldData* FieldData = polydata->GetFieldData();

	if(FieldData)
		ScannerLocationArray = FieldData->GetArray("ScannerLocation");
	
	if(ScannerLocationArray)
	{
		double ScannerLocation[3];
		ScannerLocationArray->GetTuple(0, ScannerLocation);
		ScannerLocation_ = vgl_point_3d<double>(ScannerLocation[0], ScannerLocation[1], ScannerLocation[2]);
		IsScan_ = true;
		std::clog << "Got scanner location " << ScannerLocation_ << endl;
	}
	else
	{
		//unless we specify that this is a scan file, we assume it is not
		std::clog << "No scanner location!" << endl;
		IsScan_ = false;
	}
		
	Valid_ = true;
}

void AddNormalsToPolyData(vtkSmartPointer<vtkPolyData> &polydata)
{
	std::vector<vgl_vector_3d<double> > Normals = GetNormalsFromPolydata(polydata);
	std::clog << "There are " << Normals.size() << " normals before extraction." << std::endl;
	
	//vtkPolyDataNormals* normalGenerator = vtkPolyDataNormals::New();
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkPolyDataNormals::New();
	normalGenerator->SetInput(polydata);
	
	//point normals only
	normalGenerator->SetComputePointNormals(1);
	normalGenerator->SetComputeCellNormals(0);
	
	//cell normals only
	//normalGenerator->SetComputePointNormals(0);
	//normalGenerator->SetComputeCellNormals(1);
	
	//point and cell normals
	//normalGenerator->SetComputePointNormals(1);
	//normalGenerator->SetComputeCellNormals(1);
	
	normalGenerator->SetSplitting(0);
	normalGenerator->Update();
	polydata = normalGenerator->GetOutput(); //returned by pointer/reference
	
	Normals = GetNormalsFromPolydata(polydata);
	std::clog << "There are " << Normals.size() << " normals after extraction." << std::endl;
}

std::vector<vgl_vector_3d<double> > GetNormalsFromPolydata(vtkSmartPointer<vtkPolyData> polydata)
{
	//get points
	vtkIdType idNumPointsInPolydata = polydata->GetNumberOfPoints();
	unsigned int NumPointsInPolydata = static_cast<unsigned int> (idNumPointsInPolydata);

	
	std::vector<vgl_vector_3d<double> > Normals;
	std::clog << "Getting normals...." << endl;
	
	//normals in an array
	vtkDoubleArray* NormalData1 = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));
	
	if(NormalData1)
	{ 
		unsigned int nNormals1 = polydata->GetPointData()->GetArray("Normals")->GetNumberOfTuples();
		if(nNormals1 != NumPointsInPolydata)
		{
			cout << "There are " << nNormals1 << " components in NormalData1 (double Point \"Normals\")" << endl;
			exit(-1);
		}
	
		double normal[3];
		
		for(unsigned int i = 0; i < NumPointsInPolydata; i++)
		{
			NormalData1->GetTupleValue(i, normal);
			vgl_vector_3d<double> n(normal[0], normal[1], normal[2]);
			Normals.push_back(n);
		}
		return Normals;
	}
	
	
	//actual point normals
	//vtkSmartPointer<vtkDoubleArray> NormalData2 = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetNormals());
	vtkDataArray* PointNormals = polydata->GetPointData()->GetNormals();
	
	if(PointNormals)
	{
		unsigned int nPointNormals = PointNormals->GetNumberOfTuples();
		if(nPointNormals != NumPointsInPolydata)
		{
			cout << "There are " << nPointNormals << " components in PointNormals and " << NumPointsInPolydata << " points." << endl;
			exit(-1);
		}
		
		Normals.clear();
		double normal[3];
		
		for(unsigned int i = 0; i < NumPointsInPolydata; i++)
		{
			//NormalData2->GetTupleValue(i, normal);
			PointNormals->GetTuple(i, normal);
			vgl_vector_3d<double> n(normal[0], normal[1], normal[2]);
			Normals.push_back(n);
		}
		return Normals;
	}
		
	//double normals
	vtkDataArray* CellNormals = polydata->GetCellData()->GetNormals("cellNormals");
	
	if(CellNormals)
	{
		unsigned int nCellNormals = CellNormals->GetNumberOfTuples();
		if(nCellNormals != NumPointsInPolydata)
		{
			cout << "There are " << nCellNormals << " components in CellDatal \"Cell Normals\")" << " and " << NumPointsInPolydata << "points."<< endl;
			exit(-1);
		}
	
		double normal[3];
		
		for(unsigned int i = 0; i < NumPointsInPolydata; i++)
		{
			//CellNormals->GetTupleValue(i, normal);
			CellNormals->GetTuple(i, normal);
			vgl_vector_3d<double> n(normal[0], normal[1], normal[2]);
			Normals.push_back(n);
		}

		return Normals;
	}
	
	std::clog << "There are " << Normals.size() << " Normals." << std::endl;
	
	return Normals;
}


bool ModelFile::VtpRead(const string &filename)
{
	//get all data from the file
	//vtkXMLPolyDataReader* reader = vtkXMLPolyDataReader::New();
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	std::clog << "Reading " << filename << endl;
	reader->SetFileName(filename.c_str());
	reader->Update();
	//reader->ReadAllNormals();
	//cout << "There are " << reader->GetNumberOfNormalsInFile() << " normals in the file." << endl;
	
	vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
	//vtkPolyData* polydata = reader->GetOutput();

	ConstructFromPolydata(polydata);
	
	std::clog << "Finished reading vtp file " << filename << "." << endl;
	std::clog << "NumPoints: " << NumPoints() << " NumTris: " << NumTriangles() << " NumColors: " << NumColors() << endl;
	
	return Valid_;
}

vtkSmartPointer<vtkPolyData> ModelFile::CreatePolydata(void) const
{
	//check for triangles
	bool HasTriangles;
	if(NumTriangles() > 0)
		HasTriangles = true;
	else
		HasTriangles = false;
	
	//create a polydata
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	
	vtkSmartPointer<vtkPoints> points3D = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();

	for(unsigned int i = 0; i < NumPoints(); i++)
	{	
		vtkIdType pid[1];
		vgl_point_3d<double> Point = Points_[i].getCoord();
		pid[0] = points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		Vertices->InsertNextCell(1,pid);
		
	}

	//add points/vertices
	polydata->SetPoints(points3D);
	if(!HasTriangles)
	{
		std::clog << "No triangle, setting up vertices." << std::endl;
		polydata->SetVerts(Vertices);
	}
	
	//////////// Point indices ///////////////
	std::vector<unsigned int> PointIndices(NumPoints());
	Tools::CreateIndexVector(PointIndices);
	VTKHelpers::AddArrayToPolydata<unsigned int, vtkUnsignedIntArray>(PointIndices, "PointIndex", polydata);
	
	////////////// Other Double Arrays /////////////////
	std::map<std::string, std::vector<double> >::const_iterator ArrayIter = DoubleArrays.begin();
	for( ; ArrayIter != DoubleArrays.end(); ++ArrayIter )
	{
		if(ArrayIter->second.size() == NumPoints())
			VTKHelpers::AddArrayToPolydata<double, vtkDoubleArray>(ArrayIter->second, ArrayIter->first, polydata);
	}
	
	////////////// Information /////////////////
	std::vector<double> TotalInformation;
	std::vector<double> RemainingInformation;
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		//TotalInformation.push_back(Points_[i].getTotalInformation());
		TotalInformation.push_back(Points_[i].getDoubleValue("TotalInformation"));
		//RemainingInformation.push_back(Points_[i].getRemainingInformation());
		RemainingInformation.push_back(Points_[i].getDoubleValue("RemainingInformation"));
	}
	VTKHelpers::AddArrayToPolydata<double, vtkDoubleArray>(TotalInformation, "TotalInformation", polydata);
	VTKHelpers::AddArrayToPolydata<double, vtkDoubleArray>(RemainingInformation, "RemainingInformation", polydata);
	
	
	/////////////// Colors //////////////////
	vtkSmartPointer<vtkUnsignedCharArray> Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	if(ValidColors_)
	{
		Colors->SetNumberOfComponents(3);
		Colors->SetName("Colors");
	
		for(unsigned int i = 0; i < NumPoints(); i++)
		{
			Color<unsigned char> Color = Points_[i].getColor();
			unsigned char ColorArray[3];
			CharArray(Color, ColorArray);
			Colors->InsertNextTupleValue(ColorArray);
		}
	}

	/////////// Triangles /////////////
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	if(HasTriangles)
	{
		for(unsigned int i = 0; i < NumTriangles(); i++)
		{
			vector<unsigned int> vlist = getTriangleVertices(i);
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			triangle->GetPointIds()->SetId(0,vlist[0]);
			triangle->GetPointIds()->SetId(1,vlist[1]);
			triangle->GetPointIds()->SetId(2,vlist[2]);
			triangles->InsertNextCell(triangle);
		}
	}
	
	/////////// Normals /////////////
	vtkSmartPointer<vtkFloatArray> norms = vtkSmartPointer<vtkFloatArray>::New();
	norms->SetNumberOfComponents(3);
	norms->SetNumberOfTuples(NumPoints());
	norms->SetName("Normals");
	
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		float norm[3];
		norm[0] = static_cast<float>(Points_[i].getNormal().x());
		norm[1] = static_cast<float>(Points_[i].getNormal().y());
		norm[2] = static_cast<float>(Points_[i].getNormal().z());
		norms->SetTuple(i, norm) ;
	}


	//if this is a scan, add the scanner location
	
	if(IsScan_)
	{
		std::clog << "Adding scanner location to " << Name_ << ": " << ScannerLocation_ << std::endl;
		double loc[3];
		loc[0] = ScannerLocation_.x();
		loc[1] = ScannerLocation_.y();
		loc[2] = ScannerLocation_.z();
		
		vtkSmartPointer<vtkDoubleArray> Location = vtkSmartPointer<vtkDoubleArray>::New();
		Location->SetNumberOfComponents(3);
		Location->SetName("ScannerLocation");
		Location->InsertNextTuple(loc);
		polydata->GetFieldData()->AddArray(Location);
	}
	else
	{
		std::clog << "File " << Name_ << " is not a scan! (no scanner location)" << std::endl;
	}
	
	//add triangles (if there are any)
	if(HasTriangles)
	{
		polydata->SetPolys(triangles);
	}

	//add normals
	polydata->GetPointData()->SetNormals(norms);
	
	//add colors
	if(ValidColors_)
	{
		polydata->GetPointData()->SetVectors(Colors);
	}
			
	return polydata;
}

bool ModelFile::VtpWrite(const string &Filename) const
{
	vtkSmartPointer<vtkPolyData> polydata = CreatePolydata();
	std::clog << "Number of Cells: " << polydata->GetNumberOfCells() << endl;
	
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(Filename.c_str());
	writer->SetInput(polydata);
	writer->Write();	

	std::clog << "Finished writing vtp file " << Filename << "." << endl;
	std::clog << "NumPoints: " << NumPoints() << endl;
	std::clog << "NumTris: " << NumTriangles() << endl;

	return true;//write finshed ok
}

bool ModelFile::SimpleVtpWrite(const std::string &Filename) const
{
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	
	vtkSmartPointer<vtkPoints> points3D = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();
	
	for(unsigned int i = 0; i < NumPoints(); i++)
	{	
		vtkIdType pid[1];
		vgl_point_3d<double> Point = Points_[i].getCoord();
		pid[0] = points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		Vertices->InsertNextCell(1,pid);
	}
	
	//add points/vertices
	polydata->SetPoints(points3D);
	polydata->SetVerts(Vertices);
	
	/////////////// Colors //////////////////
	vtkSmartPointer<vtkUnsignedCharArray> Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();

	Colors->SetNumberOfComponents(3);
	Colors->SetName("Colors");

	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		Color<unsigned char> Color = Points_[i].getColor();
		unsigned char ColorArray[3];
		CharArray(Color, ColorArray);
		Colors->InsertNextTupleValue(ColorArray);
	}
	
	
	//add colors
	polydata->GetPointData()->SetVectors(Colors);

	/////////// Triangles /////////////
	if(NumTriangles() > 0)
	{
		vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
		
		for(unsigned int i = 0; i < NumTriangles(); i++)
		{
			vector<unsigned int> vlist = getTriangleVertices(i);
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			triangle->GetPointIds()->SetId(0,vlist[0]);
			triangle->GetPointIds()->SetId(1,vlist[1]);
			triangle->GetPointIds()->SetId(2,vlist[2]);
			triangles->InsertNextCell(triangle);
		}
		
		//add triangles
		polydata->SetPolys(triangles);
	}
	
	/////////// Normals /////////////
	vtkSmartPointer<vtkFloatArray> norms = vtkSmartPointer<vtkFloatArray>::New();
	norms->SetNumberOfComponents(3);
	norms->SetNumberOfTuples(NumPoints());
	norms->SetName("Normals");
	
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		float norm[3];
		norm[0] = static_cast<float>(Points_[i].getNormal().x());
		norm[1] = static_cast<float>(Points_[i].getNormal().y());
		norm[2] = static_cast<float>(Points_[i].getNormal().z());
		norms->SetTuple( i, norm ) ;
	}

	//add normals
	polydata->GetPointData()->SetNormals( norms );
	
	//Write file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(Filename.c_str());
	writer->SetInput(polydata);
	writer->Write();	

	//output info
	std::clog << "Finished writing simple vtp file " << Filename << "." << endl;
	std::clog << "NumPoints: " << NumPoints() << endl;
	std::clog << "NumTris: " << NumTriangles() << endl;

	return true;//write finshed ok
}

bool ModelFile::ObjRead(const string &filename)
{
	vtkOBJReader* reader = vtkOBJReader::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	
	vtkPolyData* polydata = reader->GetOutput();
	vtkIdType NumPoints = polydata->GetNumberOfPoints();
	
	Points_.clear();
	Points_.resize(NumPoints);
	
	double point[3];
	for(int i = 0; i < NumPoints; i++)
	{
		polydata->GetPoint(i, point);
		
		vgl_point_3d<double> p = vgl_point_3d<double>(point);
		Points_[i].setCoord(p);
	}
	
	//create triangle vertex lists
	VertexLists_.clear();
	vtkIdType NumPolys = polydata->GetNumberOfPolys();
	if(NumPolys > 0)
	{
		vtkCellArray* TriangleCells = polydata->GetPolys();
		vtkIdType npts;
		vtkIdType *pts;

		VertexLists_.clear();
		TriangleCells->InitTraversal();
		while(TriangleCells->GetNextCell(npts, pts))
		{
			vector<unsigned int> List(3);
			List[0] = pts[0];
			List[1] = pts[1];
			List[2] = pts[2];
			
			VertexLists_.push_back(List);
		}	
		
		UpdateTriangles();
	}
	
	
	return true;
	
}

bool ModelFile::ObjWrite(const string &Filename) const
{
	//NOT WORKING - vtkOBJExporter does not take polydata

	vtkPoints* points3D = vtkPoints::New();
	vtkCellArray* Vertices = vtkCellArray::New();

	for ( unsigned int i = 0; i < NumPoints(); ++i )
	{	
		vtkIdType pid[1];
		vgl_point_3d<double> Point = Points_[i].getCoord();
		pid[0] = points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		Vertices->InsertNextCell(1,pid);
	}

	bool HasColors;
	if(NumColors())
		HasColors = true;
	else
		HasColors = false;

	vtkUnsignedCharArray* Colors = vtkUnsignedCharArray::New();
	if(HasColors)
	{
		Colors->SetNumberOfComponents(3);
		Colors->SetName("Colors");
	
		for ( unsigned int i = 0; i < NumColors(); ++i )
		{
			Color<unsigned char> Color = getColor(i);
			unsigned char ColorArray[3];
			ColorArray[0] = Color.getR();
			ColorArray[1] = Color.getG();
			ColorArray[2] = Color.getB();
			Colors->InsertNextTupleValue(ColorArray);
		}
	}

	bool HasTriangles;
	if(NumTriangles() > 0)
		HasTriangles = true;
	else
		HasTriangles = false;
	
	vtkCellArray* triangles = vtkCellArray::New();
	if(HasTriangles)
	{
		for(unsigned int i = 0; i < NumTriangles(); i++)
		{
			vector<unsigned int> vlist = getTriangleVertices(i);
			vtkTriangle* triangle = vtkTriangle::New();
			triangle->GetPointIds()->SetId(0,vlist[0]);
			triangle->GetPointIds()->SetId(1,vlist[1]);
			triangle->GetPointIds()->SetId(2,vlist[2]);
			triangles->InsertNextCell(triangle);
		}
	}

	vtkPolyData* polydata = vtkPolyData::New();

	polydata->SetPoints(points3D);
	polydata->SetVerts(Vertices);
	
	if(HasTriangles)
		polydata->SetPolys(triangles);

	if(HasColors)
		polydata->GetPointData()->AddArray(Colors);

	vtkOBJExporter* writer = vtkOBJExporter::New();
	//writer->SetFileName(Filename.c_str());
	writer->SetFilePrefix(Filename.c_str());
//	writer->SetInput(polydata);

	writer->Write();	

	std::clog << "Finished writing vtp file " << Filename << "." << endl;
	std::clog << "NumPoints: " << NumPoints() << endl;
	std::clog << "NumColors: " << NumColors() << endl;

	return true;//write finshed ok
}


bool ModelFile::PlyRead(const string &filename)
{
	vtkPLYReader* reader = vtkPLYReader::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	
	vtkPolyData* polydata = reader->GetOutput();
	vtkIdType NumPoints = polydata->GetNumberOfPoints();
	
	Points_.clear();
	Points_.resize(NumPoints);
	
	double point[3];
	for(int i = 0; i < NumPoints; i++)
	{
		polydata->GetPoint(i, point);
		
		vgl_point_3d<double> p = vgl_point_3d<double>(point);
		Points_[i].setCoord(p);
	}
	
	//create triangle vertex lists
	VertexLists_.clear();
	vtkIdType NumPolys = polydata->GetNumberOfPolys();
	if(NumPolys > 0)
	{
		vtkCellArray* TriangleCells = polydata->GetPolys();
		vtkIdType npts;
		vtkIdType *pts;

		VertexLists_.clear();
		TriangleCells->InitTraversal();
		while(TriangleCells->GetNextCell(npts, pts))
		{
			vector<unsigned int> List(3);
			List[0] = pts[0];
			List[1] = pts[1];
			List[2] = pts[2];
			
			VertexLists_.push_back(List);
		}	
		
		UpdateTriangles();
	}
	
	
	return true;
	
}


bool ModelFile::PlyWrite(const string &Filename) const
{
	//NOT TESTED !!!

	
	vtkPolyData* polydata = vtkPolyData::New();

	vtkPoints* Points = vtkPoints::New();
	vtkCellArray* Vertices = vtkCellArray::New();

	for(unsigned int i = 0; i < NumPoints(); ++i )
	{	
		vtkIdType pid[1];
		vgl_point_3d<double> Point = Points_[i].getCoord();
		pid[0] = Points->InsertNextPoint(Point.x(), Point.y(), Point.z());
		Vertices->InsertNextCell(1,pid);
	}

	polydata->SetPoints(Points);
	polydata->SetVerts(Vertices);
	
	if(NumTriangles() > 0)
	{
		vtkCellArray* triangles = vtkCellArray::New();
		for(unsigned int i = 0; i < NumTriangles(); i++)
		{
			vector<unsigned int> vlist = getTriangleVertices(i);
			vtkTriangle* triangle = vtkTriangle::New();
			triangle->GetPointIds()->SetId(0,vlist[0]);
			triangle->GetPointIds()->SetId(1,vlist[1]);
			triangle->GetPointIds()->SetId(2,vlist[2]);
			triangles->InsertNextCell(triangle);
		}
		
		polydata->SetPolys(triangles);
	}


	vtkPLYWriter* writer = vtkPLYWriter::New();
	writer->SetFileName(Filename.c_str());
	writer->SetInput(polydata);

	writer->Write();	

	std::clog << "Finished writing ply file " << Filename << "." << endl;
	std::clog << "NumPoints: " << NumPoints() << endl;
	std::clog << "NumTris: " << NumTriangles() << endl;

	return true;//write finshed ok
	
}


bool ModelFile::TxtWrite(const string &Filename) const
{
	ofstream fout(Filename.c_str(), ios::out);

	//write vertices
	for(unsigned int i = 0; i < NumVertices(); i++)
	{
		vgl_point_3d<double> P = Points_[i].getCoord();
		vgl_vector_3d<double> N = Points_[i].getNormal();
		//fout << P.x() << " " << P.y() << " " << P.z() << endl;
		fout << P.x() << " " << P.y() << " " << P.z() << " " << N.x() << " " << N.y() << " " << N.z() << endl;
	}
	
	fout.close();
	
	cout << "Finished writing txt file " << Filename << "." << endl
			<< *this << endl;
	
	return true;//write finshed ok
}


bool ModelFile::TxtCartesianRead(const std::string &Filename)
{
	std::vector<vgl_point_3d<double> > Points;
	
	std::ifstream fin(Filename.c_str(), ios::out);

	std::string line;
	while(getline(fin, line))
	{
		std::cout << "line: " << line << std::endl;
		double x, y, z;
		std::stringstream linestream;
		linestream << line;
		linestream >> x >> y >> z;
		std::cout << x << " " << y << " " << z << std::endl;
		Points_.push_back(vgl_point_3d<double> (x,y,z));
	}
	
	fin.close();
	
	std::clog << "Finished reading txt file " << Filename << "." << endl
			<< *this << endl;
	
	return true;
}


bool ModelFile::TxtSphericalRead(const std::string &Filename)
{
	std::vector<vgl_point_3d<double> > Points;
	
	std::ifstream fin(Filename.c_str(), ios::out);

	std::string line;
	while(getline(fin, line))
	{
		//double x, y, z;
		double p, t, r;
		std::stringstream linestream;
		linestream << line;
		linestream >> p >> t >> r;

		p = deg2rad(p);
		t = deg2rad(t);


	//kinematics
	
	double d1 = 75;
	double d2 = 53;
	double d3 = 65;
	double p2x = d1 * cos(t);
	double p2y = d1 * sin(t);
	double p2z = d2;

	double cx = p2x + d3*sin(p)*cos(t);
	double cy = p2y + d3*sin(p)*sin(t);
	double cz = p2z + d3*cos(p);

	double x = cx + r*sin(p+M_PI/2)*cos(t);
	double y = cy + r*sin(p+M_PI/2)*sin(t);
	double z = cz + r*cos(p+M_PI/2);


		/*
		if(r <= 0)
			continue;
		
		//convert spherical to cartesian
		vgl_vector_3d<double> Direction = Sphere2Rect(deg2rad(p), deg2rad(t));
		
		normalize(Direction);

		Direction *= r;

		//if( (fabs(Direction.x() > 1e5)) || (fabs(Direction.y() > 1e5)) || (fabs(Direction.z() > 1e5)) )
		//	continue;

		vgl_point_3d<double> Point = vgl_vector_to_vgl_point(Direction);

		//if( (fabs(Point.x() > 1e5)) || (fabs(Point.y() > 1e5)) || (fabs(Point.z() > 1e5)) )
		//	continue;

*/

		vgl_point_3d<double> Point = vgl_point_3d<double> (x,y,z);

		cout << Point << endl;
		
		Points_.push_back(Point);
	}
	
	fin.close();
	
	//OriginalPoints_ = Points_;
	
	std::cout << "Finished reading txt file " << Filename << "." << endl
			<< *this << endl;
	
	return true;
}


double ThetaAngle(const vgl_point_3d<double> &P, string CoordSystem)
{
	double epsilon = 1e-6;
	if(CoordSystem == "d")//dave
	{
		if(fabs(P.z()) < epsilon)
			return 0;
		else
		{
			
			//return atan2(P.getX(), P.getZ());
			return atan2(P.z(), P.x());
		}
	}
	else if(CoordSystem == "s")//scanner
	{
		if(fabs(P.x()) < epsilon)
			return 0;
		else
		{
			return atan2(P.y(), fabs(P.x()));
		}
	}
	
	assert(0);
	exit(0);
	return 0.0;
}

double PhiAngle(const vgl_point_3d<double> &P, string CoordSystem)
{
	double epsilon = 1e-6;

	if(CoordSystem == "d")//dave
	{
		if(fabs(P.z()) < epsilon)
			return 0;
		else
		{
			//return atan2(P.getY(), P.getZ());
			return atan2(P.z(), P.y());
		}
	}
	else if(CoordSystem == "s")//scanner
	{
		if(fabs(P.z()) < epsilon)
			return 0;
		else
		{
			return atan2(P.z(), fabs(P.x()));
		}
	}
	
	assert(0);
	exit(0);
	return 0.0;
}


/*
vtkPoints* ModelFile::GetVTKPoints() const
{
	vtkPoints* points3D = vtkPoints::New();
	//vtkCellArray* Vertices = vtkCellArray::New();

	for(int i = 0; i < NumPoints(); i++)
	{
		//vtkIdType pid[1];
		vgl_point_3d<double> Point = getPoint(i);
		//pid[0] = points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		//Vertices->InsertNextCell(1,pid);
	}

	return points3D;
}
*/

vtkSmartPointer<vtkPolyData> ModelFile::GetVTKPoints() const
{
	vtkSmartPointer<vtkPoints> points3D = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();

	for(unsigned int i = 0; i < NumPoints(); i++)
	{	
		vtkIdType pid[1];
		vgl_point_3d<double> Point = Points_[i].getCoord();
		pid[0] = points3D->InsertNextPoint(Point.x(), Point.y(), Point.z());
		
		Vertices->InsertNextCell(1,pid);
	}

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->SetPoints(points3D);
	pd->SetVerts(Vertices);
	
	return pd;
}



void ModelFile::ColorByValues(const vector<double> &Values)
{
	vector<Color<unsigned char> > Colors(NumPoints());
	
	//color the points by value
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		Colors[i] = ColorByValue(Values[i]);
		Points_[i].setColor(Colors[i]);
	}
	
	//setColors(Colors);
}

void ModelFile::RemoveGround(double MaxAngle)
{
	//angle in degrees
	vgl_vector_3d<double> Up(0.0, 0.0, 1.0);

	//unsigned int NumUpPoints = 0;
	vector<bool> Remove(Points_.size(), false);
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		double angle = rad2deg(AngleBetween(Points_[i].getNormal(), Up));
		
		if(fabs(angle) < MaxAngle)
		{
			//Tools::RemoveVectorElement(Points_, i);
			Remove[i] = true;
			//NumUpPoints++;
		}
	}

	vector<double> z;
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		if(!Remove[i]) //only consider Up points
		{
			z.push_back(Points_[i].getCoord().z());
		}
	}

	double MedianZ = Tools::VectorMedian(z);
	cout << "Median Z: " << MedianZ << endl;
	cout << "Average Z: " << Tools::VectorAverage(z) << endl;


	/*
	// Ground must be parallel to z = 0 (orthogonal to z axis)!
	double eps = 1e-3;
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		if(!Remove[i]) // if it is not an up point
		{
			if(Points_[i].getCoord().z() > (MedianZ + eps) ) //and it is above the "ground"
				NewPoints.push_back(Points_[i]);
		}
		
	}
	*/
	
	//Points_ = NewPoints;
	Triangles_.clear();
	VertexLists_.clear();
}

void ModelFile::ProjectOntoXY(void)
{
	vgl_vector_3d<double> Up(0.0, 0.0, 1.0);
	vector<OrientedPoint> NewPoints;
	vgl_plane_3d<double> GroundPlane(Up, vgl_point_3d<double>(0.0, 0.0, 0.0));
	
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		Ray R(Points_[i].getCoord(), Up);
		vgl_point_3d<double> Intersection;
		R.IntersectPlane(GroundPlane, Intersection);
		Points_[i].setCoord(Intersection);
	}

			
}

void ModelFile::ComputeNormals(double radius)
{
	vgl_point_3d<double> Center = getCenterOfMass();
	
#pragma omp parallel for
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		Tools::PercentComplete(i, NumPoints(), 1);
		
		vgl_point_3d<double> CurrentPoint = getCoord(i);
		vgl_sphere_3d<double> Sphere(CurrentPoint, radius);
		
		vector<vgl_point_3d<double> > PointsInSphere = geom::GetPointsInsideSphere(getCoords(), Sphere);
		
		vgl_vector_3d<double> CurrentDir = CurrentPoint - Center;
		vgl_vector_3d<double> N = geom::BestPlane(PointsInSphere).normal();
		
		
		//for "compact" models
		if(AngleBetween(CurrentDir, N) > deg2rad(90.0))
			N *= -1.0;
		
		Points_[i].setNormal(N);
	}
	
	//vector<bool> Used(Points_.size(), false);
	vector<unsigned int> Used(Points_.size(), 0);
	//FlipNormal(0, Used);
}

void ModelFile::FlipNormal(unsigned int index, vector<unsigned int> &Used)
{
	//this function operates directly on this->Points_
	cout << Tools::VectorSum(Used) << " out of " << Used.size() << endl;
	
	vector<vgl_point_3d<double> > P;
	for(unsigned int i = 0; i < Points_.size(); i++)
	{
		if(i == index)
			continue;
		
		if(!Used[i])
			P.push_back(Points_[i].getCoord());
	}
	
	vector<unsigned int> Inds = KDFuncs::KClosestPointIndices(5, Points_[index].getCoord(), P);
	
	//flip normals for neighbors
	for(unsigned int i = 0; i < Inds.size(); i++)
	{
		if(rad2deg(AngleBetween(Points_[index].getNormal(), Points_[Inds[i]].getNormal())) > 90)
		{
			Points_[Inds[i]].setNormal(-1.0*Points_[Inds[i]].getNormal());
		}
	}
	
	
	//Used[index] = true;
	Used[index] = 1;
	
	//P.clear(); //so the RAM doesn't fill up!!
	vector<vgl_point_3d<double> >().swap(P);
	
	if(Inds.size() == 0) //there are no more points
		return;
	else
	{
		for(unsigned int i = 0; i < Inds.size(); i++)
		{
			FlipNormal(Inds[i], Used);
		}
	}
}


double ModelFile::MeshResolution()
{
	//there must be triangles!
	if(NumTriangles() == 0)
	{
		cout << "No triangles! Mesh resolution is not valid!" << endl;
		return 0.0;
	}
//algorithm:	
//for each point:
//find all points that share an edge with the point
//average those edge lengths; this is the average edge length for this point

//average the average edge lengths

	vector<Edge> Edges;
	//create a list of all the edges
	for(unsigned int tri = 0; tri < NumTriangles(); tri++)
	{
		Triangle T = getTriangle(tri);
		Edge E1(T.getVertex(0), T.getVertex(1));
		Edge E2(T.getVertex(1), T.getVertex(2));
		Edge E3(T.getVertex(2), T.getVertex(0));
		//the add edge function only adds the edge if it is not already in the vector
		AddEdge(Edges,E1);
		AddEdge(Edges,E2);
		AddEdge(Edges,E3);
	}

	double TotalAverageLength = 0.0;

	for(unsigned int point = 0; point < NumPoints(); point++)
	{
		double TotalLength = 0.0;
		unsigned int NumNeighbors = 0;
		for(unsigned int edge = 0; edge < Edges.size(); edge++)
		{
			if(Edges[edge].ContainsPoint(getCoord(point)) )
			{
				TotalLength += Edges[edge].Length();
				NumNeighbors++;
			}
		}
		
		TotalAverageLength += TotalLength/static_cast<double>(NumNeighbors);
	}
	
	TotalAverageLength /= static_cast<double>(NumPoints());
	return TotalAverageLength;
}

double ModelFile::GetAverageEdgeLength()
{
	//there must be triangles!
	if(NumTriangles() == 0)
	{
		cout << "No triangles! Mesh resolution is not valid!" << endl;
		return 0.0;
	}
	
	vector<Edge> Edges;
	//create a list of all the edges
	for(unsigned int tri = 0; tri < NumTriangles(); tri++)
	{
		Triangle T = getTriangle(tri);
		Edge E1(T.getVertex(0), T.getVertex(1));
		Edge E2(T.getVertex(1), T.getVertex(2));
		Edge E3(T.getVertex(2), T.getVertex(0));
		AddEdge(Edges,E1);
		AddEdge(Edges,E2);
		AddEdge(Edges,E3);
	}

	double TotalLength = 0.0;
	for(unsigned int i = 0; i < Edges.size(); i++)
	{
		TotalLength += Edges[i].Length();
	}

	double AverageLength = TotalLength / static_cast<double>(Edges.size());
	return AverageLength;
}

void ModelFile::ColorRed()
{
	for(unsigned int i = 0; i < NumPoints(); i++)
		Points_[i].setColor(Colors::Red());
}

void ModelFile::ColorWhite()
{
	for(unsigned int i = 0; i < NumPoints(); i++)
		Points_[i].setColor(Colors::White());
}

void ModelFile::Downsample(const double MinDist)
{
	//USE MASK POINTS FILTER IN PARAVIEW INSTEAD!

	// The problem with this function is that it deletes both points in a pair that is close enough together.
	// It should be modified to delete only one of these points.
	
	std::vector<bool> PointsToDelete(NumPoints(), false);
	
	KDTree Tree(getCoords());

	unsigned int counter = 0;
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		double dist = Tree.ClosestNonZeroDistance(getCoord(i));
		if(dist < MinDist)
		{
			PointsToDelete[i] = true;
			counter++;
		}
	}

	cout << "There were " << counter << " points marked for removal." << endl;
	
	vector<OrientedPoint> NewPoints;
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		if(PointsToDelete[i] != true)
			NewPoints.push_back(this->Points_[i]);
	}

	setPoints(NewPoints);
}

void ModelFile::SnapToPoints(const vector<vgl_point_3d<double> > &Points)
{
	KDTree Tree(Points);
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		vgl_point_3d<double> CP = Tree.ClosestPoint(Points_[i].getCoord());
		Points_[i].setCoord(CP);
	}
}



bool ModelFile::AllNormalsValid()
{
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		vgl_vector_3d<double> n = Points_[i].getNormal();
		if((n.x() == 0.0) && (n.y() == 0.0) && (n.z() == 0.0))
			return false;
	}
	
	return true;
}

bool ModelFile::GetDoubleArray(const std::string &ArrayName, std::vector<double> &Array)
{
	std::map<std::string,std::vector<double> >::iterator iter;

	iter = DoubleArrays.find(ArrayName);
	if(iter == DoubleArrays.end())
		return false; //array not found!
	else
	{
		Array = iter->second;
		return true;
	}

}