
#include "VisualizableModel.h"
#include <GL/glut.h>

void VisualizableModel::DrawPoints() const
{
	std::clog << "Drawing " << NumPoints() << "points with " << NumColors() << " colors." << std::endl;
	
	glBegin(GL_POINTS);
	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		vgl_point_3d<double> P = Points_[i].getCoord();
		Color<unsigned char> cu = Points_[i].getColor();

		glColor3ub(cu.getR(), cu.getG(), cu.getB());
		glVertex3f(P.x(), P.y(), P.z());
			
	}
	glEnd();
}

void VisualizableModel::DrawTriangles(void) const
{
	if(NumTriangles() > 0)
	{
		glBegin(GL_TRIANGLES);
	
		for(unsigned int i = 0; i < NumTriangles(); i++)
		{
			Color<unsigned char> cu = Points_[i].getColor();
			glColor3ub(cu.getR(), cu.getG(), cu.getB());
				//Triangles_[i].Draw(Normals, true);
			Triangles_[i].Draw(false, true);
		}

		glEnd();
	}
	else
	{
		cout << "NO TRIANGLES!" << endl;
	}

}


void VisualizableModel::DrawQuads(void) const
{
	glBegin(GL_QUADS);

	for(unsigned int i = 0; i < NumPoints(); i++)
	{
		vgl_vector_3d<double> N = Points_[i].getNormal();
		vgl_vector_3d<double> V1 = geom::GetOrthogonalVector(N);
		vgl_vector_3d<double> V2 = cross_product(V1, N);

		glNormal3d(N.x(), N.y(), N.z());
		
		vgl_point_3d<double> P0 = Points_[i].getCoord() + .5 * V1;
		vgl_point_3d<double> P1 = Points_[i].getCoord() - .5 * V1;
		vgl_point_3d<double> P2 = Points_[i].getCoord() + .5 * V2;
		vgl_point_3d<double> P3 = Points_[i].getCoord() - .5 * V2;
		
		glVertex3d(P0.x(), P0.y(), P0.z());
		glVertex3d(P1.x(), P1.y(), P1.z());
		glVertex3d(P2.x(), P2.y(), P2.z());
		glVertex3d(P3.x(), P3.y(), P3.z());
		
	}

	glEnd();

}


void ModelFile::LookAt(const vgl_point_3d<double> &CameraPosition) const
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	vgl_point_3d<double> Closest = geom::ClosestPoint(CameraPosition, getCoords());
	vgl_point_3d<double> Farthest = geom::FarthestPoint(CameraPosition, getCoords());
	vgl_point_3d<double> Center = geom::CenterOfMass(getCoords());
	
	double znear = vgl_distance(CameraPosition, Closest);
	double zfar = vgl_distance(CameraPosition, Farthest);
	
	//gluPerspective(70, 1, 1, 10); //view angle in y direction, aspect (width to height), zNear, zFar
	gluPerspective(70, 1, znear, zfar); //view angle in y direction, aspect (width to height), zNear, zFar
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	//gluLookAt(2.0, 2.0, 10.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0); //eyex, eyey, eyez, centerofviewx, centerofviewy, centerofviewz, upx, upy, upz
	gluLookAt(CameraPosition.x(), CameraPosition.y(), CameraPosition.z(), Center.x(), Center.y(), Center.z(), 0.0, 0.0, 1.0);

}