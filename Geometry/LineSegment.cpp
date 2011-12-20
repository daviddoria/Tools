#include "LineSegment.h"

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include <vgl/vgl_point_3d.h>

void WriteLineFile(const std::vector<LineSegment> &LineSegments, const std::string &Filename)
{
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

	/*
	//setup colors
	vtkSmartPointer<vtkUnsignedCharArray> Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	Colors->SetNumberOfComponents(3);
	Colors->SetName("Colors");

	unsigned char r[3];
	Color<unsigned char> R = Red();
	r[0] = R.getR();
	r[1] = R.getG();
	r[2] = R.getB();
	Colors->InsertNextTupleValue(r);
	*/

	//add points
	for(unsigned int i = 0; i < LineSegments.size(); i++)
	{
		for(unsigned int j = 0; j < 2; j++)
		{
			double p0[3];
			p0[0] = LineSegments[i].EndPoints[j].x();
			p0[1] = LineSegments[i].EndPoints[j].y();
			p0[2] = LineSegments[i].EndPoints[j].z();
			
			pts->InsertNextPoint(p0);
		}

	}

	//add lines (connect every other point)
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

	for(unsigned int i = 0; i < LineSegments.size() * 2; i+=2)
	{
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		
		//x axis
		line->GetPointIds()->SetId(0,i);
		line->GetPointIds()->SetId(1,i+1);
		lines->InsertNextCell(line);

	}
	
	vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();

	//add the points to the dataset
	pdata->SetPoints(pts);

	//add the lines to the dataset
	pdata->SetLines(lines);

	//color the lines
	//pdata->GetCellData()->AddArray(Colors);

	//write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInput(pdata);

	writer->SetFileName(Filename.c_str());
	writer->Write();
}