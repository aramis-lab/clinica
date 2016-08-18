// Created on 17/08/2016 by Pietro Gori, Inria
//
// C++ function that takes as input three ASCII text files containing
// information about the weighted prototypes and it resturns a VTK file format
// with tubes instead than lines.
//
// Usage: WriteTube Points_filename Number_Points_Per_Curve_filename Radius_filename outfilename
//
// Input Parameters:
//	- Points_filename: text ASCII file with the coordinates of the points of the
//                     weighted prototypes.
//	- Number_Points_Per_Curve_filename: text ASCII file with the number of points
//         															per prototype
//	- Radius_filename: text ASCII file with the weights of each prototype
//	- outfilename: name of the VTK tube file
//
// Outputs:
//  - the VTK tube file

#include "vtkLine.h"
#include "vtkTubeFilter.h"
#include "vtkLineSource.h"
#include "vtkAppendPolyData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "itkVersion.h"
#include "vtkVersion.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

#include "cmath"
#include "stdio.h"
#include "iostream"
#include "vector"

typedef vnl_vector<double> VNLVectorType;
typedef vnl_matrix<double> VNLMatrixType;

int main(int argc, char** argv)
{

	if (argc < 4)
	{
		std::cerr << "Usage: " << argv[0] << " Points_filename Number_Points_Per_Curve_filename Radius_filename outfilename" << std::endl;
		return -1;
	}

	// Get Points
	std::string Points_filename = argv[1];
	std::ifstream finP(Points_filename.c_str());

	std::string lineP;
	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();

	while(std::getline(finP, lineP))
	{
		double x,y,z;
		std::stringstream linestream;
		linestream << lineP;
		linestream >> x >> y >> z;

		Points->InsertNextPoint(x, y, z);
	}

	finP.close();

	// Get Number_Points_Per_Curve
	std::vector<int> Number_Points_Per_Curve;

	std::string Number_filename = argv[2];
	std::ifstream finN(Number_filename.c_str());

	std::string lineN;

	while(std::getline(finN, lineN))
	{
		double x;
		std::stringstream linestream;
		linestream << lineN;
		linestream >> x;

		Number_Points_Per_Curve.push_back(x);
	}

	finN.close();

	// Get Radius
	std::vector<double> Radius;

	std::string Radius_filename = argv[3];
	std::ifstream finR(Radius_filename.c_str());

	std::string lineR;

	while(std::getline(finR, lineR))
	{
		double x;
		std::stringstream linestream;
		linestream << lineR;
		linestream >> x;

		Radius.push_back(x);
	}
	finR.close();

	// PARAMETERS
	char* outfilename = argv[4];
	int NPts = Points->GetNumberOfPoints ();
	int NMedoids = Number_Points_Per_Curve.size();

	if (NMedoids != Radius.size())
	{
		std::cerr << "ERROR! Size not equal!" << std::endl;
		return -1;
	}

	// Start appending Fibers
	vtkSmartPointer<vtkAppendPolyData> appendAll = vtkSmartPointer<vtkAppendPolyData>::New();

	int k = 0;

	for (int i=0;i<NMedoids;i++)
	{

		vtkSmartPointer<vtkPoints> points_temp = vtkSmartPointer<vtkPoints>::New();

		for (int j=0;j<Number_Points_Per_Curve[i];j++)
		{
	    	double Point_temp[3];
			  Points->GetPoint(k,Point_temp);
	    	points_temp->InsertNextPoint(Point_temp);
	    	k++;
		}

		// Create a line
		vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
		line->SetPoints(points_temp);

		// Create a tube
		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
		tubeFilter->SetInputConnection(line->GetOutputPort());
		tubeFilter->SetRadius(Radius[i]);
		tubeFilter->SetNumberOfSides(10);
		tubeFilter->Update();

		#if VTK_MAJOR_VERSION < 6
			appendAll->AddInput(tubeFilter->GetOutput());
		#else
			appendAll->AddInputData(tubeFilter->GetOutput());
		#endif
	}

	appendAll->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(outfilename);
	#if VTK_MAJOR_VERSION < 6
		writer->SetInput(appendAll->GetOutput());
	#else
		writer->SetInputData(appendAll->GetOutput());
	#endif
	writer->Update();
}
