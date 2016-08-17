// Created on 17/08/2016 by Pietro Gori, Inria
//
// C++ function that takes as input a fiber bundle in .VTK format and computes the
// gramiam matrix between all streamlines. Every cell of the gramiam matrix is
// the inner product between two streamlines in the framework of weighted
// currents.
//
// Usage: Gramiam FiberBundle dimension lambdaW Lambda_Start Lambda_End
//
// Input Parameters:
//	- FiberBundle: filename fiber bundle in .vtk format
//	- dimension: dimension of the points of the stremlines (i.e. 3 for 3D)
//	- lambdaW: bandwidth of the geometric kernel of usual currents
//	- Lambda_Start: bandwidth of kernel relative to the starting points
//	- Lambda_End: bandwidth of kernel relative to the ending points
// To note: Streamlines must have a consistent orientation ! For instance they
// should all have the same starting and ending ROIs (Region Of Interest)
//
// Outputs:
// 3 binary files, let N be equal to the number of streamlines
//	- graph.diag: it is a vector [Nx1] with the squared norm of each streamline.
//                Every value is saved as a char of 4 bits
//	- graph.bin: It is a vector of char. If first writes the number of Nodes
//	            (i.e. number fo streamlines) as a char of 4 bits. Then it writes
//							the cumulative degree sequence, which means that for each
//							streamline i it writes the number of streamlines that have an
//							inner product greater than 0 as a char of 8 bits. Then it writes
//							the numbers of all these streamlines as a char of 4 bits.
//	- graph.weights: A vector with the inner products different from 0 between
//									 the streamlines. They are chars of 4 bits. The squared norm
//                   of each streamline is not considered.
//
// To note, this is the style accepted in the function community.
//
// Example: Gram matrix is [2 0 2; 3 4 6; 0 0 2].
//	graph.diag contains: 2 4 2 (squared norms, diagonal)
//	graph.bin contains: 3 (number streamlines) 1 2 0 (number entries different from 0)
//											2 (last column) 0 2 (first and last columns) (no value there are only zeros)
//  graph.weights contains: 2 3 6 (the inner products different from zero)

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

	//std::cout<<"Number 1 = "<<Number_Points_Per_Curve[0]<<std::endl;

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


	//std::cout<<"Radius 1 = "<<Radius[0]<<std::endl;
    //std::cout<<"Radius 2 = "<<Radius[1]<<std::endl;

	// OUTFILENAME
	char* outfilename = argv[4];


	// PARAMETERS
	int NPts = Points->GetNumberOfPoints ();
	//std::cout<<"NPts = "<<NPts<<std::endl;
	int NMedoids = Number_Points_Per_Curve.size();
	//std::cout<<"NMedoids = "<<NMedoids<<std::endl;
    //std::cout<<"Radius size = "<<Radius.size()<<std::endl;

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
			//std::cout<<"Points temp = "<<Point_temp[0]<<", "<<Point_temp[1]<<", "<<Point_temp[2]<<std::endl;
	    	points_temp->InsertNextPoint(Point_temp);
	    	k++;
		}

		// Create a line
		vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
		line->SetPoints(points_temp);

		//vtkSmartPointer<vtkPolyDataWriter> writer1 = vtkSmartPointer<vtkPolyDataWriter>::New();
		//writer1->SetFileName("line.vtk");
		//writer1->SetInput(line->GetOutput()); // devo mettere vtkPolyData *
		//writer1->Update();

		// Create a tube
		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
		tubeFilter->SetInputConnection(line->GetOutputPort());
		tubeFilter->SetRadius(Radius[i]);
		tubeFilter->SetNumberOfSides(10);
		tubeFilter->Update();

		//vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
		//writer2->SetFileName("tube.vtk");
		//writer2->SetInput(tubeFilter->GetOutput());
		//writer2->Update();

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
