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

#include <cmath>
#include <cstdio>
#include <time.h>
#include <sys/time.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "vtkPolyData.h"
#include <vtkSmartPointer.h>
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkPoints.h"
#include "vtkPolyDataReader.h"
#include <vtkCellArray.h>
#include <vtkDataArray.h>

#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "itksys/SystemTools.hxx"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
	if (argc < 6 )
	{
		std::cerr << "Usage: " << argv[0] << " FiberBundle dimension lambdaW Lambda_Start Lambda_End" << std::endl;
		return -1;
	}

	// Code to see how many threads there in the PC
	int nThreads, tid;
	#pragma omp parallel private(tid)
	{
		tid = omp_get_thread_num();
		if (tid == 0)
		{
			nThreads = omp_get_num_threads();
			printf("Total number of threads: %d\n", nThreads);
		}
	}

	omp_set_num_threads(nThreads); // Code uses all threads available

	// Reading Parameters
	char* Bundle_name = argv[1];
	int Dimension = atoi(argv[2]);
	double lambdaW = atof(argv[3]);
	double lambdaA = atof(argv[4]);
	double lambdaB = atof(argv[5]);

	cout << "Bundle to analyse: " << Bundle_name << endl;
	cout << "Lambda geometry: " << lambdaW << "\nLambda Start (basal): " << lambdaA << "\nLambda End (cortical): " << lambdaB << endl;

	//Creating polydata
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(Bundle_name);
	reader->Update();
	vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

	// Initialisation variables
	vector<vector<pair<int,float> > > links;
	int NumFibers;
	int NumberOfPoints;
	vtkIdType *indices;
	vtkIdType numberOfPoints;
	unsigned int lineCount;
	vtkCellArray* Lines;
	VectorXi NumberPointsPerFiber;
	VectorXf WeightFibers;
	MatrixXf Points;
	vector< MatrixXf > ListPointsFibers;
	vector< MatrixXf > Centers;
	vector< MatrixXf > Tangents;
	MatrixXf First;
	MatrixXf Last;
	MatrixXf c1; // Matrix of double
	MatrixXf t1;
	VectorXf f1; // Vector of double
	VectorXf l1;
	MatrixXf c2;
	MatrixXf t2;
	VectorXf f2;
	VectorXf l2;
	RowVectorXf point;
	RowVectorXf p0;
	RowVectorXf p1;

	// Reading the polydata
	NumFibers = polyData->GetNumberOfLines();
	cout << "Number fibers: " << NumFibers << endl;

	links.resize(NumFibers);

	NumberOfPoints= polyData->GetNumberOfPoints();
	cout << "Number points: " << NumberOfPoints << endl;

	Lines = polyData->GetLines();

	NumberPointsPerFiber.resize(NumFibers);
	WeightFibers.resize(NumFibers);

	lineCount = 0;
	for (Lines->InitTraversal(); Lines->GetNextCell(numberOfPoints, indices); lineCount++)
	{
		NumberPointsPerFiber(lineCount)=numberOfPoints;
	}
	//cout << "NumberPointsPerFiber: " << NumberPointsPerFiber << endl;

	if( NumberPointsPerFiber.sum() != NumberOfPoints )
		throw runtime_error("Total Number of Points Mismatched!");

	Points.resize(NumberOfPoints,Dimension);

	for (unsigned int i = 0; i < NumberOfPoints; i++)
	{
		 double p[3];
		 polyData->GetPoint(i, p);
		 for (int dim = 0; dim < Dimension; dim++)
		 {
			float pf = (float)(p[dim]);
		 	Points(i, dim) = pf;
		 }
	}

	// List Points per Fiber
	unsigned int temp = 0;
	ListPointsFibers.resize(NumFibers);

	for (unsigned int lineCount = 0; lineCount<NumFibers; lineCount++)
	{
		ListPointsFibers[lineCount].resize(NumberPointsPerFiber(lineCount), Dimension);
		for (unsigned int i = 0; i < NumberPointsPerFiber(lineCount,0); i++)
		{
			point = Points.row(temp);
			ListPointsFibers[lineCount].row(i)=point;
			temp++;
		}
	}

	// Computing centers, tangents, first and last points
	Centers.resize(NumFibers);
	Tangents.resize(NumFibers);
	First.resize(NumFibers,Dimension);
	Last.resize(NumFibers,Dimension);

	for (unsigned int i = 0; i < NumFibers; i++)
	{
		Centers[i].resize(NumberPointsPerFiber(i)-1, Dimension);
		Tangents[i].resize(NumberPointsPerFiber(i)-1, Dimension);

		First.row(i)=ListPointsFibers[i].row(0);
		Last.row(i)=ListPointsFibers[i].row(NumberPointsPerFiber(i)-1);

		for (unsigned int j = 0; j < NumberPointsPerFiber(i)-1; j++)
		{
			p0 = ListPointsFibers[i].row(j);
			p1 = ListPointsFibers[i].row(j+1);
			Centers[i].row(j) = (p0 + p1) / 2.0;
			Tangents[i].row(j) = p1 - p0;
		}
	}

	struct timeval start, end;
	double delta;

	VectorXf Diagonal;
	Diagonal.resize(NumFibers);

	// Multi-Threading
	gettimeofday(&start, NULL);

	for(unsigned int i=0; i<NumFibers; i++)
	{
		c1.setZero(NumberPointsPerFiber(i)-1, Dimension);
		t1.setZero(NumberPointsPerFiber(i)-1, Dimension);
		f1.setZero(Dimension);
		l1.setZero(Dimension);

		c1 = Centers[i];
		t1 = Tangents[i];
		f1 = First.row(i);
		l1 = Last.row(i);

		#pragma omp parallel for private(c2,t2,f2,l2) shared(links,Diagonal,i,t1,c1,f1,l1,lambdaW,lambdaA,lambdaB)
		for(unsigned int j=i; j<NumFibers; j++)
		{
			c2.setZero(NumberPointsPerFiber(j)-1, Dimension);
			t2.setZero(NumberPointsPerFiber(j)-1, Dimension);
			f2.setZero(Dimension);
			l2.setZero(Dimension);

			c2 = Centers[j];
			t2 = Tangents[j];
			f2 = First.row(j);
			l2 = Last.row(j);

			// Computation norm usual currents
			float norm2 = 0;
			float res_tang;
			float res_center;
			for (int p=0; p<NumberPointsPerFiber(i)-1; p++)
			{
				for (int q=0; q<NumberPointsPerFiber(j)-1; q++)
				{
					res_tang = t1.row(p)*t2.row(q).transpose();
					res_center = ( c1.row(p)-c2.row(q) ) * ( (c1.row(p)-c2.row(q)).transpose() );
					norm2 = norm2+res_tang*exp(-res_center/(lambdaW*lambdaW));
				}
			}

			// Computation of the other two kernels, only if the norm of usual currents
			// is greater than 1e-7, otherwise it writes 0
			float norm2_f = 0;
			float res_f;
			float res_l;
			if (abs(norm2)>1e-7)
			{
				res_f = (f1-f2).transpose()*(f1-f2);
				res_l = (l1-l2).transpose()*(l1-l2);
				norm2_f = norm2 * exp(-res_f/(lambdaA*lambdaA) - res_l/(lambdaB*lambdaB));

				// If the norm is smaller than 1e-7, it writes 0
				if (abs(norm2_f)>1e-7)
				{

					if (i==j)
					{
						Diagonal[i]=norm2_f;
					}
					else
					{
						#pragma omp critical // This is important, otherwise there is an error of Segmentation Fault
						{
							links[i].push_back(make_pair(j,norm2_f));
							links[j].push_back(make_pair(i,norm2_f));
						}
					}
				}
			}

		} // end for j

		if (remainder(i,10000)==0)
			cout << "Iter " << i << endl;

	} // end for i

// Diagonal Element

	ofstream Diag;
  Diag.open("graph.diag", fstream::out | fstream::binary);
	float element = 0;;
	for(unsigned int i=0; i<NumFibers; i++)
	{
		element = Diagonal[i];
		Diag.write((char *)(&element),4);
	}

	Diag.close();

// Gramiam

	ofstream Gram;
  	Gram.open("graph.bin", fstream::out | fstream::binary);

	ofstream Weights;
  	Weights.open("graph.weights", fstream::out | fstream::binary);

	unsigned int s = links.size();
	if (s!=NumFibers)
		cerr << "PROBABLY THERE IS AN ERROR! Number of nodes should be equal to the number of fibers" << endl;

	// outputs number of nodes
	Gram.write((char *)(&s),4);

	// outputs cumulative degree sequence
	long tot=0;
	for (unsigned int i=0 ; i<s ; i++) {
		tot+=(long)links[i].size();
		Gram.write((char *)(&tot),8);
	}

	// outputs links
	for (unsigned int i=0 ; i<s ; i++) {
		for (unsigned int j=0 ; j<links[i].size() ; j++)
		{
			int dest = links[i][j].first;
			float weight = links[i][j].second;
			Gram.write((char *)(&dest),4);
			Weights.write((char *)(&weight),4);
		}
	}

	Gram.close();
	Weights.close();

// TIMER
	gettimeofday(&end, NULL);
	delta = double(end.tv_sec  - start.tv_sec) + double(end.tv_usec - start.tv_usec) / 1.e6;
	printf ("It took %f seconds \n",delta);

	return 0;
}
