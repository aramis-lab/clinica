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
#include <string>
#include <stdlib.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>

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

#include "itksys/SystemTools.hxx"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Eigen;


int main(int argc, char** argv)
{
	if (argc < 5 )
	{
		cerr << "Usage: " << argv[0] << " PointsFibers PointsMedoids Tau Surface_filename lambda" << endl;
		return -1;
	}

	// Time
	struct timeval start, end;
	double delta;
	gettimeofday(&start, NULL);


	// Parameters to see how many threads there are
	Eigen::setNbThreads(0);; // Otherwise Eigen use OpenMP

	int nThreads, tid;
	#pragma omp parallel private(tid)
	{
		tid = omp_get_thread_num();

		if (tid == 0)
		{
			nThreads = omp_get_num_threads();
			//printf("Total number of threads: %d\n", nThreads);
		}
	}

	omp_set_num_threads(nThreads); // in teoria non serve

	// Reading Parameters
	char* PointsFibersFile = argv[1];
	char* PointsMedoidsFile = argv[2];
	char* TauFile = argv[3];
	char* SurfaceFile = argv[4];
	double lambda = atof(argv[5]);

	//cout << "Filename Points Fibers: " << PointsFibersFile << endl;
	//cout << "Filename Index Points: " << IndexMedoidsFile << endl;
	//cout << "Filename Tau: " << TauFile << endl;
	//cout << "Filename Surface: " << SurfaceFile << endl;
	//cout << "lambda: " << lambda << endl;

// Parameters
	int NumberPointsSurface;
	unsigned int NFibers;
	MatrixXf PointsFibers;	// [NFibersx3]
	VectorXf ExactTau;
	MatrixXf PointsMedoids;	// [NMedoidsx3]
	unsigned int NMedoids;
	MatrixXf GammaB;
	double PI = 3.14159265358979323846;
	VectorXf DensityB;
	VectorXf DensityM;
	MatrixXf GammaM;
	MatrixXf PointCoordinates;
	MatrixXf PointsFibersTemp;
	MatrixXf PointsMedoidsTemp;


// Surface

	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(SurfaceFile);
	reader->Update();
	vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

	NumberPointsSurface = polyData->GetNumberOfPoints();

	PointCoordinates.setZero(NumberPointsSurface,3);
	for (unsigned int i = 0; i < NumberPointsSurface; i++)
	{
		 double p[3];
		 polyData->GetPoint(i, p);
		 for (int dim = 0; dim < 3; dim++)
		 	PointCoordinates(i, dim) = p[dim];
	}

	//int NumFaces = polyData->GetNumberOfCells();

//cout << "PointCoordinates: \n" << PointCoordinates << "\n" << endl;

// Coordinates Points Bundle

	ifstream PointsFibersStream;
	PointsFibersStream.open(PointsFibersFile,fstream::in | fstream::binary);
	if ( !PointsFibersStream.is_open() )
	{
		cerr << "Error! Could not open Fibers file!" << endl;
		return -1;
	}

	// Read number of nodes on 4 bytes
	PointsFibersStream.read((char *)&NFibers, 4);
	assert(PointsFibersStream.rdstate() == ios::goodbit);

	PointsFibers.setZero(1,NFibers*3);

	PointsFibersStream.read((char *)&PointsFibers(0,0), NFibers*3*4);
	PointsFibersStream.close();

//cout << "PointsFibers: \n" << PointsFibers << "\n" <<  endl;

	PointsFibersTemp.resize(3,NFibers);
	PointsFibers.resizeLike(PointsFibersTemp);
	PointsFibers.transposeInPlace();
	PointsFibersTemp.resize(1,1);

//cout << "PointsFibers: \n" << PointsFibers << "\n" <<  endl;

// Coordinates Points Medoids

	ifstream PointsMedoidsStream;
	PointsMedoidsStream.open(PointsMedoidsFile,fstream::in | fstream::binary);
	if ( !PointsMedoidsStream.is_open() )
	{
		cerr << "Error! Could not open Medoids file!" << endl;
		return -1;
	}

	// Read number of nodes on 4 bytes
	PointsMedoidsStream.read((char *)&NMedoids, 4);
	assert(PointsMedoidsStream.rdstate() == ios::goodbit);

	PointsMedoids.setZero(1,NMedoids*3);

	PointsMedoidsStream.read((char *)&PointsMedoids(0,0), NMedoids*3*4);
	PointsMedoidsStream.close();

//cout << "PointsMedoids: \n" << PointsMedoids << "\n" <<  endl;

	PointsMedoidsTemp.resize(3,NMedoids);
	PointsMedoids.resizeLike(PointsMedoidsTemp);
	PointsMedoids.transposeInPlace();
	PointsMedoidsTemp.resize(1,1);

//cout << "PointsMedoids: \n" << PointsMedoids << "\n" <<  endl;

// Tau
	ifstream TauStream;
	TauStream.open(TauFile,fstream::in | fstream::binary);
	if ( !TauStream.is_open() )
	{
		cerr << "Error! Could not open Tau file!" << endl;
		return -1;
	}

	ExactTau.setZero(NMedoids);
	TauStream.read((char *)&ExactTau(0), NMedoids*4);
	TauStream.close();

//cout << "ExactTau: \n" << ExactTau << "\n" <<  endl;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// CODE ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
	// Gamma Fibers
	try{
		GammaB.setZero(NumberPointsSurface,NFibers);
	}
	catch (const std::bad_alloc&){
		cerr << "Error to allocate memory! Fiber bundle too big!" << endl;
		return -1;
	}

	// Gamma Fiber Bundle
	for(unsigned int dim=0; dim<3; dim++)
	{
		MatrixXf MV;
		MV.setZero(NumberPointsSurface,NFibers);
		#pragma omp parallel for shared(MV,PointCoordinates)
		for (unsigned int i=0; i<NFibers; i++)
		{
			MV.col(i)=PointCoordinates.col(dim);
		}

		#pragma omp parallel for shared(MV,PointsFibers)
		for (unsigned int i=0; i<NumberPointsSurface; i++)
		{
			MV.row(i) = MV.row(i) - PointsFibers.col(dim).transpose();
		}

		MV=MV.array().pow(2);
		GammaB=GammaB+MV;
	}

	PointsFibers.resize(1,1);

	GammaB=-GammaB/(2*lambda*lambda);
	GammaB=GammaB.array().exp();
	GammaB=GammaB*( ( 1/sqrt(8*PI*PI*PI) ) * ( 1/(lambda*lambda*lambda) ) );

	DensityB.setZero(NumberPointsSurface);
	DensityB=(GammaB.rowwise().sum())/NFibers;

	GammaB.resize(1,1);

	// Gamma Medoids

	try{
		GammaM.setZero(NumberPointsSurface,NMedoids);
	}
	catch (const std::bad_alloc&){
		cerr << "Error to allocate memory! Too many medoids!" << endl;
		return -1;
	}

	for(unsigned int dim=0; dim<3; dim++)
	{
		MatrixXf MV;
		MV.setZero(NumberPointsSurface,NMedoids);
		#pragma omp parallel for shared(MV,PointCoordinates)
		for (unsigned int i=0; i<NMedoids; i++)
		{
			MV.col(i)=PointCoordinates.col(dim);
		}

		#pragma omp parallel for shared(MV,PointsFibers)
		for (unsigned int i=0; i<NumberPointsSurface; i++)
		{
			MV.row(i) = MV.row(i) - PointsMedoids.col(dim).transpose();
		}

		MV=MV.array().pow(2);
		GammaM=GammaM+MV;
	}

	PointCoordinates.resize(1,1);
	PointsMedoids.resize(1,1);

	GammaM=-GammaM/(2*lambda*lambda);
	GammaM=GammaM.array().exp();
	GammaM=GammaM*( ( 1/sqrt(8*PI*PI*PI) ) * ( 1/(lambda*lambda*lambda) ) );

	#pragma omp parallel for shared(GammaM,ExactTau)
	for (unsigned int i=0; i<NumberPointsSurface; i++)
	{
		//GammaM.row(i)=GammaM.row(i).array() * ExactTau.transpose().array();
		GammaM.row(i)=GammaM.row(i).cwiseProduct(ExactTau.transpose());
	}

	DensityM.setZero(NumberPointsSurface);
	DensityM=(GammaM.rowwise().sum())/(ExactTau.sum());

	GammaM.resize(1,1);

	cout << "Number fibers: " << NFibers << ", sum tau: " << ExactTau.sum() << endl;

	ExactTau.resize(1);

	ofstream DensityFibersBin;
	DensityFibersBin.open("DensityFibers.bin", fstream::out | fstream::binary);
	for (unsigned int j=0 ; j<NumberPointsSurface ; j++) {
      		float fib = DensityB[j];
      		DensityFibersBin.write((char *)(&fib),sizeof(float));
    	}
  	DensityFibersBin.close();

	ofstream DensityMedoidsBin;
	DensityMedoidsBin.open("DensityMedoids.bin", fstream::out | fstream::binary);
	for (unsigned int j=0 ; j<NumberPointsSurface ; j++) {
      		float med = DensityM[j];
      		DensityMedoidsBin.write((char *)(&med),sizeof(float));
    	}
  	DensityMedoidsBin.close();

// TIMER
	gettimeofday(&end, NULL);
	delta = double(end.tv_sec  - start.tv_sec) + double(end.tv_usec - start.tv_usec) / 1.e6;
	printf ("It took %f seconds \n",delta);

	return 0;

} // end main
