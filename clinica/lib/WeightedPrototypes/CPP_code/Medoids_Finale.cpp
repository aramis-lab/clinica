// Created on 17/08/2016 by Pietro Gori, Inria
//
// C++ function that takes as input a Gram matrix and the result of the clustering
// based on modularity and it computes the weighted prototypes as explained in
// "Parsimonious Approximation of Streamline Trajectories in White Matter Fiber Bundles",
// IEEE Transactions on Medical Imaging, 2016
//
// Usage: Medoids_Finale Links Weights Diagonal minValueTau degree_precision outlier_limit Number_Fascicles Index_fibers1 Index_fibers2 ...
//
// Input Parameters:
//	- Links: graph.bin of compute_gramiam.cpp
//	- Weights: graph.weights of compute_gramiam.cpp
//	- Diagonal: graph.diag of compute_gramiam.cpp
//	- minValueTau: We remove the prototypes that approximate less than
//                    minValueTau fibers. Default value is 1
//	- degree_precision: percentage of the norm of the bundle explained by the
//                      weighted prototypes. Default value is 0.15, which means
//                      that the weighted prototypes will explain (1-0.15)*100 %
//                      of the norm of the bundle in the framework of weighted currents
//	- outlier_limit: maximum average angle (in radians) that a streamline may
//                   have with the other streamlines in the framework of weighted
// 									 currents. Default value is 1.5359 = 88 degrees
//	- Number_Fascicles: number of fascicles (clusters)
//	- Index_fibersXXX: for each cluster (fascicle) we have a text file with the indexes
//										 of the streamlines of the bundle belonging to the cluster
//
// Outputs:
// 3 binary files
//	- Outliers_global: indexes of the streamlines of the bundle which are
//										 considered outliers
//	- Medoids_global_normalised: indexes of the streamlines of the bundle which
//															 are considered prototypes
//	- Tau_global_normalised: Weights of the prototypes
//

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

// Remove row and column z from the square matrix Gamma
MatrixXf update_GramMM (MatrixXf Gamma, int z)
{
	z=z+1; // This is the index of the row/column which I want to eliminate
	MatrixXf TempGamma;
	TempGamma.setZero(Gamma.rows()-1,Gamma.cols()-1);
	TempGamma.topLeftCorner(z-1,z-1)=Gamma.topLeftCorner(z-1,z-1);
	TempGamma.topRightCorner(z-1,Gamma.rows()-z)=Gamma.topRightCorner(z-1,Gamma.rows()-z);
	TempGamma.bottomLeftCorner(Gamma.rows()-z,z-1)=Gamma.bottomLeftCorner(Gamma.rows()-z,z-1);
	TempGamma.bottomRightCorner(Gamma.rows()-z,Gamma.rows()-z)=Gamma.bottomRightCorner(Gamma.rows()-z,Gamma.rows()-z);
	return TempGamma;
}

// Remove row z and column Medoids(z)+1 from the rectangular matrix GramiamMN
MatrixXf update_GramMN (MatrixXf GramiamMN, VectorXi Medoids, int z)
{
	int zM=z+1; // This is the index of the row which I want to eliminate
	int zN=Medoids(z)+1; // This is the index of the columns which I want to eliminate
	MatrixXf TempGamma;
	TempGamma.setZero(GramiamMN.rows()-1,GramiamMN.cols()-1);
	TempGamma.topLeftCorner(zM-1,zN-1)=GramiamMN.topLeftCorner(zM-1,zN-1);
	TempGamma.topRightCorner(zM-1,GramiamMN.cols()-zN)=GramiamMN.topRightCorner(zM-1,GramiamMN.cols()-zN);
	TempGamma.bottomLeftCorner(GramiamMN.rows()-zM,zN-1)=GramiamMN.bottomLeftCorner(GramiamMN.rows()-zM,zN-1);
	TempGamma.bottomRightCorner(GramiamMN.rows()-zM,GramiamMN.cols()-zN)=GramiamMN.bottomRightCorner(GramiamMN.rows()-zM,GramiamMN.cols()-zN);
	return TempGamma;
}

// Remove element z from vector OriginalIndex
VectorXi update_vector (VectorXi OriginalIndex, int z)
{
	z=z+1; // Number of element to delete
	VectorXi TempOriginalIndex;
	TempOriginalIndex.setZero(OriginalIndex.size()-1);
	TempOriginalIndex.head(z-1)=OriginalIndex.head(z-1); // I want the first z-1 elements
	TempOriginalIndex.tail(OriginalIndex.size()-z)=OriginalIndex.tail(OriginalIndex.size()-z); // I want the last N-z elements
	return TempOriginalIndex;
}

// Find Prototype in matrix Gamma
int find_medoid (MatrixXf Gamma)
{
	int IF = Gamma.rows();
	VectorXf Degree2;
	Degree2.setZero(IF);

	#pragma omp parallel for shared(Degree2,Gamma)
	for (unsigned int j=0; j<IF; j++)
	{
		Degree2(j) = Gamma.row(j).sum();
		Degree2(j) = Degree2(j) / sqrt(Gamma(j,j));
	}

	Degree2=Degree2.array().square();
	int z;
	double s = Degree2.maxCoeff(&z);
	return z;
}

// Reshape matrix Gamma taking only the rows and columns in Indexa and Indexb respectively
MatrixXf reshape (MatrixXf Matrix, VectorXi Indexa, VectorXi Indexb)
{
	int L1 = Indexa.size();
	int L2 = Indexb.size();
	MatrixXf Result;
	Result.setZero(L1,L2);

	#pragma omp parallel for shared(Matrix,Result,Indexa,Indexb) collapse(2)
	for (unsigned int k=0; k<L1; k++)
	{
		for (unsigned int kk=0; kk<L2; kk++)
		{
			Result(k,kk)=Matrix(Indexa(k),Indexb(kk));
		}
	}

	return Result;
}

// Compute the weight tau for each prototype
VectorXf exact_tau (MatrixXf Gamma, VectorXi Medoids)
{
	MatrixXf TempMM;
	MatrixXf TempM;
	int L = Medoids.size();
	int IF = Gamma.rows();
	VectorXf ExactTau;
	ExactTau.setZero(L);
	VectorXi OriginalIndex;
	OriginalIndex.setZero(IF);
	OriginalIndex.setLinSpaced(IF,0,IF-1);
	TempMM=reshape(Gamma,Medoids,Medoids);
	TempM=reshape(Gamma,Medoids,OriginalIndex);
	VectorXf Ones;
	Ones.setOnes(IF);
	ExactTau=TempMM.ldlt().solve(TempM*Ones);
	return ExactTau;
}

// Compute the cost function
 float Cost_function (float InitNorm2, VectorXf ExactTau, MatrixXf Gamma, VectorXi Medoids)
{
	MatrixXf TempMM;
	MatrixXf TempM;
	int L = Medoids.size();
	int IF = Gamma.rows();

	VectorXi OriginalIndex;
	OriginalIndex.setZero(IF);
	OriginalIndex.setLinSpaced(IF,0,IF-1);

	TempMM=reshape(Gamma,Medoids,Medoids);
	TempM=reshape(Gamma,Medoids,OriginalIndex);

	VectorXf Ones;
	Ones.setOnes(IF);

	float Result = InitNorm2 - 2*ExactTau.transpose()*TempM*Ones + ExactTau.transpose()*TempMM*ExactTau;

	return Result;
}


int main(int argc, char** argv)
{
	if (argc < 8 )
	{
		cerr << "Usage: " << argv[0] << " Links Weights Diagonal minValueTau degree_precision outlier_limit Number_Fascicles Index_fibers1 Index_fibers2 ..." << endl;
		return -1;
	}

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

	omp_set_num_threads(nThreads); // maybe it is not necessary...

	// Reading Parameters
	char* Links_file = argv[1];
	char* Weight_file = argv[2];
	char* Diag_file = argv[3];
	double minValueTau = atof(argv[4]);
	double degreePrecision = atof(argv[5]);
	float outlierlimit = atof(argv[6]);
	int NumberFascicles = atoi(argv[7]);

	if ((argc - 8) != NumberFascicles)
	{
		cerr << "Error! The number of fascicles files is not correct!" << endl;
		return -1;
	}

	std::vector<char*> Fibers_file(NumberFascicles);
	unsigned int temp = 8;
	for (unsigned int f=0; f<NumberFascicles; f++)
	{
		Fibers_file[f] = argv[temp];
		temp++;
	}

// Time
	struct timeval start, end;
	double delta;
	gettimeofday(&start, NULL);

// Parameters
	unsigned int nb_nodes;
	unsigned long nb_links;
	unsigned int IF;
	vector<unsigned long> degrees;
	vector<unsigned int> links;
	VectorXf weights;
	VectorXf diag;
	vector<unsigned int> indexFibers;
	VectorXi indexFibersNoOutliers;
	MatrixXf Gramiam;
	unsigned int NFibers;
	MatrixXf GramTemp;
	VectorXi Medoids;
	VectorXf Ehistory;
	VectorXi OriginalIndex;
	float InitialNorm;
	int L;
	VectorXf ExactTau;
	float E;
	MatrixXf MatrixTemp;
	vector<unsigned int> IndexOutliers;
	unsigned int NOut;
	bool AllOutliers;
	unsigned int OldSize;
	unsigned int OldSizeOutliers;

	VectorXi MedoidsGlobal;
	VectorXf ExactTauGlobal;
	VectorXi OutliersGlobal;

	// Normalisation
	MatrixXf GramiamMM;
	MatrixXf GramiamMN;
	VectorXi MedoidsGlobalNoOutlier;
	vector<unsigned int> OutliersGlobalNorm;

// Links
	ifstream InputLinks;
	InputLinks.open(Links_file,fstream::in | fstream::binary);
	if ( !InputLinks.is_open() )
	{
		cerr << "Error! Could not open Links file!" << endl;
		return -1;
	}

	// Read number of nodes on 4 bytes
	InputLinks.read((char *)&nb_nodes, 4);
	assert(InputLinks.rdstate() == ios::goodbit);

	// Read cumulative degree sequence: 8 bytes for each node
	// cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
	degrees.resize(nb_nodes);
	InputLinks.read((char *)&degrees[0], nb_nodes*8);

	// Read links: 4 bytes for each link (each link is counted twice)
	nb_links=degrees[nb_nodes-1];
	links.resize(nb_links);
	InputLinks.read((char *)(&links[0]), (long)nb_links*8);

	InputLinks.close();

//Read weights: 4 bytes for each link (each link is counted twice)
	ifstream InputWeights;
	InputWeights.open(Weight_file,fstream::in | fstream::binary);
	if ( !InputWeights.is_open() )
	{
		cerr << "Error! Could not open Weights file!" << endl;
		return -1;
	}

	weights.resize(nb_links);
	InputWeights.read((char *)&weights[0], (long)nb_links*4);
	InputWeights.close();

// Diagonal
	ifstream InputDiagonal;
	InputDiagonal.open(Diag_file,fstream::in | fstream::binary);
	if ( !InputDiagonal.is_open() )
	{
		cerr << "Error! Could not open Diagonal file!" << endl;
		return -1;
	}

	diag.resize(nb_nodes);
	InputDiagonal.read((char *)&diag[0], (long)nb_links*4);
	InputDiagonal.close();

	NFibers = nb_nodes;

	for (unsigned int f=0; f<NumberFascicles; f++)
	{
		// Clearing parameter
		IF=0;
		NOut=0;
		IndexOutliers.resize(1);
		IndexOutliers.clear();
		AllOutliers=false;

		// Index fibers
		ifstream InputFibers;
		InputFibers.open(Fibers_file[f],fstream::in | fstream::binary);
		if ( !InputFibers.is_open() )
		{
			cerr << "Error! Could not open Fiber file " << Fibers_file[f] << endl;
			return -1;
		}

		InputFibers.read((char *)&IF, 4);

		indexFibers.resize(IF);
		InputFibers.read((char *)&indexFibers[0], (long)IF*4);
		InputFibers.close();

		indexFibersNoOutliers.setZero(IF);
		#pragma omp parallel for shared(indexFibersNoOutliers,indexFibers)
		for (unsigned int k=0; k<IF; k++)
			indexFibersNoOutliers[k] = indexFibers[k];

		Gramiam.setZero(IF,IF);

		for (unsigned int k=0; k<IF; k++)
		{
			unsigned int indexa;
			unsigned int indexb;
			if ( indexFibers[k]==0 )
			{
				indexa=0;
				indexb=degrees[ indexFibers[k] ];
			}
			else // it is greater than 0
			{
				indexa=degrees[ indexFibers[k]-1 ];
				indexb=degrees[ indexFibers[k] ];
			}
			int dest;
			float w;
			#pragma omp parallel for shared(Gramiam,links,weights,indexa,indexb,indexFibers) private(dest,w) collapse(2)
			for (unsigned int j=indexa; j<indexb; j++)
			{
				for (unsigned int kk=0; kk<IF; kk++)
				{
					dest = links[j];
					w = weights[j];

					if ( dest==indexFibers[kk] )
					{
						Gramiam(k,kk)=w;
					}
				}
			}
		}

		#pragma omp parallel for shared(Gramiam,diag,indexFibers) collapse(2)
		for (unsigned int i=0; i<NFibers; i++)
		{
			for (unsigned int k=0; k<IF; k++)
			{
				if (i==indexFibers[k])
				{
					Gramiam(k,k)=diag[i];
				}
			}
		}

///////////////////////////////// DETECTION OUTLIERS ////////////////////////

		MatrixTemp=Gramiam;
		float InformationAngle = 0;
		float mean_length = 0;
		float std_length = 0;
		float bound_max = 0;
		float bound_min = 0;
		bool foundOutlier = false;

		// I have to split the two computations into two loops in order to use parallel computing
		VectorXf Diagonal = Gramiam.diagonal();
		#pragma omp parallel for shared(MatrixTemp,Diagonal)
		for (unsigned int i=0; i<IF; i++)
		{
			MatrixTemp.row(i) = MatrixTemp.row(i) / sqrt( Diagonal(i) );
		}

		#pragma omp parallel for shared(MatrixTemp)
		for (unsigned int i=0; i<IF; i++)
		{
			MatrixTemp.col(i) = MatrixTemp.col(i) / sqrt( Diagonal(i) );
		}

		// Constraints on the length
		mean_length = Diagonal.cwiseSqrt().mean();

		VectorXf ToDelete;
		ToDelete.setOnes(IF);
		ToDelete *= mean_length;

		std_length = sqrt( (Diagonal.cwiseSqrt()-ToDelete).array().square().sum() / (IF-1) );

		// Do not parallelize
		vector<unsigned int>::iterator it = IndexOutliers.begin();
		// Assuming it follows a Gaussian distribution
		bound_max = mean_length + 2.6*std_length; // it account for the greatest 1%
		bound_min = mean_length - 2.6*std_length; // it account for the smallest 1%

		for (unsigned int i=0; i<IF; i++)
		{
			InformationAngle = ( MatrixTemp.row(i).array().acos().sum() ) / IF; // Constraint on the angle
			// If the angle is almost 90 (i.e. outlierlimit) or if it is too short or too long
			if( (InformationAngle>outlierlimit) || ( sqrt(Diagonal(i))<bound_min ) || ( sqrt(Diagonal(i))>bound_max ) )
			{
				it = IndexOutliers.insert(it, i);
				foundOutlier = true;
			}
		}

		MatrixTemp.setZero(1,1);
		ToDelete.setZero(1);
		Diagonal.setZero(1);

		if (foundOutlier)
		{
			NOut = IndexOutliers.size();

			IF=IF-NOut;
			if (IF==0)
			{
				cout << "Fascicle wiht only Outliers! Threshold for outlier detection probably too high." << endl;
				AllOutliers=true;
			}

			if (!AllOutliers)
			{
				// Do not parallelize since we need to remove index in decreasing order
				unsigned int index;
				for (unsigned int j=0; j<NOut; j++)
				{
					index = IndexOutliers[j];
					indexFibersNoOutliers = update_vector(indexFibersNoOutliers, index);
					Gramiam = update_GramMM (Gramiam, index);
				}

				if ( (Gramiam.rows() != IF) || 	(Gramiam.cols() != IF) )
					cerr << "Error in updating Gramiam after Outlier detection! " << endl;
			}
		}

///////////////////////////////// FIRST ITERATION ///////////////////////////

		if (!AllOutliers)
		{

			int z = find_medoid (Gramiam);
			OriginalIndex.setZero(IF); // IF does not consider the Outliers
			OriginalIndex.setLinSpaced(IF,0,IF-1);

			Medoids.setZero(1);
			Medoids(0)=z;

			Ehistory.setZero(1);
			InitialNorm = sqrt(Gramiam.sum());
			Ehistory(0)=InitialNorm;

			ExactTau.setZero(1);
			ExactTau = exact_tau (Gramiam, Medoids);

			E = 0;
			E = Cost_function (Gramiam.sum(), ExactTau, Gramiam, Medoids);

			int found=0;
			int iter=0;

			if (E<0)
			{
				E=-sqrt(abs(E));
				found = 1;
			}
			else
				E=sqrt(E);

			Ehistory.conservativeResize(Ehistory.size()+1);
			Ehistory(Ehistory.size()-1)=E;

			GramTemp=Gramiam;

			if (E<=degreePrecision*InitialNorm)
			{
				cout << "OK! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
				cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
				cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
				found=1;
			}

			if ( Medoids.size()>round(IF/3) )
			{
				cout << "Maxinum number of fibers! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
				cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
				cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
				found=1;
			}

			if ( iter ==  (IF-1) )
			{
				cout << "Maxinum number of iterations! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
				cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
				cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
				found=1;
			}

			//////////////////////////////// other iterations //////////////////
			while (found==0)
			{
				iter=iter+1;

				if (remainder(iter,100)==0)
					cout << "Iter " << iter << endl;

				// Update Gram matrix
				VectorXf vectorZ;
				vectorZ = GramTemp.col(z);
				GramTemp = GramTemp - (1/GramTemp(z,z)) * vectorZ * vectorZ.transpose();


				// Delete medoid already chosen; keeping the original reference
				OriginalIndex = update_vector ( OriginalIndex,  z);
				GramTemp = update_GramMM ( GramTemp,  z);

				// Find new Medoid
				z = find_medoid (GramTemp);

				// Update Medoids
				Medoids.conservativeResize(Medoids.size()+1);
				Medoids(Medoids.size()-1)=OriginalIndex(z);

				// Stop criteria
				ExactTau = exact_tau (Gramiam, Medoids);

				int ind;
				double minimum = ExactTau.minCoeff(&ind);
				if (minimum<minValueTau)
				{
					int trovato=0;
					while (trovato==0)
					{
						Medoids=update_vector(Medoids, ind);
						ExactTau = exact_tau (Gramiam, Medoids);

						minimum = ExactTau.minCoeff(&ind);
						if (minimum>=minValueTau)
							trovato=1;

						if (ExactTau.size()<1)
						{
							cerr << "ERROR! minValueTau too small" << endl;
							return -1;
						}
					}
				}

				E = Cost_function (Gramiam.sum(), ExactTau, Gramiam, Medoids);
				if (E<0)
				{
					E=-sqrt(abs(E));
					found = 1;
				}
				else
					E=sqrt(E);

				Ehistory.conservativeResize(Ehistory.size()+1);
				Ehistory(Ehistory.size()-1)=E;

				if (E<=degreePrecision*InitialNorm)
				{
					cout << "OK! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
					cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
					found=1;
				}

				if ( Medoids.size()>round(IF/3) )
				{
					cout << "Maxinum number of fibers! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
					cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
					found=1;
				}

				if ( iter ==  (IF-1) )
				{
					cout << "Maxinum number of iterations! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
					cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
					found=1;
				}

			} // end while

			if ( Medoids.size() != ExactTau.size() )
			{
				std::cerr << "PROBLE! Medoids and Exact tau should be equal!" << std::endl;
				return -1;
			}

	// Saving Medoids

			if (f==0)
			{
				MedoidsGlobal.setZero(Medoids.size());
				ExactTauGlobal.setZero(Medoids.size());
				OutliersGlobal.setZero(IndexOutliers.size());
				#pragma omp parallel for shared(MedoidsGlobal,indexFibersNoOutliers,Medoids)
				for (unsigned int j=0; j<Medoids.size(); j++)
				{
					MedoidsGlobal[j]=indexFibersNoOutliers[Medoids[j]]; // In MedoidsGlobal we have the original index of the fiber
				}

				#pragma omp parallel for shared(OutliersGlobal,indexFibers,IndexOutliers)
				for (unsigned int j=0; j<IndexOutliers.size(); j++)
				{
					OutliersGlobal[j]=indexFibers[IndexOutliers[j]]; // In OutliersGlobal we have the original index of the fiber
				}

				ExactTauGlobal.tail(Medoids.size())=ExactTau;

			}
			else
			{
				OldSize = MedoidsGlobal.size();
				MedoidsGlobal.conservativeResize(OldSize+Medoids.size());
				ExactTauGlobal.conservativeResize(OldSize+Medoids.size());

				OldSizeOutliers=OutliersGlobal.size();
				OutliersGlobal.conservativeResize(OldSizeOutliers+IndexOutliers.size());
				#pragma omp parallel for shared(MedoidsGlobal,indexFibersNoOutliers,Medoids)
				for (unsigned int j=0; j<Medoids.size(); j++)
				{
					MedoidsGlobal[OldSize+j]=indexFibersNoOutliers[Medoids[j]];
				}

				#pragma omp parallel for shared(OutliersGlobal,indexFibers,IndexOutliers)
				for (unsigned int j=0; j<IndexOutliers.size(); j++)
				{
					OutliersGlobal[OldSizeOutliers+j]=indexFibers[IndexOutliers[j]];
				}

				ExactTauGlobal.tail(Medoids.size())=ExactTau;
			}
		}
		else // Fascicles with only Outliers, there is probably a problem
		{
			if (f==0)
			{
				OutliersGlobal.setZero(IndexOutliers.size());
				#pragma omp parallel for shared(OutliersGlobal,indexFibers,IndexOutliers)
				for (unsigned int j=0; j<IndexOutliers.size(); j++)
				{
					OutliersGlobal[j]=indexFibers[IndexOutliers[j]]; // In OutliersGlobal we have the original index of the fiber
				}
			}
			else
			{
				OldSizeOutliers=OutliersGlobal.size();
				OutliersGlobal.conservativeResize(OldSizeOutliers+IndexOutliers.size());

				#pragma omp parallel for shared(OutliersGlobal,indexFibers,IndexOutliers)
				for (unsigned int j=0; j<IndexOutliers.size(); j++)
				{
					OutliersGlobal[OldSizeOutliers+j]=indexFibers[IndexOutliers[j]];
				}
			}
		}

	} // end for Fascicles

	ofstream Outbin;
	Outbin.open("Outliers_global", fstream::out | fstream::binary);
	for (unsigned int j=0 ; j<OutliersGlobal.size() ; j++) {
      		int ou = OutliersGlobal[j];
      		Outbin.write((char *)(&ou),4);
    	}
  	Outbin.close();

////////////////////////////////////// EMPTYNG MEMORY /////////////////////////
Medoids.setZero(1);

IndexOutliers.resize(1);
IndexOutliers.clear();

Gramiam.setZero(1,1);
GramTemp.setZero(1,1);
Ehistory.setZero(1);
ExactTau.setZero(1);
OriginalIndex.setZero(1);

indexFibers.resize(1);
indexFibers.clear();

indexFibersNoOutliers.setZero(1);

///////////////////////////////////////////// NORMALISATION /////////////////////////////////////
// We put all prototypes together and we recompute the weights of the prototypes
// This is useful to recompute the weights of the prototypes close to two
// cluseters (fascicles)

cout << "Normalisation \n" << endl;

unsigned int NMedoids = MedoidsGlobal.size();
NOut = OutliersGlobal.size();

OutliersGlobalNorm.resize(NOut);
#pragma omp parallel for shared(OutliersGlobalNorm,OutliersGlobal)
for (unsigned int j=0; j<NOut; j++)
	OutliersGlobalNorm[j]=OutliersGlobal[j];

OutliersGlobal.setZero(1);

////////////////////////////////////// GRAMIAM COMPUTATION //////////////////

	GramiamMM.setZero(NMedoids,NMedoids);
	GramiamMN.setZero(NMedoids,NFibers);

	#pragma omp parallel for shared(GramiamMM,GramiamMN,MedoidsGlobal,degrees,links,weights)
	for (unsigned int k=0; k<NMedoids; k++)
	{
		unsigned int indexa;
		unsigned int indexb;
		if ( MedoidsGlobal[k]==0 )
		{
			indexa=0;
			indexb=degrees[ MedoidsGlobal[k] ];
		}
		else // it is greater than 0
		{
			indexa=degrees[ MedoidsGlobal[k]-1 ];
			indexb=degrees[ MedoidsGlobal[k] ];
		}

		for (unsigned int j=indexa; j<indexb; j++)
		{
			int dest = links[j];
			float w = weights[j];
			GramiamMN(k,dest)=w;

			for (unsigned int kk=0; kk<NMedoids; kk++)
			{
				if ( dest==MedoidsGlobal[kk] )
				{
					GramiamMM(k,kk)=w;
				}
			}
		}
	}

	#pragma omp parallel for shared(MedoidsGlobal,GramiamMM,GramiamMN,diag)
	for (unsigned int k=0; k<NMedoids; k++)
	{
		GramiamMM(k,k)=diag[MedoidsGlobal[k]];
		GramiamMN(k,MedoidsGlobal[k])=diag[MedoidsGlobal[k]];
	}

///////////////////////////////////// UPDATE PARAMETERS /////////////////////////
	NFibers=NFibers-NOut;
	float InitialNorm2=weights.sum()+diag.sum();

	#pragma omp parallel for shared(InitialNorm2,OutliersGlobalNorm,degrees,weights,links)
	for (unsigned int k=0; k<NOut; k++)
	{
		unsigned int indexa;
		unsigned int indexb;
		if ( OutliersGlobalNorm[k]==0 )
		{
			indexa=0;
			indexb=degrees[ OutliersGlobalNorm[k] ];
		}
		else // it is greater than 0
		{
			indexa=degrees[ OutliersGlobalNorm[k]-1 ];
			indexb=degrees[ OutliersGlobalNorm[k] ];
		}

		for (unsigned int j=indexa; j<indexb; j++)
		{
			int dest = links[j];
			float w = weights[j];
			int found = 0;

			for (unsigned int kk=0; kk<NOut; kk++)
			{
				if ( dest==OutliersGlobalNorm[kk] )
				{
					found=1;
				}
			}

			#pragma omp critical // Otherwise Segmentation fault!!
			{
				if (found == 0)
					InitialNorm2=InitialNorm2-2*w;
				else
					InitialNorm2=InitialNorm2-w;
			}
		}
	}

	for (unsigned int i=0; i<NOut; i++)
	{
		InitialNorm2=InitialNorm2-diag[OutliersGlobalNorm[i]];
	}

///////////////////////////////// Emptying memory //////////////////////////

	links.resize(1);
	links.clear();
	weights.setZero(1);
	diag.setZero(1);
	degrees.resize(1);
	degrees.clear();

////////////////////////////////////// REMOVING OUTLIERS ////////////////////
////// IT TAKES A LOT OF TIME /////// IT SHOULD BE OPTIMISED!!! /////////////

	MedoidsGlobalNoOutlier=MedoidsGlobal;

	if (NOut>0)
	{
		sort(OutliersGlobalNorm.begin(), OutliersGlobalNorm.end(), std::greater<int>());
		MatrixXf TempGamma;
		for (unsigned int j=0; j<NOut; j++)
		{
			if (remainder(j,500)==0)
				cout << "Outlier " << j << " out of "<< NOut << endl;

			unsigned int index = OutliersGlobalNorm[j];
			TempGamma.setZero(GramiamMN.rows(),GramiamMN.cols()-1);
			TempGamma.leftCols(index)=GramiamMN.leftCols(index);
			TempGamma.rightCols(GramiamMN.cols()-1-index)=GramiamMN.rightCols(GramiamMN.cols()-index-1);
			GramiamMN=TempGamma;

			#pragma omp parallel for shared(MedoidsGlobalNoOutlier)
			for (unsigned int k=0; k<NMedoids; k++)
			{
				if (MedoidsGlobalNoOutlier[k]==index)
				{
					cerr << "Problem With Outliers and Prototypes!!!!" << std::endl;
				}
				if (MedoidsGlobalNoOutlier[k]>index)
					MedoidsGlobalNoOutlier[k]=MedoidsGlobalNoOutlier[k]-1;
			}

		}
		TempGamma.setZero(1,1);
	}

	OutliersGlobalNorm.resize(1);
	OutliersGlobalNorm.clear();

	///////////////////////////////////// NORMALISATION /////////////////////////

	VectorXf Ones;
	E=0;
	Ones.setOnes(NFibers);

	if (InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal>0)
	{
		E = sqrt(InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal);
		cout << "Norm exaplined before normalisation: " << 100 - E*100/sqrt(InitialNorm2) << " % of the initial norm (without considering Outliers)" << endl;
		cout << "Norm Prototypes is " << ( sqrt(ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal) / sqrt(InitialNorm2) )*100 << " % of the initial norm" << endl;
	}
	else
	{
		E = sqrt(abs(InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal));
		cout << "Norm exaplined before normalisation: " << 100 + E*100/sqrt(InitialNorm2) << " % of the initial norm (without considering Outliers)" << endl;
		cout << "Norm Prototypes is " << ( sqrt(ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal) / sqrt(InitialNorm2) )*100 << " % of the initial norm" << endl;
	}

 	ExactTauGlobal=GramiamMM.ldlt().solve(GramiamMN*Ones);

	if (InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal>0)
	{
		E = sqrt(InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal);
		cout << "Norm exaplined after normalisation: " << 100 - E*100/sqrt(InitialNorm2) << " % of the initial norm (without considering Outliers)" << endl;
		cout << "Norm Prototypes is " << ( sqrt(ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal) / sqrt(InitialNorm2) )*100 << " % of the initial norm" << endl;
	}
	else
	{
		E = sqrt(abs(InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal));
		cout << "Norm exaplined after normalisation: " << 100 + E*100/sqrt(InitialNorm2) << " % of the initial norm (without considering Outliers)" << endl;
		cout << "Norm Prototypes is " << ( sqrt(ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal) / sqrt(InitialNorm2) )*100 << " % of the initial norm" << endl;
	}

	//////////////////////////////////// REDUNDANCY //////////////////////
	// Here we remove the prototypes that explain less than minValueTau
	// we might have two prototypes from two different fascicles that are close
	// to each other. We may need just one of them.

	int ind;
	double minimum = ExactTauGlobal.minCoeff(&ind);
	if (minimum<minValueTau)
	{
		int trovato=0;
		while (trovato==0)
		{
			GramiamMN=update_GramMN(GramiamMN,MedoidsGlobalNoOutlier,ind);
			MedoidsGlobalNoOutlier=update_vector(MedoidsGlobalNoOutlier, ind);
			MedoidsGlobal=update_vector(MedoidsGlobal, ind);
			GramiamMM=update_GramMM(GramiamMM,ind);

			ExactTauGlobal=GramiamMM.ldlt().solve(GramiamMN*Ones);

			minimum = ExactTauGlobal.minCoeff(&ind);
			if (minimum>=minValueTau)
			{
				trovato=1;
				if (InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal>0)
				{
					E = sqrt(InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal);
					cout << "Norm exaplined after normalisation and redundancy: " << 100 - E*100/sqrt(InitialNorm2) << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal) / sqrt(InitialNorm2) )*100 << " % of the initial norm" << endl;
				}
				else
				{
					E = sqrt(abs(InitialNorm2 - 2 * ExactTauGlobal.transpose()*GramiamMN*Ones + ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal));
					cout << "Norm exaplined after normalisation and redundancy: " << 100 + E*100/sqrt(InitialNorm2) << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal) / sqrt(InitialNorm2) )*100 << " % of the initial norm" << endl;
				}
			}

			if (ExactTauGlobal.size()<1)
			{
				cerr << "ERROR! minValueTau too big" << endl;
				return -1;
			}
		}
	}
	else
		cout << "No Redundancy! \n" << endl;

cout << "Number Prototypes is: " << MedoidsGlobal.size() << " out of " << NFibers+NOut << " fibers, with a reduction of " << MedoidsGlobal.size()*100/(NFibers+NOut) <<  "% \n" << endl;

	if ( MedoidsGlobal.size() != ExactTauGlobal.size() )
	{
		std::cerr << "PROBLEM! Prototypes and Exact tau should be equal!" << std::endl;
		return -1;
	}

cout << "Total number of outliers: " << NOut << " out of " << NFibers+NOut << " fibers, which means " << NOut*100/(NFibers+NOut) << "% \n" << endl;

	ofstream Medbin;
	Medbin.open("Medoids_global_normalised", fstream::out | fstream::binary);
	for (unsigned int j=0 ; j<MedoidsGlobal.size() ; j++) {
      		int me = MedoidsGlobal[j];
      		Medbin.write((char *)(&me),4);
    	}
  	Medbin.close();

	ofstream taubin;
	taubin.open("Tau_global_normalised", fstream::out | fstream::binary);
	for (unsigned int j=0 ; j<ExactTauGlobal.size() ; j++) {
      		float ta = ExactTauGlobal[j];
      		taubin.write((char *)(&ta),4);
    	}
  	taubin.close();

	gettimeofday(&end, NULL);
	delta = double(end.tv_sec  - start.tv_sec) + double(end.tv_usec - start.tv_usec) / 1.e6;
	printf ("It took %f seconds \n",delta);

	return 0;

} // end main
