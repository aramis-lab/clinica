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

MatrixXf update_GramMM (MatrixXf Gamma, int z)
{
	z=z+1; // I can interpret z as the number of element I want
	
	MatrixXf TempGamma;
	TempGamma.setZero(Gamma.rows()-1,Gamma.cols()-1);

	TempGamma.topLeftCorner(z-1,z-1)=Gamma.topLeftCorner(z-1,z-1);	
	TempGamma.topRightCorner(z-1,Gamma.rows()-z)=Gamma.topRightCorner(z-1,Gamma.rows()-z);	
	TempGamma.bottomLeftCorner(Gamma.rows()-z,z-1)=Gamma.bottomLeftCorner(Gamma.rows()-z,z-1);
	TempGamma.bottomRightCorner(Gamma.rows()-z,Gamma.rows()-z)=Gamma.bottomRightCorner(Gamma.rows()-z,Gamma.rows()-z);

//cout << "GammaMM: \n" << Gamma << "\n" << endl;
//cout << "TempGamma: \n" << TempGamma << "\n" << endl;

	return TempGamma;
}


MatrixXf update_GramMN (MatrixXf GramiamMN, VectorXi Medoids, int z)
{
	int zM=z+1; // This is the index for the rows which I want to eliminate
	int zN=Medoids(z)+1; // This is the index for the columns which I want to eliminate
	
	MatrixXf TempGamma;
	TempGamma.setZero(GramiamMN.rows()-1,GramiamMN.cols()-1);	
	TempGamma.topLeftCorner(zM-1,zN-1)=GramiamMN.topLeftCorner(zM-1,zN-1);
	TempGamma.topRightCorner(zM-1,GramiamMN.cols()-zN)=GramiamMN.topRightCorner(zM-1,GramiamMN.cols()-zN);	
	TempGamma.bottomLeftCorner(GramiamMN.rows()-zM,zN-1)=GramiamMN.bottomLeftCorner(GramiamMN.rows()-zM,zN-1);
	TempGamma.bottomRightCorner(GramiamMN.rows()-zM,GramiamMN.cols()-zN)=GramiamMN.bottomRightCorner(GramiamMN.rows()-zM,GramiamMN.cols()-zN);

//cout << "GramiamMN: \n" << GramiamMN << "\n" << endl;
//cout << "TempGamma: \n" << TempGamma << "\n" << endl;

	return TempGamma;	
}

VectorXi update_vector (VectorXi OriginalIndex, int z)
{
	z=z+1; // I can interpret z as the number of element I want
	VectorXi TempOriginalIndex;
	//cout << "OriginalIndex: \n" << OriginalIndex << "\n" << endl;		
	TempOriginalIndex.setZero(OriginalIndex.size()-1);
	TempOriginalIndex.head(z-1)=OriginalIndex.head(z-1); // I want the first z-1 elements
	//cout << "TempOriginalIndex: \n" << TempOriginalIndex << "\n" << endl;		
	TempOriginalIndex.tail(OriginalIndex.size()-z)=OriginalIndex.tail(OriginalIndex.size()-z); // I want the last N-z elements	
	//cout << "TempOriginalIndex: \n" << TempOriginalIndex << "\n" << endl;
	
	return TempOriginalIndex;
}

MatrixXf update_matrix (MatrixXf Gamma, int z)
{
	z=z+1; // I can interpret z as the number of element I want
	
	MatrixXf TempGamma;
	TempGamma.setZero(Gamma.rows()-1,Gamma.cols()-1);
	//cout << "Gamma: \n" << Gamma << "\n" << endl;
	TempGamma.topLeftCorner(z-1,z-1)=Gamma.topLeftCorner(z-1,z-1);	
	TempGamma.topRightCorner(z-1,Gamma.rows()-z)=Gamma.topRightCorner(z-1,Gamma.rows()-z);	
	TempGamma.bottomLeftCorner(Gamma.rows()-z,z-1)=Gamma.bottomLeftCorner(Gamma.rows()-z,z-1);
	TempGamma.bottomRightCorner(Gamma.rows()-z,Gamma.rows()-z)=Gamma.bottomRightCorner(Gamma.rows()-z,Gamma.rows()-z);
	//cout << "TempGamma: \n" << TempGamma << "\n" << endl;

	return TempGamma;
}


int find_medoid (MatrixXf Gamma)
{
	int IF = Gamma.rows();
	
	//MatrixXf Temp;
	//Temp.setZero(IF,IF);	

	VectorXf Degree2;
	Degree2.setZero(IF);

	#pragma omp parallel for shared(Degree2,Gamma) // add private and public variables // add barrier otherwise
	for (unsigned int j=0; j<IF; j++)	
	{		
		//Temp.row(j) = Gamma.row(j) / sqrt(Gamma(j,j));	
		Degree2(j) = Gamma.row(j).sum();	
		Degree2(j) = Degree2(j) / sqrt(Gamma(j,j));		
	}

	//cout << "Temp matrix:\n" << Temp << "\n" << endl;
	//cout << "Degree:\n" << Degree2 << "\n" << endl;

	Degree2=Degree2.array().square();  
	//cout << "Degree^2:\n" << Degree2 << "\n" << endl;

	int z;
	double s = Degree2.maxCoeff(&z); 

	//cout << "Max index: \n" << z << "\n" << endl;
	return z;
}

MatrixXf reshape (MatrixXf Matrix, VectorXi Indexa, VectorXi Indexb)
{
	int L1 = Indexa.size();
	int L2 = Indexb.size();
	MatrixXf Result;
	Result.setZero(L1,L2);

	#pragma omp parallel for shared(Matrix,Result,Indexa,Indexb) collapse(2)// add private and public variables // add barrier otherwise
	for (unsigned int k=0; k<L1; k++) 
	{
		for (unsigned int kk=0; kk<L2; kk++) 
		{
			Result(k,kk)=Matrix(Indexa(k),Indexb(kk));			
		}
	}

	return Result;
}

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
	//cout << "Original Index: \n" << OriginalIndex << "\n" << endl;
	
	TempMM=reshape(Gamma,Medoids,Medoids);
	//cout << "TempMM:\n" << TempMM << "\n" << endl;

	TempM=reshape(Gamma,Medoids,OriginalIndex);
	//cout << "TempM:\n" << TempM << "\n" << endl;

	VectorXf Ones;
	Ones.setOnes(IF);
	//cout << "Ones:\n" << Ones << "\n" << endl;	

	ExactTau=TempMM.ldlt().solve(TempM*Ones);
	//cout << "ExactTau:\n" << ExactTau << "\n" << endl;	

	return ExactTau;
}

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
		cerr << "Usage: " << argv[0] << " Links Weights Diagonal minValueTau (1) degree_precision (0.1 or 0.15) outlier_limit (0.05) Number_Fascicles Index_fibers1 Index_fibers2 ..." << endl;
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
	
	omp_set_num_threads(nThreads); // in teoria non serve
	
	// Reading Parameters
	char* Links_file = argv[1];
	char* Weight_file = argv[2];
	char* Diag_file = argv[3];	
	double minValueTau = atof(argv[4]);
	double degreePrecision = atof(argv[5]);
	float outlierlimit = atof(argv[6]);
	int NumberFascicles = atoi(argv[7]);

	//cout << "Filename Links: " << Links_file << endl;
	//cout << "Filename Weights: " << Weight_file << endl;
	//cout << "Filename Diagonal Gramiam: " << Diag_file << endl;	
	//cout << "Min value tau: " << minValueTau << endl;
	//cout << "Degree precision: " << degreePrecision << endl;
	//cout << "Outlier limit: " << outlierlimit << endl;
	//cout << "Number Fascicles: " << NumberFascicles << endl;

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
		//cout << "Filename of the fibers: " << Fibers_file[f] << endl;
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

//cout << "Reading Links" << endl;
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
//cout << "Reading Diagonal" << endl;
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
		
		// Emptying parameter
		IF=0;
		NOut=0;
		IndexOutliers.resize(1);
		IndexOutliers.clear();
		AllOutliers=false;

		// Index fibers	
//cout << "Reading Fibers: " << Fibers_file[f] << endl;

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

//cout << "Number links: " << nb_links << endl;
//cout << "Number fibers (nodes): " << nb_nodes << endl;
//cout << "Number fibers module: " << IF << endl;
		
		Gramiam.setZero(IF,IF);		

//cout << "Gramiam matrix initialised :\n" << Gramiam << "\n" << endl;

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
						//Gramiam(kk,k)=w;
					}
				}			
			}
		}

//cout << "Gramiam matrix for fibers without diagonal :\n" << Gramiam << "\n" << endl;

		#pragma omp parallel for shared(Gramiam,diag,indexFibers) collapse(2)// add private and public variables // add barrier otherwise	
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

//cout << "Gramiam matrix for fibers:\n" << Gramiam << "\n" << endl;

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
		//cout << "Diagonal:\n" << Diagonal.transpose() << "\n" << endl;
		//cout << "Diagonal.cwiseSqrt():\n" << Diagonal.cwiseSqrt().transpose() << "\n" << endl;
		mean_length = Diagonal.cwiseSqrt().mean();
		//cout << "mean_length: " << mean_length << "\n" << endl;

		VectorXf ToDelete;
		ToDelete.setOnes(IF);
		ToDelete *= mean_length;

		std_length = sqrt( (Diagonal.cwiseSqrt()-ToDelete).array().square().sum() / (IF-1) );
		//cout << "std_length: " << std_length << "\n" << endl;
	
		// Do not parallelize
		vector<unsigned int>::iterator it = IndexOutliers.begin();
		bound_max = mean_length + 2.4*std_length;
		//cout << "bound_max: " << bound_max << "\n" << endl;
		bound_min = mean_length - 2.4*std_length;
		//cout << "bound_min: " << bound_min << "\n" << endl;

		for (unsigned int i=0; i<IF; i++)	
		{
			//cout << "MatrixTemp.row(i): " << MatrixTemp.row(i) << endl;							
			//cout << "Angle (rad) MatrixTemp.row(i): " << MatrixTemp.row(i).array().acos() << endl;			
			InformationAngle = ( MatrixTemp.row(i).array().acos().sum() ) / IF;

			//cout << "Information Angle: "  << InformationAngle << endl;
			//outlierlimit=1.55;

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
//cout << NOut << " outliers found" << endl;
	
			IF=IF-NOut;
			if (IF==0)
			{
				cout << "Fascicle wiht only Outliers! Threshold for outlier detetion probably too high." << endl;
				AllOutliers=true;
				//ofstream Medbin;
				//Medbin.open("Med.bin", fstream::out | fstream::binary);		    	
			  	//Medbin.close();	  
				//return -1;	
								
			}

			if (!AllOutliers)
			{
				// Do not parallelize since we need to remove index in decreasing order 
				unsigned int index;
	//cout << "Gramiam:\n" << Gramiam << "\n" << endl;
				for (unsigned int j=0; j<NOut; j++)
				{
					index = IndexOutliers[j];	
					indexFibersNoOutliers = update_vector(indexFibersNoOutliers, index);		
					Gramiam = update_matrix (Gramiam, index);
	//cout << "index:\n" << index << "\n" << endl;
	//cout << "Gramiam:\n" << Gramiam << "\n" << endl;
				}		

				if ( (Gramiam.rows() != IF) || 	(Gramiam.cols() != IF) )
					cerr << "Error in updating Gramiam after Outlier detection! " << endl;
			}
		}	

///////////////////////////////// FIRST ITERATION ///////////////////////////	
		if (!AllOutliers)	
		{	

			int z = find_medoid (Gramiam);
			//cout << "z: \n" << z << "\n" << endl;
			//cout << "IF: \n" << IF << "\n" << endl;
		
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

			//cout << "ExactTau: \n" << ExactTau << "\n" << endl;	

			int found=0;	
			int iter=0;

			if (E<0)
			{
				E=-sqrt(abs(E));  	
				found = 1;
			}
			else
				E=sqrt(E);

			//cout << "E: \n" << E << "\n" << endl;

			Ehistory.conservativeResize(Ehistory.size()+1);
			Ehistory(Ehistory.size()-1)=E;	
			//cout << "Ehistory: \n" << Ehistory << "\n" << endl;	

			GramTemp=Gramiam;
			//cout << "GramTemp: \n" << GramTemp << "\n" << endl;

			if (E<=degreePrecision*InitialNorm)
			{	
				cout << "OK! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
				cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
				cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
				//cout << "Tau values:\n" << ExactTau << endl;
				found=1;
			}

			if ( Medoids.size()>round(IF/3) )
			{
				cout << "Maxinum number of fibers! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
				cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
				cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
				//cout << "Tau values:\n" << ExactTau << endl;
				found=1;
			}		
	
			if ( iter ==  (IF-1) )
			{		
				cout << "Maxinum number of iterations! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
				cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
				cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
				//cout << "Tau values:\n" << ExactTau << endl;
				found=1;
			} 

			//////////////////////////////// other iterations //////////////////
			while (found==0)
			{
				iter=iter+1; 
			//cout << "iter: \n" << iter << "\n" << endl;

				if (remainder(iter,100)==0)
					cout << "Iter " << iter << endl;

				// Update Gram matrix  
				VectorXf vectorZ;
				vectorZ = GramTemp.col(z);
				GramTemp = GramTemp - (1/GramTemp(z,z)) * vectorZ * vectorZ.transpose();
				//cout << "GramTemp: \n" << GramTemp << "\n" << endl;			
		
				// Delete medoid already chosen; keeping the original reference  
				//cout << "OriginalIndex: \n" << OriginalIndex << "\n" << endl;			
				OriginalIndex = update_vector ( OriginalIndex,  z);
				//cout << "OriginalIndex: \n" << OriginalIndex << "\n" << endl;	
				GramTemp = update_matrix ( GramTemp,  z);
				//cout << "GramTemp: \n" << GramTemp << "\n" << endl;	

				// Find new Medoid  
				z = find_medoid (GramTemp);
				//cout << "z: \n" << z << "\n" << endl;

				// Update Medoids
				Medoids.conservativeResize(Medoids.size()+1);
				Medoids(Medoids.size()-1)=OriginalIndex(z);	
				//cout << "Medoids: \n" << Medoids << "\n" << endl;	

				// Stop criteria    
				ExactTau = exact_tau (Gramiam, Medoids);
				//cout << "ExactTau: \n" << ExactTau << "\n" << endl;		

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

				//cout << "E: \n" << E << "\n" << endl;
	
				Ehistory.conservativeResize(Ehistory.size()+1);
				Ehistory(Ehistory.size()-1)=E;	
				//cout << "Ehistory: \n" << Ehistory << "\n" << endl;
		
				if (E<=degreePrecision*InitialNorm)
				{	
					cout << "OK! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
					cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
					//cout << "Tau values:\n" << ExactTau << endl;
					found=1;
				}
	
				if ( Medoids.size()>round(IF/3) )
				{
					cout << "Maxinum number of fibers! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
					cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
					//cout << "Tau values:\n" << ExactTau << endl;
					found=1;
				}		
		
				if ( iter ==  (IF-1) )
				{		
					cout << "Maxinum number of iterations! Explained " << 100 - E*100/InitialNorm << " % of the initial norm" << endl;
					cout << "Norm Prototypes is " << ( sqrt(ExactTau.transpose()*reshape(Gramiam,Medoids,Medoids)*ExactTau) / InitialNorm)*100 << " % of the initial norm" << endl;
					cout << "Number of medoids used: " << Medoids.size() << " out of "  << IF+NOut << " fibers, " << "with " << NOut << " outliers \n" << endl;
					//cout << "Tau values:\n" << ExactTau << endl;
					found=1;
				} 		
		
			} // end while

			if ( Medoids.size() != ExactTau.size() )
			{
				std::cerr << "PROBLE! Medoids and Exact tau should be equal!" << std::endl;
				return -1;
			}

			//cout << "Medoids: \n" << Medoids << "\n" << endl;
			//cout << "ExactTau: \n" << ExactTau << "\n" << endl;		
			//cout << "Ehistory: \n" << Ehistory << "\n" << endl;

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

				//cout << "IndexOutliers: \n" << endl;
				#pragma omp parallel for shared(OutliersGlobal,indexFibers,IndexOutliers)  
				for (unsigned int j=0; j<IndexOutliers.size(); j++)
				{
					//cout << IndexOutliers[j] << endl;
					OutliersGlobal[j]=indexFibers[IndexOutliers[j]]; // In OutliersGlobal we have the original index of the fiber					
				}
				//cout << "\n" << endl;			
		
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

				//cout << "IndexOutliers: \n" << endl;
				#pragma omp parallel for shared(OutliersGlobal,indexFibers,IndexOutliers) 
				for (unsigned int j=0; j<IndexOutliers.size(); j++)
				{
					//cout << IndexOutliers[j] << endl;
					OutliersGlobal[OldSizeOutliers+j]=indexFibers[IndexOutliers[j]];				
				}
				//cout << "\n" << endl;
		
				ExactTauGlobal.tail(Medoids.size())=ExactTau;	
			}

	/*
			std::ostringstream ossH;
			ossH << "Ehistory_" << Fibers_file[f] << std::ends;	
			ofstream historybin;
			historybin.open(ossH.str().c_str(), fstream::out | fstream::binary);
			for (unsigned int j=0 ; j<Ehistory.size() ; j++) {
		      		float eh = Ehistory[j];
		      		historybin.write((char *)(&eh),4);
		    	}
		  	historybin.close();

			std::ostringstream ossM;
			ossM << "Med_" << Fibers_file[f] << std::ends;	
			ofstream Medbin;		
			Medbin.open(ossM.str().c_str(), fstream::out | fstream::binary);
			for (unsigned int j=0 ; j<Medoids.size() ; j++) {
		      		int me = Medoids[j];
		      		Medbin.write((char *)(&me),4);
		    	}
		  	Medbin.close();	    

			std::ostringstream ossT;
			ossT << "ExactTau_" << Fibers_file[f] << std::ends;	
			ofstream taubin;
			taubin.open(ossT.str().c_str(), fstream::out | fstream::binary);
			for (unsigned int j=0 ; j<ExactTau.size() ; j++) {
		      		float ta = ExactTau[j];
		      		taubin.write((char *)(&ta),4);
		    	}
		  	taubin.close();

			std::ostringstream ossO;
			ossO << "Outliers_" << Fibers_file[f] << std::ends;	
			ofstream Outliersbin;
			Outliersbin.open(ossO.str().c_str(), fstream::out | fstream::binary);
			for (unsigned int j=0 ; j<IndexOutliers.size() ; j++) {
		      		int outl = IndexOutliers[j];
		      		Outliersbin.write((char *)(&outl),4);
		    	}
		  	Outliersbin.close();	
	*/
		}
		else
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

	//cout << "MedoidsGlobal: \n" << MedoidsGlobal << "\n" << endl;
	//cout << "ExactTauGlobal: \n" << ExactTauGlobal << "\n" << endl;
	//cout << "OutliersGlobal: \n" << OutliersGlobal << "\n" << endl;

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

//cout << "1" << endl;	
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

//cout << "Gram (MxN):\n" << GramiamMN << "\n" << endl;
//cout << "Gram (MxM):\n" << GramiamMM << "\n" << endl;

	#pragma omp parallel for shared(MedoidsGlobal,GramiamMM,GramiamMN,diag) 
	for (unsigned int k=0; k<NMedoids; k++)
	{
		GramiamMM(k,k)=diag[MedoidsGlobal[k]];	
		GramiamMN(k,MedoidsGlobal[k])=diag[MedoidsGlobal[k]];		
	}
			

//cout << "Gram (MxN):\n" << GramiamMN << "\n" << endl;
//cout << "Gram (MxM):\n" << GramiamMM << "\n" << endl;	

///////////////////////////////////// UPDATE PARAMETERS /////////////////////////	
	NFibers=NFibers-NOut;
	float InitialNorm2=weights.sum()+diag.sum();

//cout << "InitialNorm2 with Outliers :\n" << InitialNorm2 << endl;

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

			#pragma omp critical // altrimenti fa Segmentation fault !!
			{
				if (found == 0)
					InitialNorm2=InitialNorm2-2*w;
				else
					InitialNorm2=InitialNorm2-w;	
			}						
		}
	}

//cout << "InitialNorm2 without Outliers 1:\n" << InitialNorm2 << endl;
	
	for (unsigned int i=0; i<NOut; i++)	
	{
		InitialNorm2=InitialNorm2-diag[OutliersGlobalNorm[i]];

	}

//cout << "InitialNorm2 without Outliers :\n" << InitialNorm2 << endl;	

///////////////////////////////// Emptying memory //////////////////////////

	links.resize(1);
	links.clear();
	weights.setZero(1);	
	diag.setZero(1);
	degrees.resize(1);
	degrees.clear();

////////////////////////////////////// REMOVING OUTLIERS ////////////////////
////// IT TAKES A LOT OF TIME /////// IT SHOULD BE OPTIMISED!!! /////////////

//cout << "GramiamMN:\n" << GramiamMN << "\n" << endl;
//cout << "MedoidsGlobal:\n" << MedoidsGlobal << "\n" << endl;
	
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
//cout << "index: " << index  << endl;
//cout << "GramiamMN.col(index): " << GramiamMN.col(index).transpose()  << endl;		
			TempGamma.setZero(GramiamMN.rows(),GramiamMN.cols()-1);
			TempGamma.leftCols(index)=GramiamMN.leftCols(index);
//cout << "TempGamma:\n" << TempGamma << "\n" << endl;
			TempGamma.rightCols(GramiamMN.cols()-1-index)=GramiamMN.rightCols(GramiamMN.cols()-index-1);
//cout << "TempGamma:\n" << TempGamma << "\n" << endl;
			GramiamMN=TempGamma;

			#pragma omp parallel for shared(MedoidsGlobalNoOutlier) 
			for (unsigned int k=0; k<NMedoids; k++)	
			{
				if (MedoidsGlobalNoOutlier[k]==index)
				{
					cerr << "Problem With Outliers and Medoids!!!!!!!!!!!!!" << std::endl;					
				}					
				if (MedoidsGlobalNoOutlier[k]>index)
					MedoidsGlobalNoOutlier[k]=MedoidsGlobalNoOutlier[k]-1;
			}
		
		}
		TempGamma.setZero(1,1);		
	}

//cout << "GramiamMN:\n" << GramiamMN << "\n" << endl;
//cout << "MedoidsGlobalNoOutlier:\n" << MedoidsGlobalNoOutlier << "\n" << endl;

	OutliersGlobalNorm.resize(1);
	OutliersGlobalNorm.clear();

	///////////////////////////////////// NORMALISATION /////////////////////////	

	VectorXf Ones;
	E=0;
	Ones.setOnes(NFibers);

//cout << "ExactTauGlobal.transpose()*GramiamMN*Ones:\n" << ExactTauGlobal.transpose()*GramiamMN*Ones << "\n" << endl;
//cout << "ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal:\n" << ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal << "\n" << endl;

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


//cout << "E before normalisation:\n" << E << "\n" << endl;

 	ExactTauGlobal=GramiamMM.ldlt().solve(GramiamMN*Ones);

//cout << "ExactTauGlobal after normalisation:\n" << ExactTauGlobal << "\n" << endl;
//cout << "InitialNorm2:\n" << InitialNorm2 << "\n" << endl;
//cout << "ExactTauGlobal.transpose()*GramiamMN*Ones:\n" << ExactTauGlobal.transpose()*GramiamMN*Ones << "\n" << endl;
//cout << "ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal:\n" << ExactTauGlobal.transpose()*GramiamMM*ExactTauGlobal << "\n" << endl;

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

cout << "Number Medoids is: " << MedoidsGlobal.size() << " out of " << NFibers+NOut << " fibers, with a reduction of " << MedoidsGlobal.size()*100/(NFibers+NOut) <<  "% \n" << endl;

	if ( MedoidsGlobal.size() != ExactTauGlobal.size() )
	{
		std::cerr << "PROBLE! Medoids and Exact tau should be equal!" << std::endl;
		return -1;
	}

cout << "Total number of outliers: " << NOut << " out of " << NFibers+NOut << " fibers, which means " << NOut*100/(NFibers+NOut) << "% \n" << endl;

//cout << "Medoids: \n" << MedoidsGlobal << "\n" << endl;
//cout << "ExactTau: \n" << ExactTauGlobal << "\n" << endl;

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



