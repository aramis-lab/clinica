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
