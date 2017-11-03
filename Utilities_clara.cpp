/*
Library of functions for lattice gas simulations and otherwise. C++ implemntations, making them classes. 
*/

#include "Utilities.h"
using namespace std;

Utilities::Utilities(){}


fftw_complex** Utilities::twoDAllocatecomplexmatrix(int n1, int n2){
	fftw_complex **matrix;
	int loopi,loopj;
	matrix=(fftw_complex **)malloc(n1*sizeof(fftw_complex*));
	for (loopi=0;loopi<n1;loopi++){
		matrix[loopi]=(fftw_complex *)malloc(n2*sizeof(fftw_complex));
	}
	return matrix;
}



double*** Utilities::threeDAllocatedoublematrix(int n1, int n2, int n3){
	double ***matrix;
	int loopi,loopj,loopk;
	matrix=(double ***)malloc(n1*sizeof(double**));
	for (loopi=0;loopi<n1;loopi++){
		matrix[loopi]=(double **)malloc(n2*sizeof(double*));
		for(loopj=0;loopj<n2;loopj++){
			matrix[loopi][loopj]=(double *)malloc(n3*sizeof(double));
			for(loopk=0;loopk<n3;loopk++){
				matrix[loopi][loopj][loopk]=0;
			}
		}
	}
	return matrix;
}

double**** Utilities::fourDAllocatedoublematrix(int n1, int n2, int n3,int n4){
	double ****matrix;
	int loopi,loopj,loopk,loopl;
	matrix=(double ****)malloc(n1*sizeof(double***));
	for (loopl=0;loopl<n1;loopl++){
	  matrix[loopl]=(double ***)malloc(n2*sizeof(double**));
	  for (loopi=0;loopi<n2;loopi++){
		matrix[loopl][loopi]=(double **)malloc(n3*sizeof(double*));
		for(loopj=0;loopj<n3;loopj++){
			matrix[loopl][loopi][loopj]=(double *)malloc(n4*sizeof(double));
			for(loopk=0;loopk<n4;loopk++){
				matrix[loopl][loopi][loopj][loopk]=0;
			}
		}
	}
	}
	return matrix;
}

int*** Utilities::threeDAllocateintmatrix(int n1, int n2, int n3){
  int ***matrix;
  int loopi,loopj,loopk;
  matrix=(int ***)malloc(n1*sizeof(int**));
  for (loopi=0;loopi<n1;loopi++){
    matrix[loopi]=(int **)malloc(n2*sizeof(int*));
    for(loopj=0;loopj<n2;loopj++){
      matrix[loopi][loopj]=(int *)malloc(n3*sizeof(int));
      for(loopk=0;loopk<n3;loopk++){
	matrix[loopi][loopj][loopk]=0;
      }
    }
  }
  return matrix;
}
int* Utilities::oneDAllocateintmatrix(int n1){
 int *matrix;
  int loopi,loopj;
  matrix=(int *)malloc(n1*sizeof(int));
  for (loopi=0;loopi<n1;loopi++){
    matrix[loopi]=0;
  }
  return matrix;
}

double* Utilities::oneDAllocatedoublematrix(int n1){
  double *matrix;
  int loopi,loopj;
  matrix=(double *)malloc(n1*sizeof(double));
  for (loopi=0;loopi<n1;loopi++){
    matrix[loopi]=0;
  }
  return matrix;
}


int** Utilities::twoDAllocateintmatrix(int n1, int n2){
  int **matrix;
  int loopi,loopj;
  matrix=(int **)malloc(n1*sizeof(int*));
  for (loopi=0;loopi<n1;loopi++){
    matrix[loopi]=(int *)malloc(n2*sizeof(int));
    for(loopj=0;loopj<n2;loopj++){
      matrix[loopi][loopj]=0;
    }
  }
  return matrix;
}

double** Utilities::twoDAllocatedoublematrix(int n1, int n2){
  double **matrix;
  int loopi,loopj;
  matrix=(double **)malloc(n1*sizeof(double*));
  for (loopi=0;loopi<n1;loopi++){
    matrix[loopi]=(double *)malloc(n2*sizeof(double));
    for(loopj=0;loopj<n2;loopj++){
      matrix[loopi][loopj]=0;
    }
  }
  return matrix;
}


