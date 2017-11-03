/*
Library of functions for lattice gas simulations and otherwise. C++ implemntations, making them classes. 
*/
//#define CLUSTER
#ifndef UTILGAURD
#define UTILGAURD
#include <cstdlib>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#ifdef CLUSTER
#include "/home/svaikunt/local/lib/include/fftw3.h"
#endif 
#ifndef CLUSTER
#include <fftw3.h>
#endif
using namespace std;
class Utilities{
  public:
  Utilities();
  int ***threeDAllocateintmatrix(int , int , int );
  double ***threeDAllocatedoublematrix(int , int , int );
   double ****fourDAllocatedoublematrix(int , int , int,int );
  int *oneDAllocateintmatrix(int );
  double *oneDAllocatedoublematrix(int );
  int **twoDAllocateintmatrix(int , int );
  double **twoDAllocatedoublematrix(int , int );
  fftw_complex **twoDAllocatecomplexmatrix(int,int);
};
#endif 
