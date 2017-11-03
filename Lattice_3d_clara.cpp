/*
Library of functions for lattice gas simulations and otherwise. C++ implemntations, making them classes. 
*/
#define PI 3.14159
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Lattice_3d.h"
#ifdef CLUSTER
#include "/home/svaikunt/local/lib/include/fftw3.h"
#endif 
#ifndef CLUSTER
#include <fftw3.h>
#endif
#include <assert.h>
#include <math.h>
using namespace std;

Lattice_3d::~Lattice_3d(){
//	gsl_rng_free(r);
} //destructor

Lattice_3d::Lattice_3d(){}
//constructor


void Lattice_3d::initialize(int N1, int Lx1, int Ly1){
	N=N1;
	Lx=Lx1;
    Ly=Ly1;
	initializenocc();
        set_dimension(Lx,Ly,N);
        countmovie=0;
        hstats_cap=0;
        int loopi,loopj,loopk,loopl;
        for (loopi=0;loopi<Ly;loopi++){
          hstatscap.push_back(0.0);
        }
        for (loopi=0;loopi<Lx;loopi++){
          hinstant.push_back(0.0);
        }
        for (loopi=0;loopi<Lx;loopi++){
          vector< vector<int> > twoD;
          for (loopj=0;loopj<Ly;loopj++){
            vector<int> oneD;
            for (loopk=0;loopk<N;loopk++){
              oneD.push_back(0);
            }
            twoD.push_back(oneD);
          }
          nocc.push_back(twoD);
        }
 
}



void Lattice_3d::initializenocc(){
	nocc.clear();
	int loopi,loopj,loopk,loopl;
	for (loopi=0;loopi<Ly;loopi++){
	  hstatscap.push_back(0.0);
	}
	for (loopi=0;loopi<Lx;loopi++){
	  vector< vector<int> > twoD;
	  for (loopj=0;loopj<Ly;loopj++){
	    vector<int> oneD;
	    for (loopk=0;loopk<N;loopk++){
	      oneD.push_back(0);
	    }
	    twoD.push_back(oneD);
	  }
	  nocc.push_back(twoD);
	}
	
}



void Lattice_3d::computehstats_cap(){
  hstats_cap=hstats_cap+1;
  for (int loopi=0;loopi<Lx;loopi++){
    for(int loopj=0;loopj<Ly;loopj++){
      for(int loopk=0;loopk<N;loopk++){
	hstatscap[loopj]+=nocc_image_processed[loopi][loopj][loopk]*pow(Lx,-1.0);
      }
     }
   }
}


void Lattice_3d::computeinstantinterface(){
 for (int loopi=0;loopi<Lx;loopi++){
  hinstant[loopi]=0;
  for (int loopj=0;loopj<Ly;loopj++){
   for(int loopk=0;loopk<N;loopk++){
    hinstant[loopi]+=nocc_image_processed[loopi][loopj][loopk];
    }
   }
 }
}

/*
void Lattice_3d::compute1Dcircularinterface(int n){ //only works for a 2D system with a 1D interface
 //get COM of circular region
 int err=0;
 if (hinstant.size()!=n){
   hinstant.resize(n);
 }
 for (int i=0;i<n;i++){
  hinstant[i]=0;
 }
 int sum=0;
 double COMi=0.0,COMj=0.0,COMk=0.0;
 for (int loopi=0;loopi<Lx;loopi++){
   for (int loopj=0;loopj<Ly;loopj++){
      for(int loopk=0;loopk<N;loopk++){
        if(nocc_image_processed[loopi][loopj][loopk]==1){
          COMi+=loopi;
          COMj+=loopj;
          COMk+=loopk;
          sum++;
         }
        }
       }
     }
 COMi=COMi/sum;
 COMj=COMj/sum;
 COMk=COMk/sum;
 double phi=0.;
 double dphi=2*PI/double(n);
 double phi1=dphi;
 double loopi,loopj;
 int intloopi, intloopj;
 int loopk=0;
 double htemp,htemp1;
 for(int i=0;i<n;i++){
   loopi=COMi;
   loopj=COMj;
   if(i==0){
    do{
       loopi=loopi+cos(phi);
        if (loopi>Lx) loopi=loopi-Lx;
        if (loopi<0) loopi=Lx+loopi;
       intloopi=int(loopi);
       loopj=loopj+sin(phi);
         if (loopj>Ly) loopi=loopj-Ly;
         if (loopj<0) loopj=Ly+loopj;
       intloopj=int(loopj);
       } while(nocc_image_processed[intloopi][intloopj][loopk]==1);
    htemp=sqrt(pow((intloopi-COMi),2.0)+pow((intloopj-COMj),2.0));
   }
   else {
     htemp=htemp1; 
   }
   do{
       loopi=loopi+cos(phi1);
        if (loopi>Lx) loopi=loopi-Lx;
        if (loopi<0) loopi=Lx+loopi;
       intloopi=int(loopi);
       loopj=loopj+sin(phi1);
        if (loopj>Ly) loopj=loopj-Ly;
        if (loopj<0) loopj=Ly+loopj;
       intloopj=int(loopj);
       } while(nocc_image_processed[intloopi][intloopj][loopk]==1);
   htemp1=sqrt(pow((intloopi-COMi),2.0)+pow((intloopj-COMj),2.0));
   hinstant[i]=0.5*(htemp+htemp1);   
   phi=phi1;
   phi1+=dphi;
 }
}
*/


void Lattice_3d::writeoutlattice_XYZ(char *outputfile1){
  int loopi,loopj;
  if (countmovie==0){
    fileout.open(outputfile1);
    fileout<<Lx*Ly<<"\n";
    fileout<<"Iteration\n";
    for (loopi=0;loopi<Lx;loopi++){
      for(loopj=0;loopj<Ly;loopj++){
        if (nocc_image_processed[loopi][loopj][0]==1)
         fileout<<"O"<<"\t"<<loopi<<"\t"<<loopj<<"\t"<<0<<"\n";
        else 
         fileout<<"O"<<"\t"<<0<<"\t"<<0<<"\t"<<-5<<"\n";
      }
    }
   countmovie=countmovie+1;
   }
   else{
     fileout<<Lx*Ly<<"\n";
     fileout<<"Iteration\n";
    for (loopi=0;loopi<Lx;loopi++){
      for(loopj=0;loopj<Ly;loopj++){
        if (nocc_image_processed[loopi][loopj][0]==1)
         fileout<<"O"<<"\t"<<loopi<<"\t"<<loopj<<"\t"<<0<<"\n";
        else
         fileout<<"O"<<"\t"<<0<<"\t"<<0<<"\t"<<-5<<"\n"; 
      }
    }
   }
}



/*void Lattice_3d::neighborlist(int *box1,int loopk1,int loopk2,int loopk3){
  int loopk11=loopk1+(int)(floor(3*gsl_rng_uniform(r)-1));
  int loopk22=loopk2+(int)(floor(3*gsl_rng_uniform(r)-1));
  int loopk33=loopk3+(int)(floor(3*gsl_rng_uniform(r)-1));
  if (loopk11>Lx-1)
    loopk11=0;
  if (loopk11<0)
    loopk11=Lx-1;
  if (loopk22>Ly-1)
    loopk22=0;
  if (loopk22<0)
    loopk22=Ly-1;
  if (loopk33>N-1)
    loopk33=0;
  if (loopk33<0)
    loopk33=N-1;
  box1[0]=loopk11;
  box1[1]=loopk22;
  box1[2]=loopk33;
}
*/
void Lattice_3d::whichbox(int *box1, double rx,double ry,double rz){
  box1[0]=(int)(floor(rx/coarsegrain+0.5));
  box1[1]=(int)(floor(ry/coarsegrain+0.5));
  box1[2]=(int)(floor(rz/coarsegrain+0.5));
  
  
}
