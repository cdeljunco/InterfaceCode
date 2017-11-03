/* Compile uisng: g++ 2Dinterface_capillarysurfacetension.cpp Lattice_3d.cpp Cluster_analysis.cpp Utilities.cpp -lfftw3 -o whatever.out */
#define PI 3.14159
//#define CLUSTER
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
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
#include <assert.h>

#include "Lattice_3d.h"

#define ARGS 8

//I am going to write a new code now which computes the average density profile. I will do this for many different area cross sections and then estimate a \gamma_{cap} for the lattice gas. 
//Hopefully, this approach gives an answer that matches \epsilon^2/pi.
using namespace std;
void fouriertransform();
void avgabsfouriermode(fftw_complex *,double *, double *, int ); //average absolute fourier mode?
void setheightreal(fftw_complex *,double *,int ,double);
fftw_plan p;

int main( int argc,char *argv[]){
        int loopi,loopj,loopk,loopl;
	int checkmerge;
	if (argc==1){//Needs Lidx and Lidy of the orginal box to read files. Also needs Ly and Lx=L/coarsegrain 
	  cout<<"Usage Lx Ly Lidy Lidx (for input XYZ file)) Pe\n";
	  exit(1);
	}
        int N=1;
        int Lx=atoi(argv[1]);
        int Ly=atoi(argv[2]);
        double Lidx=atof(argv[3]); //What is Lid?
        double Lid=atof(argv[4]);
	char* filename=argv[5];        
        int x_in=Lx*Ly; //Internal size of the lattice                       
        Lattice_3d mylattice;  //Declares an instance of Lattice_3d. This class has a lot of functions built in 
        mylattice.initialize(N,Lx,int(Ly/2));      //Intializing. mylattice.nocc is an array into whcih all the data is written in -> all what data? 


        //hstats is going to be one dimensional vector which stores the density profile
	vector<double> hstats;
	for (loopk=0;loopk<int(Ly/2);loopk++){
	    hstats.push_back(0.0);//Just fills the vector with zero values -> push_back adds values to the end of the vector
	}
       //initialize and set 3 arrays to 0. 
        double avgabsheightfourier[Lx];
        double stdabsheightfourier[Lx];
        double hinstant[Lx];
        for(int i=0;i<Lx;i++){
          avgabsheightfourier[i]=0;
          stdabsheightfourier[i]=0;
          hinstant[i]=0;
        }
        //SETTING UP FFTW
        fftw_complex*heightreal=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Lx); //fftw_complex is a data type fftw_cpmplex[2] where the 0 element holds the real part and the 1 element holds th eimaginary part;
        fftw_complex*heightfourier=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Lx);
        p=fftw_plan_dft_1d(Lx,heightreal, heightfourier,FFTW_FORWARD,FFTW_MEASURE); //(N, in, out, FFTW_FORWARD is direction of FT, FFTW_MEASURE tells fftw to try and measure run time  of several ffts to find the fastest one for array of size N -> this overwrtes in & out so make plan before filling them!! In general, it is best practic eot make the plan first).
        //The object heightreal gets fourier transformed   

       //FILE INPUT AND OUTPUT 
       char outfilename[300],infilename[300],hfilename[300],ffilename[300];
       FILE *inputfile,*outputFourier;
       char line[300];
       sprintf(outfilename,"NewOutputMovie_%s.XYZ",filename);
       sprintf(infilename,"%s.XYZ",filename); //file with interface xyz coords
       inputfile=fopen(infilename,"r");          
       ofstream Hout;
       sprintf(hfilename,"NewLatticeOutput_%s.stats",filename); //file to write Height stats       
       Hout.open(hfilename);
 
       //setting counter;
       loopj=0;
       char junk;
       double latx,laty,latz;
       double prod=0;
       while((fgets(line, sizeof(line), inputfile)) !=NULL){//looping through all frames until eof
          mylattice.initializenocc();//Resets everything to zero at each frame
          fgets(line,sizeof(line),inputfile); // This and above fgets skip first two lines in each frame which are formatted lines for vmd
          loopj++;
          for (int i=0; i<=x_in-1; i++){
            fscanf(inputfile,"%c %lf %lf %lf\n",&junk,&latx,&laty,&latz);
            if (latz>=0 && laty<Ly/2.0){ //cell is occupied 
             mylattice.nocc[latx][laty][latz]=1;//Data from the processed file is stored in nocc
            }
          }
	  //Read data file.
	  if (loopj%1==0){
	    mylattice.uf_initialize(Lx*int(Ly/2.0)*N,mylattice.nocc);//Initializing the union find algorithm
        int i=0; 
	    do{
              for (int occ=0;occ<2;occ++){
	      cout<<occ<<"\n";
		  mylattice.hoshen_kopelman(occ);//calling hte hoshen_kopelman algortihm to find clusters of contiguous cells -> again, why twice? 
		  mylattice.check_labelling(occ);
	       }
	       checkmerge=mylattice.merge_clusters();//clusters are merged 
	       cout<<"Number of clusers merged\t"<<checkmerge<<"\t"<<"Iteration"<<loopj<<"\n";
	       mylattice.nocc_image_preprocessed=mylattice.nocc_image_processed;//images of the actual lattice nocc used by union find
	     }
	     while(checkmerge!=0);
             //The system should just have two clusters now. One of red and one of blue.
	     mylattice.computehstats_cap();
             mylattice.writeoutlattice_XYZ(outfilename);
             mylattice.computeinstantinterface();//The position of the instantaenous interface is computed 
             for (int loopx=0;loopx<Lx;loopx++){ 
               hinstant[loopx]=mylattice.hinstant[loopx];
             }
             setheightreal(heightreal,hinstant,Lx,0);

             fouriertransform(); //FT heightreal to heightfourier using plan p
             //runnign average computed
             avgabsfouriermode(heightfourier,avgabsheightfourier,stdabsheightfourier,Lx);
             prod++;
             //Write to file 
	     //mylattice.writeoutlattice(outputmov,outputmov1,mylattice.nocc,0);//The last entry here is a FLAG. FLAG=0 writes the lattice
	     //mylatticebulk.writeoutlattice(poutputmov,poutputmov1,mylattice.nocc_image_processed,0);//FLAG=1 writes the processed image. 
	    }
        }//all of the above is while loop that happens for each frame
	hstats=mylattice.hstatscap;
	//Hout<<"#\t"<<epsilon<<"\n";
	for(loopj=0;loopj<int(Ly/2.0);loopj++){
	  Hout<<loopj<<"\t"<<hstats[loopj]/mylattice.hstats_cap<<"\n";
	}
	Hout.flush();
	Hout.close();
        for(loopi=0;loopi<Lx;loopi++) {
          avgabsheightfourier[loopi] = avgabsheightfourier[loopi] * pow(prod*double(Lx),-1.0);//FT is unnormalized, so need to divide by n (=Lx)
          stdabsheightfourier[loopi] = sqrt(stdabsheightfourier[loopi] * pow(prod*double(Lx),-1.0) - pow(avgabsheightfourier[loopi],2));
        }

        //Output written out.
        int loopk1;
        sprintf(ffilename,"Fourierout_%s.data",filename);
        outputFourier=fopen(ffilename,"w");
        for (loopi=0;loopi<Lx;loopi++){
         loopk1=loopi;//positive frequencies are stored in first half
         if (loopk1>0.5*Lx)
           loopk1=Lx-loopk1;//negative frequencies are stored in second half, (-k/n same as n-k/n bc periodic) - FFTW documentation
//Need to get correct FT coefficients here! end up with symmetric array with repeated data.ls fourer
         fprintf(outputFourier,"%f\t%.6f\t%.6f\t%.2f\t%.2f\n",2*M_PI*double(loopk1)/double(Lx),avgabsheightfourier[loopi],stdabsheightfourier[loopi]/sqrt(prod),heightfourier[loopi][0],heightfourier[loopi][1]);       
         }
}
	
	/********************FREE MEMORY*************************///////
	//return 0;
	
void fouriertransform(){
  fftw_execute(p);
}


void avgabsfouriermode(fftw_complex *a_1,double *a_2, double *std, int i_1)
{
  int loopi;
  for (loopi=0;loopi<i_1;loopi++){
     a_2[loopi]=a_2[loopi]+pow(a_1[loopi][0],2.0)+pow(a_1[loopi][1],2.0); // h(k)+|h(k)|^2?  what is the value of h(k)
     std[loopi] = std[loopi] + pow(pow(a_1[loopi][0],2.0)+pow(a_1[loopi][1],2.0),2.0);

  }
}

//mylattice->setheightreal_1(mylattice->heightreal,height,mylattice->x,mylattice->y,asize,avgh);
void setheightreal(fftw_complex *a_1,double *m_1,int i_1,double temp){
        int loopi;
        //printf("\nFrame");
        for (loopi=0;loopi<i_1;loopi++){
                a_1[loopi][0]=m_1[loopi]-temp; //set real part to computed height of lattice.
                a_1[loopi][1]=0; //set imaginary part to zero since it is real
                //printf("\n%f",a_1[loopi][0]);
        }

}

