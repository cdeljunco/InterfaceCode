/* 
This is a script that takes an XYZ file containing a trajectory of a lattice with some points filled and some points empty,
and returns the fourier spectrum of height fluctuations at the interface.
The input lattice must have an interface running along the middle of the box in the x-direction, that is, dividing the top and bottom halves of the box.

For each frame of the trajectory:
- The Hoshen-Kopelman algorithm is perfomed on the lattice to obtain a single cluster of filled points with a linear interface. 
    - The resulting movie can be checked to make sure the algorithm was successful and gives a single linear interface.
- The height of the interface in units of lattice spacing is computed.
- The function h(x) is fourier-transformed to get h(k) from which |h(k)^2| (square of the real part) is calculated.
At the end, the values of |h(k)^2| are averaged. 

The output is a file containing k, (|h(k)|^2).

Compile uisng: 
g++ 2Dinterface_capillarysurfacetension_clara.cpp Lattice_3d_clara.cpp Cluster_analysis_clara.cpp Utilities_clara.cpp -lfftw3 -o whatever.out

run using ./whatever.out Lx** Ly** filename

where Lx, Ly are the dimensions of the coarse-grained lattice

**NB** the x and y directions are reversed from the simulation if you used the Script_analysis_latticenew_1.py file to generate the lattice! 
For example, if the sim box was Lx_0 = 100 wide and Ly_0 = 200 tall, and we used a coarse-graining lattice spacing of 2,
then we would use Lx = Ly_0/2 = 100 and Ly = Lx_0/2 = 50 as input parameters to this script.

filename is the name of the XYZ file containing the lattice trajectory *without* the .XYZ extension.
*/

#define PI 3.14159
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>
#include <assert.h>

#include "Lattice_3d.h"

#define ARGS 8


using namespace std;

void fouriertransform();

void avgabsfouriermode(fftw_complex *,double *, double *, int );

void setheightreal(fftw_complex *,double *,int ,double);

fftw_plan p;

int main( int argc,char *argv[]){

    int loopi,loopj,loopk,loopl;
	int checkmerge;

	if (argc==1){
	  cout << "Usage Lx Ly filename(without extension)\n" ;
	  exit(1);
	}
        
    //reading in command line arguments
    int N = 1;
    int Lx = atoi( argv[1] );
    int Ly = atoi( argv[2] );
	char* filename = argv[3];       

    cout << "Setting up lattice\n" ;
    
    int x_in = Lx * Ly; //Number of lattice points       
               
    Lattice_3d mylattice;  //Declares an instance of Lattice_3d. This class has a lot of functions built in 
    mylattice.initialize( N, Lx, int( Ly / 2 ) );      //Intializing. mylattice.nocc is an array into which all the data is written


    //hstats is going to be one dimensional vector which stores the density profile
	vector<double> hstats;
	for (loopk = 0;loopk < int( Ly / 2 ); loopk++ ){
	    hstats.push_back(0.0);//Just fills the vector with zero values -> push_back adds values to the end of the vector
	}
 
    double avgabsheightfourier[ Lx ]; //average, absolute height of the interface in fourier space
    double stdabsheightfourier[ Lx ]; // standard deviation of the absoute height of the interface in fourier space. 
    double hinstant[ Lx ]; //height at 1 moment in time

    for( int i = 0; i < Lx ; i++ ){
        avgabsheightfourier[i] = 0;
        stdabsheightfourier[i] = 0;
        hinstant[i] = 0;
        }

    cout << "Setting up FFTs\n" ;


    //SETTING UP Fast Fourier Transform 

    //fftw_complex is a data type fftw_complex[2] where the 0 element holds the real part and the 1 element holds th eimaginary part.
    // We need to allocate 1 for each column to store the heights in real space.
    fftw_complex* heightreal = ( fftw_complex* )fftw_malloc( sizeof( fftw_complex ) * Lx); 
    //And one for each column to store the heights in fourier space
    fftw_complex* heightfourier = ( fftw_complex* )fftw_malloc( sizeof( fftw_complex ) * Lx);

    /*from the fftw documentation:
    fftw_plan_dft_1d(
    number of points,
    input, 
    where to store output, 
    FFTW_FORWARD is direction of FT, 
    FFTW_MEASURE tells fftw to try and measure run time of several ffts to find the fastest one for array of size N 
    -> FFTW_MEASURE overwrites in & out so make the plan before filling them!! In general, it is best practice th make the plan first.
    */
    p=fftw_plan_dft_1d( Lx, heightreal, heightfourier, FFTW_FORWARD, FFTW_MEASURE); 
       
    //FILE INPUT AND OUTPUT 

    char outfile1name[300], outfile2name[300], infilename[300], hfilename[300], ffilename[300];
    FILE *inputfile, *outputFourier;
    char line[ 300 ];

    //PreHK movie is just a check that the lattice is being correctly read in. It should look the same as the input file.
    sprintf( outfile1name, "PreHKOutputMovie_%s.XYZ", filename);

    //NewOutputMovie is the lattice after the hoshen-kopelman and cluster merge algorithm has acted. 
    //It should contain a single, fluctuating interface with no overhangs or bubbles. 
    //MAKE SURE TO CHECK IT - PROBLEMS CAN HAPPEN HERE  
    sprintf( outfile2name, "NewOutputMovie_%s.XYZ", filename);

    //Input file with lattice
    sprintf( infilename, "%s.XYZ", filename);
    inputfile = fopen( infilename, "r");       
    if (inputfile  == NULL) {
    	perror("Failed to open input file ");
    	return 1;	
	}
    //Output file storing the spectrum of fourier amplitudes. (k, <|h(k)^2|>   
    ofstream Hout;
    sprintf( hfilename, "NewLatticeOutput_%s.stats", filename); //file to write Height stats       
    Hout.open( hfilename );
 
    loopj = 0;
    double junk;
    double latx,laty,latz;
    double prod = 0;

    cout << "Entering loop\n" ;

    while((fgets(line, sizeof(line), inputfile)) !=NULL){//looping through all frames until end of file

        cout << loopj << "\n" ;

        mylattice.initializenocc();//Resets everything to zero at each frame


        // This and above fgets skip first two lines in each frame which are formatted lines for vmd
        fgets(line,sizeof(line),inputfile); 
        loopj++;
        for (int i = 0; i <= x_in - 1; i++){ //x_in is the number of lattice points.

            fscanf( inputfile, "%lf %lf %lf %lf\n" , &junk, &latx, &laty, &latz);

            //cout << latx << "\t" << laty << "\t" << latz << "\n";

            if ( latz >= 0 && laty < Ly / 2.0){ //cell is occupied - we only need to store one half of the interface to compute fourier modes
                mylattice.nocc[ latx ][ laty ][ latz ] = 1; //nocc = 1 if lattice site is occupied, otherwise = 0
                }
            }
        cout << "Read frame " << loopj << "\n" ;
        //cout << outfile1name << "\n" ;
        //mylattice.writeoutnocclattice_XYZ( outfile1name );

        //cout << "Wrote to file\n" ;
	    if (loopj % 1 == 0){ //if (loopj % n == 0) -  do the following computation every n frames. Typically this frame selection happens in the lattice pre-processing (Script_analysis_latticeNew_1.py), in that case use 1.
    
            //Initializing the union find algorithm
	        mylattice.uf_initialize( Lx * int(Ly/2.0) * N, mylattice.nocc );
            int i = 0; 

	        do{
                for (int occ = 0; occ < 2; occ++ ){

		        mylattice.hoshen_kopelman( occ );//calling the hoshen_kopelman algortihm to find clusters of contiguous cells -> once for occupancy 0, once for occupancy 1. 
		        mylattice.check_labelling( occ );

	            }
                //get number of clusters that were merged
	            checkmerge = mylattice.merge_clusters();

	            //cout<<"Number of clusers merged\t" << checkmerge << "\t" << "Iteration" << loopj << "\n";

                //images of the actual lattice nocc used by union find
	            mylattice.nocc_image_preprocessed = mylattice.nocc_image_processed;
	     
                } while( checkmerge != 0 );
            //The system should just have two clusters now. One of red and one of blue.
	     
            //get the column heights
            mylattice.computehstats_cap();

            //Write the merged lattice to an xyz fil - that way we can check whether the algorithm gives a reasonable output
            mylattice.writeoutlattice_XYZ( outfile2name );

            //get the position of the interface
            mylattice.computeinstantinterface();//The position of the instantaenous interface is computed 
 
            for (int loopx = 0;loopx < Lx;loopx++){ 

                hinstant[ loopx ] = mylattice.hinstant[ loopx ];
                cout << hinstant[ loopx ] << "\n";

             }

            setheightreal(heightreal, hinstant, Lx, 0);

            //FT heightreal to heightfourier using plan p  
            fouriertransform(); 

            //add FT height vector to the average
            avgabsfouriermode(heightfourier, avgabsheightfourier, stdabsheightfourier, Lx);

            prod++;

	        }

        } //End of file

    cout << "Calculating averages\n" ;

    hstats=mylattice.hstatscap;

	for(loopj = 0 ; loopj < int( Ly/2.0 ) ; loopj++ ){
	  Hout << loopj << "\t" << hstats[ loopj ] / mylattice.hstats_cap << "\n";
	}


	Hout.flush();
	Hout.close();
        
    for(loopi = 0;loopi < Lx; loopi++) {
        avgabsheightfourier[loopi] = avgabsheightfourier[loopi] * pow( prod * double(Lx) , -1.0);//FT is unnormalized, so need to divide by n (=Lx)
        stdabsheightfourier[loopi] = sqrt(stdabsheightfourier[loopi] * pow( prod * double(Lx) ,-1.0) - pow(avgabsheightfourier[loopi],2));
        }

    //Output written out

    cout << "Writing results to file\n" ;

    int loopk1;
    sprintf(ffilename, "Fourierout_%s.data", filename);
    outputFourier = fopen( ffilename, "w");

    for (loopi = 0; loopi < Lx ; loopi++){
        loopk1 = loopi;//positive frequencies are stored in first half

        if (loopk1 > 0.5 * Lx){
           loopk1 = Lx - loopk1;//negative frequencies are stored in second half, (-k/n same as n-k/n bc periodic) - FFTW documentation
            }

         fprintf( outputFourier, "%f\t%.6f\t%.6f\t%.2f\t%.2f\n",2*M_PI*double( loopk1 ) / double( Lx ), avgabsheightfourier[loopi], stdabsheightfourier[loopi]/ sqrt( prod ), heightfourier[loopi][0] , heightfourier[loopi][1] );       
         
        }

}
	

void fouriertransform(){

  fftw_execute(p);

}


void avgabsfouriermode( fftw_complex *a_1, double *a_2, double *std, int i_1){

  int loopi;
  for (loopi = 0; loopi < i_1; loopi++){
     a_2[loopi] = a_2[ loopi ] + pow( a_1[ loopi ][ 0 ] , 2.0 ) + pow( a_1[ loopi ][ 1 ], 2.0); // add |h(k)|^2 of this frame to running average
     std[loopi] = std[loopi] + pow( pow( a_1[ loopi ][ 0 ], 2.0 ) + pow( a_1[ loopi ][ 1 ], 2.0), 2.0); // add to standard deviation.
  }

}

void setheightreal( fftw_complex *a_1, double *m_1, int i_1, double temp){

    int loopi;
    for (loopi = 0; loopi < i_1; loopi++ ){
        a_1[ loopi ][ 0 ] = m_1[ loopi ] - temp; //set real part to computed height of lattice.
        a_1[ loopi ][ 1 ] = 0; //set imaginary part to zero since it is real
     }

}

