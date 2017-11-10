/*
Library of functions for lattice gas simulations and otherwise. C++ implemntations, making them classes. 
*/
#define PI 3.14159
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Cluster_analysis.h"
#include <fftw3.h>
#include <assert.h>
using namespace std;

Cluster_analysis::~Cluster_analysis(){

}//destructor

Cluster_analysis::Cluster_analysis(){}
//constructor

void Cluster_analysis::set_dimension(int x,int y,int z){
  
    Lx = x;
    Ly = y;
    N = z; //always just 1 for our 2d analysis
    initialized = 1;

}


int Cluster_analysis::uf_find( int x ){

  int y = x;

  while (labels[y] != y){
    y = labels[y];
  }

  while (labels[x] != x){
    int z = labels[x];
    labels[x] = y;
    x = z;
  }

  return y;
}

int Cluster_analysis::uf_make_set(void){

  labels[0]++;
  assert(labels[0] < n_labels);
  labels.push_back( labels[0] );
  return labels[0];

}

int Cluster_analysis::uf_union(int x,int y){

  return labels[ uf_find(x) ] = uf_find(y);

}

void Cluster_analysis::uf_initialize(int max_labels, const vector< vector < vector<int> > >& nocc_pass){

    int loopi, loopj, loopk;

    if (initialized != 1){
        cout << "Cluster_analysis class not initialized properly !";
        assert( initialized == 1 );
        }

    n_labels = Lx * Ly * N; //1 label per lattice site
    labels.clear();
    labels.push_back(0);

    nocc_image.clear();
    nocc_image_liq.clear();
    nocc_image_vap.clear();
    nocc_image_preprocessed.clear();
    nocc_image_processed.clear();

    for (loopi = 0; loopi < Lx; loopi++){

	    vector< vector<int> > twoD;
	    for (loopj = 0; loopj < Ly; loopj++){
	        vector<int> oneD;
	        for (loopk = 0;loopk < N; loopk++){
	            oneD.push_back( 0 );
	        }
	        twoD.push_back( oneD );
	    }

	  nocc_image.push_back( twoD );
	  nocc_image_liq.push_back( twoD );
	  nocc_image_vap.push_back( twoD );
	  nocc_image_processed.push_back( twoD );
	  //nocc_image_preprocessed.push_back(twoD);
	  nocc_image_preprocessed = nocc_pass;//what is this?
  }	
}

void Cluster_analysis::uf_done(void){

    n_labels = 0 ;
    nocc_image.clear();
    nocc_image_liq.clear();
    nocc_image_vap.clear();
    nocc_image_preprocessed.clear();
    nocc_image_processed.clear();
    labels.clear();

}


int Cluster_analysis::max(int loopi, int loopj, int loopk){

    if (!!loopi == 1) // I think this means: if loopi is not 0.
        return loopi;

    if (!!loopj == 1)
        return loopj;

    if (!!loopk == 1)
        return loopk;

}

int Cluster_analysis::maxbind(int loopi, int loopj, int loopk){

    if (!!loopi * !!loopj == 1 ) //if neither loopi not loopj is 0.
        return uf_union( loopi, loopj);

    if (!!loopi * !!loopk == 1)
        return uf_union( loopi, loopk);

    if (!!loopk * !!loopj == 1)
        return uf_union(loopk, loopj);
}
 
int Cluster_analysis::checksite(int lattice_occ,int occ){

    if (lattice_occ == occ)
        return 1;
    else 
        return 0;

}

//This routine needs to run twice, once for occ=1 and onece for occ=0;  -> why? Not how hoshen-kopelman works...
//uf_initialize(L*L*N); has to be called before this is called

void Cluster_analysis::hoshen_kopelman( int occ ){

    //check
    if ( n_labels != Lx * Ly * N ) {
        cout<<"Initialize union find first\n";
        assert(0);
    }

    labels.clear();
    labels.push_back(0);
    nocc_image.clear();

    for (int loopi = 0; loopi < Lx; loopi++ ){
	    vector< vector<int> > twoD;

	    for (int loopj = 0; loopj < Ly; loopj++){
	        vector<int> oneD;

	        for (int loopk = 0 ; loopk < N ; loopk++){
	            oneD.push_back(0);
	        }
	        twoD.push_back( oneD );
	    }
	    nocc_image.push_back(twoD);
     }


    for (int loopi = 0; loopi < Lx; loopi++){
        for (int loopj = 0; loopj < Ly; loopj++){
            for (int loopk = 0; loopk < N; loopk++){
        
                //same as if (nocc_image_preprocessed[loopi][loopj][loopk] == occ)
	            if (checksite( nocc_image_preprocessed[loopi][loopj][loopk], occ )){

	                int lookx = ( loopi == 0 ? 0 : nocc_image[loopi-1][loopj][loopk]); // lookx = 0 if loopi = 0, otherwise lookx = nocc_mage of lattice site at loopi - 1               
                    int looky = ( loopj == 0 ? 0 : nocc_image[loopi][loopj-1][loopk]);
	                int lookz = ( loopk == 0 ? 0 : nocc_image[loopi][loopj][loopk-1]);


	                switch (!!lookx + !!looky + !!lookz){
	                
                        case 0: //if the lattice sites at x-1, y-1 and z-1 are all 0 OR if we are at an edge
	                        nocc_image[loopi][loopj][loopk] = uf_make_set(); //
	                        break;
	                    case 1: //if one of these is not 0 set nocc_image to its value
	                        nocc_image[loopi][loopj][loopk] = max( lookx, looky, lookz);
	                        break;
	                    case 2: //if 2 of these are not 0 set nocc_image to 
	                        nocc_image[loopi][loopj][loopk] = maxbind( lookx, looky, lookz);
	                        break;
	                    case 3:
	                        nocc_image[loopi][loopj][loopk] = uf_union( lookx, looky);
	                        nocc_image[loopi][loopj][loopk] = uf_union( lookz, nocc_image[loopi][loopj][loopk]);
	                        break;
	                }
	            }
            }
        }
    } 

  //cout<<"Ok then\n";
  //cout.flush();
  //relabelling
  vector <int> new_labels;
  vector <int> count_labels;
  for (int loopi=0;loopi<n_labels;loopi++){
    new_labels.push_back(0);
    count_labels.push_back(0);
  }  
  for (int loopi=0;loopi<Lx;loopi++){
    for (int loopj=0;loopj<Ly;loopj++){
      for (int loopk=0;loopk<N;loopk++){
	if(occ==nocc_image_preprocessed[loopi][loopj][loopk])
	{
	  //cout<<loopi<<loopj<<loopk<<"\n";
	  int x=uf_find(nocc_image[loopi][loopj][loopk]);
	  if (new_labels[x]==0){
	    new_labels[0]++;
	    new_labels[x]=new_labels[0];
	  }
	  nocc_image[loopi][loopj][loopk]=new_labels[x];
	  count_labels[nocc_image[loopi][loopj][loopk]]++;
	}
      }
    }
  }
  
  find_max=0;
  int counter=0;
  int counterlabel=0;
  for (vector <int>::iterator it=count_labels.begin();it !=count_labels.end();++it){
    if (*it>counterlabel){
      find_max=counter;
      counterlabel=*it;
    }
    counter++;
  }
  
  if (occ==1)
  {
    nocc_image_liq=nocc_image;
    find_max_liq=find_max;
    labels_liq=labels;
  }
  if (occ==0){
    nocc_image_vap=nocc_image;
    find_max_vap=find_max;
    labels_vap=labels;
  }
}


int Cluster_analysis::neighborcheck(int occ, int x,int loopi,int loopj,int loopk){
  if (occ==1){
    nocc_image=nocc_image_vap;
    find_max=find_max_vap;
  }
  else{
    nocc_image=nocc_image_liq;
    find_max=find_max_liq;
  }
  int Nor=(loopi==0 ? 0 : nocc_image[loopi-1][loopj][loopk]);
  int Sou=(loopi==Lx-1 ? 0 : nocc_image[loopi+1][loopj][loopk]);
  int Eas=(loopj==0 ? 0 : nocc_image[loopi][loopj-1][loopk]);
  int Wes=(loopj==Ly-1 ? 0 : nocc_image[loopi][loopj+1][loopk]);
  int Up=(loopk==0 ? 0 : nocc_image[loopi][loopj][loopk-1]);
  int Down=(loopk==N-1 ? 0 : nocc_image[loopi][loopj][loopk+1]);
  if (Nor==find_max || Sou==find_max || Eas==find_max || Wes==find_max || Up==find_max || Down==find_max)
    return 0;
  else 
    return 1;
}
int Cluster_analysis::merge_clusters(){
  vector <double> labelflip_liq;
  vector <double> labelflip_vap;
  int checkmerge=0;
  for (int loopi=0;loopi<n_labels;loopi++){
    labelflip_liq.push_back(0);
    labelflip_vap.push_back(0);
  }  
  for (int loopi=0;loopi<Lx;loopi++){
    for (int loopj=0;loopj<Ly;loopj++){
      for (int loopk=0;loopk<N;loopk++){
	if (nocc_image_preprocessed[loopi][loopj][loopk]){
	  labels=labels_liq;
	  int x=nocc_image_liq[loopi][loopj][loopk];
	  if (x!=find_max_liq){
	    if (!labelflip_liq[x]){
	      if (!neighborcheck(1,x,loopi,loopj,loopk)){
		nocc_image_processed[loopi][loopj][loopk]=0;
		labelflip_liq[x]=1;
		checkmerge++;
	      }
	    }
	    else{
	      nocc_image_processed[loopi][loopj][loopk]=0;
	    }
	  }
	  else{
	    nocc_image_processed[loopi][loopj][loopk]=1;
	  }
	}
	else{
	  labels=labels_vap;
	  int x=nocc_image_vap[loopi][loopj][loopk];
	  if (x!=find_max_vap){
	    if (!labelflip_liq[x]){
	      if (!neighborcheck(0,x,loopi,loopj,loopk)){
		nocc_image_processed[loopi][loopj][loopk]=1;
		labelflip_vap[x]=1;
		checkmerge++;
	      }
	    }
	    else{
	      nocc_image_processed[loopi][loopj][loopk]=1;
	    }
	  }
	  else{
	    nocc_image_processed[loopi][loopj][loopk]=0;
	  }
	}
      }
    }
  }
  return checkmerge;
}

//Even this has to be called twice, once for occ=1 and once for occ=0
void Cluster_analysis::check_labelling(int occ){
  if (occ==1)
    nocc_image=nocc_image_liq;
  if (occ==0)
    nocc_image=nocc_image_vap;
  //cout<<occ<<"\t"<<nocc_image[L-1][L-1][N-1]<<"\n";
  int Nor,Eas,Wes,Sou,Up,Down;
  for (int loopi=0;loopi<Lx;loopi++){
    for (int loopj=0;loopj<Ly;loopj++){
      for (int loopk=0;loopk<N;loopk++){
	if (nocc_image_preprocessed[loopi][loopj][loopk]==occ){
	  
	  Nor=(loopi==0 ? 0 : nocc_image[loopi-1][loopj][loopk]);
	  Sou=(loopi==Lx-1 ? 0 : nocc_image[loopi+1][loopj][loopk]);
	  Eas=(loopj==0 ? 0 : nocc_image[loopi][loopj-1][loopk]);
	  Wes=(loopj==Ly-1 ? 0 : nocc_image[loopi][loopj+1][loopk]);
	  Up=(loopk==0 ? 0 : nocc_image[loopi][loopj][loopk-1]);
	  Down=(loopk==N-1 ? 0 : nocc_image[loopi][loopj][loopk+1]);
	  if (Down!=0 && nocc_image[loopi][loopj][loopk]!=Down)
	    cout<<loopi<<loopj<<loopk<<"\t"<<Down<<"\t"<<nocc_image[loopi][loopj][loopk]<<"\t"<<occ<<"\n";
	  assert(Nor==0 || nocc_image[loopi][loopj][loopk]==Nor);
	  assert(Sou==0 || nocc_image[loopi][loopj][loopk]==Sou);
	  assert(Eas==0 || nocc_image[loopi][loopj][loopk]==Eas);
	  assert(Wes==0 || nocc_image[loopi][loopj][loopk]==Wes);
	  assert(Up==0 || nocc_image[loopi][loopj][loopk]==Up);
	  assert(Down==0 || nocc_image[loopi][loopj][loopk]==Down);
	}
      }
    }
  }
}
