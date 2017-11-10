/*
Library of functions for lattice gas simulations and otherwise. C++ implemntations, making them classes. 
*/
#ifndef SOSGAURD
#define SOSGAURD
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#endif

using namespace std;


class Cluster_analysis{
private:
  int Lx,Ly,N,initialized;
public:
  int n_labels;
  int find_max;
  int find_max_liq;
  int find_max_vap;
  vector <int> labels;
  vector <int> labels_liq;
  vector <int> labels_vap;
  vector< vector < vector<int> > > nocc_image_liq;
  vector< vector < vector<int> > > nocc_image_vap;
  vector< vector < vector<int> > > nocc_image_processed;
  vector< vector < vector<int> > > nocc_image;
  vector< vector < vector<int> > > nocc_image_preprocessed;
 Cluster_analysis();
 ~Cluster_analysis();
  void set_dimension(int,int,int);
  int get_dimension(void);
  int uf_find(int);
  int uf_make_set();
  int uf_union(int,int);
  void uf_initialize(int,const vector< vector < vector<int> > >& );
  void uf_done(void);
  int max(int,int,int);
  int maxbind(int, int,int);
  int checksite(int,int);
  void hoshen_kopelman(int);
  int neighborcheck(int,int,int,int,int);
  int merge_clusters();
  void check_labelling(int);
};





