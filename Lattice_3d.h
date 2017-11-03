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
#include "Utilities.h"
#include "Cluster_analysis.h"
#endif
#define FLAG_MAX 2
using namespace std;


class Lattice_3d: public Utilities, public Cluster_analysis {
private:
	int N,Lx,Ly,maxstep,size,asize,sollocprobe,nsol,hsol,cts;
	double epsilon,Te,epsilon1;
	int istheresolute;
	int countlattice;
	double v_e_HS;
	double Umbrella;
	double mu;
	double rad;
	double coarsegrain;
	double K_umbrella;
	int delta_rad;
	double deltaint;
	int extent;
	int count[FLAG_MAX];
public:
  int hstats_cap;
  int countmovie;
  ofstream freenergy,Covariance;
  gsl_rng *r;
  long int randomseed;
  int SOSTYPEDEF;
  int GAUSSIAN;
  
  vector< vector < vector<int> > > nocc;
  vector< vector < vector<int> > > nocc1;
  vector< vector < vector<int> > > nocctemp;
  double rsol[20][3];
  vector< vector< vector < vector<double> > > > cv;
  vector< vector< vector < vector<double> > > > cv1;
  vector<double> cr;
  vector<double> hstats;
  vector<double> hstatscap;
  vector<double> hinstant; 
  ofstream fileout,fileout1;
  //fftw_complex *heightreal,*heightfourier,**heightfourier_sol;
  //fftw_plan p;
  Lattice_3d();
  int box[3];
  ~Lattice_3d();
  void initialize(int, int , int);
  void initializenocc();
  void writeoutlattice_XYZ(char *);
  void runmontecarlopropagate_lattice(int ,int,double);
  double montecarlopropagate_lattice(int,double);
  void att_movesolute();
  void att_changerad();
  double energy_modify();
  double Umbrellabias();
  double volumeoccsolute();
  double volumeincell(double *,int, int,int);
  void writesolstats_lattice(int);
  double energyint_lattice(int, int,int);
  double r2(double,double,double);	
  void setnewrad(double);
  int checkbox(double *,int, int,int);
  void cube_sphere_overlap_MC();
  void cube_cube_overlap_MC();
  void whichbox(int*,double,double,double);
  int deltafromrad(double);
  void att_flip_rad();
  double energy_lattice_total();
  double energy_lattice_total(double,double);
  void computehstats();
  void computehstats_cap();
  void computeinstantinterface();
  void compute1Dcircularinterface(int);
  void reinitialize_hstatscap();
  void setnewepsilon(double);
  void neighborlist(int *,int,int,int);
};





