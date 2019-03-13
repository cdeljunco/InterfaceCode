/*
This is a script to compute <F^2>, <F>^2, <F^2>_0, and the work done per cycle given a restart file that contains a snapshot in the steady state.

run using mpirun -np 1 Lammps_drivensystem_F.out 1 input.in Pe tau rho dt

The output is 2 files: 

dateTimeString_Drivenv2.work.density''.Pe''.Tau''.fyoverfx''.dt''.stats

and
  
dateTimeString_Drivenv2.initCorrs.density''.Pe''.Tau''.fyoverfx''.dt''.stats

where dateTimeString is a string containing the date and time in YYMMDD-HHMMSS format and '' indicates where the relevant input variable shows up.

 - the 'work' output file contains 4 columns:

time | cycle number (units of tau) | running average of work per cycle value | work for this cycle

The work rate is the average of the last column divided by tau.

 - the 'initCorrs' file contains 13 columns:

time | running average <f(0)^2> | <f(0)^2> for this loop  | ignore | r.a. <f(0)f(dt)> | <f(0)f(dt)> this loop | ignore | ignore | ignore | ignore | r.a. <eta(0)f(dt)> | <eta(0)f(dt) this loop | ignore 
  
*/
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// c++_driver = simple example of how an umbrella program
//              can invoke LAMMPS as a library on some subset of procs
// Syntax: c++_driver P in.lammps
//         P = # of procs to run LAMMPS on
//             must be <= # of procs the driver code itself runs on
//         in.lammps = LAMMPS input script
// See README for compilation instructions
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <ctime>
#include "math.h"
#include "string.h"
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "update.h"
#include "library.h"
#include "modify.h"
#include "fix.h"
#include "assert.h"
using namespace LAMMPS_NS;

int main(int narg, char **arg)
{
  // setup MPI and various communicators
  // driver runs on all procs in MPI_COMM_WORLD
  // comm_lammps only has 1st P procs (could be all or any subset)

  MPI_Init(&narg,&arg);
  char outputfile[100];
  char inputline[1000];
  if (narg ==1) {
    printf("Syntax: c++_driver P in.lammps Pe Tau rho\n");
    exit(1);
  }
  long int randomseed=100025525;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,randomseed);
  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  int nprocs_lammps = atoi(arg[1]);
  //input file name is arg[2]
  double Pe=atof(arg[3]);
  double tau=atof(arg[4]);
  double rho=atof(arg[5]);
  double dt=atof(arg[6]);
  if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int tmax = 100*tau;
  
  int lammps;
  if (me < nprocs_lammps) lammps = 1;
  else lammps = MPI_UNDEFINED;
  MPI_Comm comm_lammps;
  MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
  
  // open LAMMPS input script

  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[2],"r");
    if (fp == NULL) {
      printf("ERROR: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  // run the input script thru LAMMPS one line at a time until end-of-file
  // driver proc 0 reads a line, Bcasts it to all procs
  // (could just send it to proc 0 of comm_lammps and let it Bcast)
  // all LAMMPS procs call input->one() on the line
  char *arg1 = "-screen"; //suppress output to screen
  char *arg2 = "none"; 
  //char *arg3 = "-log";
  //char *arg4 = "none";
  //char *arg5 = "-echo";
  //char *arg6 = "none";
  char *commargs[2]={arg1,arg2};//,arg3,arg4,arg5,arg6};
  LAMMPS *lmp;
  if (lammps == 1) lmp = new LAMMPS(2,commargs,comm_lammps);
  //if (lammps == 1) lmp = new LAMMPS(0,NULL,comm_lammps);
  // SEE EARLIER PART OF CODE FOR CONDITIONS WHEN lammps!=1. This happens only when input is faulty.
  

  int n,loopi,loopj;
  char line[1024];
  while (1) {
    if (me == 0) {
	//if (me==0) is true when the processor has id 0. 
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
      if (n == 0) fclose(fp);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
    if (lammps == 1) lmp->input->one(line);
  }
    
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  FILE *fileout,*fileouttest, *fileoutD, *forcesfile;
  sprintf(outputfile,"DrivenLJ.work.density%.2f.Pe%.2f.Tau%.2f.2f.dt%.4f.stats",rho,Pe,tau,dt);
  fileout=fopen(outputfile,"w");


  double *PotEng;
  double *Terun;
  double *oldPotEng;
  double *KinEng;
  double *oldKinEng;
  double *ptr;
  int natoms = static_cast<int> (lmp->atom->natoms);
  int nlocal= static_cast<int> (lmp->atom->nlocal);

  int counter =0;
  double dw=0;
  double dwavg=0;
  double dwf=0.;
  double dwfavg=0.;
  double D=0.;
  double Dave=0.;
  double gamma=100;
  double gamma2 = 100*100;
  int told = 0;
  int skip = 1;
  double fx = 0;
  double fy = 0.;
  double phase = 0.;
  double fxnc = 0.;
  double fync = 0.;
  double ex = 0.;
  double ey = 0.;
  
  lmp->input->one("run 0"); //make sure everything is updated

  double **x=lmp->atom->x;
  double **f=lmp->atom->f;
  int *type =lmp->atom->type;
  bigint *ntimestep=&lmp->update->ntimestep;
  std::cout << *ntimestep << std::endl;
  int time = *ntimestep;
  //std:: cout << time << std::endl;


//Loop to calculate work

tmax = tmax/dt; 

  for (int loopj=0;loopj<tmax;loopj++){

    if (loopj%int(tau/dt)==0){
        std::cout << "loop " << loopj << " of " << tmax << std::endl;
    }

    time = *ntimestep;

    phase = 2*M_PI*time*dt/tau;
    ex = sin(phase);
    ey = cos(phase);
    fxnc = Pe*ex;
    fync = Pe*ey;

for (int loopi=0;loopi<nlocal;loopi++){

    int typei = type[loopi];

    fx =  f[loopi][0]-(2-typei)*fxnc; //2-type = 2-1=1 for driven particles , 2-type = 2-2=0 for undriven.
    fy =  f[loopi][1]-(2-typei)*fync;

//compute work new way (dW = int F(x)dx = F_conservative(t) * dx_driven(t)

    if (typei==1){
     dwf -= (fx) * fxnc * dt / gamma +  (fy) * fync * dt / gamma ;
    }


}


//bin/print work values 
   if (loopj%((int)(tau/dt))==0 && loopj>0){

    dwfavg += dwf;
    counter += 1;
    fprintf( fileout, "%f\t%d\t%f\t%f\n", loopj * dt, counter, dwfavg / counter, dwf);
    fflush( fileout ) ;
    dwf = 0;
     
    }
 
//integrate at the end of the loop
   lmp->input->one("run 1");

}

if (lammps == 1) delete lmp;

 // close down MPI
  MPI_Finalize();
}

