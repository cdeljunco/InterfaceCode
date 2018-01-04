#!/usr/bin/python

'''
This code takes as input a lammps trajectory file, the x and y dimensions of the simulation box, and a lattice cell length l. 
The trajectory must contain a mixture of particles labeled type 1 and type 2 which have x and y positions. It must be 2d and rectangular. 
This code analyzes each frame and 
coarse-grains the simulation box and assigns to each lattice cell of length lxl 
a value of 1 if it contains mostly type 2 particles, and a value of 0 if it contains mostly type 2 particles.
The lattice is printed to an xyz file that can be read by VMD.
'''

import sys
import string
import math
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt
from numpy import linalg as LA
from string import Template
from dump import dump
#from namefiles import namefiles
  

#phi==density of lattice
#cutoff= coarse graining lenght scale 
#L= Size of box along y dimension
#Lx= Size of box along x dimension
#print_from = frame number after which to print. 

phi=string.atof(sys.argv[1]);
cutoff=string.atof(sys.argv[2]);
L=string.atoi(sys.argv[3]);
Lx=string.atoi(sys.argv[4]);
print_from = string.atoi(sys.argv[5])
filename = sys.argv[6]; #filename without extension
xlo=0;
ylo=0;
xhi=Lx;
yhi=L;


#set up coarse grained lattice 
 
xlattice = int( np.divide(xhi,cutoff) ) #number of lattice points in x-direction
ylattice = int( np.divide(yhi,cutoff) ) #number of lattice points in y-direction
latticex = np.zeros( xlattice )  #array of zeros, each representing a lattice point along the x-axis
latticey = np.zeros( ylattice ) # ditto in y 

#Input file name 

trajectory = "{:}.lammpstrj".format(filename)



print "Input file:{0:s}".format(trajectory)


#Commands from the pizza.py library, used to read in and pre-process lammpstrj files.
 
d=dump(trajectory);
d.sort()
time=d.time()

#skip: Number of frames to skip at the beginning of the trajectory  - needs to be 0 in order to properly fix center of mass
#skip1: Number of frames to skip between gathering data
#print_from : Frame number at which to start doing the coarse-grained analysis and printing to file. If we know how long the system takes to equilibrate, we can start printing only after it's equilibrated.
skip = 0;
skip1 = 1;

count = 0
phobulk = 0

#Create output files 

#xyz output file will contain the lattice. You should watch it on vmd to be sure it doesn't have any bubbles, etc.
#name = "newLattice_cutoff{:}_{:}".format(cutoff,filename) 
outputfile = "newLattice_cutoff{:}_{:}.XYZ".format(cutoff, filename)

fxyz = open( outputfile,"w" )

#Fix the center of mass. Look at displacements and subtract the total net displacement. 

for t in time:

    if count > skip and count % skip1 == 0:

        lattice = np.zeros( ( xlattice, ylattice ) ) #array of zeros, each representing a lattice point

        #elegant module in pizza.py library. Easy way to process the lammpstrj files. d.vecs() goes frame by frame.  
        idlist, typelist, xlist, ylist, zlist = d.vecs( t, "id", "type", "x", "y", "z" )      
     #The coarse-grained anaysis begins 
     #Note to Clara: You can use the x y positions of the atoms in xlist and ylist array to do your analysis. The interface should be correcly positioned in this 
     #reference frame.

        if count > print_from:
       
          for i in range(len(idlist)):
              
              #find which lattice cell the particle is in.

              xi = int( np.divide( xlist[i] ,cutoff )) 
              yi = int( np.divide( ylist[i] ,cutoff ))

              #print xi,yi

              if xi < 0:
                  xi = 0
                  print "WARNING: x lattice site out of bounds", xi, xlattice

              if xi >= xlattice:
                  xi = xlattice - 1
                  print "WARNING: x lattice site out of bounds", xi, xlattice

              if yi < 0:
                  yi = 0
                  print "WARNING: y lattice site out of bounds", yi, ylattice

              if yi >= ylattice:
                  yi = ylattice - 1
                  print "WARNING: y lattice site out of bounds", yi, ylattice

              if typelist[i] == 1: #Red particles
                  colorcount = 1
                  lattice[xi][yi] += 1

              # if particle is blue, colorcount = 0, do not increment lattice[xi][yi]
        
          fxyz.write( "%g\n"%( xlattice * ylattice )) #comment lines in xyz file
          fxyz.write(" Iteration\n " )
          latticexcm = 0
          latticeycm = 0
          counter = 0

          #Average density in the middle of the box, i.e. in the blue bulk.
          phobulk = np.mean(lattice[int(xlattice/2):int(xlattice/2)+3, int(ylattice/2):int(ylattice/2)+3])


          for i in range( xlattice ): 

              for j in range( ylattice ):

                  if np.absolute( lattice[i][j] ) >= 0.5 * phobulk: #Then it is red  - a totally red cell has phobulk red particles in it.
                      lattice[i][j] = 1
                      latticexcm += i
                      latticeycm += j
                      counter = counter + 1
                  else: 
                      lattice[i][j] = 0 #Then it is blue or empty
          
          for i in range( xlattice ):
              for j in range( ylattice ):
         
                      if lattice[i][j] == 1:
                          fxyz.write("%g\t%g\t%g\t%g\n"%(1,j,i,0)) #lattice is flipped so that the interface is horizontal.
                      else: 
                          fxyz.write("%g\t %g\t%g\t%g\n"%(1,0,0,-5)) #if the lattice cell is blue, write it out at point 0, 0, -5 - basically it just doesn't show up.
    
    print("Incrementing count")
    count += 1
 
 
    

