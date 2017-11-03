#!/usr/bin/python

'''
This code takes as input a lammps trajectory file, the x and y dimensions of the simulation box, and a lattice cell length l. 
The trajectory must contain a mixture of particles labeled type 1 and type 2 which have x and y positions. It must be 2d and rectangular. 
This code analyzes each frame, corrects any center of mass drift, 
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
from namefiles import namefiles
  

#phi==density of lattice
#cutoff= coarse graining lenght scale 
#L= Size of box along y dimension
#Lx= Size of box along x dimension 

phi=string.atof(sys.argv[1]);
cutoff=string.atof(sys.argv[2]);
L=string.atoi(sys.argv[3]);
Lx=string.atoi(sys.argv[4]);
filename=sys.argv[5]; #filename without extension
xlo=0;
ylo=0;
xhi=Lx;
yhi=L;

#skip: Number of frames to skip at the beginning of the trajectory 
#skip1 > skip: Number of frames to skip between gathering data

skip=10;
skip1=1;


#set up coarse grained lattice 
 
xlattice = int( np.divide(xhi,cutoff) ) #number of lattice points in x-direction
ylattice = int( np.divide(yhi,cutoff) ) #number of lattice points in y-direction
latticex = np.zeros( xlattice )  #array of zeros, each representing a lattice point along the x-axis
latticey = np.zeros( ylattice ) # ditto in y 

print lattice[0][0],xlattice,ylattice


#Input file name 

trajectory = "{:}.lammpstrj".format(filename)



print "Input file:{0:s}".format(trajectory)


#Commands from the pizza.py library, used to read in and pre-process lammpstrj files.
 
d=dump(trajectory);
d.sort()
time=d.time()

stats = np.zeros( len(time) )
statsx = np.zeros( len(time) )
statsy = np.zeros( len(time) )
count = 0
sumx = 0
sumy = 0
phobulk = 0

#Create output files 

#xyz output file will contain the lattice. You should watch it on vmd to be sure it doesn't have any bubbles, etc.
name="newLattice_cutoff{:}_{:}".format(cutoff,filename) 
outputfile="{:}.XYZ".format(namefiles.output(name))

fxyz=open(outputfile,"w")

#Fix the center of mass. Look at displacements and subtract the total net displacement. 

grandsumx=0
grandsumy=0

for t in time:

    if count > skip and count % skip1 == 0:

        lattice=np.zeros( ( xlattice, ylattice ) ) #array of zeros, each representing a lattice point

        #elegant module in pizza.py library. Easy way to process the lammpstrj files. d.vecs() goes frame by frame.  
        idlist, typelist, xlist, ylist, zlist = d.vecs( t, "id", "type", "x", "y", "z" )
        cmx = 0
        cmy = 0;
        counter = 0;
        type = 0;

        #Compute center of mass 
        for i in range (len(idlist)):
            if typelist[i] == 1:
                cmx += xlist[i];
                cmy += ylist[i];
                counter = counter + 1;
   
        cmx = cmx / counter;
        cmy = cmy / counter;
   
        #compute displacement of each atom. Lammps scrambles the output. The i^th atom in the frame need not have the label i. the d.sort() routine used above takes care of this 
        dispx = np.zeros( len(idlist) )
        dispy = np.zeros( len(idlist) )
        #print t, cmx,cmy

        if count > 0 :
            oldidlist, oldtypelist, oldxlist, oldylist, oldzlist = d.vecs( time[ count - skip1 ], "id", "type", "x", "y", "z")     
     
            #The arrays dispx and dispy store the frame to frame displacemnts of each atom  
            #sumx and sumy store the total frame to frame displacements in the x and y directions 
            #grandsumx and grandsumy store the total  displacement across frames. 
            #grandsumx and grandsumy are subtracted from the x and y positions of the atoms respectively. 

            sumx = 0
            sumy = 0
            countersum = 0

            for i in range (len(idlist)):
      
                if typelist[i] != oldtypelist[i]: #check that the atoms were properly sorted
                    print " Warning, atoms scrambled ", typelist[i] , " ", oldtypelist[i]
      
                dispx[i] = xlist[i] - oldxlist[i]
                dispy[i] = ylist[i] - oldylist[i]
                
                #check minimum image condition
                if dispx[i] > 0.5 * xhi:  
                    dispx[i] = dispx[i] - xhi
                if dispx[i] < -0.5 * xhi: 
                    dispx[i] = dispx[i] + xhi
   #   
                if dispy[i] > 0.5 * yhi:
                    dispy[i] = -dispy[i] + yhi
                if dispy[i] < -0.5 * yhi:
                    dispy[i] = -dispy[i] - yhi 
      
                if typelist[i] == 2 or typelist[i] == 1: #it should always be 1 or 2 ...
       
                    sumx += dispx[i]
                    sumy += dispy[i]
                    countersum = countersum + 1
                    
            sumx = sumx / countersum
            sumy = sumy / countersum
            grandsumx += sumx
            grandsumy += sumy

            #Fixing displacement.  
            for i in range (len(idlist)):

                xlist[i] = xlist[i] - grandsumx
       
                while xlist[i] < 0 or xlist[i] > = xhi:
                    if xlist[i] < 0:
                        xlist[i] + = xhi
                    if xlist[i] >= xhi:
                        xlist[i] = xlist[i] - xhi
                
                ylist[i] = ylist[i] - grandsumy
                
                while ylist[i] < 0 or ylist[i] > = yhi:
                    if ylist[i] < 0:
                        ylist[i] += yhi
                    if ylist[i] >= yhi:
                        ylist[i] = ylist[i] - yhi
     
            #Check the center of mass (this is an extra check) 
            counter = 0
            cmx = 0	
            cmy = 0

            for i in range (len(idlist)): #only count type 1 particles - this keeps the slab in the middle of the box.
                
                if typelist[i] == 1:
                    cmx += xlist[i];
                    cmy += ylist[i];
                    counter = counter+1;
                
            cmx = cmx / counter;
            cmy = cmy / counter;
            #print t,cmx,cmy 
     
            #Fixing center of mass.  

            for i in range (len(idlist)):

                xlist[i] = xlist[i] - cmx + xhi / 2.0
       
                if xlist[i] < 0:
                    xlist[i] += xhi
                if xlist[i] >= xhi:
                    xlist[i] = xlist[i] - xhi
            
                ylist[i] = ylist[i] - cmy + yhi / 2.0
                if ylist[i] < 0:
                    ylist[i] += yhi
                if ylist[i] >= yhi:
                    ylist[i] = ylist[i] - yhi

            #Check the center of mass 
            counter = 0
            cmx = 0
            cmy = 0
            for i in range (len(idlist)):

                if typelist[i] == 1:
                    cmx += xlist[i];
                    cmy += ylist[i];
                    counter = counter + 1;
     
            cmx = cmx / counter;
            cmy = cmy / counter;

            print("cmx : " , cmx, " cmy : ", cmy)

            #Now the slab shold reallly be in the middle.
 
      
     #The coarse-grained anaysis begins 
     #Note to Clara: You can use the x y positions of the atoms in xlist and ylist array to do your analysis. The interface should be correcly positioned in this 
     #reference frame.


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

        xleft = xi - 1;
        xright = xi + 1;
        ydown = yi - 1;
        yup = yi + 1;
       
        if xleft<0:
            xleft = xlattice - 1

        if ydown < 0:
            ydown = ylattice - 1

        if xright > xlattice-1 :
            xright = 0;

        if yup > ylattice - 1:
            yup = 0;
    
        colorcount = 0

        if typelist[i] == 2: #Blue particles
            colorcount = 1
            lattice[xi][yi] += colorcount

            type=type+1 #moving up in the typelist?

        else: 
            colorcount = 0 #red particles
        
        lattice[xi][yi] += colorcount

        if count > skip1:

            fxyz.write( "%g\n"%( xlattice * ylattice )) #comment lines in xyz file
            fxyz.write(" Iteration\n " )
            latticexcm = 0
            latticeycm = 0
            counter = 0

            #Average density way over to the left, i.e. in the blue bulk.
            phobulk = np.mean(lattice[0:3, 0:3]) / 9.0

            print phobulk
        
            phobulk1 = cutoff * cutoff * phi #number of particles in one lattice cell

            for i in range( xlattice ): 

                for j in range( ylattice ):

                    if np.absolute( lattice[i][j] ) >= 0.5 * phobulk: #Then it is blue - a totally blue cell has phobulk blue particles in it.
                        lattice[i][j] = 1
                        latticexcm += i
                        latticeycm += j
                        counter = counter + 1
                    else: 
                        lattice[i][j] = 0 #Then it is red or empty
     
            for i in range( xlattice ):
                for j in range( ylattice ):
   
                        if lattice[i][j] == 1:
                                fxyz.write("%g\t%g\t%g\t%g\n"%(1,j,i,0)) #lattice is flipped so that the interface is horizontal.
                        else:
                                fxyz.write("%g\t %g\t%g\t%g\n"%(1,0,0,-5)) #if the lattice cell is red, write it out at point 0, 0, -5 - basically it just doesn't show up.
    count += 1
 
 
    

