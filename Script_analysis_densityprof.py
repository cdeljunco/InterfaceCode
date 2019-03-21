#!/usr/bin/python

'''
This is a script to compute the number density of active and passive particles as a function of x.

It takes as input the name of a lammpstrj file, the dimensions of the box and returns a file with the columns: x, <rho(x)>, and <(rho(x))^2> 

'''


import sys #How to import a library
import string
import math
import numpy as np #nickname a library to type it faster
import scipy.special as sc
from string import Template #why re-import using from?
from dump import dump
from scipy.optimize import curve_fit 
  

sigma = 1.0

L = string.atoi(sys.argv[1]); #convert strign to int
Lx = string.atoi(sys.argv[2]);   
filename = sys.argv[3]; #file without extension

xlo = 0;
ylo = 0;
xhi = Lx;
yhi = L;


#Input file name

trajectory = "{:}.lammpstrj".format(filename)

print "Input file:{0:s}".format(trajectory)

#read in trajectory
d = dump(trajectory);
d.sort()
time = d.time()
dx = 0.1;
rho  = np.zeros(int(Lx / dx))
rho_2 = np.zeros(int(Lx / dx))
rho_total  = np.zeros(int(Lx / dx))
rho_total_2 = np.zeros(int(Lx / dx))
count = 0
sumx = 0
sumy = 0


#look at displacements ..and subtract the total netdisplacement. 
grandsumx = 0
grandsumy = 0


for t in time: 
           
    idlist, typelist, xlist, ylist, zlist = d.vecs(t,"id","type","x","y","z")
    
    # Part 1 : fix center of mass
    cmx = 0
    cmy = 0;
    counter = 0;
    type = 0;
    
    #1.a Compute center of mass of 1 type of particle
    for i in range ( len( idlist ) ): #Careful! Python uses indentation to mark blocks
        if typelist[i] == 1:
            cmx += xlist[i];
            cmy += ylist[i];
            counter = counter + 1;
    
    cmx = cmx / counter;
    cmy = cmy / counter;
    
    #1.b compute displacement of each atom.  
    dispx = np.zeros( len( idlist ) )
    dispy = np.zeros( len( idlist ) )
    
    if count > 1:
        oldidlist, oldtypelist, oldxlist, oldylist, oldzlist = d.vecs( time[ count - 1], "id", "type", "x", "y", "z")

        sumx = 0
        sumy = 0
        countersum = 0
        for i in range ( len( idlist ) ):
            if typelist[i] != oldtypelist[i]:
                print typelist[i], oldtypelist[i]

            dispx[i] = xlist[i] - oldxlist[i]
            dispy[i] = ylist[i] - oldylist[i]

            if dispx[i] > 0.5 * xhi:
                dispx[i] = dispx[i] - xhi
            if dispx[i] < -0.5 * xhi:
                dispx[i] = dispx[i] + xhi

            if dispy[i] > 0.5 * yhi:
                dispy[i] = -dispy[i] + yhi
            if dispy[i] < -0.5 * yhi:
                dispy[i] = -dispy[i] - yhi

            sumx = dispx[i] + sumx
            sumy = dispy[i] + sumy
            countersum = countersum + 1
        
        sumx = sumx / countersum
        sumy = sumy / countersum
        grandsumx += sumx
        grandsumy += sumy
        

        #Adjust atomd based on displacement.

        for i in range ( len( idlist ) ):
            xlist[i] = xlist[i] - grandsumx
            while xlist[i] < 0 or xlist[i] >= xhi:
                if xlist[i] < 0:
                    xlist[i] += xhi
                if xlist[i] >= xhi:
                    xlist[i] = xlist[i] - xhi

            while ylist[i] < 0 or ylist[i] >= yhi:
                ylist[i] = ylist[i] - grandsumy
                if ylist[i] < 0:
                    ylist[i] += yhi
                if ylist[i] >= yhi:
                    ylist[i] = ylist[i] - yhi

        #Check the center of mass (this is an extra check)
        counter = 0
        cmx = 0
        cmy = 0
        for i in range (len(idlist)):
            if typelist[i] == 2:
                cmx += xlist[i]
                cmy += ylist[i]
                counter = counter + 1
        cmx = cmx / counter
        cmy = cmy / counter

        #Fixing center of mass.
        for i in range (len( idlist ) ):
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


    rhotemp, binedges = np.histogram( xlist, bins = int( Lx/  dx ), weights = [2 - x for x in typelist] ) # gives a weight of 1 to particle type 1, weight of 0 to particle type 2 - counts type 1 particles. 
    rho = np.add( rhotemp/(dx * L), rho) 
    rho_2 = np.add( np.square( rhotemp / (dx * L) ) , rho_2 ) #standard dev
    
    rhotemp, binedges = np.histogram( xlist, bins = int( Lx/  dx )) # gives a weight of 1 to all particles 
    rho_total = np.add( rhotemp/(dx * L), rho_total) 
    rho_total_2 = np.add( np.square( rhotemp / (dx * L) ) , rho_total_2 ) #standard dev
      
    count = count + 1

factor = count 
             
rho = rho / factor # get number density by dividing by area = dx * L, take average over frames by dividing by factor. 
rho_2 = np.sqrt( rho_2 / factor - np.square(rho) ) / math.sqrt(factor)
rho_total = rho_total / factor # get number density by dividing by area = dx * L, take average over frames by dividing by factor. 
rho_total_2 = np.sqrt( rho_total_2 / factor - np.square(rho_total) ) / math.sqrt(factor)


x = xlo + ( dx / 2.0)

outputfile = "densityprof.wtotal.dx{:.1f}.{:}.stats".format(dx, filename)

fstats = open( outputfile, "w" )  

for i in range(len(rho)):
    fstats.write( "%f\t%f\t%f\t%f\t%f\n"%( x, rho[i], rho_2[i] , rho_total[i], rho_total_2[i] ))
    x += dx

fstats.close()

