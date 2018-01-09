#!/usr/bin/python


'''
This code calculates the average frame-to-frame displacement in the y-direction as a function of the position in the box.  First, the displacement of the CM is fixed.  Then, the average y-direction velocity of each particle type and the total is computed in its bulk region. 
cutoff= coarse graining length scale 
L= Size of box along y dimension
Lx= Size of box along x dimension 


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
  

L = string.atoi(sys.argv[1])
Lx = string.atoi(sys.argv[2])
cg = string.atof(sys.argv[3])
filename = sys.argv[4]
xlo = 0
ylo = 0
xhi = Lx
yhi = L
 
skip = 0  #Number of frames to skip at the beginning of the trajectory before beginning analysis
skip1 = 1  #skip1 (must be > skip) : Analyze every skip1'th frame. 
collect = 20000 #number of frames to collect data for - if you want to do all of them, just make it very big! 

#Input file name 

trajectory = "{:}.lammpstrj".format(filename)

print "Input file:{0:s}".format(trajectory)

d = dump(trajectory)
d.sort()
time = d.time()
v1 = np.zeros( int(np.divide( xhi, cg )) )
v2 = np.zeros( int(np.divide( xhi, cg )) )
v1_std = np.zeros( int(np.divide( xhi, cg )) )
v2_std = np.zeros( int(np.divide( xhi, cg )) )
countersum1 = np.zeros( int(np.divide (xhi, cg)) )
countersum2 = np.zeros( int(np.divide( xhi, cg )) )
count = 0 #count is for frames, counter is for particles - sorry, it's confusing
sumx = 0
sumy = 0
grandsumx = 0
grandsumy = 0

for t in time:
    #print "count: ", count
    if count > skip and count % skip1 == 0 and (( count - skip ) / skip1) < collect:

        #elegant module in pizza.py library. Easy way to process the lammpstrj files. d.vecs() goes frame by frame.  
        idlist, typelist, xlist, ylist, zlist = d.vecs( t, "id", "type", "x", "y", "z" )
        cmx = 0
        cmy = 0
        counter = 0

        #Compute center of mass 
        for i in range (len( idlist )):
            if typelist[i] == 1:
                cmx += xlist[i]
                cmy += ylist[i]
                counter = counter + 1

        cmx = cmx / counter
        cmy = cmy / counter


        #compute displacement of each atom.

        dispx = np.zeros( len( idlist ) )
        dispy = np.zeros ( len( idlist ) )

        #Create an empty array to store last x-position of particles *after the COM displacement has been corrected*. This will be used to compute the displacment.
        if count == skip1:
            oldxlist2 = np.zeros( len( idlist ) )
            oldylist2 = np.zeros( len ( idlist ) )

        if count > skip + skip1:
            oldidlist, oldtypelist, oldxlist, oldylist, oldzlist = d.vecs( time[ count - skip1 ], "id", "type", "x", "y", "z" )

        if count > skip + skip1:      
           #The arrays dispx and dispy store the frame to frame displacemnts of each atom  
           #sumx and sumy store the total frame to frame displacements in the x and y directions 
           # grandsumx and grandsumy store the total  displacement across frames. 
           # grandsumx and grandsumy are subtracted from the x and y positions of hte atoms respectively. 
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
          #Fixing displacement.  
          
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
  
        if count > skip + ( 2 * skip1 ):
        #Now, calculate the profile of velocities across the box
            x = xlo
            j = 0   
            while x < xhi: 
                if x + cg > xhi:
                    print "Warning, bin out of bounds"
                sumy1 = 0
                sumy2 = 0
                j = int( np.divide( x, cg ) )
                for i in range( len( idlist ) ):
                    if typelist[i] == 1 and xlist[i] > x and xlist[i] < x + cg:
                        dispy[i] = ylist[i] - oldylist2[i]
                        if dispy[i] > 0.5 * yhi:
                            dispy[i] = -dispy[i] + yhi
                        if dispy[i] < -0.5 * yhi:
                            dispy[i] = -dispy[i] - yhi 

                        countersum1[j] = countersum1[j] + 1
                        sumy1 = sumy1 + dispy[i]
                 
                    if typelist[i] == 2 and xlist[i] > x and xlist[i] < x + cg:
                        dispy[i] = ylist[i] - oldylist2[i]
                        if dispy[i] > 0.5 * yhi:
                            dispy[i] = -dispy[i] + yhi
                        if dispy[i] < -0.5 * yhi:
                            dispy[i] = -dispy[i] - yhi 
                        
                        countersum2[j] = countersum2[j] + 1
                        sumy2 = sumy2 + dispy[i]
                
                v1[j] = v1[j] + sumy1 #bin average velocity of each particle typ ein small slice of the lattive
                v2[j] = v2[j] + sumy2
                v1_std[j] = v1_std[j] + pow( sumy1, 2.0 )
                v2_std[j] = v2_std[j] + pow( sumy2, 2.0 )
                x = x + cg
  
        for i in range( len( idlist ) ):
            oldxlist2[i] = xlist[i]
            oldylist2[i] = ylist[i]
        
    count = count + 1


#Output

outputfile = "vprofY.skip1{:g}.{:}.stats".format( skip1, filename )
foutput = open( outputfile, "w" )

x = xlo + ( 0.5 * cg )
for j in range( len( v1 ) ):
    if countersum1[j] > 0:
        v1[j] = v1[j] / countersum1[j]
        if v1_std[j] / countersum1[j] < pow( v1[j], 2.0 ):
            print "Warning: computed negative standard deviation for type 1 particles"
        v1_std[j] = math.sqrt( ( v1_std[j] / countersum1[j] ) - pow( v1[j], 2.0 ) ) / math.sqrt( countersum1[j] )

    if countersum2[j] > 0:
        v2[j] = v2[j] / countersum2[j]
        if v2_std[j] / countersum2[j] < pow( v2[j], 2.0 ):
            print "Warning: computed negative standard deviation for type 2 particles"
        v2_std[j] = math.sqrt( (v2_std[j] / countersum2[j] ) - pow( v2[j], 2.0 ) ) / math.sqrt( countersum2[j] )

    foutput.write( "%f\t%f\t%f\t%f\t%f\n" % ( x, v1[j], v1_std[j], v2[j], v2_std[j]) )
    x = x + cg

foutput.close()
