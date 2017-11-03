#!/usr/bin/python

#this is a script to compute the number fluctuations in a rectangular probe volume at the interface that is constrained to contained an equal number of red & blue particles

import sys #How to import a library
import string
import math
import numpy as np #nickname a library to type it faster
import scipy.special as sc
from numpy import linalg as LA
from string import Template #why re-import using from?
from dump import dump

#phi==density of lattice
#npvolsx = number of probe volumes (x dimension)
#npvolsy = number of probe volumes (y dimension) 
#L= Size of box along y dimension
#Lx= Size of box along x dimension
#prad =  length of box in y-direction - variable name is leftover from circular probe volume.
sigma=1.0


phi=string.atof(sys.argv[1]);
L=string.atoi(sys.argv[2]);
Lx=string.atoi(sys.argv[3]);
px = string.atof(sys.argv[4]); #x-location of interface to measure at - should be one that has no gap at frames measured.
filename=sys.argv[5]; #prefix of file: single drive, no drive, elliptical drive...
xlo=0;
ylo=0;
xhi=Lx;
yhi=L;
xbox=3.0; #length of box in x-direction
deltaxlo=px-(xbox/2.0)
deltaxhi=px+(xbox/2.0)

#skip: Number of frames to skip as the trajectory is being read in 
#skip>skip1: Number of frames to skipped as the read in trajectory is analyzed - make sure to set so that interface is always the same
#collect = number of frames to collect for
skip=100;
skip1=1;
collect=800; #desired number of statistics - may have to play with this 

#tol: maximum ratio-1 of r/b particles allowed in probe volume at interface. Should depend on probe volume and density.
nave= np.array([4,6,10,16,20])#desired values
nprads = len(nave)
prad = np.zeros(nprads)
for i in range(nprads):
 prad[i]=nave[i]/(xbox*phi) #calculated from bulk density 
 print prad[i]
tol = 0 #exactly equal number of particles
#Input file name

trajectory="{:}.lammpstrj".format(filename)
                       
print "Input file:{0:s}".format(trajectory)

#Commands from the pizza.py librbary. VEry usefull to read lammps trajectories. Will be usefull for the future.
d=dump(trajectory);
d.sort()
time=d.time()

stats=[]
snap=[]
loc=[]
frac=[]

for i in range(nprads):
 stats.append([])
 snap.append([])
 loc.append([])
 frac.append([])
 

#stats = np.zeros((nprads,collect),dtype=np.int)
#snap = np.zeros((nprads, collect),dtype=np.int)
#loc= np.zeros((nprads,collect),dtype=np.int)
#frac = np.zeros((nprads,1),dtype=np.float)


count = 0
it = np.zeros(nprads,dtype=np.int)
frame = np.zeros(nprads,dtype=np.int)

sumx=0
sumy=0



 #look at displacements ..and subtract the total netdisplacement. 
grandsumx=0
grandsumy=0

for t in time:
        if count>skip and count%skip1==0 and np.amin(it)<collect:
           for i in range(nprads):
            if it[i]<collect:
             frame[i]=count
	#elegant module in pizza.py library. Easy way to process the lammpstrj files. d.vecs() goes frame by frame.  
           idlist,typelist,xlist,ylist,zlist=d.vecs(t,"id","type","x","y","z")
           cmx=0
           cmy=0;
           counter=0;
           type=0;
           #Compute center of mass of 1 type of particle
           for i in range (len(idlist)): #Careful! Python uses indentation to mark blocks
            if typelist[i]==1:
             cmx+=xlist[i];
             cmy+=ylist[i];
             counter=counter+1;
             #i.e. we just did this for every type 1 atom
           cmx=cmx/counter;
           cmy=cmy/counter;
           #compute displacement of each atom. Lammps scrambles the output. The i^th atom in the frame need not have the label i. the d.sort() routine takes care of this 
           dispx=np.zeros(len(idlist))
           dispy=np.zeros(len(idlist))
           #print t, cmx,cmy
           if count>1:
            oldidlist,oldtypelist,oldxlist,oldylist,oldzlist=d.vecs(time[count-skip1],"id","type","x","y","z")


           if count>1:      
               #The arrays dispx and dispy store the frame to frame displacemnts of each atom  
               #sumx and sumy store the total frame to frame displacements in the x and y directions 
               # grandsumx and grandsumy store the total  displacement across frames. 
               # grandsumx and grandsumy are subtracted from the x and y positions of hte atoms respectively. 
               sumx=0
               sumy=0
               countersum=0
               for i in range (len(idlist)):
                if typelist[i]!=oldtypelist[i]:
                  print typelist[i],oldtypelist[i]
                dispx[i]=xlist[i]-oldxlist[i]
                dispy[i]=ylist[i]-oldylist[i]
                if dispx[i]>0.5*xhi:
                 dispx[i]=dispx[i]-xhi
                if dispx[i]<-0.5*xhi: 
                 dispx[i]=dispx[i]+xhi
             #   
                if dispy[i]>0.5*yhi:
                 dispy[i]=-dispy[i]+yhi
                if dispy[i]<-0.5*yhi:
                 dispy[i]=-dispy[i]-yhi 
                if typelist[i]==2 or typelist[i]==1: #What does it mean if it is neither? Is this just a check for errors?
                 sumx=dispx[i]+sumx
                 sumy=dispy[i]+sumy
                 countersum=countersum+1
               sumx=sumx/(countersum)
               sumy=sumy/(countersum)
               grandsumx+=sumx
               grandsumy+=sumy
               #Fixing displacement.  
               for i in range (len(idlist)):
                 xlist[i]=xlist[i]-grandsumx
                 while xlist[i]<0 or xlist[i]>=xhi:
                  if xlist[i]<0:
                   xlist[i]+=xhi
                  if xlist[i]>=xhi:
                   xlist[i]=xlist[i]-xhi
                 ylist[i]=ylist[i]-grandsumy
                 if ylist[i]<0:
                  ylist[i]+=yhi
                 if ylist[i]>=yhi:
                  ylist[i]=ylist[i]-yhi

               #Check the center of mass (this is an extra check) 
               counter=0
               cmx=0	
               cmy=0
               for i in range (len(idlist)):
                if typelist[i]==1:
                 cmx+=xlist[i];
                 cmy+=ylist[i];
                 counter=counter+1;
               cmx=cmx/counter;
               cmy=cmy/counter;
               #print t,cmx,cmy 
               #Fixing center of mass.  
               for i in range (len(idlist)):
                 xlist[i]=xlist[i]-cmx+xhi/2.0
                 if xlist[i]<0:
                  xlist[i]+=xhi
                 if xlist[i]>=xhi:
                  xlist[i]=xlist[i]-xhi
                 ylist[i]=ylist[i]-cmy+yhi/2.0
                 if ylist[i]<0:
                  ylist[i]+=yhi
                 if ylist[i]>=yhi:
                  ylist[i]=ylist[i]-yhi
               #Check the center of mass 
               counter=0
               cmx=0
               cmy=0
               for i in range (len(idlist)):
                if typelist[i]==1:
                 cmx+=xlist[i];
                 cmy+=ylist[i];
                 counter=counter+1;
               cmx=cmx/counter;
               cmy=cmy/counter;

               #end of Suri's part, beginning of mine.
#For the range of y values, starting distance at the bottom of the box and going up to one box length from the top.
               for j in range(nprads): 
		       py=0
		       while py<(yhi-2*prad[j]) and it[j]<collect: 
			   n1=0
			   n2=0
			   ntot=0
			   for i in range(len(idlist)):
			       if xlist[i]>deltaxlo and xlist[i]<deltaxhi and ylist[i]>py and ylist[i]<py+prad[j]:
				   if typelist[i]==1:
				       n1=n1+1
				       ntot=ntot+1
				   if typelist[i]==2:
				       n2=n2+1
				       ntot=ntot+1
			   if ntot!=0: 	
			    frac[j].append(float(n1)/float(ntot))
	#Check if the probe volume contains equal numbers of r/b within tol (what should tol be?)
			   if abs(n2-n1)<=tol: #If so, bin statistics, and move up box by 2 box lengths(avoid correlations/overlap of probe volumes)
			       stats[j].append(int(ntot))
			       snap[j].append(frame[j])
			       loc[j].append(py) 
                               it[j]=it[j]+1
			       py=py+2*prad[j]
			   else:
			       py=py+prad[j] #If not, move up and test again. 
  
        count=count+1

for i in range(nprads):

        px=25                            
                                      
	print "stats collected for probe box length {0:g}={1:g}".format(prad[i],it[i])
	print "frame number {0:g}".format(frame[i])
	stats1=np.asarray(stats[i][0:it[i]-1],dtype=np.int)

	#Bin statistics and output.
	subbins = it[i]/10
	temp1=np.zeros(10)

	statsFinal=np.bincount(stats1)
	var=np.zeros(len(statsFinal))
	for j in range(len(statsFinal)):
	    for k in range(0,10):
		temp = np.bincount(stats1[k*subbins:(k+1)*subbins])
		norm = float(np.sum(temp))
		temp.resize(len(statsFinal))
		temp1[k]=float(temp[j]/norm)
	    var[j]=np.nanvar(temp1)

	outputfile="results/numfluc_intconstrain_boxL{0:g}x{1:g}collect{2:g}_{3:s}.stats".format(prad[i],px,collect,filename)
	
	fstats=open(outputfile,"w")
	norm=float(np.sum(statsFinal))

	for j in range(len(statsFinal)):
	    fstats.write("%d\t%d\t%f\t%f\n"%(j,statsFinal[j],float(statsFinal[j])/norm,var[j]))
	fstats.close()

	outputfile="results/numfluc_intconstrainConfigs_boxL{0:g}x{1:g}collect{2:g}_{3:s}.stats".format(prad[i],px,collect,filename)
	#Output configuration statistics here!
	fconfig=open(outputfile,"w")

	for j in range(len(loc[i])):
	  fconfig.write("%d\t%f\n"%(snap[i][j],loc[i][j]))
	fconfig.close()

	outputfile="results/numfluc_intconstrainFrac_boxL{0:g}x{1:g}collect{2:g}_{3:s}.stats".format(prad[i],px,collect,filename)
	fstats=open(outputfile,"w")

	for j in range(len(frac[i])):
	   fstats.write("%f\n"%(frac[i][j]))

	fstats.close()




