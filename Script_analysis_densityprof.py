#!/usr/bin/python

import sys #How to import a library
import string
import math
import numpy as np #nickname a library to type it faster
import scipy.special as sc
from numpy import linalg as LA
from string import Template #why re-import using from?
from dump import dump
from namefiles import namefiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
  

def func(x,a,b,c,d):
    return b*sc.erf((x-d)/a)+c

def residuals(p,y,x):
    a,b,c,d = p
    err = y - b*sc.erf((x-d)/a)+c
    return err

#This code calculates the number density profile of passive and active particles across the simulation box from a lammpstrj file. NB I think there is potential to mess up in the iteration over t check this.
sigma=1.0


L=string.atoi(sys.argv[1]); #convert strign to int
Lx=string.atoi(sys.argv[2]);   
filename=sys.argv[3]; #file without extension
pe = string.atof(sys.argv[4]);
tau = string.atof(sys.argv[5]);
phi = string.atof(sys.argv[6]);
tmax = string.atoi(sys.argv[7]);
xlo=0;
ylo=0;
xhi=Lx;
yhi=L;


#skip: Number of frames to skip as the trajectory is being read in 
#skip>skip1: Number of frames to skipped as the read in trajectory is analyzed
skip=100;
skip1=2;

#Input file name

trajectory = "{:}.lammpstrj".format(filename)

print "Input file:{0:s}".format(trajectory)

#Commands from the pizza.py librbary. VEry usefull to read lammps trajectories. Will be usefull for the future.
d=dump(trajectory);
d.sort()
time=d.time()
dens=np.zeros(Lx)
densVar=np.zeros(Lx)
frac = np.zeros(Lx)
fracVar=np.zeros(Lx)
count=0
sumx=0
sumy=0


#look at displacements ..and subtract the total netdisplacement. 
grandsumx=0
grandsumy=0


for t in time[:tmax]:
      if count>skip and count%skip1==0:
           #elegant module in pizza.py library. Easy way to process the lammpstrj files. d.vecs() goes frame by frame.  
           idlist,typelist,xlist,ylist,zlist=d.vecs(t,"id","type","x","y","z")
           a = np.array(len(typelist))
           a.fill(2)
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

               #output trajectory after CM fix to double-check

               ## checkfile.write("%d\n%f\n"%(len(idlist),t))
               ## for i in range(len(idlist)):
               ##     if typelist[i]==1:
               ##        checkfile.write("%s\t%f\t%f\t%f\n"%("O",xlist[i],ylist[i],0.0))
               ##     if typelist[i]==1:
               ##        checkfile.write("%s\t%f\t%f\t%f\n"%("H",xlist[i],ylist[i],0.0))
               #print t,cmx, cmy, sumx, sumy,grandsumx


             #End of part copied from Suri's coarse=graining code.



               denstemp,binedges = np.histogram(xlist,bins=Lx,weights=-(np.subtract(typelist,a))) # gives a weight of 1 to particle type 1, weight of 0 to particle type 2. 
               denstemp2,binedges = np.histogram(xlist,bins=Lx) #gets total particle number
               denstemp3 = np.divide(denstemp,denstemp2) #divide type 1 by total - denstemp3 
               dens = np.add(denstemp2,dens) #add total # of particles to average
               densVar = np.add(np.square(dens),densVar) #standard dev
               frac = np.add(denstemp3,frac) #add fraction of type 1 particles to average
               fracVar = np.add(np.square(denstemp3),fracVar) #standard dev
      count=count+1

factor = ((len(time[:tmax])-skip)/skip1)
             
dx=float(Lx/float(len(dens)))
dens = dens/factor
frac = frac/factor

densVar = np.sqrt((densVar/factor)-np.square(dens))/math.sqrt(factor)
fracVar = np.sqrt((fracVar/factor)-np.square(frac))/math.sqrt(factor)

x=xlo+(dx/2.0)

name="densityprof_tmax{:}_{:}".format(tmax,filename)
outputfile = "{:}.stats".format(namefiles.output(name))

fstats=open(outputfile,"w")  
for i in range(len(dens)):
 fstats.write("%f\t%f\t%f\t%f\t%f\n"%(x,dens[i],densVar[i],frac[i],fracVar[i]))
 x+=dx
fstats.close()

outputfile = "InterfaceWidth_4.stats"
of=open(outputfile,"a")

x=np.arange(xlo+(dx/2.0),xhi,dx)

#plt.plot(x,dens)
#plt.show()

p0 = [3.0, 0.5, -0.5, 25]

print "Fitting..."

z,pcov = curve_fit(func,x[:len(x)/2],frac[:len(x)/2],sigma = fracVar[:len(x)/2],absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
#print plsq[2]

#print z[0],z[1],z[2],z[3]
#plt.plot(x[:len(x)/2],func(x[:len(x)/2],z[0],z[1],z[2],z[3]))
#plt.show()

of.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(L,pe,tau,phi,tmax,z[0],perr[0],z[1],perr[1],z[2],perr[2],z[3],perr[3]))

of.close()
