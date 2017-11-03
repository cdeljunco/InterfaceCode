#!/usr/bin/python

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
  

#This code calculates the velocity profile across the system of the y-compoennt of the velocity.  First, the displacement of the CM is fixed.  Then, the average y-direction velocity of each particle type and the total is computed in its bulk region. 
#cutoff= coarse graining lenght scale 
#L= Size of box along y dimension
#Lx= Size of box along x dimension 

L=string.atoi(sys.argv[1]);
Lx=string.atoi(sys.argv[2]);
cg=string.atof(sys.argv[3]);
filename=sys.argv[4];
xlo=0;
ylo=0;
xhi=Lx;
yhi=L;
 
#skip: Number of frames to skip as the trajectory is being read in 
#skip1>skip: Number of frames to skipped as the read in trajectory is analyzed 
skip=0;
skip1=1;
collect=100;

#Input file name 


trajectory = "{:}.lammpstrj".format(filename)

print "Input file:{0:s}".format(trajectory)

#Commands from the pizza.py librbary. VEry usefull to read lammps trajectories. Will be usefull for the future. 
d=dump(trajectory);
d.sort()
time=d.time()
vytype1=np.zeros(np.divide(xhi,cg))
vytype2=np.zeros(np.divide(xhi,cg))
vytype1_std=np.zeros(np.divide(xhi,cg))
vytype2_std=np.zeros(np.divide(xhi,cg))
countersum1=np.zeros(np.divide(xhi,cg))
countersum2=np.zeros(np.divide(xhi,cg))
count=0
sumx=0
sumy=0

#check= "dump/checktrajectory_{0:s}.XYZ".format(filename)
#checkfile = open(check,"w")

#Output file names 

name="vprofUnnormalized_{:}".format(filename)
outputfile="{:}.stats".format(namefiles.output(name))

foutput=open(outputfile,"w")

#look at displacements
grandsumx=0
grandsumy=0

for t in time:
  if count>skip and count%skip1==0 and (count-skip)/skip1<collect:
   #elegant module in pizza.py library. Easy way to process the lammpstrj files. d.vecs() goes frame by frame.  
   idlist,typelist,xlist,ylist,zlist=d.vecs(t,"id","type","x","y","z")
   cmx=0
   cmy=0;
   counter=0;
   type=0;
   #Compute center of mass 
   for i in range (len(idlist)):
    if typelist[i]==1:
     cmx+=xlist[i];
     cmy+=ylist[i];
     counter=counter+1;
   cmx=cmx/counter;
   cmy=cmy/counter;
   #compute displacement of each atom. Lammps scrambles the output. The i^th atom in the frame need not have the label i. the d.sorrt() routine takes care of this 
   dispx=np.zeros(len(idlist))
   dispy=np.zeros(len(idlist))
   if count==skip1:
    oldxlist2=np.zeros(len(idlist))
    oldylist2=np.zeros(len(idlist))
   #print t, cmx,cmy
   if count>skip+skip1:
    oldidlist,oldtypelist,oldxlist,oldylist,oldzlist=d.vecs(time[count-skip1],"id","type","x","y","z")
   if count>skip+skip1:      
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
      if typelist[i]==2 or typelist[i]==1:
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
      if typelist[i]==2: #what is going on here??
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
 
      
     #print t,cmx, cmy, sumx, sumy,grandsumx

#Center of mass should now be fixed in 

#output trajectory after CM fix to double-check

#     checkfile.write("%d\n%f\n"%(len(idlist),t))
#     for i in range(len(idlist)):
#            if typelist[i]==1:
#               checkfile.write("%s\t%f\t%f\t%f\n"%("O",xlist[i],ylist[i],0.0))
#            if typelist[i]==2:
#               checkfile.write("%s\t%f\t%f\t%f\n"%("H",xlist[i],ylist[i],0.0))



#Now, calculate the profile of velocities across the box
   if count>skip+(2*skip1):

     x=xlo
     j=0   
     while x<xhi: 
      if x+cg>xhi:
       print "Warning, oob"
      sumy1=0
      sumy2=0
      j=int(np.divide(x,cg))
      for i in range(len(idlist)):
	if typelist[i]==1 and xlist[i]>x and xlist[i]<x+cg:
	 dispy[i]=ylist[i]-oldylist2[i]
         if dispy[i]>0.5*yhi:
          dispy[i]=-dispy[i]+yhi
         if dispy[i]<-0.5*yhi:
          dispy[i]=-dispy[i]-yhi 
         countersum1[j]=countersum1[j]+1
         sumy1=sumy1+dispy[i]
        if typelist[i]==2 and xlist[i]>x and xlist[i]<x+cg:
	 dispy[i]=ylist[i]-oldylist2[i]
         if dispy[i]>0.5*yhi:
          dispy[i]=-dispy[i]+yhi
         if dispy[i]<-0.5*yhi:
          dispy[i]=-dispy[i]-yhi 
         countersum2[j]=countersum2[j]+1
         sumy2=sumy2+dispy[i]
      
      vytype1[j]= vytype1[j] + sumy1 #bin average velocity of each particle typ ein small slice of the lattive
      vytype2[j]= vytype2[j] + sumy2
      vytype1_std[j]=vytype1_std[j] + math.pow(sumy1,2.0)
      vytype2_std[j]=vytype2_std[j] + math.pow(sumy2,2.0)
      x=x+cg

   for i in range(len(idlist)):
      oldxlist2[i]=xlist[i]
      oldylist2[i]=ylist[i]
      
   
  count=count+1

x=xlo+(0.5*cg)

print np.sum(vytype1)
print np.sum(vytype2)

factor=float(collect)
for j in range(len(vytype1)):
   vytype1[j]=vytype1[j]/factor
   if vytype1_std[j]/factor<pow(vytype1[j],2.0):
    print "type 1 error"
   vytype1_std[j]=math.sqrt((vytype1_std[j]/factor)-pow(vytype1[j],2.0))/math.sqrt(factor)
   vytype2[j]=vytype2[j]/factor
   if vytype2_std[j]/factor<pow(vytype2[j],2.0):
    print "type 2 error"
   vytype2_std[j]=math.sqrt((vytype2_std[j]/factor)-pow(vytype2[j],2.0))/math.sqrt(factor)
   foutput.write("%f\t%f\t%f\t%f\t%f\n"%(x,vytype1[j],vytype1_std[j],vytype2[j],vytype2_std[j]))
   x=x+cg
#checkfile.close()
foutput.close()
