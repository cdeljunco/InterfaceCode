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

#This code plots the displacement vectors (based on it's displacement from one timestep to the next)

filename=sys.argv[1];
frame=string.atoi(sys.argv[2])
xhi=string.atoi(sys.argv[3])
yhi=string.atoi(sys.argv[4])

trajectory = "{:}.lammpstrj".format(filename)

print "Input file:{0:s}".format(trajectory)

d=dump(trajectory);
d.sort()
time=d.time()

oldidlist,oldtypelist,oldxlist,oldylist,oldzlist=d.vecs(time[frame-1],"id","type","x","y","z")
idlist,typelist,xlist,ylist,zlist=d.vecs(time[frame],"id","type","x","y","z")

dispx=np.zeros(len(idlist))
dispy=np.zeros(len(idlist))

for i in range (len(idlist)):

      if typelist[i]!=oldtypelist[i]:
        print typelist[i],oldtypelist[i]

      dispx[i]=xlist[i]-oldxlist[i]
      dispy[i]=ylist[i]-oldylist[i]

      if dispx[i]>0.5*xhi:
       dispx[i]=dispx[i]-xhi
      if dispx[i]<-0.5*xhi: 
       dispx[i]=dispx[i]+xhi
      
      if dispy[i]>0.5*yhi:
       dispy[i]=-dispy[i]+yhi
      if dispy[i]<-0.5*yhi:
       dispy[i]=-dispy[i]-yhi 


plt.figure()

plt.quiver(xlist,ylist,dispx,dispy)

fig.set_size_inches(10,10)

output="vfield.png"

fig.savefig(output,dpi=500)
