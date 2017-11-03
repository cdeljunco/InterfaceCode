#!/usr/bin/python

import sys
import string
import math

lList = [80,160,320,640]

x= " "
print "L", 10*x, "Timesteps", 10*x, "Time", 10*x
for L in lList:
    print L, 10*x, 100000000*(L**2)/(40**2), 10*x, 100000*(L**2)/(40**2) 



