#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 19:24:51 2017

@author: punit
"""

import matplotlib.pyplot as plt
import numpy as np
wavepacket=np.loadtxt('time1.txt')
print(len(wavepacket))
x=np.linspace(0,400,200)
plt.plot([],'r--', label="Real")
plt.plot([],'b--', label="Imaginary")
for col in range(16):
    if(col%2==0):    
        plt.plot(x,wavepacket[50:250,col]+0.5*col,'r--')
    else:
        plt.plot(x,wavepacket[50:250,col]+0.5*(col-1),'b--')
plt.legend()
axes = plt.gca()
axes.set_yticklabels([-0.005,0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035,0.040])
#axes.set_ylim([-0.5,3])
#plt.yticks(np.arange(min(y), max(y)+1, 10.0))
plt.xlabel('Mesh Points')
plt.ylabel('Time in steps of 0.005s')
plt.title(' Time Evolution of the wavepacket')