#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 14:06:56 2017

@author: punit
"""

import matplotlib.pyplot as plt
import numpy as np
corr=np.loadtxt('correlation.txt')
time=np.loadtxt('time.txt')

plt.plot([],'r--', label="Correlation Function")
#for col in range(0,10):
#    if(col%2==0):    
plt.plot(time[0:2000],corr[0:2000,0],'r')

plt.legend()
plt.xlabel('Time')
plt.title(' Correlation Plot')