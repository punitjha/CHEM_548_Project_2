#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 13:25:05 2017

@author: punit
"""
import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('eigen1.txt')
eigenvalues2=np.loadtxt('eigen2.txt')
frank=np.loadtxt('frank.txt')
#print (eigenvalues1[0])
#print (frank[0,:])
for x in range(len(eigenvalues2)):
    eigenvalues2[x]=eigenvalues2[x]-eigenvalues1[0]
#print (eigenvalues2)
plt.bar(eigenvalues2[0:50],frank[0:50])

#plt.xcorr(eigenvalues2, frank, normed=True, usevlines=True, maxlags=10, hold=None, data=None)
#plt.title("Gaussian Histogram")
#plt.xlabel("Value")
#plt.ylabel("Frequency")
plt.show()