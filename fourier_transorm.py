#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 15:36:59 2017

@author: punit
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 14:06:56 2017

@author: punit
"""

import matplotlib.pyplot as plt
import numpy as np
ft=np.loadtxt('FT1.txt')
eigenvalues1=np.loadtxt('eigen1.txt')
eigenvalues2=np.loadtxt('eigen2.txt')
for x in range(len(eigenvalues2)):
    eigenvalues2[x]=eigenvalues2[x]-eigenvalues1[0]
z=np.linspace(0,100000,len(ft))
plt.plot(z,ft)
plt.xlabel('Energy Difference')
plt.title(' Fourier Tansform')
plt.show()


plt.figure()
eigenvalues1=np.loadtxt('eigen1.txt')
eigenvalues2=np.loadtxt('eigen2.txt')
frank=np.loadtxt('frank.txt')
#print (eigenvalues1[0])
#print (frank[0,:])
for x in range(len(eigenvalues2)):
    eigenvalues2[x]=eigenvalues2[x]-eigenvalues1[0]
#print (eigenvalues2)
plt.bar(eigenvalues2[0:32],frank[0:32])
plt.xlabel('Energy Difference')
plt.title(' Frank Condon Factors')
#plt.xcorr(eigenvalues2, frank, normed=True, usevlines=True, maxlags=10, hold=None, data=None)
#plt.title("Gaussian Histogram")
#plt.xlabel("Value")
#plt.ylabel("Frequency")
plt.show()