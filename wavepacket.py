#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:44:49 2017

@author: punit
"""

import matplotlib.pyplot as plt
import numpy as np
wavepacket=np.loadtxt('newwave.txt')
eigenvectors1=np.loadtxt('eigenvectors1.txt')
x=np.linspace(0,10,400)
plt.plot(x,wavepacket,'r--',label="|w0>")
plt.plot(x,eigenvectors1[:,0],'bo',label="|v0>")
plt.legend()
plt.xlabel('Mesh Points')
plt.title(' Wave Packet formed after Excitation and the Initial Wavepacket')