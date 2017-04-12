#!/usr/bin/python

import numpy as np

eigs1 = np.loadtxt('eigen1.txt')
eigs2 = np.loadtxt('eigen2.txt')

for i in range(len(eigs1)):
	print(eigs2[i] - eigs1[0])

