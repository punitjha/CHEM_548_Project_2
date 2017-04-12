#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 13:33:45 2017

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
corr2=np.loadtxt('corr2.txt')
time=np.loadtxt('time2.txt')

plt.plot([],'r--', label="Correlation Function")
#for col in range(0,10):
#    if(col%2==0):    
plt.plot(time[0:100],corr2[0:100,0],'r')

plt.legend()
plt.xlabel('Time')
plt.title(' Correlation Plots')