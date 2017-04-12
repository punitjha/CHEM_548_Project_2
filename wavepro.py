#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:49:54 2017

@author: punit
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:44:49 2017

@author: punit
"""

import matplotlib.pyplot as plt
import numpy as np
wavepacket=np.loadtxt('wavepro1.txt')
print(len(wavepacket))
x=np.linspace(-100,300,400)
plt.plot([],'r--', label="Real")
plt.plot([],'b--', label="Imaginary")
plt.plot(x,wavepacket[:,0]+0 ,'r--')
plt.plot(x,wavepacket[:,1]+0 ,'b--')
plt.plot(x,wavepacket[:,10]+0.5,'r--')
plt.plot(x,wavepacket[:,11]+0.5,'b--')
plt.plot(x,wavepacket[:,22]+1.0,'r--')
plt.plot(x,wavepacket[:,23]+1.0,'b--')
plt.plot(x,wavepacket[:,30]+1.5 ,'r--')
plt.plot(x,wavepacket[:,31]+1.5 ,'b--')



axes = plt.gca()
axes.set_yticklabels([-0.05,0, 0.01, 0.22,0.30 ])
plt.legend()
plt.xlabel('Mesh Points')
plt.ylabel('Time ')
plt.title(' Time Evolution of the wavepacket -- Finite difference method')