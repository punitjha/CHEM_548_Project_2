
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack


ft=np.loadtxt('FT2.txt')
print(len(ft))
N=len(ft)
ts=100/N
fs=1/ts
fs1=-fs/2
fs2=fs/2
fs3=fs/N
print(fs, fs1, fs2, fs3)
#list2 = np.arange(-500,500,fs3)
list2=np.linspace(0,550,10000)
#print(len(list2))
#plt.plot((list2/556)*100 ,ft[0:50000])
plt.plot(list2,ft[0:10000])
plt.title(' Fourier Tranform Plot')
plt.figure()
eigenvalues1=np.loadtxt('eigen1.txt')
eigenvalues2=np.loadtxt('eigen2.txt')
frank=np.loadtxt('frank.txt')
#print (eigenvalues1[0])
#print (frank[0,:])
for x in range(len(eigenvalues2)):
    eigenvalues2[x]=eigenvalues2[x]-eigenvalues1[0]
#print (eigenvalues2)
plt.bar(eigenvalues2[0:26],frank[0:26])
plt.xlabel('Energy Difference')
plt.title(' Frank Condon Factors (to the 26 bound states)')
#plt.xcorr(eigenvalues2, frank, normed=True, usevlines=True, maxlags=10, hold=None, data=None)
#plt.title("Gaussian Histogram")
#plt.xlabel("Value")
#plt.ylabel("Frequency")
plt.show()