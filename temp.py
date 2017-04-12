import matplotlib.pyplot as plt
import numpy as np

#plt.plotfile('eigenvectors.txt', delimiter=' ')
#x=np.linspace(0,100,100)
#plt.plot(x,y, label='Loaded from file!')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Interesting Graph\nCheck it out')
#plt.legend()
#plt.show()

#m = np.loadtxt("eigenvectors.txt", delimiter="\t")
eigenvectors1=np.loadtxt('eigenvectors1.txt')
eigenvalues1=np.loadtxt('eigen1.txt')
potential1=np.loadtxt('potent1.txt')
eigenvectors2=np.loadtxt('eigenvectors2.txt')
eigenvalues2=np.loadtxt('eigen2.txt')
potential2=np.loadtxt('potent2.txt')
#R=np.zeros(eigenvectors.shape[0])
#eigenvectors=np.vstack([eigenvectors,R])
for col in range(1):
    #print(ary[:,col])
    #print (arr[col])
    ppp=eigenvectors1[:,col]+eigenvalues1[col]
    x=np.linspace(0,10,400)
    plt.plot(x,ppp)
plt.plot(x,potential1,'r--',label="Potential1")
plt.legend()
for col in range(17):
    #print(ary[:,col])
    #print (arr[col])
    ppp=eigenvectors2[:,col]+eigenvalues2[col]
    plt.plot(x,ppp)
plt.plot(x,potential2,'r--',label="Potential2")
plt.legend()
plt.xlabel('Mesh Points')
plt.ylabel('Energy ')
plt.title(' Displaced Morse Potential')
plt.show()
#datafile = open('eigenvectors.txt', 'r') # Open the file to read it
#data = datafile.readlines()
#l = [ row.strip().split() for row in datafile] 
#a=5+l[1][2]
#print (a)

#datafile.close()
##print ( len(data))
##print (data[0])
#N = []
##S = []
##K = []
#i=0
#for row in data:
#    this_data = row.split()
#    N.append(this_data[2])
#    i=i+1
#    plt.plot(N)
##    S.append(float(this_data[1]))
##    K.append(float(this_data[2]))
#print (N)
#x=np.linspace(0,100,100)
#
#plt.plot(x,ary[:,0])
##plt.plot(x,K)
#plt.ylabel('some numbers')
#plt.show()

#c = np.genfromtxt('eigenvectors.txt', names=True, dtype=None)
#print( type(l))
##print (c.dtype)
#print (c['7435446036652e07'])
##x=np.linspace(0,10,100)
##plt.plot(x,c['7435446036652e07'])
##plt.show()