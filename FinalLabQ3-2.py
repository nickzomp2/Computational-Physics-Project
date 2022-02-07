#####Code for question #3a,b #####
#Essentially the same code as #2, but with a smaller domain (99 vs 100)

import numpy as np
import matplotlib.pyplot as plt
from math import tanh, e
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys
from time import clock

start=clock()

#constants
g=9.8 #gravitational acceleration
Nx=99 #number of points in x direction
Ny=99 #number of points in y direction
dx=1.0 #spacing in x
dy=1.0 #spacing in y-
dt=1.0/600.0
it=1500# number of  iteraions
tf=1.0
time=np.arange(0.0,tf,dt)
cx=dt/dx
cy=dt/dy

#empty arrays of Nx x Ny size
b=np.zeros((it+1,int(Nx/dx)+1,int(Ny/dy)+1))
i=np.copy(b)
M=np.copy(b)
N=np.copy(b)
h=np.copy(b)
D=np.copy(b)
s=np.copy(b)

#coordinate system set-up
X=np.arange(0.0,Nx+dx,dx)
Y=np.arange(0.0,Ny+dy,dy)
X,Y=np.meshgrid(X,Y)

n=0.0
m=0.0
L=np.loadtxt('longisland1.txt') #import grid of height values

#determining initial conditions & calculating depth profile function
for m in range(0,int(Nx/dx)+1): #m=y
    for n in range(0,int(Ny/dy)+1): #n=x
        b[0,m,n]=L[m,n]-min(L.flatten())+1
        h[0,m,n]=(10.0*e**((-((n-50.0)**2)/20.0)-(((m-50.0)**2)/20.0)))-min(L.flatten())+1
        s[0][m][n]=np.abs(min(L.flatten())+1)
        D[0,m,n]=h[0,m,n]+b[0,m,n]
        M[0,m,n]=0*np.sqrt(g*h[0,m,n])
        
print np.abs(min(L.flatten())+1)
        
#boundary conditions
M[0,:,0]=0.0 #equals 0 at boundaries
N[0,:,0]=0.0
M[0,:,len(X)-1]=0.0
N[0,:,len(Y)-1]=0.0
M[0,0,:]=0.0 
N[0,0,:]=0.0
M[0,len(X)-1,:]=0.0
N[0,len(Y)-1,:]=0.0

#################################################################################
#viewing initial depth profile (part of 3a)
fig1 = plt.figure(figsize=(13.66,7.68), dpi=120)
ax = fig1.gca(projection='3d')
surf=ax.plot_surface(X, Y, b[0], rstride=1, cstride=1,linewidth=0)
surf=ax.plot_surface(X, Y, s[0], rstride=10, cstride=10,linewidth=0,alpha=0.5)
ax.set_zlim(0, 30)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.view_init(elev=20, azim=150)
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
plt.title("Long Island Depth Profile")
plt.show()
################################################################################
#comment below code to run 3a#

#to run 3b, comment out above flot of fig1, between the hashtags

#####STEPPING#####

#Equations 4-6 of the lab, in code form
def Mf(t,i,j,M,N,h):
    A=0.5*cx*(((M[t,i+1,j]**2)/D[t,i+1,j])-((M[t,i-1,j]**2)/D[t,i-1,j])+(0.5*g*((h[t,i+1,j]**2)-(h[t,i-1,j]**2))))
    B=0.5*cy*((M[t,i,j+1]*N[t,i,j+1]/D[t,i,j+1])-(M[t,i,j-1]*N[t,i,j-1]/D[t,i,j-1]))
    C=g*b[0,i,j]*0.5*cx*(h[t,i+1,j]-h[t,i-1,j])
    return M[t,i,j]-A-B-C
    
def Nf(t,i,j,M,N,h):
    E=0.5*cy*(((N[t,i,j+1]**2)/D[t,i,j+1])-((N[t,i,j-1]**2)/D[t,i,j-1])+(0.5*g*((h[t,i,j+1]**2)-(h[t,i,j-1]**2))))
    F=0.5*cx*((M[t,i+1,j]*N[t,i+1,j]/D[t,i+1,j])-(M[t,i-1,j]*N[t,i-1,j]/D[t,i-1,j]))
    G=g*b[0,i,j]*0.5*cx*(h[t,i,j+1]-h[t,i,j-1])
    return N[t,i,j]-E-F-G
                
def hf(t,i,j,M,N,h):
    return h[t,i,j]-(0.25*cx*(M[t,i+1,j]-M[t,i-1,j]))-(0.25*cy*(N[t,i,j+1]-N[t,i,j-1]))-\
    (0.25*cx*(M[t+1,i+1,j]-M[t+1,i-1,j]))-(0.25*cy*(N[t+1,i,j+1]-N[t+1,i,j-1]))
    
t=0
i=0
j=0

#indexing for each t, for each x, for each y value
for t in range(0,it):
    for i in range(1,len(X)-1):
        for j in range(1,len(Y)-1):
            
            
            #reinforce boundary conditions
            M[t,:,0]=0.0 #equals 0 at boundaries
            N[t,:,0]=0.0
            M[t,:,len(Y)-1]=0.0
            N[t,:,len(Y)-1]=0.0
            M[t,0,:]=0.0 #equals 0 at boundaries
            N[t,0,:]=0.0
            M[t,len(X)-1,:]=0.0
            N[t,len(X)-1,:]=0.0
            
            #determining M,N at t+1
            M[t+1,i,j]=Mf(t,i,j,M,N,h)
            N[t+1,i,j]=Nf(t,i,j,M,N,h)
            
            
            if np.isnan(M[t+1][i][j])==True: #nan checker
                print t, i, j
                sys.exit("Error message A")
            if np.isnan(N[t+1][i][j])==True: #nan checker
                print t, i, j
                sys.exit("Error message B")
                
#calculating h at t+1 from M,N values at t,t+1
    for i in range(0,len(X)-1):

        for j in range(0,len(Y)-1):
            
           #forward/backwards difference calculations for h value at boundaries

            if i==0:
                h[t+1,i,j]=h[t,i,j]-(0.25*cx*(M[t,i+1,j]-M[t,i,j]))-(0.25*cy*(N[t,i,j+1]-N[t,i,j]))-\
    (0.25*cx*(M[t+1,i+1,j]-M[t+1,i,j]))-(0.25*cy*(N[t+1,i,j+1]-N[t+1,i,j]))
    
                h[t+1,i+1,j]=h[t+1,i,j]
                                
            if i==len(X)-2:
                h[t+1,99,j]=h[t,99,j]-(0.25*cx*(M[t,99,j]-M[t,98,j]))-(0.25*cy*(N[t,99,j]-N[t,99,j-1]))-\
    (0.25*cx*(M[t+1,99,j]-M[t+1,98,j]))-(0.25*cy*(N[t+1,99,j]-N[t+1,99,j-1]))
                
            if j==0:
                h[t+1,i,j]=h[t,i,j]-(0.25*cx*(M[t,i+1,j]-M[t,i,j]))-(0.25*cy*(N[t,i,j+1]-N[t,i,j]))-\
    (0.25*cx*(M[t+1,i+1,j]-M[t+1,i,j]))-(0.25*cy*(N[t+1,i,j+1]-N[t+1,i,j]))
                
            if j==len(Y)-2:
                h[t+1,i,j]=h[t,i,j]-(0.25*cx*(M[t,i,j]-M[t,i-1,j]))-(0.25*cy*(N[t,i,j]-N[t,i,j-1]))-\
    (0.25*cx*(M[t+1,i,j]-M[t+1,i-1,j]))-(0.25*cy*(N[t,i,j]-N[t,i,j-1]))
                
                        
            else:
                #calculating h at t+1
                h[t+1,i,j]=hf(t,i,j,M,N,h) 
                
            if np.isnan(h[t+1][i][j])==True:
                print t, i, j
                sys.exit("Error message E")
            
            #calculating D at t+1
            D[t+1,i,j]=h[t+1,i,j]+b[0,i,j] 
    
    #reinforce conditions of D at boundary
    D[t+1,0,:]=h[t+1,0,:]+b[0,0,:] 
    D[t+1,:,0]=h[t+1,:,0]+b[0,:,0] 
    D[t+1,len(X)-1,:]=h[t+1,len(X)-1,:]+b[0,len(X)-1,:] 
    D[t+1,:,len(Y)-1]=h[t+1,:,len(Y)-1]+b[0,:,len(Y)-1] 
    
    #reinforce bondary conditions of h
    h[t+1,:,99]=h[t+1,:,98]
    h[t+1,:,0]=h[t+1,:,1]
    h[t+1,0,:]=h[t+1,1,:]#
    h[t+1,98,:]=h[t+1,99,:]#
    
   
    if t%20.0==0.0:    
        print t
    
print ''


######PLOTTING#####
frame=10.0 #number of plots to skip

for f in range(0,it,1):
    if f%frame == 0.0:
        
        fig = plt.figure(figsize=(13.66,7.68), dpi=120)
        ax = fig.gca(projection='3d')
        surf=ax.plot_surface(X, Y, b[0], rstride=5, cstride=5,linewidth=0)
        surf=ax.plot_surface(X, Y, h[f], rstride=1, cstride=1,cmap=cm.coolwarm,linewidth=0)
        ax.set_zlim(0, 30)
        ax.view_init(elev=20, azim=150)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        plt.savefig('LI' + str(f) + '.png')
        plt.title('Shallow Water Wave Propagation (t=' + str(f*dt) + ')')
        plt.show()
        plt.close()
        print "fig: ", f

print ''
end=clock()
diff=end-start
print 'Run time: ', diff, ' sec'