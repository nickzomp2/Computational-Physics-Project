#####Code for question #1 #####

import numpy as np
import matplotlib.pyplot as plt
from math import tanh, e
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

g=9.8 #gravitational acceleration
Nx=100 #number of points in x direction
Ny=100 #number of points in y direction
dx=1.0 #spacing in x
dy=1.0 #spacing in y
dt=0.01 
it=1 # number of iteraions
tf=1.0 #end time
time=np.arange(0.0,tf,dt) #time array

#empty arrays of Nx x Ny size
b=np.zeros((int(tf/dt)+1,int(Nx/dx)+1,int(Ny/dy)+1))
i=np.copy(b)
h=np.copy(b)
s=np.copy(b)


#coordinate system
X=np.arange(0.0,Nx+dx,dx)
Y=np.arange(0.0,Ny+dy,dy)
X,Y=np.meshgrid(X,Y)

n=0.0
m=0.0

#determining initial conditions & calculating depth profile function (uncomment plots of choice before running)
for m in range(0,int(Nx/dx)+1): #m=y
    for n in range(0,int(Ny/dy)+1): #n=x
        b[0][m][n]=50.0*tanh((0.069*n)- 5.5)+50.0 # hyperbolic depth profile
        #b[0][m][n]=n #linear depth profile
        #b[0][m][n]=np.sqrt(100*n) #square root profile
        h[0,m,n]=(100.0*e**((-((n-40.0)**2)/20.0)-(((m-50.0)**2)/20.0)))+100 #initial Gaussian conditon of wave #1
        #h[0,m,n]=(20.0*e**((-((n-40.0)**2)/50.0)-(((m-50.0)**2)/20000.0)))+100 #initial Gaussian conditon of wave #2
        s[0][m][n]=100 #represents sea level

#3d plotting --> based off of matplotlib's demo code: surface3d_demo.py
fig1 = plt.figure(figsize=(13.66,7.68), dpi=120)
ax = fig1.gca(projection='3d')
surf=ax.plot_surface(X, Y, b[0], rstride=1, cstride=1,cmap=cm.coolwarm,linewidth=0)
#surf=ax.plot_surface(X, Y, h[0], rstride=1, cstride=1,cmap=cm.coolwarm,linewidth=0)
surf=ax.plot_surface(X, Y, s[0], rstride=10, cstride=10,linewidth=0,alpha=0.5)
ax.set_zlim(0, 200)
ax.view_init(elev=20, azim=240)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#plt.title('Gaussian Initial Wave (A=100.0, sigmax=sigmay=20.0)')
plt.show()