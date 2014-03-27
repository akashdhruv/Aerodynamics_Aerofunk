import numpy as np
import matplotlib.pyplot as plt
from math import *
import aerofunk as af
from scipy import integrate

N = 200                           
xStart,xEnd = -2.0,2.0            
yStart,yEnd = -2.0,2.0            
x = np.linspace(xStart,xEnd,N)    
y = np.linspace(yStart,yEnd,N)    
X,Y = np.meshgrid(x,y)            

Uinf = 1.0        
alphaInDegrees = 0.0       
alpha = alphaInDegrees*pi/180


uFreestream = Uinf*cos(alpha)
vFreestream = Uinf*sin(alpha)


psiFreestream = + Uinf*cos(alpha)*Y - Uinf*sin(alpha)*X

sigma = 2.5    # strength of the source sheet

uPanel    = np.empty((N,N),dtype=float)
vPanel    = np.empty((N,N),dtype=float)

# boundaries of the sheet
ymin,ymax = -1.0,1.0

# computing the velocity field
for i in range(N):
    for j in range(N):
        
        func = lambda s : X[i,j]/(X[i,j]**2+(Y[i,j]-s)**2)
        uPanel[i,j] = sigma/(2*pi)*integrate.quad(func,ymin,ymax)[0]
        
        func = lambda s:(Y[i,j]-s)/(X[i,j]**2+(Y[i,j]-s)**2)
        vPanel[i,j] = sigma/(2*pi)*integrate.quad(func,ymin,ymax)[0]
        
u=uPanel+uFreestream
v=vPanel+vFreestream

size = 8
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.axvline(0.0,(ymin-yStart)/(yEnd-yStart),(ymax-yStart)/(yEnd-yStart),\
            color='#CD2305',linewidth=4)
cont = plt.contourf(X,Y,np.sqrt(u**2+v**2),levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(cont)
cbar.set_label('U',fontsize=16);
plt.show()