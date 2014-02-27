import numpy as np
import matplotlib.pyplot as plt
from math import *
import aerofunk as af

N = 200                           
xStart,xEnd = -4.0,4.0            
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

st=1.0
xs=0.0
ys=1.0

source=af.source(st,xs,ys)
sourceimage=af.source(st,xs,-ys)

vortex=af.vortex(st,xs,ys)
vorteximage=af.vortex(-st,xs,-ys)

source.vel(X,Y)
sourceimage.vel(X,Y)

vortex.vel(X,Y)
vorteximage.vel(X,Y)

source.psi(X,Y)
sourceimage.psi(X,Y)

vortex.psi(X,Y)
vorteximage.psi(X,Y)

"""u=source.u+sourceimage.u
v=source.v+sourceimage.v
psi=source.psi+sourceimage.psi"""

u=vortex.u+vorteximage.u
v=vortex.v+vorteximage.v
psi=vortex.psi+vorteximage.psi

size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex.x,vortex.y,c='#CD2305',s=80,marker='o')
plt.scatter(vorteximage.x,vorteximage.y,c='#CD2305',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);
plt.axis('equal')
plt.show()