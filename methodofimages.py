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
xs=1.0
ys=1.0

source=af.source(st,xs,ys)
sourceimage=af.source(st,xs,-ys)

vortex1=af.vortex(st,xs,ys)
vorteximage1=af.vortex(-st,xs,-ys)

vortex2=af.vortex(-st,-xs,ys)
vorteximage2=af.vortex(st,-xs,-ys)

source.vel(X,Y)
sourceimage.vel(X,Y)

vortex1.vel(X,Y)
vorteximage1.vel(X,Y)
vortex2.vel(X,Y)
vorteximage2.vel(X,Y)

source.psi(X,Y)
sourceimage.psi(X,Y)

vortex1.psi(X,Y)
vorteximage1.psi(X,Y)
vortex2.psi(X,Y)
vorteximage2.psi(X,Y)

"""u=source.u+sourceimage.u
v=source.v+sourceimage.v
psi=source.psi+sourceimage.psi"""

u=vortex1.u+vortex2.u+vorteximage1.u+vorteximage2.u
v=vortex1.v+vortex2.v+vorteximage1.v+vorteximage2.v
psi=vortex1.psi+vortex2.psi+vorteximage1.psi+vorteximage2.psi

size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter([vortex1.x,vorteximage1.x],[vortex1.y,vorteximage1.y],c='#CD2305',s=80,marker='o')
plt.scatter([vortex2.x,vorteximage2.x],[vortex2.y,vorteximage2.y],c='#CD2305',s=80,marker='o')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);
plt.axis('equal')
plt.show()
