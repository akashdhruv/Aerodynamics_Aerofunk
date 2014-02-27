import numpy as np
import matplotlib.pyplot as plt
from math import *
import aerofunk as af

N = 200                           # Number of points in each direction
xStart,xEnd = -4.0,4.0            # x-direction boundaries
yStart,yEnd = -2.0,2.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid 

Uinf = 1.0        
alphaInDegrees = 0.0       
alpha = alphaInDegrees*pi/180


uFreestream = Uinf*cos(alpha)
vFreestream = Uinf*sin(alpha)


psiFreestream = + Uinf*cos(alpha)*Y - Uinf*sin(alpha)*X



xd,yd=0.0,0.0
mu=5.0
ga=5.0
xv,yv,=0.0,0.0

doublet=af.doublet(mu,xd,yd)
vortex=af.vortex(ga,xv,yv)

doublet.vel(X,Y)
doublet.psi(X,Y)

vortex.vel(X,Y)
vortex.psi(X,Y)

R = sqrt(mu/(2*pi*Uinf))

print ga/(4*pi*Uinf*R)


u=doublet.u+uFreestream+vortex.u
v=doublet.v+vFreestream+vortex.v
psi=doublet.psi+psiFreestream+vortex.psi

# plotting

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xd,yd,c='#CD2305',s=80,marker='o')
circle = plt.Circle((0,0),radius=R,color='#CD2305',alpha=0.5)
plt.gca().add_patch(circle)
plt.contour(X,Y,psi,linewidths=2,linestyles='--')
plt.show()
plt.axis('equal')

Cp = 1.0-(u**2+v**2)/Uinf**2

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='--')
plt.show()


theta=np.linspace(0,2*pi,50)
utheta=-2*Uinf*np.sin(theta)-ga/(2*pi*R)
Cpt= 1 - (utheta/Uinf)**2
plt.figure()
plt.plot(theta,Cpt)
plt.show()