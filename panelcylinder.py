import numpy as np
import matplotlib.pyplot as plt
from math import *
import aerofunk as af
from scipy import integrate

Uinf=1

Np=100
r=1

def I(pi,pj):
    def func(s):
        return (+(pi.xc-(pj.xa-sin(pj.beta)*s))*cos(pi.beta)\
				+(pi.yc-(pj.ya+cos(pj.beta)*s))*sin(pi.beta))\
			   /((pi.xc-(pj.xa-sin(pj.beta)*s))**2\
			   + (pi.yc-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]

def J(pi,pj):
    def func(s):
        return (-(pi.xc-(pj.xa-sin(pj.beta)*s))*sin(pi.beta)\
				+(pi.yc-(pj.ya+cos(pj.beta)*s))*cos(pi.beta))\
			   /((pi.xc-(pj.xa-sin(pj.beta)*s))**2\
			   + (pi.yc-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]


xb=r*np.cos(np.linspace(0,2*pi,Np+1))
yb=r*np.sin(np.linspace(0,2*pi,Np+1))

panel=np.empty(Np,dtype=object)

for i in range(Np):
    panel[i]=af.panel(xb[i],yb[i],xb[i+1],yb[i+1])

plt.figure()
plt.plot(xb,yb)
plt.scatter([p.xa for p in panel],[p.ya for p in panel],c='#CD2305')
plt.scatter([p.xc for p in panel],[p.yc for p in panel],c='k')
plt.axis('equal')
plt.show()

A=np.zeros((Np,Np),dtype=float)

for i in range(Np):
    for j in range(Np):
        if(i==j):
            A[i][j]=0.5
        else:
            A[i][j]=0.5/pi*I(panel[i],panel[j])

b = -Uinf*np.cos([p.beta for p in panel])

sig=np.linalg.solve(A,b)

for i in range(Np):
    panel[i].sigma=sig[i]


for i in range(Np):
    for j in range(Np):
        if (i!=j):
            A[i,j] = 0.5/pi*J(panel[i],panel[j])
        else:
            A[i][j]=0
            
B = -Uinf*np.sin([p.beta for p in panel])
sigma = np.array([p.sigma for p in panel])
Vt = np.dot(A,sigma) + B
for i in range(Np):
    panel[i].Vt = Vt[i]

for i in range(Np):
    panel[i].Cp = 1 - (panel[i].Vt/Uinf)**2
    
    

plt.figure()
plt.scatter([p.xc for p in panel],[p.Cp for p in panel])
plt.show()