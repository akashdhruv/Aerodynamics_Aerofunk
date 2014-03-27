import numpy as np
import matplotlib.pyplot as plt
from math import *
import aerofunk as af
from scipy import integrate


"""Defining function NACA to get airfoil co-ordinates"""

def NACA(naca,c,n):
    
    H=float(naca[0])/100
    p=float(naca[1])/10
    t1=float(naca[2])
    t2=float(naca[3])
    T=10*t1+t2
    T=T/100

    beta=np.linspace(0,pi,n+1)
    xc=np.zeros(np.size(beta))
    for i in range(n+1):
        xc[i]=c*(1-.5*(1-cos(beta[i])))
    
    thdis=np.zeros(np.size(xc))
    
    for i in range(n+1):
        thdis[i]=5*T*c*(0.2969*sqrt(xc[i]/c)-0.126*xc[i]/c-0.3537*(xc[i]/c)**2 +0.2843*(xc[i]/c)**3-0.1015*(xc[i]/c)**4)
    
    camberline=np.zeros(np.size(beta))
    
    if(p!=0.0 and H!=0.0):
        for i in range(n+1):
            if(xc[i] <= p*c):
                camberline[i]=(H/p**2)*xc[i]*(2*p-xc[i]/c)
            elif(xc[i] > p*c):
                camberline[i]=(H/(1-p)**2)*(c-xc[i])*(1+xc[i]/c-2*p)
    
    xu=np.zeros(np.size(xc))
    xl=np.zeros(np.size(xc))
    zu=np.zeros(np.size(xc))
    zl=np.zeros(np.size(xc))
    tht=np.zeros(np.size(xc))
                
    if(p==0 or H==0):
        xu=xc
        zu=thdis
        xl=xc
        zl=-thdis
    else:
        for i in range(n+1):
            if(xc[i] <= p*c):
                tht[i]=atan((2*H/p)*(-xc[i]/(c*p)+1))
            elif(xc[i] > p*c):
                tht[i]=atan((2*H/(1-p**2))*(p-(xc[i]/c)))
            xu[i]=xc[i]-thdis[i]*sin(tht[i])
            zu[i]=camberline[i]+thdis[i]*cos(tht[i])
            xl[i]=xc[i]+thdis[i]*sin(tht[i])
            zl[i]=camberline[i]-thdis[i]*cos(tht[i])
        
        
    X=np.zeros((n+n+1,1),dtype=float)
    Z=np.zeros((n+n+1,1),dtype=float)
    for i in range(n+1):
        X[i]=xl[i]
        Z[i]=zl[i]
    

    for i in range(n):
        X[n+1+i]=xu[n-i-1]
        Z[n+1+i]=zu[n-i-1]
        
    return X,Z


"""Class to store freestream values"""

class freestream:
    def __init__(self,uinf,alpha):
        self.uinf=uinf
        self.alpha=alpha*pi/180

"""Defining function I to get influence coeffcients"""
def I(xc,yc,pj,dxdz,dydz):
	def func(s):
		return (+(xc-(pj.xa-sin(pj.beta)*s))*dxdz+(yc-(pj.ya+cos(pj.beta)*s))*dydz)/((xc-(pj.xa-sin(pj.beta)*s))**2+(yc-(pj.ya+cos(pj.beta)*s))**2)
	return integrate.quad(lambda s:func(s),0.,pj.length)[0]

"""Function to bulid matrix A"""

def buildA(p):
    N=len(p)
    A=np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            if(i!=j):
                A[i][j]=0.5/pi*I(p[i].xc,p[i].yc,p[j],cos(p[i].beta),sin(p[i].beta))
            elif(i==j):
                A[i][j]=0.5

    return A

"""Function to build matrix rhs"""
def buildrhs(p,fs):
    N=len(p)
    Rhs=np.zeros((N),dtype=float)
    for i in range(N):
		Rhs[i] = -fs.uinf*cos(fs.alpha-p[i].beta)
    return Rhs

"""Function to get tangential velocity"""

def tvel(p,fs):
    N = len(p)
    A = np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            if (i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),cos(p[i].beta))
    Rhs = fs.uinf*np.sin([fs.alpha-pp.beta for pp in p])
    var = np.array([pp.sigma for pp in p])
    vt = np.dot(A,var)+Rhs
    for i in range(N):
		p[i].vt = vt[i]

"""Function to get Cp"""
def cp(p,fs):
	for i in range(len(p)):
		p[i].Cp = 1-(p[i].vt/fs.uinf)**2

"""_____________________________MAIN___________________________"""

N=40
c=1
naca="0012"

xp,yp=NACA(naca,c,N)

M=len(xp)

panel=np.empty(M-1,dtype=object)

for i in range(M-1):
    panel[i]=af.panel(xp[i],yp[i],xp[i+1],yp[i+1])


plt.figure
plt.plot(xp,yp)
plt.scatter([p.xc for p in panel],[p.yc for p in panel])
plt.axis("equal")
plt.show()

uinf=1.0
alpha=5.0

freestream=freestream(uinf,alpha)

A=buildA(panel)
Rhs=buildrhs(panel,freestream)

var = np.linalg.solve(A,Rhs)
for i in range(len(panel)):
	panel[i].sigma = var[i]

tvel(panel,freestream)
cp(panel,freestream)

