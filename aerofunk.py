import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import integrate



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


class doublet:
    def __init__(self,mu,x,y):
        self.mu=mu
        self.x,self.y=x,y
    def vel(self,X,Y):
        self.u = - self.mu/(2*pi)*((X-self.x)**2-(Y-self.y)**2)/((X-self.x)**2+(Y-self.y)**2)**2
        self.v = - self.mu/(2*pi)*2*(X-self.x)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)**2
    def psi(self,X,Y):
        self.psi=-self.mu/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        

class vortex:
    def __init__(self,gam,x,y):
        self.gam=gam
        self.x,self.y=x,y
    def vel(self,X,Y):
        self.u = +self.gam/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        self.v = -self.gam/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
    def psi(self,X,Y):
        self.psi = self.gam/(4*pi)*np.log((X-self.x)**2+(Y-self.y)**2)
        

class source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    def vel(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    def psi(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))

class panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa=xa
        self.ya=ya
        self.xb=xb
        self.yb=yb
        self.xc=(xa+xb)/2
        self.yc=(ya+yb)/2
        self.length=sqrt((xa-xb)**2+(ya-yb)**2)
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        self.sigma=0
        self.Vt=0
        self.Cp=0
