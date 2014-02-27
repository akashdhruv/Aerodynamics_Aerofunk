import numpy as np
import matplotlib.pyplot as plt
from math import *


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