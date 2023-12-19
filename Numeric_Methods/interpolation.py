import numpy as np
from numpy.polynomial import polynomial as P
import math
import matplotlib.pyplot as plt
import NMclasses

def LNdepN(x, f, a, b, ns, nf, typeOfGrid):
    Er = list()
    N = list()
    for n in range(ns, nf):
        h = (b-a)/n
        Xi = list()
        Yi = list()
        for i in range(n+1):
            if typeOfGrid == 'Normal':
                Xi.append(a+i*h)
            if typeOfGrid == 'Cheb':
                xi = (b+a)/2+((b-a)/2)*math.cos(((2*i+1)/(2*n+2))*math.pi)
                Xi.append(xi)
            Yi.append(f(Xi[i]))
        sum = 0
        for k in range(len(Xi)):
            p=1
            for j in range(len(Xi)):
                if j!=k :
                    c = np.poly([Xi[j]])/(Xi[k]-Xi[j])
                    p = np.polymul(p,c)
            term = p*Yi[k]
            sum = sum+term

        E = lambda x: np.polyval(sum,x)-f(x)
        Err = (E(x))
        Er.append(np.max(Err))
        N.append(n)
    return Er, N 
    
def LNdepLen(f, a, b, n, typeOfGrid):
    Er = list()
    Len = list()
    for k in range(12):
        a = a-1
        b = b+1
        x = np.linspace(a, b, 500)
        length = b-a
        h = (b-a)/n
        Xi = list()
        Yi = list()
        for i in range(n+1):
            if typeOfGrid == 'Normal':
                Xi.append(a+i*h)
            if typeOfGrid == 'Cheb':
                xi = (b+a)/2+((b-a)/2)*math.cos(((2*i+1)/(2*n+2))*math.pi)
                Xi.append(xi)
            Yi.append(f(Xi[i]))
        sum = 0
        for k in range(len(Xi)):
            p=1
            for j in range(len(Xi)):
                if j!=k :
                    c = np.poly([Xi[j]])/(Xi[k]-Xi[j])
                    p = np.polymul(p,c)
            term = p*Yi[k]
            sum = sum+term

        E = lambda x: np.polyval(sum,x)-f(x)
        Err = E(x)
        Er.append(np.max(Err))
        Len.append(length)
    return Er, Len 

def SPdepN(f, a, b, ns, nf, typeOfGrid):
    E = list()
    N = list()
    for n in range(ns, nf+1):
        spX, spY, sp = NMclasses.interpolation(f, a, b, n, typeOfGrid).splineInterpolation()
        Err = list()
        F = f(spX)
        for i in range(len(F)):
            Err.append(abs(-F[i]+spY[i]))
        E.append(np.max(Err))
        N.append(n)
    return E, N

