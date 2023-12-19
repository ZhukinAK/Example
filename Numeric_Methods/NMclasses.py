import math
import random
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 

from numpy.polynomial import polynomial as P


class Plots:
    def __init__(self, n, ars, title= None, xlabel=None, ylabel=None, legend = None):          
        self.n = n
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.ars = ars
        self.legend = legend
        
        self.fig = plt.figure(self.n)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid(True) 
        self.ax.set_title((self.title))                                     
        self.ax.set_xlabel((self.xlabel))
        self.ax.set_ylabel((self.ylabel))
        
        
        self.build()
        
    def build(self, scale = 'Normal'):
        colors = ['b', 'g', 'r', 'k', 'm']
        for number in range(len(self.ars)): 
            col = random.choice(colors)
            if self.legend != None:
                self.ax.plot(self.ars[number][0], self.ars[number][1], col, label = self.legend[number])
                if scale == 'logX':
                    self.ax.set_xscale('log')
                if scale == 'logY':
                    self.ax.set_yscale('log')
                if scale == 'loglog':
                    self.ax.set_xscale('log')
                    self.ax.set_yscale('log')
            else:
                self.ax.plot(self.ars[number][0], self.ars[number][1], col)
                if scale == 'logX':
                    self.ax.set_xscale('log')
                if scale == 'logY':
                    self.ax.set_yscale('log')
                if scale == 'loglog':
                    self.ax.set_xscale('log')
                    self.ax.set_yscale('log')
            colors.remove(col)
        self.ax.legend()
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(1))

        
class interpolation: 
    def __init__(self, f, a, b, n, typeOfGrid = None):
        self.a = a
        self.b = b
        self.f = f
        self.n = n
        self.typeOfGrid = typeOfGrid
        
    def splineInterpolation(self):
        h = (self.b-self.a)/self.n
        Xi = list()
        Yi = []
        for i in range(self.n+1):
            if self.typeOfGrid == 'Normal':
                Xi.append(self.a+i*h)
            if self.typeOfGrid == 'Cheb':
                xi = (self.b+self.a)/2+((self.b-self.a)/2)*math.cos(((2*i+1)/(2*self.n+2))*math.pi)
                Xi.append(xi)       
            Yi.append(self.f(Xi[i]))

        splines = self._splining(Xi, Yi)
        x = []
        y = []

        for i in range(1, len(splines)):
            interval = np.linspace(splines[i].xs, splines[i].xf, 20)
            interval = interval.tolist()
            x.append(interval)
            for j in range(20):
                y.append(np.polyval(splines[i].create(), x[i-1][j]))
        x = sum(x, [])
        return x, y, splines

    def _splining(self, Xi, Yi):

        splines = list()
        splines.append(spline(0,0,0,0,0,0))
        for i in range(1, self.n+1):
            splines.append(spline(0,0,0,0,0,0))
            splines[i].a = Yi[i]
            splines[i].xs = Xi[i-1]
            splines[i].xf = Xi[i]

        splines[0].c = 0   ## производная в начале промежутка
        splines[self.n].c = 0 ## в конце
        

        delta = [0 for t in range(self.n)]
        lambd = [0 for t in range(self.n)]

        for i in range(1, self.n):
            hi  = Xi[i+1] - Xi[i]
            hip = Xi[i] - Xi[i-1]
            B = hip
            C = 2*(hi+hip)
            D = hi
            R = 3*((Yi[i+1]-Yi[i])/hi - (Yi[i]-Yi[i-1])/hip)
            delta[i] = -D/(B*delta[i-1]+C)
            lambd[i] = (R - B*lambd[i-1])/(B*delta[i-1]+C)

        for i in range(self.n-1, 0, -1):
            splines[i].c = delta[i] * splines[i+1].c + lambd[i]

        for i in range(1, self.n+1):
            splines[i].d = (splines[i].c-splines[i-1].c)/(3*(Xi[i]-Xi[i-1]))
            splines[i].b = (Yi[i] - Yi[i - 1]) / (Xi[i]-Xi[i-1]) + 1/3*(Xi[i]-Xi[i-1])*splines[i-1].c + 2/3*(Xi[i]-Xi[i-1])*splines[i].c
        return splines

    def lagrange(self):
        self.n = self.n-1
        h = (self.b-self.a)/self.n
        Xi = list()
        Yi = list()
        for i in range(self.n+1):
            if self.typeOfGrid == 'Normal':
                Xi.append(self.a+i*h)
            if self.typeOfGrid == 'Cheb':
                xi = (self.b+self.a)/2+((self.b-self.a)/2)*math.cos(((2*i+1)/(2*self.n+2))*math.pi)
                Xi.append(xi)
            Yi.append(self.f(Xi[i]))
        summ = 0
        for k in range(len(Xi)):
            p=1
            for j in range(len(Xi)):
                if j!=k :
                    c = np.poly([Xi[j]])/(Xi[k]-Xi[j])
                    p = np.polymul(p,c)
            term = p*Yi[k]
            summ = summ+term
        return Xi, Yi, summ

class spline:
    def __init__(self, a, b, c, d, xs, xf):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.xs = xs
        self.xf = xf
    
    def create(self):
        p1 = np.poly([self.xf])
        p2 = np.polymul(p1,p1)
        p3 = np.polymul(p2, p1)
        p32 = np.polyadd(self.d*p3, self.c*p2)
        p321 = np.polyadd(p32, self.b*p1)
        p3210 = np.polyadd(p321, [self.a])
        return p3210
    
class Integtration:
    def __init__(self, func, a, b, n=None):
        self.func = func
        self.a = a
        self.b = b
        self.n = n
        
    def rectangles(self):
        integral = 0.0
        step = (self.b - self.a) / self.n
        X = [self.a+i*step for i in range(self.n)]
        for x in range(len(X)):
            integral += step * self.func(X[x] + step / 2)
        return integral

    def __lobattoAlg(self, start, end):
        integral = 0.0
        t_i = np.array([-1, -math.sqrt(1/5), math.sqrt(1/5), 1])
        weighs = np.array([1/6, 5/6, 5/6, 1/6])
        for i in range(len(weighs)):
            integral += self.func( (end - start)*t_i[i]/2 + (start + end)/2 )*weighs[i]
        integral *= (end - start)/2
        return integral

    def lobattoInt(self):
        integral = 0.0
        step = (self.b - self.a)/self.n
        X = [self.a + i*step for i in range(self.n+1)]
        for i in range(self.n):
            integral += self.__lobattoAlg(X[i], X[i+1])
        return integral



        