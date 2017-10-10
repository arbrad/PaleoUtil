# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 20:40:03 2017

@author: Aaron
"""

import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats
from scipy.interpolate import Rbf

class Params:
    def __init__(self):
        # grid of zones
        self.X = 10
        self.Y = 10
        # number of individuals
        self.N = 100
        # resolution, e.g., with .25, 16 subzones/zone
        self.Res = 0.25
        # sample points for creating eco-impact functions
        self.CrossExtrema = 10
        # amount of change between zone-specific impact functions
        self.CrossJitter = 0.5
        # width of Gaussian in impact functions
        self.StdDev = 1
        # density-impact % on landscape
        self.Density = -1
        # 1/Head = chance of individual accepting bad move is 1/Bad
        self.Bad = 5
        # whether to reproduce
        self.Reproduce = True

def bound(a, l, u):
    for i in range(len(a)):
        a[i] = min(u, max(l, a[i]))

def grid(params):
    X, Y, Res = params.X, params.Y, params.Res
    return np.meshgrid(np.arange(0, X+Res, Res),
                       np.arange(0, Y+Res, Res))

class Sog:
    '''Sum of Gaussians'''
    def __init__(self, params, x, y, z):
        cov = [[params.StdDev**2, 0], [0, params.StdDev**2]]
        xi, yi = grid(params)
        co = np.zeros((*xi.shape,2))
        for m in range(xi.shape[0]):
            for n in range(xi.shape[1]):
                co[m,n,0] = xi[m,n]
                co[m,n,1] = yi[m,n]
        self.f = np.zeros(xi.shape)
        for i in range(len(x)):
            d = scipy.stats.multivariate_normal([x[i], y[i]], cov)
            self.f += z[i]*d.pdf(co)
        self.f /= len(x)
    def __call__(self, xi, yi):
        return self.f

def cross(params, verb=False):
    X, Y = params.X, params.Y
    N = params.CrossExtrema
    x = np.random.rand(N) * X
    y = np.random.rand(N) * Y
    z = np.random.rand(N) - 0.5
    xi, yi = grid(params)
    rv = [[None for _ in range(Y)] for _ in range(X)]
    for i in range(X):
        for j in range(Y):
            for k in range(len(x)):
                x[k] += random.normalvariate(0, params.CrossJitter)
                y[k] += random.normalvariate(0, params.CrossJitter)
                z[k] += random.normalvariate(0, params.CrossJitter)
                bound(x, 1, X-1)
                bound(y, 1, Y-1)
                bound(z, -1, 1)
            if params.Density != 0:
                x[-1] = i + 0.5
                y[-1] = j + 0.5
                z[-1] = params.Density
            rbv = Sog(params, x, y, z) #Rbf(x, y, z, function='gaussian')
            rv[i][j] = rbv(xi, yi)
            if verb:
                zi = rv[i][j]
                plt.subplot(Y, X, Y*(Y-j-1)+i+1)
                plt.pcolor(xi, yi, zi, cmap=mplt.cm.jet)
                plt.scatter(x, y, c=z, cmap=mplt.cm.jet)
                plt.xlim(0, X)
                plt.ylim(0, Y)
                plt.colorbar()
    return rv

class Sim:
    def __init__(self, params, seed=None):
        random.seed(seed)
        self.params = params
        self.cross = cross(params)
        self.search = [(params.X/2, params.Y/2) for _ in range(params.N)]
        self.scape = np.ones(self.cross[0][0].shape)
        self.scape /= sum(sum(self.scape))
        self.grid = grid(self.params)
    def crosses(self):
        xi, yi = grid(self.params)
        pi = 0
        for j in range(self.params.Y-1,-1,-1):
            for i in range(self.params.X):
                pi += 1
                plt.subplot(self.params.Y, self.params.X, pi)
                plt.pcolor(xi, yi, self.cross[i][j], cmap=mplt.cm.jet)
    def major(self, x, y):
        return int(x), int(y)
    def minor(self, x, y):
        return int(round(y/self.params.Res, 0)), int(round(x/self.params.Res, 0))
    def incrScape(self, factor=1):
        zs = [self.scape[self.minor(x, y)] for x, y in self.search]
        for i, (x, y) in enumerate(self.search):
            X, Y = int(x), int(y)
            self.scape += zs[i]/factor/len(self.search)*self.cross[X][Y]
        self.scape -= self.scape.min()
        #for x, y in self.search:
        #    self.scape[self.minor(x, y)] *= self.params.Density
        self.scape /= sum(sum(abs(self.scape)))
    def incrOpts(self):
        mean, cov = [0, 0], [[self.params.Res**2, 0], [0, self.params.Res**2]]
        for i in range(len(self.search)):
            x, y = self.search[i]
            z0 = self.scape[self.minor(x, y)]
            dx, dy = np.random.multivariate_normal(mean, cov)
            x = min(self.params.X-.0001, max(0, x+dx))
            y = min(self.params.Y-.0001, max(0, y+dy))
            z1 = self.scape[self.minor(x, y)]
            if z1 > z0 or (not random.randint(0, self.params.Bad) and random.uniform(0, z0) < z1):
                self.search[i] = (x, y)
        if self.params.Reproduce:
            N = len(self.search)
            Z = sum(self.scape[self.minor(*co)] for co in self.search)
            for i in range(N):
                z = self.scape[self.minor(*self.search[i])]
                if random.uniform(0, Z) <= z:
                    self.search.append(self.search[i])
            if len(self.search) > N:
                x = [(self.scape[self.minor(*co)], (i, co)) for i, co in enumerate(self.search)]
                x.sort()
                x = x[-N:]
                x = [y for _, y in x]
                x.sort()
                self.search = [y for _, y in x]
    def plot(self):
        xi, yi = self.grid
        im = plt.pcolor(xi, yi, self.scape, cmap=mplt.cm.jet)
        plt.xlim(0, self.params.X)
        plt.ylim(0, self.params.Y)
        x, y = [x for x, _ in self.search], [y for _, y in self.search]
        plt.scatter(x, y, c='black', marker='.')
    def run(self, N, scape=1, plot=1, animate=True):
        if not animate:
            n = math.ceil((N/plot)**.5)
            C = n
            R = 1
            while R*C <= N/plot:
                R += 1
            pi = 0
        for n in range(N):
            s.incrScape(scape)
            s.incrOpts()
            if n%plot == 0:
                if animate:
                    plt.subplot(1, 1, 1)
                else:
                    pi += 1
                    plt.subplot(R, C, pi)
                s.plot()
                plt.pause(0.001)

def xplot(fn):
    fn = np.vectorize(fn)
    lo, hi, dpu = 0, 10, 10
    N = (hi-lo)*dpu+1
    mesh = np.linspace(lo, hi, N)
    xi, yi = np.meshgrid(mesh, mesh)
    zi = np.fromfunction(lambda x, y: fn(x/dpu, y/dpu), (N, N))
    plt.pcolor(xi, yi, zi, cmap=mplt.cm.jet)
    plt.colorbar()
