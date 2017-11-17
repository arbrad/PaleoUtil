# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 20:40:03 2017

@author: Aaron
"""

import csv
import itertools
import math
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats
from scipy.interpolate import Rbf
import sympy

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

class RQ:
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

class NoLineage:
    def __init__(self, N, K, M, G, maxS):
        self.maxS = maxS
        self.macroRate = 0.01
        g = N//G
        nodes = []
        for i in range(G):
            nodes.append(np.arange(g*i, N if i+1 == g else g*(i+1)))
        self.neighbors = [np.random.choice(nodes[min(i//g, len(nodes)-1)], K, replace=False) for i in range(N)]
        for i in range(G-1):
            self.neighbors[g*(i+1)-1][-1] = g*(i+1)
        self.interact = np.random.normal(scale=0.05, size=(M,M))
        for i in range(self.interact.shape[0]):
            self.interact[i,i] = 0
        self.ma2mi = np.random.normal(scale=0.05, size=(N,M))
        self.mi2ma = np.random.normal(scale=0.05, size=(N,M))
        self.dist = np.zeros((N, M))
        for i, n in enumerate(np.random.choice(np.arange(N), M, replace=False)):
            self.dist[n,i] = 1
        self.size = [1 for _ in range(N)]
        self.colors = 'bgrcmyk'
    def incrReproduce(self):
        N = self.dist.shape[0]
        M = self.dist.shape[1]
        gain = np.zeros(self.dist.shape)
        for n in range(N):
            dist = np.zeros(M)
            for m in range(M):
                rate = self.ma2mi[n,m]
                for l in range(M):
                    rate += self.dist[n,l] * self.interact[m,l]
                # density is bad
                rate -= 0.05/(1+math.exp(-10*(self.dist[n,m]-0.5)))
                dist[m] = (1+rate) * self.dist[n,m]
            self.dist[n,:] += dist
            np.copyto(dist, self.dist[n,:])
            s = sum(self.dist[n,:])
            if s > 1: 
                self.dist[n,:] = self.dist[n,:]/s
                dist -= self.dist[n,:]
                dist = dist*self.size[n]/len(self.neighbors[n])
                for k in self.neighbors[n]:
                    gain[k,:] += dist
            self.size[n] *= (1+self.macroRate*self.dist[n,:].dot(self.mi2ma[n,:].transpose()))
        self.dist += gain/len(self.neighbors[0])
        for n in range(self.dist.shape[0]):
            s = sum(self.dist[n,:])
            if s > 1: self.dist[n,:] = self.dist[n,:]/s
        if len(self.neighbors) < self.maxS:
            for n in range(N):
                if self.size[n] > 2:
                    self.neighbors.append(np.copy(self.neighbors[n]))
                    self.neighbors[n][random.randint(0,len(self.neighbors[n])-1)] = len(self.neighbors)-1
                    self.ma2mi = np.append(self.ma2mi, [np.copy(self.ma2mi[n,:])], axis=0)
                    self.mi2ma = np.append(self.mi2ma, [np.copy(self.mi2ma[n,:])], axis=0)
                    self.dist = np.append(self.dist, [self.dist[n,:]/2+sum(self.dist)/(2*sum(sum(self.dist)))], axis=0)
                    self.size[n] = 1
                    self.size.append(1)
                    print(n, len(self.neighbors)-1)
    def mutate(self):
        M = self.interact.shape[0]
        x = self.interact
        self.interact = np.zeros((M+1,M+1))
        self.interact[0:M,0:M] = x
        self.interact[M,:] = np.random.normal(scale=0.05, size=M+1)
        self.interact[:,M] = np.random.normal(scale=0.05, size=M+1)
        self.interact[M,M] = 0
        x = self.dist
        self.dist = np.zeros((x.shape[0],M+1))
        self.dist[:,0:M] = x
        n = random.randint(0,self.dist.shape[0]-1)
        self.dist[n,M] = 0.1
        self.dist[n,:] /= sum(self.dist[n,:])
    def plot(self):
        x = np.arange(len(self.neighbors))
        plt.subplot(2,1,1)
        bot = np.zeros(len(self.neighbors))
        for k in range(self.dist.shape[1]):
            h = self.dist[:,k].transpose()
            plt.bar(x, h, bottom=bot, color=self.colors[k%len(self.colors)])
            bot += h
        plt.subplot(2,1,2)
        plt.bar(x, self.size, color='blue')
    def run(self, N, plot=1):
        for n in range(N):
            self.incrReproduce()
            if n%plot == 0: 
                self.plot()
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

def fixpointDist(N, P, A=1, prec=1e-12):
    assert P+A <= N
    Q = N - A - P
    a = 1
    b = 1 if P else 0
    c = 1 if Q else 0
    n = a+b+c
    a, b, c = a/n, b/n, c/n
    while True:
        an = A/N * a + b
        bn = P/N * a
        cn = Q/N * a
        n = an+bn+cn
        an, bn, cn = an/n, bn/n, cn/n
        if abs(an-a) < prec and abs(bn-b) < prec and abs(cn-c) < prec:
            return an, bn, cn
        else:
            a, b, c = an, bn, cn
def solveFixpointDist(N, P, A=1):
    a, b, c = sympy.symbols('a b c')
    Q = N-A-P
    c0 = a+b+c-1
    c1 = A/N*a + b - a*(a+b)
    c2 = P/N*a - b*(a+b)
    c3 = Q/N*a - c*(a+b)
    x = sympy.solve((c0, c1, c2, c3), a, b, c)
    return [y for y in x if y[2] < 1 and all(z >= 0 for z in y)][0]
def asexFreq(N, P, A=1):
    Q = N-A-P
    a, b, c = fixpointDist(N, P, A)
    return (1-a/A, 1-b/P if P > 0 else None, 1-c/Q if Q > 0 else None)
def plotBodyDists(Nmax, maxA=1, data=False, onlyAsex=False):
    x, y0, y1, y2, y3, col = [], [], [], [], [], []
    title = ['Autozooid', 'Asexual non-auto', 'Non-asexual', 'Non-sexual']
    for N in range(2, Nmax+1):
        for A in range(1, min(maxA,N)+1):
            for P in range(N-A+1):
                Q = N-A-P
                a, b, c = fixpointDist(N, P, A)
                assert abs(a+b+c-1) < 1e-6
                for i, (d, D) in enumerate([(a,A), (b,P), (c,Q)]):
                    if not D: continue
                    asex = 1 - d/D
                    var = 2*int(A>1)+int(i==0)
                    x.append(N + (var-3/2)/10)
                    y0.append(a)     # auto
                    y1.append(b)     # asexual non-auto, producing only auto
                    y2.append(c)     # non-asexual non-auto
                    y3.append(asex)  # non-ovicell
                    col.append('bgrc'[var])
                    # blue: A = 1, ovi != auto
                    # green: A = 1, ovi is auto
                    # red: A > 1, ovi != auto
                    # cyan: A > 1, ovi is auto
    ys, dim = ([y3], 1) if onlyAsex else ([y0, y1, y2, y3], 2)
    # plot in reverse order to emphasize simplest
    xrev = x[:]
    xrev.reverse()
    col.reverse()
    for i, y in enumerate(ys):
        ax = plt.subplot(dim, dim, i+1)
        yrev = y[:]
        yrev.reverse()
        ax.scatter(xrev, yrev, c=col)
        plt.title(title[i])
        plt.ylim((0, 1))
        for j in range(len(y)):
            ax.annotate(j, (x[j], y[j]))
    if data:
        with open('data09.16.2010.txt') as f:
            raw = [(int(x['polymorph.types']), float(x['nrr'])) for 
                   x in csv.DictReader(f, delimiter='\t') if 
                   x['Phylum'] == 'Bryozoa' and x['nrr'] != 'NA' and int(x['polymorph.types']) > 1]
        ax.scatter([x+random.normalvariate(0, 0.1) for x, _ in raw], [x for _, x in raw], c='black', marker='.')
