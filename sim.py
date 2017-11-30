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

def bitstr(l):
    return ''.join('1' if x else '0' for x in l)

class BryParams:
    def __init__(self):
        self.feed = 0.1
        self.asex = 1
        self.sex = 1
        self.sexAsex = 0.5
        self.cost = 0.25
        self.copy = 0.01
        self.change = 0.1
        self.repop = 10
bryPm = BryParams()

class BrySpecies:
    def __init__(self, mutateFrom=None):
        if mutateFrom is None:
            self.topo = [[True]]
            self.sex = [True]
            self.lineage = None
        else:
            while True:
                self.topo = [x[:] for x in mutateFrom.topo]
                self.sex = mutateFrom.sex[:]
                if random.uniform(0, 1) < bryPm.copy:
                    # copy type
                    i = random.randint(0, len(self.sex)-1)
                    self.topo.append(self.topo[i][:])
                    for l in self.topo:
                        l.append(l[i])
                    self.sex.append(self.sex[i])
                elif random.uniform(0, 1) < bryPm.change:
                    # mutate attribute
                    N = len(self.sex)
                    i = random.randint(0, N*(N+1)-1)
                    if i < N:
                        self.sex[i] ^= True
                    else:
                        i -= N
                        self.topo[i%N][i//N] ^= True
                self.minimize()
                if any(self.sex):
                    break
            self.lineage = mutateFrom
        assert len(self.topo) == len(self.sex)
        assert len(self.topo) == len(self.topo[0])
    def minimize(self):
        reach = set()
        self.reach(0, reach)
        unreach = list(set(range(len(self.sex))) - reach)
        unreach.sort()
        unreach.reverse()
        for i in unreach:
            self.sex.pop(i)
            self.topo.pop(i)
            for l in self.topo:
                l.pop(i)
    def reach(self, i, reach):
        if i in reach: return
        reach.add(i)
        for j in range(len(self.topo[i])):
            if self.topo[i][j]:
                self.reach(j, reach)
    def asexFreq(self):
        if all(self.sex): return 0
        N = len(self.sex)
        CYCLE = 1
        freq = [[1/N for _ in self.sex]]
        freqPer = [0 for _ in self.sex]
        c = 0
        while True:
            c += 1
            if c > 1000:
                CYCLE += 1
                c = 0
                #print(CYCLE, self.sex, self.topo)
            for i in range(N):
                s = sum(self.topo[i])
                freqPer[i] = freq[-1][i]/s if s else 0
            nfreq = []
            for i in range(N):
                nfreq.append(sum(self.topo[j][i]*freqPer[j] for j in range(N)))
            s = sum(nfreq)
            if not s: return -1 
            nfreq = [x/s for x in nfreq]
            for j in range(len(freq)):
                if all(abs(nfreq[i]-freq[j][i]) < 1e-9 for i in range(N)):
                    s = 0
                    for k in range(j, len(freq)):
                        s += sum(x for i, x in enumerate(freq[k]) if not self.sex[i])
                    s /= len(freq)-j
                    return s
            if len(freq) == CYCLE:
                for j in range(len(freq)-1):
                    freq[j] = freq[j+1]
                freq.pop()
            freq.append(nfreq)
    def __str__(self):
        return (bitstr(self.sex) + ' ' + 
                ' '.join(bitstr(x) for x in self.topo) + ' ' +
                '%.2f'%self.asexFreq())

class BryAnimal:
    def __init__(self, colony, btype, x, y):
        self.colony = colony
        species = colony.species
        self.asex = species.topo[btype]
        self.sex = species.sex[btype]
        self.x = x
        self.y = y
    def act(self, energy, hasSpace):
        canAsex = energy >= bryPm.asex and hasSpace and any(self.asex)
        canSex = energy >= bryPm.sex and self.sex
        if canAsex and canSex:
            canSex = random.uniform(0, 1) <= bryPm.sexAsex
            canAsex = not canSex
        if canAsex:
            i = random.choice([j for j in range(len(self.asex)) if self.asex[j]])
            return bryPm.asex, i, None
        if canSex:
            return bryPm.sex, None, BrySpecies(self.colony.species)
        return 0, None, None
    def die(self, board):
        board.remove(self.x, self.y)

class BryColony:
    def __init__(self, species, x, y):
        self.species = species
        founder = BryAnimal(self, 0, x, y)
        self.animals = [founder]
        self.energy = 0
    def feed(self, unitEnergy):
        self.energy += sum(unitEnergy * (bryPm.cost if ani.sex else 1) for 
                           ani in self.animals)
    def act(self, board):
        eggs, newa = [], []
        anis = self.animals[:]
        random.shuffle(anis)
        for i, ani in enumerate(anis):
            z = board.random(ani.x, ani.y, 1)
            used, clone, egg = ani.act(self.energy / (len(self.animals)-i), z is not None)
            self.energy -= used
            if clone is not None:
                x, y = z
                clone = BryAnimal(self, clone, x, y)
                board.insert(x, y, clone)
                newa.append(clone)
            if egg is not None:
                eggs.append((egg, ani.x, ani.y, self.asexFreq()))
        self.animals.extend(newa)
        return eggs
    def size(self):
        return len(self.animals)
    def asexFreq(self):
        return 1 - sum(a.sex for a in self.animals)/len(self.animals)
    def die(self, board):
        for ani in self.animals:
            ani.die(board)
    def __str__(self):
        return '%3d %s'%(self.size(), str(self.species))

class BryBoard:
    def __init__(self, X, Y, uniform=False):
        self.cells = [[None for _ in range(X)] for _ in range(Y)]
        self.N = 0
        self.M = X*Y
        self.uniform = uniform
        if uniform: self.clear()
    def modxy(self, x, y):
        return x % len(self.cells[0]), y % len(self.cells)
    def cell(self, x, y):
        x, y = self.modxy(x, y)
        return self.cells[y][x]
    def insert(self, x, y, ani):
        x, y = self.modxy(x, y)
        assert self.cells[y][x] is None
        self.cells[y][x] = ani
        self.N += 1
        if self.uniform: 
            if (x, y) == self.empty[-1]:
                self.empty.pop()
            else:
                self.uniform = False
    def remove(self, x, y):
        x, y = self.modxy(x, y)
        assert self.cells[y][x] is not None
        self.cells[y][x] = None
        self.N -= 1
        if self.uniform: 
            self.uniform = False
    def random(self, x, y, r):
        if self.full(): return None
        if self.uniform:
            return self.empty[-1]
        rows = list(range(y-r, y+r+1))
        while len(rows):
            i = random.randint(0, len(rows)-1)
            b = rows[i] % len(self.cells)
            rows.pop(i)
            row = self.cells[b]
            space = [a % len(row) for a in range(x-r, x+r+1) if row[a % len(row)] is None]
            if not space: continue
            a = random.choice(space)
            return a, b
        return None
    def size(self):
        return self.M
    def full(self):
        return self.N == self.M
    def density(self):
        return self.N / self.M
    def clear(self):
        for y in self.cells:
            for i in range(len(y)):
                y[i] = None
        self.N = 0
        if self.uniform:
            self.empty = [(i, j) for 
                          i in range(len(self.cells)) for 
                          j in range(len(self.cells[0]))]
            random.shuffle(self.empty)
    def __str__(self):
        return '\n'.join(''.join(' ' if x is None else 'o' for x in y) for y in self.cells)
        
class BrySim:
    def __init__(self, X=51, Y=51, R=None):
        self.X = X
        self.Y = Y
        self.board = BryBoard(X, Y, uniform = R is None)
        if R is None: R = X // 2
        self.eggRad = R
        self.colonies = [BryColony(BrySpecies(), 0, 0)]
        x, y = self.board.random(0, 0, self.eggRad)
        self.board.insert(x, y, self.colonies[0].animals[0])
        self.thrAsexFreq = []
        self.empAsexFreq = []
        self.fillTime = []
        self.meanAF = []
        self.unwMeanAF = []
        self.meanCS = []
        self.stdCS = []
        self.meanNBT = []
    def run(self, steps=1000):
        fillTime = 0
        for st in range(steps):
            fillTime += 1
            for c in self.colonies:
                c.feed(bryPm.feed)
            eggs = []
            for c in self.colonies:
                eggs.extend(c.act(self.board))
            random.shuffle(eggs)
            eggs = eggs[0:min(len(eggs), bryPm.repop)]
            cleared = False
            while True:
                unique = set()
                for egg, x, y, _ in eggs:
                    z = self.board.random(x, y, self.eggRad)
                    if z is None: break
                    a, b = z
                    colony = BryColony(egg, a, b)
                    self.board.insert(a, b, colony.animals[0])
                    self.colonies.append(colony)
                    if cleared: unique.add(str(colony.species))
                if cleared: print('\n'.join(list(unique)))
                if self.board.full() and eggs:
                    cs = [c.size() for c in self.colonies]
                    nc = len(cs)
                    mcs = np.mean(cs)
                    sdcs = np.std(cs)
                    afs = [c.asexFreq() for c in self.colonies]
                    maf = sum(cs[i]*afs[i] for i in range(nc))/sum(cs)
                    unwmaf = sum(afs)/nc
                    nbt = [len(c.species.sex) for c in self.colonies]
                    mnbt = sum(nbt)/nc
                    print('Repopulating with %d (%d %d %.0f %.2f %.2f)'%(
                            len(eggs), nc, max(cs), mcs, sdcs, maf))
                    self.board.clear()
                    self.colonies.clear()
                    cleared = True
                else:
                    break
            if cleared:
                assert eggs
                af = [(egg.asexFreq(), caf) for egg, _, _, caf in eggs]
                af = [x for x in af if x[0] != -1]
                if af: 
                    self.thrAsexFreq.append(sum(x for x, _ in af)/len(af))
                    self.empAsexFreq.append(sum(x for _, x in af)/len(af))
                    self.fillTime.append(fillTime)
                    self.meanAF.append(maf)
                    self.unwMeanAF.append(unwmaf)
                    self.meanCS.append(mcs)
                    self.stdCS.append(sdcs)
                    self.meanNBT.append(mnbt)
                fillTime = 0
    def plot(self):
        fig = plt.figure()
        ax = plt.subplot(2, 2, 1)
        ax.plot(self.thrAsexFreq)
        ax.plot(self.empAsexFreq)
        ax.plot(self.meanAF)
        ax.plot(self.unwMeanAF)
        ax.set_title('Limit/empirical egg & wght/unwght mean colony asex freq')
        ax = plt.subplot(2, 2, 2)
        ax.plot(self.fillTime)
        ax.set_title('Sim cycles to fill')
        ax = plt.subplot(2, 2, 4)
        ax.plot(self.meanCS)
        ax.plot(self.stdCS)
        ax.set_title('Mean/std colony size')
        ax = plt.subplot(2, 2, 3)
        ax.plot(self.meanNBT)
        ax.set_title('Mean colony num body types')
    def size2colonies(self):
        sz2cs = {}
        for c in self.colonies:
            if c.size() not in sz2cs: sz2cs[c.size()] = []
            sz2cs[c.size()].append(c)
        css = list(sz2cs.items())
        css.sort()
        css.reverse()
        return css
    def analyze(self):
        css = self.size2colonies()
        big = css[0][0]
        for sz, cs in css:
            if sz < big/2:
                print(sz, len(cs))
            else:
                for c in cs:
                    print(str(c))
                print()
