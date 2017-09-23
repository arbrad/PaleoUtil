# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 15:48:10 2017

@author: Aaron
"""

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
from scipy.stats import binom

lncr_ = {}
def lncr(Rs):
    key = ','.join(str(x) for x in Rs)
    if key not in lncr_:
        N = sum(Rs)
        num = sum(math.log(x) for x in range(N, Rs[-1], -1))
        den = sum(sum(math.log(x) for x in range(1, x+1)) for x in Rs[:-1])
        lncr_[key] = num-den
    return lncr_[key]
def lmultinomial(Rs, Ps):
    '''Log of multinomial PDF.'''
    return lncr(Rs) + sum(Rs[i]*math.log(Ps[i]) for i in range(len(Ps)))
def multinomial(Rs, Ps):
    '''Multinomial PDF on len(Rs) [== len(Ps)] categories with Rs[i] items
    and probability Ps[i] for category i. Implicitly, N == sum(Rs).'''
    assert len(Rs) == len(Ps)
    for i in range(len(Rs)):
        if Ps[i] < 1e-5:
            if Rs[i]: return 0
            return multinomial(Rs[:i]+Rs[i+1:], Ps[:i]+Ps[i+1:])
    return math.exp(lmultinomial(Rs, Ps))

def hdi(x, pct=0.95):
    '''Compute highest density interval on single-mode distribution x. Returns
    (left, mode, right).'''
    norm = sum(x)
    i, s = max(enumerate(x), key=lambda x: x[1])
    l, r = i-1, i+1
    while s/norm < pct:
        if l < 0:
            s += x[r]
            r += 1
        elif r == len(x):
            s += x[l]
            l -= 1
        elif x[l] > x[r]:
            s += x[l]
            l -= 1
        else:
            s += x[r]
            r += 1
    return (max(0,l), i, min(len(x)-1,r))
def hdis(xs, pct=0.95):
    '''Compute HDIs over sequence of x,y-defined distributions.'''
    # compute hdis for each of cs in terms of mesh
    hdis = []
    for x, y in xs:
        l, m, r = hdi(y, pct)
        l, m, r = x[l], x[m], x[r]
        hdis.append((l, m, r))
    return hdis

def naiveAffinity(dist, absAff):
    try:
        if absAff:
            a, b = dist[0], dist[2]
        else:
            a, b = (dist[x+0]/(dist[x+0]+dist[x+1]) for x in (0,2))
        return a/(a+b) if a+b else -1
    except:
        return -1

def options(**kwargs):
    class Options:
        def __init__(self):
            self.grid = 100
            self.margin = False
            self.marginf = False
            self.marginh = False
            self.smart = False
            self.absAff = False
            self.hdip = 0.95
            self.simple = False
    o = Options()
    for k, v in kwargs.items():
        setattr(o, k, v)
    return o

def pAffinity(dist=[0,250,0,1000], plot=True, **kwargs):
    '''For dist=[H1-taxa-fossil, H1-other, H2-taxa-fossil, H2-other],
    compute the probability distribution of affinity for H1. When
    margin=False, the computation uses the sample frequency of taxa
    fossils, f, and H1 fossils, h; otherwise, it considers the joint
    a-f-h distribution and marginalizes out f and h. For small samples,
    the computed distributions differ w/ and w/o marginalizing, but
    as the sampe size increases (for the same ratios), they converge.
    grid controls the mesh granularity. Default affinity is "relative,"
    while absAff=True computes "absolute" affinity.'''
    opts = options(**kwargs)
    def g(a,f,h):
        # a = affinity for H1
        # f = frequency of taxa of interest
        # h = frequency of H1 samples compared to H2
        # Let x, y be frequencies defining the four categories:
        #  xh (1-x)h y(1-h) (1-y)(1-h)
        # i.e., x is the proportion of taxa fossils in H1, y is
        # the same in H2. We want
        #  a = x/(x+y)
        #  f = xh + y(1-h)
        # Solving for x and y in terms of a, f, and h yields the
        # following:
        if not opts.absAff:
            z = (1-a)/a
            x = f/(h+z*(1-h))
            y = z*x
        else:
            # instead, define a = xh/f
            x = a*f/h
            y = x*h/(1-h)*(1-a)/a
        pcat = [x*h, (1-x)*h, y*(1-h), (1-y)*(1-h)]
        if any(p < 0 or p > 1 for p in pcat):
            return 0
        return multinomial(dist, pcat)
    N = sum(dist)
    mesh = np.linspace(1/opts.grid,1-1/opts.grid,opts.grid-1)
    if N == 0: 
        return mesh, [1 for _ in mesh]
    marginf = marginh = opts.margin
    marginf = marginf or opts.marginf
    marginh = marginh or opts.marginh
    if opts.smart:
        # for almost all situations, marginalizing just one, especially f,
        # is sufficient
        marginf = marginf or any(x < 10 for x in dist)
        marginh = marginh or any(dist[i]+dist[i+1] < 2 for i in (0, 2))
    if opts.simple or not (marginf or marginh):
        n = dist[0]+dist[2]
        k = dist[0]
        h = sum(dist[:2])/N
        def g(a):
            z = (1-a)/a
            p = 1/(h+z*(1-h))*h
            return binom.pmf(k, n, p)
        x = [g(a) for a in mesh]
    else:
        meshf = mesh if marginf else [(dist[0]+dist[2])/N]
        meshh = mesh if marginh else [sum(dist[:2])/N]
        x = [sum(g(a,f,h) for f in meshf for h in meshh) for a in mesh]
    unit = 1/opts.grid
    s = sum(x)
    if not s:
        pass #print(dist)
    x = [e/s/unit for e in x]
    if plot: 
        _plotAffinity(mesh, x)
        print('%.2f'%naiveAffinity(dist, opts.absAff))
    return mesh, x

def _plotAffinity(mesh, x):
    p = plt.plot(mesh, x)
    l, _, r = hdi(x)
    plt.axvline(mesh[l], color=p[-1].get_color(), linestyle=':')
    plt.axvline(mesh[r], color=p[-1].get_color(), linestyle=':')
    plt.xlabel('Affinity')
    plt.yticks([])

def pAffinityPct(a, b, h, plot=True, priorA=None, priorH=None, **kwargs):
    opts = options(**kwargs)
    def g(a, h):
        z = (1-a)/a
        p = h/(h+z*(1-h))
        rv = binom.pmf(k, n, p)
        if priorA:
            rv = rv * priorA(a)
        return rv
    n = a+b
    k = a
    mesh = np.linspace(1/opts.grid,1-1/opts.grid,opts.grid-1)
    if priorH:
        x = [sum(g(a, h)*priorH(h) for h in mesh) for a in mesh]
    else:
        x = [g(a, h) for a in mesh]
    s = sum(x)
    x = [e/s*opts.grid for e in x]
    if plot: 
        _plotAffinity(mesh, x)
        na = a*(1-h)/(a*(1-h) + b*h)
        print('%.2f'%na)    
    return mesh, x

def pAffinityDiff(dist1, dist0, plot=True, p1=None, p0=None, **kwargs):
    '''For two distributions as in pAffinity, compute the 
    probability distribution of the difference of affinities, D1-D0.'''
    p0 = (pAffinity(dist0, plot=plot, **kwargs) if p0 is None else p0)[1]
    p1 = (pAffinity(dist1, plot=plot, **kwargs) if p1 is None else p1)[1]
    opts = options(**kwargs)
    x = [0 for _ in range(2*len(p0)+1)]
    mesh = np.linspace(1/opts.grid,1-1/opts.grid,opts.grid-1)
    for i, a0 in enumerate(mesh):
        for j, a1 in enumerate(mesh):
            d = a1-a0
            x[int((d+1)*opts.grid+0.5)] += p0[i]*p1[j]
    ls = np.linspace(-1+1/opts.grid,1-1/opts.grid,2*opts.grid-1)
    s = sum(x)
    x = [e/s*opts.grid for e in x]
    if plot:
        p = plt.plot(ls, x)
        l, _, r = hdi(x)
        plt.axvline(ls[l], color=p[-1].get_color(), linestyle=':')
        plt.axvline(ls[r], color=p[-1].get_color(), linestyle=':')
    return ls, x

def pAffinityChanges(dists, times=None, wiggle=0, **kwargs):
    '''Plot HDIs around modes of changes in affinity over given
    sequence of fossil distributions.'''
    assert (times is None or len(dists) == len(times))
    # compute all affinity distributions silently
    opts = options(**kwargs)
    ps = [pAffinity(d, False, **kwargs) for d in dists]
    phdis = hdis(ps, opts.hdip)
    # compute difference distributions silently using ps
    a, b = (1, 0) if times else (0, 1)
    cs = [pAffinityDiff(None, None, False, ps[i+b], ps[i+a], **kwargs) for
          i in range(len(ps)-1)]
    chdis = hdis(cs, opts.hdip)
    # plot: affinities/changes with HDIs over time
    for hs, offset in ((phdis, 0), (chdis, 0.5)):
        if times:
            b = offset+wiggle
            a = 1-b
            x = [-(a*times[i]+b*times[i+1]) for i in range(len(times)-1)]
            if len(hs) == len(times):
                x.append(-(a*times[-1] + b*(2*times[-1]-times[-2])))
        else:
            x = [x+offset+wiggle for x in range(len(hs))]
        plt.errorbar(x, 
                     [x[1] for x in hs],
                     np.array([[x[1]-x[0] for x in hs], 
                               [x[2]-x[1] for x in hs]]),
                     fmt='.')

def pAffinityAncDes(dists, anc=None, wiggle=0, **kwargs):
    '''Plot descendant vs ancestor affinities as modes with HDIs.'''
    # compute all affinity distributions silently
    opts = options(**kwargs)
    ps = [pAffinity(d, False, **kwargs) for d in dists]
    phdis = hdis(ps, opts.hdip)
    if anc:
        aps = [pAffinity(d, False, **kwargs) for d in anc]
        aphdis = hdis(aps, opts.hdip)
    else:
        aps = ps
        aphdis = phdis
    # plot: descendant affinity vs ancestor affinity
    def eb(hdi): 
        return np.array([[x[1]-x[0] for x in hdi],
                         [x[2]-x[1] for x in hdi]])
    plt.errorbar([x[1]+wiggle for x in aphdis[1:]], 
                 [x[1]+wiggle for x in phdis[:-1]], 
                 eb(aphdis[1:]), 
                 eb(phdis[:-1]), 
                 fmt=',')
    plt.plot([0,1], [0,1], c='black')

def simulateAffinity(N, start=[10,40,15,10], delta=10, **kwargs):
    dists = [np.array(start)]
    for _ in range(N-1):
        ch = np.array([random.randint(-delta, delta) for _ in dists[-1]])
        dists.append(dists[-1]+np.array(ch))
        for i in range(len(dists[-1])):
            if dists[-1][i] < 0:
                dists[-1][i] = 0
    for absAff in [False, True]:
        pAffinityChanges([x.tolist() for x in dists],
                          absAff=absAff, wiggle=0.1*absAff, **kwargs)
        na = [naiveAffinity(d, absAff) for d in dists]
        nca = [na[i+1]-na[i] for i in range(len(na)-1)]
        nca.insert(0, 0)
        print('in n-aff n-afc | h1t h1o h2t h2o')
        for i in range(len(na)):
            print('%2d % .2f % .2f | %s'%(i, na[i], nca[i],
                                         ' '.join('%3d'%x for x in dists[i])))

def testAffinity(a, h, N, runs, f=None, **kwargs):
    nok = 0
    for _ in range(runs):
        p = a*h / (a*h + (1-a)*(1-h))
        if not f:
            nt = binom.rvs(N, p)
            no = N-nt
            pa = pAffinityPct(nt, no, h, False, **kwargs)
        else:
            H0 = binom.rvs(N, h)
            H1 = N-H0
            z = (1-a)/a
            f0 = f/(h + z*(1-h))
            f1 = z*f0
            t0 = binom.rvs(H0, f0)
            t1 = binom.rvs(H1, f1)
            pa = pAffinity([t0, H0-t0, t1, H1-t1], False, **kwargs)
        hdi = hdis([pa], **kwargs)[0]            
        ok = hdi[0] <= a and a <= hdi[2]
        if runs <= 100:
            if not f:
                print('%2d %2d %.2f %.2f %.2f %s'%(nt, no, hdi[0], hdi[1], hdi[2], '*' if ok else ''))
            else:
                print('%2d %2d %2d %2d %.2f %.2f %.2f %s'%(t0, H0-t0, t1, H1-t1, hdi[0], hdi[1], hdi[2], '*' if ok else ''))
        nok += ok
    print(nok, runs)