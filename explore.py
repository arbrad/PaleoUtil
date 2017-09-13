# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 16:43:55 2017

@author: Aaron

A file for exploratory code. Fruitful explorations will be rewritten into 
clean code.
"""

from affinity import *
import bisect
import csv
from functools import reduce
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np
from pbdbDb import DB
import random

def plotFamilies(db, comps=[0,1,2], top=5):
    fam2sp = db.fieldSubsets('family', db.fieldSubset('order', 'Scleractinia'))
    fams = [(len(sp), fam) for fam, sp in fam2sp.items()]
    fams.sort()
    fams.reverse()
    sp2fam = {}
    for i, (sz, fam) in enumerate(fams):
        if i >= top: break
        print(fam, sz)
        for sp in fam2sp[fam]:
            sp2fam[sp] = fam
    f2i = {f:i for i, (_, f) in enumerate(fams)}
    def color(sp):
        if sp not in sp2fam: return 0
        fi = f2i[sp2fam[sp]]+1
        return fi #300*fi+np.mean(db.fossilRange(sp))
    db.plot(comps, color)

def plotTrait(db, field, trait, doHas, comps):
    scler = db.fieldSubset('order', 'Scleractinia')
    has = db.fieldSubset(field, trait) & scler
    hasnt = db.fieldSubset(field, trait, False) & scler
    both = has & hasnt
    which = has if doHas else hasnt
    db.plot(comps, lambda sp: sum(int(sp in x) for x in [scler,which,both]))

def traitAnalysis(field='diet', trait='photosymbiotic', has=False, db=None):
    db = db if db != None else DB('../../Paleo/Simpson/Corals/pbdb_animalia_marine_090217.csv')
    db.computePCA(['environment'])
    plotTrait(db, field, trait, has, [0,3])
    db.computePCA(['lithology*'])
    plotTrait(db, field, trait, has, [0,1])
    return db

# from Hopkins'14 Table 1
carbonate = ['limestone','carbonate','reef rocks','bafflestone','dolomite',
             'framestone','grainstone','lime mudstone','packstone','rudstone',
             'floatstone','wackestone']
clastic = ['shale','siliciclastic','volcaniclastic','claystone',
           'conglomerate','mudstone','phyllite','quartzite','sandstone',
           'siltstone','slate','schist']
shallow = ['coastal indet.','delta front','delta plain','deltaic indet.',
           'estuary/bay','foreshore','interdistributary bay','lagoonal',
           'lagoonal/restricted shallow subtidal','marginal marine indet.',
           'open shallow subtidal','fluvial-deltaic indet.','paralic indet.',
           'peritidal','prodelta','sand shoal','shallow subtidal indet.',
           'shoreface','transition zone/lower shoreface',
           'intrashelf/intraplatform reef','reef','buildup or bioherm',
           'perireef or subreef','platform/shelf-margin reef']
deep = ['basinal (carbonate)','basinal (siliceous)','basinal (siliciclastic)',
        'deep-water indet.','deep subtidal indet.','deep subtidal ramp',
        'deep subtidal shelf','offshore','offshore indet.','offshore shelf',
        'slope','submarine fan','offshore ramp','basin reef','slope/ramp reef']
def binaryAnalysis(db, comps, field, a, b):
    scler = db.fieldSubset('order', 'Scleractinia')
    asub = db.fieldSubset(field, a) 
    bsub = db.fieldSubset(field, b)
    both = asub & bsub
    db.plot(comps, 
            lambda sp: 4*int(sp in scler)+3-int(sp in bsub)-2*int(sp in asub))
    
def timeAnalysis(db, start=300):
    scler = db.fieldSubset('order', 'Scleractinia')
    db.computePCA(['environment'])
    db.plot([0,3], db.colorByTime(start=start), species=scler)
    db.plotMeanPCAOverTime([0,3], species=scler, start=start)
    db.computePCA(['lithology*'])
    db.plot([0,1], db.colorByTime(start=start), species=scler)
    db.plotMeanPCAOverTime([0,1], species=scler, start=start)

def location(r, degrees=10):
    lat, lng = 'lat', 'lng'
    if lat not in r or lng not in r:
        return
    def f(s): return math.floor(float(r[s])//degrees+0.5)
    try:
        return f(lng), f(lat)
    except:
        return

def fossilDistributions(data, field, h1_, h2_, 
                        trackLevel='accepted_name',
                        timeLevel=5,
                        lblFn=location,
                        start=541, end=0):
    # Time splits
    with open('resources/time.csv') as f:
        splits = [float(t['max_ma']) for t in csv.DictReader(f) if 
                  int(t['scale_level']) == timeLevel]
        splits.sort()
    # Expand field into multiple if necessary
    if field[-1] == '*':
        fields = [field[:-1]+'1', field[:-1]+'2']
    else:
        fields = [field]
    # Intersection of values with category.
    h1desc, h2desc = set(h1_), set(h2_)
    def inter(vals, h):
        vals = set(x.strip().strip('"') for x in vals.split(','))
        return bool(vals & h)
    # Categories and labels
    h1, h2, h1lbl, h2lbl = ([set() for _ in range(len(splits))] for _ in range(4))
    sp2i, sp2lbl = {}, {}
    # Add to category data structures
    def addTo(isH1, isH2, i, species):
        # add overlaps to both categories
        if isH1: h1[i].add(species)
        if isH2: h2[i].add(species)
        if species not in sp2i: sp2i[species] = set()
        sp2i[species].add(i)
    # Required fields
    must = [trackLevel, 'max_ma', 'min_ma'] + fields
    for r in data:
        if any(x not in r for x in must): continue
        isH1 = isH2 = False
        for f in fields:
            isH1 = isH1 or inter(r[f], h1desc)
            isH2 = isH2 or inter(r[f], h2desc)
        e, l = float(r['max_ma']), float(r['min_ma'])
        if l > start or e < end: continue
        e = min(e, start)
        l = max(l, end)
        eb = bisect.bisect_left(splits, e)
        lb = bisect.bisect_right(splits, l)
        if lb > eb:
            assert e == l
            lb = eb
        assert lb <= eb
        assert eb < len(h1)
        sp = r[trackLevel]
        lbl = lblFn(r)
        for i in range(lb, eb+1):
            addTo(isH1, isH2, i, sp)
            if isH1: h1lbl[i].add(lbl)
            if isH2: h2lbl[i].add(lbl)
        if sp not in sp2lbl: sp2lbl[sp] = set()
        if lbl is not None: sp2lbl[sp].add(lbl)
    # Remove singletons and add species throughout their ranges
    for species, intervals in sp2i.items():
        mini, maxi = min(intervals), max(intervals)
        if mini == maxi:
            # remove singleton
            for h in (h1, h2):
                if species in h[mini]:
                    h[mini].remove(species)
        else:
            # propagate last environment occurrences
            def which(j):
                return [h for h in (h1, h2) if species in h[j]]
            last = which(mini)
            for j in range(mini+1, maxi):
                curr = which(j)
                if not curr:
                    for h in last:
                        h[j].add(species)
                else:
                    last = curr
    # Prepare final list of fossil distributions
    dists = [list(x) for x in zip(h1, h2)]
    # Eliminate empty intervals
    while not any(dists[-1][i] for i in range(2)): 
        dists.pop()
        h1lbl.pop()
        h2lbl.pop()
    # Prepare times
    splits.insert(0, 0)
    times = [(splits[i]+splits[i+1])/2 for i in range(len(splits)-1)]
    times = times[:len(dists)]
    # Eliminate empty intervals
    while not any(dists[0][i] for i in range(2)):
        dists.pop(0)
        times.pop(0)
        h1lbl.pop(0)
        h2lbl.pop(0)
    # Compute frequency of h1 labels relative to all labels for each interval
    flbl = []
    for i in range(len(h1lbl)):
        x, y = len(h1lbl[i]), len(h2lbl[i])
        if x+y == 0: 
            flbl.append(None)
        else:
            flbl.append(x/(x+y))
    return times, dists, sp2lbl, flbl

def speciesInClades(data, trackLevel='accepted_name', groupLevel='order'):
    g2s = {}
    for r in data:
        if trackLevel not in r or groupLevel not in r:
            continue
        nm = r[trackLevel]
        gp = r[groupLevel]
        if gp == '': continue
        if gp not in g2s: g2s[gp] = set()
        g2s[gp].add(nm)
    return g2s

def splitDists(dists, sp):
    rv = []
    for h1, h2 in dists:
        rv.append([h1&sp, h1-sp, h2&sp, h2-sp])
    return rv

def lenDist(d):
    return [len(x) for x in d]

def plotAffinity(title, times, dists, sp2loc, color=None, **kwargs):
    assert len(dists[0]) == 4
    assert len(dists) == len(times)
    fig = plt.figure()
    def ld(d):
        return [len(x) for x in d]
    def f(dl, op, d):
        ds = [ld([op(d[i][0], d[1-i][0]), 
                  d[i][1], 
                  op(d[i][2], d[1-i][2]), 
                  d[i][3]]) for i in range(2)]
        dl.append(pAffinityDiff(ds[1], ds[0], False, **kwargs))
    def nlocs(sps):
        locs = set()
        for sp in sps:
            locs = locs | sp2loc[sp]
        return len(locs)
    def flocs(d):
        return nlocs(d[0]|d[2])/nlocs(reduce(set.union, d))
    xx = None
    def plotHdis(ax, title, ps, diff, color=None):
        nonlocal xx
        if color: color = color[:]
        hs = hdis(ps)
        if diff:
            x = [-(times[i]+times[i+1])/2 for i in range(len(times)-1)]
        else:
            x = [-x for x in times]
        for i in (0, -1):
            while hs[i][2]-hs[i][0] > 0.5:
                hs.pop(i)
                x.pop(i)
                if color: color.pop(i)
                if not hs: return
        ax.set_title(title)
        if not diff:
            xx = x[-1]-5, x[0]+5
        ax.set_xlim(xx)
        ax.axhline(0 if diff else 0.5, color='black', linestyle=':')
        ax.errorbar(x,
                    [x[1] for x in hs],
                    np.array([[x[1]-x[0] for x in hs],
                              [x[2]-x[1] for x in hs]]),
                    fmt='.', zorder=0)
        if color:
            im = ax.scatter(x, [x[1] for x in hs], c=color, cmap=mplt.cm.jet, 
                            marker='.', zorder=100)
            fig.colorbar(im)
    dall, ddel, dana, deo = [], [], [], []
    for i in range(len(dists)):
        d0 = dists[i]
        ld0 = ld(d0)
        dall.append(pAffinity(ld0, False, **kwargs))
        if i == 0: continue
        d1 = dists[i-1]
        ddel.append(pAffinityDiff(None, None, False, p1=dall[-2], p0=dall[-1], **kwargs))
        f(dana, set.intersection, [d0, d1])
        f(deo, set.difference, [d0, d1])
    plotHdis(fig.add_subplot(4, 1, 1), title+' (affinity)', dall, False, color=color)
    cc = color[:-1] if color else None
    plotHdis(fig.add_subplot(4, 1, 2), 'Change: all', ddel, True, color=cc)
    plotHdis(fig.add_subplot(4, 1, 3), 'Change: anagenetic', dana, True, color=cc)
    plotHdis(fig.add_subplot(4, 1, 4), 'Change: extinction + origination', deo, True, color=cc)

class EnvAffinity:
    def __init__(self, data, field, h1desc, h2desc,
                 trackLevel='accepted_name',
                 timeLevel=5):
        self.data = data
        self.field = field
        self.trackLevel = trackLevel
        self.times, self.dists, self.sp2loc, self.floc = fossilDistributions(
                data, field, h1desc, h2desc, trackLevel)
    def focusOn(self, groupLevel='order'):
        self.g2s = speciesInClades(self.data, self.trackLevel, groupLevel)
        gs = []
        for g in self.g2s:
            sp = self.g2s[g]
            gdists = splitDists(self.dists, sp)
            affs = []
            for d in gdists:
                ld = lenDist(d)
                if ld[0] + ld[2] >= 10:  # kind of arbitrary choice
                    affs.append(naiveAffinity(ld, False))
            if affs: 
                gs.append((math.floor(math.log(len(sp), 2)), max(affs)-min(affs), g))
        gs.sort()
        gs.reverse()
        for i in range(len(gs)):
            if gs[i][0] < 9:
                gs = gs[:i]
                break
        self.interesting = [(i, g[-1], len(self.g2s[g[-1]]), g[1]) for i, g in enumerate(gs)]
        return self.interesting
    def plot(self, idx, start=None, end=None, **kwargs):
        assert hasattr(self, 'interesting')
        assert 0 <= idx and idx < len(self.interesting)
        gp = self.interesting[idx][1]
        sp = self.g2s[gp]
        times, dists, floc = self.times[:], splitDists(self.dists, sp), self.floc[:]
        while ((start is not None and times[-1] > start) or 
                (not dists[-1][0] and not dists[-1][2])):
            times.pop()
            dists.pop()
            floc.pop()
        while ((end is not None and times[0] < end) or
                (not dists[0][0] and not dists[0][2])):
            times.pop(0)
            dists.pop(0)
            floc.pop(0)
        plotAffinity(gp+': '+self.field, times, dists, self.sp2loc, floc, **kwargs)

def affinityAnalysis(data, timeLevel=5, doColor=False, **kwargs):
    timesEnv, distsEnv, sp2locEnv, locEnv = fossilDistributions(
            data, 'environment', shallow, deep, timeLevel=timeLevel, start=250, **kwargs)
    timesLit, distsLit, sp2locLit, locLit = fossilDistributions(
            data, 'lithology*', carbonate, clastic, timeLevel=timeLevel, start=250, **kwargs)
    g2s = speciesInClades(data)
    distsEnv = splitDists(distsEnv, g2s)
    distsLit = splitDists(distsLit, g2s)
    plotAffinity('Environment: shallow/deep', timesEnv, distsEnv, sp2locEnv, locEnv if doColor else None)
    plotAffinity('Lithology: carbonate/clastic', timesLit, distsLit, sp2locLit, locLit if doColor else None)

def redQueenQM(orig=0.20, extn=0.04, slope=.001, decayO=.0025, incrE=.0025, grid=1):
    '''Model in Quental & Marshall 2013. Default parameters create Fig 3D.'''
    div = [1]
    first = True
    mesh = np.linspace(0, 1000, grid*1000)
    for time in mesh:
        o = orig - slope*div[-1]
        e = extn + slope*div[-1]
        if first and o <= e:
            first = False
            eq = time
        div.append(div[-1]*(1+(o-e)/grid))
        if div[-1] <= 1: break
        orig -= decayO/grid
        extn += incrE/grid
    mesh = mesh[:len(div)]
    plt.plot(mesh, div)
    plt.axvline(eq, linestyle=':')