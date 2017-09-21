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
import sympy

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
reef = ['intrashelf/intraplatform reef','reef','buildup or bioherm',
        'perireef or subreef','platform/shelf-margin reef','basin reef',
        'slope/ramp reef']
nonreef = ['coastal indet.','delta front','delta plain','deltaic indet.',
           'estuary/bay','foreshore','interdistributary bay','lagoonal',
           'lagoonal/restricted shallow subtidal','marginal marine indet.',
           'open shallow subtidal','fluvial-deltaic indet.','paralic indet.',
           'peritidal','prodelta','sand shoal','shallow subtidal indet.',
           'shoreface','transition zone/lower shoreface','basinal (carbonate)',
           'basinal (siliceous)','basinal (siliciclastic)',
           'deep-water indet.','deep subtidal indet.','deep subtidal ramp',
           'deep subtidal shelf','offshore','offshore indet.','offshore shelf',
           'slope','submarine fan','offshore ramp']

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

def location(r, latm=10, lngm=10):
    lat, lng = 'paleolat', 'paleolng'
    if lat not in r or lng not in r:
        return
    def f(s, degrees): return math.floor(float(r[s])//degrees+0.5)
    try:
        return f(lng, lngm), f(lat, latm)
    except:
        return

def getSplits(timeLevel=5):
    with open('resources/time.csv') as f:
        splits = [float(t['max_ma']) for t in csv.DictReader(f) if 
                  int(t['scale_level']) == timeLevel]
        splits.sort()
    return splits
def getFields(field):
    if field[-1] == '*':
        return [field[:-1]+'1', field[:-1]+'2']
    else:
        return [field]    

def fossilDistributions(data, field, h1_, h2_, 
                        trackLevel, trackRank,
                        timeLevel, timef, start, end,
                        lblFn,
                        dropSingle, correct,
                        macro):
    splits = getSplits(timeLevel)
    fields = getFields(field)
    # Intersection of values with category.
    h1desc, h2desc = set(h1_), set(h2_)
    def inter(vals, h):
        vals = set(x.strip().strip('"') for x in vals.split(','))
        return bool(vals & h)
    # Categories and labels
    h1, h2 = ([{} for _ in range(len(splits))] for _ in range(2))
    h1lbl, h2lbl = ([set() for _ in range(len(splits))] for _ in range(2))
    sp2i, sp2lbl = {}, {}
    # Add to category data structures
    def addTo(isH1, isH2, i, species):
        # add overlaps to both categories
        def f(h):
            if species not in h[i]: h[i][species] = 0
            h[i][species] += 1
        if isH1: f(h1)
        if isH2: f(h2)
        if species not in sp2i: sp2i[species] = set()
        sp2i[species].add(i)
    # Required fields
    must = [trackLevel, 'accepted_rank', 'max_ma', 'min_ma'] + fields
    for r in data:
        if any(x not in r for x in must): continue
        # Habitat
        isH1 = isH2 = False
        for f in fields:
            isH1 = isH1 or inter(r[f], h1desc)
            isH2 = isH2 or inter(r[f], h2desc)
        # The big mess o' time
        e, l = timef(float(r['max_ma'])), timef(float(r['min_ma']))
        assert l <= e
        if l > start or e < end: continue
        e = min(e, start)
        l = max(l, end)
        eb = bisect.bisect_left(splits, e)
        lb = bisect.bisect_right(splits, l)
        if lb > eb:
            # Right on the boundary, eh? Assume at the start of the later
            # interval, e.g., post-extinction.
            assert e == l
            lb = eb
        orig = e, l, eb, lb
        while correct and lb < eb:
            # Hmmm... one species over multiple intervals. If there's only
            # a relatively minor overlap, go with just the main interval.
            mid = splits[lb]
            assert l <= mid
            assert mid <= e
            a = (splits[lb-1] if lb > 0 else 0) - mid
            b = mid - (splits[lb+1] if lb+1 < len(splits) else end)
            x = (mid-l)/a
            y = (e-mid)/b
            # Compare against 1 my.
            if mid-l < 1 and x < y:
                l = mid
                lb += 1
            elif e-mid < 1 and y < x:
                e = mid
                eb -= 1
            else:
                # Big overlap, so go with it.
                break
        assert lb-1 <= orig[3]
        assert eb+1 >= orig[2]
        assert lb <= eb
        assert eb < len(h1)
        # Tracking level
        sp = r[trackLevel]
        accept = (trackLevel != 'accepted_name' or 
                  trackRank is None or 
                  r['accepted_rank'] == trackRank)
        # Label is, e.g., a spatial grid position
        lbl = lblFn(r)
        for i in range(lb, eb+1):
            if accept: addTo(isH1, isH2, i, sp)
            if isH1: h1lbl[i].add(lbl)
            if isH2: h2lbl[i].add(lbl)
        if sp not in sp2lbl: sp2lbl[sp] = set()
        if lbl is not None: sp2lbl[sp].add(lbl)
    # Remove singletons and add species throughout their ranges
    for species, intervals in sp2i.items():
        mini, maxi = min(intervals), max(intervals)
        if dropSingle and mini == maxi:
            # Meh... misses the case of a multi-interval single occurrence
            for h in (h1, h2):
                if species in h[mini]:
                    del h[mini][species]
        elif macro:
            # Propagate last environment occurrences
            def which(j):
                return [h for h in (h1, h2) if species in h[j]]
            last = which(maxi)
            for j in range(maxi-1, mini, -1):
                curr = which(j)
                if not curr:
                    for h in last:
                        if species not in h[j]: h[j][species] = 0
                        h[j][species] += 1
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
        x, y = h1lbl[i], h2lbl[i]
        c = x&y
        x, y = x-c, y-c
        x, c, y = len(x), len(c)/2, len(y)
        n, d = x+c, x+c+c+y
        if d == 0: 
            flbl.append(None)
        else:
            flbl.append(n/d)
    # Convert dict->set if macro
    if macro:
        dists = [tuple(set(x.keys()) for x in x) for x in dists]
    return times, dists, sp2lbl, flbl

def speciesInClades(data, trackLevel='accepted_name', groupLevel='order',
                    restrictLevel=None, restrict=None):
    g2s = {}
    for r in data:
        if trackLevel not in r or groupLevel not in r:
            continue
        if restrictLevel and restrictLevel not in r:
            continue
        if restrictLevel:
            rg = r[restrictLevel]
            if ((type(restrict) == set and rg not in restrict) or 
                (type(restrict) == str and rg != restrict)):
                continue
        nm = r[trackLevel]
        gp = r[groupLevel]
        if nm == '' or gp == '': continue
        if gp not in g2s: g2s[gp] = set()
        g2s[gp].add(nm)
    return g2s

def splitDists(dists, sp):
    def inter(x, s):
        y = {}
        for k in (set(x.keys()) & s):
            y[k] = x[k]
        return y
    def diff(x, s):
        y = {}
        for k in (set(x.keys()) - s):
            y[k] = x[k]
        return y
    rv = []
    for h1, h2 in dists:
        if type(h1) == set:
            rv.append([h1&sp, h1-sp, h2&sp, h2-sp])
        else:
            rv.append([inter(h1,sp), diff(h1,sp), inter(h2,sp), diff(h2,sp)])
    return rv

def lenDist(d):
    def f(x):
        if type(x) == set:
            return len(x)
        elif type(x) == dict:
            return sum(x for _, x in x.items())
        assert False
    return [f(x) for x in d]
def sdInter(x, y):
    assert type(x) == type(y)
    if type(x) == set:
        return x & y
    assert type(x) == dict
    rv = {}
    for sp in (set(x.keys()) & set(y.keys())):
        rv[sp] = x[sp]
    return rv
def sdDiff(x, y):
    assert type(x) == type(y)
    if type(x) == set:
        return x - y
    assert type(x) == dict
    rv = {}
    for sp in (set(x.keys()) - set(y.keys())):
        rv[sp] = x[sp]
    return rv

def plotAffinity(title, times, dists, sp2loc, color=None, ax=None, hf=None, xx=None, std=False, **kwargs):
    assert len(dists[0]) == 4
    assert len(dists) == len(times)
    if not ax: fig = plt.figure()
    def ld(d): return lenDist(d)
    def f(dl, op, d, ext=None, orig=None):
        ds = [ld([op(d[i][0], d[1-i][0]), 
                  d[i][1], 
                  op(d[i][2], d[1-i][2]), 
                  d[i][3]]) for i in range(2)]
        p0, p1 = (pAffinity(d, False, **kwargs) for d in ds)
        if ext is not None: ext.append(p0)
        if orig is not None: orig.append(p1)
        dl.append(pAffinityDiff(None, None, False, p1=p1, p0=p0, **kwargs))
    def flocs(d):
        return nlocs(d[0]|d[2])/nlocs(reduce(set.union, d))
    def plotHdis(ax, title, ps, diff, color=None, shift=0):
        nonlocal xx
        if color: color = color[:]
        hs = hdis(ps)
        if diff:
            x = [-(times[i]+times[i+1])/2 for i in range(len(times)-1)]
        elif shift:
            x = [-x for x in (times[1:] if shift == 1 else times[:-1])]
            if color: color = color[1:] if shift == 1 else color[:-1]
            assert len(x) == len(ps)
        else:
            x = [-x for x in times]
        # Prune down to interesting interval
        for i in (0, -1):
            while hs[i][2]-hs[i][0] > 0.6:
                hs.pop(i)
                x.pop(i)
                if color: color.pop(i)
                if hf: hf.pop(i)
                if not hs: return
        ax.set_title(title)
        if not diff and not xx:
            xx = x[-1]-5, x[0]+5
        ax.set_xlim(xx)
        if std: ax.set_ylim((-1, 1) if diff else (0, 1))
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
    dall, ddel, dana, deo, dext, dorig = ([] for _ in range(6))
    for i in range(len(dists)):
        d0 = dists[i]
        ld0 = ld(d0)
        if hf:
            dall.append(pAffinityPct(ld0[0], ld0[2], hf[i], False, **kwargs))
        else:
            dall.append(pAffinity(ld0, False, **kwargs))
        if i == 0: continue
        d1 = dists[i-1]
        ddel.append(pAffinityDiff(None, None, False, p1=dall[-2], p0=dall[-1], **kwargs))
        f(dana, sdInter, [d0, d1])
        f(deo, sdDiff, [d0, d1], dext, dorig)
    if ax:
        plotHdis(ax, title+' (affinity)', dall, False, color=color)
    else:
        plotHdis(fig.add_subplot(3, 2, 1), title+' (affinity)', dall, False, color=color)
        cc = color[:-1] if color else None
        plotHdis(fig.add_subplot(3, 2, 2), 'Change: all', ddel, True, color=cc)
        plotHdis(fig.add_subplot(3, 2, 4), 'Change: anagenetic', dana, True, color=cc)
        plotHdis(fig.add_subplot(3, 2, 3), 'Affinity: (pseudo-)extinction', dext, False, color=color, shift=1)
        plotHdis(fig.add_subplot(3, 2, 5), 'Affinity: origination', dorig, False, color=color, shift=-1)
        plotHdis(fig.add_subplot(3, 2, 6), 'Change: (pseudo-)extinction + origination', deo, True, color=cc)
        fig.tight_layout()

def plotAffinityByRange(title, times, dists, sp2loc, start=541, end=0, hf=None, **kwargs):
    def nlocs(sps):
        locs = set()
        for sp in sps:
            locs = locs | sp2loc[sp]
        return len(locs)
    r2h = {}
    for i, d in enumerate(dists):
        if times[i] > start or times[i] < end: continue
        r = nlocs(d[0]|d[2])
        ld = lenDist(d)
        if hf:
            p = pAffinityPct(ld[0], ld[2], hf[i], False, **kwargs)
        else:
            p = pAffinity(ld, False, **kwargs)
        hdi = hdis([p])[0]
        if r not in r2h: r2h[r] = []
        r2h[r].append((hdi, i))
    xy = []
    for r, hs in r2h.items():
        for j, (h, i) in enumerate(hs):
            xy.append((r+j/len(hs), h, times[i]))
    xy.sort()
    x, y, c = [list(x) for x in zip(*xy)]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title)
    ax.axhline(0.5, color='black', linestyle=':')
    ax.errorbar(x,
                [x[1] for x in y],
                np.array([[x[1]-x[0] for x in y],
                          [x[2]-x[1] for x in y]]),
                fmt='.')
    im = ax.scatter(x, [x[1] for x in y], c=c, cmap=mplt.cm.jet, 
                    marker='.', zorder=100)
    fig.colorbar(im)

def readData(file='../../Paleo/Simpson/Corals/pbdb_animalia_marine_090217.csv'):
    return DB(db=file).data

class EnvAffinity:
    def __init__(self, data, field, h1desc, h2desc,
                 trackLevel='accepted_name', trackRank=None,
                 timeLevel=5, timef=lambda x: x,
                 start=541, end=0,
                 dropSingle=False, correct=False,
                 macro=True,
                 latDegrees=5, lngDegrees=5):
        self.data = data
        self.field = field
        self.trackLevel = trackLevel
        self.macro = macro
        loc = lambda r: location(r, latDegrees, lngDegrees)
        self.times, self.dists, self.sp2loc, self.floc = fossilDistributions(
                data, field, h1desc, h2desc, trackLevel, trackRank, timeLevel, 
                timef, start, end, loc, dropSingle, correct, macro)
    def focusOn(self, groupLevel='phylum', threshold=9, **kwargs):
        self.g2s = speciesInClades(self.data, self.trackLevel, groupLevel, **kwargs)
        gs = []
        for g in self.g2s:
            sp = self.g2s[g]
            gdists = splitDists(self.dists, sp)
            affs = []
            for d in gdists:
                ld = lenDist(d)
                if ld[0] + ld[2] >= threshold:
                    affs.append(naiveAffinity(ld, False))
            if affs: 
                gs.append((math.floor(math.log(len(sp), 2)), max(affs)-min(affs), g))
        gs.sort()
        gs.reverse()
        for i in range(len(gs)):
            if gs[i][0] < threshold:
                gs = gs[:i]
                break
        self.interesting = [(i, g[-1], len(self.g2s[g[-1]]), g[1]) for i, g in enumerate(gs)]
        return self.interesting
    def gt(self, rng=0.5):
        return [x for x in self.interesting if x[-1] > rng]
    def group(self, idx):
        assert hasattr(self, 'interesting')
        if type(idx) == str:
            gp = idx
        else:
            assert 0 <= idx and idx < len(self.interesting)
            gp = self.interesting[idx][1]
        assert gp in self.g2s
        return gp
    def trim(self, dists, start=541, end=0):
        times, floc = self.times[:], self.floc[:]
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
        return times, dists, floc
    def plotByTime(self, idx, start=None, end=None, color=False, ax=None, **kwargs):
        gp = self.group(idx)
        sp = self.g2s[gp]
        times, dists, floc = self.trim(splitDists(self.dists, sp), start, end)
        h = None if self.macro else [0.5 for _ in times]
        plotAffinity(gp+': '+self.field, times, dists, self.sp2loc, floc if color else None, ax, hf=h, **kwargs)
    def plotByRange(self, idx, **kwargs):
        gp = self.group(idx)
        sp = self.g2s[gp]
        times, dists, floc = self.trim(splitDists(self.dists, sp))
        h = None if self.macro else [0.5 for _ in times]
        plotAffinityByRange(gp+': '+self.field, times, dists, self.sp2loc, hf=h, **kwargs)
    def plotGenera(self, taxon, level='family', ax=None):
        assert not self.macro
        g2s = speciesInClades(self.data, self.trackLevel, 'genus', restrictLevel=level, restrict=taxon)
        s2g = {}
        for g, sps in g2s.items():
            for sp in sps:
                s2g[sp] = g
        g2c = {g:[0,0] for g in g2s}
        for d in self.dists:
            for sp, g in s2g.items():
                for i in range(2):
                    if sp in d[i]:
                        g2c[g][i] += d[i][sp]
        ps = []
        for g, (t0, t1) in g2c.items():
            if not (t0 or t1): continue
            ps.append(pAffinityPct(t0, t1, 0.5, False))
        if not ps: return
        mesh = ps[0][0]
        x = np.zeros((len(mesh),))
        for _, p in ps:
            x += np.array(p)
        s = sum(x)/(len(mesh)+1)
        x = [e/s for e in x]
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title(taxon+' ('+level+'): genera-level affinities')
        ax.plot(mesh, x)
    def plotGeneraFor(self, taxa):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for level, taxon in taxa:
            self.plotGenera(taxon, level, ax)
    def plots(self, idx, **kwargs):
        self.plotByTime(idx, **kwargs)
        self.plotByRange(idx, **kwargs)

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
    
def bivalve(data, biDb=None, smart=True, file=None):
    if not biDb:
        biDb = DB(db='../../Paleo/Simpson/Corals/pbdb_bivalvia_091317.csv').data
    if file:
        file = open(file, 'w')
    hasComp = [r for r in biDb if r['composition'] != '']
    cal = set([r['accepted_name'] for r in hasComp if r['composition'] == 'aragonite'])
    ara = set([r['accepted_name'] for r in hasComp if r['composition'] != 'aragonite'])
    both = cal | ara
    ncarb, nclas = 0, 0
    c2c, a2c = {}, {}
    ca, cl = set(carbonate), set(clastic)
    for r in data:
        nm = r['accepted_name']
        lith = r['lithology1']+','+r['lithology2']
        lith = lith.split(',')
        lith = set([x.strip().strip('"') for x in lith])
        carb = bool(lith & ca)
        clas = bool(lith & cl)
        ncarb += carb
        nclas += clas
        if nm not in both: continue
        g = r['genus']
        if g == '': continue
        if nm in cal:
            if g not in c2c: c2c[g] = [0, 0]
            if carb: c2c[g][0] += 1
            if clas: c2c[g][1] += 1
        if nm in ara:
            if g not in a2c: a2c[g] = [0, 0]
            if carb: a2c[g][0] += 1
            if clas: a2c[g][1] += 1
    print(ncarb, nclas)
    for title, g2c in (('Aragonite', a2c), ('Calcite', c2c)):
        n1, nm, n2, t1, t2 = 0, 0, 0, 0, 0
        gs = []
        for g, cnts in g2c.items():
            d = [cnts[0], ncarb-cnts[0], cnts[1], nclas-cnts[1]]
            t1 += cnts[0]
            t2 += cnts[1]
            try:
                h = hdis([pAffinity(d, False, smart=smart)])[0]
            except:
                # floating point issues
                continue
            mixed = False
            if h[0] > 0.5: 
                n1 += 1
            elif h[2] < 0.5: 
                n2 += 1
            else:
                nm += 1
                mixed = True
            if not mixed:
                gs.append(((h[0] > 0.5, h[2] < 0.5, h[2]-h[1], g), 
                           '  %-32s: %.2f %.2f %.2f'%(g, h[0], h[1], h[2])))
        gs.sort()
        print('%s: carbonate: %d, mixed: %d, clastic: %d'%(title, n1, nm, n2), file=file)
        h = hdis([pAffinity([t1, ncarb-t1, t2, nclas-t2], False, smart=smart)])[0]
        print('  Affinity: %.2f %.2f %.2f'%(h[0], h[1], h[2]), file=file)
        for _, s in gs:
            print(s, file=file)
    if file: file.close()

def scenario3(data):
    lit = EnvAffinity(data, 'lithology*', carbonate, clastic)
    fig = plt.figure()
    lit.focusOn('class')
    ax = fig.add_subplot(3,1,1)
    lit.plotByTime('Gastropoda', color=False, ax=ax, xx=[-541,0])
    ax = fig.add_subplot(3,1,3)
    lit.plotByTime('Bivalvia', color=False, ax=ax, xx=[-541,0])
    lit.focusOn('phylum')
    ax = fig.add_subplot(3,1,2)
    lit.plotByTime('Brachiopoda', color=False, ax=ax, xx=[-541,0])
    fig.tight_layout()

def affinityIsMode():
    a, h, m, n = sympy.symbols("a h m n")
    p = (a*h)/(a*h + (1-a)*(1-h))
    x = p**m * (1-p)**n
    return sympy.solve(sympy.Eq(sympy.diff(x, a)))
def affinityIsMode3():
    a, b, h, g, m, n, l = sympy.symbols("a b h g m n l")
    p = (a*h)/(a*h + b*g + (1-a-b)*(1-h-g))
    q = (b*g)/(a*h + b*g + (1-a-b)*(1-h-g))
    x = p**m * q**n * (1-p-q)**l
    dxa = sympy.diff(x, a).simplify()
    dxb = sympy.diff(x, b).simplify()
    sola = m/h / (m/h + n/g + l/(1-h-g))
    solb = n/g / (m/h + n/g + l/(1-h-g))
    return [dx.subs([(a, sola), (b, solb)]).simplify() for dx in (dxa, dxb)]
def affinityIsNotMode():
    a, h, f, g, m, n, o, p = sympy.symbols("a h f g m n o p")
    z = (1-a)/a
    f0 = f/(h+z*(1-h))
    f1 = z*f0
    #f0, f1 = f, g
    x = (f0*h)**m * ((1-f0)*h)**n * (f1*(1-h))**o * ((1-f1)*(1-h))**p
    return sympy.solve(sympy.Eq(sympy.diff(x, a)))