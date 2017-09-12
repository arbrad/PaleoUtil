# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 16:43:55 2017

@author: Aaron

A file for exploratory code. Fruitful explorations will be rewritten into 
clean code.
"""

from functools import reduce
from pbdbDb import DB
from affinity import pAffinityDiff, hdis
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

def fossilDistributions(data, field, h1_, h2_, 
                        groupLevel='order', group='Scleractinia',
                        trackLevel='accepted_name',
                        start=250, end=0, timeLevel=5,
                        lblFn=lambda x: 0):
    with open('resources/time.csv') as f:
        splits = [float(t['max_ma']) for t in csv.DictReader(f) if 
                  int(t['scale_level']) == timeLevel]
        splits.sort()
    h1, h2 = set(h1_), set(h2_)
    # Intersection of values with category.
    def inter(vals, h):
        vals = set(x.strip().strip('"') for x in vals.split(','))
        return bool(vals & h)
    # Categories: h1/h2 taxon/other
    h1t, h1o, h2t, h2o, h1lbl, h2lbl = ([set() for _ in range(len(splits)+1)] for _ in range(6))
    h_ = [h1o, h2o, h1t, h2t]
    sp2i, sp2lbl = {}, {}
    # Add to category data structures
    def addTo(isScler, isH1, isH2, i, species):
        # add overlaps to both categories
        if isH1: h_[2*isScler][i].add(species)
        if isH2: h_[2*isScler+1][i].add(species)
        if species not in sp2i: sp2i[species] = set()
        sp2i[species].add(i)
    # Process fields
    if field[-1] == '*':
        fields = [field[:-1]+'1', field[:-1]+'2']
    else:
        fields = [field]
    # Required fields
    must = [trackLevel, groupLevel, 'max_ma', 'min_ma'] + fields
    for r in data:
        if any(x not in r for x in must): continue
        isGroup = r[groupLevel] == group
        isH1 = isH2 = False
        for f in fields:
            isH1 = isH1 or inter(r[f], h1)
            isH2 = isH2 or inter(r[f], h2)
        e, l = float(r['max_ma']), float(r['min_ma'])
        if l > start: continue
        if e < end: continue
        eb = bisect.bisect_left(splits, e)
        lb = bisect.bisect_right(splits, l)
        assert lb <= eb
        assert eb < len(h1t)
        sp = r[trackLevel]
        lbl = lblFn(r)
        for i in range(lb, eb+1):
            addTo(isGroup, isH1, isH2, i, sp)
            if isH1: h1lbl[i].add(lbl)
            if isH2: h2lbl[i].add(lbl)
        if sp not in sp2lbl: sp2lbl[sp] = set()
        if lbl is not None: sp2lbl[sp].add(lbl)
    # Remove singletons and add species throughout their ranges
    for species, intervals in sp2i.items():
        mini, maxi = min(intervals), max(intervals)
        for h in h_:
            if species in h[mini]:
                if mini == maxi:
                    h[mini].remove(species)
                else:
                    for j in range(mini+1, maxi):
                        h[j].add(species)
    # Prepare final list of fossil distributions
    dists = [list(x) for x in zip(h1t, h1o, h2t, h2o)]
    while not (dists[-1][0] or dists[-1][2]): 
        dists.pop()
        h1lbl.pop()
        h2lbl.pop()
    # Prepare times
    splits.insert(0, 0)
    times = [(splits[i]+splits[i+1])/2 for i in range(len(splits)-1)]
    times = times[0:len(dists)]
    while not (dists[0][0] or dists[0][2]):
        dists.pop(0)
        times.pop(0)
        h1lbl.pop(0)
        h2lbl.pop(0)
    return times, dists, sp2lbl, [len(h1lbl[i])/(len(h1lbl[i])+len(h2lbl[i])) for i in range(len(h1lbl))]
def affinityAnalysis(data, timeLevel=5, ret=False, degrees=10, doColor=False, **kwargs):
    def location(r):
        lat, lng = 'lat', 'lng'
        if lat not in r or lng not in r:
            return
        def f(s): return math.floor(float(r[s])//degrees+0.5)
        try:
            return f(lng), f(lat)
        except:
            return
    timesEnv, distsEnv, sp2locEnv, locEnv = fossilDistributions(
            data, 'environment', shallow, deep, timeLevel=timeLevel, lblFn=location, **kwargs)
    timesLit, distsLit, sp2locLit, locLit = fossilDistributions(
            data, 'lithology*', carbonate, clastic, timeLevel=timeLevel, lblFn=location, **kwargs)
    def aa(title, dists, times, sp2loc, loc):
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
            return nlocs(d[0]|d[2]) #/nlocs(reduce(set.union, d))
        def plotHdis(ax, title, ps, diff, color=None):
            hs = hdis(ps)
            if diff:
                x = [-(times[i]+times[i+1])/2 for i in range(len(times)-1)]
            else:
                x = [-x for x in times]
            ax.set_title(title)
            ax.set_xlim([-times[-1]-5, -times[0]+5])
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
        dall, ddel, dana, deo, color = [], [], [], [], []
        for i in range(len(dists)):
            d0 = dists[i]
            ld0 = ld(d0)
            dall.append(pAffinity(ld0, False, **kwargs))
            if doColor: color.append(loc[i])
            if i == 0: continue
            d1 = dists[i-1]
            ddel.append(pAffinityDiff(None, None, False, p1=dall[-2], p0=dall[-1], **kwargs))
            f(dana, set.intersection, [d0, d1])
            f(deo, set.difference, [d0, d1])
        plotHdis(fig.add_subplot(4, 1, 1), title+' (affinity)', dall, False, color=color)
        plotHdis(fig.add_subplot(4, 1, 2), 'Change: all', ddel, True, color=color[:-1])
        plotHdis(fig.add_subplot(4, 1, 3), 'Change: anagenetic', dana, True, color=color[:-1])
        plotHdis(fig.add_subplot(4, 1, 4), 'Change: extinction + origination', deo, True, color=color[:-1])
    aa('Environment: shallow/deep', distsEnv, timesEnv, sp2locEnv, locEnv)
    aa('Lithology: carbonate/clastic', distsLit, timesLit, sp2locLit, locLit)
    if ret: return timesEnv, distsEnv, timesLit, distsLit

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