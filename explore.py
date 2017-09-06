# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 16:43:55 2017

@author: Aaron

A file for exploratory code. Fruitful explorations will be rewritten into 
clean code.
"""

from functools import reduce
import operator as op
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

def fossilDistributions(data, field, h1_, h2_, start=260, timeLevel=5):
    with open('resources/time.csv') as f:
        splits = [float(t['max_ma']) for t in csv.DictReader(f) if 
                  int(t['scale_level']) == timeLevel]
        splits.sort()
    h1, h2 = set(h1_), set(h2_)
    def inter(vals, h):
        vals = set(x.strip().strip('"') for x in vals.split(','))
        return bool(vals & h)
    h1t = [0 for _ in range(len(splits))]
    h1o, h2t, h2o = (h1t[:] for _ in range(3))
    h_ = [h1o, h2o, h1t, h2t]
    def addTo(isScler, isH1, isH2, i):
        h_[2*isScler+isH2][i] += 1
        # add overlaps to both categories
        if isH1: h_[2*isScler][i] += 1
    if field[-1] == '*':
        fields = [field[:-1]+'1', field[:-1]+'2']
    else:
        fields = [field]
    must = ['order', 'max_ma', 'min_ma'] + fields
    for r in data:
        if any(x not in r for x in must): continue
        isScler = r['order'] == 'Scleractinia'
        isH1 = isH2 = False
        for f in fields:
            isH1 = isH1 or inter(r[f], h1)
            isH2 = isH2 or inter(r[f], h2)
        e, l = float(r['max_ma']), float(r['min_ma'])
        if e > start: continue
        eb = bisect.bisect_left(splits, e)
        lb = bisect.bisect_right(splits, l)
        assert lb <= eb and eb < len(h1t)
        for i in range(lb, eb+1):
            addTo(isScler, isH1, isH2, i)
    dists = list(zip(h1t, h1o, h2t, h2o))
    while not (dists[-1][0] or dists[-1][2]): dists.pop()
    splits.insert(0, 0)
    times = [(splits[i]+splits[i+1])/2 for i in range(len(splits)-1)]
    times = times[0:len(dists)]
    return times, dists
def affinityAnalysis(data, smart=False, grid=100, hdip=0.95, timeLevel=5, 
                     ret=False, absAff=False):
    timesEnv, distsEnv = fossilDistributions(data, 'environment', shallow, deep, 
                                             timeLevel=timeLevel)
    timesLit, distsLit = fossilDistributions(data, 'lithology*', carbonate, clastic,
                                             timeLevel=timeLevel)
    plt.figure()
    plt.title('Environment: shallow/deep')
    pAffinityChanges(distsEnv, times=timesEnv, smart=smart, grid=grid, hdip=hdip,
                     absAff=absAff)
    plt.figure()
    plt.title('Lithology: carbonate/clastic')
    pAffinityChanges(distsLit, times=timesLit, smart=smart, grid=grid, hdip=hdip,
                     absAff=absAff)
    if ret: return timesEnv, distsEnv, timesLit, distsLit

lncr_ = {}
def lncr(N, Rs):
    key = str(N) + ':' + ','.join(str(x) for x in Rs)
    if key not in lncr_:
        num = sum(math.log(x) for x in range(N, Rs[-1], -1))
        den = sum(sum(math.log(x) for x in range(1, x+1)) for x in Rs[:-1])
        lncr_[key] = num-den
    return lncr_[key]
def lmultinomial(N, Rs, Ps):
    return lncr(N, Rs) + sum(Rs[i]*math.log(Ps[i]) for i in range(len(Ps)))
def multinomial(N, Rs, Ps):
    assert sum(Rs) == N
    assert len(Rs) == len(Ps)
    for i in range(len(Rs)):
        if Ps[i] < 1e-5:
            if Rs[i]: return 0
            return multinomial(N, Rs[:i]+Rs[i+1:], Ps[:i]+Ps[i+1:])
    return math.exp(lmultinomial(N, Rs, Ps))

def hdi(x, pct=0.95):
    '''Compute highest density interval on distribution x. Returns
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
    for x in xs:
        l, m, r = hdi(x[1], pct)
        l, m, r = x[0][l], x[0][m], x[0][r]
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

def pAffinity(dist=[0,250,0,1000], grid=100, margin=False, plot=True, 
              smart=False, absAff=False):
    '''For dist=[H1-taxa-fossil, H1-other, H2-taxa-fossil, H2-other],
    compute the probability distribution of affinity for H1. When
    margin=False, the computation uses the sample frequency of taxa
    fossils, f, and H1 fossils, h; otherwise, it considers the joint
    a-f-h distribution and marginalizes out f and h. For small samples,
    the computed distributions differ w/ and w/o marginalizing, but
    as the sampe size increases (for the same ratios), they converge.
    grid controls the mesh granularity.'''
    N = sum(dist)
    if N == 0: return
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
        if not absAff:
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
        return multinomial(N, dist, pcat)
    mesh = np.linspace(1/grid,1-1/grid,grid-1)
    if smart:
        margin = any(x < 10 for x in dist)
    if margin:
        # Marginalize, which can matter for small sample sizes
        x = [sum(g(a,f,h) for f in mesh for h in mesh) for a in mesh]
    else:
        # Use sample means for f, h
        f = (dist[0]+dist[2])/N
        h = sum(dist[:2])/N
        x = [g(a,f,h) for a in mesh]
    unit = 1/grid
    s = sum(x)
    if not s:
        print(dist)
    x = [e/s/unit for e in x]
    if plot:
        p = plt.plot(mesh, x)
        l, _, r = hdi(x)
        plt.axvline(mesh[l], color=p[-1].get_color(), linestyle=':')
        plt.axvline(mesh[r], color=p[-1].get_color(), linestyle=':')
        half = math.floor(0.5*(len(x)-1))
        print('%.2f %.2f'%(unit*sum(x[:1+half]),
                           naiveAffinity(dist, absAff)))
    return mesh, x

def pAffinityChange(dist0, dist1, grid=100, margin=False,
                    plot=True, p0=None, p1=None, absAff=False):
    '''For two distributions as in pAffinity, compute the 
    probability distribution of the change in affinity for H1.'''
    p0 = pAffinity(dist0, grid, margin, plot=True, absAff=absAff)[1] if p0 is None else p0
    p1 = pAffinity(dist1, grid, margin, plot=True, absAff=absAff)[1] if p1 is None else p1
    x = [0 for _ in range(2*len(p0)+1)]
    mesh = np.linspace(1/grid,1-1/grid,grid-1)
    for i, a0 in enumerate(mesh):
        for j, a1 in enumerate(mesh):
            d = a1-a0
            x[int((d+1)*grid+0.5)] += p0[i]*p1[j]/grid
    ls = np.linspace(-1+1/grid,1-1/grid,2*grid-1)
    if plot:
        p = plt.plot(ls, x)
        l, _, r = hdi(x)
        plt.axvline(ls[l], color=p[-1].get_color(), linestyle=':')
        plt.axvline(ls[r], color=p[-1].get_color(), linestyle=':')
    return ls, x

def pAffinityChanges(dists, hdip=0.95, grid=100, margin=False, times=None, 
                     smart=False, absAff=False):
    '''Plot HDIs around modes of changes in affinity over given
    sequence of fossil distributions.'''
    # compute all affinity distributions silently
    ps = [pAffinity(d, grid, margin, False, smart=smart, absAff=absAff) for d in dists]
    phdis = hdis(ps, hdip)
    # compute difference distributions silently using ps
    a, b = (1, 0) if times else (0, 1)
    cs = [pAffinityChange(None, None, grid, margin, False, 
                          ps[i+a][1], ps[i+b][1]) for
          i in range(len(ps)-1)]
    chdis = hdis(cs, hdip)
    for hs, offset in ((phdis, 0), (chdis, 0.5)):
        if times:
            if offset:
                x = [-(times[i]+times[i+1])/2 for i in range(len(times)-1)]
            else:
                x = [-x for x in times]
        else:
            x = [x+offset for x in range(len(hs))]
        plt.errorbar(x, 
                     [x[1] for x in hs],
                     np.array([[x[1]-x[0] for x in hs], 
                               [x[2]-x[1] for x in hs]]),
                     fmt='.')

def simulateAffinity(N, start=[10,40,15,10], delta=10, grid=100, margin=False):
    dists = [np.array(start)]
    for _ in range(N-1):
        ch = np.array([random.randint(-delta, delta) for _ in dists[-1]])
        dists.append(dists[-1]+np.array(ch))
        for i in range(len(dists[-1])):
            if dists[-1][i] < 0:
                dists[-1][i] = 0
    pAffinityChanges([x.tolist() for x in dists])
    na = [naiveAffinity(d) for d in dists]
    nca = [na[i+1]-na[i] for i in range(len(na)-1)]
    nca.insert(0, 0)
    print('in n-af n-afc | h1t h1o h2t h2o')
    for i in range(len(na)):
        print('%2d % .2f % .2f | %s'%(i, na[i], nca[i],
                                     ' '.join('%3d'%x for x in dists[i])))
