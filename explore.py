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
        if not Ps[i]:
            if Rs[i]: return 0
            return multinomial(N, Rs[:i]+Rs[i+1:], Ps[:i]+Ps[i+1:])
    return math.exp(lmultinomial(N, Rs, Ps))

def pAffinity(dist=[0,250,0,1000], grid=100, margin=False, plot=True):
    '''For dist=[H1-taxa-fossil, H1-other, H2-taxa-fossil, H2-other],
    compute the probability distribution of affinity for H1. When
    margin=False, the computation uses the sample frequency of taxa
    fossils, f, and H1 fossils, h. When margin=True, the distribution 
    is computed by marginalizing out f and h from a joint a-f-h
    distribution. grid controls the mesh granularity.'''
    N = sum(dist)
    def g(a,f,h):
        # a = affinity for H1
        # f = frequency of taxa of interest
        # h = frequency of H1 samples compared to H2
        return multinomial(N, dist,
                           [a*f*h,
                            (1-a*f)*h,
                            (1-a)*f*(1-h),
                            (1-(1-a)*f)*(1-h)])
    if margin:
        # Naive mesh, will switch to MCMC
        mesh = np.linspace(0,1,grid+1)
        x = [sum(g(a,f,h) for f in mesh for h in mesh) for a in mesh]
    else:
        # Use sample means for f, h
        f = (dist[0]+dist[2])/N
        h = sum(dist[:2])/N
        x = [g(a,f,h) for a in np.linspace(0,1,grid+1)]
    unit = 1/grid
    s = sum(x)
    x = [e/s/unit for e in x]
    if plot:
        plt.plot(np.linspace(0,1,len(x)), x)
        tenth = math.floor(0.1*(len(x)-1))
        print('%.2f %.2f'%(unit*sum(x[:1+tenth]),unit*sum(x[-tenth-1:])))
    return x

def pAffinityChange(dist0, dist1, grid=100, margin=False):
    '''For two distributions as in pAffinity, cmpute the 
    probability distribution of the change in affinity for H1.'''
    p0 = pAffinity(dist0, grid, margin, plot=True)
    p1 = pAffinity(dist1, grid, margin, plot=True)
    x = [0 for _ in range(2*len(p0)-1)]
    mesh = np.linspace(0,1,grid+1)
    for i, a0 in enumerate(mesh):
        for j, a1 in enumerate(mesh):
            d = a1-a0
            x[int((d+1)*grid+0.5)] += p0[i]*p1[j]/grid
    plt.plot(np.linspace(-1,1,2*grid+1), x)
    print([d[0]/d[1]/(d[0]/d[1]+d[2]/d[3]) for d in [dist0,dist1]])