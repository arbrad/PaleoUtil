# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 16:43:55 2017

@author: Aaron

A file for exploratory code. Fruitful explorations will be rewritten into 
clean code.
"""

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
