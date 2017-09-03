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

def plotPhoto(db, comps):
    scler = db.fieldSubset('order', 'Scleractinia')
    photo = db.fieldSubset('diet', 'photosymbiotic') & scler
    notphoto = db.fieldSubset('diet', 'photosymbiotic', False) & scler
    both = photo & notphoto
    db.plot(comps, lambda sp: sum(int(sp in x) for x in [scler,notphoto,both]))

def photoAnalysis():
    db = DB('../../Paleo/Simpson/Corals/pbdb_animalia_marine_090217.csv')
    db.computePCA(['environment'])
    plotPhoto(db, [0,3])
    db.computePCA(['lithology*'])
    plotPhoto(db, [0,1])