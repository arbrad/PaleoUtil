# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 14:07:54 2017

@author: Aaron
"""

import bisect
import csv
import itertools
import math
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np
import re
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

class DB:
    '''Holds and analyzes PBDB dataset.'''
    
    # geneeral setup
    # resource directory
    rd = 'resources/'
    # fields to exclude from dataset
    excludeFields = re.compile('descript|comments|basis')
    # descriptors to exclude
    excludeValues = set(['','not reported','planktonic','buildup or bioherm',
                         'suspension feeder','microcarnivore',
                         'carbonate indet.','marine indet.','red','orange',
                         'yellow','green','blue','brown','black','gray',
                         'red or brown','white'])
    
    def __init__(self, db,
                 time='time.csv', timeLevel=4,
                 taxaName='accepted_name',
                 data=None):
        self.taxaName = taxaName
        # dataset
        if data:
            self.data = data
        else:
            with open(db, encoding='utf-8') as f:
                self.data = []
                for r in csv.DictReader(f):
                    if self.name(r) != '':
                        rr = {f:v for f,v in r.items() if not DB.excludeFields.search(f)}
                        self.data.append(rr)
        # time bins
        self.timeLevel = timeLevel
        with open(DB.rd+time) as f:
            my = [float(t['max_ma']) for t in csv.DictReader(f) if 
                  int(t['scale_level']) == timeLevel]
            my.sort()
            del my[-1]
            self.timeSplits = my
        # time ranges by species
        self.range = {}
        for r in self.data:
            sp = self.name(r)
            if sp not in self.range: self.range[sp] = (0, math.inf)
            ce,cl = self.range[sp]
            e,l = float(r['max_ma']), float(r['min_ma'])
            self.range[sp] = (max(ce,e), min(cl,l))
        # time bins by species
        self.bins = {}
        for sp,(e,l) in self.range.items():
            eb = bisect.bisect_right(self.timeSplits, e)
            lb = bisect.bisect_left(self.timeSplits, e)
            assert lb <= eb
            self.bins[sp] = list(range(lb,eb+1))
    
    # Typical splitting function
    def splitByLocation(self, degrees=30, modifier=''):
        latstr, lngstr = modifier+'lat', modifier+'lng'
        def f(r):
            try:
                lat, lng = float(r[latstr]), float(r[lngstr])
                if lat > 90: lat -= 180
                if lat <= -90: lat += 180
                if lng > 90: lng -= 180
                if lng <= -90: lng += 180
                mlat, mlng = math.ceil(lat/degrees), math.ceil(lng/degrees)
                return (mlat, mlng)
            except:
                return 'No location'
        return f

    def split(self, splitFn):
        '''Split this DB into many according to the splitFn.'''
        part = {}
        for r in self.data:
            k = splitFn(r)
            if k not in part: part[k] = []
            part[k].append(r)
        dbs = {}
        for k, db in part.items():
            dbs[k] = DB(db=None, timeLevel=self.timeLevel,
                        taxaName=self.taxaName, data=db)
        return dbs

    def period(self, binId):
        '''Return range corresponding to the bin ID.'''
        e = 4500 if binId == len(self.timeSplits) else self.timeSplits[binId]
        l = 0 if binId == 0 else self.timeSplits[binId-1]
        return (e, l)

    def speciesByTime(self):
        '''Return a list of species occurring in each time bin.'''
        bins = [[] for _ in range(len(self.timeSplits)+1)]
        for sp, binIds in self.bins.items():
            for binId in binIds:
                bins[binId].append(sp)
        return bins

    def computePCA(self, fields=['environment','lithology*']):
        '''Fit PCA to dataset, using values from fields as Boolean dimensions.'''
        # observation dimensions: descriptors used in fields across dataset
        self.dims = list(self.valuesInFieldsOfData(fields))
        self.dims.sort()
        self.dim2i = {self.dims[i]:i for i in range(len(self.dims))}
        # data: frequency of observations per species
        sp = set()
        for r in self.data:
            sp.add(self.name(r))
        self.species = list(sp)
        self.species.sort()
        self.sp2i = {self.species[i]:i for i in range(len(self.species))}
        # construct observation matrix
        self.Obs = np.zeros((len(self.species), len(self.dims)))
        nocc = [0 for _ in self.species]
        for r in self.data:
            si = self.sp2i[self.name(r)]
            nocc[si] += 1
            for dim in self.valuesInFields(fields, r):
                self.Obs[si,self.dim2i[dim]] += 1
        for r in range(self.Obs.shape[0]):
            for c in range(self.Obs.shape[1]):
                self.Obs[r,c] /= nocc[r]
        # compute PCA
        self.pca = PCA(whiten=True)
        self.pca.fit(self.Obs)
        self.components()

    def components(self, compThresh=0.3, contrThresh=0.02):
        '''Report contributing principal components.'''
        for r in range(self.pca.components_.shape[0]):
            if self.pca.explained_variance_ratio_[r] < contrThresh:
                break
            dims = []
            for i in self.pca.components_[r,:].argsort():
                if abs(self.pca.components_[r,i]) < compThresh:
                    continue
                dims.append('%s:%.2f'%(self.dims[i],self.pca.components_[r,i]))
            print('%2d %.2f %s'%(r,self.pca.explained_variance_ratio_[r],', '.join(dims)))

    # Typical coloring function constructors
    def colorByTime(self, extreme=min):
        return lambda sp: extreme(self.bins[sp])
    def colorBySpecies(self, species):
        return lambda sp: int(sp in species)

    # Typical species-subset constructor for plotting
    def fieldSubset(self, fields, value):
        '''Subset function constructor.'''
        sp = set()
        for r in self.data:
            if value in self.valuesInFields(fields, r):
                sp.add(self.name(r))
        return sp 
    
    def plot(self, components=[0,1], colorOf=None, species=None, dimThresh=0.3):
        '''Plot projection of (subset of) dataset onto principal components.'''
        assert len(components) < 5
        sfm, sfn = [(1,1), (2,2), (2,3)][len(components)-2]
        # select rows: subset if species list is given
        rows = list(range(len(self.species)))
        if species: rows = [self.sp2i[sp] for sp in species]
        # project
        P = self.pca.transform(self.Obs[rows,:])
        # color labels
        color = None
        if colorOf:
            color = []
            for sp in self.species:
                if not species or sp in species:
                    color.append(colorOf(sp))
        # draw each plot
        fig = plt.figure()
        for i, comps in enumerate(itertools.combinations(components,2)):
            a, b = P[:,comps[0]], P[:,comps[1]]
            ax = fig.add_subplot(sfm, sfn, i+1)
            # emphasize higher-valued colors
            if color:
                x = list(zip(color,a,b))
                x.sort()
                c,a,b = tuple(list(y) for y in zip(*x))
            else:
                c = None
            # plot
            im = ax.scatter(a,b,c=c,cmap=mplt.cm.nipy_spectral,marker='.')
            if color: plt.colorbar(im)
            # add primary contributing dimensions as vectors
            for j in range(len(self.dims)):
                if any(abs(self.pca.components_[a,j]) >= dimThresh for a in comps):
                    u,v = (self.pca.components_[a,j] for a in comps)
                    ax.quiver(0,0,u,v,angles='xy',scale_units='xy',scale=1,width=.005,color='red')
                    ax.annotate(self.dims[j], (u,v), color='red')

    def name(self, r):
        if self.taxaName not in r: return ''
        return r[self.taxaName]

    def expand(self, field):
        if field[-1] == '*':
            return [field[:-1]+str(s) for s in range(1,3)]
        else:
            return [field]

    def values(self, data, exclude=None):
        if exclude == None: exclude = DB.excludeValues
        rv = []
        for x in data.split(','):
            x = x.strip().strip('"')
            if x not in exclude:
                rv.append(x)
        return rv
    def valuesInFieldsOfData(self, fields):
        values = set()
        for r in self.data:
            for v in self.valuesInFields(fields, r):
                values.add(v)
        return values
    def valuesInFields(self, fields, r):
        if type(fields) != list: fields = [fields]
        values = set()
        for ff in fields:
            for f in self.expand(ff):
                for v in self.values(r[f]):
                    values.add(v)
        return values

def testPCA():
    db = DB('../../Paleo/Simpson/Corals/scler.082417a1.csv')
    db.computePCA()
    db.components()
    db.plot([0,1,2], db.colorByTime(), db.fieldSubset('order','Scleractinia'))
    return db

def speciesOverTime(db, degrees=None, modifier=''):
    '''Construct a plot of species-over-time by region of the world.
    Set modifer='paleo' to use paleocoordinates instead of modern
    coordinates. By default, the world is treated as one region.'''
    if degrees:
        dbs = db.split(db.splitByLocation(degrees=degrees, modifier=modifier))
    else:
        dbs = {(0,0):db}
    if 'No location' in dbs:
        print('No locations:', len(dbs['No location'].data))
        del dbs['No location']
    keys = set(dbs.keys())
    if degrees:
        loci = list(range(math.ceil(-89/degrees), math.ceil(90/degrees)+1))
        for x in loci:
            for y in loci:
                k = (x, y)
                if k not in keys:
                    keys.add(k)
    keys = list(keys)
    keys.sort(key=lambda x: (-x[0], x[1]))
    n = round(len(keys)**0.5)
    fig = plt.figure()
    for i, loc in enumerate(keys):
        if loc not in dbs: continue
        db = dbs[loc]
        ax = fig.add_subplot(n, n, i+1)
        z = db.speciesByTime()
        x = [-db.period(i)[1] for i in range(len(z))]
        y = [len(s) for s in z]
        ax.plot(x, y, '.')
        lr = LinearRegression()
        fz = 0
        for i in range(len(y)-1,-1,-1):
            if y[i]:
                fz = i+1
                break
        x, y = x[:fz], y[:fz]
        lr.fit(np.array([[v] for v in x]), np.array(y))
        a, b = lr.coef_, lr.intercept_
        ax.plot([x[0],x[-1]], [b,b+a*x[-1]])
        ax.set_title(str(loc))
        
