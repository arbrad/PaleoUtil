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
    
    # general setup
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
    
    def __init__(self, 
                 data=None, db=None,
                 time='time.csv', timeLevel=5,
                 taxaName='accepted_name'):
        assert data or db
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
        for sp, (e,l) in self.range.items():
            eb = bisect.bisect_left(self.timeSplits, e)
            lb = bisect.bisect_right(self.timeSplits, l)
            assert lb <= eb
            assert lb == 0 or self.timeSplits[lb-1] <= l
            assert eb == len(self.timeSplits) or self.timeSplits[eb] >= e
            self.bins[sp] = (lb, eb)
        # remove singletons: species that appear in only one bin
        for sp, _ in self.range.items():
            l, e = self.bins[sp]
            if l == e: del self.bins[sp]
    
    # Typical splitting function
    def splitByLocation(self, degrees=30, modifier=''):
        latstr, lngstr = modifier+'lat', modifier+'lng'
        def f(r):
            try:
                lat, lng = float(r[latstr]), float(r[lngstr])
                if lat > 180: lat -= 360
                if lat <= -180: lat += 360
                if lng > 180: lng -= 360
                if lng <= -180: lng += 360
                mlat, mlng = math.ceil(lat//degrees-0.5), math.ceil(lng//degrees-0.5)
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
            dbs[k] = DB(data=db, timeLevel=self.timeLevel, taxaName=self.taxaName)
        return dbs

    def species(self):
        '''Iterator over species in dataset.'''
        for sp in self.range:
            yield sp

    def fossilRange(self, species):
        '''Fossil range of species.'''
        return self.range[species]

    def interval(self, binId):
        '''Range corresponding to bin ID.'''
        e = 4500 if binId == len(self.timeSplits) else self.timeSplits[binId]
        l = 0 if binId == 0 else self.timeSplits[binId-1]
        return (e, l)

    def speciesByTime(self):
        '''List of species occurring in each time bin.'''
        bins = [[] for _ in range(len(self.timeSplits)+1)]
        for sp, (l,e) in self.bins.items():
            for binId in range(l, e+1):
                bins[binId].append(sp)
        return bins

    def computePCA(self, fields=['lithology*'], exclude=None):
        '''Fit PCA to dataset, using values from fields as Boolean dimensions.'''
        # observation dimensions: descriptors used in fields across dataset
        self.dims = list(self.valuesInFieldsOfData(fields, exclude))
        self.dims.sort()
        self.dim2i = {self.dims[i]:i for i in range(len(self.dims))}
        # data: frequency of observations per species
        sp = set()
        for r in self.data:
            sp.add(self.name(r))
        self.species_ = list(sp)
        self.species_.sort()
        self.sp2i = {self.species_[i]:i for i in range(len(self.species_))}
        # construct observation matrix
        self.Obs = np.zeros((len(self.species_), len(self.dims)))
        nocc = [0 for _ in self.species_]
        for r in self.data:
            si = self.sp2i[self.name(r)]
            nocc[si] += 1
            for dim in self.valuesInFields(fields, r, exclude):
                self.Obs[si,self.dim2i[dim]] += 1
        for r in range(self.Obs.shape[0]):
            for c in range(self.Obs.shape[1]):
                self.Obs[r,c] /= nocc[r]
                assert self.Obs[r,c] >= 0
                assert self.Obs[r,c] <= 1
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
    def colorByTime(self, summary=np.mean, start=4000):
        def color(sp):
            try:
                rng = self.fossilRange(sp)
                if rng[1] > start: return start
                return summary(rng)
            except:
                return -1
        return color
    def colorBySpecies(self, species):
        return lambda sp: int(sp in species)

    # Typical species-subset constructor for plotting
    def fieldSubset(self, fields, value, pol=True):
        '''Subset function constructor.'''
        if not hasattr(self, 'fieldSubset_'):
            self.fieldSubset_ = {}
        key = str(fields) + ':' + str(value) + ':' + str(pol)
        if key in self.fieldSubset_:
            return self.fieldSubset_[key]
        sp = set()
        if type(value) != str:
            value = set(value)
        for r in self.data:
            if ((type(value) == set and 
                 bool(value & self.valuesInFields(fields, r)) == pol)
                or (value in self.valuesInFields(fields, r)) == pol):
                sp.add(self.name(r))
        self.fieldSubset_[key] = sp
        return sp 
    def fieldSubsets(self, fields, species=None):
        '''Constructs value->set(species) map.'''
        vmap = {}
        for r in self.data:
            sp = self.name(r)
            if sp == '' or (species and sp not in species): continue
            for v in self.valuesInFields(fields, r):
                if v not in vmap: vmap[v] = set()
                vmap[v].add(sp)
        return vmap
    
    def plot(self, components=[0,1], colorOf=None, species=None, dimThresh=0.3):
        '''Plot projection of (subset of) dataset onto principal components.'''
        rows = self.speciesRows(species)
        P = self.pca.transform(self.Obs[rows,:])
        # color labels
        color = None
        if colorOf:
            color = []
            for sp in self.species_:
                if not species or sp in species:
                    color.append(colorOf(sp))
        # define scatter function
        def points(x, y):
            return P[:,x], P[:,y], None, color
        self.plotPCA(components, dimThresh, points)

    def plotPCA(self, components, dimThresh, points):            
        # draw each plot
        assert len(components) < 5
        sfm, sfn = [(1,1), (2,2), (2,3)][len(components)-2]
        fig = plt.figure()
        for i, comps in enumerate(itertools.combinations(components,2)):
            a, b, covs, c = points(comps[0], comps[1])
            ax = fig.add_subplot(sfm, sfn, i+1)
            # emphasize higher-valued colors
            if c:
                hasCovs = covs != None
                if not hasCovs: covs = [None for _ in a]
                x = list(zip(c,a,b,covs))
                x.sort()
                c,a,b,covs = tuple(list(y) for y in zip(*x))
                if not hasCovs: covs = None
            # plot
            im = ax.scatter(a,b,c=c,cmap=mplt.cm.jet,marker='.')
            if covs:
                # plot 95% confidence ellipses
                for i, cov in enumerate(covs):
                    if cov is None: continue
                    evals, evecs = np.linalg.eig(cov)
                    e0, e1 = evals[0], evals[1]
                    if not e0 or not e1: continue
                    v = evecs[:, 0 if e0 > e1 else 1]
                    # Chi-squared, 2 df, 95% conf
                    w = 2*(5.991*e0)**0.5
                    h = 2*(5.991*e1)**0.5
                    alpha = 180/math.pi*math.atan2(v[1],v[0])
                    el = mplt.patches.Ellipse((a[i],b[i]), w, h, alpha, 
                                              fill=False, ec=im.to_rgba(c[i]))
                    ax.add_patch(el)
            if c: plt.colorbar(im)
            # add primary contributing dimensions as vectors
            for j in range(len(self.dims)):
                if any(abs(self.pca.components_[a,j]) >= dimThresh for a in comps):
                    u,v = (self.pca.components_[a,j] for a in comps)
                    ax.quiver(0,0,u,v,angles='xy',scale_units='xy',scale=1,width=.005,color='red')
                    ax.annotate(self.dims[j], (u,v), color='red')
            ax.set_title(str(comps))

    def plotMeanPCAOverTime(self, components=[0,1], species=None, dimThresh=0.3,
                            start=4000, conf=True):
        '''Plot mean projection of species per time bin.'''
        def points(x, y):
            a, b, covs, c = [], [], [] if conf else None, []
            for i, spst in enumerate(self.speciesByTime()):
                if self.interval(i)[1] > start: continue
                sps = set(spst)
                if species: sps &= species
                if not sps: continue
                rows = self.speciesRows(sps)
                P = self.pca.transform(self.Obs[rows,:])
                a.append(np.mean(P[:,x]))
                b.append(np.mean(P[:,y]))
                if conf:
                    if len(rows) > 1:
                        covs.append(np.cov(P[:,[x,y]].transpose(), bias=True))
                    else:
                        covs.append(None)
                c.append(self.interval(i)[1])
            return a, b, covs, c
        self.plotPCA(components, dimThresh, points)                

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
    def valuesInFieldsOfData(self, fields, exclude=None):
        values = set()
        for r in self.data:
            for v in self.valuesInFields(fields, r, exclude):
                values.add(v)
        return values
    def valuesInFields(self, fields, r, exclude=None):
        if type(fields) != list: fields = [fields]
        values = set()
        for ff in fields:
            for f in self.expand(ff):
                for v in self.values(r[f], exclude):
                    values.add(v)
        return values

    def speciesRows(self, species):
        if species == None: return list(range(len(self.species_)))
        x = [self.sp2i[sp] for sp in species]
        x.sort()
        return x

def testPCA():
    db = DB('../../Paleo/Simpson/Corals/scler.082417a1.csv')
    db.computePCA()
    db.components()
    db.plot([0,1,2], db.colorByTime(), db.fieldSubset('order','Scleractinia'))
    return db

def speciesOverTime(db, degrees=None, modifier='', byBin=False, cutoffThresh=.05):
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
        loci = list(range(math.ceil(-179/degrees), math.ceil(180/degrees)+1))
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
        if byBin:
            # count by time bins, but normalize to species/my
            z = db.speciesByTime()
            x = [-db.interval(i)[1] for i in range(len(z))]
            # normalize to species/my, with a minimum of 1 my per interval
            y = [len(s)/max(1,db.interval(i)[0]-db.interval(i)[1]) for 
                 i, s in enumerate(z)]
        else:
            # use fossil ranges explicitly, by my
            x = [-my for my in range(2000)]
            y = [0 for _ in x]
            for sp in db.species():
                e, l = db.fossilRange(sp)
                for i in range(math.floor(l), math.floor(e)+1):
                    y[i] += 1
            # cut graph to to start at first significant time
            thresh = cutoffThresh * max(y)
            fnz = len(y)-1
            while y[fnz] < thresh: fnz -= 1
            x, y = x[:fnz+1], y[:fnz+1]
        ax.plot(x, y, '.')
        lr = LinearRegression()
        fz = 0
        for i in range(len(y)-1,-1,-1):
            if y[i]:
                fz = i+1
                break
        x, y = x[:fz], y[:fz]
        if len(x) > 1:
            weights = [1 for _ in x]
            if byBin:
                # weight by bin interval, taking care of final long one
                intervals = [db.interval(i)[0]-db.interval(i)[1] for i in range(len(y))]
                maxW = max(intervals[:-1])
                intervals[-1] = min(maxW, intervals[-1])
                weights = intervals
            lr.fit(np.array([[v] for v in x]), np.array(y), weights)
            a, b = lr.coef_, lr.intercept_
            ax.plot([x[0],x[-1]], [b,b+a*x[-1]])
        ax.set_title(str(loc))
        
