# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 23:35:00 2017

@author: Aaron

Examples of how to use the mess o' code.
"""

from explore import *
import scipy.stats

def run():
    print('Reading data')
    data = readData()  # put the path to your PBDB CSV file here
    print('Kept', len(data), 'entries')
    
    # 1. Generate lithology diversity-affinity time series plots for Mollusca.
    
    # correct=True narrows suspicious time ranges
    # dropSingle (default: False) drops species appearing in only one time interval
    print('Creating EnvAffinity: lithology* diversity')
    lit = EnvAffinity(data, 'lithology*', carbonate, clastic, correct=True)
    print('Plotting for Mollusca')
    lit.focusOn('phylum')  # prints list of "interesting" phyla
    lit.plotByTime('Mollusca', start=450)
    # look at Gastropoda; example of restriction to taxa (optional)
    print('Plotting for Gastropoda')
    lit.focusOn('class', restrictLevel='phylum', restrict='Mollusca')
    lit.plotByTime('Gastropoda', start=450)
    # Change-change plot for Gastropoda
    print('Plotting change-change')
    lit.plotChangeVsChange('Gastropoda')
    
    # 2. Occurrence-affinity time series.
    print('Creating EnvAffinity: lithology* occurrence')
    occ = EnvAffinity(data, 'lithology*', carbonate, clastic, correct=True, macro=False)
    print('Plotting for Mollusca')
    occ.focusOn('phylum')
    occ.plotByTime('Mollusca', start=450)
    
    # 3. Generate lithology genus-level affinities for various taxa.
    print('Creating EnvAffinity: lithology* occurrence')
    # timeLevel=1 creates (with time.csv) just one bin for 541-0
    gen = EnvAffinity(data, 'lithology*', carbonate, clastic, macro=False, timeLevel=1)
    print('Plotting genera-level affinities for Bryozoa')
    gen.plotGenera('Bryozoa', 'phylum')
    print('Plotting genera-level affinities for Gastropoda (blue), Bivalvia (orange), Bryozoa (green), Scleractinia (red)')
    gen.plotGeneraFor([('class', 'Gastropoda'), ('class', 'Bivalvia'), ('phylum', 'Bryozoa'), ('order', 'Scleractinia')])
    
    # Other combos: shallow/deep, reef/nonreef

def generateFigures(data=None, fdir='../../Paleo/Affinity/', **kwargs):
    saveOpts = {'bbox_inches':'tight', 'dpi':600}
    plt.figure()
    for a, b in [(1,1), (1,2), (0,2), (9,12), (9,24)]:
        pAffinityPct(a, b, 0.5, **kwargs)
    plt.savefig(fdir+'scenario1.png', **saveOpts)
    plt.figure()
    p0 = pAffinityPct(10, 7, 0.5, True, **kwargs)
    p1 = pAffinityPct(5, 17, 0.5, True, **kwargs)
    pAffinityDiff(None, None, True, p1, p0, **kwargs)
    plt.xlabel('(Change in) Affinity')
    plt.savefig(fdir+'scenario2.png', **saveOpts)
    plt.figure()
    for a, b in [(1,4), (10,40)]:
        pAffinityPct(a, b, 0.3, **kwargs)
    plt.savefig(fdir+'scenario4.png', **saveOpts)
    plt.figure()
    grid = int(kwargs['grid']) if 'grid' in kwargs else 100
    mesh = np.linspace(1/grid, 1-1/grid, grid-1)
    priorA = lambda x: scipy.stats.beta.pdf(x, 4, 4)
    priorH = lambda x: scipy.stats.beta.pdf(x, 12, 12)
    plt.plot(mesh, [priorH(x) for x in mesh])
    plt.plot(mesh, [priorA(x) for x in mesh])
    pAffinityPct(3, 5, 0.5, True, **kwargs)
    pAffinityPct(3, 5, 0.5, True, priorA=priorA, priorH=priorH, **kwargs)
    pAffinityPct(3, 5, 0.5, True, priorH=priorH, **kwargs)
    plt.savefig(fdir+'prior.png', **saveOpts)
    plt.figure()
    pAffinityPct(1, 4, 0.7)
    for c, d in [(13,2), (27,8), (55,20)]:
        pAffinity([1, c, 4, d], margin=True)
    plt.xlim([0, 0.5])
    plt.savefig(fdir+'scenario5.png', **saveOpts)
    if data:
        generateGeneraPlots(data, fdir, saveOpts)

def generateGeneraPlots(data, fdir, saveOpts):
    lit = EnvAffinity(data, 'lithology*', carbonate, clastic, macro=False, timeLevel=1)
    lit.plotGeneraFor([('class', 'Bivalvia'), ('phylum', 'Bryozoa'), ('order', 'Scleractinia')], 
                      multiPlot=True,
                      figsize=(4,8))
    plt.savefig(fdir+'scenario3.png', **saveOpts)
    lit.plotGeneraFor([('genus', 'Abra (Abra)'), ('genus', 'Daonella'), ('genus', 'Actinostreon'), ('genus', 'Toucasia')],
                      grid=1000)
    plt.savefig(fdir+'scenario3_ex.png', **saveOpts)