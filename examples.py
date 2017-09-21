# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 23:35:00 2017

@author: Aaron

Examples of how to use the mess o' code.
"""

from explore import *

def run():
    print('Reading data')
    data = readData()  # put the path to your PBDB CSV file here
    print('Kept', len(data), 'entries')
    
    # 1. Generate lithology diversity-affinity time series plots for Mollusca.
    
    # correct=True narrows suspicious time ranges
    # dropSingle (default: False) drops species appearing in only one time interval
    print('Creating EnvAffinity: lithology* diversity')
    lit = EnvAffinity(data, 'lithology*', carbonate, clastic, correct=True)
    print('Plotting for Mollusca and Brachiopoda')
    lit.focusOn('phylum')  # prints list of "interesting" phyla
    lit.plotByTime('Mollusca', start=450)
    lit.plotByTime('Brachiopoda', start=450)
    # look at Gastropoda and Bivalvia; example of restriction to taxa (optional)
    print('Plotting for Gastropoda and Bivalvia')
    lit.focusOn('class', restrictLevel='phylum', restrict='Mollusca')
    lit.plotByTime('Gastropoda', start=450)
    lit.plotByTime('Bivalvia', start=450)
    
    # 2. Generate lithology genus-level affinities for various taxa.
    print('Creating EnvAffinity: lithology* occurrence')
    gen = EnvAffinity(data, 'lithology*', carbonate, clastic, correct=True, macro=False)
    print('Plotting genera-level affinities for Gastropoda (blue), Bivalvia (orange), Bryozoa (green), Scleractinia (red)')
    gen.plotGeneraFor([('class', 'Gastropoda'), ('class', 'Bivalvia'), ('phylum', 'Bryozoa'), ('order', 'Scleractinia')])
    
    # Other combos: shallow/deep, reef/nonreef
