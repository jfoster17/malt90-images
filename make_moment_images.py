#!/usr/bin/env python
# encoding: utf-8
"""
"""

import sys
import os
import MALT90Source
import pickle
import matplotlib.pylab as plt
from matplotlib.ticker import MaxNLocator
from astropy.table import Table
import numpy as np
import malt90_catalog as mcat

lines = ["hnco413","c2h","sio","h41a",
         "hc13ccn","hnco404","ch3cn","hc3n",
         "h13cop","hn13c","13cs","13c34s",
         "hcop","hnc","n2hp","hcn",
         ]
        
mainlines = ["hcop","hnc","n2hp","hcn"]


def generate_pickle():
    all_cores = []
    t = mcat.read_latest("../results/malt90catalog.cat")

    #print(t)
    #f = open("../malt90_lineinfo.cat")
    for i,source in enumerate(t['agid']):
        if i < 10000:
            print(source)
            vel = float(t['velocity'][i])
            sourcename = t['malt90_map'][i]
            a_lon = float(t['ag_long'][i])
            a_lat = float(t['ag_lat'][i])

            distances = np.array(np.sqrt((t['ag_long']-a_lon)**2+(t['ag_lat']-a_lat)**2))
            min_distance = 3/60.

            nearby = t[(distances < min_distance)]
            #print(nearby)

            source_lons = nearby['ag_long']
            source_lats = nearby['ag_lat']
            source_ids = nearby['agid']

            all_cores.append(MALT90Source.MALT90Source(source,sourcename,vel,t[i],source_lons,source_lats,source_ids))

    return(all_cores)



def make_moment_images(all_cores):
    tagall = "AllLines"

    format = ".png"
    for core in all_cores:
        core.make_moment_images(lines,tagall,format=format,mom='mom0',outdir="/Volumes/Screwdriver/malt90/results/mom0/NewLines")




def main():
    #Of course unpickling does not work for astropy.table
    #Fixed in dev version, but hard to merge myself. 
    #Just generate the list every time.
    all_cores = generate_pickle()
    make_moment_images(all_cores)
    #do_spectra_plots(all_cores)


if __name__ == '__main__':
    main()
