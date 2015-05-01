#!/usr/bin/env python
# encoding: utf-8
"""
Make nice looking images of the 0th moment maps

This is a wrapper script that calls the functions
in MALT90Source which actually know how to display
the moment maps. In this wrapper we just load in
a catalog and generate images for all the entries.

Main other dependencies to check/edit:
MALT90Source
malt90_catalog
malt_params

Generally you will want to run make_agal_moments.py
before running this script to update the actual
moment maps before vizualizing them.

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
import malt_params as malt

lines = ["hnco413","c2h","sio","h41a",
         "hc13ccn","hnco404","ch3cn","hc3n",
         "h13cop","hn13c","13cs","13c34s",
         "hcop","hnc","n2hp","hcn",
         ]
        
mainlines = ["hcop","hnc","n2hp","hcn"]


def generate_pickle():
    all_cores = []
    t = mcat.read_latest(malt.base+"/results/malt90catalog.cat")

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

            source_lons = nearby['ag_long']
            source_lats = nearby['ag_lat']
            source_ids = nearby['agid']

            all_cores.append(MALT90Source.MALT90Source(source,sourcename,vel,t[i],source_lons,source_lats,source_ids))

    return(all_cores)



def make_moment_images(all_cores):
    tagall = "AllLines"
    format = ".png"
    for core in all_cores:
        core.make_moment_images(lines,tagall,format=format,mom='mom0',outdir=malt.base+"/results/mom0/NewLines")


def main():
    #Pickling does not work for old versions of astropy.table
    #Just generate the list every time for consistency.
    all_cores = generate_pickle()
    make_moment_images(all_cores)


if __name__ == '__main__':
    main()
