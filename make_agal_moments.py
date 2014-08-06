#!/usr/bin/env python
# encoding: utf-8
"""
make_agal_moments

Make moment maps using the catalog from Jill/Scott

Make a moment map for each ATLASGAL source. Simpler
than anything else.

6/5/14: New ATLASGAL names/catalog requires a remake

Using the consenvels.dat (Thu Jun 19 01:20:25 2014) 
from Scott. Need to cross-correlate with the main list
to get AGAL IDs


"""

import sys,os
import numpy as np

import malt_params as malt
import malt90_catalog as mcat
import moment_map as moment_map

from astropy.table import Table
def main():
    #t = mcat.read("/Users/jonathanfoster/Desktop/Current/Malt90/malt90_lineinfo.cat")
    t = mcat.readnew("/Users/jonathanfoster/Desktop/Current/Malt90/ATLASGALsources_coveredby_MALT90.list")
    vv = Table.read('consenvels.dat',names=["source","novbest","vbest","veldev"],format="ascii")
    for i,source in enumerate(t['ag_id']):
        number = float(source[2:])
        agalname = t['ag_name'][i]
        goodrow = vv[vv['source'] == agalname]
        malt90_map = t['malt90_map_filename'][i]
        malt90_pos = t['malt90_pos'][i]
        
        if (number > 0) and ((malt90_pos == "YCM") or (malt90_pos == "YEM")): #To allow me to specify
            print("Doing source: "+source)
            print("With MALT90 Map: "+malt90_map)
            print("This is AGAL: "+agalname)
            vel = goodrow['vbest']
            if vel < -300 or vel > 300:
                vel = 0
            if not goodrow:
                vel = 0
            moment_map.do_source(malt90_map,malt.lines,
                                     outname=source,vel=vel)

if __name__ == '__main__':
    main()
