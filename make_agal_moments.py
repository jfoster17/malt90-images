#!/usr/bin/env python
# encoding: utf-8
"""
make_agal_moments

Make moment maps using the catalog from Jill/Scott

Make a moment map for each ATLASGAL source. Simpler
than 

"""

import sys,os
import numpy as np

import malt_params as malt
import malt90_catalog as mcat
import moment_map as moment_map

def main():
    t = mcat.read("/Users/jonathanfoster/Desktop/Current/Malt90/malt90_lineinfo.cat")
    for i,source in enumerate(t['ag_id']):
        number = float(source[2:])
        if number > 808: #To allow me to specify
            malt90_map = t['malt90_map_filename'][i]
            print("Doing source: "+source)
            print("With MALT90 Map: "+malt90_map)
            vel = t['consensus_velocity'][i]
            if vel < -300 or vel > 300:
                vel = 0
            moment_map.do_source(malt90_map,malt.lines,
                                     outname=source,vel=vel)

if __name__ == '__main__':
    main()
