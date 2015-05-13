#!/usr/bin/env python
# encoding: utf-8
"""
Make moment maps of ATLASGAL sources

usage: make_agal_moments [-h] [-i n] 

optional arguments:
-h       show this help message and exit
-i n     start with ATLASGAL id number n (useful if crashes)

A given MALT90 map may contain multiple ATLASGAL sources.
Furthermore, each ATLASGAL source may be associated 
with multiple velocity components, so we have to make 
moment maps around each velocity component. This script
will generate moment maps in malt90_images/mommaps/ and
this directory will be organized like the main MALT90 
moment maps directory, but will have an entry for every
ATLASGAL velocity component, named with their AGALID.
"""

import sys,os,getopt

import malt_params as malt
import malt90_catalog as mcat
import moment_map as moment_map

def main():
    start_id = 0

    try:
        opts,args = getopt.getopt(sys.argv[1:],"i:h")
    except getopt.GetoptError,err:
        print(str(err))
        print(__doc__)
        sys.exit(2)
    for o,a in opts:
        if o == "-i":
            start_id = int(a)
        elif o == "-h":
            print(__doc__)
            sys.exit(1)
        else:
            assert False, "unhandled option"
            print(__doc__)
            sys.exit(2)
    
    t = mcat.read_latest(malt.base+"/results/malt90catalog.cat")
    for i,source in enumerate(t['agid']):
        number = float(source[2:])
        agalname = t['source'][i]
        malt90_map = t['malt90_map'][i]
        vel = float(t['velocity'][i])
        
        if (number > start_id): #To allow me to specify a start
            print("Doing source: "+source)
            print("With MALT90 Map: "+malt90_map)
            print("This is AGAL: "+agalname)
            if vel < -300 or vel > 300: #Set junk velocities to 0
                vel = 0
            moment_map.do_source(malt90_map,malt.lines,
                                     outname=source,vel=vel)

if __name__ == '__main__':
    main()
