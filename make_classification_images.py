#!/usr/bin/env python
# encoding: utf-8
"""
Make Spitzer/Herschel images of MALT90 sources.

usage: make_classification_images [-h] [-s] [-d] [-i n] [-f m]

optional arguments:
-h       show this help message and exit
-s       ONLY make Spitzer images (default is to make both)
-d       ONLY make Herschel dust temp/column images (default both)
-i n     start with ATLASGAL id number n (useful if crashes)
-f m     end with ATLASGAL id number m (useful if crashes)

This script uses the complicated cut-out/cube/image generation
logic and tracking system built into the MALT90 object code. 
At this point, each MALT90 map has the appropriate cut-outs
already prepared from the large Spitzer and Herschel mosaics. 
This script takes those cutouts and makes three-color images
with the appropriate annotations and files names for each
ATLASGAL image. Spitzer images are use for classification, hence
the slightly confusing name.

Due to memory leak issues, it is possible that this script
will crash when run on the full catalog. It is possible to
specify a range of numbers (AGAL IDs) to run to circumvent this.

This script also contains some old code to generate a table
of sources that have bad Spitzer/MIPS data for some reason
in the appropriate formate to grab the WISE data. This process
is sort of done by hand, and is currently not working well.

"""


import os,sys,getopt
import malt90_catalog as mcat
import AgalSource
import AgalSourceWise
import numpy as np
import gc as garbagecollection
import malt_params as malt



def main():
    bad_mips = [] #List of sources which have bad MIPS/Spitzer
                  #Currently not being used?
    
    do_spitzer = True
    do_herschel = True
    
    try:
        opts,args = getopt.getopt(sys.argv[1:],"i:f:sdh")
    except getopt.GetoptError,err:
        print(str(err))
        print(__doc__)
        sys.exit(2)
    for o,a in opts:
        if o == "-i":
            start_id = int(a)
        elif o == "-f":
            final_id = int(a)
        elif o == "-s":
            do_spitzer = True
            do_herschel = False
        elif o == "-d":
            do_spitzer = False
            do_herschel = True
        elif o == "-h":
            print(__doc__)
            sys.exit(1)
        else:
            assert False, "unhandled option"
            print(__doc__)
            sys.exit(2)
    
    
    #Read in catalog
    t = mcat.read_latest("../results/malt90catalog.cat")
    for i,source in enumerate(t['agid']):
        num = float(source[2:])
        sourcename = t['malt90_map'][i]
        if (num > start_id) and (num < final_id):
            print("="*10+source+"="*10)
            sourcename = t['malt90_map'][i]
            print(sourcename)
            
            a_lon = t['ag_long'][i]
            a_lat = t['ag_lat'][i]
            distances = np.array(np.sqrt((t['ag_long']-a_lon)**2
                                 +(t['ag_lat']-a_lat)**2))
            min_distance = 3/60.
            nearby = t[(distances < min_distance)]
            source_lons = nearby['ag_long']
            source_lats = nearby['ag_lat']
            source_ids = nearby['agid']
            cat_row = t[i] 
            
            aim = malt.image_dir+source+"_Aim.png"
            gim = malt.image_dir+source+"_Gim.png"
            mim = malt.image_dir+source+"_Mim.png"
            
            #This if/else must have had a point at some point in the past
            if ((not os.path.isfile(aim)) or (not os.path.isfile(gim)) or (not os.path.isfile(mim))):
                do_images = True                
            else:
                do_images = True
            
            if (num in bad_mips) and do_images: 
                agal = AgalSourceWise.AgalSourceWise(source,cat_row,
                                             source_lons,source_lats,source_ids)
                if do_spitzer:
                    agal.make_classification_images()
                if do_herschel:
                    agal.make_herschel_result_images()
                del agal
            else:
                if do_images:
                    agal = AgalSource.AgalSource(source,cat_row,
                                             source_lons,source_lats,source_ids)
                    if do_spitzer:
                        agal.make_classification_images()
                    if do_herschel:
                        agal.make_herschel_result_images()
                else:
                    pass
                
            print("="*26)
            garbagecollection.collect()
            
def make_wise_lookup_table(t):
    """
    Old code to make a table for getting WISE data. Should not be needed.
    """
    import atpy
    ras = []
    decs = []
    
    #bad_hand = [4926]
    bad_hand = bad_mips
    
    for i,source in enumerate(t['ag_id']):
        num = float(source[2:])
        sourcename = t['malt90_map_filename'][i]
        a_lon = t['ag_long'][i]
        a_lat = t['ag_lat'][i]
        cat_row = t[i] 
        #if (a_lat > 1.) or (a_lat < -1) or (num in bad_hand) :
        if (num in bad_hand) :
            print(num)
            agal = AgalSource.AgalSource(source,cat_row,
                                         None,None,None)
            ras.append(agal.malt90source.apos_ra)
            decs.append(agal.malt90source.bpos_dec)
    tbl = atpy.Table()
    tbl.add_column('ra',np.array(ras),unit='degrees')
    tbl.add_column('dec',np.array(decs),unit='degrees')
    tbl.write("wise-table-for-lookup-more.tbl")

    

if __name__ == '__main__':
    main()
