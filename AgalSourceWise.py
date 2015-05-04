import os,sys
import Malt90Source
import Malt90SourceBadMips
import aplpy
import gc as garbagecollection
import matplotlib
import pyfits
import AgalSource

class AgalSourceWise(AgalSource.AgalSource):
    """
        Does WISE instead of MIPS
    """
    
    def __init__(self,atlasgal_id,cat_row,near_lons,near_lats,near_ids):
        self.agal_id = atlasgal_id
        self.agal_name = cat_row['ag_name']
        #self.velocity = float(cat_row['consensus_velocity'])
        
        self.glon = float(cat_row['ag_long'])
        self.glat = float(cat_row['ag_lat'])
        
        self.near_lons = near_lons
        self.near_lats = near_lats
        
        self.awin = 38/60.
        self.bwin = 38/60.
        
        self.malt90_map = cat_row['malt90_map_filename']
        self.malt90source = Malt90SourceBadMips.Malt90SourceBadMips(self.malt90_map,
                                                      self.agal_id)
