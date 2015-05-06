import os,sys
import Malt90Source
import aplpy
import gc as garbagecollection
import matplotlib
import pyfits
import malt_params as malt

class AgalSource:
    """
    A class to store information about each ATLASGAL source in MALT90
    
    Each Agal object has a Malt90 object based on the identification
    made in the linecatalog.
    
    """
    
    def __init__(self,atlasgal_id,cat_row,near_lons,near_lats,near_ids):
        self.agal_id = atlasgal_id
        self.agal_name = cat_row['source']
        #self.velocity = float(cat_row['consensus_velocity'])
        
        self.glon = float(cat_row['ag_long'])
        self.glat = float(cat_row['ag_lat'])
        
        self.near_lons = near_lons
        self.near_lats = near_lats
        
        self.awin = 38/60.
        self.bwin = 38/60.
        
        self.malt90_map = cat_row['malt90_map']
        self.malt90source = Malt90Source.Malt90Source(self.malt90_map,
                                                      self.agal_id)
        
    def make_herschel_result_images(self):
        try:
            self.malt90source.update_herschel_results()
        except:
        #Hack to deal with files missing Herschel coverage
        #if self.malt90source.atlasgal_mos == "None":
            d,h = pyfits.getdata(self.malt90source.filenames["M24c"],header=True)
            d[:] = 0
            pyfits.writeto(self.malt90source.filenames["RH_T"],d,h,clobber=True)
            pyfits.writeto(self.malt90source.filenames["RH_N"],d,h,clobber=True)

        #try:
        #    self.malt90source.update_atlasgal()
        #except ValueError:
        #Hack to deal with files missing ATLASGAL coverage
        #if self.malt90source.atlasgal_mos == "None":
        #    d,h = pyfits.getdata(self.malt90source.filenames["M24c"],header=True)
        #    d[:] = 0
        #    pyfits.writeto(self.malt90source.filenames["Agal"],d,h,clobber=True)

        
        Tim,Nim = self.malt90source.get_herschel_result_images()
        
        images = ["RH_Tim","RH_Nim"]

        #Add ATLASGAL-specific info        
        for im,name in zip([Tim,Nim],images):
        
            im.show_markers(self.near_lons,self.near_lats,edgecolor='blue',facecolor='none',marker='+',s=80,lw=1)
            im.show_markers(self.glon,self.glat,edgecolor='red',facecolor='none',marker='+',s=80,lw=1)
            im.show_ellipses(self.glon,self.glat,0.01,0.01,angle=0,facecolor='none',edgecolor='green',lw=2)
            out_filename = malt.image_dir+self.agal_id+"_"+name+".png"
            im.save(out_filename)
            print("Saving to: "+out_filename)
            im.close()
        matplotlib.pyplot.clf()
        matplotlib.pyplot.close() 
        matplotlib.pyplot.close() 
        
        matplotlib.pyplot.close() 
    
        garbagecollection.collect()
        
    def make_classification_images(self):
        #This checks the database, so is not wasteful
        self.malt90source.update_spitzer()
        #try:
        self.malt90source.update_atlasgal()
        #except ValueError:
        #Hack to deal with files missing ATLASGAL coverage
        #if self.malt90source.atlasgal_mos == "None":
        ##    print("Failed to update ATLASGAL")
        #    d,h = pyfits.getdata(self.malt90source.filenames["M24c"],header=True)
        #    d[:] = 0
        #    pyfits.writeto(self.malt90source.filenames["Agal"],d,h,clobber=True)
        
        #This returns references to each of the images
        Mim,Gim,Aim = self.malt90source.get_classification_images()
        images = ["Mim","Gim","Aim"]

        #Add ATLASGAL-specific info        
        for im,name in zip([Mim,Gim,Aim],images):
        
            im.show_markers(self.near_lons,self.near_lats,edgecolor='blue',facecolor='none',marker='+',s=80,lw=1)
            im.show_markers(self.glon,self.glat,edgecolor='red',facecolor='none',marker='+',s=80,lw=1)
            im.show_ellipses(self.glon,self.glat,0.01,0.01,angle=0,facecolor='none',edgecolor='green',lw=2)
            out_filename = malt.image_dir+self.agal_id+"_"+name+".png"
            im.save(out_filename)
            print("Saving to: "+out_filename)
            im.close()
        matplotlib.pyplot.clf()
        matplotlib.pyplot.close() 
        matplotlib.pyplot.close() 
        
        matplotlib.pyplot.close() 
    
        garbagecollection.collect()
        
        
        