import os,sys
import matplotlib
matplotlib.use('Agg')
import reproject_map
import astLib.astCoords as ac
import astLib.astWCS as astWCS
import pyfits
import numpy as np
import aplpy
import numdisplay.zscale as zscale
import montage
import gc as garbagecollection
import matplotlib
import bisect
import glob
import shutil
#Somehow the new agpy breaks mpfit install, rendering this problematic
#just need to install my own new version of mpfit, probably.
#import herschel_utilities
from astroquery.magpis import Magpis
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import malt_params as malt


class Malt90Source:
    """A class to store information about a Malt90Source"""

    def __init__(self,name,agal_id):
        self.name = name
        self.safe_name = name.replace('.','_')
        self.id   = agal_id
        try:
            self.glon = float(self.name[1:7])
            self.glat = float(self.name[8:14])
        except ValueError: #This is a not-observed source
            pass

        radec = ac.convertCoords("GALACTIC","J2000",self.glon,self.glat,2000.)
        self.ra  = radec[0]
        self.dec = radec[1]

        self.comments = self.get_comments()


        ### Define the cutouts ###
        self.cutouts = ["M24c","IR4c","IR3c","IR2c","IR1c",
                        "TMkc","TMhc","TMjc",
                        "H500c","H350c","H250c","H170c","H70c",
                        "Agal","RH_T","RH_N",
                        ]

        ### Define the standard cubes ###
        self.cubes = {"Mcube":["M24c","IR4c","IR1c"],
                      "Gcube":["IR4c","IR2c","IR1c"],
                      "Tcube":["TMkc","TMhc","TMjc"],
                      "Hcube":["H500c","H350c","H250c"]
                      }

        ### Define the standard images ###
        self.images = {"Mim":["Mcube"],
                       "Gim":["Gcube"],
                       "Tim":["Tcube"],
                       "Him":["Hcube"],
                       "Aim":["Agal"],
                       "RH_Tim":["RH_T"],
                       "RH_Nim":["RH_N"],
                    }

        self.storage_dir = malt.base+"/results/"
        self.auxiliary_log = self.storage_dir+"aux_data.log"

        ### Define the filename conventions ###
        self.set_filenames()
        #This seems possibly redundant
        self.apos = self.glon
        self.bpos = self.glat
        try:
            posradec = ac.convertCoords("GALACTIC","J2000",self.apos,self.bpos,2000.)
            self.apos_ra = posradec[0]
            self.bpos_dec = posradec[1]
        except IndexError:
            print("Failed to convert coordinates")
            self.apos = self.glon
            self.bpos = self.glat
        self.data_base = malt.base+"/sources/"
        self.line_name_trans = {"n2hp":r"$\mathrm{N_2H^+}$","hcop":r"$\mathrm{HCO^+}$","hcn":r"$\mathrm{HCN}$","hnc":r"$\mathrm{HNC}$",
                                "13c34s":r"$\mathrm{^{13}C^{34}S}$","13cs":r"$\mathrm{^{13}CS}$","c2h":r"$\mathrm{C_2H}$","ch3cn":r"$\mathrm{CH_3CN}$",
                                "h13cop":r"$\mathrm{H^{13}CO^{+}}$","h41a":r"$\mathrm{H41\alpha}$","hc3n":r"$\mathrm{HC_3N}$",
                                "hc13ccn":r"$\mathrm{HC^{13}CCN}$","hn13c":r"$\mathrm{HN^{13}C}$","hnco404":r"$\mathrm{HNCO 4_{0,4}}$",
                                "hnco413":r"$\mathrm{HNCO 4_{1,3}}$","sio":r"$\mathrm{SiO}$"}
        self.line_spectra    = {"n2hp":[],"hcop":[],"hcn":[],"hnc":[],
                                "13c34s":[],"13cs":[],"c2h":[],"ch3cn":[],
                                "h13cop":[],"h41a":[],"hc3n":[],
                                "hc13ccn":[],"hn13c":[],"hnco404":[],
                                "hnco413":[],"sio":[]}
                   
        self.moment_name_trans = {"snr0":r"$\mathrm{SNR}$"}
    
    def get_moment_map_path(self,molecule,moment):
        mol = "_"+molecule+"_"
        a = os.path.join(self.data_base,self.name,self.name+mol+"mommaps",self.name+mol+moment+".fits")
        return(a)

    def get_cube_path(self,molecule):
        mol = "_"+molecule+"_"
        a = os.path.join(self.data_base,self.name,self.name+mol+"MEAN.fits")
        return(a)
    
    def get_herschel_path(self,wavelength):
        a = os.path.join(self.storage_dir,"cutouts",self.safe_name+"_H"+wavelength+"c.fits")
        return(a)
    
    def do_twomass_images(self):
        self.do_image("Tim")
        
        
    def get_moment_at_location(self,molecule,moment,xx,yy):
        moment_map = self.get_moment_map_path(molecule,moment)
        data = pyfits.getdata(moment_map)
        val = data[yy,xx]
        return(val)
    
    def get_max_snr_in_center(self,molecule):
        snr = self.get_moment_map_path(molecule,'snr0')
        data = np.nan_to_num(pyfits.getdata(snr))
        max_val = np.max(data[8:18,8:18])
        return(max_val)
    
    def do_snr_images(self):
        mainlines = {"n2hp":"yellow","hcop":"cyan","hcn":"orange","hnc":"magenta"}
        otherlines = {"13c34s":"aliceblue","13cs":"white","c2h":"coral","ch3cn":"bisque","h13cop":"aqua","h41a":"crimson","hc3n":"chartreuse","hc13ccn":"deeppink",
                     "hn13c":"fuchsia","hnco404":"gray","hnco413":"powderblue","sio":"gold"}
        self.identify_spitzer_mosaics()
        levels = np.arange(0,2,3) #This is bad now
        for line in mainlines.keys():
        #    #print(line)
            snr = self.get_moment_map_path(line,"snr0")
            self.do_image("Mim",contour_file=snr,contour_levels=levels,contour_color=mainlines[line])
            self.do_image("Gim",contour_file=snr,contour_levels=levels,contour_color=mainlines[line])
        for line in otherlines.keys():
            #print(line)
            snr = self.get_moment_map_path(line,"snr0")
            self.do_image("Mim",contour_file=snr,contour_levels=levels,contour_color=otherlines[line])
            self.do_image("Gim",contour_file=snr,contour_levels=levels,contour_color=otherlines[line])
    
    def get_classification_images(self):
        self.identify_spitzer_mosaics()
        self.identify_atlasgal_mosaics()
        
        #Set these levels to be appropriate. Somehow
        levels = np.arange(0.2,4,0.5)
        atlasgal = self.filenames["Agal"]
        Mim = self.do_image("Mim",contour_file=atlasgal,contour_levels=levels,contour_color='yellow')
        Gim = self.do_image("Gim",contour_file=atlasgal,contour_levels=levels,contour_color='yellow')
        Aim = self.do_grayscale("Aim",contour_file=atlasgal,contour_levels=levels,contour_color='yellow')
        return(Mim,Gim,Aim)
    
    def get_herschel_result_images(self):
        self.identify_herschel_result_mosaics_new()
        self.identify_atlasgal_mosaics()
        
        levels = np.arange(0.2,4,0.5)
        atlasgal = self.filenames["Agal"]
        Tim = self.do_hcolorscale("RH_Tim",contour_file=atlasgal,contour_levels=levels,contour_color='yellow',vmin=10,vmax=30,color="hot")
        Nim = self.do_hcolorscale("RH_Nim",contour_file=atlasgal,contour_levels=levels,contour_color='yellow',vmin=0.01,vmax=0.3,color="gist_yarg",do_log=True)
        return(Tim,Nim)
    
    def set_filenames(self):
        """ Define filename conventions """
        cutouts = self.storage_dir+"cutouts"
        cubes   = self.storage_dir+"cubes"
        images  = self.storage_dir+"images"
        self.filenames = {}
        for obj,datapath in zip([self.cutouts,self.cubes,self.images],[cutouts,cubes,images]):
            for dd in obj:
                self.filenames[dd] = datapath+"/"+self.safe_name+"_"+dd+".fits"
        pass


    def get_comments(self):
        """Get freeform comments about a source from text database."""
        pass


    def update_atlasgal(self):
        """Make Atlasgal sub-mosiacs and images"""
        cutouts = ["Agal"]
        cubes = []
        images  = ["Aim"]
        self.identify_atlasgal_mosaics()
        self.update_files(cutouts,cubes,images,"Atlasgal")

    def update_herschel_results(self):
        """Make Herschel result sub-mosiacs and images
        RH = "Results Herschel"
        
        """
        cutouts = ["RH_T","RH_N"]
        cubes = []
        images  = ["RH_Tim","RH_Nim"]
        self.identify_herschel_result_mosaics_new()
        self.update_files(cutouts,cubes,images,"HerschelResult")

    def update_spitzer(self):
        """Make Spitzer sub-mosiacs and images"""
        cutouts = ["M24c","IR1c","IR2c","IR3c","IR4c"]
        cubes   = ["Mcube","Gcube"]
        images  = ["Mim","Gim"]
        self.identify_spitzer_mosaics()
        self.update_files(cutouts,cubes,images,"Spitzer")

    def update_herschel(self):
        """Make Herschel sub-mosiacs and images"""
        cutouts = ["H500c","H350c","H250c"]#,"H170c","H70c"]
        cubes   = ["Hcube"]
        images  = ["Him"]
        self.identify_herschel_mosaics()
        self.update_files(cutouts,cubes,images,"Herschel")

    def update_twomass(self):
        """Make 2MASS sub-mosiacs and images"""
        cutouts = ["TMkc","TMhc","TMjc"]
        cubes   = ["Tcube"]
        images  = ["Tim"]
        self.update_files(cutouts,cubes,images,"Twomass")

    def get_undone(self,test_keys):
        """Figure out which of test_keys needs to be processed"""
        undone_keys = []
        curr_version,data_version = self.read_auxiliary_log()
        for test_key in test_keys:
            #print(self.name)
            #print(data_version[test_key])
            #print(curr_version[test_key])
            if data_version[test_key] < curr_version[test_key]:
                undone_keys.append(test_key)
        return(undone_keys)

    def update_files(self,cutouts,cubes,images,survey_type):
        """ Generically update a tree of cutouts, cubes and images """
        ### Do cut-outs ###
        undone_cutouts = self.get_undone(cutouts)
        for undone_cutout in undone_cutouts:
            try:
                print("Doing "+undone_cutout)
                self.do_cutout(undone_cutout)
            #except IOError:
            except SyntaxError:
                print("Syntax Error")
                #pass
                #Maybe I should restore catching this to make it more robust?
                #print(sys.exc_info()[0])

        ### Do image cubes ###
        #undone_cubes should be a set so that append->add without fuss about dupes
        undone_cubes = self.get_undone(cubes)
        for cube in self.cubes:
            for undone_cutout in undone_cutouts:
                if undone_cutout in cube:
                    undone_cubes.append(cube) #Redo cubes based on new cut-outs
        for undone_cube in undone_cubes:
            try:
                self.do_cube(undone_cube)
            except:
                pass

        ### Do images ###
        undone_images = self.get_undone(images)
        for image in self.images:
            for undone_cube in undone_cubes:
                if undone_cube in image:
                    undone_images.append(image) #Redo images based on new images
        for undone_image in undone_images:
            try:
                self.do_image(undone_image)
            except:
                print("Failed to make image")
                pass
    def do_grayscale(self,undone_image,contour_file=None,contour_levels=5,contour_color=None):
        """
        Make a grayscale image. For now, this means ATLASGAL.
        """
        version = 4
        success = False
        
        print("Making an image for "+undone_image)
        print("Using files... "+str(self.images[undone_image]))
        print("Full files are: "+self.filenames[self.images[undone_image][0]])
        
        print(contour_file)
        
        gc = aplpy.FITSFigure(self.filenames[self.images[undone_image][0]])
        gc.show_grayscale()
        gc.tick_labels.set_xformat("dd.dd")  
        gc.tick_labels.set_yformat("dd.dd")
        if contour_file:
            #gc.show_contour(contour_file,levels=contour_levels,colors="black",linewidths = 3.,smooth=3.,zorder=1,alpha=0.5)
            gc.show_contour(contour_file,levels=contour_levels,colors=contour_color,linewidths = 1.,zorder=1)
        success = True
        if success:
            self.report_success(undone_image,version)
        return(gc)
        
    def do_hcolorscale(self,undone_image,contour_file=None,contour_levels=5,contour_color=None,vmin=0,vmax=100,color="gist_heat",do_log=False):
        """
        Make a colorscale image for Herschel.
        """
        version = 4
        success = False
    
        print("Making an image for "+undone_image)
        print("Using files... "+str(self.images[undone_image]))
        print("Full files are: "+self.filenames[self.images[undone_image][0]])
    
    
        print(contour_file)
    
        gc = aplpy.FITSFigure(self.filenames[self.images[undone_image][0]])
        if do_log:
            stretch = 'log'
        else:
            stretch = 'linear'
        gc.show_colorscale(cmap=color,vmin=vmin,vmax=vmax,stretch=stretch)
        gc.add_colorbar()
        gc.colorbar.show()
        gc.tick_labels.set_xformat("dd.dd")  
        gc.tick_labels.set_yformat("dd.dd")
        if contour_file:
            #gc.show_contour(contour_file,levels=contour_levels,colors="black",linewidths = 3.,smooth=3.,zorder=1,alpha=0.5)
            gc.show_contour(contour_file,levels=contour_levels,colors=contour_color,linewidths = 1.,zorder=1)
        success = True
        if success:
            self.report_success(undone_image,version)
        return(gc)
        
        

    def do_image(self,undone_image,contour_file=None,contour_levels=5,contour_color=None):
        """Make an image"""
        version = 4
        success = False
        print("Making an image for "+undone_image)
        print("Using files... "+str(self.images[undone_image]))
        
        cube_data = pyfits.getdata(self.filenames[self.images[undone_image][0]])
        z1_k,z2_k = zscale.zscale(np.nan_to_num(cube_data[0,...]))
        z1_h,z2_h = zscale.zscale(np.nan_to_num(cube_data[1,...]))
        z1_j,z2_j = zscale.zscale(np.nan_to_num(cube_data[2,...]))
        
        aplpy.make_rgb_image(self.filenames[self.images[undone_image][0]],
                            self.filenames[self.images[undone_image][0]].replace('.fits','_2d.png'),
                            vmin_r = z1_k,vmax_r = z2_k,vmin_g = z1_h, vmax_g = z2_h, vmin_b = z1_j, vmax_b = z2_j)
        gc = aplpy.FITSFigure(self.filenames[self.images[undone_image][0]].replace('.fits','_2d.fits'))
        gc.show_rgb(self.filenames[self.images[undone_image][0]].replace('.fits','_2d.png'))
        gc.tick_labels.set_xformat("dd.dd")  
        gc.tick_labels.set_yformat("dd.dd")
        #out_filename = self.filenames[undone_image].replace('.fits','.pdf')
        
        if contour_file:
            #gc.show_contour(contour_file,levels=contour_levels,colors="black",linewidths = 3.,smooth=3.,alpha=0.5)
            gc.show_contour(contour_file,levels=contour_levels,colors=contour_color,linewidths = 1.)
            #chunks = contour_file.split('_')
            #print(chunks)
            #linename = chunks[3]
            #moment_name = chunks[4][0:-5]
            #tag = "_"+linename+"_"+moment_name
            #print(tag)
            #out_filename = out_filename.replace('.pdf',tag+'.pdf')
            #self.add_labels(gc,undone_image,linename,moment_name = moment_name,color=contour_color)
        else:
            pass
            #self.add_labels(gc,undone_image)
            
        #gc.show_rectangles(self.apos,self.bpos,0.0625,0.0625,ec='w',lw=2,facecolor='none',zorder=1)
        #gc.show_ellipses(self.apos,self.bpos,self.awin*2,self.bwin*2,angle=self.twin,facecolor='none',ec='green',lw=2)
        
        #print(self.awin,self.bwin)
        success = True
        if success:
            self.report_success(undone_image,version)
        return(gc)
    
    def set_label_text(self,undone_image):   
        if undone_image.startswith("M"):
            label1 = r"$\mathrm{GLIMPSE}\/3.6\mu \mathrm{m}$"
            label2 = r"$\mathrm{GLIMPSE}\/8.0\mu \mathrm{m}$"
            label3 = r"$\mathrm{MIPS}\/24\mu \mathrm{m}$"
        elif undone_image.startswith("G"):
            label1 = r"$\mathrm{GLIMPSE}\/3.6\mu \mathrm{m}$"
            label2 = r"$\mathrm{GLIMPSE}\/4.5\mu \mathrm{m}$"
            label3 = r"$\mathrm{GLIMPSE}\/8.0\mu \mathrm{m}$"
        elif undone_image.startswith("T"):
            label1 = r"$\mathrm{2MASS}\/J$"
            label2 = r"$\mathrm{2MASS}\/H$"
            label3 = r"$\mathrm{2MASS}\/K$"
            
        return(label1,label2,label3)
    
    def add_labels(self,gc,undone_image,linename=None,moment_name = None,color=None):
        label1,label2,label3 = self.set_label_text(undone_image)
        if linename:
            gc.show_rectangles(self.apos-0.031,self.bpos+0.031,0.017,0.014,facecolor='white',zorder=5)
            gc.show_rectangles(self.apos-0.031,self.bpos+0.026,0.017,0.005,facecolor='black',zorder=5)
            gc.add_label(0.79,0.83,self.line_name_trans[linename]+" "+self.moment_name_trans[moment_name]+" ---",
                        color=color,relative="True", weight="extra bold",horizontalalignment="left",zorder=10)
            #contour_positions = np.zeros([2,2])
            #contour_positions[0,0] = self.apos-0.038
            #contour_positions[0,1] = self.apos-0.034
            ##contour_positions[1,0] = self.bpos+0.0265
            #contour_positions[1,1] = self.bpos+0.0265
            #gc.show_lines([contour_positions],color=color,zorder=10)
        else:
            gc.show_rectangles(self.apos-0.031,self.bpos+0.031,0.017,0.014,facecolor='white')
            
        gc.add_label(0.79,0.95,label1,color='blue',relative="True", weight="extra bold",horizontalalignment="left",zorder=10)
        gc.add_label(0.79,0.91,label2,color='green',relative="True", weight="extra bold",horizontalalignment="left",zorder=10)
        gc.add_label(0.79,0.87,label3,color='red',relative="True", weight="extra bold",horizontalalignment="left",zorder=10)
            

    def do_cube(self,undone_cube):
        """Make a cube"""
        version = 2
        success = False
        print("Making a cube for "+undone_cube)
        print("Using files... "+str(self.cubes[undone_cube]))
        aplpy.make_rgb_cube([self.filenames[self.cubes[undone_cube][0]],self.filenames[self.cubes[undone_cube][1]],
                             self.filenames[self.cubes[undone_cube][2]]],self.filenames[undone_cube],system="GAL")
        print("Saving to: "+self.filenames[undone_cube])
        success = True
        if success:
            self.report_success(undone_cube,version)

    def get_all_spectra(self,option):
        """Get the spectra for all lines as specified by option."""
        if option == "higal500":
            herschel_info = self.get_higal_peak_info("500")
            for line in self.line_spectra.keys():
                self.line_spectra[line] = self.get_spectra_at_higal_peak(line,herschel_info)

    def get_higal_peak_info(self,wavelength):
        """Obtain the position and column density at the peak of the (wavelength) HiGAL map"""
        datapath = self.get_cube_path("n2hp") #We assume that the WCS is the same for all lines
        d,h = pyfits.getdata(datapath,header=True)
        malt90_WCS = astWCS.WCS(h,mode="pyfits")
        hipath = self.get_herschel_path("500")
        dd_500,hh_500 = pyfits.getdata(hipath,header=True)
        herschel_WCS_500 = astWCS.WCS(hh_500,mode="pyfits")
        maximum = 0
        bestxx = 10
        bestyy = 10
        for xx in range(8,16):
            for yy in range(8,16):
                #print(xx,yy,dd[yy,xx])
                if dd_500[yy,xx] > maximum:
                    maximum = dd_500[yy,xx]
                    bestxx = xx
                    bestyy = yy
        #print(bestxx,bestyy)
        radec = herschel_WCS_500.pix2wcs(bestxx,bestyy)
        xy = malt90_WCS.wcs2pix(radec[0],radec[1])
        #xy = [14,12]
        #print(xy)
        n2hp = self.get_moment_at_location("n2hp","mom0",xy[1],xy[0])
        hcop = self.get_moment_at_location("hcop","mom0",xy[1],xy[0])
        ratio = hcop/n2hp
        print(ratio)
        
        hipath = self.get_herschel_path("250")
        dd_250,hh_250 = pyfits.getdata(hipath,header=True)
        herschel_WCS_250 = astWCS.WCS(hh_250,mode="pyfits")
        xy_250 = herschel_WCS_250.wcs2pix(radec[0],radec[1])
        
        S_500 = maximum
        S_250 = dd_250[yy,xx] #This is a very simplistic (i.e. wrong) 250 micron flux
        
        S_500 = S_500 #Convert to common units
        S_250 = S_250
        
        dust_T = 10.
        column = 1000.
        #kappa_500 = herschel_utilities.get_kappa(500)
        #dust_T = herschel_utilities.find_temperature(S_250/(18.1**2),S_500/(36.9**2))
        #column = herschel_utilities.find_column(dust_T,kappa_500,S_500*1000.,36.9,0.5) #1000 takes Jy/beam to mJy/beam
        #print(dust_T)
        #print(column)
        #print(S_500)
        info = {"dust_T":dust_T,"column":column,"xy":xy,"S_500":S_500}
        return(info)

    def get_spectra_at_higal_peak(self,line,herschel_info):
        """Obtain the spectrum for a given line at the specified location"""
        datapath = self.get_cube_path(line)
        d,h = pyfits.getdata(datapath,header=True)
        xy = herschel_info["xy"]
        self.get_velocity()
        cen_chan = int(((self.velocity*1e3)-h['CRVAL3'])/h['CDELT3']+h['CRPIX3'])
        #print(line)
        #print(self.velocity)
        #print(cen_chan)
        fullsize = 302
        min_chan = max(0,cen_chan-int(fullsize/2.))
        max_chan = min(4096,cen_chan+int(fullsize/2.))
        spectra1 = d[min_chan:max_chan,round(xy[1]),round(xy[0])]
        spectra2 = d[min_chan:max_chan,round(xy[1])-1,round(xy[0])]
        spectra3 = d[min_chan:max_chan,round(xy[1])+1,round(xy[0])]
        spectra4 = d[min_chan:max_chan,round(xy[1]),round(xy[0])-1]
        spectra5 = d[min_chan:max_chan,round(xy[1]),round(xy[0])+1]
        spectra = (spectra1 + spectra2 + spectra3 + spectra4 + spectra5)/5.

        Av = herschel_info["S_500"]
        #spectra = spectra/(Av/50.)
        
        fullspec = np.zeros(fullsize)
        if len(spectra) < fullsize:
            n = fullsize-len(spectra)
            #print(n)
            #print(len(spectra))
            #print(len(fullspec[n:]))
            if min_chan == 0:
                fullspec[n:] = spectra 
            if max_chan == 4096:
                pass
                #spectra = np.append(spectra,np.zeros(n))
        else: 
            fullspec = spectra
                
        return(fullspec)


    def do_cutout(self,undone_cutout):
        """Make a cutout"""
        success = True
        version = 2
        print("Making a cut-out for: "+undone_cutout)
        cutsize = 0.08

        #print(self.filenames)

        if undone_cutout.startswith("I"):
            #Spitzer/GLIMPSE
            channel = undone_cutout[2]
            mosaic = "self.glimpse_mosaic_"+channel
            try:
                reproject_map.do_reprojection(eval(mosaic),"GALACTIC","GALACTIC",size=cutsize,center=(self.apos,self.bpos),
                                                outfile=self.filenames[undone_cutout])
            except montage.status.MontageError:
            #    success=False
            #    print()
                return(False)
            os.remove(self.filenames[undone_cutout].replace('.fits','_area.fits'))
        elif undone_cutout.startswith("M"):
            #Spitzer/MIPSGAL
            try:
                #reproject_map.do_reprojection(self.mips_mosaic,"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                #                            outfile=self.filenames[undone_cutout])
                reproject_map.do_reprojection(self.mips_mosaic,"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                            outfile=self.filenames[undone_cutout])
            except:
                print("Failed to reproject MIPS!!!")
                return(False)
            os.remove(self.filenames[undone_cutout].replace('.fits','_area.fits'))
            d,h = pyfits.getdata(self.filenames[undone_cutout],header=True)
            d[np.where(d!=d)] = 2000.0
            pyfits.writeto(self.filenames[undone_cutout],d,h,clobber=True)
        elif undone_cutout.startswith("T"):
            #2MASS
            coord_string = str(self.ra)+","+str(self.dec)
            headname = self.filenames[undone_cutout].replace(".fits",".hdr")
            montage.mHdr(coord_string,cutsize,headname,system="galactic")
            band = undone_cutout[2] #This gets correct filter from name
            print("2MASS "+band+" "+self.filenames[undone_cutout]+" "+headname)
            montage.mExec("2MASS",band,output_image=self.filenames[undone_cutout],region_header=headname)
            os.remove(headname)
            
        elif undone_cutout.startswith("R"):
            maptype = undone_cutout[3:]
            print("Actually doing reproject")
            print(undone_cutout)
            print(maptype)
            mosaic1 = "self.RH_mos_"+maptype+"_try1"
            mosaic2 = "self.RH_mos_"+maptype+"_try2"
            
            try:
                reproject_map.do_reprojection(eval(mosaic1),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                            outfile=self.filenames[undone_cutout])
            except:
                reproject_map.do_reprojection(eval(mosaic2),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                outfile=self.filenames[undone_cutout])
            os.remove(self.filenames[undone_cutout].replace('.fits','_area.fits'))
            
        elif undone_cutout.startswith("A"):
            #ATLASGAL does NOT use existing mosaics.
            #Instead, we fetch from Magpis
            print("Doing Atlasgal")
            #print(self.filenames)
            #ATLASGAL
            image = Magpis.get_images(coordinates.Galactic(float(self.glon), float(self.glat),
                    unit=(u.deg,u.deg)), image_size=cutsize*u.deg, survey='atlasgal')
            print(image)
            fits.writeto(self.filenames[undone_cutout],
                         image[0].data,image[0].header,clobber=True)
                
                
                #reproject_map.do_reprojection(self.mips_mosaic,"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                #                            outfile=self.filenames[undone_cutout])
                #print(self.atlasgal_mos)
                #reproject_map.do_reprojection(self.atlasgal_mos,"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                #                            outfile=self.filenames[undone_cutout])
            #except IndexError:
            #    print("Failed to fetch data for ATLASGAL!!!")
            #    return(False)
            #os.remove(self.filenames[undone_cutout].replace('.fits','_area.fits'))
            
        elif undone_cutout.startswith("H"):
            #Herschel HiGal
            print("Trying to make a cut-out")
            channel = undone_cutout[1:4]
            try:
                os.mkdir("temp")
                mosaic = "self.H"+channel+"c_try1_file1"
                file1 = self.filenames[undone_cutout].replace("cutouts","temp").replace(".fits","1.fits")
                reproject_map.do_reprojection(eval(mosaic),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                outfile=file1,list_o_files=True,hdu=1)
                mosaic = "self.H"+channel+"c_try1_file2"
                file2 = self.filenames[undone_cutout].replace("cutouts","temp").replace(".fits","2.fits")
                reproject_map.do_reprojection(eval(mosaic),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                outfile=file2,list_o_files=True,hdu=1)
                #montage.mosaic(self.storage_dir+"temp/",self.storage_dir+"temp2/")
                #Now we assume these files are the exact same size and projection
                dd1,hh1 = pyfits.getdata(file1,header=True)
                dd2,hh2 = pyfits.getdata(file2,header=True)
                try:
                    dd = np.average(np.dstack((dd1,dd2)),axis=2)
                except ValueError:
                    dd = dd1 #Remove this hack!!
                pyfits.writeto(self.filenames[undone_cutout],dd,hh1,clobber=True)
            except (montage.status.MontageError, IndexError):
                print("Failure!!")
                try:
                    mosaic = "self.H"+channel+"c_try2_file1"
                    file1 = self.filenames[undone_cutout].replace("cutouts","temp").replace(".fits","1.fits")
                    reproject_map.do_reprojection(eval(mosaic),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                    outfile=file1,list_o_files=True,hdu=1)
                    mosaic = "self.H"+channel+"c_try2_file2"
                    file2 = self.filenames[undone_cutout].replace("cutouts","temp").replace(".fits","2.fits")
                    reproject_map.do_reprojection(eval(mosaic),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                    outfile=file2,list_o_files=True,hdu=1)
                    #montage.mosaic(self.storage_dir+"temp/",self.storage_dir+"temp2/")
                    #Now we assume these files are the exact same size and projection
                    dd1,hh1 = pyfits.getdata(file1,header=True)
                    dd2,hh2 = pyfits.getdata(file2,header=True)
                    try:
                        dd = np.average(np.dstack((dd1,dd2)),axis=2)
                    except ValueError:
                        dd = dd1 #Remove this hack!!
                    pyfits.writeto(self.filenames[undone_cutout],dd,hh1,clobber=True)
                except montage.status.MontageError:
                    success = False
            finally:
                try:
                    shutil.rmtree("temp")
                except OSError:
                    pass
            
            #try to reproject try 1 500 
            #If this fails, try 2
            #Then coadd/median combine the re_projected images
            
        #Look-up source mosaic / get 2MASS image / get multiple files
            #Lead character in cutout can tell us the basics
        #Extract using montage
        #Do any necessary post-processing.
            #MIPS fix NANs to 2000
            #Herschel rationalize header?
        #Save to the appropriate place
        print("Saving to: "+self.filenames[undone_cutout])
        if success:
            self.report_success(undone_cutout,version)

    def read_auxiliary_log(self):
        """ Read info from the log """
        f = open(self.auxiliary_log,"r")
        for i,line in enumerate(f):
            #print(i)
            if i == 0:
                keys = line.split()
                #print(keys)
            elif i == 1:
                current_versions = line.split()
                #print(current_versions)
                #print(keys)
                pair_list = []
                for key,vers in zip(keys,current_versions):
                    pair_list.append((key,vers))
                current_best_version = dict(pair_list)
            else:
                data_items = line.split()
                #print(data_items)
                if data_items[0] == self.name:
                    pair_list = []
                    for key,vers in zip(keys,data_items):
                        pair_list.append((key,vers))
                    data_version = dict(pair_list)
        f.close()
        return(current_best_version,data_version)

    def report_success(self,entry,version):
        """ Record a successfully updated file in the log"""
        os.rename(self.auxiliary_log,self.auxiliary_log+"~")
        source = open(self.auxiliary_log+"~","r")
        destination = open(self.auxiliary_log,"w")
        for i,line in enumerate(source):
            if i == 0:
                keys = line.split()
                destination.write(line)
            else:
                data_items = line.split()
                if data_items[0] == self.name:
                    pair_list = []
                    for key,vers in zip(keys,data_items):
                        pair_list.append((key,vers))
                    data_version = dict(pair_list)
                    data_version[entry] = version
                    for key in keys:
                        destination.write(str(data_version[key])+"\t")
                    destination.write("\n")
                else:
                    destination.write(line)
        source.close()
        destination.close()

    def identify_herschel_mosaics(self):
        """Identify the Herscehl mosaics for each source"""
        herschel_path = malt.herschel_path
        list_all_mosaics = glob.glob(os.path.join(herschel_path,"Glon*"))
        list_all_mosaics.sort()
        #print(list_all_mosaics)
        coords = [float(os.path.basename(a)[4:7]) for a in list_all_mosaics]
        coords = list(set(coords))
        coords.sort()
        #print(coords)
        
        def find_le(a, x):
            'Find rightmost value less than or equal to x'
            i = bisect.bisect_right(a, x)
            if i:
                return a[i-1]
            raise ValueError
            
        def find_ge(a, x):
            'Find leftmost item greater than or equal to x'
            i = bisect.bisect_left(a, x)
            if i != len(a):
                return a[i]
            raise ValueError   
        self.h_try1 = find_le(coords,self.glon)
        try:
            self.h_try2 = find_ge(coords,self.glon)
        except ValueError:
            self.h_try2 = self.h_try1
        #print(self.h_try1)
        oddly_named = [303,305,323]
        if (self.h_try1 in oddly_named) or (self.h_try2 in oddly_named):
            suffix = ["one","two"]
        else:
            suffix = ["Nom","Orth"]
            
        #Lots of redundant code here
        #500 micron
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try1)).zfill(3)+"-"+suffix[0],"level2","PLW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H500c_try1_file1 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try1)).zfill(3)+"-"+suffix[1],"level2","PLW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H500c_try1_file2 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try2)).zfill(3)+"-"+suffix[0],"level2","PLW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H500c_try2_file1 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try2)).zfill(3)+"-"+suffix[1],"level2","PLW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H500c_try2_file2 = currfile
        
        #350 micron
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try1)).zfill(3)+"-"+suffix[0],"level2","PMW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H350c_try1_file1 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try1)).zfill(3)+"-"+suffix[1],"level2","PMW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H350c_try1_file2 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try2)).zfill(3)+"-"+suffix[0],"level2","PMW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H350c_try2_file1 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try2)).zfill(3)+"-"+suffix[1],"level2","PMW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H350c_try2_file2 = currfile
        
        #250 micron
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try1)).zfill(3)+"-"+suffix[0],"level2","PSW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H250c_try1_file1 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try1)).zfill(3)+"-"+suffix[1],"level2","PSW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H250c_try1_file2 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try2)).zfill(3)+"-"+suffix[0],"level2","PSW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H250c_try2_file1 = currfile
        dir1 = os.path.join(herschel_path,"Glon"+str(int(self.h_try2)).zfill(3)+"-"+suffix[1],"level2","PSW")
        currfile = glob.glob(os.path.join(dir1,"*.fits"))
        self.H250c_try2_file2 = currfile
        
        
        print(self.H500c_try1_file1)
        print(self.H500c_try1_file2)
        print(self.H500c_try2_file1)
        print(self.H500c_try2_file2)
        print(self.H350c_try1_file1)
        print(self.H350c_try1_file2)
        print(self.H350c_try2_file1)
        print(self.H350c_try2_file2)
        print(self.H250c_try1_file1)
        print(self.H250c_try1_file2)
        print(self.H250c_try2_file1)
        print(self.H250c_try2_file2)

    def identify_atlasgal_mosaics(self):
        """ Select the ATLASGAL mosaic for each source."""
        agal_path = malt.agal_path
        if (self.glon) > -1 and (self.glon < 10):
            self.atlasgal_mos = agal_path+"ATLASGAL_main.fits"
        elif (self.glon) > 330 and (self.glon < 360):
            self.atlasgal_mos = agal_path+"ATLASGAL_main.fits"
        elif (self.glon > 10) and (self.glon < 30):
            #Find the correct positive file
            if (self.glon < 12):
                self.atlasgal_mos = agal_path+"ATLASGAL.10.5.fits"
            elif (self.glon < 15) and (self.glon > 12):
                self.atlasgal_mos = agal_path+"ATLASGAL.13.5.fits"
            elif (self.glon < 18) and (self.glon > 15):
                self.atlasgal_mos = agal_path+"ATLASGAL.16.5.fits"
            elif (self.glon < 21) and (self.glon > 18):
                self.atlasgal_mos = agal_path+"ATLASGAL.19.5.fits"
        elif (self.glon > 307.5) and (self.glon < 322):
            self.atlasgal_mos = agal_path+"ATLASGAL-m38-m52.fits"
        elif (self.glon > 299.6) and (self.glon < 305):
            self.atlasgal_mos = agal_path+"ATLASGAL-m55-m60.fits"
        else:
            self.atlasgal_mos = None
        print("Start attempt to identify Agal")
        print(self.glon)
        print(self.atlasgal_mos)
        print("End attempt to identify Agal")
        
    def identify_herschel_result_mosaics_new(self):
        """Select the Herschel result (T and logN) mosaic
        
        The new versionf of T and logN from Yannet are
        now just one large mosaic covering every(?) source.
        Keep the same structure for consistency.
        
        """
        h_path = malt.h_path
        
        self.RH_mos_T_try1 = os.path.join(h_path,"Tdust_ait.fits")
        self.RH_mos_T_try2 = os.path.join(h_path,"Tdust_ait.fits")
        self.RH_mos_N_try1 = os.path.join(h_path,"N_ait.fits")
        self.RH_mos_N_try2 = os.path.join(h_path,"N_ait.fits")
        
        
    def identify_spitzer_mosaics(self):
        """ Selects the GLIMPSE/MIPSGAL mosaic for each source. """
        glimpse_path = malt.glimpse_path
        mips_path = malt.mips_path
        glimpseII_path = malt.glimpseII_path
        mipsgalII_path = malt.mipsgalII_path
        
        if self.glon > 180:
            self.l = self.glon-360.
        else:
            self.l = self.glon
        self.b = self.glat

        if self.l < 0:
            l_lookup = 360. + self.l
        else:
            l_lookup = self.l
        frame = str(int(round(l_lookup/3)*3)).zfill(3)
        #print(frame)
        if frame == "360":
            frame = "000"
        if abs(self.l) < 4.5:
            glimpse_mosaic_base = glimpseII_path+"v3.5_3.1x4.5_mosaics/"+"GLM_"+str(frame)+"00+0000_mosaic_I"
            self.glimpse_mosaic_1 = glimpse_mosaic_base+"1.fits"
            self.glimpse_mosaic_2 = glimpse_mosaic_base+"2.fits"
            self.glimpse_mosaic_3 = glimpse_mosaic_base+"3.fits"
            self.glimpse_mosaic_4 = glimpse_mosaic_base+"4.fits"
        elif abs(self.l) < 7.5:
            glimpse_mosaic_base = glimpseII_path+"v3.5_3.1x3.45_mosaics/"+"GLM_"+str(frame)+"00+0000_mosaic_I"
            self.glimpse_mosaic_1 = glimpse_mosaic_base+"1.fits"
            self.glimpse_mosaic_2 = glimpse_mosaic_base+"2.fits"
            self.glimpse_mosaic_3 = glimpse_mosaic_base+"3.fits"
            self.glimpse_mosaic_4 = glimpse_mosaic_base+"4.fits"
        elif abs(self.l) < 10.5: #We're in GLIMPSE II mosaics
            glimpse_mosaic_base = glimpseII_path+"v3.5_3.1x2.4_mosaics/"+"GLM_"+str(frame)+"00+0000_mosaic_I"
            self.glimpse_mosaic_1 = glimpse_mosaic_base+"1.fits"
            self.glimpse_mosaic_2 = glimpse_mosaic_base+"2.fits"
            self.glimpse_mosaic_3 = glimpse_mosaic_base+"3.fits"
            self.glimpse_mosaic_4 = glimpse_mosaic_base+"4.fits"
        else:
            glimpse_mosaic_base = glimpse_path+"GLM_"+str(frame)+"00+0000_mosaic_I"
            glimpse4_mosaic = glimpse_path+"mosaicv3_ch4_"+str(frame)+".fits"
            self.glimpse_mosaic_1 = glimpse_mosaic_base+"1.fits"
            self.glimpse_mosaic_2 = glimpse_mosaic_base+"2.fits"
            self.glimpse_mosaic_3 = glimpse_mosaic_base+"3.fits"
            self.glimpse_mosaic_4 = glimpse4_mosaic
        #print(self.glimpse_mosaic_4)
        mframe = str(int(round(l_lookup))).zfill(3)
        if mframe == "360":
            mframe = "000"

        if self.b < 0:
            bname = 'n'
        else:
            bname = 'p'
        #Normal
        if abs(self.b) < 1.1:
            mips_base = mipsgalII_path+"MG"+mframe+"0"+bname+'005_024.fits'
        #Hack for high points
        else:
            mips_base = mipsgalII_path+"MG"+mframe+"0"+bname+'015_024.fits'

        self.mips_mosaic = mips_base

        self.IR1c = self.name+"_I1.fits"
        self.IR2c = self.name+"_I2.fits"
        self.IR3c = self.name+"_I3.fits"
        self.IR4c = self.name+"_I4.fits"
        self.M24c = self.name+"_M1.fits"

        print(self.mips_mosaic)