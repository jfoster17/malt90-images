import os
import aplpy
import numpy as np
import pylab
import math
import pyfits
import idl_stats
import collections
import copy
import interpolate_image
from astropy.table import Table, Row

class MALT90Source():
    """Base class for MALT90 sources. Can be inherited to deal with e.g. Pilot, year1, year2"""

    def __init__(self,atlasgal_id,goodname,velocity,cat_row,source_lons,source_lats,source_ids):
        self.name = goodname
        #self.glon = float(goodname[1:8])
        #self.glat = float(goodname[8:15])
        
        self.glon = float(cat_row['ag_long'])
        self.glat = float(cat_row['ag_lat'])
        self.agal_name = cat_row['source']
        self.malt90_name = cat_row['malt90_map']
        
        self.agal_id = atlasgal_id
        self.aname = self.agal_id
        self.source_lons = source_lons
        self.source_lats = source_lats
        self.source_ids = source_ids
        try:
            self.velocities = [float(velocity)]
        except ValueError:
            self.velocities = [float(velocity.split(',')[0]),float(velocity.split(',')[1])]
        self.cat_row = cat_row

        self.max_intensities = {'hcop':cat_row['hcop_ii'],
                           'hnc':cat_row['hnc_ii'],
                           'n2hp':cat_row['n2hp_ii'],
                           'hcn':cat_row['hcn_ii'],
                           'h13cop':cat_row['h13cop_ii'],
                           'hn13c':cat_row['hn13c_ii'],
                           '13cs':cat_row['13cs_ii'],
                           '13c34s':cat_row['13c34s_ii'],
                           'hc13ccn':cat_row['hc13ccn_ii'],
                           'hnco404':cat_row['hnco404_ii'],
                           'ch3cn':cat_row['ch3cn_ii'],
                           'hc3n':cat_row['hc3n_ii'],
                           'hnco413':cat_row['hnco413_ii'],
                           'c2h':cat_row['c2h_ii'],
                           'sio':cat_row['sio_ii'],
                           'h41a':cat_row['h41a_ii']}


        self.data_dir = '/Volumes/Screwdriver/malt90'
        self.result_dir = '/Volumes/Screwdriver/malt90/results/'
        self.momkeys = ["mom0","err0","mom1","err1","mom2","err2","npix","snr0"]
        self.peak = collections.defaultdict(dict)
        self.peak_pos = {}
        lineint  = {"hcop":0.,"hcn":0.,"hnc":0.,"n2hp":0.,
                                "13c34s":0.,"13cs":0.,"c2h":0.,"ch3cn":0.,
                                "h13cn":0.,"h13cop":0.,"h41alpha":0.,"hc3n":0.,
                                "hc13ccn":0.,"hnco404":0.,"hnco413":0.,"sio":0.}

        linepos = {"hcop":[],"hcn":[],"hnc":[],"n2hp":[],
                                "13c34s":[],"13cs":[],"c2h":[],"ch3cn":[],
                                "h13cn":[],"h13cop":[],"h41alpha":[],"hc3n":[],
                                "hc13ccn":[],"hnco404":[],"hnco413":[],"sio":[]}
        #spectra = {"hcop":np.empty([512,2]),"hcn":np.empty([512,2]),
         #                       "hnc":np.empty([512,2]),"n2hp":np.empty([512,2])}

       # self.center_spectra = copy.deepcopy(spectra)


    def get_path_to_cube(self,line):
        """Get path to a cube"""
        filename = os.path.join(self.data_dir,"gridzilla",line,self.name+"_"+line+"_"+"MEAN.fits")
        return(filename)


    def do_image(self,gc,fname,levels):
        gc.tick_labels.set_yformat(format="d.dd")
        gc.tick_labels.set_xformat(format="ddd.dd")
        interpolated_fname = interpolate_image.interpolate(fname)
        try:
            gc.show_contour(interpolated_fname,cmap='gist_yarg',filled=True,levels=levels,extend='both',smooth=3)
        except AttributeError:
            pass
        gc.show_contour(interpolated_fname,colors='black',levels=levels,extend='both',smooth=3)
        #gc.clabel(levels,inline=1)
        gc.ticks.set_color('black')


    def get_display_range(self,line=None, mom='mom0',secondvel=False,fixed=False,vmin=None,vmax=None):
        """Return an appropriate estimate for displaying data"""
        if mom == 'mom0':
            mmin = 0
            mmax = self.max_intensities[line]
            #print(self.agal_id)
            #print(mmax)
            if mmax == -999:
                mmax = 0
            return(mmin,mmax)

    def make_moment_images(self,linelist,tag,format='.png',mom='mom0',secondvel=False,outdir="Images"):
        latex_dict = {'hcop':r'HCO$^+$','n2hp':r'N$_2$H$^+$','hnc':r'HNC','hcn':r'HCN','c2h':r'C$_2$H',
                                'sio':r'SiO','h13cop':r'H$^{13}$CO$^{+}$','h13cn':r'H$^{13}$CN','hnco413':r'HNCO 4$_{1,3}$',
                                'hnco404':r'HNCO 4$_{0,4}$','hc3n':r'HC$_{3}$N','13c34s':r'$^{13}$C$^{34}$S',
                                'hc13ccn':r'HC$^{13}$CCN','13cs':r'$^{13}$CS','h41alpha':r'H41$\alpha$','ch3cn':r'CH$_{3}$CN',
                                'hn13c':r'HN$^{13}$C','h41a':r'H41$\alpha$'}
        n_images = len(linelist)
        n = math.sqrt(n_images)
        if n_images == 4:
            mainline = True
        else:
            mainline = False
        xoff = 0.15
        yoff = 0.10
        xcsize = 0.05
        
        vmin = str(round(self.velocities[0]-8.5,1))
        vmax = str(round(self.velocities[0]+8.5,1))

        fullimage = pylab.figure(figsize = (9,9))
        xsize = (1-(xoff+3*xcsize))/(n)
        ysize = xsize

        for i,line in enumerate(linelist):
            #print(line)
            fname = self.get_path(line,mom,secondvel)
            #print(fname)
            
            mmin, mmax = self.get_display_range(line=line,mom=mom,secondvel=secondvel,fixed=False)

            if mmin == mmax:
                levels = np.linspace(1000,1010,num=9) #Do not show contours
            else:
                midpoint = (mmin+mmax)/2
                lev1 = np.linspace(mmin,midpoint,num=6,endpoint=False)
                lev2 = np.linspace(midpoint,mmax,num=3)
                levels = np.array([lev1[0],lev1[1],lev1[2],lev1[3],lev1[4],lev1[5],lev2[0],lev2[1],lev2[2]])
                #levels = np.array([1.0,2.0,3.0,4.5,6.0,7.5,10,13,16])
            
            gc = aplpy.FITSFigure(fname,figure=fullimage,subplot=[xoff+(i%n)*xsize,yoff+(i//n)*ysize,xsize,ysize],auto_refresh=False)

            self.do_image(gc,fname,levels)
            if i%n != 0:
                gc.axis_labels.hide_y()
                gc.tick_labels.hide_y()
            if i//n != 0:
                gc.axis_labels.hide_x()
                gc.tick_labels.hide_x()
            gc.axis_labels.set_xtext("")
            gc.axis_labels.set_ytext("")
            if mainline:
                gc.add_label(0.03,0.93,latex_dict[line],fontsize=16,relative=True,color='black',backgroundcolor=(1,1,1,0.7),horizontalalignment='left',alpha=0.9)
            else:
                gc.add_label(0.05,0.85,latex_dict[line],fontsize=16,relative=True,color='black',backgroundcolor=(1,1,1,0.7),horizontalalignment='left',alpha=0.9)
                gc.add_label(0.05,0.70,str(round(mmax,1))+" K km/s",size=6,relative=True,color='black',backgroundcolor=(1,1,1,0.7),horizontalalignment='left',alpha=0.9)
        
            if not mainline:
                gc.ticks.set_xspacing(0.03)
            if mainline:
                gc.ticks.set_xspacing(0.02)
            gc.show_markers(self.source_lons,self.source_lats,edgecolor='blue',facecolor='none',marker='+',s=20)
            gc.show_markers(self.glon,self.glat,edgecolor='red',facecolor='none',marker='+',s=20)

        #gc = aplpy.FITSFigure('ColorSwatch.fits',figure=fullimage,subplot=[0.90,0.1,0.05,0.35],auto_refresh=False)
        #gc.show_colorscale(cmap='gist_yarg',vmin=-1.5,vmax=8.5)
        #gc.axis_labels.hide_x()
        #gc.axis_labels.set_ytext('')
        #if mom == 'snr0':
        #    label_text = "SNR"
        #if mom == 'mom0':
        #    label_text = "K km/s"
        #if mom == 'mom1':
        #    label_text = "km/s"
        #if mom == 'mom2':
        #    label_text = "km/s"
            
        #gc.add_label(0.5,1.1,label_text,relative=True)

        #gc.tick_labels.hide()
        #gc.ticks.set_length(0)
        #for i in range(9):
        #    gc.add_label(-0.5,(i)*0.125,np.around(levels[i],1),relative=True)

        fullimage.text(0.15,0.87,self.agal_id+" : "+self.agal_name+" : "+self.malt90_name,fontsize = 14)#+"  Levels = "+str(np.around(levels,2))+" (K km/s)")
        fullimage.text(0.15,0.83,"Velocity Range: "+vmin+" to "+vmax+" km/s")
        pylab.figtext(0.5,0.03,"Galactic Longitude [degrees]",ha='center')
        pylab.figtext(0.03,0.5,"Galactic Latitude [degrees]",rotation='vertical',va='center')

        if format == ".png":
            filename = self.agal_id
        else:
            filename = self.agal_id.replace('.','_')

        if secondvel:
            fullimage.savefig(outdir+"/"+filename+"_"+mom+"_v2"+tag+format)
        else:
            fullimage.savefig(outdir+"/"+filename+"_"+mom+"_"+tag+format)
        pylab.close()

    def get_path(self,line='hcop',mmap='emap',secondvel = False):
        """Get path to a specific moment map"""
        if secondvel:
            filename = os.path.join(self.data_dir,"sources",self.name,self.name+"_v2_"+line+"_mommaps",self.name+"_v2_"+line+"_"+mmap+".fits")
        else:
            filename = os.path.join(self.data_dir,"mommaps",line,self.aname+"_"+line+"_mommaps",self.aname+"_"+line+"_"+mmap+".fits")

        return(filename)
