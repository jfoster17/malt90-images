import os
import reproject_map
import astLib.astCoords as ac
import pyfits
import numpy as np
import aplpy
import numdisplay.zscale as zscale
import montage
import Malt90Source

class Malt90SourceBadMips(Malt90Source.Malt90Source):
    """A class to store information and do stuff for about a Malt90Source with bad MIPS"""

    def set_label_text(self,undone_image):   
        if undone_image.startswith("M"):
            label1 = r"$\mathrm{GLIMPSE}\/3.6\mu \mathrm{m}$"
            label2 = r"$\mathrm{GLIMPSE}\/8.0\mu \mathrm{m}$"
            label3 = r"$\mathrm{WISE}\/22\mu \mathrm{m}$"
        elif undone_image.startswith("G"):
            label1 = r"$\mathrm{GLIMPSE}\/3.6\mu \mathrm{m}$"
            label2 = r"$\mathrm{GLIMPSE}\/4.5\mu \mathrm{m}$"
            label3 = r"$\mathrm{GLIMPSE}\/8.0\mu \mathrm{m}$"
        return(label1,label2,label3)

    def identify_spitzer_mosaics(self):
        """ Selects the GLIMPSE/MIPSGAL image to display. """
        glimpse_path = "/Volumes/Mako3/glimpsev3/"
        mips_path = "/Volumes/Mako3/mipsgal/"
        glimpseII_path = "/Volumes/Data/GLIMPSEII/"
        mipsgalII_path = "/Volumes/Data/MIPSGALII/"
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

        if self.b < 0:
            bname = 'n'
        else:
            bname = 'p'
        
        all_wise = os.listdir("WISE_Files")
        search_string = "ra"+str(self.apos_ra)[0:6]
        print(search_string)
        for wise_file in all_wise:
            if search_string in wise_file:
                if "w4" in wise_file:
                    self.mips_mosaic = "WISE_Files/"+wise_file
                if "w3" in wise_file:
                    self.glimpse_mosaic_4 = "WISE_Files/"+wise_file
                if "w2" in wise_file:
                    self.glimpse_mosaic_2 = "WISE_Files/"+wise_file
                    self.glimpse_mosaic_3 = "WISE_Files/"+wise_file
                if "w1" in wise_file:
                    self.glimpse_mosaic_1 = "WISE_Files/"+wise_file
        print(self.mips_mosaic)
        #Very special case
        #if (self.id == "AG6645") or (self.id == "AG6647") or (self.id == "AG6673") or (self.id == "AG6678") or (self.id == "AG6679"):
        #    for wise_file in all_wise:
        #        if search_string in wise_file:
        #            if "w4" in wise_file:
        #                self.mips_mosaic = "WISE_Files/"+wise_file
        #            if "w3" in wise_file:
        #                self.glimpse_mosaic_4 = "WISE_Files/"+wise_file
        #            if "w2" in wise_file:
        #                self.glimpse_mosaic_2 = "WISE_Files/"+wise_file
        #                self.glimpse_mosaic_3 = "WISE_Files/"+wise_file
        #            if "w1" in wise_file:
        #                self.glimpse_mosaic_1 = "WISE_Files/"+wise_file
        #print("IRAC 1")
        #print(self.IR1c)    
        self.IR1c = self.name+"_I1.fits"
        self.IR2c = self.name+"_I2.fits"
        self.IR3c = self.name+"_I3.fits"
        self.IR4c = self.name+"_I4.fits"
        self.M24c = self.name+"_M1.fits"


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
            print(eval(mosaic))
            reproject_map.do_reprojection(eval(mosaic),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                            outfile=self.filenames[undone_cutout])
            
            try:
                reproject_map.do_reprojection(eval(mosaic),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                outfile=self.filenames[undone_cutout])
            except montage.status.MontageError:
                try:
                    reproject_map.do_reprojection(eval(mosaic),"GALACTIC","GALACTIC",size=cutsize,center=(self.apos,self.bpos),
                                                    outfile=self.filenames[undone_cutout])
                except montage.status.MontageError:
                    return(False)
                    
            #    success=False
            #    print()
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
        elif undone_cutout.startswith("R"):
            maptype = undone_cutout[3:]
            mosaic1 = "self.RH_mos_"+maptype+"_try1"
            mosaic2 = "self.RH_mos_"+maptype+"_try2"
            
            try:
                reproject_map.do_reprojection(eval(mosaic1),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                            outfile=self.filenames[undone_cutout])
            except montage.status.MontageError:
                reproject_map.do_reprojection(eval(mosaic2),"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                                outfile=self.filenames[undone_cutout])
            os.remove(self.filenames[undone_cutout].replace('.fits','_area.fits'))
        
        elif undone_cutout.startswith("T"):
            #2MASS
            coord_string = str(self.ra)+","+str(self.dec)
            headname = self.filenames[undone_cutout].replace(".fits",".hdr")
            montage.mHdr(coord_string,cutsize,headname,system="galactic")
            band = undone_cutout[2] #This gets correct filter from name
            print("2MASS "+band+" "+self.filenames[undone_cutout]+" "+headname)
            montage.mExec("2MASS",band,output_image=self.filenames[undone_cutout],region_header=headname)
            os.remove(headname)
        elif undone_cutout.startswith("A"):
            #print(self.filenames)
            #ATLASGAL
            try:
                #reproject_map.do_reprojection(self.mips_mosaic,"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                #                            outfile=self.filenames[undone_cutout])
                print(self.atlasgal_mos)
                reproject_map.do_reprojection(self.atlasgal_mos,"J2000","GALACTIC",size=cutsize,center=(self.apos_ra,self.bpos_dec),
                                            outfile=self.filenames[undone_cutout])
            except IndexError:
                print("Failed to reproject ATLASGAL!!!")
                return(False)
            os.remove(self.filenames[undone_cutout].replace('.fits','_area.fits'))
            
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

