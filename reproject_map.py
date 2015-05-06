#!/usr/bin/env python
# encoding: utf-8
"""
Re-project a map from equitorial to galactic coordinates.

Optionally trim to make it square.
Optrionally trim to a fixed size.
"""

try:
    import montage
except ImportError:
    import montage_wrapper as montage
import os
import pyfits
import astLib.astCoords as ac
import astLib.astWCS as aw
import subprocess

def main():
    pass

def do_reprojection(filename,input_sys,output_sys,size=None,center=None,outfile=None,list_o_files=False,hdu=0):
    system_lookup= {"GALACTIC":"galactic","J2000":"equatorial"}
    filename_lookup = {"GALACTIC":"GAL","J2000":"EQ"}

    if list_o_files:
        filename = filename[0]
    d,h = pyfits.getdata(filename,hdu,header=True)
    if not size:
        size = get_size(h)
        
    #Center coords always in J2000
    if not center:
        x,y = get_center_position(h,input_sys,"J2000")
    else:
        x,y = ac.convertCoords(input_sys,"J2000",center[0],center[1],2000.)
    pix_size = get_pixel_size(h)
    #print(x,y)
    print(pix_size)
    coord_string = str(x)+" "+str(y)
    #print(coord_string)
    #Note that older versions on Montage do not work because of normal split
    #dividing the coordinates. Fix is to use shlex.split instead.
    montage.mHdr(coord_string,size,"Test.hdr",system=system_lookup[output_sys],pix_size=pix_size*3600)
    #print(yo)
    montage.mSubimage(filename,"Test.fits",x,y,size*2,hdu=hdu)
    if not outfile: #Try a reasonable guess for output file
        outfile = filename.replace(".fits","_"+filename_lookup[output_sys]+".fits")
    montage.mProject("Test.fits",outfile,"Test.hdr")
    #Yay, this works.

def get_center_position(h,input_sys,output_sys):
    """Get position at center of map, convert to other system."""
    WCS = aw.WCS(h,mode="pyfits")
    xcen = int(h["NAXIS1"]/2.)
    ycen = int(h["NAXIS2"]/2.)

    original_coords = WCS.pix2wcs(xcen,ycen)
    modified_coords = ac.convertCoords(input_sys,output_sys,original_coords[0],original_coords[1],2000.)
    modified_coords = original_coords
    return(modified_coords[0],modified_coords[1])

def get_size(h):
    """Establish a reasonable size for the output image."""
    guess_one = abs(h["NAXIS1"]*h["CDELT1"])
    WCS = aw.WCS(h,mode="pyfits")
    guess_two = WCS.getFullSizeSkyDeg()
    #print(guess_one)
    #print(guess_two)
    return(guess_one)

def get_pixel_size(h):
    try:
        pixel_size = h["CDELT1"]
    except KeyError:
        pixel_size = h["CD1_1"]
        
    return(abs(pixel_size))

def make_header():
    pass

if __name__ == '__main__':
    main()
