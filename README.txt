# Ye Olde Guide to Making MALT90 Images #

For historical reasons, the scripts to generate MALT90 images are currently 
split into two separate tasks.

Task 1: Classification Images (Spitzer/WISE and Herschl dust temp/column)
Task 2: Moment Map creation and imaging

The top-level scripts are all called make_* and all have some help text 
available by running, for instance:
 
    make_agal_moments -h

These scripts generally work by calling methods in other objects. There 
are two main kinds of objects:
MALT90MomentSource.py -- a fairly simple class that just makes moment 
                         map images
MALT90Source.py       -- the base class for doing complicated stuff. 
                         There is a bunch of legacy code in here for 
                         grander plans that never happened, which 
                         makes it a little complicated. The main 
                         complication is that the code to make FITS 
                         cut-outs and combine these into three-color
                         cubes and images is very general to try and
                         work with lots of auxiliary data sources. Also,
                         to be efficient, there's a whole logging setup
                         to try and avoid duplication of expensive steps 
                         when making image. None of this should concern 
                         you, but it's all in there anyway.
                         
                         Malt90SourceBadMips, AgalSource, and 
                         AgalSourceWISE are all light wrappers around
                         this class with a few things changed/added.

## Requirements ##

The Montage program (http://montage.ipac.caltech.edu/index.html)

The following python libraries

* numpy
* matplotlib
* scipy
* pyfits
* astropy
* aplpy
* montage-wrapper
* astroquery
* astLib

## Location ##

Everything is currently set up in 

/DATA/MALT_1/MALT90/malt90_images/

with the following structure:
gridzilla/           #Symbolic link to the data cubes
sources/             #Symbolic link to the data organized by source
herschel-atlasgal/   #The Herschel dust temperature and column maps
mommaps/             #Moment maps generated for each velocity component
results/             #Contains both data and results; important parts:
        malt90catalog.cat          #current catalog
        aux_data.log               #log to track image generation
        images/                    #Spitzer/Classiciation and Herschel images
        mom0/NewLines              #Moment 0 pretty images
        spectra_images/            #Currently not used
        WISE_Files/                #WISE data to use where Spitzer is bad
scripts/              #The scripts to do all the image generation

The scripts can also be found on my github account:
https://github.com/jfoster17/malt90-images

## Prepwork ##

There are two main files in the scripts/ directory that you may need to 
edit when using a new version of the catalog. 

The most important is malt90_catalog.py, which has some ugly code 
for parsing the current catalog into an astopy table. Due to limitations of 
the astropy table code, it's currently hard-coded to read with fixed-format 
column numbers. The basic problem is that the current catalog has repeating 
column names (SNR for every line, for instance), and this causes problems with 
the astropy table read. A more flexible parser would be a worthwhile thing, I 
just didn't bother. Define a new reader-function if you need to and point  
the read_latest() function to it.

The other configuration file is malt_params.py. This has a lot of file 
locations hand-coded and is imported by many other scripts. Unless you need 
to go back and regenerate very early steps in the process, this file 
should be properly configured to work on Draco as is.


## Task 1: Classification/Herschel Images ##

### Updating Files ###

There is a whole long complicated procedure to create the necessary 
FITS cutouts of the Spitzer files from the raw mosaics. This is 
tied into a system for tracking these steps, which is important because 
they can be slow. HOWEVER, you should not have to do anything with 
this code sub-system because all the files have already been made for 
each MALT90 map (and each ATLASGAL/velocity source by necessity 
points to an existing MALT90 map). That is, there are a whole bunch 
of FITS cut-outs for the MALT90 map centered at G000.006+00.156, 
but it doesn't matter how many velocity components or ATLASGAL sources 
you decide are inside this map. The code looks up which MALT90 map 
is relevant, grabs those FITS cut-outs, and only names the output 
files after the ATLASGAL/catalog source. 

### Image Generation ###
#### make_classification_images.py ####

For eatch catalog source, we create an AgalSource (which is just 
a thin wrapper around a Malt90Source) that knows how to find the 
appropriate Spitzer images for the MALT90 map, as well as the 
ATLASGAL image for that map and the Herschel Tdust/Ncol data. These 
data files are then displayed with appropriate annotations and saved 
in results/images/ as PNG files.

On my machine, this script sometimes runs into memory issues despite 
my attempts at making sure that the garbage collection is working 
correctly. Possibly this is due to matplotlib objects not 
getting released/destroyed properly? The command line options allow 
you to specify a range of ATLASGAL ID numbers, so you could do them 
in batches of <1000 to circumvent this problem if it crops up on 
DRACO (I have not extensively tested this).

There is some old cruft code in make_classification_images.py 
related to figuring out which MALT90 sources had problematic or 
nonexistant Spitzer images. Since you should never have to deal 
with this, you don't have to worry about it.

## Task 2: Moment Map Creation and Imaging ##

### Moment Map Creation ###
#### make_agal_moments.py ####

We make moment maps for each entry in the source catalog, since the
velocity range over which we integrate is listed in the catalog.

The main script for this is make_agal_moments.py. Without any command line 
line options, this will run through the current version of the catalog and 
generate moment maps for each entry, saving them in mommaps/ in a directory 
that looks like this:
mommaps/
	13c34s/
		AG001_13c34s_mommaps/
		.../
	13cs/
	.../

### Moment Map Imaging ###
#### make_moment_images.py #### 

The zeroth moment maps are imaged for all the MALT90 transitions using the 
catalog measured integrated intensities to set the scale (and to suppress 
the image when it would just be noise). You really want to run 
make_agal_moments.py after changing the catalog before running this. To be 
safe, I also recommend deleting everything inside 
    results/mommaps/NewLines/
before running the script. 