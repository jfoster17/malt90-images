#This file stores variables for use
#throughout the MALT90 reduction pipeline
#Updated for remaking moment maps on Butch
#where not all the pipeline will work.


lines = ["hnco413","c2h","sio","h41a",
         "hc13ccn","hnco404","ch3cn","hc3n",
         "h13cop","hn13c","13cs","13c34s",
         "hcop","hnc","n2hp","hcn",
         ]

#Redundant with data_dir

base = '/Volumes/Screwdriver/malt90'

#reduce_malt
vnum = {"rename":"1.5","ldata":"1.6","gzilla":"1.6","arrange":"1.6","mommaps":"1.7"}
#sd = '/nfs/atapplic/malt/reduce/'
sd = '/epp/atapplic/malt/malt90-analysis-code/reduce/' #Seems to be new location
data_dir = base

#ReduceLog
log_location = base+'reduction_log.txt'
lock    = base+'lock.txt'
source_dir = base+'raw_data/' #This is rather badly named
#classification_file = sd+'MALT90_classifications_2012.txt'

#preprocess
rename_dir = '/DATA/MALT_1/MALT90/data/renamed/'
#data_dir = '/DATA/MALT_1/MALT90/raw_data/'

#cal
cal_dir = base+"cal/"
ver_dir = base+"results/verification/"

#These are paths used by do_by_hand.py
byhand_rename_dir = base+'data/byhand/renamed/'
byhand_data_dir   = base+'data/byhand/'

#ATLASGAL information        
agal_path = "/Volumes/Data1/Herschel-ATLASGAL/"

#Herschel information
h_path = "/Volumes/Data1/Herschel-ATLASGAL/"

#Spitzer information
glimpse_path = "/Volumes/Mako3/glimpsev3/"
mips_path = "/Volumes/Mako3/mipsgal/"
glimpseII_path = "/Volumes/Data1/GLIMPSEII/"
mipsgalII_path = "/Volumes/Data1/MIPSGALII/"


#OLD Herschel information
#Uncomment and update if we need to remake Herschel images from raw files.
#herschel_path = "/Volumes/Mako3/higal_public_reg/"

