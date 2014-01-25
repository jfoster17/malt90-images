#This file stores variables for use
#throughout the MALT90 reduction pipeline
#Updated for remaking moment maps on Butch
#where not all the pipeline will work.

#Redundant with data_dir

base = '/Volumes/Scratch/malt90'

#reduce_malt
vnum = {"rename":"1.5","ldata":"1.6","gzilla":"1.6","arrange":"1.6","mommaps":"1.7"}
#sd = '/nfs/atapplic/malt/reduce/'
sd = '/epp/atapplic/malt/malt90-analysis-code/reduce/' #Seems to be new location
data_dir = base+

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
