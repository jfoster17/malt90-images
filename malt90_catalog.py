"""
Just read in the malt90_lineinfo catalog

"""
from astropy.table import Table

def read_yeah_this_changed_again(infile):
    t = Table.read(infile,format="ascii.fixed_width",
                ##         id,s, L  lo la mi mn  mm  mp, mv,hco,hnc,n2h,hcn,h13 ,hn13,13cs,hc13,404, ch3 ,sio ,hc3n, c2h,413,13c34,h41
                col_starts=[0,7 ,30,36,46,60,69, 91,107,115,213,429,645,861,1077,1293,1509,1725,1941,2157,2373,2589,2805,3021,3237,3453],
                col_ends  =[6,29,31,47,57,67,90,106,114,125,222,438,654,870,1086,1302,1518,1734,1950,2166,2382,2598,2814,3030,3246,3462],
                header_start=None,
                data_start = 0,
                guess=False,
                names=["agid","source","L","ag_long","ag_lat","m90id","malt90_name",
                       "malt90_map","m90_p","velocity",
                       "hcop_ii","hnc_ii","n2hp_ii","hcn_ii",
                       "h13cop_ii","hn13c_ii","13cs_ii","hc13ccn_ii",
                       "hnco404_ii","ch3cn_ii","sio_ii","hc3n_ii",
                       "c2h_ii","hnco413_ii","13c34s_ii","h41a_ii"],
                )
    print(t[0])

    return(t)

def readlatest(infile):
    t = Table.read(infile,format="ascii.fixed_width",
                col_starts=[0,7 ,30,32,43,51,59, 81,103,110,209,421,633,845,1057,1269,1481,1693,1905,2117,2329,2541,2753,2956,3168,3389],
                col_ends  =[6,29,31,42,50,58,80,102,109,120,218,430,642,854,1066,1278,1490,1702,1914,2126,2338,2550,2762,2965,3177,3398],
                header_start=None,
                data_start = 0,
                guess=False,
                names=["agid","source","L","ag_long","ag_lat","m90id","malt90_name",
                       "malt90_map","m90_p","velocity",
                       "hcop_ii","hnc_ii","n2hp_ii","hcn_ii",
                       "h13cop_ii","hn13c_ii","13cs_ii","hc13ccn_ii",
                       "hnco404_ii","ch3cn_ii","sio_ii","hc3n_ii",
                       "c2h_ii","hnco413_ii","13c34s_ii","h41a_ii"],
                )
    #print(t[0:30])
    return(t)


def readnew(infile):
    t = Table.read(infile,format = "ascii",
         names = ('ag_id', 'ag_name', 'ag_long', 'ag_lat', 'unkown1', 'unknown2', 'unknown3', 'malt90_map_filename',
                 'unknown4', 'unknown5', 'malt90_pos', 'malt90_name','malt90_id')
                )
    return(t)


def read(infile):
    t = Table.read(infile,format = "ascii",
         names = ('ag_id', 'ag_name', 'ag_long', 'ag_lat', 'malt90_id', 'malt90_name',
                  'malt90_map_filename', 'malt90_pos', 'classification',
                  'classification_string', 'consensus_velocity',
                  'hcop_name', 'hcop_status', 'hcop_whichfit', 'hcop_snmom', 
                  'hcop_snfit', 'hcop_rms', 'hcop_fitprob', 'hcop_Tastar', 
                  'hcop_dTastar', 'hcop_Vlsr', 'hcop_dVlsr', 'hcop_fwhm', 
                  'hcop_dfwhm','hcop_ii', 'hcop_dii', 'hcop_detection',
                  'hnc_name', 'hnc_status', 'hnc_whichfit', 'hnc_snmom', 
                  'hnc_snfit', 'hnc_rms', 'hnc_fitprob', 'hnc_Tastar', 
                  'hnc_dTastar', 'hnc_Vlsr', 'hnc_dVlsr', 'hnc_fwhm', 
                  'hnc_dfwhm','hnc_ii', 'hnc_dii', 'hnc_detection',
                  'n2hp_name', 'n2hp_status', 'n2hp_whichfit', 'n2hp_snmom', 
                  'n2hp_snfit', 'n2hp_rms', 'n2hp_fitprob', 'n2hp_Tastar', 
                  'n2hp_dTastar', 'n2hp_Vlsr', 'n2hp_dVlsr','n2hp_fwhm', 
                  'n2hp_dfwhm','n2hp_ii', 'n2hp_dii', 'n2hp_detection',
                  'hcn_name', 'hcn_status', 'hcn_whichfit', 'hcn_snmom', 
                  'hcn_snfit', 'hcn_rms', 'hcn_fitprob', 'hcn_Tastar', 
                  'hcn_dTastar', 'hcn_Vlsr', 'hcn_dVlsr', 'hcn_fwhm', 
                  'hcn_dfwhm','hcn_ii', 'hcn_dii', 'hcn_detection',
                  'h13cop_name', 'h13cop_status', 'h13cop_whichfit', 
                  'h13cop_snmom', 'h13cop_snfit', 'h13cop_rms',
                  'h13cop_fitprob', 'h13cop_Tastar', 'h13cop_dTastar', 
                  'h13cop_Vlsr', 'h13cop_dVlsr', 'h13cop_fwhm', 'h13cop_dfwhm',
                  'h13cop_ii', 'h13cop_dii', 'h13cop_detection',
                  'hn13c_name', 'hn13c_status', 'hn13c_whichfit', 'hn13c_snmom', 
                  'hn13c_snfit', 'hn13c_rms', 'hn13c_fitprob', 'hn13c_Tastar', 
                  'hn13c_dTastar', 'hn13c_Vlsr', 'hn13c_dVlsr', 'hn13c_fwhm', 
                  'hn13c_dfwhm','hn13c_ii', 'hn13c_dii', 'hn13c_detection',
                  '13cs_name', '13cs_status', '13cs_whichfit', '13cs_snmom', 
                  '13cs_snfit', '13cs_rms', '13cs_fitprob', '13cs_Tastar', 
                  '13cs_dTastar', '13cs_Vlsr', '13cs_dVlsr',
                  '13cs_fwhm', '13cs_dfwhm','13cs_ii', '13cs_dii', '13cs_detection',
                  'hc13ccn_name', 'hc13ccn_status', 'hc13ccn_whichfit', 
                  'hc13ccn_snmom', 'hc13ccn_snfit', 'hc13ccn_rms',
                  'hc13ccn_fitprob', 'hc13ccn_Tastar', 'hc13ccn_dTastar', 
                  'hc13ccn_Vlsr', 'hc13ccn_dVlsr', 'hc13ccn_fwhm', 'hc13ccn_dfwhm',
                  'hc13ccn_ii', 'hc13ccn_dii', 'hc13ccn_detection',
                  'hnco404_name', 'hnco404_status', 'hnco404_whichfit', 
                  'hnco404_snmom', 'hnco404_snfit', 'hnco404_rms',
                  'hnco404_fitprob', 'hnco404_Tastar', 'hnco404_dTastar', 
                  'hnco404_Vlsr', 'hnco404_dVlsr', 'hnco404_fwhm', 
                  'hnco404_dfwhm','hnco404_ii', 'hnco404_dii', 'hnco404_detection',
                  'ch3cn_name', 'ch3cn_status', 'ch3cn_whichfit', 'ch3cn_snmom', 
                  'ch3cn_snfit', 'ch3cn_rms', 'ch3cn_fitprob', 'ch3cn_Tastar', 
                  'ch3cn_dTastar', 'ch3cn_Vlsr', 'ch3cn_dVlsr', 'ch3cn_fwhm', 
                  'ch3cn_dfwhm','ch3cn_ii', 'ch3cn_dii', 'ch3cn_detection',
                  'hc3n_name', 'hc3n_status', 'hc3n_whichfit', 'hc3n_snmom', 
                  'hc3n_snfit', 'hc3n_rms', 'hc3n_fitprob', 'hc3n_Tastar', 
                  'hc3n_dTastar', 'hc3n_Vlsr', 'hc3n_dVlsr',
                  'hc3n_fwhm', 'hc3n_dfwhm','hc3n_ii', 'hc3n_dii', 'hc3n_detection',
                  'hnco413_name', 'hnco413_status', 'hnco413_whichfit', 
                  'hnco413_snmom', 'hnco413_snfit', 'hnco413_rms',
                  'hnco413_fitprob', 'hnco413_Tastar', 'hnco413_dTastar', 
                  'hnco413_Vlsr', 'hnco413_dVlsr', 'hnco413_fwhm', 'hnco413_dfwhm',
                  'hnco413_ii', 'hnco413_dii', 'hnco413_detection',
                  'c2h_name', 'c2h_status', 'c2h_whichfit', 'c2h_snmom', 
                  'c2h_snfit','c2h_rms', 'c2h_fitprob', 'c2h_Tastar', 'c2h_dTastar', 
                  'c2h_Vlsr', 'c2h_dVlsr', 'c2h_fwhm', 'c2h_dfwhm','c2h_ii', 
                  'c2h_dii', 'c2h_detection',
                  '13c34s_name', '13c34s_status', '13c34s_whichfit', '13c34s_snmom', 
                  '13c34s_snfit', '13c34s_rms', '13c34s_fitprob', '13c34s_Tastar', 
                  '13c34s_dTastar', '13c34s_Vlsr', '13c34s_dVlsr', '13c34s_fwhm', 
                  '13c34s_dfwhm','13c34s_ii', '13c34s_dii', '13c34s_detection',
                  'sio_name', 'sio_status', 'sio_whichfit', 'sio_snmom', 
                  'sio_snfit','sio_rms', 'sio_fitprob', 'sio_Tastar', 'sio_dTastar', 
                  'sio_Vlsr', 'sio_dVlsr', 'sio_fwhm', 'sio_dfwhm','sio_ii', 
                  'sio_dii', 'sio_detection', 'h41a_name', 'h41a_status', 
                  'h41a_whichfit', 'h41a_snmom', 'h41a_snfit', 'h41a_rms',
                  'h41a_fitprob', 'h41a_Tastar', 'h41a_dTastar', 'h41a_Vlsr', 
                  'h41a_dVlsr',
                  'h41a_fwhm', 'h41a_dfwhm','h41a_ii', 'h41a_dii', 'h41a_detection'
                      ))
    return(t)
