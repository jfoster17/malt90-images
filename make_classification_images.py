"""

6/5/2014: Attempts to do this for new images.

8/4/2014: Once more updating for new (hopefully final) catalog

11/19/2014: Yet again the catalog has changed.

"""


import os,sys
import malt90_catalog as mcat
import AgalSource
import AgalSourceWise
import numpy as np
import gc as garbagecollection
import malt_params as malt


#bad_mips = [780.0, 2584.0, 3325.0, 3329.0, 4354.0,
#    4357.0, 4361.0, 4362.0, 4363.0, 4378.0, 4379.0,
#    4381.0, 4767.0, 4769.0, 5041.0, 5043.0, 5045.0,
#    5047.0, 5050.0, 5257.0, 5505.0, 6640.0, 6645.0,
#    6647.0, 6650.0, 6652.0, 6657.0, 6662.0, 6669.0,
#    6670.0, 6673.0, 6674.0, 6678.0, 6679.0, 6680.0,
#    6691.0, 6695.0, 6706.0, 6711.0, 6712.0, 6715.0,
#    6738.0, 6758.0, 6776.0, 6802.0, 6812.0, 6874.0,
#    7017.0, 7023.0, 7024.0, 7035.0, 7048.0, 7095.0,
#    7150.0, 7153.0, 7201.0]

bad_mips = [4737.0, 
4926.0,
5609.0,
6112.0,
6113.0,
6350.0,
6591.0,
6649.0,
7047.0,
7188.0,
7589.0,
8158.0,
8160.0,
8161.0,
8405.0,
8499.0,
7226.0,
7227.0,
7229.0,
7230.0,
7232.0,
7233.0,
7234.0,
7235.0,
7236.0,
]

bad_mips = []

   #780.0, 2584.0, 3325.0, 3329.0, 4354.0,
   #4357.0, 4361.0, 4362.0, 4363.0, 4378.0, 4379.0,
   #4381.0, 4767.0, 4769.0, 5041.0, 5043.0, 5045.0,
   #5047.0, 5050.0, 5257.0, 5505.0, 6640.0, 6645.0,
   #6647.0, 6650.0, 6652.0, 6657.0, 6662.0, 6669.0,
   #6670.0, 6673.0, 6674.0, 6678.0, 6679.0, 6680.0,
   #6691.0, 6695.0, 6706.0, 6711.0, 6712.0, 6715.0,
   #6738.0, 6758.0, 6776.0, 6802.0, 6812.0, 6874.0,
   #7017.0, 7023.0, 7024.0, 7035.0, 7048.0, 7095.0,
   #7150.0, 7153.0, 7201.0, 4926.0]

#6645 -> 


def make_wise_lookup_table(t):
    import atpy
    ras = []
    decs = []
    
    #bad_hand = [4926]
    bad_hand = bad_mips
    
    for i,source in enumerate(t['ag_id']):
        num = float(source[2:])
        sourcename = t['malt90_map_filename'][i]
        a_lon = t['ag_long'][i]
        a_lat = t['ag_lat'][i]
        cat_row = t[i] 
        #if (a_lat > 1.) or (a_lat < -1) or (num in bad_hand) :
        if (num in bad_hand) :
            print(num)
            agal = AgalSource.AgalSource(source,cat_row,
                                         None,None,None)
            ras.append(agal.malt90source.apos_ra)
            decs.append(agal.malt90source.bpos_dec)
    tbl = atpy.Table()
    tbl.add_column('ra',np.array(ras),unit='degrees')
    tbl.add_column('dec',np.array(decs),unit='degrees')
    tbl.write("wise-table-for-lookup-more.tbl")
                                    

def main():
    #Read in catalog
    t = mcat.read_latest("../results/malt90catalog.cat")
    #print(t['agid'])
    #return(0)
    #t = mcat.read("malt90_lineinfo_sm.cat")
    #make_wise_lookup_table(t)
    for i,source in enumerate(t['agid']):
        num = float(source[2:])
        #print(num)
        #This is to do just the WISE sources
        #if (num in bad_mips) :

        #This deals with a couple oddball problem cases
        
        #badnums = [6681,6811,7017]
        #if (num in badnums):
        sourcename = t['malt90_map'][i]
        #m90_p = t['malt90_p'][i]
        
        #bad_agals.txt lits 7188 and onwards
        
        #malt90s = t[t['malt90_map_filename'] != 'notObs']
        #print(len(malt90s))
        
        if (num > 0): #and (num < 1950):
            print("="*10+source+"="*10)
            #vel = t['consensus_velocity'][i]
            sourcename = t['malt90_map'][i]
            print(sourcename)
            
            a_lon = t['ag_long'][i]
            a_lat = t['ag_lat'][i]
            distances = np.array(np.sqrt((t['ag_long']-a_lon)**2
                                 +(t['ag_lat']-a_lat)**2))
            min_distance = 3/60.
            nearby = t[(distances < min_distance)]
            source_lons = nearby['ag_long']
            source_lats = nearby['ag_lat']
            source_ids = nearby['agid']
            cat_row = t[i] 
            
            aim = malt.image_dir+source+"_Aim.png"
            gim = malt.image_dir+source+"_Gim.png"
            mim = malt.image_dir+source+"_Mim.png"
            
            if ((not os.path.isfile(aim)) or (not os.path.isfile(gim)) or (not os.path.isfile(mim))):
                do_images = True                
            else:
                do_images = True
                
            
            if (num in bad_mips) and do_images: 
                agal = AgalSourceWise.AgalSourceWise(source,cat_row,
                                             source_lons,source_lats,source_ids)
                agal.make_classification_images()
                agal.make_herschel_result_images()
                del agal
            else:
                #try:
                if do_images:
                    #try:
                    agal = AgalSource.AgalSource(source,cat_row,
                                             source_lons,source_lats,source_ids)
                    agal.make_classification_images()
                    agal.make_herschel_result_images()
                    #except :
                    #    f1 = open('newest_bad_agals.txt','a')
                    #    output = str(num)+"\n"
                    #    f1.write(output)
                    #    f1.close()
                        
                else:
                    pass
                #agal.make_herschel_result_images()
                #except:
                #    f1 = open('bad_agals_herschel.txt','a')
                #    output = str(num)+"\n"
                #    f1.write(output)
                #    f1.close()
                    
                #del agal
                
            print("="*26)
            garbagecollection.collect()

    

if __name__ == '__main__':
    main()
