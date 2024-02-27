"""
Adapted for code repository on 2024-02-26

description: Data fro Supplementary Figures 8 and 9 - 
             print out of hazard frequencies and intensities of present and future
             results of SI Tables are frequency and intensity changes, e.g. frequency
             future - frequency present.

@author: simonameiler
"""

import sys
import numpy as np

#Load Climada modules
from climada.util.constants import SYSTEM_DIR # loads default directory paths for data
from climada.hazard import TropCyclone

def main(region, period):
    
    
    ###########################################################################
    ############### A: define constants, functions, paths #####################
    ###########################################################################
    
    # define paths
    haz_dir = SYSTEM_DIR/"hazard/future"
    
    res = 300
    region = str(region) # AP, IO, SH, WP
    period = str(period) # fut1. fut2
    

    # load hazard
    # make list
    h1_min, h1_max = (1, 6) # model
    h2_min, h2_max = (1, 3) # scenario
    h3_min, h3_max = (1, 2) # cat
    h4_min, h4_max = (1, 2) # wind
    
    model_key = {1: 'CESM2', 
                 2: 'CNRM-CM6-1', 
                 3: 'EC-Earth3', 
                 4: 'IPSL-CM6A-LR', 
                 5: 'MIROC6', 
                # 6: 'MPI-ESM1-2-HR', 
                 6: 'UKESM1-0-LL'}

    ssp_haz_key = {1: 'ssp245', 
                   2: 'ssp370',
                   3: 'ssp585'}
    
    cat_key = {1: 'CRH',
               2: 'SD'}
    
    wind_model_key = {1: 'H08',
                      2: 'ER11'}
    
    # present climate
    tc_haz_base_dict = {}
    for h1 in range(h1_min, h1_max+1):
        for h2 in range(h2_min, h2_max+1):
            for h3 in range(h3_min, h3_max+1):
                for h4 in range(h4_min, h4_max+1):
                    haz_base_str = f"TC_{region}_0300as_CHAZ_{model_key[h1]}_base_{ssp_haz_key[h2]}_80ens_{cat_key[h3]}_{wind_model_key[h4]}.hdf5"
                    tc_haz_base = TropCyclone.from_hdf5(haz_dir.joinpath(haz_base_str))
                    ev_filt = np.where(tc_haz_base.intensity.sum(axis=1)>0)[0].tolist()
                    haz_base = tc_haz_base.select(event_id = ev_filt)
                    haz_base.check()
                    tc_haz_base_dict[str(model_key[h1])+'_'+str(ssp_haz_key[h2])+'_'+str(cat_key[h3])+'_'+str(wind_model_key[h4])] = haz_base
    
    # future climate
    tc_haz_fut_dict = {}
    for h1 in range(h1_min, h1_max+1):
        for h2 in range(h2_min, h2_max+1):
            for h3 in range(h3_min, h3_max+1):
                for h4 in range(h4_min, h4_max+1):
                    haz_fut_str = f"TC_{region}_0300as_CHAZ_{model_key[h1]}_{period}_{ssp_haz_key[h2]}_80ens_{cat_key[h3]}_{wind_model_key[h4]}.hdf5"
                    tc_haz_fut = TropCyclone.from_hdf5(haz_dir.joinpath(haz_fut_str))
                    ev_filt = np.where(tc_haz_fut.intensity.sum(axis=1)>0)[0].tolist()
                    haz_fut = tc_haz_fut.select(event_id = ev_filt)
                    haz_fut.check()
                    tc_haz_fut_dict[str(model_key[h1])+'_'+str(ssp_haz_key[h2])+'_'+str(cat_key[h3])+'_'+str(wind_model_key[h4])] = haz_fut
    
    # store frequency information - present
    for h, haz in tc_haz_base_dict.items():
        print(h, haz.frequency.sum())    
    
    # store frequency information - future
    for h, haz in tc_haz_fut_dict.items():
        print(h, haz.frequency.sum())

    # store intensity information - present
    for h, haz in tc_haz_base_dict.items():  
        print(h, np.max(haz.intensity, axis=1).mean())
        
    # store intensity information - future
    for h, haz in tc_haz_fut_dict.items():
        print(h, np.max(haz.intensity, axis=1).mean())
    
if __name__ == "__main__":
    main(*sys.argv[1:])
