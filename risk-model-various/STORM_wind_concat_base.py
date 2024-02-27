"""
Adapted for code repository on 2024-02-26

description: Load STORM windfields from genesis basins, concatenate to study regions,
            apply frequency bias correction, save. Present climate.

@author: simonameiler
"""

import sys
import numpy as np
import copy as cp

# import CLIMADA modules:
from climada.hazard import TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################
#i_wind = ['H08', 'ER11']
# i_ens = range(10)
# i_basin = ['EP', 'NA', 'NI', 'SI', 'SP', 'WP']
def main(i_wind):

    i_wind = str(i_wind)
    
    haz_dir = SYSTEM_DIR/"hazard"/"STORM_present"
    
    regions = {'AP': ['EP', 'NA'],
               'IO': 'NI',
               'SH': ['SI', 'SP'],
               'WP': 'WP'}
    
    freq_corr_STORM = 1/1000
    
    # load all STORM hazard files and append to list
    for reg in regions.keys():
        for i_ens in range(10):
            if reg == 'IO':
                i_basin = regions[reg]
                tc_hazard = TropCyclone.from_hdf5(haz_dir/f"TC_{i_basin}_{i_ens}_0300as_STORM_{i_wind}.hdf5")
                tc_hazard.frequency = np.ones(tc_hazard.size)*freq_corr_STORM
                tc_hazard.write_hdf5(haz_dir/f"TC_{reg}_{i_ens}_0300as_STORM_{i_wind}.hdf5")
            
            else:
                reg_list = regions[reg]
                STORM_hazard = []
                for i_basin in reg_list:
                    tc_hazard = TropCyclone.from_hdf5(haz_dir/f"TC_{i_basin}_{i_ens}_0300as_STORM_{i_wind}.hdf5")
                    STORM_hazard.append(tc_hazard)
                        
                STORM_master = cp.deepcopy(STORM_hazard[0])
                for haz in range(1,len(STORM_hazard)):
                    STORM_master.append(STORM_hazard[haz])            
                
                STORM_master.frequency = np.ones(STORM_master.size)*freq_corr_STORM
                STORM_master.write_hdf5(haz_dir/f"TC_{reg}_{i_ens}_0300as_STORM_{i_wind}.hdf5")


if __name__ == "__main__":
    main(*sys.argv[1:]) 
    

