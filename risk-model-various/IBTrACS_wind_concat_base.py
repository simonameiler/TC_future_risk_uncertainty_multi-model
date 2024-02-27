"""
Adapted for code repository on 2024-02-26

description: Load single year TC hazard sets from probabilistic IBTrACS, concatenate, 
            split by basin and save.
            
@author: simonameiler
"""

import sys
import numpy as np
import copy as cp

# import CLIMADA modules:
from climada.hazard import TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################
#windmodel = ['H08', 'ER11']

def main(windmodel):

    windmodel = str(windmodel)
    
    res = 300
    
    haz_dir = SYSTEM_DIR/"hazard"/"future"
    
    # boundaries of (sub-)basins (lonmin, lonmax, latmin, latmax)
    BASIN_BOUNDS = {
        # North Atlantic/Eastern Pacific Basin
        'AP': [-180.0, 0.0, 0.0, 65.0],
    
        # Indian Ocean Basin
        'IO': [30.0, 100.0, 0.0, 40.0],
    
        # Southern Hemisphere Basin
        'SH': [-180.0, 180.0, -60.0, 0.0],
    
        # Western Pacific Basin
        'WP': [100.0, 180.0, 0.0, 65.0],
    }
    
    reg_id = {'AP': 5000, 'IO': 5001, 'SH': 5002, 'WP': 5003}
    
    # function to split hazard set up per basin
    def basin_split_haz(hazard, basin):
        """ Split CHAZ global hazard up into ocean basins of choice """
        tc_haz_split = TropCyclone()
        # get basin bounds
        x_min, x_max, y_min, y_max = BASIN_BOUNDS[str(basin)]
        basin_idx = (hazard.centroids.lat > y_min) & (
                     hazard.centroids.lat < y_max) & (
                     hazard.centroids.lon > x_min) & (
                     hazard.centroids.lon < x_max)
        hazard.centroids.region_id[basin_idx] = reg_id[basin]
        tc_haz_split = hazard.select(reg_id=reg_id[basin]) 
        return tc_haz_split

    # create tropcyclone object for 1980 and append the subsequent years
    ib_haz = TropCyclone.from_hdf5(f"TC_global_0{res}as_IBTrACS_prob_present_{windmodel}_1980.hdf5")
    for year in range(1981, 2021):
        haz_str = f"TC_global_0{res}as_IBTrACS_prob_present_{windmodel}_{year}.hdf5"
        tc_haz = TropCyclone.from_hdf5(haz_dir/haz_str)
        ib_haz.append(tc_haz)
        
    ib_haz.write_hdf5(haz_dir.joinpath(f'TC_global_0{res}as_IBTrACS_prob_present_{windmodel}'))
    
    # call basin split function and save results
    for bsn in BASIN_BOUNDS:
        tc_haz_basin = TropCyclone()
        tc_haz_basin = basin_split_haz(ib_haz, bsn)
        tc_haz_basin.write_hdf5(haz_dir.joinpath(
            f"TC_{bsn}_0{res}as_IBTrACS_prob_present_{windmodel}.hdf5"))

if __name__ == "__main__":
    main(*sys.argv[1:]) 
    

