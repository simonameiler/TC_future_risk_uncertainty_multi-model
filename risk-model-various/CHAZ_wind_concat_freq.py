#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adapted for code repository on 2024-02-26

description: Load TC windfield subsets from CHAZ model, concatenate, apply frequency
            bias correction and save.
            
@author: simonameiler
"""

import sys
import os
import re
import numpy as np
import copy as cp
from pathlib import Path

# import CLIMADA modules:
from climada.hazard import TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################

def main(model, scenario, cat, wind, period):

    model = str(model) # CESM2, ...
    scenario = str(scenario) # ssp245, ssp370, ssp585
    cat = str(cat) # CRH, SD
    wind = str(wind) # H08, ER11
    period = str(period) # base, fut1, fut2
    
    haz_in = SYSTEM_DIR/"hazard"/"future"/"CHAZ2"
    haz_out = SYSTEM_DIR/"hazard"/"future"
    
    
    # prepare frequency bias correction
    yrly_freq_IB_lit = {'AP': 25.3,
                        'IO': 2.0,
                        'SH': 21.6,
                        'WP': 22.5,
                        'global': 71.4}
    
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
    
    # function for frequency bias correction of CHAZ tracks
    def freq_bias_corr(hazard, bsn, years):
        num_tracks = hazard.intensity.max(axis=1).getnnz(axis=0)[0]
        freq_IB = yrly_freq_IB_lit[bsn]
        cor = freq_IB/(num_tracks/years)
        freq_corr = cor/years
        hazard.frequency = np.ones(hazard.size)*freq_corr
        return hazard
    
    
    # load all CHAZ hazard files and append to list
    file_list = []
    file1 = f'TC_global_0300as_CHAZ_{model}_{period}_{scenario}_2ens00'
    file2 = f'_{cat}_{wind}_'
    pattern = re.compile(file1 + r'(\d+)' + file2 + r'(\d+)' + r'\.hdf5')
    
    for filename in os.listdir(haz_in):
        match = re.match(pattern, filename)
        if match:
            groups = match.groups()
            if len(groups) == 2:
                file_list.append(filename)

    # make hazard object from the list of files
    CHAZ_hazard = TropCyclone()
    for fl in file_list:
        tc_hazard = TropCyclone.from_hdf5(haz_in/fl)
        CHAZ_hazard.append(tc_hazard)
    
    #apply frequency bias correction and save
    CHAZ_hazard = freq_bias_corr(CHAZ_hazard, 'global', 1600)
    CHAZ_hazard.write_hdf5(haz_out.joinpath(f'TC_global_0300as_CHAZ_{model}_{period}_{scenario}_80ens_{cat}_{wind}.hdf5'))
    
    # call basin split function and frequency bias correction and save results
    for bsn in BASIN_BOUNDS:
        tc_haz_basin = TropCyclone()
        tc_haz_basin = basin_split_haz(CHAZ_hazard, bsn)
        # frequency correction 
        tc_haz_basin = freq_bias_corr(tc_haz_basin, bsn, 1600) #years 8000 because 400 ensembles * 20 years
        tc_haz_basin.write_hdf5(haz_out.joinpath(
            f"TC_{bsn}_0300as_CHAZ_{model}_{period}_{scenario}_80ens_{cat}_{wind}.hdf5")) 

if __name__ == "__main__":
    main(*sys.argv[1:])
