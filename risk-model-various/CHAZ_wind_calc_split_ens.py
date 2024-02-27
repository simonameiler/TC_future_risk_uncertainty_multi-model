#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adapted for code repository on 2024-02-26

description: Load TC tracks from the CHAZ model and calculate the 2D windfield
            after Holland (2008) and Emanuel and Rotunno (2011).
"""

import sys
import copy
import numpy as np
from pathos.pools import ProcessPool as Pool

# import CLIMADA modules:
from climada.hazard import Centroids, TCTracks, TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################

def main(i_file, model, scenario, cat, wind, period):

    i_file = int(i_file) # 0..9
    model = str(model) # CESM2, ...
    scenario = str(scenario) # ssp245, ssp370, ssp585
    cat = str(cat) # CRH, SD
    wind = str(wind) # H08, ER11
    period = str(period) # base, fut1, fut2
    
    year_range_dict = {'base': (1995, 2014),
                       'fut1': (2041, 2061),
                       'fut2': (2081, 2100)
                       }
    year_range = year_range_dict[period]
    ens_nums = np.arange(0, 40, 5).tolist()
    
    chaz_dir = SYSTEM_DIR.joinpath('tracks','CHAZ', 'CMIP6', scenario)
    cent_str = SYSTEM_DIR.joinpath("earth_centroids_0300as_global.hdf5")
    haz_dir = SYSTEM_DIR/"hazard"/"future"/"CHAZ2"
    
    # Initiate CHAZ tracks
    def init_CHAZ_tracks_ens(model, i_file, cat, year_range, ens_nums):
        fname = chaz_dir.joinpath(f"{model}_Global_2100_2ens00{i_file}_{cat}_compressed.nc")
        tracks_CHAZ = TCTracks.from_simulations_chaz(fname, year_range=year_range, ensemble_nums=ens_nums)
        #tracks_CHAZ.equal_timestep(time_step_h=1)
        return tracks_CHAZ
    
    # call functions
    tc_tracks = init_CHAZ_tracks_ens(model, i_file, cat, year_range, ens_nums)

    # load centroids from this source
    cent = Centroids.from_hdf5(cent_str)
    cent_tracks = cent.select(extent=tc_tracks.get_extent(5))

    pool = Pool()
    k = 1000
    for n in range(0, tc_tracks.size, k):
        tracks = copy.deepcopy(tc_tracks)
        tracks.data = tracks.data[n:n+k]
        tracks.equal_timestep(time_step_h=1.)
        tc = TropCyclone.from_tracks(tracks, centroids=cent_tracks, pool=pool)
        haz_str = f"TC_global_0300as_CHAZ_{model}_{period}_{scenario}_2ens00{i_file}_{cat}_{wind}_{n}.hdf5"
        tc.write_hdf5(haz_dir.joinpath(haz_str))
    pool.close()
    pool.join()

if __name__ == "__main__":
    main(*sys.argv[1:])
