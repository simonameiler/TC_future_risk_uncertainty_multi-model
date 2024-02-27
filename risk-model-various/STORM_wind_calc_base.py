"""
Adapted for code repository on 2024-02-26

description: Load TC tracks from the STORM model and calculate the 2D windfield
            after Holland (2008) and Emanuel and Rotunno (ER11). Present climate model simulations.

@author: simonameiler
"""

import sys
import os
import numpy as np

# import CLIMADA modules:
from climada.hazard import Centroids, TCTracks, TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################
# i_ens = range(10)
# i_basin = ['EP', 'NA', 'NI', 'SI', 'SP', 'WP']
def main(i_basin, i_ens, windmodel):
    
    i_basin = str(i_basin)
    windmodel = str(windmodel)
    
    storm_dir = SYSTEM_DIR.joinpath('tracks','STORM', 'present')
    haz_dir = SYSTEM_DIR/"hazard"/"STORM_present"
    haz_str = f"TC_{i_basin}_{i_ens}_0300as_STORM_{windmodel}.hdf5"
    
    cent_str = SYSTEM_DIR.joinpath("earth_centroids_0300as_global.hdf5")
    
    def init_STORM_tracks(i_basin, i_ens):
        """ Load STORM tracks for the present climate."""
        fname = f"STORM_DATA_IBTRACS_{i_basin}_1000_YEARS_{i_ens}.txt"
        tracks_STORM = TCTracks.from_simulations_storm(os.path.join(storm_dir, fname))
        tracks_STORM.equal_timestep(time_step_h=1.)
        return tracks_STORM
    
    # call functions
    tc_tracks = TCTracks()
    tc_tracks = init_STORM_tracks(i_basin, i_ens)
    
    # load centroids from this source
    cent = Centroids.from_hdf5(cent_str)

    tc_hazard = TropCyclone.from_tracks(tc_tracks, centroids=cent, model=windmodel)
    tc_hazard.write_hdf5(haz_dir.joinpath(haz_str))
    tc_hazard.check()

if __name__ == "__main__":
    main(*sys.argv[1:])

