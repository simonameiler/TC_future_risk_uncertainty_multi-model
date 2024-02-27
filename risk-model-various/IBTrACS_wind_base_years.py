"""
Adapted for code repository on 2024-02-26

description: Load TC tracks from probabilistic IBTrACS and calculate the 2D windfield
            after Holland (2008) and Emanuel and Rotunno (2011).
            
@author: simonameiler
"""

import sys

# import CLIMADA modules:
from climada.hazard import Centroids, TCTracks, TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################

def main(windmodel, year):

    windmodel = str(windmodel)
    year = int(year)
        
    res = 300
    
    IB_synth_dir = SYSTEM_DIR/"tracks"/"IBTrACS_prob"/f"{year}"
    haz_dir = SYSTEM_DIR.joinpath('hazard','future')
    haz_str = f"TC_global_0{res}as_IBTrACS_prob_present_{windmodel}_{year}.hdf5"
    
    cent_str = SYSTEM_DIR.joinpath("earth_centroids_0300as_global.hdf5")

    # load global, probabilistic IBTrACS
    tracks = TCTracks.from_netcdf(IB_synth_dir)
    
    # load centroids from this source
    cent = Centroids.from_hdf5(cent_str)

    tc_hazard = TropCyclone.from_tracks(tracks, centroids=cent, model=windmodel)
    tc_hazard.write_hdf5(haz_dir.joinpath(haz_str))


if __name__ == "__main__":
    main(*sys.argv[1:])