"""
Adapted for code repository on 2024-02-26

description: Load TC tracks from the MIT model and calculate the 2D windfield
            after Holland (2008) and Emanuel and Rotunno (2011).

@author: simonameiler
"""

import sys
import numpy as np
from scipy.io import loadmat
import datetime as dt

# import CLIMADA modules:
from climada.hazard import Centroids, TCTracks, TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################

def main(region, model, scenario, wind_model):
    
    region = str(region)
    model = str(model)
    scenario = str(scenario)
    wind_model = str(wind_model)
        
    res = 300
    yrs_total = 20
    
    tracks_dir = SYSTEM_DIR.joinpath('tracks','Kerry','future')
    fname = tracks_dir.joinpath(f"Meiler_{region}_{model}_{scenario}.mat")
    haz_dir = SYSTEM_DIR.joinpath('hazard','future')
    haz_str = f"TC_{region}_0{res}as_MIT_{model}_{scenario}_{wind_model}.hdf5"
    
    cent_str = SYSTEM_DIR.joinpath("earth_centroids_0300as_global.hdf5")
    
    # Initiate MIT tracks
    def init_MIT_tracks(region):
        if region == 'SH':
            tracks_MIT = TCTracks.from_simulations_emanuel(fname, hemisphere='S')
        else:
            tracks_MIT = TCTracks.from_simulations_emanuel(fname, hemisphere='N')
        tracks_MIT.equal_timestep(time_step_h=1)
        return tracks_MIT
    
    # call functions
    tc_tracks = init_MIT_tracks(region)
    
    # load centroids from this source
    cent = Centroids.from_hdf5(cent_str)

    tc_hazard = TropCyclone.from_tracks(tc_tracks, centroids=cent, model=wind_model)
    
    # apply frequency correction according to the freq scalar provided with the
    # event sets
    fname = tracks_dir.joinpath(f"Meiler_{region}_{model}_{scenario}.mat")
    freq_year = loadmat(fname)['freqyear'][0].tolist()
    
    event_year = np.array([
        dt.datetime.fromordinal(d).year 
        for d in tc_hazard.date.astype(int)])

    for i, yr in enumerate(np.unique(event_year)):
        yr_mask = (event_year == yr)
        yr_event_count = yr_mask.sum()
        tc_hazard.frequency[yr_mask] = (
        np.ones(yr_event_count) * freq_year[i] /
        (yr_event_count * yrs_total)
        )    
    tc_hazard.write_hdf5(haz_dir.joinpath(haz_str))
    tc_hazard.check()

if __name__ == "__main__":
    main(*sys.argv[1:])