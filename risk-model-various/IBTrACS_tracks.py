"""
Adapted for code repository on 2024-02-26

description: Load TC tracks from IBTrACS record and generate synthetic TC set.
            
@author: simonameiler
"""

import sys
import copy as cp

# import CLIMADA modules:
from climada.hazard import TCTracks
from climada.util.constants import SYSTEM_DIR

############################################################################
    
def main(year):

    year = int(year)
    IB_tracks_dir = SYSTEM_DIR/"tracks"/"IBTrACS"
    IB_synth_dir = SYSTEM_DIR/"tracks"/"IBTrACS_p"/f"{year}"
    
    # Load IBTrACS and generate probabilistic set
    tracks = TCTracks.from_ibtracs_netcdf(year_range=(year,year+1))
    tracks_IB = TCTracks()
    for i in range(0,6):
        filterdict = {'category': i}
        tracks_IB.data.extend(tracks.subset(filterdict).data)
    # post processing, increase time steps for smoother wind field:
    tracks_IB.equal_timestep(time_step_h=1., land_params=False)
    tracks_IB.write_netcdf(IB_tracks_dir)
    
    # generate probabilistic tracks
    IB_tracks_synth = cp.deepcopy(tracks_IB)
    IB_tracks_synth.data = [x for x in IB_tracks_synth.data if x.time.size > 1]
    IB_tracks_synth.calc_perturbed_trajectories(nb_synth_tracks=24)
    IB_tracks_synth.write_netcdf(IB_synth_dir)

if __name__ == "__main__":
    main(*sys.argv[1:])