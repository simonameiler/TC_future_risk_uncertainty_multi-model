"""
Adapted for code repository on 2024-02-26

description: Generate future TC hazard sets for probabilistic IBTrACS for various
            RCPs, regions, periods and both windmodels
            
@author: simonameiler
"""

import sys
from scipy.sparse import csr_matrix

# import CLIMADA modules:
from climada.hazard import Centroids, TropCyclone
from climada.util.constants import SYSTEM_DIR

############################################################################

def main(region, scenario, future, windmodel):
    
    res = 300
    reg = str(region)
    rcp = int(scenario)
    future = int(future)
    windmodel = str(windmodel)
    
    haz_dir = SYSTEM_DIR.joinpath('hazard','future')
    haz_str = f"TC_{reg}_0{res}as_IBTrACS_prob_present_{windmodel}.hdf5"
    
    tc_hazard = TropCyclone.from_hdf5(haz_dir.joinpath(haz_str))

    tc_hazard_cc = tc_hazard.apply_climate_scenario_knu(
                ref_year=future, rcp_scenario=rcp)
    haz_fut_str = f"TC_{reg}_0{res}as_IBTrACS_prob_{rcp}_{future}_{windmodel}.hdf5"
    tc_hazard_cc.write_hdf5(haz_dir.joinpath(haz_fut_str))

if __name__ == "__main__":
    main(*sys.argv[1:])
