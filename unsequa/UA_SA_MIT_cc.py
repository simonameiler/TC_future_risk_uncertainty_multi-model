"""
Adapted for code repository on 2024-02-26

description: Uncertainty and sensitivity analysis MIT - climate change contribution,
             exposure fixed at present-day baseline

@author: simonameiler
"""

import sys
import scipy as sp
import numpy as np
import logging

#Load Climada modules
from climada.util.constants import SYSTEM_DIR # loads default directory paths for data
from climada.hazard import Hazard

import climada.util.coordinates as u_coord
from climada.entity import Exposures
from climada.engine.unsequa import InputVar, CalcDeltaImpact
from climada.entity.impact_funcs.trop_cyclone import ImpfSetTropCyclone

def main(region, period, N_samples):
    
    LOGGER = logging.getLogger(__name__)
    
    ###########################################################################
    ############### A: define constants, functions, paths #####################
    ###########################################################################
    
    # define paths
    haz_dir = SYSTEM_DIR/"hazard/future"
    unsequa_dir = SYSTEM_DIR/"unsequa"
    
    res = 300
    ref_year = 2005
    N_samples = int(N_samples)
    region = str(region) # AP, IO, SH, WP
    period = str(period) # 20thcal (1995-2014), cal (2041-2060), _2cal (2081-2100)
    
    # translate period in the hazard object naming (K. Emanuel) to ref_year for
    # exposure scaling
    year_from_per = {'20thcal': 2005,
                     'cal': 2050,
                     '_2cal': 2090}
    
    ###########################################################################
    ########## B: load and define hazard, exposure, impf_sets #################
    ###########################################################################
    
    # load hazard
    # make list
    h1_min, h1_max = (1, 9)
    h2_min, h2_max = (1, 3)
    h3_min, h3_max = (1, 2)
    
    model_key = {1: 'cesm2', 
                 2: 'cnrm6', 
                 3: 'ecearth6', 
                 4: 'fgoals', 
                 5: 'ipsl6', 
                 6: 'miroc6', 
                 7: 'mpi6', 
                 8: 'mri6', 
                 9: 'ukmo6'}

    ssp_haz_key = {1: 'ssp245', 
                   2: 'ssp370',
                   3: 'ssp585'}
    
    wind_model_key = {1: '',
                      2: 'ER11_'}
    
    # present climate
    tc_haz_base_dict = {}
    for h1 in range(h1_min, h1_max+1):
        for h3 in range(h3_min, h3_max+1):
            haz_base_str = f"TC_{region}_0{res}as_{wind_model_key[h3]}MIT_{model_key[h1]}_20thcal.hdf5"
            tc_haz_base = Hazard.from_hdf5(haz_dir.joinpath(haz_base_str))
            tc_haz_base.check()
            tc_haz_base_dict[str(wind_model_key[h3])+str(model_key[h1])] = tc_haz_base
    
    # future climate
    tc_haz_fut_dict = {}
    for h1 in range(h1_min, h1_max+1):
        for h2 in range(h2_min, h2_max+1):
            for h3 in range(h3_min, h3_max+1):
                haz_fut_str = f"TC_{region}_0{res}as_{wind_model_key[h3]}MIT_{model_key[h1]}_{ssp_haz_key[h2]}{period}.hdf5"
                tc_haz_fut = Hazard.from_hdf5(haz_dir.joinpath(haz_fut_str))
                tc_haz_fut.check()
                tc_haz_fut_dict[str(wind_model_key[h3])+str(model_key[h1])+'_'+str(ssp_haz_key[h2])] = tc_haz_fut

        
    # load exposure
    # present
    e3_min, e3_max = (1, 9)
    
    mn_key = {1: [0.5, 0.5],
              2: [0.5, 1.0],
              3: [0.5, 1.5], 
              4: [1.0, 0.5],
              5: [1.0, 1.0],
              6: [1.0, 1.5],
              7: [1.5, 0.5],
              8: [1.5, 1.0],
              9: [1.5, 1.5]}
    
    # set exposure region_id to code regions for impact function mapping
    iso3n_per_region = impf_id_per_region = ImpfSetTropCyclone.get_countries_per_region()[2]
    code_regions = {'NA1': 1, 'NA2': 2, 'NI': 3, 'OC': 4, 'SI': 5, 'WP1': 6, \
                    'WP2': 7, 'WP3': 8, 'WP4': 9, 'ROW': 10}
    
    exp_base_dict = {}
    for e3 in range(e3_min, e3_max+1):
        [m, n] = mn_key[e3]
        ent_str = f"litpop_0300as_{ref_year}_{region}_{m}-{n}.hdf5"
        exp_base = Exposures.from_hdf5(SYSTEM_DIR.joinpath(ent_str))
        exp_base.assign_centroids(tc_haz_base)
        exp_base.value_unit = 'USD'
        
        # match exposure with correspoding impact function
        for calibration_region in impf_id_per_region:
            for country_iso3n in iso3n_per_region[calibration_region]:
                exp_base.gdf.loc[exp_base.gdf.region_id==country_iso3n, 'impf_TC'] = code_regions[calibration_region]
    
        # get iso3alpha codes from exposure region_ids
        exp_natids = np.unique(exp_base.gdf.region_id).tolist()
        #exp_iso = u_coord.country_to_iso(exp_natids, representation="alpha3")
    
        # add iso_code to exposure gdf
        for natid in exp_natids:
            exp_base.gdf.loc[exp_base.gdf.region_id==natid, 'iso_code'] = \
            u_coord.country_to_iso(natid, representation="alpha3")
    
        exp_base_dict[str(mn_key[e3])] = exp_base
    
    
    ###########################################################################
    ############## C: define input variables and parameters ###################
    ###########################################################################
    
    # hazard
    # present    
    def haz_base_func(gc_model, wind_model, HE_base, tc_haz_base_dict=tc_haz_base_dict):
        haz_base = tc_haz_base_dict[
            str(wind_model_key[wind_model])+str(model_key[gc_model])]

        rng = np.random.RandomState(int(HE_base))
        rnd_ids = [rng.choice(np.arange(500*n+1, 500*(n+1)+1), int(n_ev), replace=False) for n in range(20)]
        event_ids = np.concatenate(rnd_ids).tolist()
        haz_base = haz_base.select(event_id=event_ids)
        
        return haz_base
        
    
    haz_base_distr = {"gc_model": sp.stats.randint(low=h1_min, high=h1_max+1),
                      "wind_model": sp.stats.randint(low=h3_min, high=h3_max+1),
                      "HE_base": sp.stats.randint(0, 2**32 - 1)}
    
    n_ev = 400 #80% of all events are selected for samples of every year
    haz_base_iv = InputVar(haz_base_func, haz_base_distr)
    
    # future
    def haz_fut_func(gc_model, ssp_haz, wind_model, HE_fut, tc_haz_fut_dict=tc_haz_fut_dict):
        haz_fut = tc_haz_fut_dict[
            str(wind_model_key[wind_model])+str(model_key[gc_model])+'_'+str(ssp_haz_key[ssp_haz])]
        
        rng = np.random.RandomState(int(HE_fut))
        rnd_ids = [rng.choice(np.arange(500*n+1, 500*(n+1)+1), int(n_ev), replace=False) for n in range(20)]
        event_ids = np.concatenate(rnd_ids).tolist()
        haz_fut = haz_fut.select(event_id=event_ids)
        
        return haz_fut
    
    haz_fut_distr = {"gc_model": sp.stats.randint(low=h1_min, high=h1_max+1),
                     "ssp_haz": sp.stats.randint(low=h2_min, high=h2_max+1),
                     "wind_model": sp.stats.randint(low=h3_min, high=h3_max+1),
                     "HE_fut": sp.stats.randint(0, 2**32 - 1)}
    
    haz_fut_iv = InputVar(haz_fut_func, haz_fut_distr)

    # exposure
    # present   
    def exp_base_func(mn_exp, exp_base_dict=exp_base_dict):
        return exp_base_dict[str(mn_key[mn_exp])]
    
    # derive distribution from GDP growth scaling factors; see exp_scale_future.py
    exp_base_distr = {"mn_exp": sp.stats.randint(low=e3_min, high=e3_max+1)}
    
    exp_base_iv = InputVar(exp_base_func, exp_base_distr)

    # impact function set - same input parameters for present and future
    def impf_func(v_half):
        impf_set_MED = ImpfSetTropCyclone.from_calibrated_regional_ImpfSet(
            calibration_approach='EDR', q=v_half)
        return impf_set_MED

    impf_distr = {"v_half": sp.stats.uniform(.25, .5)}
    impf_iv = InputVar(impf_func, impf_distr)

    ###########################################################################
    ####################### D: calculate uncertainty ##########################
    ###########################################################################
    # just placeholders for now ... needs more work
    
    calc_imp = CalcDeltaImpact(exp_base_iv, impf_iv, haz_base_iv, 
                               exp_base_iv, impf_iv, haz_fut_iv)
    
    output_imp = calc_imp.make_sample(N=N_samples)
    output_imp.get_samples_df()
    output_imp = calc_imp.uncertainty(output_imp, calc_at_event=True)
    
    ###########################################################################
    ####################### D: calculate sensitivity ##########################
    ###########################################################################
    
    output_imp = calc_imp.sensitivity(output_imp)

    output_imp.to_hdf5(unsequa_dir.joinpath(
        f"unsequa_TC_{year_from_per[period]}_{region}_0{res}as_MIT_{N_samples}_cc.hdf5"))
    
if __name__ == "__main__":
    main(*sys.argv[1:])
