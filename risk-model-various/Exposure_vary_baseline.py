"""
Adapted for code repository on 2024-02-26

description: Exposure - Create nine different Exposure layer after different
             Lit and Pop exponent combinations for the reference year 2005

@author: simonameiler
"""

import logging

# import CLIMADA modules:
from climada.util.constants import SYSTEM_DIR
from climada.entity import LitPop

LOGGER = logging.getLogger(__name__)

res = 300
ref_year = 2005

# ISO3 codes of all countries retrieved from the SSP Database (OECD model)
# with projected GDP growth

all_cntrs = ['ABW', 'AFG', 'AGO', 'ALB', 'ARE', 'ARG', 'ARM', 'AUS', 'AUT', 'AZE', 
             'BDI', 'BEL', 'BEN', 'BFA', 'BGD', 'BGR', 'BHR', 'BHS', 'BIH', 'BLR', 
             'BLZ', 'BOL', 'BRA', 'BRB', 'BRN', 'BTN', 'BWA', 'CAF', 'CAN', 'CHE', 
             'CHL', 'CHN', 'CIV', 'CMR', 'COD', 'COG', 'COL', 'COM', 'CPV', 'CRI', 
             'CUB', 'CYP', 'CZE', 'DEU', 'DJI', 'DNK', 'DOM', 'DZA', 'ECU', 'EGY', 
             'ERI', 'ESP', 'EST', 'ETH', 'FIN', 'FJI', 'FRA', 'GAB', 'GBR', 'GEO', 
             'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'GRC', 'GTM', 'GUY', 'HKG', 'HND', 
             'HRV', 'HTI', 'HUN', 'IDN', 'IND', 'IRL', 'IRN', 'IRQ', 'ISL', 'ISR', 
             'ITA', 'JAM', 'JOR', 'JPN', 'KAZ', 'KEN', 'KGZ', 'KHM', 'KOR', 'KWT', 
             'LAO', 'LBN', 'LBR', 'LBY', 'LCA', 'LKA', 'LSO', 'LTU', 'LUX', 'LVA', 
             'MAC', 'MAR', 'MDA', 'MDG', 'MDV', 'MEX', 'MKD', 'MLI', 'MLT', 'MMR', 
             'MNE', 'MNG', 'MOZ', 'MRT', 'MUS', 'MWI', 'MYS', 'NAM', 'NCL', 'NER', 
             'NGA', 'NIC', 'NLD', 'NOR', 'NPL', 'NZL', 'OMN', 'PAK', 'PAN', 'PER', 
             'PHL', 'PNG', 'POL', 'PRI', 'PRT', 'PRY', 'PSE', 'PYF', 'QAT', 'ROU', 
             'RUS', 'RWA', 'SAU', 'SDN', 'SEN', 'SGP', 'SLB', 'SLE', 'SLV', 'SOM', 
             'SRB', 'STP', 'SUR', 'SVK', 'SVN', 'SWE', 'SWZ', 'SYR', 'TCD', 'TGO', 
             'THA', 'TJK', 'TKM', 'TLS', 'TON', 'TTO', 'TUN', 'TUR', 'TWN', 'TZA', 
             'UGA', 'UKR', 'URY', 'USA', 'UZB', 'VCT', 'VEN', 'VNM', 'VUT', 'WSM', 
             'YEM', 'ZAF', 'ZMB', 'ZWE']

# make exposure from OECD countries
choice_mn = [[0.5, 0.5], [0.5, 1.0], [0.5, 1.5], 
             [1.0, 0.5], [1.0, 1.0], [1.0, 1.5], 
             [1.5, 0.5], [1.5, 1.0], [1.5, 1.5]]

for [m, n] in choice_mn:
    exp = LitPop.from_countries(all_cntrs, res_arcsec=res, reference_year=ref_year, exponents=(m,n))
    exp.set_geometry_points()
    exp.set_lat_lon()
    exp_str = f"litpop_0{res}as_{ref_year}_global_{m}-{n}.hdf5"
    exp.write_hdf5(SYSTEM_DIR.joinpath(exp_str))
    
    # split exposure up for the four regions
    # boundaries of (sub-)basins (lonmin, lonmax, latmin, latmax)
    basin_bounds = {
        # North Atlantic/Eastern Pacific Basin
        'AP': [-180.0, -30.0, 0.0, 65.0],
    
        # Indian Ocean Basin
        'IO': [30.0, 100.0, 0.0, 40.0],
    
        # Southern Hemisphere Basin
        'SH': [-180.0, 180.0, -60.0, 0.0],
    
        # Western Pacific Basin
        'WP': [100.0, 180.0, 0.0, 65.0],
    }
    
    for reg, bounds in basin_bounds.items():
        x_min, x_max, y_min, y_max = bounds
        exp_reg = LitPop()
        exp_reg.gdf = exp.gdf.cx[x_min:x_max, y_min:y_max]
        exp_reg_str = SYSTEM_DIR.joinpath(f"litpop_0{res}as_{ref_year}_{reg}_{m}-{n}.hdf5")
        exp_reg.write_hdf5(exp_reg_str)