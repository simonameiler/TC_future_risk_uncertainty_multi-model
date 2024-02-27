"""
Adapted for code repository on 2024-02-26

description: generation of annual growth rates relative to 2020, from SSP public database

@author: evelynm; adapted @simonameiler
"""
import pandas as pd
import numpy as np

gdp_file = '/Users/simonameiler/Documents/WCR/Active_Research_Projects/TC_future/data/iamc_db.csv'

"""
This input csv file was downloaded from the SSP public database (v2.0)
Keywan Riahi, Detlef P. van Vuuren, Elmar Kriegler, Jae Edmonds, Brian C. O’Neill, Shinichiro Fujimori, Nico Bauer, Katherine Calvin, Rob Dellink, Oliver Fricko, Wolfgang Lutz, Alexander Popp, Jesus Crespo Cuaresma, Samir KC, Marian Leimbach, Leiwen Jiang, Tom Kram, Shilpa Rao, Johannes Emmerling, Kristie Ebi, Tomoko Hasegawa, Petr Havlík, Florian Humpenöder, Lara Aleluia Da Silva, Steve Smith, Elke Stehfest, Valentina Bosetti, Jiyong Eom, David Gernaat, Toshihiko Masui, Joeri Rogelj, Jessica Strefler, Laurent Drouet, Volker Krey, Gunnar Luderer, Mathijs Harmsen, Kiyoshi Takahashi, Lavinia Baumstark, Jonathan C. Doelman, Mikiko Kainuma, Zbigniew Klimont, Giacomo Marangoni, Hermann Lotze-Campen, Michael Obersteiner, Andrzej Tabeau, Massimo Tavoni.
The Shared Socioeconomic Pathways and their energy, land use, and greenhouse gas emissions implications: An overview, Global Environmental Change, Volume 42, Pages 153-168, 2017,
ISSN 0959-3780, DOI:110.1016/j.gloenvcha.2016.05.009

Selection: 1. Region - all countries, 2. Scenarios - GDP - IASA/OECD/PIK, 3. Variable - GDP (growth Total)
"""

df_growth = pd.read_csv(gdp_file, header=0)

def calc_growthrate(df_growth):
    growth = df_growth[['Model', 'Scenario', 'Region']]
    proj_avail = np.array(df_growth.columns[5:-1].astype(int).values)
    years = np.arange(2021,2100)
    growth['2020'] = 1
    for year in years:
        proj_sel = proj_avail[((proj_avail-year)<=0) & ((proj_avail-year)>-5)]
        growth[str(year)] = growth.apply(
            lambda row: row[str(year-1)]*(1+df_growth.iloc[row.name][str(proj_sel[0])]/100),
            axis=1)
    
    return growth
        
growth = calc_growthrate(df_growth)
growth.to_csv('/Users/simonameiler/Documents/WCR/Active_Research_Projects/TC_future/data/ssps_gdp_annual.csv')
