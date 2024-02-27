"""
Adapted for code repository on 2024-02-26

description: Supplementary Figure 4 & 5 - plotting of sensitivity indices - absolute values
             Supplementary Table 3 - largest sensitivity indices - absolute values
             line 84: change salib='S1' (first-order) to 'ST' for total-order sensitivity index

@author: simonameiler
"""

import numpy as np
import logging
import seaborn as sns
import matplotlib.pyplot as plt

#Load Climada modules
from climada.engine.unsequa import UncOutput
from climada.util.constants import SYSTEM_DIR


LOGGER = logging.getLogger(__name__)

###########################################################################
############### A: define constants, functions, paths #####################
###########################################################################

# define paths
unsequa_dir = SYSTEM_DIR/"unsequa"
res_dir = SYSTEM_DIR/"results"

res = 300
region = ['AP', 'IO', 'SH', 'WP']
period = [2050, 2090]
#N_samples = 2**10

model = ['MIT', 'CHAZ', 'STORM', 'IBTrACS']
end_dict = {'MIT': '_abs',
            'STORM': '_abs_v3',
            'IBTrACS': '_abs_v3',
            'CHAZ': '_abs'}

samples_dict = {'MIT': 2**10,
                'STORM': 2**10,
                'IBTrACS': 2**10,
                'CHAZ': 2**10}

# make dictionary of unsequa output
output_dict= {}
for reg in region:
    for per in period:  
        for mdl in model:
            if mdl == 'STORM':
                unsequa_str = f"unsequa_TC_2050_{reg}_0{res}as_{mdl}_{samples_dict[mdl]}{end_dict[mdl]}.hdf5"
                output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
                output_dict[str(reg)+'_'+str(per)+'_'+str(mdl)] = output_imp
            else:            
                unsequa_str = f"unsequa_TC_{per}_{reg}_0{res}as_{mdl}_{samples_dict[mdl]}{end_dict[mdl]}.hdf5"
                output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
                output_dict[str(reg)+'_'+str(per)+'_'+str(mdl)] = output_imp

#%% get largest si per dict
import pandas as pd

# Preallocate a list to collect data
data_list = []

for run, si_dict in output_dict.items():
    si_df = si_dict.get_largest_si('S1')
    st_df = si_dict.get_largest_si('ST')
    reg, year, model = run.split('_')
    
    si_EAD = si_df.param[si_df.metric=='aai_agg'].iloc[0]
    si_rp100 = si_df.param[si_df.metric=='rp100'].iloc[0]
    st_EAD = st_df.param[si_df.metric=='aai_agg'].iloc[0]
    st_rp100 = st_df.param[si_df.metric=='rp100'].iloc[0]
    
    data_list.append([reg, year, model, si_EAD, si_rp100, st_EAD, st_rp100])

# Convert list of lists to DataFrame
df = pd.DataFrame(data_list, columns=['region', 'year', 'model', 'si_EAD', 'si_rp100', 'st_EAD', 'st_rp100'])
df.to_excel(res_dir.joinpath('largest_si_abs.xlsx'), index=False)

#%%
salib = 'ST'
model = ['MIT', 'CHAZ', 'STORM', 'IBTrACS']
lst = ['mn_exp', 'ssp_exp', 'gdp_model', 'HE_fut', 'tcgi_var', 'ssp_haz', 'gc_model', 'wind_model', 'v_half']

# Consider both periods 2050 and 2090
periods = [2050, 2090]
sens1_df_dicts = {per: {} for per in periods}

for per in periods:
    sens1_df_dict = {}
    
    for mdl in model:
        mdl_dict = {}
        
        for reg in region:
            df_S1 = output_dict[f"{reg}_{per}_{mdl}"].get_sensitivity(salib_si=salib)
            
            df = df_S1[["param","aai_agg","rp100"]]
            df = df.set_index("param")
            df.rename(index={'rcp': 'ssp_haz'}, inplace=True)
            df.rename(index={'ensemble_fut': 'HE_fut'}, inplace=True)
            df.rename(index={'ensemble_pres': 'HE_base'}, inplace=True)
            
            df_reord = df.reindex(lst)
            
            mdl_dict[str(reg)] = df_reord
            
        sens1_df_dict[mdl] = mdl_dict
    
    sens1_df_dicts[per] = sens1_df_dict
    
#%%
y = np.arange(len(lst))
colPalette_m = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']

from matplotlib.legend_handler import HandlerBase

class SingleColorHandler(HandlerBase):
    def __init__(self, color, **kwargs):
        super().__init__(**kwargs)
        self.color = color

    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        patch = plt.Rectangle([x0, y0], width, height, facecolor=self.color,
                              edgecolor='k', transform=trans)
        return [patch]

small_shift = 0.15
bar_height = 0.25

output_btt = ['rp100', 'aai_agg']

region_btt = region[::-1]

sens_names = [
    'Exposure urban/rural weighting',
    'SSP exposures',
    'GDP model',
    'Event subsampling future',
    'TCGI moisture variable',
    'SSP hazard',
    'GCM',
    'Wind model',
    'Vulnerability function midpoint']

models = ['MIT', 'CHAZ', 'STORM', 'IBTrACS_p']

# plot
fig, axes = plt.subplots(figsize=(14, 7), nrows=2, ncols=4, sharey=True, sharex=True)
fig.subplots_adjust(bottom=0.20, hspace=0.35, wspace=0.15)

# Specify the models per row, skipping the model for axes[1,2]
models_per_row = [
    ['MIT', 'CHAZ', 'STORM', 'IBTrACS'],
    ['MIT', 'CHAZ', None, 'IBTrACS']
]

# Plot data for both 2050 and 2090
for idx, per in enumerate(periods):
    for ax, mdl in zip(axes[idx], models_per_row[idx]):
        if mdl is None:  # Skip the model for axes[1,2]
            ax.axis('off')  # Turn off the axis
            continue

        sens_dict = sens1_df_dicts[per][mdl]
        for r, reg in enumerate(region_btt):
            sens_dict[reg].replace(np.nan, 0, inplace=True)
            sens_dict[reg][sens_dict[reg] < 0] = 0
            val_tot = 0
            for i, val in enumerate(sens_dict[reg][output_btt[0]]):
                ax.barh(r-small_shift, width=val, height=bar_height, left=val_tot, color=colPalette_m[i])
                val_tot += val
            val_tot = 0
            for i, val in enumerate(sens_dict[reg][output_btt[1]]):
                ax.barh(r+small_shift , width=val, height=bar_height, left=val_tot, color=colPalette_m[i])
                val_tot += val
            axes[0,r].set_title(models[r])
            
        if ax == axes[idx][0]:
            ax.set_ylabel(f'{per}')

for ax in axes[1]:
    ax.set_xlabel(salib)

axes[0,0].annotate('EAD', xy=(0.3, 3.05))
axes[0,0].annotate('RP 100', xy=(0.3, 2.75))

# Adding subfigure labels
labels = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', '', 'g)']
for i, ax in enumerate(axes.ravel()):
    ax.text(-0.1, 1.02, labels[i], transform=ax.transAxes, fontsize=12)

plt.setp(axes[0,0], yticks=range(len(region)), yticklabels=region_btt)
plt.setp(axes[1,0], yticks=range(len(region)), yticklabels=region_btt)

# Legend
# Legend
handles = [plt.Rectangle((0,0), 1, 1, color='none') for _ in colPalette_m]
hmap = {handle: SingleColorHandler(color) for handle, color in zip(handles, colPalette_m)}
axes[1, 0].legend(handles=handles, labels=sens_names, handler_map=hmap, ncol=4, loc='center left',
                  bbox_to_anchor=(-0.15, -0.45), fancybox=False, shadow=False, fontsize='large')

save_fig_str = f"SA_all-models_{salib}_abs.png"
plt.savefig(res_dir.joinpath(save_fig_str), dpi=300, facecolor='w',
            edgecolor='w', orientation='portrait',
            format='png', bbox_inches='tight', pad_inches=0.1)