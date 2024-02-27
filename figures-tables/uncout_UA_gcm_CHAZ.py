"""
Adapted for code repository on 2024-02-26

description: Supplementary Figure 6 - plotting of TC risk change (EAD) for each GCM from CHAZ
             Supplementary Figure 7 - plotting of TC risk change (rp100) for each GCM from CHAZ
             line 70: change metric="EAD" to "rp100"

            
@author: simonameiler
"""

import numpy as np
import copy as cp
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#Load Climada modules
from climada.util.constants import SYSTEM_DIR # loads default directory paths for data
from climada.engine.unsequa import UncOutput

    
LOGGER = logging.getLogger(__name__)

###########################################################################
############### A: define constants, functions, paths #####################
###########################################################################

# define paths
unsequa_dir = SYSTEM_DIR/"unsequa"
res_dir = SYSTEM_DIR/"results"

res = 300
ref_year = 2005
region = ['AP', 'IO', 'SH', 'WP']
period = [2050, 2090]
N_samples = 2**10

# make dictionary of unsequa output
output_dict= {}
for reg in region:
    for per in period:
        unsequa_str = f"unsequa_TC_{per}_{reg}_0{res}as_CHAZ_{N_samples}_main.hdf5"
        output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
        output_dict[str(reg)+'_'+str(per)] = output_imp

# make dataframe of all results over regions and periods
# samples_df is identical for all UncOutput objects in the output_dict
# idea: extend the samples_df with the results of interest
# first, get indexes where hazard and exposure SSPs match
ssp245_idx = (output_imp.samples_df.ssp_haz == 1.0) & (output_imp.samples_df.ssp_exp == 2.0)
ssp370_idx = (output_imp.samples_df.ssp_haz == 2.0) & (output_imp.samples_df.ssp_exp == 3.0)
ssp585_idx = (output_imp.samples_df.ssp_haz == 3.0) & (output_imp.samples_df.ssp_exp == 5.0)
ssp_idx = ssp245_idx + ssp370_idx + ssp585_idx

output_df = cp.deepcopy(output_imp.samples_df)
for reg in region:
    for per in period:
        output_df[str(reg)+'_'+str(per)+'_EAD_unc'] = output_dict[
            str(reg)+'_'+str(per)].aai_agg_unc_df
        output_df[str(reg)+'_'+str(per)+'_rp100_unc'] = output_dict[
            str(reg)+'_'+str(per)].freq_curve_unc_df.rp100
        output_df[str(reg)+'_'+str(per)+'_EAD_ssp_unc'] = output_dict[
            str(reg)+'_'+str(per)].aai_agg_unc_df[ssp_idx]
        output_df[str(reg)+'_'+str(per)+'_rp100_ssp_unc'] = output_dict[
            str(reg)+'_'+str(per)].freq_curve_unc_df.rp100[ssp_idx]

#%%
metric = "EAD"

labels_dict = {(0,0): 'a)',
               (0,1): 'b)',
               (1,0): 'c)',
               (1,1): 'd)',
               (2,0): 'e)',
               (2,1): 'f)',
               (3,0): 'g)',
               (3,1): 'h)'}

TCR_CHAZ = {1.0: 2.0,
            2.0: 2.22,
            3.0: 2.3,
            4.0: 2.35,
            5.0: 1.55,
            6.0: 2.77}

models_TCR = [2.0, 2.22, 2.30, 2.35, 1.55, 2.77]
models_srtd = ['miroc6', 'cesm2', 'cnrm6', 'ecearth', 'ipsl6', 'ukmo6']

output_df.gc_model.replace(TCR_CHAZ, inplace=True)
output_df.sort_values('gc_model', inplace=True)

# options for secondary y-axis plots - sort TCR values
TCR_list = list(TCR_CHAZ.values())
TCR_list.sort()
plt_points = np.arange(0,6)

# okay, now make this pretty
fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(12,12), sharex=True, sharey=False)
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
# # Set your custom color palette
customPalette_l = sns.hls_palette(n_colors=2, s=0.6)
customPalette_d = sns.hls_palette(n_colors=2, l=0.4, s=1.)
for r, reg in enumerate(region):
    for p, per in enumerate(period):
        sns.stripplot(
            data=output_df, x="gc_model", y=f"{reg}_{per}_{metric}_unc", hue="tcgi_var",
            marker=".", dodge=True, alpha=.75, zorder=1, legend=True, palette=customPalette_d,
            ax=ax[r,p])

        ax[r,0].text(-0.25, 0.5, reg, transform=ax[r,0].transAxes,
                        fontsize=12, rotation=0)
        ax[r,p].text(-0.1, 1.05, labels_dict[r,p], transform=ax[r,p].transAxes,
                     fontsize=12)
        ax[0,0].set_title('2050')
        ax[0,1].set_title('2090')
        ax[r,p].get_legend().remove()
        ax[r,p].get_yaxis().set_visible(True)
        ax[r,p].set(xlabel='GCMs', ylabel=f'\u0394 {metric} (%)')
        ax[r,p].set_xticklabels(models_srtd, rotation=90)
        
        secax_y = ax[r,p].twinx()
        secax_y.plot(plt_points, TCR_list, marker='*', color='k', markersize=12, ls='')
        secax_y.set_ylabel('TCR')
        sns.despine()
handles, labels = ax[1,1].get_legend_handles_labels()
handles2 = Line2D([0], [0], marker='*', color='k', label='TCR', linestyle = 'None',
                          markerfacecolor='k', markersize=12)
handles.append(handles2)
ax[1,1].legend(handles=handles, labels=['CRH', 'SD', 'TCR'], loc="upper left", bbox_to_anchor=(1.1, 0.25), handletextpad=0)

# save_fig_str = f"UA_TC_risk_CHAZ_{metric}.png"
# plt.savefig(res_dir.joinpath(save_fig_str), dpi=300, facecolor='w',
#             edgecolor='w', orientation='portrait', 
#             format='png', bbox_inches='tight', pad_inches=0.1)


