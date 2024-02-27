"""
Adapted for code repository on 2024-02-26

description: Figure 1 - plotting of TC risk change drivers EAD
             Supplementary Figure 1 - plotting of TC risk change drivers rp100
             line 102: change metric="EAD" to "rp100"
             Saving values describing the boxplots of Figure 1, SI Fig. 1 in more detail.
            
@author: simonameiler
"""

import pandas as pd
import logging
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
deltas = ["main", "cc", "soc"]
delta_dict = {"main": "total",
              "cc": "CC",
              "soc": "SOC"}

samples_dict = {'MIT': 2**10,
                'STORM': 2**10,
                'IBTrACS': 2**10,
                'CHAZ': 2**10}

model = ['MIT', 'CHAZ', 'STORM', 'IBTrACS']

# make dictionary of unsequa output
output_dict= {}
for reg in region:
    for delta in deltas:
        for per in period:
            for mdl in model:
                if mdl == "IBTrACS":
                    unsequa_str = f"unsequa_TC_{per}_{reg}_0{res}as_{mdl}_{samples_dict[mdl]}_{delta}_v3.hdf5"
                    output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
                    output_dict[str(reg)+'_'+str(per)+'_'+str(delta_dict[delta])+'_'+str(mdl)] = output_imp
                elif mdl == 'STORM':
                    unsequa_str = f"unsequa_TC_2050_{reg}_0{res}as_{mdl}_{samples_dict[mdl]}_{delta}_v3.hdf5"
                    output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
                    output_dict[str(reg)+'_'+str(per)+'_'+str(delta_dict[delta])+'_'+str(mdl)] = output_imp
                else:
                    unsequa_str = f"unsequa_TC_{per}_{reg}_0{res}as_{mdl}_{samples_dict[mdl]}_{delta}.hdf5"
                    output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
                    output_dict[str(reg)+'_'+str(per)+'_'+str(delta_dict[delta])+'_'+str(mdl)] = output_imp
                    
#%% prep dataframe 

delta_df_dict = {}

for reg in region:
    delta_df_mdl = {}  # Initialize for each region
    for mdl in model:
        delta_df_list = []
        for delta in deltas:
            for per in period:
                df = pd.DataFrame()
                df['EAD'] = output_dict[str(reg)+'_'+str(per)+'_'+str(delta_dict[delta])+'_'+str(mdl)].aai_agg_unc_df
                df['rp100'] = output_dict[str(reg)+'_'+str(per)+'_'+str(delta_dict[delta])+'_'+str(mdl)].freq_curve_unc_df.rp100
                df['delta'] = delta_dict[delta]
                df['year'] = per
                delta_df_list.append(df)

        delta_df = pd.concat(delta_df_list)
        
        # These operations should be outside the deltas and period loops but inside the model loop
        df_cc = delta_df[delta_df.delta == 'CC']
        df_soc = delta_df[delta_df.delta == 'SOC']
        df_sum = df_cc + df_soc
        df_sum.delta = 'sum'
        df_sum.year[df_sum.year == 4100] = 2050
        df_sum.year[df_sum.year == 4180] = 2090
        delta_df = pd.concat([delta_df, df_sum])
        
        delta_df_mdl[mdl] = delta_df

    delta_df_dict[reg] = delta_df_mdl


#%% plot 
metric = "EAD"

# Function to concatenate dataframes with a new source column
def concat_dfs(dfs, labels):
    return pd.concat([df.assign(source=label) for df, label in zip(dfs, labels)], ignore_index=True)

# Updated to 2 columns
fig, ax = plt.subplots(ncols=2, nrows=4, figsize=(8,8), sharex=True, sharey=True)

labels = iter('abcdefgh')  # create an iterator for subplot labels
#custom_colors = ['lightblue', 'orange', 'green', 'purple']
custom_colors = ['#2b83ba', '#abdda4', '#fdae61', '#d7191c']

for r, reg in enumerate(['AP', 'IO', 'SH', 'WP']):
    dfs = [delta_df_dict[reg]['MIT'], delta_df_dict[reg]['CHAZ'], delta_df_dict[reg]['STORM'], delta_df_dict[reg]['IBTrACS']]
    data_labels = ['MIT', 'CHAZ', 'STORM', 'IBTrACS']
    combined_df = concat_dfs(dfs, data_labels)
    
    # Set all values of the source "STORM" in the year 2090 to NaN
    combined_df.loc[(combined_df['source'] == 'STORM') & (combined_df['year'] == 2090), metric] = float('nan')
    
    for c, year in enumerate([2050, 2090]):
        # Filter for each year
        year_df = combined_df[combined_df['year'] == year]

        sns.boxplot(data=year_df, x="delta", hue="source", y=metric, width=0.4, palette=custom_colors,
                    order=["CC", "SOC", "sum", "total"], showfliers=False, ax=ax[r, c], dodge=True)
        
        ax[r, c].get_legend().remove()
        ax[r, c].set_yscale('symlog')
        sns.despine()
        ax[r, c].set(xlabel="", ylabel="")
        ax[0, c].set_title(f"{year}")
        if c == 0:
            ax[r, c].set(ylabel=f"\u0394 {metric} (%)")
            ax[r, c].text(-0.25, 0.5, reg, transform=ax[r, c].transAxes, ha='center', va='center', fontsize=12)  # Add 'reg' info
        # else:
        #     ax[r, c].text(1.15, 0.5, reg, transform=ax[r, c].transAxes, ha='center', va='center', fontsize=12, fontweight='bold')  # Add 'reg' info
        ax[r, c].axhline(0, ls='dotted', color='k')

        # Add subplot labels
        label = next(labels)  # get the next label from the iterator
        ax[r, c].text(-0.1, 1.05, label+')', transform=ax[r, c].transAxes, fontsize=12)

# Create custom legend outside of the plots
fig_labels = ['MIT', 'CHAZ', 'STORM', 'IBTrACS_p']
handles = [mpatches.Patch(color=custom_colors[i], label=fig_labels[i]) for i in range(len(fig_labels))]
fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.15, 0.5))

plt.tight_layout()

# save_fig_str = f"delta_all_{metric}.png"
# plt.savefig(res_dir.joinpath(save_fig_str), dpi=300, facecolor='w', 
#             edgecolor='w', orientation='portrait', 
#             format='png', bbox_inches='tight', pad_inches=0.1) 

#%% store boxplot statistics
all_stats = []

for r, reg in enumerate(['AP', 'IO', 'SH', 'WP']):
    dfs = [delta_df_dict[reg]['MIT'], delta_df_dict[reg]['CHAZ'], delta_df_dict[reg]['STORM'], delta_df_dict[reg]['IBTrACS']]
    data_labels = ['MIT', 'CHAZ', 'STORM', 'IBTrACS']
    combined_df = concat_dfs(dfs, data_labels)

    for c, year in enumerate([2050, 2090]):
        year_df = combined_df[combined_df['year'] == year]
        for delta_value in ["CC", "SOC", "sum", "total"]:
            for source in data_labels:
                df_subset = year_df[(year_df['delta'] == delta_value) & (year_df['source'] == source)]
                stats = df_subset[metric].describe(percentiles=[.25, .5, .75])
                stats['reg'] = reg
                stats['year'] = year
                stats['delta'] = delta_value
                stats['source'] = source
                all_stats.append(stats)

# Convert to DataFrame
stats_df = pd.DataFrame(all_stats)

# 2. Rename the 'source' column to 'model'
stats_df = stats_df.rename(columns={'source': 'model'})

stats_df = stats_df.drop(columns=['count'])


# 1. Format floats
float_cols = ['mean', 'std', 'min', '25%', '50%', '75%', 'max']
for col in float_cols:
    stats_df[col] = stats_df[col].apply(lambda x: round(x, 2))
    
cols_order = ['reg', 'delta', 'year', 'model', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
stats_df = stats_df[cols_order]

# set all STORM values of 2090 to 'N/A'
stats_df.loc[(stats_df['model'] == 'STORM') & (stats_df['year'] == 2090), ['25%', '50%', '75%', 'std', 'min', 'max', 'mean']] = 'N/A'

# Save to Excel
# stats_df.to_excel(res_dir.joinpath(f'boxplot_statistics_all-models_{metric}.xlsx'), index=False)
