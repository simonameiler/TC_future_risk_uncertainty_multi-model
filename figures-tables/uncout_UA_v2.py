"""
Adapted for code repository on 2024-02-26

description: Figure 2 - plotting kernel density distribution EAD
             Supplementary Figure 2 - plotting kernel density distribution rp100
             Supplementary Table 1 - saving values describing the locations of max density
            
@author: simonameiler
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd
from climada.engine.unsequa import UncOutput
from climada.util.constants import SYSTEM_DIR
import matplotlib.patches as patches
from matplotlib.lines import Line2D

unsequa_dir = SYSTEM_DIR/"unsequa"
res_dir = SYSTEM_DIR/"results"

res = 300
region = ['AP', 'IO', 'SH', 'WP']
fut_year = 2050

model = ['MIT', 'CHAZ', 'STORM', 'IBTrACS']
end_dict = {'MIT': '_main',
            'STORM': '_main_v3',
            'IBTrACS': '_main_v3',
            'CHAZ': '_main'}

samples_dict = {'MIT': 2**10,
                'STORM': 2**10,
                'IBTrACS': 2**10,
                'CHAZ': 2**10}

def construct_output_df(year, models, regions, samples_dict, end_dict, res, unsequa_dir):
    """Construct output dataframe for a given year."""
    output_dict= {}
    for reg in regions:
        for mdl in models:
            # Skip STORM for 2090
            if year == 2090 and mdl == "STORM":
                continue

            unsequa_str = f"unsequa_TC_{year}_{reg}_0{res}as_{mdl}_{samples_dict[mdl]}{end_dict[mdl]}.hdf5"
            output_imp = UncOutput.from_hdf5(unsequa_dir.joinpath(unsequa_str))
            output_dict[str(reg)+'_'+str(mdl)] = output_imp

    output_df = pd.DataFrame()
    for reg in regions:
        for mdl in models:
            # Check if the model-year combination is in the output_dict before accessing
            if str(reg)+'_'+str(mdl) in output_dict:
                output_df[str(reg)+'_'+str(mdl)+'_EAD_unc'] = output_dict[str(reg)+'_'+str(mdl)].aai_agg_unc_df
                output_df[str(reg)+'_'+str(mdl)+'_rp100_unc'] = output_dict[str(reg)+'_'+str(mdl)].freq_curve_unc_df.rp100
            
    return output_df


# Fetch and construct the output_df here
years = [2050, 2090]
output_dfs = {}
for year in years:
    output_dfs[year] = construct_output_df(year, model, region, samples_dict, end_dict, res, unsequa_dir)

color_dict = {
    'MIT': '#2b83ba',
    'CHAZ': '#abdda4',
    'STORM': '#fdae61',
    'IBTrACS': '#d7191c'
}

labels_dict = {
    (0,0): 'A', (0,1): 'B', (0,2): 'C', (0,3): 'D',
    (1,0): 'E', (1,1): 'F', (1,2): 'G', (1,3): 'H'
}

# Set font sizes and line weights
plt.rcParams.update({
    'font.size': 8,            # Default font size
    'axes.titlesize': 9,       # Font size for figure part labels (A, B, C, etc.)
    'axes.labelsize': 8,       # Font size for axis labels
    'xtick.labelsize': 6,      # Font size for x-tick labels
    'ytick.labelsize': 6,      # Font size for y-tick labels
    'legend.fontsize': 7.5,    # Font size for legend
    'lines.linewidth': 0.28,   # Line weight
})

# Plotting
for mtrc in ['EAD', 'rp100']:
    # Create a new figure for each metric
    fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(7.25, 3.6), sharex=False, sharey=True)
    
    # Max density table for each metric
    max_density_table = pd.DataFrame(columns=['Year', 'Column', 'Region', 'Max Density Point'])
    
    # Inner loop for each year
    for idx, year in enumerate([2050, 2090]):
        
        output_df = construct_output_df(year, model, region, samples_dict, end_dict, res, unsequa_dir)
        
        for r, reg in enumerate(region):
            cols = [f'{reg}_{mdl}_{mtrc}_unc' for mdl in model if f'{reg}_{mdl}_{mtrc}_unc' in output_df.columns]
            
            # Color palette based on available models
            custom_colors = [color_dict[mdl] for mdl in model if f'{reg}_{mdl}_{mtrc}_unc' in output_df.columns]
            
            df_plot = output_df.loc[:, cols]
            sns.kdeplot(df_plot, fill=False, legend=False, ax=ax[idx,r], palette=custom_colors, linewidth=1.0)
            
            sns.despine()
            ax[0,r].set_xlim(-3,10)
            ax[1,r].set_xlim(-3,25)
            ax[0,r].text(0.5, 1.05, reg, transform=ax[0,r].transAxes)
            ax[idx,r].text(-0.15, 1.05, labels_dict[idx,r], transform=ax[idx,r].transAxes, 
                          fontsize=9, fontweight='bold')
            #ax[:,r].get_legend().remove()
            
            # Compute kernel density estimate
            for col, color in zip(cols, custom_colors):
                kde = gaussian_kde(df_plot[col].dropna())
                x_vals = np.linspace(-3, 10, 1000)
                y_vals = kde(x_vals)
                
                # Find the maximum density point
                max_density_point = x_vals[np.argmax(y_vals)]
                
                # Append to the max_density_table for each metric
                max_density_table = max_density_table.append({
                    'Year': year, 'Column': col, 'Region': reg, 
                    'Max Density Point': max_density_point}, ignore_index=True)
        
            ax[idx,1].set_xlabel(f'{mtrc} change {year} (%)')
        #ax[1,1].set_xlabel(f'{mtrc} change {year} (%)')
            ax[idx,1].xaxis.set_label_coords(1.25, -.2)
        #ax[1,1].xaxis.set_label_coords(1.25, -.15)
        
    # custom_handles = [
    #      patches.Rectangle((0, 0), 1, 1, color='#2b83ba', alpha=1.),
    #      patches.Rectangle((0, 0), 1, 1, color='#abdda4', alpha=1.),
    #      patches.Rectangle((0, 0), 1, 1, color='#fdae61', alpha=1.),
    #      patches.Rectangle((0, 0), 1, 1, color='#d7191c', alpha=1.)
    #  ]
    custom_handles = [
    Line2D([0], [0], color='#2b83ba', lw=1.0),
    Line2D([0], [0], color='#abdda4', lw=1.0),
    Line2D([0], [0], color='#fdae61', lw=1.0),
    Line2D([0], [0], color='#d7191c', lw=1.0)
        ]
    
    custom_labels = ['MIT', 'CHAZ', 'STORM', 'IBTrACS_p']

    ax[1,3].legend(handles=custom_handles, labels=custom_labels, loc="upper left", bbox_to_anchor=(1.1, 1.5), fontsize=6.5)
    plt.tight_layout()
    
    # # # Save max_density_table for the metric
    # max_density_table.to_csv(res_dir.joinpath(f'max_density_points_{mtrc}.csv'), index=False)
    
    # # Save the figure for the metric
    save_fig_str = f"UA_all-models_{mtrc}_v2.png"
    plt.savefig(res_dir.joinpath(save_fig_str), dpi=300, facecolor='w', edgecolor='w', orientation='portrait', format='png', bbox_inches='tight', pad_inches=0.1)
