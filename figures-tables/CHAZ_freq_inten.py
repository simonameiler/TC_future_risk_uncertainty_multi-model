"""
Adapted for code repository on 2024-02-26

description: Supplementary Figure 8 - CHAZ intensity change for different GCMs and TCGIs
             Supplementary Figure 9 - CHAZ frequency change for different GCMs and TCGIs
             plots data from data/CHAZ_freq.xlsx, data/CHAZ_int.xlsx; calulated by CHAZ_freq_inten_check.py
            
@author: simonameiler
"""

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

#Load Climada modules
from climada.util.constants import SYSTEM_DIR

res_dir = SYSTEM_DIR/"results"

TCR_CHAZ = {'CESM2': 2.0,
            'CNRM-CM6-1': 2.22,
            'EC-Earth3': 2.3,
            'IPSL-CM6A-LR': 2.35,
            'MIROC6': 1.55,
            'UKESM1-0-LL': 2.77}

#%%
df = pd.read_excel(res_dir.joinpath("CHAZ_freq.xlsx"))
# Sort models by TCR values
sorted_models = sorted(TCR_CHAZ.keys(), key=lambda x: TCR_CHAZ[x])
prefixes = ['AP', 'IO', 'SH', 'WP']

plt.rcParams.update({
    'font.size': 8,            # Default font size
    'axes.titlesize': 9,       # Font size for figure part labels (A, B, C, etc.)
    'axes.labelsize': 8,       # Font size for axis labels
    'xtick.labelsize': 6,      # Font size for x-tick labels
    'ytick.labelsize': 6,      # Font size for y-tick labels
    'legend.fontsize': 7.5,    # Font size for legend
    'lines.linewidth': 0.75,   # Line weight
})

# Create subplots: 4 rows, 2 columns
fig, axes = plt.subplots(len(prefixes), 2, figsize=(7.25, len(prefixes) * 2.5), sharex='col', sharey='row')

# Loop through each prefix and plot the corresponding columns
for i, prefix in enumerate(prefixes):
    delta1_col = f"{prefix}_delta1"
    delta2_col = f"{prefix}_delta2"

    # Plotting for delta1 column of the current prefix
    sns.boxplot(data=df, x='model', y=delta1_col, hue='TCGI', ax=axes[i, 0], order=sorted_models, width=0.4)
    axes[0, 0].set_title('2050')
    axes[i, 0].axhline(y=0, color='black', linestyle='dotted')  # Adding dotted horizontal line
    axes[i, 0].set_ylabel('Frequency Change (-)')
    axes[i, 0].set_xlabel('')
    #axes[i, 0].text(0.02, 0.95, prefix, transform=axes[i, 0].transAxes, verticalalignment='top')  # Adding region label

    # Secondary y-axis for delta1 column
    ax2_delta1 = axes[i, 0].twinx()
    ax2_delta1.scatter(sorted_models, [TCR_CHAZ[model] for model in sorted_models], marker='*', color='black')
    ax2_delta2 = axes[i, 1].twinx()
    ax2_delta2.scatter(sorted_models, [TCR_CHAZ[model] for model in sorted_models], marker='*', color='black')
    ax2_delta2.set_ylabel('TCR')

    # Plotting for delta2 column of the current prefix
    sns.boxplot(data=df, x='model', y=delta2_col, hue='TCGI', ax=axes[i, 1], order=sorted_models, width=0.4)
    axes[0, 1].set_title('2090')
    axes[i, 1].axhline(y=0, color='black', linestyle='dotted')  # Adding dotted horizontal line
    axes[i, 1].set_ylabel('')
    axes[i, 1].set_xlabel('')
    #axes[i, 1].text(0.02, 0.95, prefix, transform=axes[i, 1].transAxes, verticalalignment='top')  # Adding region label
    axes[i, 0].tick_params(axis='x', rotation=90)
    axes[i, 1].tick_params(axis='x', rotation=90)
    # Remove the legends for each subplot
    axes[i, 0].get_legend().remove()
    axes[i, 1].get_legend().remove()
    # Adding subfigure labels
    labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    axes[i,0].text(-0.3, 0.5, prefix, transform=axes[i,0].transAxes,
                        fontsize=9)
    for j, ax in enumerate(axes.ravel()):
        ax.text(-0.05, 1.05, labels[j], transform=ax.transAxes, fontsize=9, fontweight='bold')

    sns.despine()

handles, labels = axes[1,1].get_legend_handles_labels()
handles2 = Line2D([0], [0], marker='*', color='k', label='TCR', linestyle = 'None',
                          markerfacecolor='k')
handles.append(handles2)
#axes[1,1].legend(handles=handles, labels=['CRH', 'SD', 'TCR'], bbox_to_anchor=(1.9, -0.1), handletextpad=0.2)
fig.legend(handles=handles, labels=['CRH', 'SD', 'TCR'], loc='center right', bbox_to_anchor=(1.05, 0.5), fontsize=6.5)

# Adjusting the layout
#plt.tight_layout()

save_fig_str = "CHAZ_frequency.png"
plt.savefig(res_dir.joinpath(save_fig_str), dpi=300, facecolor='w',
            edgecolor='w', orientation='portrait', 
            format='png', bbox_inches='tight', pad_inches=0.1)

#%% repeat the story for intensity changes
df = pd.read_excel(res_dir.joinpath("CHAZ_int.xlsx"))
# Sort models by TCR values
sorted_models = sorted(TCR_CHAZ.keys(), key=lambda x: TCR_CHAZ[x])
prefixes = ['AP', 'IO', 'SH', 'WP']
# Create subplots: 4 rows, 2 columns
fig, axes = plt.subplots(len(prefixes), 2, figsize=(7.25, len(prefixes) * 2.5), sharex='col', sharey='row')

# Loop through each prefix and plot the corresponding columns
for i, prefix in enumerate(prefixes):
    delta1_col = f"{prefix}_delta1"
    delta2_col = f"{prefix}_delta2"

    # Plotting for delta1 column of the current prefix
    sns.boxplot(data=df, x='model', y=delta1_col, hue='TCGI', ax=axes[i, 0], order=sorted_models, width=0.4)
    axes[0, 0].set_title('2050')
    axes[i, 0].axhline(y=0, color='black', linestyle='dotted')  # Adding dotted horizontal line
    axes[i, 0].set_ylabel('Intensity Change (m/s)')
    axes[i, 0].set_xlabel('')
    #axes[i, 0].text(0.02, 0.95, prefix, transform=axes[i, 0].transAxes, verticalalignment='top')  # Adding region label

    # Secondary y-axis for delta1 column
    ax2_delta1 = axes[i, 0].twinx()
    ax2_delta1.scatter(sorted_models, [TCR_CHAZ[model] for model in sorted_models], marker='*', color='black')
    ax2_delta2 = axes[i, 1].twinx()
    ax2_delta2.scatter(sorted_models, [TCR_CHAZ[model] for model in sorted_models], marker='*', color='black')
    ax2_delta2.set_ylabel('TCR')

    # Plotting for delta2 column of the current prefix
    sns.boxplot(data=df, x='model', y=delta2_col, hue='TCGI', ax=axes[i, 1], order=sorted_models, width=0.4)
    axes[0, 1].set_title('2090')
    axes[i, 1].axhline(y=0, color='black', linestyle='dotted')  # Adding dotted horizontal line
    axes[i, 1].set_ylabel('')
    axes[i, 1].set_xlabel('')
    #axes[i, 1].text(0.02, 0.95, prefix, transform=axes[i, 1].transAxes, verticalalignment='top')  # Adding region label
    axes[i, 0].tick_params(axis='x', rotation=90)
    axes[i, 1].tick_params(axis='x', rotation=90)
    # Remove the legends for each subplot
    axes[i, 0].get_legend().remove()
    axes[i, 1].get_legend().remove()
    # Adding subfigure labels
    labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    axes[i,0].text(-0.3, 0.5, prefix, transform=axes[i,0].transAxes,
                        fontsize=9)
    for j, ax in enumerate(axes.ravel()):
        ax.text(-0.05, 1.05, labels[j], transform=ax.transAxes, fontsize=9, fontweight='bold')

    sns.despine()

handles, labels = axes[1,1].get_legend_handles_labels()
handles2 = Line2D([0], [0], marker='*', color='k', label='TCR', linestyle = 'None',
                          markerfacecolor='k')
handles.append(handles2)
#axes[1,1].legend(handles=handles, labels=['CRH', 'SD', 'TCR'], bbox_to_anchor=(1.9, -0.1), handletextpad=0.2)
fig.legend(handles=handles, labels=['CRH', 'SD', 'TCR'], loc='center right', bbox_to_anchor=(1.05, 0.5), fontsize=6.5)

# Adjusting the layout
#plt.tight_layout()

save_fig_str = "CHAZ_intensity.png"
plt.savefig(res_dir.joinpath(save_fig_str), dpi=300, facecolor='w',
            edgecolor='w', orientation='portrait', 
            format='png', bbox_inches='tight', pad_inches=0.1)