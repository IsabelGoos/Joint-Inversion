import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import sys
import os
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path, 'scripts'))
from plot_tools import *

path_project  = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
folder_data   = os.path.join(path_project, 'data/')
folder_output = os.path.join(path_project, 'BinsStudy/plots/')
files  = ['DUNE-prem64-OutputExpTruth.root',
          'DUNE-prem128-OutputExpTruth.root',
          ]
channels = ['h_tracks_all', 'h_showers_all']
xticks   = [1, 2, 4, 6, 8, 10, 20, 40]
n_ct_bins  = 300
n_enu_bins = 99
min_bin_count = np.zeros((n_ct_bins, n_enu_bins))
max_bin_count = np.zeros((n_ct_bins, n_enu_bins))

for file in files[:1]: # Only test the first file for now
    for channel in channels[:1]: # Only test the first channel for now
        # create histogram
        histo        = SyntheticData(folder_data, file, channel)
        
        x, y, z      = histo.get_histo()

        enu_regrouped, ct_regrouped, histo_regrouped = histo.regroup_bins_DG(n_pois_norm = 2)
    
        print("*********************************************************************************************************************")
        print(" ")
        print(np.shape(enu_regrouped))
        print(np.shape(ct_regrouped))
        print(np.shape(histo_regrouped))

        fig, ax = plt.subplots(figsize=(10, 6))

        pcm = ax.pcolormesh(enu_regrouped,    # x edges: shape (n_enu + 1,)
                            ct_regrouped,     # y edges: shape (n_ct  + 1,)
                            histo_regrouped,  # data:    shape (n_ct, n_enu)
                            norm=mcolors.LogNorm(),
                            cmap='viridis',
                            shading='flat')   # 'flat' = one color per cell, requires edges

        cbar = fig.colorbar(pcm, ax=ax)
        cbar.set_label(channel, fontsize=12)

        ax.set_xscale('log')
        ax.set_xlabel('Neutrino Energy [GeV]', fontsize=12)
        ax.set_ylabel(r'$\cos\theta$', fontsize=12)
        ax.set_title(f"Regrouped — {channel}\nshape: {histo_regrouped.shape}")

        plt.tight_layout()
        plt.show()







