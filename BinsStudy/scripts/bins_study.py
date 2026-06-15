import numpy as np
import matplotlib.pyplot as plt

import sys
import os
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path, 'scripts'))
from plot_tools import *

path_project  = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
folder_data   = os.path.join(path_project, 'data/')
folder_output = os.path.join(path_project, 'BinsStudy/plots/')
files  = ['DUNE-prem64-OutputExpTruth.root',
          'HK-prem64-OutputExpTruth.root',
          'ORCA-prem64-OutputExpTruth.root',
          'PERFECT-prem64-OutputExpTruth.root',
          'DUNE-prem128-OutputExpTruth.root',
          'HK-prem128-OutputExpTruth.root',
          'ORCA-prem128-OutputExpTruth.root',
          'PERFECT-prem128-OutputExpTruth.root']
channels = ['h_tracks_all', 'h_showers_all']
xticks   = [1, 2, 4, 6, 8, 10, 20, 40]
n_ct_bins  = 300
n_enu_bins = 99
min_bin_count = np.zeros((n_ct_bins, n_enu_bins))
max_bin_count = np.zeros((n_ct_bins, n_enu_bins))

for file in files:
    for channel in channels:
        # create histogram
        histo        = SyntheticData(folder_data, file, channel)
        x, y, z      = histo.get_histo()
        # plot detected events
        figure       = plot_histo(x, y, z, cbartxt=channel, xlog=True, xticks=xticks, xticklabels=xticks)
        filename_out = f"{file.replace('OutputExpTruth.root','')}{channel}"
        save_fig(figure, folder_output, filename_out)
        plt.close(figure)
        # run over different amount of bins (energy and cos(theta)) and find smallest and highest bin count
        for i in range(n_ct_bins):
            if i % 10 == 0: print(i)
            for j in range(n_enu_bins):
                n_ct  = i+1
                n_enu = j+1
                #x_rebinned, y_rebinned, histo_rebinned = histo.rebin_histo(n_ct, n_enu)
                x_rebinned, y_rebinned, histo_rebinned = histo.rebin_histo(n_ct, n_enu, firsts=False)
                min_bin_count[i, j] = np.log10(np.maximum(np.amin(histo_rebinned), 1e-12))
                max_bin_count[i, j] = np.log10(np.maximum(np.amax(histo_rebinned), 1e-12))
        # plot the smallest bin count
        figure = plot_binstudy(min_bin_count)
        filename_out = f"{file.replace('OutputExpTruth.root','')}{channel}-binstudy"
        save_fig(figure, folder_output, filename_out)
        plt.close(figure)









