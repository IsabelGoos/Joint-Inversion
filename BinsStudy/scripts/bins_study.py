import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plot_tools import *

folder = "/Users/igoos/Desktop/projects/Joint-Inversion/BinsStudy/data/"
files  = ["DUNE-prem128-OutputExpTruth.root",
          "HK-prem128-OutputExpTruth.root",
          "ORCA-prem128-OutputExpTruth.root",
          "PERFECT-prem128-OutputExpTruth.root"]
channels = ['h_tracks_all', 'h_showers_all']
folder_out = "/Users/igoos/Desktop/projects/Joint-Inversion/BinsStudy/plots/"
xticks = [1, 2, 4, 6, 8, 10, 20, 40]

for k in range(len(files)):
    for l in range(len(channels)):
        print()
        histo = SyntheticData(folder, files[k], channels[l])
        x, y, z = histo.get_histo()
        figure = plot_histo(x,y,z, cbartxt=channels[l], xlog=True, xticks=xticks, xticklabels=xticks)
        filename_out = files[k].replace("OutputExpTruth.root", channels[l])
        save_fig(figure, folder_out, filename_out)
        plt.close(figure)
    
        n_ct_bins  = 300
        n_enu_bins = 99
        min_bin_count = np.zeros((n_ct_bins, n_enu_bins))
        max_bin_count = np.zeros((n_ct_bins, n_enu_bins))
        for i in range(n_ct_bins):
            for j in range(n_enu_bins):
                print(i)
                n_ct  = i+1
                n_enu = j+1
                x_rebinned, y_rebinned, histo_rebinned = histo.rebin_histo(n_ct, n_enu)
                min_bin_count[i,j] = np.log10(np.amin(histo_rebinned))
                max_bin_count[i,j] = np.log10(np.amax(histo_rebinned))

        figure = plot_binstudy(min_bin_count)
        filename_out = files[k].replace("OutputExpTruth.root", str(channels[l]) + "-binstudy")
        save_fig(figure, folder_out, filename_out)
        plt.close(figure)











