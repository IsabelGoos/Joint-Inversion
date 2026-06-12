import numpy as np
import matplotlib.pyplot as plt 
import uproot

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.math_utils import *
from utils.plot_utils import *

class SyntheticData:
    def __init__(self, folder, filename, channel):
        """ 
        See readme.txt for the possible channels. 
        """
        self.folder   = folder
        self.filename = filename
        self.channel  = channel

    def get_rootdata(self):
        """ 
        Read the data given in ROOT format by EarthProbe. 
        """
        file = uproot.open(self.folder + self.filename)
        data = file[self.channel]
        return data.to_numpy()
    
    def get_histo(self):
        """ 
        Return the histogram, whatever the channel. 
        """
        data = self.get_rootdata()
        # x, y, z(, b) = number, energy, costheta(, Bjorken-Y)
        x, y, z = data[:3] 
        var_z   = np.squeeze(x).T 
        return y, z, var_z

    def rebin_histo(self, ct_rebin=20, enu_rebin=5):
        """
        Rebin a 2D histogram uniformly by summing bin contents.
        """
        x, y, z = self.get_histo()
        ct_bins, enu_bins = z.shape
        # trim to divisible sizes
        ct_bins_new   = (ct_bins  // ct_rebin)  * ct_rebin
        enu_bins_new  = (enu_bins // enu_rebin) * enu_rebin
        histo_trimmed = z[:ct_bins_new, :enu_bins_new]
        # -> to keep the first bins
        x_trimmed = x[:enu_bins_new] 
        y_trimmed = y[:ct_bins_new]
        # -> to keep the last bins
        #x_trimmed = x[-enu_bins_new:] 
        #y_trimmed = y[-ct_bins_new:]
        # reshape and sum
        histo_rebinned = histo_trimmed.reshape(ct_bins_new  // ct_rebin,  ct_rebin,
                                               enu_bins_new // enu_rebin, enu_rebin)
        histo_rebinned = histo_rebinned.sum(axis=(0, 2))
        x_rebinned = x_trimmed[::(enu_bins // enu_rebin)]
        y_rebinned = y_trimmed[::(ct_bins  // ct_rebin)]
        return x_rebinned, y_rebinned, histo_rebinned

    def regroup_bins(self, n_pois_norm=25):
        """ WORK IN PROGRESSSSSSSSSS
        Rebin a 2D histogram in the sense of arXiv:2408.07015
        """
        x, y, z = self.rebin_histo(ct_rebin=80, enu_rebin=100)
        ct_bins, enu_bins = z.shape
        histo_rebinned = []
        x_rebinned = [0]
        current = np.zeros(ct_bins)
        for i in len(enu_bins):
            current += z[:, i]
            if np.amin(current) >= n_pois_norm:
                histo_rebinned.append(current)
                current = np.zeros(ct_bins)

        # handle leftover last bin
        if np.amin(current) > 0:
            histo_rebinned.append(current.copy())
        x_rebinned = x
        y_rebinned = y
        return histo_rebinned


















# cosas rescatadas, ya pronto las tirare XD

# 80-300
# 100-600

######### Figure 1 paper
#import seaborn as sns
#CMAP1 = 'cividis'
#CMAP2 = sns.color_palette('vlag', as_cmap=True)
#folder = "/Users/igoos/Desktop/projects/Joint-Inversion/BinsStudy/data/" # adapt
#filename_mat = "DUNE-prem64-OutputExpTruth.root" # adapt maybe
##channel = 'h_tracks_numu_cc'
#channel = 'h_osc_numu_nutau'
#cbartxt = r"P$_{\nu_\mu \rightarrow \nu_\tau}$"
#histo = SyntheticData(folder, filename_mat, channel)
#figura, x, y, z = histo.plot_histo(CMAP1, cbartxt, cbrticks=None, Delta=True, xlim=True)
##plt.show()
#folder_out = "/Users/igoos/Desktop/projects/Joint-Inversion/BinsStudy/plots/"
#filename_out = "osc1"
#save_fig(figura, folder_out, filename_out)

         #ax.set_xticks([1, 3, 6, 9, 12, 15])        
          #ax.set_xticklabels([1, 3, 6, 9, 12, 15])
          #ax.set_xticks([1, 2, 4, 6, 8, 10, 20, 30])
          #ax.set_xticklabels([1, 2, 4, 6, 8, 10, 20, 30])

#         n_obs = len(data)
#          # h_tracks_*, h_showers_*, h_int_* -> nvars = 4
#          # h_osc_*                          -> nvars = 3
