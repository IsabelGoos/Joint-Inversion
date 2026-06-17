import numpy as np
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

    def rebin_histo(self, ct_rebin=20, enu_rebin=5, firsts=True):
        """
        Rebin a 2D histogram uniformly by summing bin contents.
        """
        x, y, z = self.get_histo()
        ct_bins, enu_bins = z.shape
        # trim to divisible sizes
        ct_bins_new   = (ct_bins  // ct_rebin)  * ct_rebin
        enu_bins_new  = (enu_bins // enu_rebin) * enu_rebin
        if firsts:
            # -> to keep the first ct_rebin and enu_rebin bins
            x_trimmed = x[:(enu_bins_new+1)] 
            y_trimmed = y[:(ct_bins_new+1)]
            histo_trimmed = z[:ct_bins_new, :enu_bins_new] 
        else:
            # -> to keep the last ct_rebin and enu_rebin bins
            x_trimmed = x[-(enu_bins_new+1):] 
            y_trimmed = y[-(ct_bins_new+1):]
            histo_trimmed = z[-ct_bins_new:, -enu_bins_new:]
        # reshape
        x_rebinned = x_trimmed[::(enu_bins // enu_rebin)]
        y_rebinned = y_trimmed[::(ct_bins  // ct_rebin)]
        histo_rebinned = histo_trimmed.reshape(ct_bins_new  // ct_rebin,  ct_rebin,
                                               enu_bins_new // enu_rebin, enu_rebin)
        # sum
        histo_rebinned = histo_rebinned.sum(axis=(0, 2))
        return x_rebinned, y_rebinned, histo_rebinned

    def regroup_bins_DG(self, n_pois_norm=25, ct_rebin=80):
        """ 
        Rebin a 2D histogram in the sense of arXiv:2408.07015
        """
        x, y, z = self.rebin_histo(ct_rebin=ct_rebin, enu_rebin=100)
        ct_bins, enu_bins = z.shape
        # Stores the columns of the regrouped histogram (each column corresponds to a merged energy bin).
        histo_regrouped = []
        # Contains the current event counts of all cos(theta) bins, up to the i-th energy bin.
        counts_allct_uptoithe = np.zeros(ct_bins)
        # Upper edges of the regrouped energy bins.
        enu_upper_edges = [x[0]]
        for i in range(enu_bins):
            counts_allct_uptoithe += z[:, i]
            # Regrouping condition
            if (np.amin(counts_allct_uptoithe) >= n_pois_norm):
                histo_regrouped.append(counts_allct_uptoithe.copy())
                enu_upper_edges.append(x[i+1]) 
                counts_allct_uptoithe = np.zeros(ct_bins)
        # Handle leftover last bin 
        if (np.amin(counts_allct_uptoithe) > 0) & (np.amin(counts_allct_uptoithe) < n_pois_norm):
            histo_regrouped.append(counts_allct_uptoithe.copy())
            enu_upper_edges.append(x[-1])
        # Create the regrouped histogram 
        x_regrouped = np.array(enu_upper_edges)
        y_regrouped = y
        histo_regrouped = np.column_stack(histo_regrouped)
        return x_regrouped, y_regrouped, histo_regrouped

    def get_chi2(self, z_expected, n_pois_norm=25, ct_rebin=80):
        """ 
        Compute chi2 values with respect to an expected histogram.
        """
        x_observed, y_observed, z_observed = self.regroup_bins_DG(n_pois_norm=n_pois_norm, ct_rebin=ct_rebin)
        chi2_vals = chi2(z_expected, z_observed)
        return chi2_vals




















