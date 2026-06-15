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
            x_trimmed = x[:enu_bins_new] 
            y_trimmed = y[:ct_bins_new]
            histo_trimmed = z[:ct_bins_new, :enu_bins_new] 
        else:
            # -> to keep the last ct_rebin and enu_rebin bins
            x_trimmed = x[-enu_bins_new:] 
            y_trimmed = y[-ct_bins_new:]
            histo_trimmed = z[-ct_bins_new:, -enu_bins_new:]
        # reshape
        x_rebinned = x_trimmed[::(enu_bins // enu_rebin)]
        y_rebinned = y_trimmed[::(ct_bins  // ct_rebin)]
        histo_rebinned = histo_trimmed.reshape(ct_bins_new  // ct_rebin,  ct_rebin,
                                               enu_bins_new // enu_rebin, enu_rebin)
        # sum
        histo_rebinned = histo_rebinned.sum(axis=(0, 2))
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


















