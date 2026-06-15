import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from constants import *

def save_fig(figure, folder_out, filename_out):
    figure.savefig(f"{folder_out}{filename_out}", dpi=DPI, bbox_inches='tight', transparent=True)
    return None

def truncate_cmap(cmap_name, minval=0.0, maxval=1.0, n=256):
    cmap   = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(minval, maxval, n)) 
    return colors.LinearSegmentedColormap.from_list(f"trunc_{cmap_name}", colors)

def plot_histo(x, y, z, cmap='cividis', cbartxt=False, cbrticks=False, xlim=False, xlog=False, xticks=False, xticklabels=False, Delta=False):
    """ 
    Create figure of histogram. 
    """
    # figure
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    im      = ax.pcolormesh(x, y, z, cmap=cmap)
    # colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', location='top')
    if cbartxt:  cbar.set_label(cbartxt, labelpad=10)
    if cbrticks: cbar.set_ticks(cbrticks)
    if (len(x.shape) == 2): im.set_clim(0, 1) # Probabilities go from 0 to 1
    # labels
    ax.set_xlabel(r"$E_\nu$ [GeV]")
    ax.set_ylabel(r"$\cos \theta_z$")
    # ticks
    if xlim: ax.set_xlim(xlim[0], xlim[1])
    if xlog: ax.set_xscale('log')
    if xticks: ax.set_xticks(xticks)
    if xticklabels: ax.set_xticklabels(xticklabels)
    # twin axis
    if Delta: # Delta = epicentral distance
        ax2 = ax.twinx()
        ax2.set_ylabel(r"$\Delta$ [$^\circ$]")
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks([epi2costh(0),  epi2costh(30),  epi2costh(60),
                        epi2costh(90), epi2costh(120), epi2costh(180)])
        ax2.set_yticklabels(["0", "30", "60", "90", "120", "180"])
    return fig

def plot_binstudy(z, n_pois_norm=25):
    # figure
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.5, 7.0), constrained_layout=True)
    (n_ct_bins, n_enu_bins) = np.shape(z)
    X, Y = np.meshgrid(np.arange(n_enu_bins+1), np.arange(n_ct_bins+1), indexing='ij')
    threshold = np.log10(n_pois_norm)
    z_smaller = np.ma.masked_greater(z, threshold)
    z_greater = np.ma.masked_less_equal(z, threshold)
    blue_cmap = truncate_cmap('Blues', 0.4, 1.0)
    red_cmap  = truncate_cmap('Reds',  0.0, 0.4)
    im1 = ax.pcolormesh(X, Y, z_greater.T, cmap=blue_cmap)
    im2 = ax.pcolormesh(X, Y, z_smaller.T, cmap=red_cmap)
    # colorbar
    cbar1 = fig.colorbar(im1, ax=ax, orientation='horizontal', location='top')
    cbar2 = fig.colorbar(im2, ax=ax, orientation='horizontal', location='top')
    cbar1.set_label(r"$\log_{10}$ (Lowest bin count)", labelpad=10)
    cbar2.set_label("", labelpad=20)
    # labels
    ax.set_xlabel(r"Number of energy bins")
    ax.set_ylabel(r"Number of $\cos(\theta)$ bins")
    return fig

