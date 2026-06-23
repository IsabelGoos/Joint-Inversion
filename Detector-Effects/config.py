import numpy as np
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

import os



# Some paths (needs to be adapted by the user)
FOLDER = os.getcwd()
FOLDER = os.path.dirname(FOLDER)
PLOTS_PATH       = FOLDER + "/plots/"
EARTHMODELS_PATH = FOLDER + "/data/earthmodels/"
RESULTS_PATH     = FOLDER + "/data/results/"



# Physics constants

# Properties of the Earth
R_EARTH   = 6371.0   # km
R_LMANTLE = 5711.0   # km
R_CMB     = 3480.0   # km
R_ICORE   = 1230.0   # km
M_EARTH = 5.9722e24  # kg          (+- 0.00060e24 kg,      0.00010 %)
I_EARTH = 8.01736e31 # kg*km^2     (+- 0.00097e31 kg*km^2, 0.00012 %)

# Properties related to particle physics
M_P = 1.67262e-24 # proton mass in g
N_A = 6.02214e23  # Avogadro number in mol^-1
F_CONV = 1 / (M_P * N_A)



# Settings for plots

# Basic settings
DPI = 300

# Colors

# Straightforward definitions
COLORS_BLUE = sns.color_palette('Blues', as_cmap=True)
COLORS_BLUE_11  =  cm.Blues(np.linspace(0.3, 0.9, 11))
COLORS_BLUE_20  =  cm.Blues(np.linspace(0.3, 0.9, 20))
COLORS_BLUE_50  =  cm.Blues(np.linspace(0.3, 0.9, 50))
COLORS_GREEN_11 = cm.Greens(np.linspace(0.3, 0.9, 11))
CMAP1 = 'cividis'
CMAP2 = sns.color_palette('vlag', as_cmap=True)
CB_9  = ['#377eb8', '#ff7f00', '#4daf4a',
         '#f781bf', '#a65628', '#984ea3',
         '#999999', '#e41a1c', '#dede00']
PLASMA_50  = cm.plasma(np.linspace(0, 1, 50))
PLASMA_100 = cm.plasma(np.linspace(0, 1, 100))

# Custommade definitions
START_COLOR   = '#8EB956'
END_COLOR     = '#A03135'
CMAP_EARTH    = LinearSegmentedColormap.from_list('discrete_red_green', [START_COLOR, END_COLOR])
CMAP_EARTH_6  = [CMAP_EARTH(i / 5) for i in range(6)]
