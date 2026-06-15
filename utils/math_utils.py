import numpy as np

def epi2costh(epi):
    """ 
    Return cos(theta) for a neutrino with epicentral distance epi (degrees). 
    """
    epi_radians  = epi * np.pi / 180.0
    theta_intern = np.pi - (epi_radians + np.pi) / 2.0
    theta        = np.pi - theta_intern
    return np.cos(theta)
