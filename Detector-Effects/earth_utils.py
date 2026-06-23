import numpy as  np
import scipy.integrate as integrate
from utils import config as C
from utils import math_utils as M



# This function discretizes the PREM and GERM distributions into "N_shells" shells of equal depth.
# For each shell, it computes average values of density "rho" and electron yield "zoa".
# Note that when averaging the global constraints might (slightly) not be satisfied anymore.

def original_PREM(N_shells=64):

    radii = np.round(np.arange(C.R_EARTH/N_shells, C.R_EARTH+1, C.R_EARTH/N_shells), 2)
    rhos  = np.zeros(N_shells)
    zoas  = np.zeros(N_shells)

    for i in range(N_shells):
        radius1 = radii[i-1] if i>0 else 0.0
        radius2 = radii[i]
        rho_integral = integrate.quad(lambda x: M.earth_rho(x), radius1, radius2)
        zoa_integral = integrate.quad(lambda x: M.earth_zoa(x), radius1, radius2)
        rhos[i] = np.round(rho_integral[0]/(C.R_EARTH/N_shells), 4)
        zoas[i] = np.round(zoa_integral[0]/(C.R_EARTH/N_shells), 4)
        
    return radii, rhos, zoas
    


# This function modifies layers or shells by the percentage pc*100%, 
# starting from the profiles produced by the function original_PREM(N_shells)
# and keeping the mass and moment of inertia correct.
# "shells" is True if we want to analyze individual shells (all the N_shells shells),
# otherwise "core" is True if we want to study the core,
# otherwise "mtz" is True if we want to study the MTZ.
# Exactly one of these has to be True.
# "test" is True if we want to check that the modiified profiles respect the global constraints.

def modified_PREM(N_shells=64, pc=0.03, shells=True, core=False, mtz=False, icore=False, ocore=False, test=False):
    
    # Only for testing
    if (test == True):
        abs_err_M_before = []
        abs_err_I_before = []
        abs_err_M_after = []
        abs_err_I_after = []

    radii, rhos, zoas = original_PREM(N_shells)
    rhos_mod  = np.zeros(N_shells) # mod for modified
    true_vals = np.array([C.M_EARTH, C.I_EARTH]) # true mass and moment of inertia
    coeffs_deep    = []
    coeffs_shallow = []
    
    if (shells == True):
        n_shells_to_analyze = N_shells
    else:
        n_shells_to_analyze = 1

    for i in range(n_shells_to_analyze):

        # Define modified PREM only changing the i'th shell or layer of interest
        for j in range(N_shells):
            rhos_mod[j] = rhos[j]
        if (shells == True): 
            rhos_mod[i] = rhos[i] * (1.0 + pc)
        elif (core == True):
            for j in range(N_shells):
                if (radii[j]<=C.R_CMB + (0.5 * C.R_EARTH/N_shells)):
                    rhos_mod[j] = rhos[j] * (1.0 + pc)
        elif (mtz == True):
            for j in range(N_shells):
                if ( (C.R_EARTH-410.0) >= radii[j] ) and ( radii[j] >= (C.R_EARTH-660.0) ):
                    rhos_mod[j] = rhos[j] * (1.0 + pc)
        elif (icore == True):
            for j in range(N_shells):
                if (radii[j] <= C.R_ICORE + (0.5 * C.R_EARTH/N_shells)):
                    rhos_mod[j] = rhos[j] * (1.0 + pc)
        elif (ocore == True):
            for j in range(N_shells):
                if ( C.R_ICORE + (0.5 * C.R_EARTH/N_shells) < radii[j] ) and ( radii[j]<=C.R_CMB + (0.5 * C.R_EARTH/N_shells) ):
                    rhos_mod[j] = rhos[j] * (1.0 + pc)

        # mod_vals -> computation of M (in the first row) and I (in the second row)
        # for deep (first column) and shallow (second column) layers
        # deep    -> core
        # shallow -> mantle + crust
        mod_vals = np.zeros((2, 2))
        M_core, I_core = 0.0, 0.0
        for j in range(N_shells):
            radius1 = radii[j-1] if j>0 else 0.0
            radius2 = radii[j]
            if (shells == True) or (mtz == True): # in these cases we rescale the core and (mantle+crust) separately
                if (radius2 <= C.R_CMB + (0.5 * C.R_EARTH/N_shells)): # deep shells
                    mod_vals[0, 0] = mod_vals[0, 0] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                    mod_vals[1, 0] = mod_vals[1, 0] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
                else: # shallow shells
                    mod_vals[0, 1] = mod_vals[0, 1] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                    mod_vals[1, 1] = mod_vals[1, 1] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
            else: # in this case we don't touch the core, we rescale the mantle and the crust separately
                if (radius2 <= C.R_CMB + (0.5 * C.R_EARTH/N_shells)): # shells that will not be modified
                    M_core = M_core + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                    I_core = I_core + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
                elif (radius2 <= C.R_LMANTLE): # deep shells that are rescaled
                    mod_vals[0, 0] = mod_vals[0, 0] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                    mod_vals[1, 0] = mod_vals[1, 0] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
                else: # shallow shells that are modified
                    mod_vals[0, 1] = mod_vals[0, 1] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                    mod_vals[1, 1] = mod_vals[1, 1] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I

        # Multiply by 1.0e12 to go from g/cm^3 to kg/km^3
        if (core == True):
            M_core = 1.0e12 * M_core * 4 * np.pi / 3  # kg
            I_core = 1.0e12 * I_core * 8 * np.pi / 15 # kg*km^2
        mod_vals[0, 0] = 1.0e12 * mod_vals[0, 0] * 4 * np.pi / 3  # kg
        mod_vals[0, 1] = 1.0e12 * mod_vals[0, 1] * 4 * np.pi / 3  # kg
        mod_vals[1, 0] = 1.0e12 * mod_vals[1, 0] * 8 * np.pi / 15 # kg*km^2
        mod_vals[1, 1] = 1.0e12 * mod_vals[1, 1] * 8 * np.pi / 15 # kg*km^2

        # Compute absolute errors, only for testing
        if (test == True):
            if (shells == True) or (mtz == True):
                abs_err_M_before = np.append(abs_err_M_before, ((mod_vals[0, 0]+mod_vals[0, 1])-C.M_EARTH))
                abs_err_I_before = np.append(abs_err_I_before, ((mod_vals[1, 0]+mod_vals[1, 1])-C.I_EARTH))
            else:
                abs_err_M_before = np.append(abs_err_M_before, ((mod_vals[0, 0]+mod_vals[0, 1])-(C.M_EARTH-M_core)))
                abs_err_I_before = np.append(abs_err_I_before, ((mod_vals[1, 0]+mod_vals[1, 1])-(C.I_EARTH-I_core)))

        # Compute correction coefficients
        mod_vals_inv = np.linalg.inv(mod_vals)
        if (core == True): 
            true_vals = true_vals - np.array([M_core, I_core])
        coeffs = np.matmul(mod_vals_inv, true_vals)
        coeffs_deep    = np.append(coeffs_deep,    coeffs[0])
        coeffs_shallow = np.append(coeffs_shallow, coeffs[1])

        # Compute corrected modified PREM
        for j in range(N_shells):
            if (shells == True) or (mtz == True):
                if (radii[j] <= C.R_CMB + (0.5 * C.R_EARTH/N_shells)): # deep shells
                    rhos_mod[j] = rhos_mod[j] * coeffs[0]
                else: # shallow shells
                    rhos_mod[j] = rhos_mod[j] * coeffs[1]
            else:
                if (radii[j] > C.R_CMB + (0.5 * C.R_EARTH/N_shells)) and (radii[j] <= C.R_LMANTLE): # deep shells 
                    rhos_mod[j] = rhos_mod[j] * coeffs[0]
                elif (radii[j] > C.R_LMANTLE):  # shallow shells
                    rhos_mod[j] = rhos_mod[j] * coeffs[1]

        # Only for testing
        if (test == True):
            mod_vals_test = np.zeros((2, 2))
            for j in range(N_shells):
                radius1 = radii[j-1] if j>0 else 0.0
                radius2 = radii[j]
                if (shells == True) or (mtz == True):
                    if (radius2 <= C.R_CMB + (0.5 * C.R_EARTH/N_shells)):
                        mod_vals_test[0, 0] = mod_vals_test[0, 0] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                        mod_vals_test[1, 0] = mod_vals_test[1, 0] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
                    else:
                        mod_vals_test[0, 1] = mod_vals_test[0, 1] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                        mod_vals_test[1, 1] = mod_vals_test[1, 1] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
                else:
                    if (radii[j] > C.R_CMB + (0.5 * C.R_EARTH/N_shells)) and (radii[j] <= C.R_LMANTLE):
                        mod_vals_test[0, 0] = mod_vals_test[0, 0] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                        mod_vals_test[1, 0] = mod_vals_test[1, 0] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I
                    elif (radii[j] > C.R_LMANTLE):
                        mod_vals_test[0, 1] = mod_vals_test[0, 1] + rhos_mod[j] * ((radius2**3) - (radius1**3)) # M
                        mod_vals_test[1, 1] = mod_vals_test[1, 1] + rhos_mod[j] * ((radius2**5) - (radius1**5)) # I

            mod_vals_test[0, 0] = 1.0e12 * mod_vals_test[0, 0] * 4 * np.pi / 3  # kg
            mod_vals_test[0, 1] = 1.0e12 * mod_vals_test[0, 1] * 4 * np.pi / 3  # kg
            mod_vals_test[1, 0] = 1.0e12 * mod_vals_test[1, 0] * 8 * np.pi / 15 # kg*km^2
            mod_vals_test[1, 1] = 1.0e12 * mod_vals_test[1, 1] * 8 * np.pi / 15 # kg*km^2
            if (core == True):
                abs_err_M_after = np.append(abs_err_M_after, ((mod_vals_test[0, 0]+mod_vals_test[0, 1])-(C.M_EARTH-M_core)))
                abs_err_I_after = np.append(abs_err_I_after, ((mod_vals_test[1, 0]+mod_vals_test[1, 1])-(C.I_EARTH-I_core)))
            else:
                abs_err_M_after = np.append(abs_err_M_after, ((mod_vals_test[0, 0]+mod_vals_test[0, 1])-C.M_EARTH))
                abs_err_I_after = np.append(abs_err_I_after, ((mod_vals_test[1, 0]+mod_vals_test[1, 1])-C.I_EARTH))
 
    if (test == True):
        output = radii, rhos_mod, zoas, coeffs_deep, coeffs_shallow, abs_err_M_before, abs_err_I_before, abs_err_M_after, abs_err_I_after
    else:
        output = radii, rhos_mod, zoas, coeffs_deep, coeffs_shallow
    return output

