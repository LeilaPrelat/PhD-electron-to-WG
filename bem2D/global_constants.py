#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
global constants 
"""
import numpy as np

#%%

def constants():
    """
    Returns
    -------
    hb :  hbar in eV*s
    c :  light velocity in micron/seg
    alpha : fine structure
    me_c2_eV : electron mass * c**2 in eV
    """
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    c =  2.99792458*10**(14)      ### light velocity in micron/seg
    alpha = 1/137.0359            ### fine structure

    ### Ee to velocity : Eduardo Tesis page 43 #### 
    me_c2_eV = 510998.95069  ## eV
    #E0 = me_c2_eV*( np.sqrt( 1/(1-beta**2) ) -1) 
    
    return hb,c,alpha,me_c2_eV

#%%
