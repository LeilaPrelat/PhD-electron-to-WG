#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

convert n, k to epsilon(omega)
plot the permittivities
"""
from scipy.interpolate import interp1d
#from scipy.interpolate import CubicSpline
import numpy as np
import os

#pwd = os.getcwd()
pwd = os.path.dirname(__file__) 
 
 
plot_epsilon = 0

#%%

 
def epsilon(hbw):
    """    
    Parameters
    ----------
    hbw : energy in eV
    delta: extra loss for the imaginary part in eV
    Returns
    -------
    permittivity of Si/Ge from 
    https://refractiveindex.info/?shelf=main&book=Si&page=Aspnes  
    """
    os.chdir(pwd)
    
    table = np.loadtxt('Si_extra_loss.txt',delimiter=' ',skiprows=2)
    table = np.transpose(table)
    eV_list, epsi_real, epsi_imag = table
 
    f_real = interp1d(eV_list, epsi_real)
    f_imag = interp1d(eV_list, epsi_imag)
    
    min_Elist, max_Elist = np.min(eV_list), np.max(eV_list)
    
    #delta = 1e-1
    if hbw<min_Elist:
        rta = f_real(min_Elist) + 1j*(f_imag(min_Elist) )
        
    elif hbw>max_Elist:
        rta = f_real(max_Elist) + 1j*(f_imag(max_Elist) )
    else:
        rta = f_real(hbw) + 1j*(f_imag(hbw) )
    
    # if f_imag(hbw) < 0.1: 
    #     f_imag(hbw) = 0.1
    
    return rta
    
if plot_epsilon == 1:
    
    tamfig = [4, 3]
    tamletra = 13
    tamtitle  = 10
    tamnum = tamletra
    tamlegend = tamletra
    labelpady = 3
    labelpadx = 2
    pad = 2.5
    mk = 1
    ms = 1
    lw = 1.5
    hp = 0.5
    length_marker = 1.5
    dpi = 500
 
    import matplotlib.pyplot as plt
    
    material = 'Si'

    N = int(1e3)
    ev_max = 10
    listx = np.linspace(1, ev_max, N)
    listy_re = []
    listy_im = []
    for x in listx:
        epsi = epsilon(x)
        listy_re.append(np.real(epsi))
        listy_im.append(np.imag(epsi))
    
    ticksx = np.arange(1,ev_max+1,1)
    ticksy = np.arange(-15,45,10)
    title = 'Permittivity of %s' %(material)
#    title = ''
    labelx=r'Electron energy $\hbar\omega$ (eV)'
    labely=r'Permittivity $\epsilon_{\text{%s}}(\omega)$' %(material)
    
    plt.figure(title,figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.plot(listx,listy_re,'-',label = r'Re($\epsilon$)')
    plt.plot(listx,listy_im,'-',label = r'Im($\epsilon$)')
    plt.xticks(ticksx)
    plt.yticks(ticksy)
    plt.xscale('log')
    plt.tick_params(labelsize = tamnum,direction="in", width=1, length = 1.5,pad = 1)
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
    plt.savefig('permittivity_%s.png' %(material), format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)
    # plt.tight_layout()
    plt.show()
    