#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 16:03:14 2025

@author: lprelat
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import matplotlib as mpl
import os
import sys
 
def find_theta_V0(w, wp_norm_w, d_norm_w, h_norm_w, N, epsilon, x0_norm_w, Ee_electron_keV, bmin_norm_w):
    """
    Parameters
    ----------
    w : width of the middle waveguide in nm
    wp_norm_w : wp/w, wp width of the side waveguides
    d_norm_w : d/w, d distance between the side waveguides and the central one
    h_norm_w : h/w, h height of the waveguides
    N : discretization points
    epsilon : epsilon of the substrate
    x0_norm_w : xi/w = xf/w and nx = 1
    Ee_electron_keV : energy of electron in keV
    bmin_norm_w  : bmin/w where bmin starts from 0 (we substract the h)

    Returns
    -------
    gives the solutions of theta_mrad, V0 (eV)
    for the corresponding bmin/w
    """
    
    path_basic = os.getcwd()
    path_data = os.path.join(path_basic, 'potential_data')
    os.makedirs(path_data, exist_ok=True)
    os.chdir(path_data)

         
    me_c2_eV = 510998.95069
    Ee_electron = Ee_electron_keV * 1e3
    
    name_potential = 'potential_wp%.1f_d%.1f_h%.1f_N%i_xe%.1f.txt' %(wp_norm_w,d_norm_w,h_norm_w,N,x0_norm_w)
    try: 
        tabla_sh = np.loadtxt(name_potential,skiprows=1,delimiter=' ') 
    except FileNotFoundError:
        raise('Data',tabla_sh,'not found. Run run_potential_norm_V0.sh' )
        
    tabla_sh2 = np.transpose(tabla_sh)
    list_z_norm_w_sh, listV_normV0_sh = tabla_sh2
     
    V_interp = CubicSpline(list_z_norm_w_sh, listV_normV0_sh)

    def function_to_be_zero(value_z_norm_w, theta, V0):
        beta = np.sqrt(1 - (1 + Ee_electron / me_c2_eV) ** -2)
        gamma_e = 1 / np.sqrt(1 - beta ** 2)
        aux_function = beta * np.sin(theta)
        function = aux_function ** 2 * me_c2_eV * gamma_e / 2
        V_norm_V0 = V_interp(value_z_norm_w)
        partA = function / V0
        partB = V_norm_V0
        return partA + partB

    # if Ee_electron_keV == 200:
    #     Nvals = 300
    #     theta_mrad_vals = np.linspace(0.01, 2.25, Nvals)
    #     V0_vals = np.linspace(0.01, 1, Nvals)
    # elif Ee_electron_keV == 100:
    #     Nvals = 300
    #     theta_mrad_vals = np.linspace(0.01, 3, Nvals)
    #     V0_vals = np.linspace(0.01, 1, Nvals)
    # else:
    #     print("Unsupported Ee_electron_keV value.")
    #     sys.exit(1)
        
    Nvals = 300   
    theta_mrad_vals = np.linspace(0.01, 3, Nvals)
    V0_vals = np.linspace(0.01, 1, Nvals)

    zmin = h_norm_w
    zmax = np.max(list_z_norm_w_sh)

    X_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
    for i, V0 in enumerate(V0_vals):
        for j, theta_mrad in enumerate(theta_mrad_vals):
            theta = theta_mrad * 1e-3
            try:
                x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta, V0))
                X_vals[i, j] = x_root - h_norm_w
            except ValueError:
                X_vals[i, j] = np.nan

    X_vals_transpose = np.transpose(X_vals)
    
    tamfig2 = [5, 3]
    tamletra = 13
    tamtitle  = tamletra - 5
    tamnum = tamletra
    labelpady = 3
    labelpadx = 2
    pad = 2.5
    dpi = 500
    
    limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
    cmap = plt.cm.RdBu   # define the colormap

    # another cbar with constant a*1e-3 (real units, in microns)
    cte_cbar2 = w*1e-3
    # Create contour plot
    contour_levels = [bmin_norm_w]

    vmin1 , vmax1 = np.nanmin(X_vals), np.nanmax(X_vals)
    bounds1 =  [vmin1,0.1,bmin_norm_w,0.5,1]
    ## bounds before was np.log10(0.1) was too big 
    bounds1 =   np.logspace(np.log10(0.05), np.log10(vmax1) , 10) 
    # bounds1 =   np.logspace(np.log10(vmin1), np.log10(100) , 10) 
    norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)
    norm2 = mpl.colors.BoundaryNorm(bounds1*cte_cbar2, cmap.N) ## second colorbar with real units 

    contours = plt.contour(V0_vals, theta_mrad_vals, X_vals_transpose, levels=contour_levels,lw = 0, colors='white', linestyles='dashed' )

    tol = 1*1e-1 ## zero of the function 
    allsegs = contours.collections[0].get_paths()  # Assuming only one contour collection
    
    # Check the contour data is not misleading : not corresponding to NaN values  
    valid_paths = []
    for path in allsegs:
        vertices = path.vertices 
    
        for v0, theta_mrad in vertices:
     
            theta = theta_mrad*1e-3
            value_z_norm_w = bmin_norm_w + h_norm_w
            value_z = function_to_be_zero(value_z_norm_w, theta, v0)
            #print(v0,theta_mrad,value_z)
            if not np.isnan(value_z) and np.abs(value_z)<tol: ## filter the extracted V0,theta to only the ones that give real results 
                valid_paths.append([v0, theta_mrad])
    
    title = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f' % (wp_norm_w,h_norm_w,d_norm_w,x0_norm_w)
    if len(valid_paths)>1:
        
        list_V0_fix_bmin = []
        list_theta_fix_bmin = [] 
        for xy in valid_paths:
            x,y = xy
            list_theta_fix_bmin.append(y)
            list_V0_fix_bmin.append(x)
            
        #print(np.min(list_V0_fix_bmin))
        
    
        # Subtract the V0 and \theta values of corresponding points on different levels
        # Sorte the lists argsort()
        index_sorted = np.argsort(list_V0_fix_bmin)
        listx_sorted = np.sort(list_V0_fix_bmin)
        listy_sorted = [] 
        for index in index_sorted:
            listy_sorted.append(list_theta_fix_bmin[index])
        
        header = title + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV)
        label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,x0_norm_w) ## all normalize to width w
        label_Ee = '_Ee%ikeV' % (Ee_electron_keV)
        ## saving data for plotting ##
        np.savetxt('bmin_for_plot' + label_Ee + label_txt + '.txt', X_vals_transpose, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('V0_for_plot_for_bmin%.2f' %(bmin_norm_w) + label_Ee + label_txt + '.txt', V0_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('theta_mrad_for_plot_for_bmin%.2f' %(bmin_norm_w) + label_Ee + label_txt + '.txt', theta_mrad_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    
        ## saving solutions (V0,theta) such that bmin = bmin_vals0
        name_sol = 'V0_theta_mrad_sols_for_bmin%.2f' %(bmin_norm_w) + label_Ee + label_txt + '.txt'
        table_solutions = np.transpose([listx_sorted, listy_sorted])
        np.savetxt(name_sol, table_solutions, fmt='%.10f', delimiter='\t', header = 'V0 (eV)        theta (mrad), ' + header, encoding=None)
        
        print("Processing complete. Data saved.")
    
    else:
        print("No valid (V0, theta) solutions found. Skipping plot.")

    
    plt.figure(figsize=tamfig2)
    plt.title(title,fontsize=tamtitle-1)
    #plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
    im_show = plt.imshow(X_vals_transpose, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' , norm = norm1  ) 
        #plt.clabel(contours, fmt='%.2f', colors='green', fontsize=tamletra, manual=[(0.5, 5) ])  # Label contours
    cbar = plt.colorbar(im_show, fraction=0.046, pad=0.14 , format = '%.2f') 
    im_show2 = plt.imshow(np.array(X_vals_transpose)*cte_cbar2, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm=norm2  )  ## second colorbar with real units 
    cbar2 = plt.colorbar(im_show2, fraction=0.046, pad=0.04, orientation = 'vertical' , format = '%.2f')
    cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra-1)
    cbar.ax.tick_params(labelsize = tamnum, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
    cbar2.ax.tick_params(labelsize = tamnum, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
    cbar2.ax.set_title(r'$b_{\text{min}}$ ($\mu$m)',fontsize=tamletra-1)
    plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
    #plt.plot(listx_sorted,listy_sorted,'--',color = 'green')
 
    plt.text(0.55, 0.12, r"$E_{\text{e}} = %i$ keV" %(Ee_electron_keV),color = 'black',fontsize = tamletra)
    plt.text(0.55, 0.42, r"$x = %i$ nm" %(x0_norm_w*w),color = 'black',fontsize = tamletra)
    plt.savefig('zmin' + label_Ee + label_txt + '.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
    plt.savefig('zmin' + label_Ee + label_txt + '.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
    plt.show()
        
    print("Plot saved.")

if __name__ == "__main__":
    if len(sys.argv) != 10:
        print("Usage: python find_theta_V0.py w wp/w d/w h/w N epsilon x0/w Ee_electron_keV bmin/w")
        sys.exit(1)

    args = [float(arg) for arg in sys.argv[1:]]

    labels = ["w (nm)", "wp/w", "d/w", "h/w", "N", "epsilon", "x0/w", "Ee_electron_keV", "bmin/w"]
    print("Input Parameters:")
    for label, value in zip(labels, args):
        print(f"{label}: {value}")

    find_theta_V0(*args)
