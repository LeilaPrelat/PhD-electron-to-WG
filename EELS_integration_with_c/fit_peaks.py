## for integration over mode: fitting the peak as a lorentzian and integrate it

import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
 
#%%

def find_width_of_peak_old(listx,listy,index_mode,title1,plot_figure):
    """
    Parameters
    ----------
    listx : x array
    listy : y array
    index_mode : index of the mode
    title1 : info of the parameters for the plot
    plot_figure : if ==1 plots the fit of the curve by the Lorentzian
    Returns
    -------
    width of the peak as a array 
    x_left, x_right from the fitting
    of listy with a lorentzian
    """
    list_energy0 = listx
    list_P_integrated_over_z = listy
    # Lorentzian function
    def lorentzian(x, A, x0, gamma, offset):
        return A * gamma**2 / ((x - x0)**2 + gamma**2) + offset

     
    ######### find peaks and sort them by the highest to minimum (then we identify each mode according to its amplitude)
    peaks, _ = find_peaks(list_P_integrated_over_z, height=0 )
    listy_peaks = []
    listx_peaks = []
     
    for peak in peaks: 
        listy_peaks.append(list_P_integrated_over_z[peak])
        listx_peaks.append(list_energy0[peak])

    ## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
    sorted_index = np.argsort(listy_peaks)
    # listy_peaks_sorted = np.sort(listy_peaks) 
    # listx_peaks_sorted = list_energy0[peaks[sorted_index]]
    list_index_peaks = peaks[sorted_index]

    ind_max = list_index_peaks[index_mode]

    ind0 = int(ind_max*0.5)
    ind1 = int(ind_max*1.5)
    x_data = list_energy0[ind0:ind1] ## we fit around the maximum
    y_data = list_P_integrated_over_z[ind0:ind1]
    
    ############ good estimation for gamma #########################################
    # Get the half-maximum value
    half_max = (max(y_data) + min(y_data)) / 2
    
    # Find indices where the curve crosses the half-maximum
    above_half_max = y_data > half_max
    crossings = np.where(np.diff(above_half_max.astype(int)) != 0)[0]
    
    if len(crossings) >= 2:
        # Estimate FWHM from the distance between crossings
        x1 = x_data[crossings[0]]
        x2 = x_data[crossings[-1]]
        fwhm_estimate = abs(x2 - x1)
        gamma_estimate = fwhm_estimate / 2  # Lorentzian has FWHM = 2 * gamma
    else:
        # Fallback if we don't have clear half-max crossings
        gamma_estimate = 0.1 * (x_data[-1] - x_data[0])  # 10% of window as fallback
    
    ###############################################################################
    
    # Initial parameter guesses: [A, x0, gamma, offset]
    initial_guess = [max(y_data), list_energy0[ind_max], gamma_estimate, min(y_data)]
    
    fit_success = False  # Flag to track whether fitting was successful

    try: 
        # Perform the fit
        popt, pcov = curve_fit(lorentzian, x_data, y_data, p0=initial_guess)
        A_fit, x0_fit, gamma_fit, offset_fit = popt
        
        # Threshold at 1% of amplitude
        threshold = offset_fit + 0.01 * A_fit
        
        # Solve for x where Lorentzian equals threshold
        delta = np.sqrt((A_fit * gamma_fit**2) / (threshold - offset_fit) - gamma_fit**2)
        x_left = x0_fit - delta
        x_right = x0_fit + delta
        
        fit_success = True  # Fit completed successfully
        

    except RuntimeError as error:
        print(error)
        x_left, x_right = list_energy0[ind0],list_energy0[ind1]
    
    step2=0.001 #integration over freq
    x_integrate = np.arange(x_left, x_right + step2, step2)
    function_lorentzian = np.sum(lorentzian(x_integrate, A_fit, x0_fit, gamma_fit, offset_fit))*step2
 
    # Plot regardless of fitting result
    if plot_figure == 1:
        labelx = 'Electron energy loss $\hbar\omega$ (eV)'
        labely = r'$\Gamma(\omega)$ (1/eV)'
    
        tamfig = [4.5, 3.5]
        tamletra = 13
        tamtitle = tamletra - 3
        tamnum = tamletra
        labelpady = 3
        labelpadx = 2
        pad = 2.5
        dpi = 500
    
        #title2 = r'$V_0$ = %.3f eV, $\theta$ = %.2f mrad' % (V0, theta_mrad)
    
        plt.figure(figsize=tamfig)
        plt.title(title1, fontsize=tamtitle)
        plt.xlabel(labelx, fontsize=tamletra, labelpad=labelpadx)
        plt.ylabel(labely, fontsize=tamletra, labelpad=labelpady)
        plt.plot(list_energy0, list_P_integrated_over_z, '.-')
    
        if fit_success:
            plt.plot(x_integrate, lorentzian(x_integrate, *popt), 'r-',label = 'Lorentzian fit')
            plt.plot([x_left], [0], "x", color='black')
            plt.plot([x_right], [0], "x", color='black')
    
        plt.xticks(np.arange(list_energy0[0], list_energy0[-1]+0.05, 0.05))
        plt.tick_params(labelsize=tamnum, length=2, width=1, direction="in", which='both', pad=pad)
        plt.legend(loc = 'best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.2, handlelength=1)
        plt.savefig('EELS_integrated_over_z_fit_width.png', bbox_inches='tight', pad_inches=0.01, format='png', dpi=dpi)
        plt.savefig('EELS_integrated_over_z_fit_width.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf', dpi=dpi)
    
    

    # bottom_width = x_right - x_left
    

    return x_left, x_right, x0_fit, function_lorentzian, fit_success





def find_width_of_peak(listx1, listy1, index_mode, title1, plot_figure,
                       listx2=None, listy2=None, title2=None):
    """
    Parameters
    ----------
    listx1, listy1 : arrays
        First dataset (x, y).
    index_mode : int
        Index of the mode (peak).
    title1 : str
        Info of the parameters for the plot (first dataset).
    plot_figure : int
        If ==1, plots the fit(s) of the curve(s) by the Lorentzian.
    listx2, listy2 : arrays, optional
        Second dataset (x, y) to compare with the first.
    title2 : str, optional
        Label/title for the second dataset.
    Returns
    -------
    results : dict
        Keys: 'dataset1', 'dataset2' (if provided).
        Each contains (x_left, x_right, x0_fit, lorentzian_integral, fit_success).
    """

    def lorentzian(x, A, x0, gamma, offset):
        return A * gamma**2 / ((x - x0)**2 + gamma**2) + offset

    def process_dataset(listx, listy, index_mode):
        # --- Find peaks ---
        peaks, _ = find_peaks(listy, height=0)
        listy_peaks = listy[peaks]
        sorted_index = np.argsort(listy_peaks)
        list_index_peaks = peaks[sorted_index]

        ind_max = list_index_peaks[index_mode]
        ind0 = int(ind_max*0.5)
        ind1 = int(ind_max*1.5)
        x_data = listx[ind0:ind1]
        y_data = listy[ind0:ind1]

        # --- Estimate gamma from half max ---
        half_max = (max(y_data) + min(y_data)) / 2
        above_half_max = y_data > half_max
        crossings = np.where(np.diff(above_half_max.astype(int)) != 0)[0]

        if len(crossings) >= 2:
            x1 = x_data[crossings[0]]
            x2 = x_data[crossings[-1]]
            fwhm_estimate = abs(x2 - x1)
            gamma_estimate = fwhm_estimate / 2
        else:
            gamma_estimate = 0.1 * (x_data[-1] - x_data[0])

        initial_guess = [max(y_data), listx[ind_max], gamma_estimate, min(y_data)]
        fit_success = False

        try:
            popt, _ = curve_fit(lorentzian, x_data, y_data, p0=initial_guess)
            A_fit, x0_fit, gamma_fit, offset_fit = popt
            threshold = offset_fit + 0.01 * A_fit
            delta = np.sqrt((A_fit * gamma_fit**2) / (threshold - offset_fit) - gamma_fit**2)
            x_left = x0_fit - delta
            x_right = x0_fit + delta
            fit_success = True
        except RuntimeError:
            x_left, x_right = listx[ind0], listx[ind1]
            x0_fit, A_fit, gamma_fit, offset_fit = listx[ind_max], 0, 1, 0

        step2 = 0.001
        x_integrate = np.arange(x_left, x_right + step2, step2)
        lorentz_integral = np.sum(lorentzian(x_integrate, A_fit, x0_fit, gamma_fit, offset_fit)) * step2

        return (x_left, x_right, x0_fit, gamma_fit, A_fit, lorentz_integral, fit_success, (x_integrate, popt))

    # --- Process dataset 1 ---
    results = {}
    results['dataset1'] = process_dataset(listx1, listy1, index_mode)
 
    x_left, x_right, x0_fit, gamma_fit, A_fit, lorentzian_integral, fit_success, (x_integrate, popt1) = results['dataset1']

    # --- Process dataset 2 if provided ---
    if listx2 is not None and listy2 is not None:
        results['dataset2'] = process_dataset(listx2, listy2, index_mode)
        x_left2, x_right2, x0_fit2, gamma_fit2, A_fit2, lorentzian_integral2, fit_success2, (x_integrate2, popt2) = results['dataset2']

        
    # --- Plot ---
    if plot_figure == 1:
        
        tamfig = [4,3]
        tamletra = 14
        tamtitle  = 8
        tamnum = tamletra
        tamlegend = tamletra
        labelpady = 2
        labelpadx = 3
        pad = 2.5
        dpi = 300
        
        plt.figure(figsize=tamfig)
        plt.title(title1 if title2 is None else f"{title1} vs {title2}", fontsize=tamtitle)
        plt.xlabel(r'Electron energy loss $\hbar\omega$ (eV)', fontsize=tamletra)
        plt.ylabel(r'$\Gamma(\omega,x)$ (1/eV)', fontsize=tamletra)

        # Plot dataset 1
        x_left, x_right, x0_fit, gamma_fit, A_fit, _, fit_success, (x_integrate, popt1) = results['dataset1']

        # Plot dataset 2 (if any)
        if 'dataset2' in results:
            x_left, x_right, x0_fit, gamma_fit, A_fit, _, fit_success, (x_integrate, popt2) = results['dataset2']
            plt.plot(listx2, listy2, '--',color = 'blue', label=r"$\phi$ linear")
            if fit_success:
                plt.plot(x_integrate, lorentzian(x_integrate, *popt2), 'r--')
                energy_compare = 0.8 ## energy used to compared both cases (approx and no approx)
                plt.plot(energy_compare, lorentzian(energy_compare, *popt1),'o',color = 'black')
                #plt.plot([x_left, x_right], [0, 0], "x", color='blue')
        
        plt.plot(listx1, listy1, '.-',color = 'blue', label=r"$\phi$ from BEM")
        if fit_success:
            plt.plot(x_integrate, lorentzian(x_integrate, *popt1), 'r-', label='Fit')
            #plt.plot([x_left, x_right], [0, 0], "x", color='black')
            
        plt.xticks(np.arange(listx1[0], listx1[-1]+0.05, 0.05))
        plt.tick_params(labelsize=tamnum, length=2, width=1, direction="in", which='both', pad=pad)
        plt.legend(loc = 'upper left',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.2, handlelength=1)
        plt.savefig('EELS_integrated_over_z_fit_width_comparison.png',bbox_inches='tight',pad_inches = 0.01, dpi=dpi)
        plt.savefig('EELS_integrated_over_z_fit_width_comparison.pdf',bbox_inches='tight',pad_inches = 0.01, dpi=dpi)

    return results
















