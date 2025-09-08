#code d.h : 
#- Models a grounded plane + biased rectangular wire with rounded corners.
#- Uses the boundary element method to compute induced surface charge.
#- Computes electrostatic potential/V0 at any point in space (x,y)
# where: 
#bb = 2, aspect ratio  (side length along y (labelled as b) divided by side length along x).
#ss = 0.1, rounding radius.
#dd = 1, distance to the plane.
#N = 200, discretization points.
### All lengths are given in units of the x side-length, label as "a"

## Define the EELS integrated over z and over the peaks
## Run the code "EELS_rewrite_name_files.py" inside bem2D folder to separate the files in each frequency

import subprocess
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
 
path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem2D_files_along_z')

#%%

# def P_integrand_over_z(z_norm_w, energy, W, h, Eelectron_keV, theta, V0, list_z_norm_w, listV_normV0, mode):
#     """
#     Parameters
#     ----------
#     z_norm_w : z/w coordinate normalized to w
#     energy_eV: photon energy hb*omega in eV
#     W : width along x (nm)
#     h : heigh (nm)
#     Eelectron_keV : energy of the electron in keV
#     theta : angle of incidence in mrad
#     V0 : voltage of surface of the wg in eV
#     list_z_norm_w : list of z/width from c++ code.this array depends on wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1
#     listV_normV0 : list of V(z)/V0 from c++ code. this array depends on wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1
#     mode: mode of the EELS. we sort them by amplituted 
#     Returns
#     -------
#     P(omega) integrand over z/a 
#     see "Eq_P_integral_over_z.png"
#     """
#     epsi1 = 1
#     me_c2_eV = 510998.95069  ## me*c**2 in eV
#     Eelectron = Eelectron_keV*1e3
#     beta = np.sqrt( 1- (1 + Eelectron/me_c2_eV)**(-2) )  ## beta = v/c
#     gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
#     aux_function = beta*np.sin(theta)
#     h_norm_W = h/W
    
#     V_interp =  interp1d(list_z_norm_w, listV_normV0)
#     V_norm_V0 = V_interp(z_norm_w)    
 
#     N = 700
#     info = 'w%inm_h%inm_N%i_Ee%ikeV_mode%i' %(W,h,N,Eelectron_keV,mode)
#     name_file = "fitted_parameters_%s.txt" %(info)
    
#     table = np.loadtxt(name_file,  delimiter='\t', skiprows = 1, encoding=None)
#     table2 = np.transpose(table)
#     list_energy, list_A_fit, list_k_par_times_W_fit = table2   
 
 
#     list_integrand  = []
    
#     for j in range(len(list_energy)): 
#         k_par_times_w = list_k_par_times_W_fit[j]
#         amp = list_A_fit[j]
        
#         fit = amp*np.exp(-2*k_par_times_w*(z_norm_w-h_norm_W)) 
        
#         ## we force the potential to be negative to have the electron attraction 
#         ## WARNING: some values of z will make the square root in the denominator complex!! #issue04 --> add a modulus inside the sqrt
#         arg_denominator = np.abs(aux_function**2 + 2*V_norm_V0*V0/(me_c2_eV*gamma_e))
#         denominator = np.sqrt(arg_denominator)
#         # print(arg_denominator)
        
#         function = 2*beta*fit/denominator
#         list_integrand.append(function)
        
    
#     EELS_interp =  interp1d(list_energy, list_integrand)
    
#     return EELS_interp(energy)  


## integration over z as a sum 
def P_integrated_over_z(z_min_val, energy, W, h, Eelectron_keV, theta, V0, list_z_norm_w, listV_normV0, mode):
    """
    Parameters
    ----------
    z_min_val : minimum value of the integration over z/W to infty
    theta : angle of the electron (initial conditions) in radians
    V0 : potential of the gate in e*Volts
    Eelectron : energy of the electron in keV
    bb : aspect ratio between width and height
    dd : distance to the plane normalized to width
    energy_eV: photon energy hb*omega in eV
    a : width along x (nm)
    N : number of points in bem2d
    list_z_norm_a : list of z/width from c++ code 
    listV_normV0: list of V(z)/V0 from c++ code. this array depends on wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1
    Returns
    -------
    P(omega) integrated over z/a 
    see "Eq_P_integral_over_z.png"
    """


    epsi1 = 1
    me_c2_eV = 510998.95069  ## me*c**2 in eV
    Eelectron = Eelectron_keV*1e3
    beta = np.sqrt( 1- (1 + Eelectron/me_c2_eV)**(-2) )  ## beta = v/c
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    aux_function = beta*np.sin(theta)
    h_norm_W = h/W

 
    N = 500
    info = 'w%inm_h%inm_N%i_Ee%ikeV_mode%i' %(W,h,N,Eelectron_keV,mode)
    name_file = "fitted_parameters_%s.txt" %(info)
    
    table = np.loadtxt(name_file,  delimiter='\t', skiprows = 1, encoding=None)
    table2 = np.transpose(table)
    list_energy, list_A_fit, list_k_par_times_W_fit = table2   
    
    A_interp =  interp1d(list_energy, list_A_fit)
    k_par_times_W_interp = interp1d(list_energy, list_k_par_times_W_fit)
    
 
    V_interp =  interp1d(list_z_norm_w, listV_normV0)
    
    z_max_val = 600/W ## EELS in bem2d goes up to z = 600 nm
    
    Nz = int(1e3)
    z_vals = np.linspace(z_min_val,z_max_val,Nz)
    
    k_par_times_w =  k_par_times_W_interp(energy)
    amp = A_interp(energy)    
    fit = amp * np.exp(-2 * k_par_times_w * (z_vals - h_norm_W))
    ## EELS ##########
    
    V_norm_V0 = V_interp(z_vals)
    arg_denominator = np.abs(aux_function**2 + 2 * V_norm_V0 * V0 / (me_c2_eV * gamma_e))
    denominator = np.sqrt(arg_denominator)
    
    integrand = 2 * beta * fit / denominator

# --- integrate
    deltaz = z_vals[1] - z_vals[0]
    integral = np.sum(integrand) * deltaz * W 
    ## we add "a" because the integral is going to be over z/a
    
    # integral = np.trapz(list_integrand, list_z_norm_w)*W
    
    
    return integral

def find_width_of_peak(listx,listy,index_mode,title1,plot_figure):
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

    ind0 = int(ind_max*0.8)
    ind1 = int(ind_max*1.1)
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
        

    except RuntimeError or IndexError as error:
        print(error)
        x_left, x_right = list_energy0[ind0],list_energy0[ind1]
    
    step2=0.001 #integration over freq
    x_integrate = np.arange(x_left, x_right + step2, step2)
    function_lorentzian = np.sum(lorentzian(x_integrate, A_fit, x0_fit, gamma_fit, offset_fit))*step2
    
    # Plot regardless of fitting result
    if plot_figure == 1:
        labelx = 'Electron energy loss $\hbar\omega$ (eV)'
        labely = 'EELS per electron (1/eV)'
    
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
            plt.plot(x_data, lorentzian(x_data, *popt), 'r-')
            plt.plot([x_left], [0], "x", color='black')
            plt.plot([x_right], [0], "x", color='black')
    
        # plt.xticks(np.arange(0.5, 3.5, 0.5))
        plt.tick_params(labelsize=tamnum, length=2, width=1, direction="in", which='both', pad=pad)
        plt.savefig('EELS_integrated_over_z_fit_width.png', bbox_inches='tight', pad_inches=0.01, format='png', dpi=dpi)
    
    

    # bottom_width = x_right - x_left
    

    return x_left, x_right, x0_fit, function_lorentzian, fit_success


 