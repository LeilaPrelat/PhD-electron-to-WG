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

## fixed x         wp_norm_w,d_norm_w,h_norm_w,N,epsilon,0,h_norm_w,h_norm_w,2
def run_e_out_EELS(wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1): 
    x1 = x0
    nx = 1
    nz = 300
    
    ### all the parameters are normalized to width!!! and BEM files NOT ###
    
    os.chdir(path_basic)
    # Run the C++ program name dy.out and save a list of y and V
    #absolute_path_windows = r"E:\Desktop\Leila\EELs_omega_vs_theta_new\PhD-electron-to-WG-main\potential\dy.exe"
    
    ## WINDOWS ####
    ## for rocket --> parallelization (rochet). we need the absolute path for windows and we need to compile it in windows
    ## terminal as well: g++ .\dy.cpp -o dy.out
    
    cmd = ["./e.out", str(wp_norm_w), str(d_norm_w), str(h_norm_w), str(N), str(epsilon), str(x0), str(x1), str(nx), str(z0), str(z1), str(nz)]
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = result.stdout.strip().split('\n')
    
    list_z_norm_w = []
    list_V_norm_U = []
    for line in lines:
        try:
            zvalue, Vvalue = map(float, line.strip().split())
            list_z_norm_w.append(zvalue)
            list_V_norm_U.append(Vvalue)  
        except ValueError:
            continue
        
        
    # zz, V0_over_U_list = run_e_out_EELS(wp_norm_w,d_norm_w,h_norm_w,N,epsilon,0,h_norm_w,h_norm_w,2)
    # V0_over_U = np.abs(V0_over_U_list[0]) ## normalized to V0


    return list_z_norm_w, list_V_norm_U


# dy_cache = {}

# def dy_cached(wp,d,h,N,epsilon,x0,z0,z1,nz):
#     if d not in dy_cache:
#         print('running dy.out for d/W = %i' %(d))
#         dy_cache[d] = run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz)
#     return dy_cache[d]

def EELS_from_BEM_interpolated(energy_eV,w,h,xe,Eelectron_keV,N):
    """
    Parameters
    ----------
    energy_eV: photon energy hb*omega in eV
    w : width along x (nm)
    h : thickness along z in nm
    xe : position of the electron respect to the width 
    Eelectron_keV : energy of the electron in keV
    N : number of points in bem2d
    Returns
    -------
    EELS with ss fixed to 10% of width
    """
    os.chdir(path_data)
    # energy_label = 'energy%.3feV' %(energy_eV)
    # pp = energy_label.replace('.','')
    ze = 50  ## mistake but the EELS is along z, so ze is not fixed
    
    if xe == 0.4: ## 40%
        labelxe ='_xe40W'
    elif xe == 0.8:  ## 80%
        labelxe ='_xe80W'
    else: ## center 
        labelxe ='_xe0'
    
    # name2 = 'EELS_along_z_N%i_a%inm_h%inm_ze%inm_Ee%ikeV_%s'%(N,a,h,ze,Eelectron_keV,labelxe) + '_' + pp + '.dat'
    name2 = 'EELS_along_z_Si_N%i_a%inm_h%inm_ze%inm_Ee%ikeV%s_energy_%.7f'%(N,w,h,ze,Eelectron_keV,labelxe,energy_eV) + '.txt'
    tabla1 = np.loadtxt(name2, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
    
    listeV = tabla1_2[0]
    # listq = tabla1_2[1]
    # listx = tabla1_2[2]
    listz = tabla1_2[3]     ##
    listEELS = tabla1_2[4]
    
    listeV_correct = []
    listEELS_correct = []
    listz_correct = []
    for k in range(len(listEELS)):
        if listEELS[k]>0: 
            
            listeV_correct.append(listeV[k])
            listEELS_correct.append(listEELS[k])
            listz_correct.append(listz[k])
        # else:
        #     print(k,listEELS[k])
    
    listz_norm_w_BEM = np.array(listz_correct)/w   ## we interpolate using z/W, same axis as the potential code
    
    # listz_norm_a_BEM = np.array(listz_norm_a_BEM).flatten() ## force arrays to be 1D
    # listEELS_correct = np.array(listEELS_correct).flatten() ## force arrays to be 1D
    
    # print("First 10 z:", listz_norm_a_BEM[:10])
    # print("First 10 EELS:", listEELS_correct[:10])
    
    
    EELS_interp =  interp1d(listz_norm_w_BEM, listEELS_correct)
    
    
    # print(len(listz_norm_a_BEM),len(listEELS_correct), energy_eV)
 
    return listz_norm_w_BEM, EELS_interp

 
def P_integrand_over_z(value_z_norm_w, energy_eV, w, h, xe, Eelectron_keV, N, theta, V0, list_z_norm_w, listV_normU, V0_over_U):
    """
    Parameters
    ----------
    value_z_norm_w : z/w coordinate normalized to w
    energy_eV: photon energy hb*omega in eV
    w : width along x (nm)
    h : heigh (nm)
    xe : position of the electron respect to the width 
    Eelectron_keV : energy of the electron in keV
    N : number of points in bem2d
    list_z_norm_w : list of z/width from c++ code.this array depends on wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1
    listV_normU : list of V(z)/V0 from c++ code. this array depends on wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1
    V0_over_U : V(x=0,z=h)=V0/U we use it to normalized listV_normU and obtain listV_normV0
    Returns
    -------
    P(omega) integrand over z/a 
    see "Eq_P_integral_over_z.png"
    """
    epsi1 = 1
    me_c2_eV = 510998.95069  ## me*c**2 in eV
    Eelectron = Eelectron_keV*1e3
    beta = np.sqrt( 1- (1 + Eelectron/me_c2_eV)**(-2) )  ## beta = v/c
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    aux_function = beta*np.sin(theta)
    # print(np.sin(theta))
    # aux_function = beta*theta
    
    # list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N)
    # list_z_norm_a_interp = np.linspace(1.1,2,int(N*Nint)) 
    # print(len(list_z_norm_a),len(listV_normV0))
    V_interp =  interp1d(list_z_norm_w, listV_normU)
    ## NOW THE POTENTIAL IS NEGATIVE BY DEFINITION ## 
    V_norm_U = V_interp(value_z_norm_w)
    
    # print(V0_over_U)
    
    V_norm_V0 = V_norm_U/V0_over_U
    
    
    # h_norm_w = h/w
    # V_norm_U_0 = V_interp(h_norm_w)
    
    
    V_norm_V0_zmin = V_interp(0.7)/V0_over_U
    
 
    ## we force the potential to be negative to have the electron attraction 
    ## WARNING: some values of z will make the square root in the denominator complex!! #issue04 --> add a modulus inside the sqrt
    arg_denominator = np.abs(aux_function**2 + 2*V_norm_V0*V0/(me_c2_eV*gamma_e))
    denominator = np.sqrt(arg_denominator)
    # print(arg_denominator)
 
    listz_norm_a_BEM, EELS_interp = EELS_from_BEM_interpolated(energy_eV,w,h,xe,Eelectron_keV,N)

    EELS_value = EELS_interp(value_z_norm_w) 
    
    function = 2*beta*EELS_value/denominator
    # print(arg_denominator, denominator, value_z_norm_w, EELS_value, V_norm_V0, function)
    # print(V_norm_V0,EELS_value)
    
    return function ## we add "a" because the integral is going to be over z/a


## integration over z as a sum 
def P_integrated_over_z(z_min_val, energy_eV, w, h, xe, Eelectron_keV, N, theta, V0, list_z_norm_w, listV_normU, V0_over_U):
    """
    Parameters
    ----------
    z_min_val : minimum value of the integration over z/W to infty
    theta : angle of the electron (initial conditions) in radians
    V0 : potential of the gate in e*Volts
    xe : position of the electron respect to the width 
    Eelectron : energy of the electron in keV
    bb : aspect ratio between width and height
    dd : distance to the plane normalized to width
    energy_eV: photon energy hb*omega in eV
    a : width along x (nm)
    N : number of points in bem2d
    list_z_norm_a : list of z/width from c++ code 
    listV_normU : list of V(z)/V0 from c++ code. this array depends on wp_norm_w,d_norm_w,h_norm_w,N,epsilon,x0,z0,z1
    V0_over_U : V(x=0,z=h)=V0/U we use it to normalized listV_normU and obtain listV_normV0
    Returns
    -------
    P(omega) integrated over z/a 
    see "Eq_P_integral_over_z.png"
    """
    Nint = 10
 
 

    listz_norm_w_BEM, EELS_interp = EELS_from_BEM_interpolated(energy_eV,w,h,xe,Eelectron_keV,N)
    listz_norm_w_BEM_interp = np.linspace(z_min_val, np.max(listz_norm_w_BEM),int(N*Nint)) 
    
    
    P_array = P_integrand_over_z(listz_norm_w_BEM_interp, energy_eV, w, h, xe, Eelectron_keV, N, theta, V0, list_z_norm_w, listV_normU, V0_over_U )
    list_P = np.sum(P_array)
    # list_P = 0
    # for value_z_norm_a in listz_norm_a_BEM_interp: ## the integration should be from z_min 
    #     P_value = P_integrand_over_z(value_z_norm_a, theta, V0, Ee_electron, bb, ss, dd, energy_eV, a, N, list_z_norm_a, listV_normV0)
    #    #if np.isnan(P_value)==True: ## the square root defined in the function integrand is imaginary
    #    #     list_P = list_P + 0
    #    # else:
    #     list_P = list_P + P_value
        # print(P_value)
        
        
    delta_z_norm_w = listz_norm_w_BEM_interp[1] - listz_norm_w_BEM_interp[0]

    return list_P*delta_z_norm_w*w, listz_norm_w_BEM_interp, P_array


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


 