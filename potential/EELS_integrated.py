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

## plot the z min as a function of theta, V0
import subprocess
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem_files_EELS')

#%%

def run_dy_out(bb,ss,dd, N):
    os.chdir(path_basic)
    # Run the C++ program name dy.out and save a list of y and V
    absolute_path_windows = r"E:\Desktop\Leila\EELs_omega_vs_theta\PhD-electron-to-WG-main\PhD-electron-to-WG-main\potential\./dy.out" ## for rocket --> parallelization
    cmd = [  absolute_path_windows , str(bb), str(ss), str(dd), str(N)] 
#    cmd = ["./dy.out", str(bb), str(ss), str(dd), str(N)]
 
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = result.stdout.strip().split('\n')
    list_z = []
    list_V = []
    for line in lines:
        try:
            zvalue, Vvalue = map(float, line.strip().split())
            list_z.append(zvalue)
            list_V.append(Vvalue)  
        except ValueError:
            continue

    return list_z, list_V

dy_cache = {}

def dy_cached(bb,ss,dd, N):
    if dd not in dy_cache:
        print('running dy.out for dd = %i' %(dd))
        dy_cache[dd] = run_dy_out(bb,ss,dd, N)
    return dy_cache[dd]

def EELS_from_BEM_interpolated(energy_eV,a,h,N):
    """
    Parameters
    ----------
    energy_eV: photon energy hb*omega in eV
    a : width along x (nm)
    h : thickness along z in nm
    N : number of points in bem2d
    Returns
    -------
    EELS
    """
    os.chdir(path_data)
    energy_label = 'energy%.3feV' %(energy_eV)
    pp = energy_label.replace('.','')
    name1 = 'EELS_along_z_N%i_a%inm_h%inm'%(N,a,h) + '_' + pp + '.dat'
    tabla1 = np.loadtxt(name1, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
    
    # listeV = tabla1_2[0]
    # listq = tabla1_2[1]
    # listx = tabla1_2[2]
    listz = tabla1_2[3]     ##
    listEELS = tabla1_2[4]
    
    listz_norm_a_BEM = np.array(listz)/a   ## we interpolate using z/W, same axis as the potential code
    EELS_interp =  interp1d(listz_norm_a_BEM, listEELS)
 
    return listz_norm_a_BEM, EELS_interp
    
def P_integrand_over_z(value_z_norm_a, theta, V0, Ee_electron, bb, ss, dd, energy_eV, a, N, list_z_norm_a, listV_normV0 ):
    """
    Parameters
    ----------
    value_z_norm_a : z/a coordinate normalized to a
    theta : angle of the electron (initial conditions) in radians
    V0 : potential of the gate in e*Volts
    Ee_electron : energy of the electron in eV, this determinates beta
    bb : aspect ratio between width and height
    ss : rounding radius normalized to width 
    dd : distance to the plane normalized to width
    energy_eV: photon energy hb*omega in eV
    a : width along x (nm)
    N : number of points in bem2d
    list_z_norm_a : list of z/width from c++ code 
    listV_normV0 : list of V(z)/V0 from c++ code
    Returns
    -------
    P(omega) integrand over z/a 
    see "Eq_P_integral_over_z.png"
    """
    epsi1 = 1
    me_c2_eV = 510998.95069  ## me*c**2 in eV
    beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    aux_function = beta*np.sin(theta)
    # aux_function = beta*theta
    
    # list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N)
    # list_z_norm_a_interp = np.linspace(1.1,2,int(N*Nint)) 
    V_interp =  interp1d(list_z_norm_a, listV_normV0)
    
    V_norm_V0 = V_interp(value_z_norm_a)  ## V(z/a)/V0
    denominator = np.sqrt(aux_function**2 + 2*V_norm_V0*V0/(me_c2_eV*gamma_e))

    h = bb*a
    listz_norm_a_BEM, EELS_interp = EELS_from_BEM_interpolated(energy_eV,a,h,N)

    EELS_value = EELS_interp(value_z_norm_a) 
    
    function = 2*beta*EELS_value/denominator
    
    # print(V_norm_V0,EELS_value)
    
    return function*a ## we add "a" because the integral is going to be over z/a


## integration over z as a sum 
def P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy_eV, a, N , list_z_norm_a, listV_normV0):
    """
    Parameters
    ----------
    z_min_val : minimum value of the integration over z in z/W
    theta : angle of the electron (initial conditions) in radians
    V0 : potential of the gate in e*Volts
    Ee_electron : energy of the electron in eV, this determinates beta
    bb : aspect ratio between width and height
    ss : rounding radius normalized to width 
    dd : distance to the plane normalized to width
    energy_eV: photon energy hb*omega in eV
    a : width along x (nm)
    N : number of points in bem2d
    list_z_norm_a : list of z/width from c++ code 
    listV_normV0 : list of V(z)/V0 from c++ code
    Returns
    -------
    P(omega) integrated over z/a 
    see "Eq_P_integral_over_z.png"
    """
    Nint = 10
    h = bb*a
    listz_norm_a_BEM, EELS_interp = EELS_from_BEM_interpolated(energy_eV,a,h,N)
    listz_norm_a_BEM_interp = np.linspace(z_min_val, np.max(listz_norm_a_BEM),int(N*Nint)) 
    
    list_P = 0
    for value_z_norm_a in listz_norm_a_BEM_interp: ## the integration should be from z_min 
        P_value = P_integrand_over_z(value_z_norm_a, theta, V0, Ee_electron, bb, ss, dd, energy_eV, a, N, list_z_norm_a, listV_normV0)
        list_P = list_P + P_value
    delta_z_norm_a = listz_norm_a_BEM_interp[1] - listz_norm_a_BEM_interp[0]
    
    return list_P*delta_z_norm_a

 
def find_width_of_peak(listx,listy):
    """
    Parameters
    ----------
    listx : x array
    listy : y array
    Returns
    -------
    width of the peak as a array 
    x_left, x_right
    """
    
    print('5-Load the EELS integrated over z/W as a function of energy and find the width of the highest peak')
    
    Nint = 10
    ######## interpolation of the data ######################################
    x_int = np.linspace(np.min(listx),np.max(listx),int(len(listx)*Nint))
    cs = interp1d(listx, listy)
    y_int = cs(x_int)
    
    y_int2 = [] ## filter small negative numbers created with the interpolation
    for yy in y_int: 
        if yy<0: 
            y_int2.append(0)
        else:
            y_int2.append(yy)
        
    height = np.max(y_int2)/3 ## real modes --> interpolation creates fake peaks with small intensity, a way
    # to get rid of them is to specify a big height, then only the real peaks (without interpolation) will remain 
    ######## find peaks ####################################################
    peaks, _ = find_peaks(listy, height=height )
    
    ######## find width : we move left and right from the peak until the function starts increasing again###########################
    x_left_peak_index = []
    x_left_peak_value = []
    for peak in peaks:
        i = peak
        # Move left as long as signal is decreasing
        while i > 0 and listy[i-1] <= listy[i] and listy[i]>0:
            i -= 1
        x_left_peak_index.append(i)
        x_left_peak_value.append(listx[i])
    
    
    x_right_peak_index = []
    x_right_peak_value = []

    for peak in peaks:
        i = peak
        while i < len(listy) - 1 and listy[i+1] <= listy[i] and listy[i]>0:
            i += 1
        x_right_peak_index.append(i)
        x_right_peak_value.append(listx[i])
        
    
    # results_full = peak_widths(y_int2, peaks, rel_height=3)
    # results_full_xmin_index = []
    # for xmin in results_full[2]:
    #     results_full_xmin_index.append(int(xmin))
    # results_full_xmax_index = []
    # for xmax in results_full[3]:
    #     results_full_xmax_index.append(int(xmax))

    # results_full_xmin_index = np.array(results_full_xmin_index)
    # results_full_xmax_index = np.array(results_full_xmax_index)
    ## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
    sorted_listy = np.sort(listy[peaks])  
    sorted_index = np.argsort(listy[peaks]) 
    sorted_listx = listx[peaks[sorted_index]]
    # sorted_xmin = x_int[results_full_xmin_index[sorted_index]]
    # sorted_xmax = x_int[results_full_xmax_index[sorted_index]]
    
    # witdh = [] 
    # for j in range(len(sorted_xmin)):
    #     delta_x = sorted_xmax[j] - sorted_xmin[j]
    #     witdh.append(delta_x)
    
    labelx='Photon energy $\hbar\omega$ (eV)'
    labely='EELS per electron (1/eV)'
 
    tamfig = [4.5,3.5]
    tamletra = 13
    tamtitle  = tamletra - 3
    tamnum = tamletra
    labelpady = 3
    labelpadx = 2
    pad = 2.5
    dpi = 500
    
    print('5b-Find the wmin and wmax of the mode = mode')
    plot_figure = 1
    if plot_figure == 1:
        plt.figure(figsize=tamfig)
        plt.title('',fontsize=tamtitle)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
        plt.plot(x_int, y_int2 ,'--',lw = 1.5 )
        plt.plot(listx, listy ,'-' )
        plt.plot(sorted_listx, sorted_listy, "x")
        plt.plot(x_left_peak_value, np.ones(len(x_left_peak_value))*0, "x",color = 'blue')
        plt.plot(x_right_peak_value, np.ones(len(x_right_peak_value))*0, "x",color = 'red')
        plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
    
    return x_left_peak_value, x_right_peak_value
