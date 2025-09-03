#code e.h : 
#- Models 3 biased rectangular wires (center of width W and 2 on the side of width W', all of them the same height h).
#- Uses the boundary element method to compute induced surface charge.
#- Computes electrostatic potential/V0 at any point in space (x,y)
# where: 
#wp = 0.5 w'/w
#d = 0.2, d/w distance between the side wires normalized to w
#h = 0.6, aspect ratio height/w
#N = 200, discretization points.
#epsilon = permittivity = 9
#N = 200, discretization points.
#x0 = 0, x0/w
#x1 = 2 , x1/w
#nx = 200
#z0 = h  z0/w ## start outside the waveguide
#z1 = 2  z1/w 
#nz = 200
## x=0 is in the middle of the center waveguide
## z=0 is at the interface

### All lengths are given in units of the width of the center waveguide, label as "w"

## plot the z min as a function of theta, V0. solutions only above the wg 
## return the V0,theta for a fix value of z 

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,interp1d
from scipy.optimize import curve_fit
import os
 
#%%

#print('Run the c++ code to get the points of the potential')

def run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz): ## fixed x 
    x1 = x0
    nx = 1
    # Run the C++ program name dy.out and save a list of x, z, V
    cmd = ["./e.out", str(wp), str(d), str(h), str(N), str(epsilon), str(x0), str(x1), str(nx), str(z0), str(z1), str(nz)]
    
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

# Define linear model
def linear_func(z, a, b):
    return a*z + b

#%%

# ========= main function =========

def fit_potential_linear(w, wp_norm_w, d_norm_w, h_norm_w, N, epsilon, x0_norm_w,
                         b_norm_w, delta, plot_figure = True):
    """
    Parameters
    ----------
    w : width of the middle waveguide in nm
    wp_norm_w : wp/w, wp width of the side waveguides
    d_norm_w : d/w, d distance between the side waveguides and the central one
    h_norm_w : h/w, h height of the waveguides
    N : discretization points
    epsilon : epsilon of the substrate, should be = 9 (DC of sapphire)
    x0_norm_w : xi/w = xf/w and nx = 1
    b_norm_w  : defined such that the fit is around (zdelta - delta, zdelta + delta) with
    zdelta = b_norm_w + h/w
    delta: length of the z values for the fitting is 2*delta

    Returns
    -------
    Runs the C++ solver, normalizes potential, performs a linear fit near zdelta,
    saves slope and offset to a text file named with hdelta.
    """
    tamfig = [4,3]
    tamletra = 14
    tamtitle  = 8
    tamnum = tamletra
    tamlegend = tamletra
    labelpady = 2
    labelpadx = 3
    pad = 2.5
    dpi = 300
    
    # ---- simulation parameters ----
    z0 = h_norm_w*1.001    
    z1 = 1.2    
    nz = 500
    title = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f' % (wp_norm_w,h_norm_w,d_norm_w,x0_norm_w)
    labelx = r'$z/W$'
    labely = r'$\phi(x,z)/V_0$'
    
    # ---- run solver ----
    list_z_norm_a, listV_normU = run_e_out1(wp_norm_w, d_norm_w, h_norm_w, N, epsilon, x0_norm_w, z0, z1, nz)
    zz, V0_over_U_list = run_e_out1(wp_norm_w, d_norm_w, h_norm_w, N, epsilon, 0, h_norm_w, h_norm_w, 2)
    V0_over_U = np.abs(V0_over_U_list[0]) 
    listV_normV0 = listV_normU / V0_over_U

    # ---- interpolation ----
    V_interp = interp1d(list_z_norm_a, listV_normV0, kind='linear')

    # ---- fit window ----
    zdelta = h_norm_w + b_norm_w

    if zdelta - delta < z0:
        listx_for_fit = np.linspace(z0, zdelta + delta, int(N))
    else:
        listx_for_fit = np.linspace(zdelta - delta, zdelta + delta, int(N))

    listy_for_fit = V_interp(listx_for_fit)
    # ---- linear fitting ----
    popt, _ = curve_fit(linear_func, listx_for_fit, listy_for_fit)
    a, b = popt
    after_fitting = linear_func(listx_for_fit, a, b)

    # ---- saving data ----
    path_basic = os.getcwd()
    path_data = os.path.join(path_basic, 'potential_data')
    os.makedirs(path_data, exist_ok=True)
    os.chdir(path_data)

    # ---- save results ----
    results = np.array([[a, b]])
    header_info = f"""# Linear fit of V(z)/V0
    # Parameters:
    # wp = {wp_norm_w}, h/W = {h_norm_w}, d/W = {d_norm_w}, epsilon = {epsilon}
    # x0 = {x0_norm_w}, w = {w} nm, b_norm_w = {b_norm_w}    
    # Fit window delta = {delta}
    # Columns: slope(a), offset(b) such as V(z)/V0 = a*(z/w) + b
    """
    label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,x0_norm_w) ## all normalize to width w
    filename = f"linear_fit{label_txt}_bmin%.2f.txt" %(b_norm_w)
    np.savetxt(filename, results, header=header_info)
    print(f"Saved fit results to {filename}")
    
    # ---- plotting ----
    if plot_figure == True:
        plt.figure(figsize=tamfig)
        plt.title(title,fontsize=tamtitle)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
        plt.plot(list_z_norm_a, listV_normV0,'.-',label = r'$\phi$ from BEM' )
        plt.plot(listx_for_fit, after_fitting,'--',label = r'Linear fit' )
        plt.plot([zdelta],[linear_func(zdelta,a,b)],'o',color = 'black',label = r'$z_{\text{min}}$')
        #plt.xticks(np.arange(0.5,2.25,0.25),["0.5","","1.0","","1.5","","2.0"])
        plt.xticks(np.arange(0.5,1.2+0.1,0.1),['', '0.6', '', '0.8', '', '1.0' , '', '1.2'])
        plt.yticks(np.arange(-1,0,0.1),[r"$-1.0$" , "", r"$-0.8$", "", r"$-0.6$", "", r"$-0.4$", "", r"$-0.2$", ""])
        plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
        plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
        plt.tight_layout()
        plt.savefig(f'potential_fitted_b_norm_w{b_norm_w}.png',
                    bbox_inches='tight',pad_inches = 0.01, dpi=dpi)
        
        plt.savefig(f'potential_fitted_b_norm_w{b_norm_w}.pdf',
                    bbox_inches='tight',pad_inches = 0.01, dpi=dpi)
        
        plt.show()
    
    return a, b



if __name__ == "__main__":
    import sys
    
    # Read command line arguments
    if len(sys.argv) != 10:
        print("Usage: python fit_potential_as_linear.py w wp_norm_w d_norm_w h_norm_w N epsilon x0_norm_w b_norm_w delta")
        sys.exit(1)

    w = float(sys.argv[1])
    wp_norm_w = float(sys.argv[2])
    d_norm_w = float(sys.argv[3])
    h_norm_w = float(sys.argv[4])
    N = int(sys.argv[5])
    epsilon = float(sys.argv[6])
    x0_norm_w = float(sys.argv[7])
    b_norm_w = float(sys.argv[8])
    delta = float(sys.argv[9])

    args = [float(arg) for arg in sys.argv[1:]]

    labels = ["w (nm)", "wp/w", "d/w", "h/w", "N", "epsilon", "x0/w", "bmin/w", "delta"]
    print("Input Parameters:")
    for label, value in zip(labels, args):
        print(f"{label}: {value}")

    # Call your function
    fit_potential_linear(*args)