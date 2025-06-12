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

## plot the z min as a function of theta, V0. solutions only above the wg 
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq,root
import matplotlib as mpl

tamfig = [4, 3]
tamletra = 13
tamtitle  = tamletra - 5
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

#%%

print('IMPORTANT: change the lines to numero ymax=bb/2+15*dd numero ymin=-bb/2-dd inside the dy.cpp and compile it')
print('1-Run the c++ code to plot the potential of a rectangular waveguide')

ss = 0.1
dd = 10
N = 200
me_c2_eV = 510998.95069  ## me*c**2 in eV
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)

Nvals = 200
V0_vals = np.linspace(0.1, 10, Nvals)
theta_mrad = 5
theta = theta_mrad*1e-3
 
print('2-Define the function to be zero--> zmin (electron position) for a fixed V0 and different angles')
print('3-zmin as a function of V0 for theta=5mrad, for zmin (electron position) above the WG')
## minimum z
labelx = r'$V_0$ (eV)'
labely = r'$b_{\text{min}}/W$'

title1 = r'$\theta$ = %i mrad, $E_{\text{e}}$ = %i keV, $s/W$ =%.1f, $d/W$ = %i'%(theta_mrad,Ee_electron_keV,ss,dd)
listbb = [0.5,0.75,1]

plt.figure(figsize=tamfig)
plt.title(title1 ,fontsize=tamtitle)    
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$b_{\text{min}}/W$',fontsize=tamletra,labelpad =labelpady)

for bb in listbb: 
    def run_dy_out(bb,ss,dd, N):
        # Run the C++ program name dy.out and save a list of y and V
        cmd = ["./dy.out", str(bb), str(ss), str(dd), str(N)]
     
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
    
    list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)
    
    def function_to_be_zero(value_z_norm_a, theta, V0,bb,ss,dd):
        """
        Parameters
        ----------
        value_z_norm_a : z/a
        theta : angle of the electron (initial conditions) in radians
        V0 : potential of the gate in e*Volts
        bb : aspect ratio 
        ss : rounding radius normalized to x-length
        dd : distance to the plane normalized to x-length
        Returns
        -------
        V(z)/V0 for z min --> see notes of 2025-05-27 units_potential 
        """
        epsi1 = 1
        gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
        aux_function = beta*np.sin(theta)
        # aux_function = beta*theta
        function = aux_function**2*me_c2_eV*gamma_e/2
    
        # print(partA/partB)
    
        V_interp =  CubicSpline(list_z_norm_a, listV_normV0)
        V_norm_V0 = V_interp(value_z_norm_a)
        
        partA = function/V0
        partB = V_norm_V0
        
        return partA - partB ## add a minus to the potential of the WG --> has to be negative 
    
    zmin = bb/2 ## above the WG
    zmax = np.max(list_z_norm_a)
    # Storage for roots
    X_vals_2 = np.zeros(len(V0_vals))
    f_vals_2 = np.zeros(len(V0_vals))
    # Loop over parameter values
    for i, V0 in enumerate(V0_vals):
    
        try:
            # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
            x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta,V0,bb,ss,dd))
            X_vals_2[i] = x_root-bb/2 ## we subtract the height of the waveguide
            value_function = function_to_be_zero(x_root, theta, V0,bb,ss,dd)
            f_vals_2[i] = value_function
        except ValueError:
            # No root found in the interval
            X_vals_2[i] = np.nan
            f_vals_2[i] = np.nan
                    

    plt.plot( V0_vals, X_vals_2, label = r'$h/W$ = %.2f' %(bb))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 

plt.xticks(np.arange(0,12,2))
plt.savefig('minimum_z_vs_V0_dd%i' %(dd) + label_Ee + '.png' , format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()
   