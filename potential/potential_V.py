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
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq,root

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

list_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
 '#7f7f7f', '#bcbd22', '#17becf']

#%%

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
    

bb = 2
ss = 0.1
dd = 10
N = 200
list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)

xmin = -bb/2
xmax = bb/2
aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)

labelx = r'$z/a$'
labely = r'$V/V_0$'

plt.figure(figsize=tamfig)
plt.plot(list_z_norm_a, listV_normV0 ,'.-' )
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('Vy.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

Nint = 1
## the electron is between the gate and the rectangular waveguide. NO
## the electron after the WG 
list_z_norm_a_interp = np.linspace(np.min(list_z_norm_a), np.max(list_z_norm_a),int(N*Nint)) 
V_interp =  CubicSpline(list_z_norm_a, listV_normV0)

#%%

me_c2_eV = 510998.95069  ## me*c**2 in eV
ind = 2
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
V0 = Ee_electron/2
## minimum z
def function_to_be_zero(value_z_norm_a, theta, V0):
    """
    Parameters
    ----------
    value_z_norm_a : z/a
    theta : angle of the electron (initial conditions) in radians
    V0 : potential of the gate in e*Volts
    Returns
    -------
    V(z)/V0 for z min --> see notes of 2025-05-27 units_potential 
    """
    epsi1 = 1
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    aux_function = beta*np.tan(theta)
    # aux_function = beta*theta
    function = aux_function**2*me_c2_eV*gamma_e/2
    
    V_norm_V0 = V_interp(value_z_norm_a)
    
    partA = function/V0
    partB = V_norm_V0
    # print(partA/partB)
    return partA - partB ## add a minus to the potential of the WG --> has to be negative 



list_theta_mrad = [0.1,0.25,0.5,1,5,10]
plt.figure(figsize=tamfig)
k = 0
for theta_mrad in list_theta_mrad:
    theta = theta_mrad*1e3
    list_function_to_be_zero = [] 
    for value_z_norm_a in list_z_norm_a_interp:
        list_function_to_be_zero.append(function_to_be_zero(value_z_norm_a,theta,V0))
    
    sol = root(function_to_be_zero, [bb/2], method='hybr', args=(theta,V0))
    plt.plot(list_z_norm_a_interp, list_function_to_be_zero ,'.-' ,color = list_colors[k],label = r'$\theta$ = %.2f mrad' %(theta_mrad))
    print(sol.x)
    plt.plot([sol.x],[0],'x',color = list_colors[k])
    k = k + 1
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\Delta$',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('minimum_z_vs_angle' + label_Ee + '.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%

Nvals = 50
# Define ranges for V0 and theta
V0_vals = np.linspace(Ee_electron/5, Ee_electron, Nvals)
theta_mrad_vals = np.linspace(0.1, 10, Nvals)

zmin = bb/2 ## up from the WG
zmax = np.max(list_z_norm_a)
# Storage for roots
X_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
f_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
# Loop over parameter values
for i, V0 in enumerate(V0_vals):
    for j, theta_mrad in enumerate(theta_mrad_vals):
        theta = theta_mrad*1e3
        try:
            # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
            x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta,V0))
            X_vals[i, j] = x_root
            value_function = function_to_be_zero(x_root, theta, V0)
            f_vals[i, j] = value_function
        except ValueError:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan




plt.figure(figsize=tamfig)
V0_grid, theta_grid = np.meshgrid(np.array(V0_vals)*1e-3, theta_mrad_vals, indexing='ij')
plt.contourf(V0_grid, theta_grid, X_vals, levels=50, cmap="viridis")
plt.colorbar(label=r'$z_{min}/a$')
plt.xlabel(r'$V_0$ (keV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('minimum_z' + label_Ee + '.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

