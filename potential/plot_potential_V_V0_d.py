#code d.h : 
#- Models a grounded plane + biased rectangular wire with rounded corners.
#- Uses the boundary element method to compute induced surface charge.
#- Computes electrostatic potential/V0 at any point in space (x,y)
# where: 
#bb = 2, aspect ratio  (side length along y (labelled as b) divided by side length along x).
#ss = 0.1, rounding radius normalized to height.
#dd = 1, distance to the plane normalized to height.
#N = 200, discretization points.
### All lengths are given in units of the x side-length, label as "a"

## plot the z min as a function of d, V0
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq,root
import matplotlib as mpl

tamfig = [4, 3]
tamletra = 13
tamtitle  = tamletra - 4
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

print('IMPORTANT: change the lines to numero ymax=6; numero ymin=-6 inside the dy.cpp and compile it, if theta = 2mrad')
print('1-Run the c++ code to plot the potential of a rectangular waveguide')

bb = 2
ss = 0.1
N = 100

me_c2_eV = 510998.95069  ## me*c**2 in eV
ind = 2
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
V0 = Ee_electron/1e4 
V0 = 1

theta_mrad = 2
theta = theta_mrad*1e-3

labelx = r'$z/a$'
title1 = r'$h/a$ = %i, $s/a$ =%.1f, $\theta$ = %i mrad, $E_{\text{e}}$ = %i keV' % (bb,ss,theta_mrad,Ee_electron_keV)
if theta_mrad == 2:
    list_dd = [5,10,50,100]
else:
    list_dd = [25,100,150,200]
# list_dd = [5 ]
from mycolorpy import colorlist as mcp
color1 = mcp.gen_color(cmap="hot",n=len(list_dd)+2)
 
def run_dy_out( bb,ss,dd, N  ):
    # Run the C++ program name dy.out and save a list of y and V
    cmd = ["./dy.out", str(bb), str(ss), str(dd), str(N)]
 
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = result.stdout.strip().split('\n')
    list_z_norm_a = []
    listV_normV0 = []
    for line in lines:
        try:
            zvalue, Vvalue = map(float, line.strip().split())
            list_z_norm_a.append(zvalue)
            listV_normV0.append(Vvalue)  
        except ValueError:
            continue

                              
    return list_z_norm_a,listV_normV0
 

def function_to_be_zero(value_z_norm_a,V_norm_V0,V0):
    """
    Parameters
    ----------
    value_z_norm_a : z/a
    V_norm_V0 :  V(z/a)/V0
    V0 : potential of the gate in e*Volts
    Returns
    -------
    V(z)/V0 for z min --> see notes of 2025-05-27 units_potential 
    """
    epsi1 = 1
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    aux_function = beta*np.sin(theta)
    # aux_function = beta*theta
    function = aux_function**2*me_c2_eV*gamma_e/2
    
 
    partA = function / V0
    partB = V_norm_V0
        
    return partA - partB ## add a minus to the potential of the WG --> has to be negative 
    

zmin = -6
zmax = 6
plt.figure(figsize=tamfig)
plt.title(title1 + ', $V_0$ = %.1f eV' %(V0),fontsize=tamtitle)
k = 0
for dd in list_dd:
    print(dd)
    list_z_norm_a, listV_normV0 = run_dy_out( bb,ss,dd, N  )
    V_norm_V0 = CubicSpline(list_z_norm_a, listV_normV0)
    list_function_to_be_zero = [] 
    for value_z_norm_a in list_z_norm_a:
        V_norm_V0_value = V_norm_V0(value_z_norm_a)
        list_function_to_be_zero.append(function_to_be_zero(value_z_norm_a,V_norm_V0_value,V0))
    
#    sol1 = root(lambda z: function_to_be_zero(z, V_norm_V0_value, V0), [-bb], method='hybr')  ## electron between PEC and WG
#    sol2 = root(lambda z: function_to_be_zero(z, V_norm_V0_value, V0), [bb], method='hybr')   ## electron above the WG
    
    sol1 = brentq(lambda z: function_to_be_zero(z, V_norm_V0(z), V0), zmin, -bb/2)  ## electron between PEC and WG
    sol2 = brentq(lambda z: function_to_be_zero(z, V_norm_V0(z), V0), bb/2, zmax)   ## electron above the WG
    
    
    plt.plot(list_z_norm_a, list_function_to_be_zero ,'.-' ,color = color1[k],label = r'$d/a$ = %i' %(dd))
    print(sol1, sol2)
    

    plt.plot([sol1],[0],'x',color = 'black')

    plt.plot([sol2],[0],'x',color = 'black')
    k = k + 1
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)

plt.ylabel(r'$\Delta$',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('minimum_z_vs_dd' + label_Ee + '_V0_%ieV.png' %(V0), format='png',bbox_inches='tight',pad_inches = 0.01, dpi=dpi)  
plt.show()

#%%
print('3-Color map of zmin as a function of V0 and dd, for zmin (electron position) above the WG')

Nvals = 100
# Define ranges for V0 and theta
V0_vals = np.linspace(0.5, 1.5, Nvals)
if theta_mrad == 2:
    dd_vals = np.linspace(5, 100, Nvals)
else:
    dd_vals = np.linspace(25, 200, Nvals)
limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(dd_vals) , np.max(dd_vals)]

# Storage for roots
X_vals = np.zeros((len(V0_vals), len(dd_vals)))
f_vals = np.zeros((len(V0_vals), len(dd_vals)))

list_z_norm_a_tot = []
listV_normV0_tot = []
for dd in dd_vals:
    list_z_norm_a, listV_normV0  = run_dy_out(bb, ss, dd, N)
    list_z_norm_a_tot.append(list_z_norm_a)
    listV_normV0_tot.append(listV_normV0)

# Loop over parameter values
for i, V0 in enumerate(V0_vals):
    for j, dd in enumerate(dd_vals):
        V_norm_V0 = CubicSpline(list_z_norm_a_tot[j], listV_normV0_tot[j])
        try:
            # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
            x_root =  brentq(lambda z: function_to_be_zero(z, V_norm_V0(z), V0), bb/2, zmax) 
            X_vals[i, j] = x_root-bb/2 ## we subtract the height of the waveguide
            value_function = function_to_be_zero(x_root, V_norm_V0(x_root), V0)
            f_vals[i, j] = value_function
        except ValueError:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan
 
n_color = 10
vmin1 , vmax1 = np.nanmin(X_vals), np.nanmax(X_vals)
cmap = plt.cm.viridis  # define the colormap
bounds1 =   np.linspace(vmin1, vmax1 , n_color) 

bounds1 =  [vmin1,0.25,0.5,0.75,2,5]
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
im_show = plt.imshow(np.transpose(X_vals), extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'  ,norm = norm1 ) 
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%.2f') 
cbar.ax.set_title(r'$b^{(2)}_{min}/a$',fontsize=tamletra)
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$d/a$',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin2' + label_Ee +  '_theta_%imrad.png' %(theta_mrad), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

#%%

print('4-Color map of zmin as a function of V0 and dd, for zmin (electron position) between the WG and the PEC')
 
 
# Storage for roots
X_vals2 = np.zeros((len(V0_vals), len(dd_vals)))
f_vals2 = np.zeros((len(V0_vals), len(dd_vals)))
# Loop over parameter values
for i, V0 in enumerate(V0_vals):
    for j, dd in enumerate(dd_vals):
        V_norm_V0 = CubicSpline(list_z_norm_a_tot[j], listV_normV0_tot[j])
        try:
            # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
            x_root =  brentq(lambda z: function_to_be_zero(z, V_norm_V0(z), V0), zmin, -bb/2) 
            X_vals2[i, j] = x_root-bb/2 ## we subtract the height of the waveguide
            value_function = function_to_be_zero(x_root, V_norm_V0(x_root), V0)
            f_vals2[i, j] = value_function
        except ValueError:
            # No root found in the interval
            X_vals2[i, j] = np.nan
            f_vals2[i, j] = np.nan


vmin2 , vmax2 = np.nanmin(X_vals2), np.nanmax(X_vals2)
bounds2 =   np.linspace(vmin2, vmax2 , n_color) 
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap.N)
plt.figure(figsize=tamfig)
im_show = plt.imshow(np.transpose(X_vals2), extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'    ) 
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%.2f') 
cbar.ax.set_title(r'$b^{(1)}_{min}/a$',fontsize=tamletra)
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$d/a$',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin1' + label_Ee + '_theta_%imrad.png' %(theta_mrad), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()


