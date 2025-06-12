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
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq,root
import matplotlib as mpl

tamfig = [4, 3]
tamletra = 13
tamtitle  = tamletra - 3
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
    

bb = 0.25
ss = 0.1
dd = 500
N = 200
list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)

xmin = -bb/2
xmax = bb/2
aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)

labelx = r'$z/a$'
labely = r'$V/V_0$'

title1 = r'$h/a$ = %.2f, $s/a$ =%.1f, $d/a$ = %i' % (bb,ss,dd)
plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
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

print('2-Define the function to be zero--> zmin (electron position) for a fixed V0 and different angles')

me_c2_eV = 510998.95069  ## me*c**2 in eV
ind = 2
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
V0 = Ee_electron/1e4 
V0 = 5
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
    aux_function = beta*np.sin(theta)
    # aux_function = beta*theta
    function = aux_function**2*me_c2_eV*gamma_e/2
    
    V_norm_V0 = V_interp(value_z_norm_a)
    
    partA = function/V0
    partB = V_norm_V0
    # print(partA/partB)
    return partA - partB ## add a minus to the potential of the WG --> has to be negative 


aux_ejey = np.linspace(-np.max(listV_normV0),np.min(listV_normV0),10)

list_theta_mrad = [0.1,0.25,0.5,1,2]
from mycolorpy import colorlist as mcp
color1 = mcp.gen_color(cmap="hot",n=len(list_theta_mrad)+2)

plt.figure(figsize=tamfig)
plt.title(title1 + r', $V_0$ = %.1f eV, $E_{\text{e}}$ = %i keV' %(V0,Ee_electron_keV),fontsize=tamtitle)
k = 0
for theta_mrad in list_theta_mrad:
    theta = theta_mrad*1e-3
    list_function_to_be_zero = [] 
    for value_z_norm_a in list_z_norm_a_interp:
        list_function_to_be_zero.append(function_to_be_zero(value_z_norm_a,theta,V0))
    
    sol1 = root(function_to_be_zero, [-bb/2], method='hybr', args=(theta,V0))  ## electron between PEC and WG
    sol2 = root(function_to_be_zero, [bb/2], method='hybr', args=(theta,V0))   ## electron above the WG
    
    
    plt.plot(list_z_norm_a_interp, list_function_to_be_zero ,'.-' ,color = color1[k],label = r'$\theta$ = %.2f mrad' %(theta_mrad))
    print(sol1.x, sol2.x)
    
    if sol1.fun <= 1e-15:
        plt.plot([sol1.x],[0],'x',color = 'black')
    if sol2.fun <= 1e-15:
        plt.plot([sol2.x],[0],'x',color = 'black')
    k = k + 1
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\Delta$',fontsize=tamletra,labelpad =labelpady)
plt.plot(np.ones(10)*xmin,aux_ejey,'--',lw = 0.8,color = 'black')
plt.plot(np.ones(10)*xmax,aux_ejey,'--',lw = 0.8,color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('minimum_z_vs_angle_dd%i_hh%.2f' %(dd,bb) + label_Ee + '_V0_%ieV.png' %(V0), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%
print('3-Color map of zmin as a function of V0 and theta, for zmin (electron position) above the WG')

Nvals = 200
# Define ranges for V0 and theta
V0_vals = np.linspace(0.005, 5, Nvals)
theta_mrad_vals = np.linspace(0.1, 5, Nvals)
limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]

zmin = bb/2 ## above the WG
zmax = np.max(list_z_norm_a)
# Storage for roots
X_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
f_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
# Loop over parameter values
for i, V0 in enumerate(V0_vals):
    for j, theta_mrad in enumerate(theta_mrad_vals):
        theta = theta_mrad*1e-3
        try:
            # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
            x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta,V0))
            X_vals[i, j] = x_root-bb/2 ## we subtract the height of the waveguide
            value_function = function_to_be_zero(x_root, theta, V0)
            f_vals[i, j] = value_function
        except ValueError:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan



#%%
n_color = 10
vmin1 , vmax1 = np.nanmin(X_vals), np.nanmax(X_vals)
cmap = plt.cm.viridis  # define the colormap
bounds1 =   np.linspace(vmin1, vmax1 , n_color) 
if dd == 20:
    bounds1 =  [vmin1,0.25,0.5,5,100,200,int(vmax1)]
    # bounds1 =  [vmin1,0.25,0.5,5,int(vmax1)]
else:
    bounds1 =  [vmin1,0.25,0.5,5,50,int(vmax1)]
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)

plt.figure(figsize=tamfig)
im_show = plt.imshow(np.transpose(X_vals), extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'  ,norm = norm1  ) 
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%.2f') 
cbar.ax.set_title(r'$b^{(2)}_{min}/a$',fontsize=tamletra)
plt.xticks(np.arange(0,6,1))
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin2' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

#%%

print('4-Color map of zmin as a function of V0 and theta, for zmin (electron position) between the WG and the PEC')

theta_mrad_vals = np.linspace(0.1, 1, Nvals)
limits2 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
zmin = -bb/2-dd ## above the WG
zmax = np.max(list_z_norm_a)
# Storage for roots
X_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
f_vals = np.zeros((len(V0_vals), len(theta_mrad_vals)))
# Loop over parameter values
for i, V0 in enumerate(V0_vals):
    for j, theta_mrad in enumerate(theta_mrad_vals):
        theta = theta_mrad*1e-3
        try:
            # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
            x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta,V0))
            X_vals[i, j] = x_root-bb/2 ## we subtract the height of the waveguide
            value_function = function_to_be_zero(x_root, theta, V0)
            f_vals[i, j] = value_function
        except ValueError:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan

vmin2 , vmax2 = np.nanmin(X_vals), np.nanmax(X_vals)
bounds2 =   np.linspace(vmin2, vmax2 , n_color) 
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap.N)
plt.figure(figsize=tamfig)
im_show = plt.imshow(np.transpose(X_vals), extent = limits2, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm = norm2  ) 
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%.2f') 
cbar.ax.set_title(r'$b^{(1)}_{min}/a$',fontsize=tamletra)
plt.xticks(np.arange(0,6,1))
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin1' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()


