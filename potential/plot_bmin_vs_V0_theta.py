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
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import matplotlib as mpl
import os
 
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

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'zmin_solutions')
load_data = 0

#%%

os.chdir(path_basic)

print('Run the c++ code to get the points of the potential')

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

#%%

print('Define paramaters')

me_c2_eV = 510998.95069  ## me*c**2 in eV
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)    

wp = 0.5 ## w'/w
d = 0.2 ## d/w distance between the side wires normalized to w
h = 0.4   ## aspect ratio height/w
N = 400   ## discretization points.
epsilon = 9 ## permittivity = 9
N = 200   ## discretization points.
x0 = 0     ## x0/w --> position of the electron. options [0,0.2,0.4,0.6]. 
            ## x0=/0.6/0.8 does not give solutions because potential stops being negative: check the plots of V vs x0
z0 = h*1.001    ## z0/w ## start outside the waveguide
z1 = 2    ## z1/w 
nz = 200
# x=0 is in the middle of the center waveguide
# z=0 is at the interface
 
w = 500 ## 100 400 600   

Nvals = 500
bmin_vals0 = 50/w  ## 50 nm
# bmin_vals0 = 60/a
# bmin_vals0 = 80/a
bmin_val = bmin_vals0*w
Nvals = 300
if Ee_electron_keV == 200:
    theta_mrad_vals = np.linspace(0.01, 2.25, Nvals)
    V0_vals = np.linspace(0.01, 1, Nvals)
elif Ee_electron_keV == 100:
    theta_mrad_vals = np.linspace(0.01, 3, Nvals)
    V0_vals = np.linspace(0.01, 1, Nvals)
    
V0_0 = V0_vals[0]
theta_mrad0 = theta_mrad_vals[0]

# bmin_vals0 = 0.1
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp,h,d,x0) ## all normalize to width w

tol = 1e-6 ## zero of the function 
tol = 5*1e-1 ## zero of the function 

density = 50 ## add more points to the interface of (nan-not nan) values 

#%%

labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)
title = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f' % (wp,h,d,x0)

#%%

list_z_norm_a, listV_normV0 = run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz)
Nint = 1
## the electron is between the gate and the rectangular waveguide. NO
## the electron after the WG 
list_z_norm_a_interp = np.linspace(np.min(list_z_norm_a), np.max(list_z_norm_a),int(N*Nint)) 
V_interp =  CubicSpline(list_z_norm_a, listV_normV0)

#%%

print('2-Define the function to be zero --> zmin (electron position) for a fixed V0 and different angles')
                                
zz, V0_over_U_list = run_e_out1(wp,d,h,N,epsilon,0,h,h,2)
V0_over_U = np.abs(V0_over_U_list[0]) ## normalized to V0
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
    beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    aux_function = beta*np.sin(theta)
    # aux_function = beta*theta
    function = aux_function**2*me_c2_eV*gamma_e/2
    
    V_norm_U = V_interp(value_z_norm_a) 
    # print(V0_over_U)
    
    V_norm_V0 = V_norm_U/V0_over_U
    
    partA = function/V0
    partB = V_norm_V0
    # print(partA/partB)
    return partA + partB 

#%%

print('3-Color map of zmin as a function of V0 and theta, for zmin (electron position) above the WG')

# theta_mrad_vals = np.linspace(0.1, 1, Nvals)
limits1 = [np.min(V0_vals) , np.max(V0_vals), np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]

zmin = h ## above the WG
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
            X_vals[i, j] = x_root-h ## we subtract the height of the waveguide
            value_function = function_to_be_zero(x_root, theta, V0)
            f_vals[i, j] = value_function
        except ValueError or x_root == h:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan

#%%
X_vals_transpose = np.transpose(X_vals)

tamfig2 = [5, 3]
limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
cmap = plt.cm.RdBu   # define the colormap

# another cbar with constant a*1e-3 (real units, in microns)
cte_cbar2 = w*1e-3
# Create contour plot
contour_levels = [bmin_vals0]

vmin1 , vmax1 = np.nanmin(X_vals), np.nanmax(X_vals)
bounds1 =  [vmin1,0.1,bmin_vals0,0.5,1]
## bounds before was np.log10(0.1) was too big 
bounds1 =   np.logspace(np.log10(0.05), np.log10(vmax1) , 10) 
# bounds1 =   np.logspace(np.log10(vmin1), np.log10(100) , 10) 
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)
norm2 = mpl.colors.BoundaryNorm(bounds1*cte_cbar2, cmap.N) ## second colorbar with real units 

contours = plt.contour(V0_vals, theta_mrad_vals, X_vals_transpose, levels=contour_levels,lw = 0, colors='white', linestyles='dashed' )

vmin3 , vmax3 = np.nanmin(f_vals), np.nanmax(f_vals)
bounds3 =  [vmin1,0.1,bmin_vals0,0.5,1]
bounds3 =   np.linspace(vmin3, vmax3 , 10) 
norm3 = mpl.colors.BoundaryNorm(bounds3, cmap.N)

#%%

os.chdir(path_data)
print('3-Find the (V0,theta) values for a fixed b_min = %.2f' %(bmin_vals0))

if load_data == 0:
    # Access the contour data
    allsegs = contours.collections[0].get_paths()  # Assuming only one contour collection
    
    # Check the contour data is not misleading : not corresponding to NaN values  
    valid_paths = []
    for path in allsegs:
        vertices = path.vertices 
    
        for v0, theta_mrad in vertices:
     
            theta = theta_mrad*1e-3
            value_z_norm_a = bmin_vals0 + h
            value_z = function_to_be_zero(value_z_norm_a, theta, v0)
            #print(v0,theta_mrad,value_z)
            if not np.isnan(value_z) and np.abs(value_z)<tol: ## filter the extracted V0,theta to only the ones that give real results 
                valid_paths.append([v0, theta_mrad])
     
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
        
        
        np.savetxt('bmin' + label_Ee + label_txt + '.txt', X_vals_transpose, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('V0_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', V0_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('theta_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', theta_mrad_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    
        ## solutions (V0,theta) such that bmin = bmin_vals0
        np.savetxt('V0_sol_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', listx_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('theta_sol_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', listy_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)


    print('4-Add more points near the interface of nan-non nan values')
    
    from scipy.ndimage import binary_dilation
    from scipy.interpolate import griddata
    
    Z = X_vals_transpose
    X = V0_vals
    Y = theta_mrad_vals
    # Step 2: Detect NaN interface
    if X.ndim == 1 or Y.ndim == 1:
        X, Y = np.meshgrid(X, Y)
    nan_mask = np.isnan(Z)
    nan_edge = binary_dilation(nan_mask, iterations=1) & ~nan_mask
    interface_coords = np.array([X[nan_edge], Y[nan_edge]]).T
    
    # Step 3: Create refined points around interface
    # Use a small radius to add more points near the edge
    def add_local_refined_points(coords, radius=0.3, density=10):
        points = []
        for x0, y0 in coords:
            # Local fine mesh around each interface point
            x_local = np.linspace(x0 - radius, x0 + radius, density)
            y_local = np.linspace(y0 - radius, y0 + radius, density)
            Xl, Yl = np.meshgrid(x_local, y_local)
            points.append(np.column_stack([Xl.ravel(), Yl.ravel()]))
        return np.vstack(points)
    
    refined_pts = add_local_refined_points(interface_coords, radius=0.2, density=density)
    
    # Step 4: Evaluate function and apply NaN mask
    Xr, Yr = refined_pts[:, 0], refined_pts[:, 1]
    
    
    # Preallocate result arrays
    X_vals_refined = np.zeros(len(refined_pts))
    f_vals_refined = np.zeros(len(refined_pts))
    
    # Loop over refined points
    for i, (V0, theta_mrad) in enumerate(zip(Xr, Yr)):
        theta = theta_mrad * 1e-3
        try:
            x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta, V0))
            X_vals_refined[i] = x_root - h  # or whatever your transformation is
            f_vals_refined[i] = function_to_be_zero(x_root, theta, V0)
        except ValueError:
            X_vals_refined[i] = np.nan
            f_vals_refined[i] = np.nan
    
    Zr = X_vals_refined
    # Step 5: Combine original and refined points
    X_all = np.concatenate([X.ravel(), Xr])
    Y_all = np.concatenate([Y.ravel(), Yr])
    Z_all = np.concatenate([Z.ravel(), Zr])
    
    Xg, Yg = np.meshgrid(V0_vals, theta_mrad_vals)
    # Interpolate scattered data onto grid
    Zg = griddata((Xr, Yr), Zr, (Xg, Yg), method='linear')  # or 'cubic'

else:
    X_vals_transpose = np.loadtxt('bmin' + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
    V0_vals = np.loadtxt('V0_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
    theta_mrad_vals = np.loadtxt('theta_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
    ## solutions (V0,theta) such that bmin = bmin_vals0
    listx_sorted = np.loadtxt('V0_sol_for_bmin%.2f'%(bmin_vals0)  + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
    listy_sorted = np.loadtxt('theta_sol_mrad_for_bmin%.2f'%(bmin_vals0) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
     
#%%
print('5-Plot bmin(V0,theta)')

Nind = len(listy_sorted)
ind0_bmin =  int(Nind*0)
ind1_bmin =  int(Nind/40)
ind2_bmin =  int(Nind/1.005)


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
plt.plot(listx_sorted,listy_sorted,'--',color = 'green')
plt.savefig('zmin' + label_Ee + label_txt + '.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.savefig('zmin' + label_Ee + label_txt + '.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%

from matplotlib.colors import SymLogNorm
norm3 = SymLogNorm(linthresh=1e-15, linscale=1.0, vmin=vmin3, vmax=vmax3)
f_vals_transpose = np.transpose(f_vals)

print('6-Verification of solution: check that the value gives a zero of the function')
N0 = 200

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
x0_0 = bmin_vals0 + h
list_x0 = np.linspace(x0_0*0.95,x0_0*1.05,N0)
for ind_bmin in [ind0_bmin]:
    list_function_to_be_zero_fix = []
    V0_0 = listx_sorted[ind_bmin]
    theta_mrad_0 = listy_sorted[ind_bmin]  
    theta_0 = theta_mrad_0*1e-3
    for x0_v in list_x0: 
        function = function_to_be_zero(x0_v, theta_0, V0_0)
        list_function_to_be_zero_fix.append(function)

    plt.plot(list_x0, list_function_to_be_zero_fix ,'.-', label = 'index = %i'%(ind_bmin) )
    plt.plot([x0_0], [function_to_be_zero(x0_0, theta_0, V0_0)] ,'o' )
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$f(z/W)$',fontsize=tamletra,labelpad =labelpady)
# plt.yscale('log')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('function' + label_Ee + 'ind%i_bmin%.2f' %(ind_bmin,bmin_vals0) + label_txt + '.txt', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%