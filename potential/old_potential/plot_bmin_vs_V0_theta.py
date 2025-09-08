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

### All lengths are given in units of the x side-length, label as "a"


## plot the z min as a function of theta, V0. solutions only above the wg 
## return the V0,theta for a fix value of z 


import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import matplotlib as mpl
import os
# from scipy.interpolate import RegularGridInterpolator
# import numpy.ma as ma

scale_log = 0

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

add_more_points_near_Nan = 0

#%%

print('IMPORTANT: change the lines to numero ymax=bb/2+15*dd numero ymin=-bb/2-dd inside the dy.cpp and compile it')
print('1-Run the c++ code to plot the potential of a rectangular waveguide')

def run_dy_out(bb,ss,dd, N):
    # Run the C++ program name dy.out and save a list of x, z, V
    #cmd = ["./d.out", str(bb), str(ss), str(dd), str(N)]
    if scale_log == 1:
        cmd = ["./dy_log.out", str(bb), str(ss), str(dd), str(N)]
    else:
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

#%%

print('Define paramaters')

me_c2_eV = 510998.95069  ## me*c**2 in eV
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)    

a=500 ## 100 400 600  

bb = 0.6 ## 0.1 0.25 0.5 0.75 1  ## bb = 0.2 goes lower than 0.1 and lower than 0.4 (not linear)
ss = 0.1 ## 10% of the width
dd = 0.001/(a*1e-3)   ## = d/a  ### a = 400 nm : 12.5 for d = 5 microns / 25 for d = 10 microns / 50 for d = 20 microns 
           ### a = 300 nm : 3.33 for d = 1 microns
           ### a = 400 nm : 2.5 for d = 1 microns
dd = 20

wp = 0.5 ## w'/w
N = 200   ## discretization points.
epsilon=9 ## permittivity = 9
N = 200   ## discretization points.
x0 = 0    ## x0/w
x1 = 2    ##x1/w
nx = 200 
z0 = bb    ## z0/w ## start outside the waveguide
z1 = 2    ## z1/w 
nz = 200
# x=0 is in the middle of the center waveguide
# z=0 is at the interface


N = 500
xmin = -bb/2
xmax = bb/2

if dd <= 100: 
    Nvals = 400
else:
    Nvals = 500

if dd<= 100:
    # Define ranges for V0 and theta
    V0_vals = np.linspace(0.005, 5, Nvals)
    V0_vals = np.linspace(0.1, 3, Nvals)
    theta_mrad_vals = np.linspace(0.01, 5, Nvals)
    
    theta_mrad_vals = np.linspace(0.1, 2.5, Nvals)
    V0_vals = np.linspace(0.01, 1, Nvals)
    
    
    theta_mrad_vals = np.linspace(0.1, 5, Nvals)
    V0_vals = np.linspace(0.01, 2.5, Nvals)
    
    
    V0_0 = 0.2
    theta_mrad0 = 1.5
    
else:
    theta_mrad_vals = np.linspace(2, 6, Nvals)
    V0_vals = np.linspace(1, 6, Nvals)
    V0_0 = 3
    theta_mrad0 = 4
    
if bb == 1: 
    bmin_vals0 = 0.15
    # bmin_vals0 = 0.1
elif bb == 0.75: 
    bmin_vals0 = 0.15    ## same bmin = 60 nm


bmin_vals0 = 50/a  ## 50 nm
# bmin_vals0 = 60/a
bmin_vals0 = 80/a

Nvals = 300
theta_mrad_vals = np.linspace(0.01, 2.5, Nvals)
V0_vals = np.linspace(0.01, 0.6, Nvals)

    # bmin_vals0 = 0.1
label_txt = '_dd%.2f_hh%.2f.txt' %(dd,bb)

tol = 1e-6 ## zero of the function 
tol = 1e-1 ## zero of the function 

density = 50 ## add more points to the interface of (nan-not nan) values 

#%%

labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)
title0 = r'W = %i nm, $h/W$ = %.2f, $s/W$ = %.1f' % (a,bb,ss)
title1 = title0  + ', '  +  r'$V_0$ = %.2f eV' %(V0_0)

#%%

list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)
Nint = 1
## the electron is between the gate and the rectangular waveguide. NO
## the electron after the WG 
list_z_norm_a_interp = np.linspace(np.min(list_z_norm_a), np.max(list_z_norm_a),int(N*Nint)) 
V_interp =  CubicSpline(list_z_norm_a, listV_normV0)

#%%

print('2-Define the function to be zero--> zmin (electron position) for a fixed V0 and different angles')

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
    
    V_norm_V0 = V_interp(value_z_norm_a)
    
    partA = function/V0
    partB = V_norm_V0
    # print(partA/partB)
    return partA - partB ## add a minus to the potential of the WG --> has to be negative 

def analytical_function_bmin( theta, V0):
    """
    Parameters
    ----------
    theta : angle of the electron (initial conditions) in radians
    V0 : potential of the gate in e*Volts
    Returns
    -------
    bmin/width: See figure from folder 2025/06/05
    """
    epsi1 = 1
    beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
    gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
    
    eta=(1+1/dd)/4
    eta=1
    
    factor =  me_c2_eV*gamma_e/2
    factor_theta = beta*np.sin(theta)
    factor_log = np.log(bb*0.5*eta/dd)/V0
    
    function = factor*factor_log*factor_theta**2
    
    bmin = dd*np.exp(function)/eta-bb/2
    return bmin

#%%

print('3-Color map of zmin as a function of V0 and theta, for zmin (electron position) above the WG')

# theta_mrad_vals = np.linspace(0.1, 1, Nvals)
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
        except ValueError or x_root == bb/2:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan

#%%
X_vals_transpose = np.transpose(X_vals)

tamfig2 = [5, 3]
limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
cmap = plt.cm.RdBu   # define the colormap

# another cbar with constant a*1e-3 (real units, in microns)
cte_cbar2 = a*1e-3
# Create contour plot
contour_levels = [bmin_vals0]

vmin1 , vmax1 = np.nanmin(X_vals), np.nanmax(X_vals)
bounds1 =  [vmin1,0.1,bmin_vals0,0.5,1]
bounds1 =   np.logspace(np.log10(0.1), np.log10(100) , 10) 
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

# Access the contour data
allsegs = contours.collections[0].get_paths()  # Assuming only one contour collection

# Check the contour data is not misleading : not corresponding to NaN values  
valid_paths = []
for path in allsegs:
    vertices = path.vertices 

    for v0, theta_mrad in vertices:
 
        theta = theta_mrad*1e-3
        value_z_norm_a = bmin_vals0 + bb/2
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
    
    header = title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV)
    
    
    np.savetxt('bmin' + label_Ee + label_txt, X_vals_transpose, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('V0_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt, V0_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('theta_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt, theta_mrad_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)

    ## solutions (V0,theta) such that bmin = bmin_vals0
    np.savetxt('V0_sol_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt, listx_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('theta_sol_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt, listy_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)


#%%

if add_more_points_near_Nan == 1: 
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
            X_vals_refined[i] = x_root - bb / 2  # or whatever your transformation is
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
    
    Zfinal = Zg
else:
    Zfinal = X_vals_transpose

#%%

Nind = len(listy_sorted)
ind0_bmin =  int(Nind*0)
ind1_bmin =  int(Nind/40)
ind2_bmin =  int(Nind/1.005)


plt.figure(figsize=tamfig2)
#plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
# im_show = plt.imshow(X_vals_transpose, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' , norm = norm1  ) 
im_show = plt.imshow(Zfinal, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' , norm = norm1  ) 
    #plt.clabel(contours, fmt='%.2f', colors='green', fontsize=tamletra, manual=[(0.5, 5) ])  # Label contours
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.18 , format = '%.2f') 
im_show2 = plt.imshow(np.array(X_vals_transpose)*cte_cbar2, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm=norm2  )  ## second colorbar with real units 
cbar2 = plt.colorbar(im_show2, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.set_title(r'$b_{\text{min}}$ ($\mu$m)',fontsize=tamletra)
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.plot(listx_sorted,listy_sorted,'--',color = 'green')
# plt.plot(listx_sorted[ind0_bmin],listy_sorted[ind0_bmin],'o',color = 'purple')
# plt.plot(listx_sorted[ind1_bmin],listy_sorted[ind1_bmin],'o',color = 'purple')
# plt.plot(listx_sorted[ind2_bmin],listy_sorted[ind2_bmin],'o',color = 'purple')
plt.savefig('zmin' + label_Ee + '_bmin%.2f_dd%.2f_hh%.2f.png' %(bmin_vals0,dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.savefig('zmin' + label_Ee + '_bmin%.2f_dd%.2f_hh%.2f.pdf' %(bmin_vals0,dd,bb), format='pdf',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%
from matplotlib.colors import SymLogNorm
norm3 = SymLogNorm(linthresh=1e-15, linscale=1.0, vmin=vmin3, vmax=vmax3)
f_vals_transpose = np.transpose(f_vals)

N0 = 200

plt.figure(figsize=tamfig)
x0_0 = bmin_vals0 + bb/2
list_x0 = np.linspace(x0_0*0.95,x0_0*1.05,N0)
for ind_bmin in [ind0_bmin]:
    list_function_to_be_zero_fix = []
    V0_0 = listx_sorted[ind_bmin]
    theta_mrad_0 = listy_sorted[ind_bmin]  
    theta_0 = theta_mrad_0*1e-3
    for x0 in list_x0: 
        function = function_to_be_zero(x0, theta_0, V0_0)
        list_function_to_be_zero_fix.append(function)

    plt.plot(list_x0, list_function_to_be_zero_fix ,'.-', label = 'index = %i'%(ind_bmin) )
    plt.plot([x0_0], [function_to_be_zero(x0_0, theta_0, V0_0)] ,'o' )
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$f(z/W)$',fontsize=tamletra,labelpad =labelpady)
# plt.yscale('log')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('function' + label_Ee + 'ind%i_bmin%.2f_dd%.2f_hh%.2f.png' %(ind_bmin,bmin_vals0,dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%
"""

V0 = Ee_electron/1e4 
V0 = 0.1
aux_ejey = np.linspace(-np.max(listV_normV0),np.min(listV_normV0),10)

list_theta_mrad = [0.1,0.25,0.5,0.75,1,2]

from mycolorpy import colorlist as mcp
color1 = mcp.gen_color(cmap="hot",n=len(list_theta_mrad)+2)

plt.figure(figsize=tamfig)
plt.title(title0 + r', $V_0$ = %.1f eV, $E_{\text{e}}$ = %i keV' %(V0,Ee_electron_keV),fontsize=tamtitle)
k = 0
for theta_mrad in list_theta_mrad:
    theta = theta_mrad*1e-3
    list_function_to_be_zero = [] 
    for value_z_norm_a in list_z_norm_a_interp:
        list_function_to_be_zero.append(function_to_be_zero(value_z_norm_a,theta,V0))
    
 
    sol2 = root(function_to_be_zero, [bb/2], method='hybr', args=(theta,V0))   ## electron above the WG
    
    
    plt.plot(list_z_norm_a_interp, list_function_to_be_zero ,'.-' ,color = color1[k],label = r'$\theta$ = %.2f mrad' %(theta_mrad))
    print( sol2.x)
 
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
"""

#%%