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

Nvals = 400
me_c2_eV = 510998.95069  ## me*c**2 in eV
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)    

bb = 1
if bb== 1:
    ss = 20/300 ## W = 300
else:
    ss = 20/400 ## W = 400
dd = 5
N = 200
xmin = -bb/2
xmax = bb/2


labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)

list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)

title0 = r'$h/W$ = %.2f, $s/W$ =%.1f' % (bb,ss)
plt.figure(figsize=tamfig)
plt.title(title0,fontsize=tamtitle)
plt.plot(list_z_norm_a, listV_normV0 ,'.-',label = r'$d/W$ = %i' %(dd) )
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
#plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('Vy.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

Nint = 1
## the electron is between the gate and the rectangular waveguide. NO
## the electron after the WG 
list_z_norm_a_interp = np.linspace(np.min(list_z_norm_a), np.max(list_z_norm_a),int(N*Nint)) 
V_interp =  CubicSpline(list_z_norm_a, listV_normV0)

title1 = r'$h/W$ = %.2f, $s/W$ =%.1f, $d/W$ = %i' % (bb,ss,dd)

#%%

print('2-Define the function to be zero--> zmin (electron position) for a fixed V0 and different angles')

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
V0 = Ee_electron/1e4 
V0 = 0.1
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

aux_ejey = np.linspace(-np.max(listV_normV0),np.min(listV_normV0),10)

list_theta_mrad = [0.1,0.25,0.5,0.75,1,2]
theta_mrad0 = list_theta_mrad[2]

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

#%%

print('3-Color map of zmin as a function of V0 and theta, for zmin (electron position) above the WG')

# Define ranges for V0 and theta
V0_vals = np.linspace(0.005, 5, Nvals)
V0_vals = np.linspace(0.1, 3, Nvals)
V0_vals = np.linspace(0.01, 0.25, Nvals)
theta_mrad_vals = np.linspace(0.01, 5, Nvals)
theta_mrad_vals = np.linspace(0.001, 1, Nvals)
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
        except ValueError:
            # No root found in the interval
            X_vals[i, j] = np.nan
            f_vals[i, j] = np.nan

#%%

tamfig2 = [5, 3]
n_color = 10
vmin1 , vmax1 = np.nanmin(X_vals), np.nanmax(X_vals)
cmap = plt.cm.RdBu  # define the colormap
bounds1 =   np.linspace(vmin1, vmax1 , n_color) 
if dd == 20:
    bounds1 =  [vmin1,0.25,0.5,5,100,200,int(vmax1)]
    # bounds1 =  [vmin1,0.25,0.5,5,int(vmax1)]
else:
    bounds1 =  [vmin1,0.25,0.5,5,50,int(vmax1)]
    bounds1 =  [vmin1,0.1,0.25,0.5,1]
    bounds1 =   np.logspace(np.log10(0.1), np.log10(100) , 10) 

width=300
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)
norm2 = mpl.colors.BoundaryNorm(bounds1*width, cmap.N)
# Mark specific values (e.g., contours at Z = 0.2, 0.4, 0.6)
contour_levels = [0.15,   0.5]

plt.figure(figsize=tamfig2)
#plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
im_show = plt.imshow(np.transpose(X_vals), extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm=norm1  ) 
# contours = plt.contour(V0_vals, theta_mrad_vals, np.transpose(X_vals), levels=contour_levels, colors='green', linestyles='dashed')
# plt.clabel(contours, fmt='%.2f', colors='green',fontsize=tamletra)  # Label contours

cbar = plt.colorbar(im_show, fraction=0.046, pad=0.15   ,format = '%.2f') 
im_show2 = plt.imshow(np.transpose(X_vals)*width, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm=norm2  ) 
cbar2 = plt.colorbar(im_show2, fraction=0.046, pad=0.04, orientation = 'vertical')

cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra-1)
cbar.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.set_title(r'$b_{\text{min}}$ (nm)',fontsize=tamletra-1)
# plt.xticks(np.arange(0,np.max(V0_vals)+0.5,0.5))
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin2' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

header = title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV)

np.savetxt('bmin' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), np.transpose(X_vals), fmt='%.10f', delimiter='\t', header = header, encoding=None)
np.savetxt('V0' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), V0_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)
np.savetxt('theta_mrad' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), theta_mrad_vals, fmt='%.10f', delimiter='\t', header = header, encoding=None)

#%%

X_vals_analytical = np.zeros((len(V0_vals), len(theta_mrad_vals)))

for i, V0 in enumerate(V0_vals):
    for j, theta_mrad in enumerate(theta_mrad_vals):
        theta = theta_mrad*1e-3
        value_analytical = analytical_function_bmin(theta, V0)
        X_vals_analytical[i, j] = value_analytical
        
 #%%      

""" 
plt.figure(figsize=tamfig)
plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
im_show = plt.imshow(np.transpose(X_vals_analytical), extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'    ) 
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%.2f') 
cbar.ax.set_title(r'$b^{ana}_{min}/a$',fontsize=tamletra)
plt.xticks(np.arange(0,6,1))
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
#plt.savefig('zmin_analytical' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()
"""
 #%%  
print('3-zmin as a function of V0 for theta=5mrad, for zmin (electron position) above the WG')

zmin = bb/2 ## above the WG
zmax = np.max(list_z_norm_a)
# Storage for roots
X_vals_2 = np.zeros(len(V0_vals))
f_vals_2 = np.zeros(len(V0_vals))
# Loop over parameter values
for i, V0 in enumerate(V0_vals):
    theta_mrad = theta_mrad0
    theta = theta_mrad*1e-3
    try:
        # Solve f(z,theta, V0) = 0 in some interval z in [zmin, zmax]
        x_root = brentq(function_to_be_zero, zmin, zmax, args=(theta,V0))
        X_vals_2[i] = x_root-bb/2 ## we subtract the height of the waveguide
        value_function = function_to_be_zero(x_root, theta, V0)
        f_vals_2[i] = value_function
    except ValueError:
        # No root found in the interval
        X_vals_2[i] = np.nan
        f_vals_2[i] = np.nan
           
        #%%  
plt.figure(figsize=tamfig)
plt.title(title1 + r', $\theta$ = %.2f mrad, $E_{\text{e}}$ = %i keV' %(theta_mrad,Ee_electron_keV),fontsize=tamtitle)    
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$b_{\text{min}}/W$',fontsize=tamletra,labelpad =labelpady)
plt.plot( V0_vals, X_vals_2)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-5,frameon=0,handletextpad=0.2, handlelength=1) 
#plt.savefig('minimum_z_vs_V0_dd%i_hh%.2f' %(dd,bb) + label_Ee + '_V0_%ieV.png' %(V0), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()
   
        
        
        
        