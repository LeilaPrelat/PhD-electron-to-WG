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

## Overlap the potential with the Delta E_perp

import subprocess
import numpy as np
import matplotlib.pyplot as plt
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

scale_log = 0
med = 0

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'Eperp_error')

#%%

print('IMPORTANT: change the lines to numero ymax=bb/2+15*dd numero ymin=-bb/2-dd inside the dy.cpp and compile it')
print('1-Run the c++ code to plot the potential of a rectangular waveguide')

def run_dy_out(bb,ss,dd, N):
    # Run the C++ program name dy.out and save a list of y and V
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

# Transverse energy spread ΔE⊥ ≈ 2 * E0 * θ * Δθ
def delta_E_perp_gaussian(E0_eV, theta, delta_theta):
    return 2 * E0_eV * theta * delta_theta  # result in eV

def E_perp_gaussian(E0_eV, theta):
    return E0_eV*(theta**2)  # result in eV

#%%

print('Define paramaters')

hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
c =  2.99792458*10**(14)      ### light velocity in micron/seg
me_c2_eV = 510998.95069       ### me*c**2 in eV

Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)    

a=500
h=300
bb = h/a ## ratio between width and height 
ss = 0.1
dd = 1/(a*1e-3)   ## distance to the plane

list_dd = [25,50]
list_dd = [12.5,25,50]
list_dd = [2.5]
N = 400
xmin = -bb/2
xmax = bb/2

theta0 = 0.5*1e-3          ## mrad
delta_theta = 0.05*1e-3     ## mrad
theta1 = 1*1e-3            ## mrad
theta1 = 5*1e-3            ## mrad
# V0_0 = 0.25              ## eV

transverse_energy_delta0 = delta_E_perp_gaussian(Ee_electron,theta0,delta_theta) 
transverse_energy0 = E_perp_gaussian(Ee_electron,theta0)

transverse_energy_delta1 = delta_E_perp_gaussian(Ee_electron,theta1,delta_theta) 
transverse_energy1 = E_perp_gaussian(Ee_electron,theta1)

V0_max = np.max([transverse_energy0,transverse_energy1])

V0_0 = np.round(V0_max*1.05,2)
# V0_0 = 2.7

# epsi1 = 1                                              ## medium where the electron travels 
# beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
# gamma_e = 1/np.sqrt(1-epsi1*beta**2)                   ## gamma lorentz \gamma_{\rm e}
# q0 = me_c2_eV*beta*gamma_e/(hb*c)
# lambdda = 2*np.pi/q0     ## microns 

#%%

print('Plot V(z) overlap with \Delta E_perp: to know if the electron is stable in the trajectory')

labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)
title0 = r'W = %i nm, $h$ = %i nm, $s$ = %i nm, $E_{\text{e}}$ = %i keV' % (a,bb*a,ss*a,Ee_electron_keV)
title1 = title0  + ', '  +  r'$V_0$ = %.3f eV' %(V0_0)

labelx = r'$b$ (nm)'
labely = r'$V$ (eV)'

listV_normV0_tot = []
listderV_normV0_tot = []
list_b_tot = []

for dd_var in list_dd:
    list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd_var,N)
    listV_normV0_tot.append(listV_normV0)

    list_b_tot.append((np.array(list_z_norm_a)-bb/2)*a)
    
    dz = (list_z_norm_a[1] - list_z_norm_a[0])*a
    listderV_normV0 = np.gradient(listV_normV0, dz)
    listderV_normV0_tot.append(listderV_normV0)
    

med = med + 1

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
for j in range(len(list_dd)):
    dd_var = list_dd[j]
    real_d_microns = dd_var*a*1e-3
    plt.plot(np.array(list_b_tot[j]), np.array(listV_normV0_tot[j])*V0_0 ,'.-' ,label = r'$d$ = %i $\mu$m' %(real_d_microns))
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
# ejey0 = np.ones(len(list_b_tot[0]))*(transverse_energy0)
ejey1 = np.ones(len(list_b_tot[0]))*(transverse_energy1)
#plt.plot(np.ones(len(listV_normV0_tot[0]))*bmin_vals0*a,- np.array(listV_normV0_tot[0])*V0_0,'--',color='black',label = r'$b_{\text{min}}$ = %i nm' %(bmin_vals0*a))

plt.plot(list_b_tot[0], ejey1,'--',color='black',label = r'$\theta = %.2f$ mrad' %(theta1*1e3))
#plt.plot(list_b_tot[0], ejey0,'-.',color='black',label = r'$E_\perp$ with $\theta = %.1f$ mrad' %(theta0*1e3))

listV = np.array(listV_normV0_tot[0])*V0_0
bottom_grey = transverse_energy1 - transverse_energy_delta1
top_grey = transverse_energy1 + transverse_energy_delta1
# Find the index of the closest value
idx_x_left = (np.abs(listV - bottom_grey)).argmin()
closest_x_left = list_b_tot[0][idx_x_left]

idx_x_right = (np.abs(listV - top_grey)).argmin()
closest_x_right = list_b_tot[0][idx_x_right]

# eje0_minus = np.array(ejey0) - transverse_energy_delta0
# eje0_max = np.array(ejey0) + transverse_energy_delta0

eje1_minus = np.array(ejey1) - transverse_energy_delta1
eje1_max = np.array(ejey1) + transverse_energy_delta1

max1 = np.max(np.array(listV_normV0_tot[0])*V0_0)
min1 = np.min(np.array(listV_normV0_tot[0])*V0_0)
delta = 0.12
# delta = 0.5

plt.fill_between(list_b_tot[0], eje1_minus, eje1_max, color = 'lightgrey',label = r'$\Delta \theta$ = %.2f mrad' %(delta_theta*1e3))
#plt.fill_between(list_b_tot[0], eje0_minus, eje0_max, color = 'lightgrey',label = r'$\Delta E_\perp$ = %.2f eV' %(transverse_energy_delta0))

plt.plot([],[], color = 'white',label = r'$b_{\text{max}}$ = %.2f nm' %(closest_x_left))
plt.plot([],[], color = 'white',label = r'$b_{\text{min}}$ = %.2f nm' %(closest_x_right))
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
#plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
# plt.yscale('log')
# plt.xscale('log')
#plt.xticks(np.arange(0,300,50))
plt.ylim(min1-delta,max1+delta)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = [0.5,0.65],markerscale=2,fontsize=tamlegend-4,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.savefig('Vy_%i.png'%(med), format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

#%%

"""
print('Plot dV(z)/dz for different d')

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
for j in range(len(list_dd)):
    dd_var = list_dd[j]
    plt.plot(list_b_tot[j], np.array(listderV_normV0_tot[j])*V0_0 ,'.-',label = r'$d/W$ = %i' %(dd_var) )
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
# plt.plot(list_z_norm_a, np.ones(len(list_z_norm_a))*bmin_vals0,'--',color='black',label = r'$b_{\text{min}}$ = %i' %(bmin_vals0*a))
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\partial_z V(z)$ (eV/nm)',fontsize=tamletra,labelpad =labelpady)
# plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
#plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('der_Vy_zoom.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()
"""

#%%

