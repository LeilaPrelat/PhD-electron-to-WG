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

os.chdir(path_basic)
    
print('Run the c++ code to get the points of the potential')

def run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz): ## fixed x 
    x1 = x0
    nx = 1
    # Run the C++ program name dy.out and save a list of x, z, V
    cmd = ["./e.out", str(wp), str(d), str(h), str(N), str(epsilon), str(x0), str(x1), str(nx), str(z0), str(z1), str(nz)]
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = result.stdout.strip().split('\n')
    
    list_z_norm_w = []
    list_V_norm_u = []
    for line in lines:
        try:
            zvalue, Vvalue = map(float, line.strip().split())
            list_z_norm_w.append(zvalue)
            list_V_norm_u.append(Vvalue)  
        except ValueError:
            continue
        
    return list_z_norm_w, list_V_norm_u



# Transverse energy spread ΔE⊥ ≈ 2 * E0 * θ * Δθ
def delta_E_perp_gaussian(E0_eV, theta, delta_theta):
    return 2 * E0_eV * theta * delta_theta  # result in eV

def E_perp_gaussian(E0_eV, theta):
    return -E0_eV*(theta**2)  # result in eV

#%%

print('Define paramaters')

me_c2_eV = 510998.95069  ## me*c**2 in eV
Ee_electron_keV = 100
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)    

wp = 0.5 ## w'/w
d = 0.2  ## d/w distance between the side wires normalized to w
h = 0.6   ## aspect ratio height/w
N = 400   ## discretization points.
epsilon=9 ## permittivity = 9
N = 200   ## discretization points.
x0 = 0.4     ## x0/w --> position of the electron. options [0,0.2,0.4,0.6]. 
            ## x0=/0.6/0.8 does not give solutions because potential stops being negative: check the plots of V vs x0
z0 = h*1.001    ## z0/w ## start outside the waveguide
z1 = 2    ## z1/w 
nz = 200
# x=0 is in the middle of the center waveguide
# z=0 is at the interface


delta_theta = 0.05*1e-3     ## mrad
theta1 = 1*1e-3            ## mrad
theta1 = 3*1e-3            ## mrad
# V0_0 = 0.25              ## eV

transverse_energy_delta1 = delta_E_perp_gaussian(Ee_electron,theta1,delta_theta) 
transverse_energy1 = E_perp_gaussian(Ee_electron,theta1)

V0_max = transverse_energy1

V0_0 = np.abs(np.round(V0_max*1.5,2))
print(V0_0)
# V0_0 = 1
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f_Ee%ikeV' %(wp,h,d,x0,Ee_electron_keV) ## all normalize to width w

 

#%%

print('Plot V(z) overlap with \Delta E_perp: to know if the electron is stable in the trajectory')

labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)
title1 = r'$W_{\text{p}}/W$ = %.1f, $h/W$ = %.1f, $d/W$ = %.1f' %(wp,h,d)
title2 = r'x/W = %.1f, $U$ = %.1f eV, $E_{\text{e}}$ = %i keV' %(x0,V0_0,Ee_electron_keV)
title = title1  + '\n'  +  title2

labelx = r'$b$ (nm)'
labely = r'$V$ (eV)'

zz, V0_over_U_list = run_e_out1(wp,d,h,N,epsilon,0,h,h,2)
V0_over_U = np.abs(V0_over_U_list[0]) ## normalized to V0
list_z_norm_a, list_V_norm_u = run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz)
list_b_norm_a = np.array(list_z_norm_a)-h

listV_normV0 = np.array(list_V_norm_u)/V0_over_U

med = med + 1

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.plot(np.array(list_b_norm_a), np.array(listV_normV0)*V0_0 ,'.-')
#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
# ejey0 = np.ones(len(list_b_tot[0]))*(transverse_energy0)
ejey1 = np.ones(len(list_b_norm_a))*(transverse_energy1)
#plt.plot(np.ones(len(listV_normV0_tot[0]))*bmin_vals0*a,- np.array(listV_normV0_tot[0])*V0_0,'--',color='black',label = r'$b_{\text{min}}$ = %i nm' %(bmin_vals0*a))

plt.plot(list_b_norm_a, ejey1,'--',color='black',label = r'$\theta = %.2f$ mrad' %(theta1*1e3))
#plt.plot(list_b_tot[0], ejey0,'-.',color='black',label = r'$E_\perp$ with $\theta = %.1f$ mrad' %(theta0*1e3))
listV = np.array(listV_normV0)*V0_0
bottom_grey = transverse_energy1 - transverse_energy_delta1
top_grey = transverse_energy1 + transverse_energy_delta1
# Find the index of the closest value
idx_x_left = (np.abs(listV - bottom_grey)).argmin()
closest_x_left = list_b_norm_a[idx_x_left]

idx_x_right = (np.abs(listV - top_grey)).argmin()
closest_x_right = list_b_norm_a[idx_x_right]

# eje0_minus = np.array(ejey0) - transverse_energy_delta0
# eje0_max = np.array(ejey0) + transverse_energy_delta0

eje1_minus = np.array(ejey1) - transverse_energy_delta1
eje1_max = np.array(ejey1) + transverse_energy_delta1

max1 = np.max(np.array(listV_normV0)*V0_0)
min1 = np.min(np.array(listV_normV0)*V0_0)
delta = 0.12
# delta = 0.5

plt.fill_between(list_b_norm_a, eje1_minus, eje1_max, color = 'lightgrey',label = r'$\Delta \theta$ = %.2f mrad' %(delta_theta*1e3))
#plt.fill_between(list_b_tot[0], eje0_minus, eje0_max, color = 'lightgrey',label = r'$\Delta E_\perp$ = %.2f eV' %(transverse_energy_delta0))

plt.plot([],[], color = 'white',label = r'$b_{\text{max}}$ = %.2f nm' %(closest_x_left))
plt.plot([],[], color = 'white',label = r'$b_{\text{min}}$ = %.2f nm' %(closest_x_right))
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
#plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
# plt.yscale('log')
# plt.xscale('log')
plt.xticks(np.arange(0,1.7,0.3))
plt.ylim(min1-delta,max1+delta)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-4,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
os.chdir(path_data)
plt.savefig('Vz' + label_txt + '.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

#%%
 

