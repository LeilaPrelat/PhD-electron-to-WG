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

## plot V/V0 for different d/W and h/W

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
color1 = mcp.gen_color(cmap="hot",n=7)

tamfig = [4, 3]
tamletra = 13
tamtitle  = tamletra - 2
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


bb0 = 1 ## ratio between width and height 
dd0 = 0.2

list_dd = [0.2,0.5,1,2,5]
list_bb = [0.2,0.5,1,2,5]

N = 400

ss = 0.1 ## 10% of W
 
#%%
    
print('Plot V(z) for different d/W for a fixed h/W = %.2f' %(bb0))

labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)
title0 = r'h/W = %.2f, $s/W$ = %.1f' % (bb0,ss)

listV_normV0_tot = []
list_b_tot = []

for dd_var in list_dd:
    list_z_norm_a, listV_normV0 = run_dy_out(bb0,ss,dd_var,N)
    listV_normV0_tot.append(listV_normV0)
    list_b_tot.append((np.array(list_z_norm_a)-bb0/2))

#%%

plt.figure(figsize=tamfig)
plt.title(title0,fontsize=tamtitle)
for j in range(len(list_dd)):
    dd_var = list_dd[j]
    plt.plot(np.array(list_b_tot[j]), np.array(listV_normV0_tot[j]) ,'.-',color = color1[j] ,label = r'$d/W$ = %.1f' %(dd_var))
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
#plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
# plt.yscale('log')
# plt.xscale('log')
plt.xticks(np.arange(0,0.3,0.05))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamlegend-2,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.savefig('Vy_vs_d.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

#%%

print('Plot V(z) for different h/W for a fixed d/W = %.2f' %(dd0))

labelx = r'$z/W$'
labely = r'$V/V_0$'

# aux_ejey = np.linspace(np.min(listV_normV0),np.max(listV_normV0),10)
title0 = r'd/W = %.2f, $s/W$ = %.1f' % (dd0,ss)

listV_normV0_tot2 = []
list_b_tot2 = []

for bb_var in list_bb:
    list_z_norm_a2, listV_normV02 = run_dy_out(bb_var,ss,dd0,N)
    listV_normV0_tot2.append(listV_normV02)
    list_b_tot2.append((np.array(list_z_norm_a2)-bb_var/2))

#%%

plt.figure(figsize=tamfig)
plt.title(title0,fontsize=tamtitle)
for j in range(len(list_bb)):
    bb_var = list_bb[j]
    plt.plot(np.array(list_b_tot2[j]), np.array(listV_normV0_tot2[j]) ,'.-',color = color1[j] ,label = r'$h/W$ = %.1f' %(bb_var))
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.plot(np.ones(10)*xmin,aux_ejey,'--',color = 'black')
#plt.plot(np.ones(10)*xmax,aux_ejey,'--',color = 'black')
# plt.yscale('log')
# plt.xscale('log')
plt.xticks(np.arange(0,0.3,0.05))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamlegend-2,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.savefig('Vy_vs_ratio.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

#%%

