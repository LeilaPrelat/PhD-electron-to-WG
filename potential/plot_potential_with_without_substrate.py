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

## compared the numerical with and without a substrate for V(z)/V0
## with z>z_min = bb/2 (outside the conductor)

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
from scipy.special import ellipk, ellipe
from scipy.optimize import curve_fit

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


print('2-Run the c++ code to plot the potential of a rectangular waveguide with substrate')

def run_e_out(bb,ss,dd, N, eps):
    # Run the C++ program name e.out and save a list of y and V
    cmd = ["./e.out", str(bb), str(ss), str(dd), str(N), str(eps)]
 
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

bb = 0.5
ss = 0.1
N = 50
epsilon = 4
list_dd = [10,100,500,1000]
color1 = mcp.gen_color(cmap="hot",n=len(list_dd)+2)
labelx = r'$z/W$'
labely = r'$V/V_0$'
title1 = r'$h/W$ = %.2f, $s/W$ =%.1f, $\epsilon$ = %i' % (bb,ss,epsilon)

#%%

listV_normV0_numerical1_tot = [] ## without a substrate (the first one we made)
listV_normV0_numerical2_tot = [] ## with a substrate (the second one we made)
list_z_norm_a1_tot = []
list_z_norm_a2_tot = []
for dd in list_dd: 
    print(dd)
    list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)
    list_z_norm_a_2, listV_normV0_2 = run_e_out(bb,ss,dd,N,epsilon)
 
    listV_normV0_numerical1_tot.append(listV_normV0)
    listV_normV0_numerical2_tot.append(listV_normV0_2)
    
    list_z_norm_a1_tot.append(list_z_norm_a)
    list_z_norm_a2_tot.append(list_z_norm_a_2)
    
#%%

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
k = 0
for dd in list_dd: 
    plt.plot(list_z_norm_a1_tot[k], listV_normV0_numerical1_tot[k] ,'--',lw = 0.8 ,label = r'$d/W$ = %i' %(dd), color = color1[k])
    plt.plot(list_z_norm_a2_tot[k], listV_normV0_numerical2_tot[k] ,'-',lw = 0.8 , color = color1[k])
    k = k + 1

#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-3,frameon=0,handletextpad=0.2, handlelength=1) 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xscale('log')
plt.savefig('potential_with_epsilon_vs_dd_forbb%.2f_epsilon%i.png' %(bb,epsilon), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%
 
# def fit_ratio(z,bb,dd,A):
    
#     return A*bb*dd*z

# popt, pcov = curve_fit(fit_ratio, xdata, ydata)
 




