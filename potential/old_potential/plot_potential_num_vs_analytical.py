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

## compared the numerical and analytical (wire) expressions for V(z)/V0
## with z>z_min = bb/2 (outside the conductor)

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
from scipy.special import ellipk, ellipe
from scipy.optimize import curve_fit

def derivative_K_m(m):
    """
    Calculates the derivative of the complete elliptic integral of the first kind
    with respect to its parameter m.

    Args:
        m (float): The parameter of the elliptic integral (0 <= m < 1).

    Returns:
        float: The derivative of K(m) with respect to m.
    """

    K = ellipk(m)
    E = ellipe(m)
    return (E - (1 - m**2) * K) / (2 * m * (1 - m**2)**0.5)

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

def V_analytical(z_norm_a,dd,bb): ## V/V0
    image_dist = 2*dd  ## distance to the image conductor 
    zref = bb/2
    factor = 1+(bb/dd)*np.pi
    # if z_norm_a <= 1:
    #     factor = (1+1/dd)/4 ## for bb = 0.1 this works at large z but dont work for bb = 0.5
    # else:
    #     factor = np.exp((1+1/dd)*z_norm_a)
    
    factor = (1+1/dd)/4 ## for bb = 0.1 this works at large z but dont work for bb = 0.5
    # factor = bb/(2*np.pi) ## thick so bb big
    # factor = 1/2          ## thin  so bb small
    
    arg = 1/np.sqrt(1+bb**2)
    K =  ellipk(arg)
    derK = derivative_K_m(arg)
    
    factor = np.exp(-np.pi*0.5*derK/K)*4/np.pi
 
    #factor_bb = 13.37*bb ## part of the function that depends on the ratio between height and witdh, label as bb
    
    factor = 1
    
    Vref = np.log(zref*factor/dd) ## value inside conductor . we force the potential to be 1, reference value
    

    function = np.log(z_norm_a*factor/dd)/Vref ## normalize to the value of the conductor 
    
    delta = (1 + bb**2)**(-1/2)/np.pi
    
    numerator = 4*(z_norm_a - bb/2)**2 + delta**2
    
    
    return function

bb = 0.75
ss = 0.1
N = 50
list_dd = [100,500,1000,2000]
color1 = mcp.gen_color(cmap="hot",n=len(list_dd)+2)
labelx = r'$z/W$'
labely = r'$V/V_0$'
title1 = r'$h/W$ = %.2f, $s/W$ =%.1f' % (bb,ss)

#%%

listV_normV0_analytical_tot = []
listV_normV0_numerical_tot = []
for dd in list_dd: 
    print(dd)
    list_z_norm_a, listV_normV0 = run_dy_out(bb,ss,dd,N)
    listV_normV0_analytical = [] 
    for zz in list_z_norm_a:
        listV_normV0_analytical.append(V_analytical(zz,dd,bb))
    listV_normV0_analytical_tot.append(listV_normV0_analytical)
    listV_normV0_numerical_tot.append(listV_normV0)
    
#%%

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
k = 0
for dd in list_dd: 
    plt.plot(list_z_norm_a, listV_normV0_analytical_tot[k] ,'--',lw = 0.8 ,label = r'$d/W$ = %i' %(dd), color = color1[k])
    plt.plot(list_z_norm_a, listV_normV0_numerical_tot[k] ,'-',lw = 0.8 , color = color1[k])
    k = k + 1

#plt.plot(listx[cut:-1], np.array(listy_analytical[cut:-1]),'--',color = 'red')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-3,frameon=0,handletextpad=0.2, handlelength=1) 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xscale('log')
plt.savefig('potential_vs_dd_forbb%.2f.png' %(bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%

ratio1 = np.array(listV_normV0_analytical_tot[-1]/np.array(listV_normV0_numerical_tot[-1]))
ratio2 = np.array(listV_normV0_analytical_tot[-2]/np.array(listV_normV0_numerical_tot[-2]))
ratio3 = np.array(listV_normV0_analytical_tot[-3]/np.array(listV_normV0_numerical_tot[-3]))


labely = r'$V_{analytical}/V_{num}$'
plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_z_norm_a, np.array(ratio1) ,'--',lw = 0.8 , color = color1[1],label = r'$d/W$ = %i' %(list_dd[1])  )
plt.plot(list_z_norm_a, np.array(ratio2) ,'--',lw = 0.8 , color = color1[2] ,label = r'$d/W$ = %i' %(list_dd[2]) )
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend-3,frameon=0,handletextpad=0.2, handlelength=1) 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xscale('log')
plt.savefig('ratio_vs_dd_forbb%.2f.png' %(bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
# plt.plot(list_z_norm_a, np.array(ratio3) ,'--',lw = 0.8  )
plt.show()

#%%




# def fit_ratio(z,bb,dd,A):
    
#     return A*bb*dd*z

# popt, pcov = curve_fit(fit_ratio, xdata, ydata)
 




