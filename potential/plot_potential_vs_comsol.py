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
 
## plot V/V0(x,z) for different d/W, h/W, wp/W and different x/W

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
import os

color1 = mcp.gen_color(cmap="hot",n=6)

tamfig = [4, 3]
tamletra = 13
tamtitle  = tamletra - 1
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
path_data = os.path.join(path_basic, 'COMSOL')
 
#%%

print('Run the c++ code to get the points of the potential')

def run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz): ## fixed x 
    os.chdir(path_basic)
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

W = 500 ## nm
Wp = 250 ## nm
G = 100  ## nm. same as my "d" without normalized by width
H = 100 ## height

list_x0 = [0,100,200,400]   ## nm
# list_x0 = [100,200,400]   ## nm
# list_x0 = [0,100,200]   ## nm
list_Vair = [5,10]   ## boundary condition placed at 5 microns or 10 microns
Vair = list_Vair[0]  ## boundary condition placed at 5 microns or 10 microns

os.chdir(path_data)
list_comsol_z_norm_tot = []
list_Vcomsol_tot = []

# x0=0
# data0 = np.loadtxt('x%i_W%i_G%i_Wp%i_V%ium_Eps-10000.txt' %(x0,W,G,Wp,Vair), comments='%')
# z_nm = data0[:, 0]*1e9
# dependent_u = data0[:, 1]
# z_nm_norm = np.array(z_nm)/W
# list_comsol_z_norm_tot.append(z_nm_norm)
# list_Vcomsol_tot.append(dependent_u)

# list_x0 = [0]   ## nm


labelx = r'$z/W$'
labely = r'$\phi(x,z)/V_0$'

print('Import data from Comsol')
for x0 in list_x0:
    data = np.loadtxt('x%i_W%i_G%i_Wp%i_V%ium.txt' %(x0,W,G,Wp,Vair), comments='%')
    z_nm = data[:, 0]*1e9
    dependent_u = data[:, 1]
    z_nm_norm = np.array(z_nm)/W
    list_comsol_z_norm_tot.append(z_nm_norm)
    list_Vcomsol_tot.append(dependent_u)
    

#%%  

wp = Wp/W
d = G/W
h = H/W ## ratio
list_x0_norm = np.array(list_x0)/W

print('1-Plot V(z) for different x positions for fixed wp/W = %.2f, d/W = %.2f, h/W = %.2f' %(wp,d,h))
 
epsilon = 11.7 ## permittivity = 9
N  = 200       ## discretization points.
z0 = h*0.05    ## z0/w ## start outside the waveguide
z1 = z_nm_norm[-1]    ## z1/w 

# x=0 is in the middle of the center waveguide
# z=0 is at the interface 

title1 = r'$W$ = %i nm, $z_{\text{air}}$ = %i $\mu$m' %(W,Vair) + '\n' + r'$W_{\text{p}}/W$ = %.2f, $d/W$ = %.2f, $h/W$ = %.2f' %(wp,d,h)

listV_normU_tot_1 = []
list_z_norm_w_tot1 = []
# list_normV0 = []
j=0
for x0 in list_x0_norm:
    nz = len(list_Vcomsol_tot[j])
    list_z_norm_w, listV_normV0 = run_e_out1(wp,d,h,N,epsilon,x0,z0,z1,nz)
    listV_normU_tot_1.append(np.array(listV_normV0)) ## np.array(listV_normV0)*2+0.3
    
    # V0_over_U = listV_normV0[0]
    # list_normV0.append(normV0)
    list_z_norm_w_tot1.append(list_z_norm_w)
    j=j+1
    
zz, V0_over_U_list = run_e_out1(wp,d,h,N,epsilon,0,h,h,2)
V0_over_U = np.abs(V0_over_U_list[0]) ## normalized to V0

#%%

## mapping my results with Saad's by adding + (-V0_over_U+2)/2
plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
for j in range(len(list_x0_norm)):
    x0 = list_x0_norm[j]
    # norm_V0 = np.abs(list_normV0[j])
#    plt.plot(np.array(list_z_norm_w_tot1[j]), np.array(listV_normU_tot_1[j])/V0_over_U+1 ,'-',color = color1[j] ,label = r'$x/W$ = %.1f' %(x0))
#    plt.plot(np.array(list_comsol_z_norm_tot[j]), (np.array(list_Vcomsol_tot[j])+1)/2 ,'--',color = color1[j])
    plt.plot(np.array(list_z_norm_w_tot1[j]), np.array(listV_normU_tot_1[j]) + (-V0_over_U+2)/2 ,'-',color = color1[j] ,label = r'$x/W$ = %.1f' %(x0))
    plt.plot(np.array(list_comsol_z_norm_tot[j]), np.array(list_Vcomsol_tot[j]) ,'--',color = color1[j])
    
    
plt.plot([],[],'--',color='black',label='COMSOL')
plt.plot([],[],'-',color='black',label='BEM')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamlegend-2,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
os.chdir(path_data)
# plt.grid(1)
plt.savefig('Vz_vs_x_Vair%ium.png' %(Vair), format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vz_vs_x_Vair%ium.pdf' %(Vair), format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

del x0

#%%

labely = r'$V_{\text{comsol}}/V_{\text{code}}$'

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
for j in range(len(list_x0_norm)):
    x0 = list_x0_norm[j]
    # norm_V0 = np.abs(list_normV0[j])
    ratio =  np.array(list_Vcomsol_tot[j])/np.array(listV_normU_tot_1[j]/V0_over_U)
    plt.plot(np.array(list_comsol_z_norm_tot[j]), ratio ,'-',color = color1[j],label = r'$x/W$ = %.1f' %(x0))
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamlegend-2,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
os.chdir(path_data)
# plt.grid(1)
plt.savefig('Vz_vs_x_ratio_Vair%ium.png' %(Vair), format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()