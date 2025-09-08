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
import matplotlib.tri as tri
import matplotlib as mpl
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import curve_fit
import os

color1 = mcp.gen_color(cmap="hot",n=8)

tamfig = [4, 3]
tamfigure2 = [3.6,3]
tamfig2 = [4.4, 3]
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
 
path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'potential_vs_diff_parameters')

# Define linear model
def linear_func(z_norm_w, a, b):
    return a*z_norm_w + b

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


def run_e_out2(wp,d,h,N,epsilon,x0,x1,nx,z0,z1,nz):
    os.chdir(path_basic)
    # Run the C++ program name dy.out and save a list of x, z, V
    cmd = ["./e.out", str(wp), str(d), str(h), str(N), str(epsilon), str(x0), str(x1), str(nx), str(z0), str(z1), str(nz)]
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = result.stdout.strip().split('\n')
 
    list_x = []
    list_z = []
    list_V = []
    for line in lines:
        try:
            xvalue, zvalue, Vvalue = map(float, line.strip().split())
            list_x.append(xvalue)
            list_z.append(zvalue)
            list_V.append(Vvalue)  
        except ValueError:
            continue

    return list_x, list_z, list_V


def run_e_out0(wp,d,h,N,epsilon,z0,x0,x1,nx): ## fixed x 
    os.chdir(path_basic)
    z1 = z0
    nz = 1
    # Run the C++ program name dy.out and save a list of x, z, V
    cmd = ["./e.out", str(wp), str(d), str(h), str(N), str(epsilon), str(x0), str(x1), str(nx), str(z0), str(z1), str(nz)]
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = result.stdout.strip().split('\n')
    
    
    
    list_x = []
    list_V = []
    for line in lines:
        try:
            xvalue, Vvalue = map(float, line.strip().split())
            list_x.append(xvalue)
            list_V.append(Vvalue )  ##+ 0.778
        except ValueError:
            continue

    return list_x, list_V




#%%

print('Define paramaters')

wp_0 = 0.5 ## w'/w
wp_0 = 0.5 ## w'/w
d_0 = 0.2  ## d/w distance between the side wires normalized to w
h_0 = 0.5   ## aspect ratio height/w
epsilon = 9 ## permittivity = 9
N_0 = 200   ## discretization points.
x0_0 = 0    ## x0/w
z0 = h_0*1.001    ## z0/w ## start outside the waveguide
z1 = 3    ## z1/w 
nz = 200
# x=0 is in the middle of the center waveguide
# z=0 is at the interface 

list_x0 = [0,0.4,0.6,0.8,1]   ## norm to w
list_z0 = [h_0*0.5, h_0*1.001,h_0*1.5,h_0*2,h_0*3]   ## norm to w
list_z0  = [0.5,0.8,1,1.5]
list_z0  = [ h_0*1.001, h_0*1.5, h_0*2, h_0*3]
z_values = [  h_0*1.001, h_0*1.5, h_0*2, h_0*3]  # example values

list_wp = [0.2,0.4,0.6,0.8]  ## norm to w
list_d = [0.2,1,2,3]  ## norm to w
list_h = [0.2,0.4,0.6,0.8,1]  ## norm to w


# list_N = [200,600,800]

ticks_y = np.arange(0.,1.25,0.25)
ticks_y_label = ['0.00' , '0.25', '0.50' , '0.75','1.00'  ]
ylim_inf1, ylim_sup1 = -0.15, 1.15
xlim_inf1, xlim_sup1 = 0.05,3.15

ticks_y2 = np.arange(-0.5,1.25,0.25)
ticks_y2_label = [r'$-0.5$' , '', '0.0' , '', '0.5' , '','1.0' ]
ylim_inf2, ylim_sup2 = -0.65, 1.15

labelx = r'$z/W$'
labely = r'$\phi(x,z)/V_0$'

# Build a normalization across the full range of x0 values
from matplotlib import cm
from matplotlib.colors import Normalize
norm = Normalize(vmin=0, vmax=3.5)
cmap = cm.hot  # or "plasma", "turbo", etc.

# collect all values from all lists
all_values = np.concatenate([list_x0, list_wp, list_d, list_h, list_z0])
vmax=np.max(all_values)+1
vmin=np.min(all_values)
# define a single normalization based on global min/max
norm = Normalize(vmin=vmin, vmax=vmax)

tolf = 3
vmax1=np.max(list_x0)
vmin1=np.min(list_x0)*tolf
norm1 = Normalize(vmin=vmin1, vmax=vmax1)

vmax2=np.max(list_wp)
vmin2=np.min(list_wp)*tolf
norm2 = Normalize(vmin=vmin2, vmax=vmax2)

vmax3=np.max(list_d)
vmin3=np.min(list_d)*tolf
norm3 = Normalize(vmin=vmin3, vmax=vmax3)

vmax4=np.max(list_h)
vmin4=np.min(list_h)*tolf
norm4 = Normalize(vmin=vmin4, vmax=vmax4)


vmax5=np.max(list_z0)
vmin5=np.min(list_z0)*tolf
norm5 = Normalize(vmin=vmin5, vmax=vmax5)

norm1 = norm
norm2 = norm
norm3 = norm
norm4 = norm
norm5 = norm

#%%
    
print('1-Plot V(z) for different x positions for fixed wp/W = %.2f, d/W = %.2f, h/W = %.2f' %(wp_0,d_0,h_0))

title1 = r'$W_{\text{p}}/W$ = %.2f, $d/W$ = %.2f, $h/W$ = %.2f' %(wp_0,d_0,h_0)

listV_normV0_tot_1 = []
list_z_norm_w_tot1 = []
for x0 in list_x0:
    list_z_norm_w, listV_normU = run_e_out1(wp_0,d_0,h_0,N_0,epsilon,x0,z0,z1,nz)
    
    z_limit = (h_0/2)*1.001
    V0_norm_U = run_e_out1(wp_0,d_0,h_0,N_0,epsilon,0,z_limit,z_limit,2)[1][0]
    
    listV_normV0 = np.array(listV_normU)/V0_norm_U
    listV_normV0_tot_1.append(listV_normV0)
    
    
    list_z_norm_w_tot1.append(list_z_norm_w)

    del listV_normU, listV_normV0

#%%

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
for j, x0 in enumerate(list_x0):
    color = cmap(norm1(x0))  # map x0 to a consistent color
    x0 = list_x0[j]
    plt.plot(np.array(list_z_norm_w_tot1[j]), np.array(listV_normV0_tot_1[j]) ,'.-',color = color ,label = '%.1f' %(x0))
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.yticks(ticks_y2,ticks_y2_label)
plt.ylim(ylim_inf2, ylim_sup2)
os.chdir(path_data)
plt.savefig('Vz_vs_x.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vz_vs_x.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

del x0 

#%%

print('2-Plot V(z) for different wp/W for fixed x/W = %.2f, d/W = %.2f, h/W = %.2f' %(x0_0,d_0,h_0))

title2 = r'$x/W$ = %.2f, $d/W$ = %.2f, $h/W$ = %.2f' %(x0_0,d_0,h_0)

listV_normV0_tot2 = []
list_z_norm_w_tot2 = []
for wp in list_wp:
    list_z_norm_w, listV_normU = run_e_out1(wp,d_0,h_0,N_0,epsilon,x0_0,z0,z1,nz)
    
    
    z_limit = (h_0/2)*1.001
    V0_norm_U = run_e_out1(wp,d_0,h_0,N_0,epsilon,0,z_limit,z_limit,2)[1][0]
    
    listV_normV0 = np.array(listV_normU)/V0_norm_U
    
    listV_normV0_tot2.append(listV_normV0)
    list_z_norm_w_tot2.append(list_z_norm_w)
    
    del listV_normU, listV_normV0


#%%

plt.figure(figsize=tamfig)
plt.title(title2,fontsize=tamtitle)
for j, x0 in enumerate(list_wp):
    color = cmap(norm2(x0))  # map x0 to a consistent color
    wp = list_wp[j]
    plt.plot(np.array(list_z_norm_w_tot2[j]), np.array(listV_normV0_tot2[j]) ,'.-',color = color ,label = r'%.1f' %(wp))
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.yticks(ticks_y2,ticks_y2_label)
plt.ylim(ylim_inf2, ylim_sup2)
plt.xlim(xlim_inf1, xlim_sup1)
os.chdir(path_data)
plt.savefig('Vz_vs_wp.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vz_vs_wp.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

del wp

#%%

print('3-Plot V(z) for different d/W for fixed wp/W = %.2f, x = %.2f, h/W = %.2f' %(wp_0,x0_0,h_0))

title3 = r'$W_{\text{p}}/W$ = %.2f, $x/W$ = %.2f, $h/W$ = %.2f' %(wp_0,x0_0,h_0)

listV_normV0_tot3 = []
list_z_norm_w_tot3 = []
list_slope = []
list_offset = []
listV_normV0_tot3_fit = []

for d in list_d:
    list_z_norm_w, listV_normU = run_e_out1(wp_0,d,h_0,N_0,epsilon,x0_0,z0,z1,nz)
    
    z_limit = (h_0/2)*1.001
    V0_norm_U = run_e_out1(wp_0,d,h_0,N_0,epsilon,0,z_limit,z_limit,2)[1][0]
    
    listV_normV0 = np.array(listV_normU)/V0_norm_U
    
    list_z_norm_w = np.array(list_z_norm_w)
    listV_normV0_tot3.append(listV_normV0)
    list_z_norm_w_tot3.append(list_z_norm_w)
    
    
    ind_max = int(len(list_z_norm_w)/5)
    # ---- linear fitting ----
    popt, _ = curve_fit(linear_func, list_z_norm_w[0:ind_max], listV_normV0[0:ind_max])
    slope, offset = popt
    after_fitting = linear_func(list_z_norm_w, slope, offset)
    
    list_slope.append(slope)
    list_offset.append(offset)
    listV_normV0_tot3_fit.append(after_fitting)
    
    del listV_normU, listV_normV0

#%%

plt.figure(figsize=tamfig)
plt.title(title3,fontsize=tamtitle)
for j, x0 in enumerate(list_d):
    color = cmap(norm3(x0))  # map x0 to a consistent color
    d = list_d[j]
    plt.plot(np.array(list_z_norm_w_tot3[j]), np.array(listV_normV0_tot3[j]) ,'.-',color = color ,label = r'%.1f' %(d))
    plt.plot(np.array(list_z_norm_w_tot3[j]), np.array(listV_normV0_tot3_fit[j]) ,'--',color = color)
    
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.yticks(ticks_y,ticks_y_label)
plt.ylim(ylim_inf1, ylim_sup1)
os.chdir(path_data)
plt.savefig('Vz_vs_d.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vz_vs_d.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

del d

#%%

print('4-Plot V(z) for different h/W for fixed wp/W = %.2f, x = %.2f, d/W = %.2f' %(wp_0,x0_0,d_0))

title4 = r'$W_{\text{p}}/W$ = %.2f, $x/W$ = %.2f, $d/W$ = %.2f' %(wp_0,x0_0,d_0)

listV_normV0_tot4 = []
list_z_norm_w_tot4 = []


list_slope2 = []
list_offset2 = []
listV_normV0_tot4_fit = []

for h in list_h:
    z0 = h*1.001    ## z0/w ## start outside the waveguide
    # print(h,z0)
    list_z_norm_w, listV_normU = run_e_out1(wp_0,d_0,h,N_0,epsilon,x0_0,z0,z1,nz)
    
    z_limit = (h/2)*1.001
    V0_norm_U = run_e_out1(wp_0,d_0,h,N_0,epsilon,0,z_limit,z_limit,2)[1][0]
    
    listV_normV0 = np.array(listV_normU)/V0_norm_U
    
    list_z_norm_w = np.array(list_z_norm_w)
    listV_normV0_tot4.append(listV_normV0)
    list_z_norm_w_tot4.append(list_z_norm_w)
    
    
    ind_max = int(len(list_z_norm_w)/5)
    # ---- linear fitting ----
    popt, _ = curve_fit(linear_func, list_z_norm_w[0:ind_max], listV_normV0[0:ind_max])
    slope2, offset2 = popt
    after_fitting2 = linear_func(list_z_norm_w, slope2, offset2)
    
    
    list_slope2.append(slope2)
    list_offset2.append(offset2)
    listV_normV0_tot4_fit.append(after_fitting2)

#%%

plt.figure(figsize=tamfig)
plt.title(title4,fontsize=tamtitle)
for j, x0 in enumerate(list_h):
    color = cmap(norm4(x0))  # map x0 to a consistent color
    h = list_h[j]
    # print(x0,color)
    plt.plot(np.array(list_z_norm_w_tot4[j]), np.array(listV_normV0_tot4[j]) ,'.-',color = color ,label = r'%.1f' %(h))
    plt.plot(np.array(list_z_norm_w_tot4[j]), np.array(listV_normV0_tot4_fit[j]) ,'--',color = color)
 
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
plt.yticks(ticks_y,ticks_y_label)
plt.ylim(ylim_inf1, ylim_sup1)
plt.xlim(xlim_inf1, xlim_sup1)
os.chdir(path_data)
plt.savefig('Vz_vs_h.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vz_vs_h.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

del h

#%%

# print('5-Plot V(z) for different N for fixed wp/W = %.2f, d/W = %.2f, h/W = %.2f, x/W = %.2f' %(wp_0,d_0,h_0,x0_0))

# title1 = r'$W_{\text{p}}/W$ = %.2f, $d/W$ = %.2f, $h/W$ = %.2f, x/W = %.2f' %(wp_0,d_0,h_0,x0_0)

# listV_normV0_tot_0 = []
# list_z_norm_w_tot0 = []
# for N in list_N:
#     list_z_norm_w, listV_normV0 = run_e_out1(wp_0,d_0,h_0,N,epsilon,x0_0,z0,z1,nz)
#     listV_normV0_tot_0.append(listV_normV0)
#     list_z_norm_w_tot0.append(list_z_norm_w)

#%%

# plt.figure(figsize=tamfig)
# plt.title(title4,fontsize=tamtitle)
# for j,x0 in enumerate(list_N):
#     color = cmap(norm(x0))  # map x0 to a consistent color
#     N = list_N[j]
#     plt.plot(np.array(list_z_norm_w_tot0[j]), np.array(listV_normV0_tot_0[j]) ,'.-',color = color ,label = r'$N$ = %i' %(N))
 
# plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
# plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.xticks(np.arange(h_0,z1+0.5,0.5))
# plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
# plt.yticks(ticks_y,ticks_y_label)
# plt.ylim(ylim_inf1, ylim_sup1)
# os.chdir(path_data)
# plt.savefig('Vz_vs_N.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
# plt.savefig('Vz_vs_N.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
# plt.show()

# del N 

#%%

cmap2 = plt.cm.RdBu   # define the colormap

print('6-Plot V(x,z) for fixed wp/W = %.2f, d/W = %.2f, h/W = %.2f' %(wp_0,d_0,h_0))

xleft = -1/2 - d_0 - wp_0/2 ## position of left waveguide 
xright = - xleft
x0 = np.round(xleft*2)
x1 = np.round(xright*2)   ## x1/w
nx = nz 
z0 = h_0*1.001    ## z0/w ## start outside the waveguide
z1 = 3    ## z1/w 
list_x_norm_w_tot5, list_z_norm_w_tot5, listV_normV0_tot5 = run_e_out2(wp_0,d_0,h_0,N_0,epsilon,x0,x1,nx,z0,z1,nz)

z_limit = (h_0/2)*1.001
V0_norm_U = run_e_out1(wp_0,d_0,h_0,N_0,epsilon,0,z_limit,z_limit,2)[1][0]
listV_normV0_tot5_norm = np.array(listV_normV0_tot5)/V0_norm_U

#%%

labelx = r'$x/W$'
labely = r'$z/W$'
labelz = r'$\phi(x,z)/V_0$'

title5 = r'$W_{\text{p}}/W$ = %.2f, $d/W$ = %.2f, $h/W$ = %.2f' %(wp_0,d_0,h_0)

limits1 = [np.min(list_x_norm_w_tot5) , np.max(list_x_norm_w_tot5),np.min(list_z_norm_w_tot5) , np.max(list_z_norm_w_tot5)]
# Create triangulation from scattered (x, z) points
triang = tri.Triangulation(list_x_norm_w_tot5, list_z_norm_w_tot5)
ejey_z = np.linspace(z0,z1,10)

vmin1 , vmax1 = np.nanmin(listV_normV0_tot5_norm), np.nanmax(listV_normV0_tot5_norm)
vabs = max(abs(vmin1), abs(vmax1))
# bounds1 =   np.linspace(vmin1 , vmax1 , 10) 

# Set tighter bounds near zero
bounds1 = np.concatenate([
    np.linspace(-vabs, -0.1, 5, endpoint=True),
    [-0.05, -0.01, 0, 0.01, 0.05],
    np.linspace(0.1, vabs, 5, endpoint=True)
])


# bounds1 =   np.logspace(np.log10(vmin1), np.log10(100) , 10) 

norm1 = mpl.colors.BoundaryNorm(bounds1, cmap2.N)

fig = plt.figure(figsize=tamfig)

ax = fig.add_axes([0.15, 0.15, 0.7, 0.7])  # [left, bottom, width, height]
cax = fig.add_axes([0.15, 0.88, 0.7, 0.05])  # [left, bottom, width, height]

# plt.title(title5,fontsize=tamtitle-4)
tpc = ax.tripcolor(triang, listV_normV0_tot5_norm, shading='flat', cmap=cmap2, norm=norm1)
cbar = fig.colorbar(tpc,cax=cax, orientation="horizontal", format = '%.1f') 
cbar.ax.set_title(labelz,fontsize=tamletra-1)
cbar.ax.tick_params(labelsize = tamnum-1, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar.ax.xaxis.set_ticks_position("top")
cbar.ax.xaxis.set_label_position("top")

# Move y-axis labels and ticks to the right
ax.yaxis.tick_right()         # move tick labels to right
ax.yaxis.set_label_position("right")

# But keep tick marks on the left
ax.tick_params(axis="y", left=True, right=True)

ax.set_xlabel(labelx,fontsize=tamletra-1,labelpad =labelpadx)
ax.set_ylabel(labely,fontsize=tamletra-1,labelpad =labelpady)
ax.tick_params(labelsize = tamnum-1, length = 3 , width=1, direction="in",which = 'both', pad = pad)
 

# plt.plot([xleft],[h_0],'o',color =  'purple')
# plt.plot([xright],[h_0],'o',color =  'purple')
ax.set_xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[r"$-2$","",r"$-1$","","0","","1","","2"])
# plt.yticks(ticks_y,ticks_y_label)

ax.set_xlim(np.min(list_x_norm_w_tot5) , np.max(list_x_norm_w_tot5))

ax.set_ylim(np.min(list_z_norm_w_tot5) , np.max(list_z_norm_w_tot5))
# plt.plot(np.ones(10)*(-1/2),ejey_z,'--',color =  'black')
# plt.plot(np.ones(10)*(1/2),ejey_z,'--',color =  'black')
os.chdir(path_data)
plt.savefig('Vxz.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vxz.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=200)  
plt.show()

#%%

print('7-Plot V(z) for different z positions for fixed wp/W = %.2f, d/W = %.2f, h/W = %.2f' %(wp_0,d_0,h_0))

title1 = r'$W_{\text{p}}/W$ = %.2f, $d/W$ = %.2f, $h/W$ = %.2f' %(wp_0,d_0,h_0)

listV_normV0_tot_6 = []
list_x_norm_w_tot6 = []      
for z0 in list_z0:
    list_x_norm_w, listV_normV0 = run_e_out0(wp_0,d_0,h_0,N_0,epsilon,z0,x0,x1,nx)
    
    
    z_limit = (h_0/2)*1.001
    V0_norm_U = run_e_out1(wp_0,d_0,h_0,N_0,epsilon,x0_0,z_limit,z_limit,2)[1][0]
    
    
    listV_normV0_tot_6.append(np.array(listV_normV0)/V0_norm_U)
    list_x_norm_w_tot6.append(list_x_norm_w)

#%%

plt.figure(figsize=tamfigure2)
plt.title(title1,fontsize=tamtitle)
for j,x0 in enumerate(list_z0):
    color = cmap(norm5(x0))  # map x0 to a consistent color
    # print(x0,color)
    z0 = list_z0[j]
    plt.plot(np.array(list_x_norm_w_tot6[j]), np.array(listV_normV0_tot_6[j]) ,'.-',color = color ,label = r'%.1f' %(z0))
 
plt.xlabel(r'$x/W$',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# plt.xticks(np.arange(h_0,z1+0.5,0.5))
plt.xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[r"$-2$","",r"$-1$","","0","","1","","2"])
plt.yticks(ticks_y2,ticks_y2_label)
plt.ylim(ylim_inf2, ylim_sup2)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
os.chdir(path_data)
plt.xlim(-2,2)
plt.savefig('Vx_vs_z.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('Vx_vs_z.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.show()

del z0

#%%

print('8-Plot the derivative of V(x,z) respect to x for fixed wp/W = %.2f, d/W = %.2f, h/W = %.2f' %(wp_0,d_0,h_0))

from scipy.interpolate import CloughTocher2DInterpolator

# Ensure inputs are numpy arrays
x = np.array(list_x_norm_w_tot5)
z = np.array(list_z_norm_w_tot5)
v = np.array(listV_normV0_tot5)

# Build interpolator
interp_V = LinearNDInterpolator(list(zip(x, z)), v)
interp_V2 = CloughTocher2DInterpolator(list(zip(x, z)), v)
# Choose evaluation points (can be original points or a regular grid)
dx = 1e-8  # small step for finite difference

# Derivative w.r.t. x using central difference
dV_dx = (interp_V(x + dx, z) - interp_V(x - dx, z)) / (2 * dx)
dV_dx_2 = (interp_V2(x + dx, z) - interp_V2(x - dx, z)) / (2 * dx)

# Create mask to exclude NaN values
valid_mask = ~np.isnan(dV_dx)
valid_mask2 = ~np.isnan(dV_dx_2)

# Filter arrays for plotting and analysis
x_valid = x[valid_mask]
z_valid = z[valid_mask]
dV_dx_valid = dV_dx[valid_mask]
triang2 = tri.Triangulation(x_valid, z_valid)

x_valid_2 = x[valid_mask2]
z_valid_2 = z[valid_mask2]
dV_dx_valid_2 = dV_dx_2[valid_mask2]
triang2_2 = tri.Triangulation(x_valid_2, z_valid_2)

# Build interpolator of derivative 
interp_dV_dx = LinearNDInterpolator(list(zip(x_valid, z_valid)), dV_dx_valid)
interp_dV_dx_2 = LinearNDInterpolator(list(zip(x_valid_2, z_valid_2)), dV_dx_valid_2)
 
 #%%

labelz = r'$\partial_x \phi/U$'

limits1 = [np.min(list_x_norm_w_tot5) , np.max(list_x_norm_w_tot5),np.min(list_z_norm_w_tot5) , np.max(list_z_norm_w_tot5)]
# Create triangulation from scattered (x, z) points
 
vmin1 , vmax1 = np.nanmin(dV_dx_valid_2), np.nanmax(dV_dx_valid_2)
# bounds1 =   np.linspace(vmin1 , vmax1 , 10) 
# Force zero to be one of the boundaries
vabs = max(abs(vmin1), abs(vmax1))
# bounds1 =   np.linspace(vmin1 , vmax1 , 10) 

# Set tighter bounds near zero
bounds1 = np.concatenate([
    np.linspace(-vabs, -0.1, 5, endpoint=True),
    [-0.05, -0.01, 0, 0.01, 0.05],
    np.linspace(0.1, vabs, 5, endpoint=True)
])


# bounds1 =   np.logspace(np.log10(vmin1), np.log10(100) , 10) 
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap2.N)

plt.figure(figsize=tamfig2)
# plt.title(title5,fontsize=tamtitle-4)
tpc2 = plt.tripcolor(triang2, dV_dx_valid, shading='flat', cmap=cmap2, norm=norm1)

cbar = plt.colorbar(tpc2, fraction=0.046, pad=0.04 , format = '%.2f') 
cbar.ax.set_title(labelz,fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 3 , width=1, direction="in",which = 'both', pad = pad)
# plt.plot([xleft],[h_0],'o',color =  'purple')
# plt.plot([xright],[h_0],'o',color =  'purple')
plt.xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[r"$-2$","",r"$-1$","","0","","1","","2"])
plt.xlim(-2,2)
plt.xlim(np.min(list_x_norm_w_tot5) , np.max(list_x_norm_w_tot5))
plt.ylim(np.min(list_z_norm_w_tot5) , np.max(list_z_norm_w_tot5))
# plt.plot(np.ones(10)*(-1/2),ejey_z,'--',color =  'black')
# plt.plot(np.ones(10)*(1/2),ejey_z,'--',color =  'black')
os.chdir(path_data)
plt.savefig('dVxz_dx.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.savefig('dVxz_dx.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.09, dpi=200)  
plt.show()

 #%%
 

 
# Define range of x-values where to evaluate the cut
x_cut = x
plt.figure(figsize=tamfigure2)
plt.title(title5,fontsize=tamtitle)
# Create constant z array

# z_values = [0.5,0.8,1,1.5]  # example values

for j,x0 in enumerate(z_values):  
    zc = z_values[j]
    color = cmap(norm5(x0))  # map x0 to a consistent color
    z_cut = np.full_like(x_cut, zc)
    dV_dx_cut = (interp_V(x_cut + dx, z_cut) - interp_V(x_cut - dx, z_cut)) / (2 * dx)
    plt.plot(x_cut, dV_dx_cut, '-',color = color, label=fr'$z/W$ = {zc:.1f}')

    
plt.xlabel(r'$x/W$',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\partial_x \phi/U$',fontsize=tamletra,labelpad =labelpady)
plt.xticks([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[r"$-2$","",r"$-1$","","0","","1","","2"])
# plt.yticks(ticks_y,ticks_y_label)
plt.xlim(-2,2)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc ='best',markerscale=2,fontsize=tamletra,frameon=0,handletextpad=0.1, handlelength=1.3,labelspacing = 0.2) 
os.chdir(path_data)
plt.savefig('dVxz_dx_cut_along_x.png', format='png',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
plt.savefig('dVxz_dx_cut_along_x.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.02, dpi=dpi)  
# plt.grid(True)
plt.show()

