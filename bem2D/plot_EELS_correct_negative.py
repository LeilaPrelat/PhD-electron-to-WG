
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the EELS of a rectangular waveguide
on top of substrate
with epsilon= 12.25 + i*0.2
for different L and N to see
when the EELS numerically converges
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

 
path_basic = os.getcwd()
os.chdir(path_basic)

from global_constants import constants 
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

path_data = os.path.join(path_basic, 'convergency_negative_EELS')

list_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

#%%

tamfig = [4.7,3.5]
tamfig2 = [4.7,3.5]
tamfig3 = [4.1,3.5]
tamletra = 14
tamtitle  = 6
tamnum = tamletra
tamlegend = tamletra
labelpady = 2
labelpadx = 3
pad = 3
mk = 1
ms = 2.5
hp = 0.3
length_marker = 0
dpi = 500
lw = 0.8

deltax,deltay = 3,8

values = tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,deltax,deltay

# //     N   = the total number of parametrization points is N+50
# //     w   = total width of wg     || x
# //     h   = total thickness of wg  || y
# //     ze  = impact parameter is h/2+ze. ze distance between electron and wg surface  

W=500 # nm
h=200 ## nm

Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
me_c2_eV = 510998.95069  ## me*c**2 in eV
beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
ze = 50 

labelx=r'$L$ (nm)'
labelx2 = r"Number of discret. points, $N_{\text{t}}$ ($10^3$)"
labely=r'd$\Gamma$/d$y$ (1/eV$\cdot$nm)'
energy = 0.35
title = r'$W$ = %i nm, $h$ = %i nm, $\hbar\omega$ = %.2f eV, $z_{\rm e}$ = %i nm, $E_{\rm e}$ = %i keV' %(W,h,energy,ze,Ee_electron_keV)
os.chdir(path_data)

#%%

def total_number_points(L,N):
    R = L + W/2
    tot_perimeter = N/W*(2*W+2*L+np.pi*R+2*h)
    return tot_perimeter*1e-3 ## units of 10^3

#%%

print('Plot EELS vs L,N')

list_LN = [ [700, 500], [1000, 800], [1000, 1000], [1500, 800], [1500, 1000], [1000, 1500],
          [1500, 1500], [2000,1000], [3000,1000], [1500,5000]]

# --- assign colors by N ---
N_values = [N for _, N in list_LN]
unique_N = sorted(set(N_values))   # unique N sorted
cmap = cm.get_cmap("viridis", len(unique_N))  # gradient colormap
color_map = {N: cmap(i) for i, N in enumerate(unique_N)}  # dict N→color

# Track which N has already been labeled
plotted_labels = set()

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
L_values = []
list_EELS = []
for pairs_LN in list_LN:
    L, N = pairs_LN
 
    name3 = 'EELS_comsol2_N%i_W%inm_h%inm_L%inm_Ee%ikeV_energy%.2f.dat' %(N,W,h,L,Ee_electron_keV,energy) ## virtual geometry
 
    try:
        tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
        tabla1_2 = np.transpose(tabla1)
    except (ValueError, FileNotFoundError):
        tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
        tabla1_2 = np.transpose(tabla1)
      
    eV,x,y,z,EELS = tabla1_2
 
  # Add label only the first time this N appears
    if N not in plotted_labels:
        plt.plot(L, EELS, "o", color=color_map[N], label=f"N = {N}")
        plotted_labels.add(N)
    else:
        plt.plot(L, EELS, "o", color=color_map[N])

    L_values.append(L)
    list_EELS.append(EELS)

# deduplicate while preserving order
unique_L = []
for L in L_values:
    if L not in unique_L:
        unique_L.append(L)

# make labels from the unique values
L_labels = [f"{L:.0f}" for L in unique_L]

plt.xticks(unique_L,L_labels)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)

plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.yscale('log')
plt.ylim(np.min(list_EELS)-1e-5,1e-5)
plt.legend(loc = 'best',markerscale=1,fontsize=tamnum,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_negative_convergency_vs_L_energy%.2feV.png'%(energy),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_negative_convergency_vs_L_energy%.2feV.pdf'%(energy),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()

#%%

from sklearn.linear_model import LinearRegression

print('Plot EELS vs number of points + linear fit')

 # --- assign colors by N ---
L_values = [L for L,_ in list_LN]
unique_L = sorted(set(L_values))   # unique N sorted
cmap = cm.get_cmap("viridis", len(unique_L))  # gradient colormap
color_map = {L: cmap(i) for i, L in enumerate(unique_L)}  # dict N→color
# Track which N has already been labeled
plotted_labels = set()

plt.figure(figsize=tamfig2)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx2,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)

list_EELS = []
list_tot_numbers = []
for pairs_LN in list_LN:
    L, N = pairs_LN
 
    name3 = 'EELS_comsol2_N%i_W%inm_h%inm_L%inm_Ee%ikeV_energy%.2f.dat' %(N,W,h,L,Ee_electron_keV,energy) ## virtual geometry
 
    try:
        tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
        tabla1_2 = np.transpose(tabla1)
    except (ValueError, FileNotFoundError):
        tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
        tabla1_2 = np.transpose(tabla1)
      
    eV,x,y,z,EELS = tabla1_2
 
    tot_num_points = total_number_points(L,N)
    list_EELS.append(EELS)
    list_tot_numbers.append(tot_num_points)
    # Add label only the first time this N appears
    if L not in plotted_labels:
        plt.plot(tot_num_points, EELS, "o", color=color_map[L], label="L = %.1f $\mu$m"%(L*1e-3))
        plotted_labels.add(L)
    else:
        plt.plot(tot_num_points, EELS, "o", color=color_map[L])


# Create and fit the model
model1 = LinearRegression()
x_for_fit = np.array(list_tot_numbers).reshape(-1, 1) 

print("Fit the EELS with log(Ntot), then EELS = A*log(Ntot) + B")
model1.fit(np.log(x_for_fit), list_EELS)
# Extract slope and intercept
slope1 = model1.coef_[0]
intercept1 = model1.intercept_
print("Exponential fit of EELS vs N")
print(f"Slope: {slope1:.4f}, Intercept: {intercept1:.4f}")

x_max1 = np.exp(-intercept1/slope1)
x_fit1 = np.linspace(np.min(list_tot_numbers), x_max1, 100).reshape(-1, 1)
y_fit1 = []
for x in x_fit1: 
    y_fit1.append(slope1*np.log(x) + intercept1) ## t = e^B*Ntot^A

print("Fit the EELS with Ntot, then EELS = A*Ntot + B")
model2 = LinearRegression()
model2.fit(x_for_fit, list_EELS)
# Extract slope and intercept
slope2 = model2.coef_[0]
intercept2 = model2.intercept_
print("Exponential fit of EELS vs N")
print(f"Slope: {slope2:.4f}, Intercept: {intercept2:.4f}")
y_fit2 = model2.predict(x_fit1)
x_max2 = -intercept2/slope2

plt.plot(x_fit1,y_fit1,'r--', label="Exp. fit")
plt.plot(x_fit1,y_fit2,'r-', label="Linear fit")
plt.plot([x_max1],[0],'x',color = 'black')
plt.plot([x_max2],[0],'x',color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = [0.54,-0.02],markerscale=1,fontsize=tamlegend,frameon=0,labelspacing =0.2,handletextpad=0.2, handlelength=1)
# plt.xscale('log')
plt.xticks(np.arange(10,70,10))
plt.ylim(np.min(list_EELS)-1e-5,1e-5)
plt.tight_layout()
plt.savefig('EELS_negative_vs_tot_points_energy%.2feV.png'%(energy),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_negative_vs_tot_points_energy%.2feV.pdf'%(energy),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()

 #%%
 
print("Import the time of running and plot EELS vs time")
 
time_label = "Time(min)"

# Read file
df = pd.read_csv("results.txt", sep="\t")
# Convert 'Time(min)' column to numeric, turn errors into NaN
df[time_label] = pd.to_numeric(df[time_label], errors="coerce")
# Drop rows where time is missing (NaN)
df_clean = df.dropna(subset=[time_label]).copy()
# Convert to integer if desired
df_clean[time_label] = df_clean[time_label].astype(int)
print("Cleaned data:")
print(df_clean)
# Access as numpy arrays if needed
list_time = df_clean[time_label].to_numpy() ##in min
list_time_in_hr = np.array(list_time)/60
 
# Access as numpy arrays if needed
L_from_data = df_clean["L(nm)"].to_numpy()
N_from_data = df_clean["N"].to_numpy()
list_tot_numbers_from_data = []
for j in range(len(N_from_data)):
    L = L_from_data[j]
    N = N_from_data[j]
    tot_num_points = total_number_points(L,N)
    list_tot_numbers_from_data.append(tot_num_points)


# Create and fit the model
model = LinearRegression()
x_for_fit = np.array(list_tot_numbers_from_data).reshape(-1, 1)  
print("Fit the log(time) with log(Ntot), then log(t) = A*log(Ntot) + B --> t = e^B*Ntot^A")
model.fit(np.log(x_for_fit), np.log(list_time_in_hr))
# Extract slope and intercept
slope = model.coef_[0]
intercept = model.intercept_
print(f"Slope: {slope:.4f}, Intercept: {intercept:.4f}")
y_fit = []
for x in x_fit1: 
    y_fit.append(np.exp(intercept)*x**(slope)) ## t = e^B*Ntot^A
    #y_fit.append(slope*np.log(x) + intercept)

y_fit_max1 = np.exp(intercept)*x_max1**(slope)
y_fit_max2 = np.exp(intercept)*x_max2**(slope)

plt.figure(figsize=tamfig3)
plt.title(title,fontsize=tamtitle)
plt.ylabel("Time (hr)",fontsize=tamletra,labelpad =labelpadx)
plt.xlabel(labelx2,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_tot_numbers_from_data,list_time_in_hr,"o")
plt.plot(x_fit1,y_fit,'r-', label="Linear fit")
plt.plot([x_max1],[y_fit_max1],'x',color = 'black')
plt.plot([x_max2],[y_fit_max2],'x',color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)

plt.xscale('log')
plt.yscale("log")
# Define custom tick positions
ticks = np.arange(10, 70, 10)
plt.xticks(ticks, labels=[str(t) for t in ticks])  # force labels
#plt.xticks([2*1e4,4*1e4,6*1e4], ["20000","40000","60000"])
plt.legend(loc = "best",markerscale=1,fontsize=tamlegend,frameon=0,labelspacing =0.2,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_negative_vs_time_energy%.2feV.png'%(energy),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_negative_vs_time_energy%.2feV.pdf'%(energy),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()