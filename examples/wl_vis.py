import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
np.set_printoptions(suppress=True)

colors = {
    "soft_blue": "#4F81BD",
    "muted_orange": "#F7A800",
    "soft_red": "#D14F5D",
    "gentle_green": "#6DAF69",
    "soft_purple": "#9B59B6",
    "dusty_pink": "#E17D85",
    "slate_blue": "#6A7F99",
    "warm_tan": "#D3B69B",
    "earthy_brown": "#8E735B",
    "moss_green": "#6A8A3B",
    "rich_gold": "#B88A2A",
    "charcoal_gray": "#4A4A48"
}

colors = {
    "steel_blue": "#1F77B4",
    "light_steel_blue": "#AEC7E8",
    "orange": "#FF7F0E",
    "light_orange": "#FFBB78",
    "forest_green": "#2CA02C",
    "light_green": "#98DF8A",
    "firebrick_red": "#D62728",
    "soft_red": "#FF9896",
    "lavender": "#9467BD",
    "light_lavender": "#C5B0D5",
    "brown": "#8C564B",
    "tan": "#C49C94",
    "orchid": "#E377C2",
    "light_orchid": "#F7B6D2",
    "gray": "#7F7F7F",
    "light_gray": "#C7C7C7",
    "yellow_green": "#BCBD22",
    "light_yellow_green": "#DBDB8D",
    "turquoise": "#17BECF",
    "light_turquoise": "#9EDAE5"
}

indexed_colors_values = {index: value for index, (key, value) in enumerate(colors.items())}
# Convert hex to RGB, then to HSV
rgb_colors = [mcolors.hex2color(color) for color in colors.values()]
hsv_colors = [mcolors.rgb_to_hsv(rgb) for rgb in rgb_colors]
# Sort colors based on the Hue value (hsv_colors[i][0])
sorted_colors_by_hue = [indexed_colors_values[i] for i in np.argsort([hsv[0] for hsv in hsv_colors])]
custom_cmap = ListedColormap(sorted_colors_by_hue)
custom_cmap = ListedColormap(colors.values())
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=custom_cmap.colors)

subfolders = [ f.name for f in os.scandir(os.getcwd()) if f.is_dir() ]
print("Available directories:")
print(subfolders)
directory = input("Input directory to pull data from: ")

filename = "{}/wl_dos_bins.dat".format(directory)
bin_edges = nc.Dataset(filename)
bin_edges = np.array(bin_edges["grid data"][:], dtype=np.float64)

filename = "{}/wl_dos.dat".format(directory)
wl_logdos = nc.Dataset(filename)
wl_logdos = np.array(wl_logdos["grid data"][:], dtype=np.float64)
#print(wl_logdos)

filename = "{}/radial_densities/rho_of_E.dat".format(directory)  
rho_of_E = nc.Dataset(filename)

rho = rho_of_E.variables['rho data'][:]
U_data = rho_of_E.variables['U data'][:]

elements = ['Al', 'Ti', 'Cr', 'Mo']
n_species = 4
concentrations = 1.0/n_species*np.ones(n_species)

wl_logdos = wl_logdos - np.max(wl_logdos)
wl_logdos = np.log(np.exp(wl_logdos)/np.sum(np.exp(wl_logdos)))

filename = "{}/wl_hist.dat".format(directory)
wl_hist = nc.Dataset(filename)
wl_hist = np.array(wl_hist["grid data"][:], dtype=np.float64)

kb_ev = 8.167333262e-5
ev_to_ry = 13.605693122
kb_ry = kb_ev/ev_to_ry
n_atoms = 432
ev_to_mev = 1000
J = 0.001
N = 6

# unit conversion
bin_edges = bin_edges*ev_to_ry # converts from Ryd to eV

orig_temp = 3000
start_temp = 400
end_temp = 3000
step_size = (end_temp-start_temp)/(75)
temperatures = np.arange(start_temp, end_temp, step_size)

plt.plot(np.exp(wl_logdos))
#plt.show()
plt.close()
# ------------------------------------------------

bin_width = bin_edges[1] - bin_edges[0]

prob = np.zeros(len(bin_edges)-1)
beta = 1.0/(kb_ry*orig_temp)
for ibin, edge in enumerate(bin_edges[:-1]):
    bin_energy = edge + 0.5*bin_width
    prob[ibin] = np.exp(wl_logdos[ibin]-beta*bin_energy)
plt.plot(prob)
#plt.show()
plt.close()

wcs = np.zeros((len(U_data),n_species,n_species,2))
for i in range(len(U_data)):
  wcoft = np.zeros((n_species,n_species,2))
  for j in range(len(elements)):
      for k in range(len(elements)):
          wcoft[j,k,0] = 1.0-1.0/8.0*rho[i,1,j,k]/concentrations[j]
          wcoft[j,k,1] = 1.0-1.0/6.0*rho[i,2,j,k]/concentrations[j]
  wcs[i] = wcoft

elements = ""  # Initialize an empty string to store species data
with open("{}/input.txt".format(directory), 'r') as file:
    for line in file:
        if 'species_name' in line:  # Check if the line contains the species_name variable
            parts = line.split("=")  # Split the line by '=' and get everything after it
            if len(parts) > 1:
                elements = parts[1].strip().split()  # Get the data after '=' and strip spaces
            break  # Exit the loop after finding species_name

pairs = np.zeros([int(len(elements)*(len(elements)+1)/2),2], dtype=np.int16)
k = 0
for i in range(len(elements)):
    for j in range(i, len(elements)):
        pairs[k] = [i, j]
        k += 1

# Setup plots
#fig, [ax1, ax2, ax3, ax4] = plt.subplots(1, 4, figsize=(24, 7), constrained_layout=True)
#fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(24, 9), constrained_layout=True)

fig_height = 5; fig_width=fig_height*1.68
fig1, ax1 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
fig2, ax2 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
fig3, ax3 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
ax4 = ax3.twinx()
fig4, ax5 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)

ax1.set_xlabel(r'Energy $U$ (meV/atom)')
ax1.set_ylabel(r'Probability P($U$)')
ax1.set_title(r'Energy Histograms')

ax2.set_xlabel(r'Temperature (K)')
ax2.set_ylabel(r'$\langle U \rangle$ (meV/atom)')
ax2.set_title(r'Mean Energy $\langle U \rangle$')  

ax3.set_xlabel(r'Temperature (K)')
ax3.set_ylabel(r'$C_v$ ($k_B$/atom)')
ax3.set_title(r'Heat Capacity $C_v$') 

ax4.set_ylabel(r'$\alpha^{pq}_1$')

ax5.set_xlabel(r'Temperature (K)')
ax5.set_ylabel(r'Entropy ($k_B$/atom)')
ax5.set_title(r'Entropy $S$')

# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
gibbs_energies  = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))
asr_orders      = np.zeros((len(pairs),len(temperatures)))

# Loop over temperatures of interest
for itemp, new_temp in enumerate(temperatures):

    beta = 1.0/(kb_ev*new_temp)

    # Reweighted histogram
    prob = np.zeros(len(bin_edges)-1)
    
    # Reweight histogram
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        prob[ibin] = wl_logdos[ibin]-beta*bin_energy
   
    prob = prob-np.max(prob)/2
    prob = np.exp(prob)

    # Normalise
    prob = prob/(np.sum(prob))
    
    # Only plot every 5th histogram to avoid crowding the axes
    if itemp%5 == 0:
        strlabel = "T={:8.2f} K".format(new_temp)
        ax1.stairs(prob, bin_edges/n_atoms*ev_to_mev, label=strlabel, fill=True)
    
    # Mean energy
    mean_energy = np.dot(bin_edges[:-1]+0.5*bin_width, prob)
    mean_energies[itemp] = mean_energy

    for ipair, pair in enumerate(pairs):
      asr_order = np.dot(wcs[:,pair[0],pair[1],0], prob)
      asr_orders[ipair,itemp] = asr_order

    # Compute heat capacity using the histogram
    msq_dev = np.zeros(len(bin_edges)-1)
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
        
    heat_caps[itemp] = np.dot(msq_dev, prob)*bin_width/(kb_ev*new_temp**2)

gibbs_energies[0] = mean_energies[0]
for itemp in range(1,len(temperatures)):

  beta_i = 1.0/(kb_ev*temperatures[itemp])
  beta_j = 1.0/(kb_ev*temperatures[itemp-1])

  gibbs_energies[itemp] = beta_j*gibbs_energies[itemp-1]/beta_i+((mean_energies[itemp-1]+mean_energies[itemp])*(beta_i-beta_j))/(2*beta_i)

for itemp, new_temp in enumerate(temperatures):
  entropies[itemp] = (mean_energies[itemp]-gibbs_energies[itemp])/new_temp
# Complete plots using data computed above
local_max_indices = np.where((heat_caps[1:-1] > heat_caps[:-2]) & (heat_caps[1:-1] > heat_caps[2:]))[0] + 1
local_min_indices = np.where((heat_caps[1:-1] < heat_caps[:-2]) & (heat_caps[1:-1] < heat_caps[2:]))[0] + 1

# unit conversion
mean_energies = mean_energies/n_atoms*ev_to_mev
heat_caps = heat_caps/kb_ev/n_atoms
entropies = entropies/kb_ev/n_atoms

ax2.plot(temperatures, mean_energies, '-o', markersize=4)
ax3.plot(temperatures, heat_caps, '-o', markersize=4)
for ipair, pair in enumerate(pairs):
  ax4.plot(temperatures, asr_orders[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
ax5.plot(temperatures, entropies, '-o', markersize=4)

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax4.legend(loc='center left', bbox_to_anchor=(1.15,0.5))

ax3.tick_params(direction="in")
ax4.tick_params(direction="in")

ax2_diff = np.max(mean_energies) - np.min(mean_energies)
ax3_diff = np.max(heat_caps) - np.min(heat_caps)
ax5_diff = np.max(entropies) - np.min(entropies)

for x in local_max_indices:
  ax2.vlines(x=temperatures[x], ymin=np.min(mean_energies)-0.1*ax2_diff, ymax=np.max(mean_energies), color=colors["firebrick_red"], linestyle='--')
  ax3.vlines(x=temperatures[x], ymin=np.min(heat_caps)-0.1*ax3_diff, ymax=np.max(heat_caps), color=colors["firebrick_red"], linestyle='--')
  ax5.vlines(x=temperatures[x], ymin=np.min(entropies)-0.1*ax5_diff, ymax=np.max(entropies), color=colors["firebrick_red"], linestyle='--')

ax2.set_ylim(np.min(mean_energies)-0.025*np.abs(ax2_diff), np.max(mean_energies)+0.025*np.abs(ax2_diff))
ax3.set_ylim(np.min(heat_caps)-0.025*np.abs(ax3_diff), np.max(heat_caps)+0.025*np.abs(ax3_diff))
ax5.set_ylim(np.min(entropies)-0.025*np.abs(ax5_diff), np.max(entropies)+0.025*np.abs(ax5_diff))

fig1.savefig('figures/{}_energy_histogram.pdf'.format(''.join(elements)), bbox_inches='tight')
fig2.savefig('figures/{}_mean_energy.pdf'.format(''.join(elements)), bbox_inches='tight')
fig3.savefig('figures/{}_heat_capacity.pdf'.format(''.join(elements)), bbox_inches='tight')
fig4.savefig('figures/{}_entropy.pdf'.format(''.join(elements)), bbox_inches='tight')
