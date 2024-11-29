import numpy as np
import matplotlib.pyplot as plt
import math as m
import netCDF4 as nc
import os
np.set_printoptions(suppress=True)

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
print(wl_logdos)

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
n_atoms = 128
eV_to_meV = 1000
J = 0.001
N = 6

orig_temp = 2000
start_temp = 400
end_temp = 3000
step_size = (end_temp-start_temp)/(75)
temperatures = np.arange(start_temp, end_temp, step_size)

plt.plot(np.exp(wl_logdos))
plt.show()
# ------------------------------------------------

bin_width = bin_edges[1] - bin_edges[0]

prob = np.zeros(len(bin_edges)-1)
beta = 1.0/(kb_ry*orig_temp)
for ibin, edge in enumerate(bin_edges[:-1]):
    bin_energy = edge + 0.5*bin_width
    prob[ibin] = np.exp(wl_logdos[ibin]-beta*bin_energy)
plt.plot(prob)
plt.show()

wcs = np.zeros((len(U_data),n_species,n_species,2))
for i in range(len(U_data)):
  wcoft = np.zeros((n_species,n_species,2))
  for j in range(len(elements)):
      for k in range(len(elements)):
          wcoft[j,k,0] = 1.0-1.0/8.0*rho[i,1,j,k]/concentrations[j]
          wcoft[j,k,1] = 1.0-1.0/6.0*rho[i,2,j,k]/concentrations[j]
  wcs[i] = wcoft

labels = ['Al', 'Ti', 'Cr', 'Mo']

pairs = [[0,0], [1,1], [2,2], [3,3], [0,1], [0,2], [0,3], [1,2], [1,3], [2,3]]

# Setup plots
#fig, [ax1, ax2, ax3, ax4] = plt.subplots(1, 4, figsize=(24, 7), constrained_layout=True)
fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(24, 9), constrained_layout=True)

ax1.set_xlabel(r'Energy $(R_y)$')
ax1.set_ylabel(r'P(U)')
ax1.set_title(r'Energy Histograms')

ax2.set_xlabel(r'Temperature (K)')
ax2.set_ylabel(r'$\langle U \rangle$ $(R_y)$')
ax2.set_title(r'Mean Energy $\langle U \rangle$')  

ax3.set_xlabel(r'Temperature (K)')
ax3.set_ylabel(r'$C_v$ $(R_yK^{-1})$')
ax3.set_title(r'Heat Capacity') 

#ax4.set_xlabel('temperature T (K)')
#ax4.set_ylabel('S (kB^-1)')
#ax4.set_title('Entropy vs temperature')

ax5 = ax3.twinx()
ax5.set_ylabel(r'$\alpha^{pq}_1$')  # we already handled the x-label with ax1

# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))
asr_orders      = np.zeros((len(pairs),len(temperatures)))

# Loop over temperatures of interest
for itemp, new_temp in enumerate(temperatures):

    beta = 1.0/(kb_ry*new_temp)

    # Reweighted histogram
    prob = np.zeros(len(bin_edges)-1)
    
    # Reweight histogram
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        prob[ibin] = np.exp(wl_logdos[ibin]-beta*bin_energy)  

    # Normalise
    prob = prob/(np.sum(prob))

    # Only plot every 5th histogram to avoid crowding the axes
    if itemp%5 == 0:
        strlabel = "T={:8.2f}".format(new_temp)
        ax1.bar(bin_edges[:-1], prob, width=bin_width, align='edge', label=strlabel)
    
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
        
    heat_caps[itemp] = np.dot(msq_dev, prob)*bin_width/(kb_ry*new_temp**2)
    heat_caps[itemp] = heat_caps[itemp]

    entropy = prob*np.log(prob)
    entropy[np.isnan(entropy)] = 0
    entropies[itemp] = -np.sum(entropy)

# Complete plots using data computed above
local_max_indices = np.where((heat_caps[1:-1] > heat_caps[:-2]) & (heat_caps[1:-1] > heat_caps[2:]))[0] + 1
for x in local_max_indices:
  ax2.vlines(x=temperatures[x], ymin=np.min(mean_energies)-0.2*np.abs(np.min(mean_energies)), ymax=np.max(mean_energies), color='r', linestyle='--')
  ax3.vlines(x=temperatures[x], ymin=np.min(heat_caps)-0.2*np.abs(np.min(heat_caps)), ymax=np.max(heat_caps), color='r', linestyle='--')

ax1.legend()
ax2.plot(temperatures, mean_energies, '-o', markersize=4)
ax3.plot(temperatures, heat_caps, '-o', markersize=4)
#ax4.plot(temperatures, entropies, '-o', markersize=4)
for ipair, pair in enumerate(pairs):
  ax5.plot(temperatures, asr_orders[ipair], label = labels[pair[0]] + '-' + labels[pair[1]])

ax5.legend(loc='upper right')
ax3.tick_params(direction="in")
ax5.tick_params(direction="in")

ax2_diff = np.max(mean_energies) - np.min(mean_energies)
ax3_diff = np.max(heat_caps) - np.min(heat_caps)
ax2.set_ylim(np.min(mean_energies)-0.025*np.abs(ax2_diff), np.max(mean_energies)+0.025*np.abs(ax2_diff))
ax3.set_ylim(np.min(heat_caps)-0.025*np.abs(ax3_diff), np.max(heat_caps)+0.025*np.abs(ax3_diff))

plt.tight_layout()

plt.show()

