import numpy as np
import matplotlib.pyplot as plt
import math as m
import netCDF4 as nc
np.set_printoptions(suppress=True)

directory = input("Input directory to pull data from: ")

filename = "{}/wl_dos_bins.dat".format(directory)
bin_edges = nc.Dataset(filename)
bin_edges = np.array(bin_edges["grid data"][:], dtype=np.float64)

filename = "{}/wl_dos.dat".format(directory)
wl_logdos = nc.Dataset(filename)
wl_logdos = np.array(wl_logdos["grid data"][:], dtype=np.float64)

wl_logdos = np.log(np.exp(wl_logdos)/np.sum(np.exp(wl_logdos)))
print(wl_logdos)

filename = "{}/wl_hist.dat".format(directory)
wl_hist = nc.Dataset(filename)
wl_hist = np.array(wl_hist["grid data"][:], dtype=np.float64)

kb_ev = 8.167333262e-5
ev_to_ry = 13.605693122
kb_ry = kb_ev/ev_to_ry
n_atoms = 256
eV_to_meV = 1000
J = 0.001
N = 6

orig_temp = 400
start_temp = 100
end_temp = 1000
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
#exit()

# Setup plots
fig, [ax1, ax2, ax3, ax4] = plt.subplots(1, 4, figsize=(24, 7), constrained_layout=True)

ax1.set_xlabel('energy U (Rydbergs)')
ax1.set_ylabel('P(U)')
ax1.set_title("MUCA energy histograms")

ax2.set_xlabel('temperature T (K)')
ax2.set_ylabel('<U> (Ry)')
ax2.set_title('Mean energy vs temperature')  

ax3.set_xlabel('temperature T (K)')
ax3.set_ylabel('Cv (Ry/K)')
ax3.set_title('Heat capacity vs temperature') 

ax4.set_xlabel('temperature T (K)')
ax4.set_ylabel('S (kB^-1)')
ax4.set_title('Entropy vs temperature')

# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
mean_errors     = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))

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
    if itemp%4 == 0:
        strlabel = "T={:8.2f}".format(new_temp)
        ax1.bar(bin_edges[:-1], prob, width=bin_width, align='edge', label=strlabel)
    
    # Mean energy
    mean_energy = np.dot(bin_edges[:-1]+0.5*bin_width, prob)
    mean_energies[itemp] = mean_energy

    # Compute heat capacity using the histogram
    msq_dev = np.zeros(len(bin_edges)-1)
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
        
    heat_caps[itemp] = np.dot(msq_dev, prob)*bin_width/(kb_ry*new_temp**2)

    entropy = prob*np.log(prob)
    entropy[np.isnan(entropy)] = 0
    entropies[itemp] = -np.sum(entropy)

# Complete plots using data computed above
ax1.legend()
ax2.errorbar(temperatures, mean_energies,yerr=mean_errors, fmt='-o', markersize=4)
ax3.plot(temperatures, heat_caps, '-o', label='samples', markersize=4)
ax4.plot(temperatures, entropies, '-o', label='samples', markersize=4)
plt.show()
