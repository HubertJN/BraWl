import numpy as np
import matplotlib.pyplot as plt
import math as m
import netCDF4 as nc
np.set_printoptions(suppress=True)

filename = "dos_bins.dat"
bin_edges = nc.Dataset(filename)
bin_edges = np.array(bin_edges["grid data"][:], dtype=np.float64)

filename = "dos_probability.dat"
statP = nc.Dataset(filename)
statP = np.array(statP["grid data"][:], dtype=np.float64)

filename = "bin_probability.dat"
visited_bins = nc.Dataset(filename)
visited_bins = np.array(visited_bins["grid data"][:], dtype=np.float64)

filename = "energy_bias_all.dat"
energy_bias = nc.Dataset(filename)
energy_bias = np.array(energy_bias["grid data"][:], dtype=np.float64)

kb_ev = 8.167333262e-5
ev_to_ry = 13.605693122
kb_ry = kb_ev/ev_to_ry
n_atoms = 256
eV_to_meV = 1000
J = 0.001
N = 6

orig_temp = 400
start_temp = 200
end_temp = 800
step_size = (end_temp-start_temp)/(75)
temperatures = np.arange(start_temp, end_temp, step_size)

# Reduce histogram such that all bins are non zero
bin_width = bin_edges[1] - bin_edges[0]

bin_energies = np.zeros(len(bin_edges[:-1]))
for i, edge in enumerate(bin_edges[:-1]):
    bin_energies[i] = edge + 0.5*bin_width

bin_width = bin_edges[1] - bin_edges[0]

plt.plot(bin_energies, statP, label="Distribution")
for i in range(energy_bias.shape[0]):
    plt.plot(bin_energies, energy_bias[i], label="Reweight: {}".format(i+1))
plt.legend(loc="upper right")
plt.xlabel("Energy (Rydberg)")
plt.ylabel("Bias")
plt.title("Energy bias")
plt.show()

print(visited_bins)
flatness = np.min(visited_bins[visited_bins > 1e-8])/np.mean(visited_bins)
plt.plot(bin_energies, visited_bins)
plt.title("Flatness: {:.3f}".format(flatness))
plt.xlabel("Energy (Rydberg)")
plt.ylabel("Number of visits")
plt.show()

bins = 128
bin_width_con = (bin_edges[-1] - bin_edges[0])/(bins)
bin_edges_con = np.arange(bin_edges[0], bin_edges[-1]+bin_width_con, bin_width_con)

statP_con = np.zeros(len(bin_edges_con)-1)

for i, edge in enumerate(bin_edges[:-1]):
    bin_energy = edge + 0.5*bin_width
    j = 0
    while 1:
        if bin_energy >= bin_edges_con[j] and bin_energy < bin_edges_con[j+1]:
            statP_con[j] += statP[i]
            break
        else:
            j += 1

statP = statP_con/sum(statP_con)
bin_edges = bin_edges_con
print(statP)
# ------------------------------------------------

beta_o = 1.0/(kb_ry*orig_temp) # beta in eV
bin_width = bin_edges[1] - bin_edges[0]

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

#ax2.axvline(og_temp, linestyle='--')
#ax3.axhline(0.0, linestyle='--')
#ax3.axvline(og_temp, linestyle='--')
#ax4.axhline(0.0, linestyle='--')
#ax4.axvline(og_temp, linestyle='--')
# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
mean_errors     = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))

# Loop over temperatures of interest
for itemp, new_temp in enumerate(temperatures):

    beta_n = 1.0/(kb_ry*new_temp)
    partition_function = 0

    # Reweighted histogram
    prob = np.zeros(len(bin_edges)-1)
    
    # Reweight histogram
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        weight = np.exp((beta_o - beta_n)*bin_energy)
        prob[ibin] = statP[ibin]*weight

    # Normalise
    prob = prob/(np.sum(prob*bin_width))

    # Only plot every 5th histogram to avoid crowding the axes
    if itemp%4 == 0:
        strlabel = "T={:8.2f}".format(new_temp)
        ax1.bar(bin_edges[:-1], prob, width=bin_width, align='edge', label=strlabel)
    
    # Mean energy
    mean_energy = np.dot(bin_edges[:-1]+0.5*bin_width, prob)/np.sum(prob)
    mean_energies[itemp] = mean_energy

    # Compute heat capacity using the histogram
    msq_dev = np.zeros(len(bin_edges)-1)
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
        
    heat_caps[itemp] = np.dot(msq_dev, prob)*bin_width/(kb_ry*new_temp**2)

    entropy = (prob*bin_width)*np.log(prob*bin_width)
    entropy[np.isnan(entropy)] = 0
    entropies[itemp] = -np.sum(entropy)

# Complete plots using data computed above
ax1.legend()
ax2.errorbar(temperatures, mean_energies,yerr=mean_errors, fmt='-o', markersize=4)
ax3.plot(temperatures, heat_caps, '-o', label='samples', markersize=4)
ax4.plot(temperatures, entropies, '-o', label='samples', markersize=4)
plt.show()