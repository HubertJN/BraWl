Example 1: Simulated Annealing of the CrCoNi system.

This example demonstrates simulated annealing of the CrCoNi system. Atom-atom interactions are taken from:
C. D. Woodgate, J. B. Staunton, Phys. Rev. B 105, 115124 (2023)
DOI: https://doi.org/10.1103/PhysRevB.105.115124

The file 'input.txt' contains all relevant parameters for the simulation, while 'NiCoCr_V_ijs.txt' contains the atom-atom interaction parameters.

Run 'mpirun -n 10 ~/codes/bontewarlo/bontewarlo.run' to launch the simulation. Once it has finished, analyse the results with 01_analyse.py
