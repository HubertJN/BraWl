# Bonte Warlo

A code for performing lattice-based Monte Carlo simulations of alloys with an internal energy given by a Bragg-Williams Hamiltonian. The code will be periodically updated so it is best to check the GitHub repository.

## Citations
ANY publications/presentations/further work resulting from the use of this software should cite the original publication for which it was developed:
C. D. Woodgate, J. B. Staunton, Phys. Rev. B 105, 115124 (2023)
DOI: [10.1103/PhysRevB.105.115124](https://doi.org/10.1103/PhysRevB.105.115124)

In addition, if you use the bcc implementation, you should cite the original paper for which that lattice type was implemented:
C. D. Woodgate, J. B. Staunton, Phys. Rev. Mater. 107, 013801 (2023)
DOI: [10.1103/PhysRevMaterials.7.013801](https://doi.org/10.1103/PhysRevMaterials.7.013801)

## Compilation
At the moment the code is only tested with gfortran and OpenMPI. Put the code in a directory like `~/codes/bontewarlo`. It is my intention to test other compilers in future: watch this space!

On the Warwick SCRTP system, which uses environment modules, run
```
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 netCDF-Fortran/4.6.0
```
then you should be able to build the code with
```
make compiler=gfortran
```

## Running the code
If you navigate to the `examples` subdirectory, you should find two examples demonstrating the code's usage which can be run inside those directories.

Most of the options specified in the input file are fairly self-explanatory. The least obvious is the `mode' option. Because it is my intention to include a 2D (and potentially 1D) option in future, the first digit indicates the number of spatial dimensions for the simulation. Then the last two digits the mode. At present, the implemented options are:
- 01: Simulated Annealing. Uses the Metropolis Monte Carlo algorithm with Kawasaki dynamics to perform simulated annealing on a system in an initially random configuration.
- 02: Draw Decorellated Samples. Optionally performs simulated annealing then draws samples of the grid N Monte Carlo steps apart. Good for generating supercell configurations for use other methods. *E.g.* this recent reference where the code was used to generate training/test data for a machine-learned interatomic potential: [arXiv:230908689](https://doi.org/10.48550/arXiv.2309.08689).

## Contributing
Any/all contributions are welcome via pull requests. 

## License
This software is released under the MIT license. See the file LICENSE.txt for details.

## Questions, Brickbats, and Bouquets
[Drop me an email](mailto:christopher.woodgate@physics.org) if you have any questions, comments, or feedback!
