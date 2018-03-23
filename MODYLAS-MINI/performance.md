A guideline for performance measurement of the MODYLAS mini app
===============================================================

For high performance, the mini app is recommended to be evaluated under the following conditions:

- Execution on parallel systems with 3-D torus interconnect, such as the K computer, FX10, and Blue Gene. To obtain the best performance, logical 3D configuration of MPI processes, which is determined by .mdconf input file, must correspond to the physical 3D topology of compute nodes. On the K computer and FX10 systems, batch job scripts can be used to configure the node topology of submitted jobs. 

- Optimization of the length of the cubic sub-cells. Sub-cell length must be greater than half of the Lennard-Jones cutoff length. To obtain the best performance, it should be as close as possible to half of the cutoff length.

- Proper distribution of the atom numbers per CPU core. On the K computer and other systems with similar architecture specifications, the number of atoms per core should be greater than 20 to attain MPI parallelization efficiency as reported in reference article [JCTC,9,3201(2013)].

- Switching of MPI communication codes for FMM upper/lower levels. Timing of this switching is determined by 'ULswitch' variable in .mdconf input file.  Its recommended value is 0, 1, 2, 3 for ncell=8, 16, 32, 64.
