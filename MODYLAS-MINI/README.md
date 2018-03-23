MODYLAS-MINI
============

* version: 1.0.0 (based on MODYLAS 1.0.1)
* date: 2014/10/02
* contact: miniapp@riken.jp


About MODYLAS
-------------

This mini app is based on a general-purpose molecular dynamics (MD) simulation program, MODYLAS, which has been developed in Nagoya University and Institute for Molecular Science.
MODYLAS evaluates long range forces by Fast Multipole Method (FMM), 
and is appropriate for large-scale systems using massively parallel supercomputers.
Refer to the [MODYLAS Web site](http://www.modylas.org/) for details.

Contact point of orignal MODYLAS programs is: 
  Dr. Yoshimichi Andoh <yoshimichi.andoh@apchem.nagoya-u.ac.jp>


Compilation
-----------

To build and run the program, Fortran90 and C compilers that support OpenMP, MPI library and GNU Make are required.

 1. Obtain the package and extract its contents.

 2. Go to the src directory and edit the file "make_setting" according to
    your compilers and libraries.
    There are several example "make_setting.*" files:
    * make_setting.intel : For using Intel compilers
    * make_setting.gcc   : For using GCC compilers
    * make_setting.fx10  : For K and FX10-based systems

 3. Run make command in the src directory.
    After successful compilation, an executable file named `modylas_mini` is
    generated in the src directory.


Testing
-------

A test shell script is provided in the test directory.
To run the test interactively, simply run "./go.sh" in the the directory, or run "make test" in the src directory.
This test script runs the program `modylas_mini` with 8 MPI processes and 2 OpenMP threads per process, and compares the computed results with reference data.


How to run the program
----------------------

### Input Files

Three files are needed to run the program.

* Configuration file

    A text file specifying runtime control parameters.
    File name must be `_session_.mdconf`.

* Force field file

    A binary file containing force field and other parameters required for MD simulations.
    File name must be `_session_.mdff.bin`

* Coordinate file

    A binary file containg initial coordinate and velocity data of atoms.
    File name must be `_session_.mdxyz.bin`

For each file name, `_session_` is a "session name", an arbitrary character string to distingush runtime conditions.


### Running

To run the program, a "session name" must be given as a command line parameter.
In the following example, we run the session named `wat111`, with 64 MPI processes and 8 OpenMP threads per process.

    $ export OMP_NUM_THREADS=8
    $ mpiexec -n 64 path_to_src_directory/modylas_mini wat111   


### Output files

* Monitor file

    Output monitoring physical variables to a text file named `_session_.mdmtr`.

* Trajectory file

    Output coordinates and velocities of atoms to a binary file named `_session_.mdtrj.bin`.

* Restart file

    Output restart data to a binary file named `_session_.mdrestart.bin`.

For each output file, the first time step of output, and output step interval can be specified in the configuration file.


### Notes

#### Incompatibility of configuration file with the original MODYLAS

This mini app's configuration file `_session_.mdconf` and original MODYLAS's one `_session_.mddef` are *not* compatible.
You cannot use `_session_.mdconf` file to run the original MODYLAS program, and vice versa.

#### Restrictions on the number of MPI processes 

The number of MPI processes must be:

  - less than or equal to the number of cells (= `ncell`**3),
  - greater than 8,
  - and power of 2.

#### Running on big endian architecture

In the sample input files contained in this package, binary data are stored in little endian format.
Therefore, to read these files on a big endian architecture, you must use endian conversion functionality provided by the system.
For example, on K computer and FX10-based systems, you need to set a environmental variable `FORT90L='-Wl,-T'` before running the program.


Configuration File Format
-------------------------

Runtime control parameters are specified in the configuration file `_session_.mdconf`.
Each line in the file, a single parameter is assigned in the form  `name=value`.
White spaces at both ends of `name` and `value` are ignored.
Blank lines and comments, beginning at the # caracter to the end of the line, are also ignored.

### Output

    mntr_start(int):       monitor file output start time step (default 0)
    mntr_interval(int):    monitor file output step intergval (default 1)
    trj_start(int):        trajectory file output start time step (default 0)
    trj_interval(int):     trajectory file output step intergval (default 1)
    restart_start(int):    restart file output start time step (default 0)
    restart_interval(int): restart file output step intergval (default 1)

### Initial condition

    maxwell_velocities(bool): reset atomic velocities using Maxwell distribution (default false)
    temperature(real):        temperature(K) (default 300.0)
    randomseed(int):          seed of pseudorandom numbers (default 1235)

### Time step

    dt(real):                time step length(sec) (mandatory)
    step(int):               # of step in MD calculation (mandatory)
    nstep_skip_middle(int):  see "Multiple time step" bellow (default 1)
    nstep_skip_long(int):    see "Multiple time step" bellow (default 1)

### Domain decomposition

    manual_division(bool):  flag for manual devision (default false)
    nxdiv(int):             # of MPI processes in the x-direction
    nydiv(int):             # of MPI processes in the y-direction
    nzdiv(int):             # of MPI processes in the z-direction

    nxdiv, nydiv and nzdiv are manadtory for manual devision.

### FMM

    ncell(int):             # of cells in each direction (mandatory)
    nmax(int):              order of multipole expansion (default 4)
    ULswitch(int):          see "Communications" below (default 1)

### SHAKE/RATTLE

    shake_max_iteration(int): max. number of iterations (default 100)
    shake_tolerance(real):    convergence tolerance (default 1.0e-8)

### Other parameters

    cutoff(real):             cut-off distance(A) (mandatory)
    ewald_surface_term(bool): flag for Ewald surface term (default false)


Sample Input data
-----------------

Three sample data sets are provided in the data directory.
In each data set directory, besides input data files, a sample job script for K computer and example output files are also contained.

### wat111

 - number of atoms:    19,530
 - cell division:      8x8x8
 - FMM tree levels:    3

### wat222

 - number of atoms:    156,240
 - cell division:      16x16x16
 - FMM tree levels:    4

### wat444

 - number of atoms:    1,249,920
 - cell division:      32x32x32
 - FMM tree levels:    5


To get high performance runs on these datadaets, refer to notes in "performance.md".

If you want to run the program on other data sets, consult the developers of original MODYLAS.


Differences from the original MODYLAS program
---------------------------------------------

* Reduced the size of the code by limiting the functionality to MD calculations of water molecule systems in the microcanonical (NEV) ensembles.

* Other code clean-up.


Target exa-scale problem setting
--------------------------------

The current model setting targeted at exa-scale simulations is as follows:

 - 10^9 atoms
 - 10^9 time steps

The target performance of this simulation is to compute the above model within 150 hours.


Further details
---------------


### Multi time step

related parameters:

    config. parameter    valiable name
    --------------------+---------------------------
    steps                md_condition__howmany_steps
    nstep_skip_middle    maxMTm
    nstep_skip_long      maxMTl
    dt                   dt

force calculations:

  - for every steps

    calculate short range force (bond, angle, ...).
    (nothing to do for mini app)

  - every `maxMTm` steps

    calculate middle range forces (direct interaction part of FMM, FdW forces).

  - every `maxMTm`*`maxMTl` steps

    calculate long range forces (multipole interaction part of FMM).


Note that number of steps `step`(`md_condition__howmany_steps`) and time step length `dt` are based on the cycle of long range force calculations.
In the unit of short range force calculation cycle:

    nuber of steps   = md_condition__howmany_steps * (maxMTm * maxMTl)
    time step length = dt / (maxMTm * maxMTl) 


### Communications

There are four major parts of the inter-node communications.

- COMM_DIRECT:

    communications of atomic coodinate data in halo cells between neighbor nodes, before direct force calculations at every time step, implemented in the subroutine `comm_direct_3`.

- MIGRATION:

    communications of the coodinate and velocity data of the atoms migrate to the cells which belonging to other process node, for every `maxMTm`*`maxMTl` step, implemented in the subroutine `comm_bound`.

- COMM_FMM:

    communications of the multipole moments of the cells in each level, for every `maxMTm`*`maxMTl` step, implemented in the subroutine `comm_fmm_local_multi` (for lower levels) and `comm_fmm_local_top` (for higher levels).

- ENE_REDUCTON:

    all-reduce communications for thermodynamical variables, for every time step.

In the current implementation of the program, inter-node communications are  optimaized for the 3-D torus networks.
To achieve the best performance on those systems, refer to notes in "performance.md".

Note that the MODYLAS-MINI does not rely on the action-reaction law to calculate direct forces.
Therefore, the communications of force data after direct force calculations are not necessary.


References
----------

*  Andoh, Y. et al., "MODYLAS: A Highly Parallelized General-Purpose
   Molecular Dynamics Simulation Program for Large-Scale Systems with
   Long-Range Forces Calculated by Fast Multipole Method (FMM) and
   Highly Scalable Fine-Grained New Parallel Processing Algorithms",
   J. Chem. Theory Comput., 2013, 9 (7), pp 3201-3209,
   DOI: 10.1021/ct400203a.
