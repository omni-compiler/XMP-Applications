#$Id: make.inc,v 1.1 2009/12/02 10:24:37 ishikawa Exp $

REVISION = -D_REVISION_='"$(shell cat .svnrevision)"'

#
# specify Fortran90 MPI compiler FC, compiler options FFLAGS and linker options LDFLAGS
#

# example: GNU gfortran, OpenMP enabled.
FC      = xmpf90
CC      = mpicc
FFLAGS  = -cpp -omp -O2 -ffixed-line-length-132 -march=native -g -fopenmp -mcmodel=medium -funderscoring
LDFLAGS = $(FFLAGS)

INCLUDE = -I./
