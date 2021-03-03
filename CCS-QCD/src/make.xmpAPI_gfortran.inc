#$Id: make.inc,v 1.1 2009/12/02 10:24:37 ishikawa Exp $

REVISION = -D_REVISION_='"$(shell cat .svnrevision)"'

#
# specify Fortran90 MPI compiler FC, compiler options FFLAGS and linker options LDFLAGS
#

# example: GNU gfortran, OpenMP enabled.
OMNI_HOME = $(shell xmpcc --show-env | grep OMNI_HOME | sed -e 's/OMNI_HOME=//' -e 's/"//g')
OMNI_INC = -I$(OMNI_HOME)/include
OMNI_LIB = $(shell xmpcc --show-env | grep OMNI_LINKER_OPT | sed -e 's/OMNI_LINKER_OPT=//' -e 's/"//g')

FC      = mpif90
CC      = mpicc
FFLAGS  = -cpp -omp -O2 -ffixed-line-length-132 -march=native -g -fopenmp -mcmodel=medium -funderscoring 
LDFLAGS = $(FFLAGS) $(OMNI_LIB)

INCLUDE = -I./ $(OMNI_INC)
