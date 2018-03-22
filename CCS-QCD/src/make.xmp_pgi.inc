#$Id: make.inc,v 1.1 2009/12/02 10:24:37 ishikawa Exp $

REVISION = -D_REVISION_='"$(shell cat .svnrevision)"'

#
# specify Fortran90 MPI compiler FC, compiler options FFLAGS and linker options LDFLAGS
#

FC        = xmpf90
CC        = mpicc

FFLAGS   += -cpp -omp -O2 -Mpreprocess -m64 -Mfree -mp
CFLAGS    = -O2 -Mpreprocess -m64 -mp

ACCFLAGS  = -U_OPENACC
FFLAGS   += $(ACCFLAGS)

LDFLAGS   = $(FFLAGS)

MPIINC    = /opt/pgi/linux86-64/2017/mpi/openmpi-2.1.2/include
INCLUDE   = -I./ -I$(MPIINC)
