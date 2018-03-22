#$Id: make.inc,v 1.1 2009/12/02 10:24:37 ishikawa Exp $

REVISION = -D_REVISION_='"$(shell cat .svnrevision)"'

#
# specify Fortran90 MPI compiler FC, compiler options FFLAGS and linker options LDFLAGS
#

FC     = xmpf90
FFLAGS = -cpp -omp -Kfast,openmp -KXFILL -Kprefetch_sequential=soft

CC     = mpifccpx

LDFLAGS = $(FFLAGS)

#LIBS = 

MPIINC = /opt/FJSVtclang/GM-1.2.0-23/include/mpi/fujitsu
INCLUDE = -I./ -I$(MPIINC)

#.F90.o :

#%.mod : %.F90
#	$(FC) $(FFLAGS) $(REVISION) $(INCLUDE) -c $<

#%.o : %.F90
#	$(FC) $(FFLAGS) $(REVISION) $(INCLUDE) -c $<
