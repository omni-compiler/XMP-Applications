#$Id: make.inc,v 1.1 2009/12/02 10:24:37 ishikawa Exp $

REVISION = -D_REVISION_='"$(shell cat .svnrevision)"'

#
# specify Fortran90 MPI compiler FC, compiler options FFLAGS and linker options LDFLAGS
#
USE_XMP_API = yes
OMNI_HOME = $(shell xmpcc --show-env | grep OMNI_HOME | sed -e 's/OMNI_HOME=//' -e 's/"//g')
OMNI_INC = -I$(OMNI_HOME)/include
OMNI_LIB = $(shell xmpcc --show-env | grep OMNI_LINKER_OPT | sed -e 's/OMNI_LINKER_OPT=//' -e 's/"//g')

FC     = mpifrtpx
FFLAGS = -cpp -Kfast,openmp -KXFILL -Kprefetch_sequential=soft

CC     = mpifccpx

LDFLAGS = $(FFLAGS)
LIBS = $(OMNI_LIB)

MPIINC = /opt/FJSVxtclanga/tcsds-1.2.30a/include/mpi/fujitsu
INCLUDE = -I./ $(OMNI_INC)

#.F90.o :

#%.mod : %.F90
#	$(FC) $(FFLAGS) $(REVISION) $(INCLUDE) -c $<

#%.o : %.F90
#	$(FC) $(FFLAGS) $(REVISION) $(INCLUDE) -c $<
