# Makefile for util_lib
#                     03 Sep 2007 Written  by Takahito Nakajima
#                     28 Sep 2012 Modified by Michio Katouda
#
include ../../config/makeconfig.xmp

ifeq ($(USE_XMP_API),yes)
UTIL_INITMPI=xmpAPI_util_initmpi
UTIL_FINMPI=xmpAPI_util_finmpi
else
UTIL_INITMPI=util_initmpi
UTIL_FINMPI=util_finmpi
endif


objects_common = \
util_linout.o \
util_matout.o \
util_transchar.o \
util_tr1to2.o \
util_lowtrmzero.o \

objects_mpi = \
util_abortmpi.o \
$(UTIL_INITMPI).o \
$(UTIL_FINMPI).o \

objec77 = \

objects = $(objects_common) $(objects_mpi)

all: objcomp

objcomp : $(objects) $(objec77)
	$(RM) libutil.a
	$(AR) $(ARFLAGS) libutil.a $(objects) $(objec77)
	$(RANLIB) libutil.a

clean :
	$(RM) *.o *.mod *.xmod *.a

veryclean :
	$(RM) *.o *.mod *.xmod *.a

.SUFFIXES:.f90
.f90.o: 
	$(F90C) $(F90FLAGS) $(OMPFLAGS) -I$(NTQC_TOP)/src/parallel_mod -c $<

.SUFFIXES:.F90
.F90.o:
	$(F90C) $(F90FLAGS) $(OMPFLAGS) -I$(NTQC_TOP)/src/parallel_mod -c $<

.SUFFIXES:.f
.f.o: 
	$(F77C) $(F77FLAGS) -c $<
