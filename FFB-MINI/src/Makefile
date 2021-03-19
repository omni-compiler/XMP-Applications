LD =
LDFLAGS =

include ./make_setting

ifndef LD
LD = $(FC)
endif
ifndef LDFLAGS
LDFLAGS = $(FFLAGS)
endif

LES3X.MPI = ../bin/les3x.mpi
FFB_MINI = ../bin/ffb_mini

VERSION = 1.0.0

FFLAGS += -DFFB_MINI_VERSION=\"$(VERSION)\"

all: $(LES3X.MPI) $(FFB_MINI)

OBJS = \
  caluel.o clrcrs.o \
  csin3x.o datcnv.o dgnscl.o e2plst.o e2pmtr.o \
  errchk.o fild3x.o fld3x2.o \
  icalel.o int3dx.o \
  lesrop.o lesrpx.o lumpex.o \
  match4.o matgau.o mkcrs.o  neibr2.o \
  reordr.o \
  rfname.o sethex.o srfexx.o subcnv.o \
  mfname.o \
  miniapp_util.o

OBJS0 = \
  les3x.o bcgs3x.o bcgsxe.o calax3.o \
  callap.o elm3dx.o grad3x.o lessfx.o \
  nodlex.o pres3e.o rcmelm.o vel3d1.o vel3d2.o
OBJS1 = \
  xmpAPI_les3x.o xmpAPI_bcgs3x.o xmpAPI_bcgsxe.o xmpAPI_calax3.o \
  xmpAPI_callap.o xmpAPI_elm3dx.o xmpAPI_grad3x.o xmpAPI_lessfx.o \
  xmpAPI_nodlex.o xmpAPI_pres3e.o xmpAPI_rcmelm.o xmpAPI_vel3d1.o xmpAPI_vel3d2.o

ifeq ($(USE_XMP_API), yes)
  OBJS += $(OBJS1)
  FBMOBJ = xmpAPI_ffb_mini_main.o
else
  OBJS += $(OBJS0)
  FBMOBJ = ffb_mini_main.o
endif

ifeq (, $(findstring -DNO_METIS, $(FFLAGS)))
OBJS += metis_wrapper.o
endif

ifeq (, $(findstring -DNO_REFINER, $(FFLAGS)))
OBJS += xmpAPI_lrfnms.o xmpAPI_extrfn.o
endif

param.h: param.h.in
	./make_param_h.sh $< > $@

ifneq (, $(findstring -DPROF_MAPROF, $(FFLAGS)))

MAPROF_DIR = ma_prof/src
MAPROF_LIB = $(MAPROF_DIR)/libmaprof_f.a

FFLAGS += -I$(MAPROF_DIR)
LDFLAGS += -L$(MAPROF_DIR)
LIBS += -lmaprof_f

$(OBJS): $(MAPROF_LIB)

export
$(MAPROF_LIB):
	$(MAKE) -C $(MAPROF_DIR) f_mpi MAPROF_F="FC FFLAGS" MAX_SECTIONS=30

endif

LIB_GF2 = gf2/libgf2.a
LIB_DD_MPI = dd_mpi/libdd_mpi.a
#LIB_DD_MPI = dd_mpi/dd_mpi.o dd_mpi/ddcom4.o

$(LIB_GF2):
	$(MAKE) -C gf2

$(LIB_DD_MPI):
	$(MAKE) -C dd_mpi

$(LES3X.MPI): $(LIB_GF2) $(LIB_DD_MPI) $(OBJS) les3x_main.o
	$(LD) $(LDFLAGS) -o $@ $(OBJS) les3x_main.o $(LIB_GF2) $(LIB_DD_MPI) $(LIBS)

$(FFB_MINI): $(LIB_GF2) $(LIB_DD_MPI) $(OBJS) $(FBMOBJ) makemesh.o
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(FBMOBJ) makemesh.o $(LIB_GF2) $(LIB_DD_MPI) $(LIBS)

.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .c .o

.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<
.F.o:
	$(FC) $(FFLAGS) -c $<
.F90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -rf *.o *.mod *.xmod *.lst param.h
	$(MAKE) -C gf2 clean
	$(MAKE) -C dd_mpi clean
ifneq (, $(findstring -DPROF_MAPROF, $(FFLAGS)))
	$(MAKE) -C $(MAPROF_DIR) clean
endif

distclean: clean
	rm -rf ../bin/*
	$(MAKE) -C ../test clean

test: $(FFB_MINI)
	$(MAKE) -C ../test

test_k: $(FFB_MINI)
	$(MAKE) -C ../test test_k

test_fx10: $(FFB_MINI)
	$(MAKE) -C ../test test_fx10


xmpAPI_ffb_mini_main.o: param.h makemesh.o
ffb_mini_main.o: param.h makemesh.o

xmpAPI_les3x.o: timing.h
xmpAPI_vel3d1.o: timing.h
xmpAPI_bcgs3x.o: timing.h
xmpAPI_pres3e.o: timing.h
#bcgs3x.o: timing.h
xmpAPI_rcmelm.o: timing.h
xmpAPI_bcgsxe.o: timing.h
xmpAPI_callap.o: timing.h
xmpAPI_calax3.o: timing.h
xmpAPI_grad3x.o: timing.h
