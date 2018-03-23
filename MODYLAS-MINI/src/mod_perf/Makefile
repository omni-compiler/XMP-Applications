#CC = mpicc
#CFLAGS =

#FC = mpif90
#FFLAGS =

LIBS = -lm

C_OBJS = perf.o
F_OBJS = mod_perf.o perf.o

ifeq (, $(findstring -DMAX_SECTIONS=, $(CFLAGS)))
	CFLAGS += -DMAX_SECTIONS=20
endif

c_mpi f_mpi: CFLAGS += -DUSE_MPI
f_mpi: FFLAGS += -DUSE_MPI

all: c

c c_mpi: $(C_OBJS)

c_test: c
	$(CC) $(CFLAGS) $(LDFLAG) -o $@ c_test.c $(C_OBJS) $(LIBS)

c_mpi_test: c_mpi
	$(CC) $(CFLAGS) -DUSE_MPI $(LDFLAG) -o $@ c_test.c $(C_OBJS) $(LIBS)

f f_mpi: $(F_OBJS)

f_test: f
	$(FC) $(FFLAGS) $(LDFLAG) -o $@ f_test.F90 $(F_OBJS) $(LIBS)

f_mpi_test: f_mpi
	$(FC) $(FFLAGS) -DUSE_MPI $(LDFLAG) -o $@ f_test.F90 $(F_OBJS) $(LIBS)


.SUFFIXES:
.SUFFIXES: .f90 .F90 .c .o .mod

.f90.o:
	$(FC) $(FFLAGS) -c $<

.F90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.mod *.i90 c_test c_mpi_test f_test f_mpi_test


# DO NOT DELETE
c_test.o: perf.h
perf.o: perf.h
f_test.o: mod_perf.o

