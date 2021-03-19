include ./make_setting

# comm options
FFLAGS += -DCOMM_CUBE
FFLAGS += -DSYNC_COM

# timing option
FFLAGS += -DPROF_MAPROF

# timing option (fapp on K/FX10)
#FFLAGS += -DPROF_FAPP

PROGRAM = modylas_mini
VERSION = 1.0.0

FFLAGS += -DMODYLAS_MINI_VERSION=\"$(VERSION)\"

CFLAGS = -DUSE_MPI

ifeq ($(USE_XMP_API), yes)
FMODULES = xmpAPI_fmodules
APP_F90 = xmpAPI_app_f90
PARSE_INPUT = xmpAPI_parse_input
COMM = xmpAPI_comm
COMM_3 = xmpAPI_comm_3
COMM_FMM = xmpAPI_comm_fmm
NVE_INTEGRATE = xmpAPI_nve_integrate
K_ENERGY = xmpAPI_k_energy
DOMAIN_DIV = xmpAPI_domain_div
MPITOOL = xmpAPI_mpitool
else
FMODULES = fmodules
APP_F90 = app_f90
PARSE_INPUT = parse_input
COMM = comm
COMM_3 = comm_3
COMM_FMM = comm_fmm
NVE_INTEGRATE = nve_integrate
K_ENERGY = k_energy
DOMAIN_DIV = domain_div
MPITOOL = mpitool
endif

OBJS = $(FMODULES).o main_f90.o $(APP_F90).o \
  md_charmm_f90.o md_fmm_f90.o md_direct_f90.o \
  $(NVE_INTEGRATE).o $(K_ENERGY).o cell_edge.o \
  mt19937ar.o $(MPITOOL).o $(DOMAIN_DIV).o assign2subcell.o \
  shake_rattle_roll.o \
	$(PARSE_INPUT).o $(COMM).o $(COMM_3).o $(COMM_FMM).o \
	opening.o closing.o ConfigReader.o ConfigRead.o

all: $(PROGRAM)

ifneq (, $(findstring -DPROF_MAPROF, $(FFLAGS)))

MAPROF_DIR = ma_prof/src
MAPROF_LIB = $(MAPROF_DIR)/libmaprof_f.a

FFLAGS += -I$(MAPROF_DIR)
LDFLAGS += -L$(MAPROF_DIR)
LIBS += -lmaprof_f

export
$(MAPROF_LIB):
	$(MAKE) -C $(MAPROF_DIR) f_mpi MAPROF_F="FC FFLAGS"

$(PROGRAM): $(MAPROF_LIB) $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

else

$(PROGRAM): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

endif


.SUFFIXES:
.SUFFIXES: .f .f90 .c .o

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.mod *.xmod
ifneq (, $(findstring -DPROF_MAPROF, $(FFLAGS)))
	$(MAKE) -C $(MAPROF_DIR) clean
endif

distclean: clean
	rm -f $(PROGRAM) ../test/log ../test/wat111.mdmntr \
	    ../test/wat111.mdtrj.bin ../test/wat111.restart.bin \
		../test/maprof_output.yaml

test: $(PROGRAM)
	cd ../test && ./go.sh


# DO NOT DELETE THIS LINE - used by make depend
$(APP_F90).o: atommass.mod g_main.mod md_condition.mod md_const.mod md_file.mod
$(APP_F90).o: md_fmm_domdiv_flg.mod md_monitors.mod md_periodic.mod md_segment.mod
$(APP_F90).o: mpivar.mod nhc.mod param.mod shakerattleroll.mod trj_mpi.mod
$(APP_F90).o: trj_org.mod unitcell.mod cutoffradius.mod
assign2subcell.o: atommass.mod md_fmm.mod md_fmm_domdiv_flg.mod md_periodic.mod
assign2subcell.o: md_segment.mod mpivar.mod param.mod trj_mpi.mod trj_org.mod
assign2subcell.o: unitcell.mod
cell_edge.o: atommass.mod md_periodic.mod md_segment.mod param.mod trj_org.mod
cell_edge.o: unitcell.mod
$(COMM).o: comm_base.mod comm_bd.mod md_fmm.mod md_fmm_domdiv_flg.mod
$(COMM).o: md_forces.mod md_monitors.mod md_periodic.mod md_segment.mod mpivar.mod
$(COMM).o: param.mod trj_mpi.mod trj_org.mod unitcell.mod ompvar.mod
$(COMM_3).o: comm_base.mod comm_d3.mod md_fmm.mod md_fmm_domdiv_flg.mod
$(COMM_3).o: md_forces.mod md_monitors.mod md_periodic.mod md_segment.mod
$(COMM_3).o: mpivar.mod trj_mpi.mod trj_org.mod unitcell.mod
$(COMM_FMM).o: comm_base.mod md_fmm.mod md_fmm_domdiv_flg.mod mpivar.mod
$(PARSE_INPUT).o: atommass.mod md_charmm_lj.mod md_condition.mod md_const.mod
$(PARSE_INPUT).o: md_coulomb.mod md_file.mod md_fmm_domdiv_flg.mod
$(PARSE_INPUT).o: md_multiplestep.mod md_periodic.mod md_segment.mod mpivar.mod
$(PARSE_INPUT).o: nhc.mod param.mod shakerattleroll.mod
$(PARSE_INPUT).o: g_main.mod cutoffradius.mod md_fmm.mod md_ewald.mod
$(PARSE_INPUT).o: trj_org.mod trj_mpi.mod unitcell.mod
$(DOMAIN_DIV).o: cutoffradius.mod md_condition.mod md_fmm.mod md_fmm_domdiv_flg.mod
$(DOMAIN_DIV).o: md_forces.mod md_multiplestep.mod md_periodic.mod md_segment.mod
$(DOMAIN_DIV).o: mpivar.mod shakerattleroll.mod trj_mpi.mod trj_org.mod
$(DOMAIN_DIV).o: unitcell.mod ompvar.mod
$(K_ENERGY).o: atommass.mod md_const.mod md_fmm.mod md_fmm_domdiv_flg.mod
$(K_ENERGY).o: ompvar.mod mpivar.mod param.mod trj_mpi.mod trj_org.mod
main_f90.o: md_condition.mod md_fmm_domdiv_flg.mod md_forces.mod
main_f90.o: md_monitors.mod md_multiplestep.mod md_segment.mod mpivar.mod
main_f90.o: shakerattleroll.mod trj_mpi.mod version.mod
md_charmm_f90.o: md_ewald.mod md_fmm_domdiv_flg.mod
md_charmm_f90.o: md_forces.mod md_monitors.mod md_segment.mod
md_charmm_f90.o: trj_mpi.mod ompvar.mod
md_direct_f90.o: cutoffradius.mod md_charmm_lj.mod md_const.mod
md_direct_f90.o: md_coulomb.mod md_fmm.mod md_fmm_domdiv_flg.mod md_forces.mod
md_direct_f90.o: md_monitors.mod md_periodic.mod md_segment.mod ompvar.mod
md_direct_f90.o: md_void.mod mpivar.mod param.mod trj_mpi.mod trj_org.mod
md_fmm_f90.o: comm_base.mod md_const.mod md_coulomb.mod md_fmm.mod
md_fmm_f90.o: md_fmm_domdiv_flg.mod md_forces.mod md_monitors.mod
md_fmm_f90.o: md_periodic.mod md_segment.mod mod_wk_fmmewald.mod mpivar.mod
md_fmm_f90.o: param.mod trj_mpi.mod trj_org.mod unitcell.mod ompvar.mod
$(MPITOOL).o: mpivar.mod
$(NVE_INTEGRATE).o: atommass.mod md_condition.mod md_const.mod
$(NVE_INTEGRATE).o: md_fmm.mod md_fmm_domdiv_flg.mod md_forces.mod md_monitors.mod
$(NVE_INTEGRATE).o: md_multiplestep.mod md_periodic.mod md_segment.mod mpivar.mod
$(NVE_INTEGRATE).o: param.mod shakerattleroll.mod trj_mpi.mod unitcell.mod
shake_rattle_roll.o: atommass.mod md_fmm.mod md_fmm_domdiv_flg.mod
shake_rattle_roll.o: md_segment.mod mpivar.mod param.mod pshake.mod
shake_rattle_roll.o: pshake_init.mod shakerattleroll.mod trj_mpi.mod trj_org.mod
atommass.mod: .//$(FMODULES).o
comm_base.mod: .//$(FMODULES).o
comm_bd.mod: .//$(FMODULES).o
comm_d3.mod: .//$(FMODULES).o
cutoffradius.mod: .//$(FMODULES).o
g_main.mod: .//$(FMODULES).o
get_wtime.mod: .//$(FMODULES).o
md_charmm_lj.mod: .//$(FMODULES).o
md_condition.mod: .//$(FMODULES).o
md_const.mod: .//$(FMODULES).o
md_coulomb.mod: .//$(FMODULES).o
md_ewald.mod: .//$(FMODULES).o
md_file.mod: .//$(FMODULES).o
md_fmm.mod: .//$(FMODULES).o
md_fmm_domdiv_flg.mod: .//$(FMODULES).o
md_forces.mod: .//$(FMODULES).o
md_monitors.mod: .//$(FMODULES).o
md_multiplestep.mod: .//$(FMODULES).o
md_periodic.mod: .//$(FMODULES).o
md_segment.mod: .//$(FMODULES).o
md_void.mod: .//$(FMODULES).o
mod_wk_fmmewald.mod: .//$(FMODULES).o
mod_wk_k_ene.mod: .//$(FMODULES).o
mpivar.mod: .//$(FMODULES).o
ompvar.mod: .//$(FMODULES).o
nhc.mod: .//$(FMODULES).o
param.mod: .//$(FMODULES).o
pshake.mod: .//$(FMODULES).o
pshake_init.mod: .//$(FMODULES).o
shakerattleroll.mod: .//$(FMODULES).o
trj_mpi.mod: .//$(FMODULES).o
trj_org.mod: .//$(FMODULES).o
unitcell.mod: .//$(FMODULES).o
version.mod: .//$(FMODULES).o

opening.o: version.mod md_file.mod mpivar.mod trj_org.mod ompvar.mod
opening.o: md_fmm.mod md_fmm_domdiv_flg.mod cutoffradius.mod
opening.o: md_condition.mod md_multiplestep.mod md_ewald.mod
closing.o: mpivar.mod
main_f90.o opening.o $(NVE_INTEGRATE).o md_fmm_f90.o: timing.h90
ifneq (, $(findstring -DPROF_MAPROF, $(FFLAGS)))
timing.h90: $(MAPROF_LIB)
endif

ConfigReader.o: ConfigReader.h
ConfigReader.o: params.h
configread.mod: ConfigRead.o
$(PARSE_INPUT).o: configread.mod
