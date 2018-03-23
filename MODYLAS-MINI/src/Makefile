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

OBJS = fmodules.o main_f90.o app_f90.o \
  md_charmm_f90.o md_fmm_f90.o md_direct_f90.o \
  nve_integrate.o k_energy.o cell_edge.o \
  mt19937ar.o mpitool.o domain_div.o assign2subcell.o \
  shake_rattle_roll.o \
  parse_input.o \
  comm.o comm_3.o comm_fmm.o \
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
app_f90.o: atommass.mod g_main.mod md_condition.mod md_const.mod md_file.mod
app_f90.o: md_fmm_domdiv_flg.mod md_monitors.mod md_periodic.mod md_segment.mod
app_f90.o: mpivar.mod nhc.mod param.mod shakerattleroll.mod trj_mpi.mod
app_f90.o: trj_org.mod unitcell.mod cutoffradius.mod
assign2subcell.o: atommass.mod md_fmm.mod md_fmm_domdiv_flg.mod md_periodic.mod
assign2subcell.o: md_segment.mod mpivar.mod param.mod trj_mpi.mod trj_org.mod
assign2subcell.o: unitcell.mod
cell_edge.o: atommass.mod md_periodic.mod md_segment.mod param.mod trj_org.mod
cell_edge.o: unitcell.mod
comm.o: comm_base.mod comm_bd.mod md_fmm.mod md_fmm_domdiv_flg.mod
comm.o: md_forces.mod md_monitors.mod md_periodic.mod md_segment.mod mpivar.mod
comm.o: param.mod trj_mpi.mod trj_org.mod unitcell.mod ompvar.mod
comm_3.o: comm_base.mod comm_d3.mod md_fmm.mod md_fmm_domdiv_flg.mod
comm_3.o: md_forces.mod md_monitors.mod md_periodic.mod md_segment.mod
comm_3.o: mpivar.mod trj_mpi.mod trj_org.mod unitcell.mod
comm_fmm.o: comm_base.mod md_fmm.mod md_fmm_domdiv_flg.mod mpivar.mod
domain_div.o: cutoffradius.mod md_condition.mod md_fmm.mod md_fmm_domdiv_flg.mod
domain_div.o: md_forces.mod md_multiplestep.mod md_periodic.mod md_segment.mod
domain_div.o: mpivar.mod shakerattleroll.mod trj_mpi.mod trj_org.mod
domain_div.o: unitcell.mod ompvar.mod
k_energy.o: atommass.mod md_const.mod md_fmm.mod md_fmm_domdiv_flg.mod
k_energy.o: ompvar.mod mpivar.mod param.mod trj_mpi.mod trj_org.mod
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
mpitool.o: mpivar.mod
nve_integrate.o: atommass.mod md_condition.mod md_const.mod
nve_integrate.o: md_fmm.mod md_fmm_domdiv_flg.mod md_forces.mod md_monitors.mod
nve_integrate.o: md_multiplestep.mod md_periodic.mod md_segment.mod mpivar.mod
nve_integrate.o: param.mod shakerattleroll.mod trj_mpi.mod unitcell.mod
parse_input.o: atommass.mod md_charmm_lj.mod md_condition.mod md_const.mod
parse_input.o: md_coulomb.mod md_file.mod md_fmm_domdiv_flg.mod
parse_input.o: md_multiplestep.mod md_periodic.mod md_segment.mod mpivar.mod
parse_input.o: nhc.mod param.mod shakerattleroll.mod
parse_input.o: g_main.mod cutoffradius.mod md_fmm.mod md_ewald.mod
parse_input.o: trj_org.mod trj_mpi.mod unitcell.mod
shake_rattle_roll.o: atommass.mod md_fmm.mod md_fmm_domdiv_flg.mod
shake_rattle_roll.o: md_segment.mod mpivar.mod param.mod pshake.mod
shake_rattle_roll.o: pshake_init.mod shakerattleroll.mod trj_mpi.mod trj_org.mod
atommass.mod: .//fmodules.o
comm_base.mod: .//fmodules.o
comm_bd.mod: .//fmodules.o
comm_d3.mod: .//fmodules.o
cutoffradius.mod: .//fmodules.o
g_main.mod: .//fmodules.o
get_wtime.mod: .//fmodules.o
md_charmm_lj.mod: .//fmodules.o
md_condition.mod: .//fmodules.o
md_const.mod: .//fmodules.o
md_coulomb.mod: .//fmodules.o
md_ewald.mod: .//fmodules.o
md_file.mod: .//fmodules.o
md_fmm.mod: .//fmodules.o
md_fmm_domdiv_flg.mod: .//fmodules.o
md_forces.mod: .//fmodules.o
md_monitors.mod: .//fmodules.o
md_multiplestep.mod: .//fmodules.o
md_periodic.mod: .//fmodules.o
md_segment.mod: .//fmodules.o
md_void.mod: .//fmodules.o
mod_wk_fmmewald.mod: .//fmodules.o
mod_wk_k_ene.mod: .//fmodules.o
mpivar.mod: .//fmodules.o
ompvar.mod: .//fmodules.o
nhc.mod: .//fmodules.o
param.mod: .//fmodules.o
pshake.mod: .//fmodules.o
pshake_init.mod: .//fmodules.o
shakerattleroll.mod: .//fmodules.o
trj_mpi.mod: .//fmodules.o
trj_org.mod: .//fmodules.o
unitcell.mod: .//fmodules.o
version.mod: .//fmodules.o

opening.o: version.mod md_file.mod mpivar.mod trj_org.mod ompvar.mod
opening.o: md_fmm.mod md_fmm_domdiv_flg.mod cutoffradius.mod
opening.o: md_condition.mod md_multiplestep.mod md_ewald.mod
closing.o: mpivar.mod
main_f90.o opening.o nve_integrate.o md_fmm_f90.o: timing.h90
ifneq (, $(findstring -DPROF_MAPROF, $(FFLAGS)))
timing.h90: $(MAPROF_LIB)
endif

ConfigReader.o: ConfigReader.h
ConfigReader.o: params.h
configread.mod: ConfigRead.o
parse_input.o: configread.mod
