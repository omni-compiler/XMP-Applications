!----------------------------------------------------------------------
! Copyright (C) 2003-2014 Kensuke Iwahashi, Noriyuki Yoshii,
!                         Atsushi Yamada, Yoshimichi Andoh,
!                         Kazushi Fujimoto, Hidekazu Kojima,
!                         Fumiyasu Mizutani, and Susumu Okazaki
! All Rights Reserved.
!
! Copyright (C) 20013-2014 RIKEN AICS
! All Rights Reserved.
!
! This MD program has been developed at Nagoya University, and
! Institute for Molecular Science, National Institutes of Natural
! Sciences.
! And this work was supported by
!    Next-Generation Supercomputer Project, and
!    NAREGI Nanoscience Project,
! Ministry of Education, Culture, Sports, Science and Technology,
! Japan.
!
! This program is NOT a free software and distributed under the
! license described in the LICENSE.
! All rights are reserved by the authors of this program.
!
! The authors do NOT warrant or assume any legal liability or
! responsibility for the accuracy or completeness.
!----------------------------------------------------------------------
subroutine opening
    use version, only: MODYLAS_VERSION
    use md_file, only: session_name
    use mpivar, only: myrank, nprocs
    use ompvar, only: nomp
    implicit none

    character(8) :: date
    character(10) :: time
    include 'mpif.h'
    character(MPI_MAX_PROCESSOR_NAME) :: name
    integer :: name_len, ierr

    if (myrank /= 0) return

    call date_and_time(date, time)

    write(*, '(a)') repeat('*', 50)

    write(*, '(2x, a)') 'MODYLAS_MINI ' // MODYLAS_MINI_VERSION &
              // ' (based on MODYLAS ' // trim(MODYLAS_VERSION) // ')'
    write(*, '(4x, "date: ", a, "/", a, "/", a)') date(1:4), date(5:6), date(7:8)
    write(*, '(4x, "time: ", a, ":", a, ":", a)') time(1:2), time(3:4), time(5:6)
    call MPI_Get_processor_name(name, name_len, ierr)
    write(*, '(4x, "host: ", a)') trim(name)
    write(*, '(4x, "nodes: ", i0)') nprocs
    write(*, '(4x, "threads: ", i0)') nomp
    write(*, '(4x, "session: ", a)') trim(session_name)

    write(*, '(a)') repeat('*', 50)

end subroutine opening
!----------------------------------------------------------------------
subroutine print_params
    use mpivar
    use trj_org, only: n
    use md_fmm, only: ncell, nlevel, nmax, lgflg
    use md_fmm_domdiv_flg, only: nxdiv, nydiv, nzdiv
    use cutoffradius, only: cutrad
    use md_condition, only: md_condition__howmany_steps, dt
    use md_multiplestep, only: maxMTm, maxMTl
    use md_ewald, only: ewald_sterm
    use md_file, only: session_name
#include "timing.h90"

    implicit none

    if (myrank == 0) then
      write(*,'(/,a)') '**** parameters'
      write(*, '("  number of atoms:    ", i0)') n
      write(*, '("  number of cells:    ", i0, " x ", i0, " x ", i0)')  &
          ncell, ncell, ncell
      write(*, '("  number of nodes:    ", i0, " x ", i0, " x ", i0)')  &
          nxdiv, nydiv, nzdiv
      write(*, '("  FMM levels:         ", i0)') nlevel
      write(*, '("  expansion digree:   ", i0)') nmax
      write(*, '("  ULseitch:           ", i0)') lgflg
      write(*, '("  cut-off length[A]:  ", f0.3)') cutrad*1.0d10
      write(*, '("  Ewald surface term: ", l1)') ewald_sterm
      write(*, '("  nsteps:             ", i0)') md_condition__howmany_steps
      write(*, '("  dt[sec]:           ", 1pe10.3)') dt
      write(*, '("  nstep_skip_middle:  ", i0)') maxMTm
      write(*, '("  nstep_skip_long:    ", i0)') maxMTl
    end if

#ifdef PROF_MAPROF
    call maprof_setup("MODYLAS MINI", MODYLAS_MINI_VERSION)
    call maprof_add_section("Main_Loop", TM_MAIN_LOOP)
    call maprof_add_section("Energy_Direct", TM_ENERGY_DIRECT)
    call maprof_add_section("Comm_Direct", TM_COMM_DIRECT)
    call maprof_add_section("Migration", TM_MIGRATION)
    call maprof_add_section("Ene_Reduction", TM_ENE_REDUCTION)
    call maprof_add_section("FMM", TM_FMM)
    call maprof_add_section("Particle2M", TM_P2M)
    call maprof_add_section("M2M", TM_M2M)
    call maprof_add_section("M2L+L2L", TM_M2L_L2L)
    call maprof_add_section("L2Particle", TM_L2P)
    call maprof_add_section("Comm_FMM", TM_COMM_FMM)
    call maprof_profile_add_problem_size("n_atom", n)
    call maprof_profile_add_problem_size("ncell", ncell)
    call maprof_profile_add_problem_size("nxdiv", nxdiv)
    call maprof_profile_add_problem_size("nydiv", nydiv)
    call maprof_profile_add_problem_size("nzdiv", nzdiv)
    call maprof_profile_add_str("session", trim(session_name))
    call maprof_profile_add_int("ULseitch", lgflg)
    call maprof_profile_add_int("nstep", md_condition__howmany_steps)
#endif

end subroutine print_params
