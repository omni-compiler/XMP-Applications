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
!======================================================================
!===                                                                ===
!===     MODYLAS --- MOlecular DYnamics simulation software         ===
!===                 ^^        ^^                                   ===
!===                 for LArge Systems                              ===
!===                     ^^    ^                                    ===
!===   ! SI unit is used if no explicit specification               ===
!===                                                                ===
!======================================================================
      program modylas
      use version
      use md_condition
      use md_fmm_domdiv_flg
      use md_forces
      use md_monitors
      use md_multiplestep
      use md_segment
      use mpivar
      use trj_mpi
      use shakerattleroll
#include "timing.h90"
      implicit none
      integer(4) :: i0,k0
      call mpistart()
      call init_openmp

      call parse_args()
      call opening()
      call parse_input()
      call initialize_application()
      call print_params()
c
c     prepare initial state for MD calculation
C
      call cell_edge()
      call calc_ia2c()
      call atom2cell() 
 
      if(totnconst > 0) call update_shake_local()
 
      call comm_direct_3()
      call apply_pbc()
 
      call md_calculate_forces_charmm_short(
     -            fshort ,eneshort)
      call md_calculate_forces_charmm_middle(
     -            fmiddle,enemiddle)
      call md_calculate_forces_charmm_long(
     -            flong  ,enelong)
      wk_p_energy=eneshort+enemiddle+enelong
!$omp parallel do default(shared)
!$omp& private(k0,i0)
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        wk_f(:,i0)=fshort(:,i0)+fmiddle(:,i0)*maxMTm
     &                         +flong(:,i0)  *maxMTm*maxMTl
      enddo
      enddo

c
c     MD main loop
c
      if (myrank == 0) write(*,'(/,a)') '**** start main loop'
      TIME_START(TM_MAIN_LOOP)
      do while (mdstep<md_condition__howmany_steps)

        call nve_integrate()

        TIME_START(TM_ENE_REDUCTION)
        call calc_hamiltonian_nve() 
        TIME_STOP(TM_ENE_REDUCTION)

        mdstep=mdstep+1

        TIME_START(TM_OUTPUT)
        call record_current_state()
        TIME_STOP(TM_OUTPUT)

      enddo
      TIME_STOP(TM_MAIN_LOOP)

      if (myrank == 0) call cleanup()

      if (myrank == 0) write(*,'(/,a)') '**** modylas normally ended!'

      call closing()

      call mpiend()
      stop
      end
c----------------------------------------------------------------------
      subroutine parse_args
      use mpivar
      use md_file
      implicit none
      character(LEN=1024) :: session
      integer(4) :: narg, iargc
#ifdef __GNUC__
      intrinsic iargc
#endif
      include 'mpif.h'
      integer(4) :: ierr

      if (myrank == 0) then
        narg = iargc()
        if (narg /= 1) then
          write(*, '(a)') 'usage: modylas_mini session-name'
        end if
        call getarg(1, session)
      end if
      call mpi_bcast(session,1024,mpi_character,mpiout,
     &               mpi_comm_world,ierr)
      session_name = trim(session)

      return
      end
