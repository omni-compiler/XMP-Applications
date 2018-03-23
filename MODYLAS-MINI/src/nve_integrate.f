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
c----------------------------------------------------------------------
      subroutine nve_integrate()
c----------------------------------------------------------------------
      use atommass
      use shakerattleroll
      use md_condition
      use md_const
      use md_forces
      use md_monitors
      use trj_mpi
      use md_fmm
      use md_fmm_domdiv_flg
      use param
      use md_multiplestep
      use md_segment
      use mpivar
      use unitcell
#include "timing.h90"
      implicit none
      integer(4) :: i0,ipar,k0
      integer(4) :: MTm, MTl
      real(8) :: dtL, dthL, scaleM, scaleL

      dtL =dt/maxMTm/maxMTl
      dthL=0.5d0 * dtL

      DO MTl=1,maxMTl   !! == Long-range force ==
      DO MTm=1,maxMTm   !! == Middle-range force ==

      if(totnconst > 0) then
!$omp parallel default(shared)
!$omp& private(k0,i0)
!$omp do
        do k0=1,nselfseg
        do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          xyzstr(:,i0) = wkxyz(:,i0)
        enddo ! i0
        enddo ! k0
!$omp end do
!$omp end parallel
      endif

!$omp parallel default(shared)
!$omp& private(k0,i0,ipar)
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ipar=paranum(m2i(i0))
        wkv(1,i0)=wkv(1,i0)+wk_f(1,i0)*dthL*r_mass(ipar)
        wkv(2,i0)=wkv(2,i0)+wk_f(2,i0)*dthL*r_mass(ipar)
        wkv(3,i0)=wkv(3,i0)+wk_f(3,i0)*dthL*r_mass(ipar)
        wkxyz(1,i0) = wkxyz(1,i0) + wkv(1,i0)*dtL
        wkxyz(2,i0) = wkxyz(2,i0) + wkv(2,i0)*dtL
        wkxyz(3,i0) = wkxyz(3,i0) + wkv(3,i0)*dtL
      enddo ! i0
      enddo ! k0
!$omp end do
!$omp end parallel

      if (totnconst > 0) then
        TIME_START(TM_SHAKE)
        call shake_roll(dtL)
        TIME_STOP(TM_SHAKE)
      endif

      if(MTl==maxMTl.and.MTm==maxMTm)then
        call update_wsegc()
        TIME_START(TM_MIGRATION)
        call comm_bound() 
        TIME_STOP(TM_MIGRATION)
        TIME_START(TM_SHAKE)
        if (totnconst > 0) call update_shake_local()
        TIME_STOP(TM_SHAKE)
      endif

      TIME_START(TM_COMM_DIRECT)
      call comm_direct_3()
      TIME_STOP(TM_COMM_DIRECT)
      call apply_pbc()
!!    ^^^ short ^^^
      call md_calculate_forces_charmm_short(
     -            fshort,eneshort)
!!    ^^^ middle ^^^
      if(MTm==maxMTm)then
        scaleM=1d0
        TIME_START(TM_ENERGY_DIRECT)
        call md_calculate_forces_charmm_middle(
     -              fmiddle,enemiddle)
        TIME_STOP(TM_ENERGY_DIRECT)
      else
        scaleM=0d0
      endif
!!    ^^^ long ^^^
      if(MTl==maxMTl.and.MTm==maxMTm)then
        scaleL=1d0
        TIME_START(TM_FMM)
        call md_calculate_forces_charmm_long(
     -              flong,enelong)
        TIME_STOP(TM_FMM)
      else
        scaleL=0d0
      endif

!!    ^^^ sum short, middle, and long ^^^
      wk_p_energy=eneshort+enemiddle+enelong

!$omp parallel default(shared)
!$omp& private(k0,i0)
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
         wk_f(:,i0)=fshort(:,i0)+scaleM*fmiddle(:,i0)*maxMTm
     &                          +scaleL*flong(:,i0)  *maxMTm*maxMTl
      enddo
      enddo
!$omp end do
!$omp end parallel

!$omp parallel default(shared)
!$omp& private(k0,i0,ipar)
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ipar=paranum(m2i(i0))
        wkv(1,i0) = wkv(1,i0) + wk_f(1,i0)*dthL*r_mass(ipar)
        wkv(2,i0) = wkv(2,i0) + wk_f(2,i0)*dthL*r_mass(ipar)
        wkv(3,i0) = wkv(3,i0) + wk_f(3,i0)*dthL*r_mass(ipar)
      enddo ! i0
      enddo ! k0
!$omp end do
!$omp end parallel

      if (totnconst .gt. 0) then
        TIME_START(TM_RATTLE)
        call rattle_roll(dtL)
        TIME_STOP(TM_RATTLE)
      endif

      ENDDO  !  MT middle
      ENDDO  !  MT long

      return
      end
c----------------------------------------------------------------------
c    
c     calculate thermodynamic values
c
c----------------------------------------------------------------------
      subroutine calc_hamiltonian_nve()
c----------------------------------------------------------------------
      use md_const
      use md_monitors
      use md_periodic
      use md_condition
      use mpivar
      use unitcell
      implicit none
      real(8) :: totke
      include 'mpif.h'
      integer(4) :: ierr

!     ^^^ reduce potential energy ^^^
!coarray      call mpi_allreduce(wk_p_energy,p_energy,1,
!coarray     &       mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      p_energy = wk_p_energy
      call co_sum(p_energy)
!!

!     ^^^ reduce kinetic energy ^^^
      call k_energy_scaler(totke)

!     ^^^ temperature ^^^
      k_energy    = totke
      temperature = 2.0d0*totke*degree_of_freedom_inverse*rvkbolz
      t_energy    = p_energy + k_energy
      hamiltonian = t_energy

      return
      end
