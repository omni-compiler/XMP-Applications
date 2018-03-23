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
c======================================================================
c===================   implementation of function   ===================
c======================================================================
      subroutine md_calculate_forces_charmm_short(ftmp,enetmp)
c======================================================================
      use trj_mpi
      use md_forces
      use md_monitors
      use md_fmm_domdiv_flg
      use md_segment
      use ompvar
      implicit none
      integer(4) :: i0,k0
      integer(4) :: iam
      real(8) :: ftmp(3,nadirect)
      real(8) :: enetmp
!$    include 'omp_lib.h'

      call init_comm_buffer()

!$omp parallel default(none)
!$omp& private(k0,i0,iam)
!$omp& shared(lsegtop,lseg_natoms,nselfseg)
!$omp& shared(nomp,ftmp,w3_f)
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ftmp(1:3,i0)=0d0
        do iam = 0,nomp-1
          ftmp(1,i0) = ftmp(1,i0) + w3_f(1,i0,iam)
          ftmp(2,i0) = ftmp(2,i0) + w3_f(2,i0,iam)
          ftmp(3,i0) = ftmp(3,i0) + w3_f(3,i0,iam)
        end do ! iam
      end do ! i0
      end do ! k0
!$omp end do
!$omp end parallel

      enetmp=wk_p_energy

      return
      end
c======================================================================
      subroutine md_calculate_forces_charmm_middle(ftmp,enetmp)
c======================================================================
      use trj_mpi
      use md_forces
      use md_monitors
      use md_fmm_domdiv_flg
      use md_segment
      use ompvar
      implicit none
      integer(4) :: i0,k0
      integer(4) :: iam
      real(8) :: ftmp(3,nadirect)
      real(8) :: enetmp
!$    include 'omp_lib.h'

      call init_comm_buffer()

      call energy_direct

      call remove_void123lj()
      call remove_void123cl()

!$omp parallel default(none)
!$omp& private(k0,i0,iam)
!$omp& shared(lsegtop,lseg_natoms,nselfseg)
!$omp& shared(nomp,ftmp,w3_f)
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ftmp(1:3,i0)=0d0
        do iam = 0,nomp-1
          ftmp(1,i0) = ftmp(1,i0) + w3_f(1,i0,iam)
          ftmp(2,i0) = ftmp(2,i0) + w3_f(2,i0,iam)
          ftmp(3,i0) = ftmp(3,i0) + w3_f(3,i0,iam)
        end do ! iam
      end do ! i0
      end do ! k0
!$omp end do
!$omp end parallel

      enetmp=wk_p_energy

      return
      end
c======================================================================
      subroutine md_calculate_forces_charmm_long(ftmp,enetmp)
c======================================================================
      use trj_mpi
      use md_forces
      use md_monitors
      use md_fmm_domdiv_flg
      use md_segment
      use md_ewald
      use ompvar
      implicit none
      integer(4) :: i,i0,k0
      integer(4) :: iam
      real(8) :: ftmp(3,nadirect)
      real(8) :: enetmp
!$    include 'omp_lib.h'

      call init_comm_buffer()

      call calc_mm
      call energy_fmm()
      if(ewald_sterm)then !! sterm is on (default in FMM)
        continue  
      else                !! sterm is off
        call calc_system_dipole
        call remove_ewaldsurfaceterm
      endif

!$omp parallel default(none)
!$omp& private(k0,i0,iam)
!$omp& shared(lsegtop,lseg_natoms,nselfseg)
!$omp& shared(nomp,ftmp,w3_f)
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ftmp(1:3,i0)=0d0
        do iam = 0,nomp-1
          ftmp(1,i0) = ftmp(1,i0) + w3_f(1,i0,iam)
          ftmp(2,i0) = ftmp(2,i0) + w3_f(2,i0,iam)
          ftmp(3,i0) = ftmp(3,i0) + w3_f(3,i0,iam)
        end do ! iam
      end do ! i0
      end do ! k0
!$omp end do
!$omp end parallel

      enetmp=wk_p_energy

      return
      end
