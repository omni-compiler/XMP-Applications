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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     kinetic energy calculation of particles
c
c     (kinetic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine k_energy_scaler(k_ene_sum)
      use atommass
      use trj_org
      use trj_mpi
      use md_const
      use md_fmm
      use md_fmm_domdiv_flg
      use param
      use mpivar
      use ompvar
      implicit none
      include 'mpif.h'
      integer(4) :: ierr
!$    include 'omp_lib.h'

      integer(4) :: ii,i0,ipar,iam
      integer(4) :: i,j,k
      integer(4) :: icx0,icy0,icz0,icxyz0
      real(8) :: wk_k_ene(0:nomp-1)
      real(8) :: k_ene_sum
      real(8) :: wk_ksum
      real(8) :: wmass

      iam = 0

!$omp parallel default(none)
!$omp& private(iam,ii,i0,j,k,ipar,wmass)
!$omp& private(icx0,icy0,icz0,icxyz0)
!$omp& shared(wkv,mass,wk_k_ene,m2i,paranum)
!$omp& shared(tag,na_per_cell,lxdiv,lydiv,lzdiv)
!$    iam = omp_get_thread_num()
      wk_k_ene(iam) = 0.0d0
!$omp do
      do ii=1,lxdiv*lydiv*lzdiv
      icz0=mod(ii-1,lzdiv)   +3   
      icy0=mod(ii-1,lzdiv*lydiv)
      icy0=icy0/lzdiv      +3   
      icx0=(ii-1)/(lzdiv*lydiv)+3
      do i0=tag(icz0,icy0,icx0),
     &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
        ipar=paranum(m2i(i0))
        wmass=mass(ipar)
        wk_k_ene(iam)=wk_k_ene(iam)+wmass*(wkv(1,i0)*wkv(1,i0)
     &                                    +wkv(2,i0)*wkv(2,i0)
     &                                    +wkv(3,i0)*wkv(3,i0))
      enddo ! i0
      enddo ! ii
!$omp end do
!$omp end parallel

      wk_ksum=0d0
      do iam=0,nomp-1 
      wk_ksum=wk_ksum+0.5d0*wk_k_ene(iam)
      enddo

!coarray      call mpi_allreduce(wk_ksum,k_ene_sum,1,
!coarray     &     mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      k_ene_sum = wk_ksum
      call co_sum(k_ene_sum)
!!

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
