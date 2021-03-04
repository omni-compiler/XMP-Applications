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
c
      subroutine mpistart
      use mpivar
      implicit none
      include 'mpif.h'
      integer(4) :: ierr

!coarray      call mpi_init(ierr)
!coarray      call mpi_comm_size(mpi_comm_world,nprocs,ierr)
!coarray      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      nprocs = num_images()
      myrank = this_image()-1   
!!
      mpiout=0
      return 
      end

      subroutine mpiend
      implicit none
!coarray      include 'mpif.h'
!coarray      integer ierr
!coarray      call mpi_finalize(ierr)
      return
      end

      subroutine mpistop
      implicit none
      include 'mpif.h'
      integer ierr
      call mpi_abort(mpi_comm_world,ierr)
      return
      end
