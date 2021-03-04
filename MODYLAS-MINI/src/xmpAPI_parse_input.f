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
!------------------------------------------------------------------------
      subroutine parse_input()
      use mpivar
      use md_file
      implicit none

!--  mdconf
      call read_mdconf()
      if (myrank == 0) write(*,'(/,a)') 
     &  '**** '//trim(session_name) //'.mdconf read end successfully!'

!--  mdxyz.bin
      call read_mdxyzbin()
      if (myrank == 0) write(*,'(/,a)')
     &  '**** '//trim(session_name)//'.mdxyz.bin read end successfully!'

!--  mdff.bin
      call read_mdffbin()
      if (myrank == 0) write(*,'(/,a)')
     &  '**** '//trim(session_name)//'.mdff.bin read end successfully!'

      end
!------------------------------------------------------------------------
      subroutine read_mdconf
      use md_file
      use md_const
      use md_condition
      use md_fmm_domdiv_flg
      use md_multiplestep
      use nhc
      use param
      use md_segment
      use md_periodic
      use shakerattleroll
      use g_main
      use cutoffradius
      use md_fmm
      use md_ewald
      use ConfigRead
      implicit none

      call ConfigRead_parse(trim(session_name)//'.mdconf')

      trj_start = ConfigRead_get_int('trj_start')
      trj_interval = ConfigRead_get_int('trj_interval')
      restart_start = ConfigRead_get_int('restart_start')
      restart_interval = ConfigRead_get_int('restart_interval')
      mntr_start = ConfigRead_get_int('mntr_start')
      mntr_interval = ConfigRead_get_int('mntr_interval')
      randomseed = ConfigRead_get_int('randomseed')
      dt = ConfigRead_get_double('dt')
      md_condition__howmany_steps = ConfigRead_get_int('step')
      reset_maxwell = ConfigRead_get_bool('maxwell_velocities')
      maxwell_temperature  = ConfigRead_get_double('temperature')
      maxMTm = ConfigRead_get_int('nstep_skip_middle')
      maxMTl = ConfigRead_get_int('nstep_skip_long')
      mpi_manual_division_flg = ConfigRead_get_bool('manual_division')
      if (mpi_manual_division_flg) then
        nxdiv = ConfigRead_get_int('nxdiv')
        nydiv = ConfigRead_get_int('nydiv')
        nzdiv = ConfigRead_get_int('nzdiv')
      end if
      maxshakecycle = ConfigRead_get_int('shake_max_iteration')
      shake_tolerance = ConfigRead_get_double('shake_tolerance')
      cutrad = 1.0d-10 * ConfigRead_get_double('cutoff')
      cutrad2 = cutrad * cutrad
      ndirect = ConfigRead_get_int('ndirect')
      nmax = ConfigRead_get_int('nmax')
      lgflg = ConfigRead_get_int('ULswitch')
      ewald_sterm = ConfigRead_get_bool('ewald_surface_term')
      ncell = ConfigRead_get_int('ncell')

      call fmod_set_ncell()

      end
!------------------------------------------------------------------------
      subroutine read_mdxyzbin
      use trj_org
      use trj_mpi
      use nhc
      use md_file
      use g_main
      use mpivar
      use unitcell
      implicit none
      integer(4) :: io, n_nhc
!KF
      integer(4) :: ierr
      include 'mpif.h'
!KF end
      if(myrank==0) then
        open(f_restart_bin, file=trim(session_name)// '.mdxyz.bin', 
     &         iostat=io, status='old', 
     &         access='sequential',form='unformatted')
        if (io /= 0) then
          call abort_with_message_a('Cannot read mdxyz.bin file.')
        endif
        read(f_restart_bin,end=100) n
      endif
!coarray      call MPI_Bcast(n, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
      call co_broadcast(n,source_image=1)
!!

      if (.not. allocated(xyz)) then
         allocate(xyz(3,n))
         allocate(v(3,n))
         allocate(i2m(n))
      endif
      if(myrank.eq.0) then
c       Read coordinates and velocities of atoms.
c       Trajectory must be arranged in the same order as input.
        read(f_restart_bin) xyz(1:3,:), v(1:3,:)
c       Read positions and velocities of thermostats.
        read(f_restart_bin) n_nhc
        if (n_nhc .ne. nnhc)  call abort_with_message_a('nnhc')
        read(f_restart_bin) rss, vss
c       Read positions and velocities of barostats.
        read(f_restart_bin) n_nhc
        if (n_nhc .ne. nnhc)  call abort_with_message_a('nnhc')
        read(f_restart_bin) rssb, vssb
c       Read cell parameters (length and angles).
        read(f_restart_bin) cellx,celly,cellz, alpha,beta,gamma,vboxg
      endif ! myrank==0
  100 continue
      close(f_restart_bin)
!KF
!coarray      call MPI_Bcast(xyz  ,  3*n, MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(v    ,  3*n, MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(n_nhc,  1  , MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(rss  ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(vss  ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(rssb ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(vssb ,  5  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(cellx,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(celly,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(cellz,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(alpha,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(beta ,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(gamma,  1  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(vboxg,  9  , MPI_REAL8   , 0, MPI_COMM_WORLD, ierr)
      call co_broadcast(   xyz(1:3,1:n), source_image=1 )
      call co_broadcast(     v(1:3,1:n), source_image=1 )
      call co_broadcast(          n_nhc, source_image=1 )
      call co_broadcast(       rss(1:5), source_image=1 )
      call co_broadcast(       vss(1:5), source_image=1 )
      call co_broadcast(      rssb(1:5), source_image=1 )
      call co_broadcast(      vssb(1:5), source_image=1 )
      call co_broadcast(          cellx, source_image=1 )
      call co_broadcast(          celly, source_image=1 )
      call co_broadcast(          cellz, source_image=1 )
      call co_broadcast(          alpha, source_image=1 )
      call co_broadcast(           beta, source_image=1 )
      call co_broadcast(          gamma, source_image=1 )
      call co_broadcast( vboxg(1:3,1:3), source_image=1 )
!!
!KF end
      cellxh = 0.5d0 * cellx
      cellyh = 0.5d0 * celly
      cellzh = 0.5d0 * cellz
      cellvol = cellx*celly*cellx

      return
      end
!----------------------------------------------------------------------
      subroutine read_mdffbin
      use md_file
      use mpivar
      use trj_org
      use param
      use md_periodic
      use md_segment
      use atommass
      use md_coulomb
      use md_charmm_lj
      use shakerattleroll
      use pshake_init
      use md_void
      implicit none
      include 'mpif.h'
      integer(4) :: force_field
      integer(4) :: num
      integer(4),allocatable :: dummy(:)
      integer(4) :: i, j, k, ierr, io

      if(myrank==0) then
        open(f_mdff, file=trim(session_name)//'.mdff.bin', iostat=io,
     &       status='old', access='sequential', form='unformatted')
        if(io.ne.0) then
          call abort_with_message_a('Cannot open mdff.bin file.')
        endif
      endif

      if(myrank==0) then
        read(f_mdff) force_field
        if (force_field /= 100) then
           call abort_with_message_a('Non-CHARMM potential is used!')
        endif
      endif

!-- parameter_number
      call read_bcast_int(npara, 1, 1, f_mdff)
!*    if (myrank == 0) write(*,*) 'npara=', npara

      allocate(paranum(n))
      call read_bcast_int(paranum, 1, n, f_mdff)

!-- read_segment
      call read_bcast_int(nsegments, 1, 1, f_mdff)

      allocate(seg_natoms(nsegments),segtop(nsegments))
      call read_bcast_int(seg_natoms, 1, nsegments, f_mdff)
      call read_bcast_int(segtop, 1, nsegments, f_mdff)

      allocate(seg_cx(nsegments),seg_cy(nsegments),seg_cz(nsegments))

!-- molecule (not used)
      if (myrank==0) then
        read(f_mdff) num
        allocate(dummy(num))
        read(f_mdff) dummy
        read(f_mdff) dummy
        deallocate(dummy)
      end if
      
!-- mass
      allocate(mass(npara))
      allocate(r_mass(npara))
      call read_bcast_real(mass, 1, npara, f_mdff)

      r_mass(:) = 1.0d0 / mass(:)

!-- charge
      allocate(chgv(npara))
      call read_bcast_real(chgv, 1, npara, f_mdff)

!-- LJ
      allocate(epsilon_sqrt(npara),R_half(npara))
      call read_bcast_real(epsilon_sqrt, 1, npara, f_mdff)
      call read_bcast_real(R_half, 1, npara, f_mdff)

!-- shake
      call read_bcast_int(totnconst, 1, 1, f_mdff)
      if (totnconst==0)  then
        call abort_with_message_a('totnconst = 0 !')
      endif

      call read_bcast_int(rngrp_ps, 1, 1, f_mdff)

      allocate(ShakeGroupLeader(npara))
      call read_bcast_int(ShakeGroupLeader, npara, 1, f_mdff)

      allocate(nconstraints(rngrp_ps))
      call read_bcast_int(nconstraints, rngrp_ps, 1, f_mdff)

      allocate(atom1S(rngrp_ps,10), atom2S(rngrp_ps,10))
      allocate(slength(rngrp_ps,10))
      atom1S(:,:) = -1
      atom2S(:,:) = -1
      slength(:,:) = -1.0d0
      if(myrank==0)  then
        do i = 1, rngrp_ps
          do j = 1, nconstraints(i)
            read(f_mdff) atom1S(i,j), atom2S(i,j), slength(i,j)
          enddo
        enddo
      endif
!coarray      call MPI_Bcast(atom1S, rngrp_ps*10, MPI_INTEGER4, 0, 
!coarray     &               MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(atom2S, rngrp_ps*10, MPI_INTEGER4, 0, 
!coarray     &               MPI_COMM_WORLD, ierr)
!coarray      call MPI_Bcast(slength, rngrp_ps*10, MPI_REAL8, 0, 
!coarray     &               MPI_COMM_WORLD, ierr)
      call co_broadcast(  atom1S, source_image=1 )
      call co_broadcast(  atom2S, source_image=1 )
      call co_broadcast( slength, source_image=1 )
!!

!-- void
      call read_bcast_int(nvoid, 1, 1, f_mdff)
!*    if (myrank == 0) write(*,*) 'input,nvoid=',nvoid
      if(nvoid==0) return

      call read_bcast_int(num, 1, 1, f_mdff)

      allocate(void_n(npara))
      allocate(void_atom1(npara,num),void_atom2(npara,num))
      call read_bcast_int(void_n, 1, npara, f_mdff)
      call read_bcast_int(void_atom1, 1, npara*num, f_mdff)
      call read_bcast_int(void_atom2, 1, npara*num, f_mdff)

!--
      if(myrank==0) close(f_mdff)

      return
      end
!-------------------------------------------------------------------------
      subroutine read_bcast_int(data, n, m, in)
      use mpivar
      implicit none
      include 'mpif.h'
      integer(4) :: n, m, data(m, n), in
      integer(4) :: i, ierr
      do i = 1, n
        if (myrank == 0) read(in) data(1:m,i)
      end do
      call MPI_Bcast(data, n*m, MPI_INTEGER4, 0,  MPI_COMM_WORLD, ierr)
      end
!-------------------------------------------------------------------------
      subroutine read_bcast_real(data, n, m, in)
      use mpivar
      implicit none
      include 'mpif.h'
      integer(4) :: n, m, in
      real(8) :: data(m, n)
      integer(4) :: i, ierr
      do i = 1, n
        if (myrank == 0) read(in) data(1:m,i)
      end do
!coarray      call MPI_Bcast(data, n*m, MPI_REAL8, 0,  MPI_COMM_WORLD, ierr)
      call co_broadcast( data(1:m,1:n), source_image=1 )
!!
      end
