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
! MODULEs
c----------------------------------------------------------------------
      module version
      character(20) :: MODYLAS_version="1.0.1"
      end module
c----------------------------------------------------------------------
      module atommass
      real(8),allocatable :: mass(:),r_mass(:)
      real(8) :: totalmass
      end module
c----------------------------------------------------------------------
      module trj_org
      real(8),allocatable :: xyz(:,:)
      real(8),allocatable :: v(:,:)
      integer(4) :: n=0
      end module
c----------------------------------------------------------------------
      module trj_mpi
      integer(4) ::  trj_mpi_img_dims(1)
      ! real(8),allocatable :: wkxyz(:,:)[:]
      real(8), POINTER :: wkxyz(:,:) => null ()
      integer(8) :: wkxyz_desc

      real(8),allocatable :: wkv(:,:)
      ! integer(4),allocatable :: i2m(:), m2i(:)[:]
      integer(4),allocatable :: i2m(:)
      integer(4), POINTER :: m2i(:) => null ()
      integer(8) :: m2i_desc

      !integer(4),allocatable :: tag(:,:,:),na_per_cell(:,:,:)[:]
      integer(4),allocatable :: tag(:,:,:)
      integer(4), POINTER :: na_per_cell(:,:,:) => null ()
      integer(8) :: na_per_cell_desc

      integer(4) :: na1cell,na5cell,nadirect 
      integer(4) :: naline,narea 
      real(8),parameter :: na1cellmargin=2.00d0 ! ya  !!100% margin
      end module
c----------------------------------------------------------------------
      module unitcell
      real(8) :: cellx=0.0d0, celly=0.0d0, cellz=0.0d0
      real(8) :: cellxh=0.0d0, cellyh=0.0d0, cellzh=0.0d0
      real(8) :: alpha=90.0d0, beta=90.0d0, gamma=90.0d0
      real(8) :: cellvol=0.0d0
      end module
c----------------------------------------------------------------------
      module param
      implicit none
      integer(4) :: npara
      integer(4),allocatable :: paranum(:)
      end module
c----------------------------------------------------------------------
      module nhc
      integer(4),parameter :: mnhc=5,mnys=5,nnhc=5
      real(8) :: rss(mnhc),vss(mnhc)
      real(8) :: rssb(mnhc),vssb(mnhc)
      real(8) :: box(3,3),vboxg(3,3)
      end module
c----------------------------------------------------------------------
      module cutoffradius
      real(8) :: cutrad, cutrad2
      end module
c----------------------------------------------------------------------
      module shakerattleroll
      real(8),allocatable :: xyzstr(:,:)
      integer(4) :: totnconst=0
      integer(4) :: totnconstL,maxibn,l0max
      integer(4),allocatable :: ibseL(:),shkijL(:,:,:)
      real(8),allocatable :: rmassL(:,:,:)
      real(8),allocatable :: dij2L(:,:)
      integer(4), allocatable :: ShakeGroupLeader(:)
      integer(4), allocatable :: nconstraints(:)
      integer(4), allocatable :: atom1S(:,:), atom2S(:,:)
      real(8), allocatable :: slength(:,:)
      integer(4),allocatable :: kshake(:)
      real(8) :: shake_tolerance
      integer(4) :: maxshakecycle
      end module
c----------------------------------------------------------------------
      module pshake
      integer(4),parameter:: n_cnst_max=6
      integer(4),parameter:: n_atom_max=4
      integer(4):: n_type_ps
      integer(4),allocatable:: type_psL(:)
      integer(4),allocatable:: gn2ips(:)
      integer(4),allocatable:: couple_ps(:,:,:), n_atom_ps(:)
      real(8),allocatable:: a_0(:,:,:), a_0_sym(:,:,:)
      real(8),allocatable:: rdij2_ps(:,:), mass_ps(:,:)
      end module
!----------------------------------------------------------------------
      module pshake_init
      integer(4):: rngrp_ps
      real(8),allocatable::r_init_ps(:,:,:)
      integer(4),allocatable:: n_cnst_ps(:)
      end module
!----------------------------------------------------------------------
      module md_charmm_lj
      real(8),allocatable :: epsilon_sqrt(:), R_half(:)
      end module
c----------------------------------------------------------------------
      module md_void
      integer(4) :: nvoid
      integer(4),allocatable :: void_atom1(:,:),void_atom2(:,:)
      integer(4),allocatable :: void_n(:)
      end module
c----------------------------------------------------------------------
      module md_fmm
      integer(4) :: ndirect=2, ncell=32, nlevel=5
      integer(4) :: ndcellmargin=0
      integer(4) :: nimage=96
      integer(4),allocatable,dimension(:,:) :: lddir
      integer(4) :: nload,nchunk
      real(8),allocatable,dimension(:,:) :: fa,pre
      integer(4),allocatable :: ndseg_fmmn(:,:,:)
      integer(4),allocatable,dimension(:,:) :: m_wk
      integer(4) :: nmerge_fmm
      integer(4) :: max_nsegments_per_cell
      integer(4) :: nmax=4
      integer(4) :: mdg 
      complex(8),allocatable,dimension(:,:) :: premm,preml,prell
      complex(8),allocatable,dimension(:,:,:,:,:,:) :: shmm,shml,shll
      complex(8),allocatable,dimension(:) :: winput, woutput, wewald
      complex(8),allocatable,dimension(:,:) :: shew
      real(8) :: margin_mnpc=5.0d0
      integer(4) :: lgflg=1 ! 0,1,2,3, when switch local2global
      real(8)::sysdpl(3),wk_sysdpl(3)
      end module
c----------------------------------------------------------------------
      module md_const
      real(8),parameter :: md_PI=3.1415926535897932384626433832795029d0
      real(8),parameter :: md_PI_sqrt=1.77245385090551602729816748334d0
      real(8),parameter :: md_r_PI_sqrt=1.0d0/md_PI_sqrt
      real(8),parameter :: md_AVOGADRO=6.02214199d+23
      real(8),parameter :: md_ELEMENTARY_CHARGE=1.602176462D-19
      real(8),parameter :: md_ATOMIC_MASS_UNIT=1.0d-3/md_AVOGADRO
      real(8),parameter :: md_VACUUM_DIELECTRIC_CONSTANT=8.854187817D-12
      real(8),parameter :: md_BOLTZMANN=1.3806503D-23
      real(8),parameter :: md_DEGREE=md_PI/180.0d0
      real(8),parameter :: md_CALORIE=4.184d0
      real(8),parameter :: md_E_CONVERT=md_AVOGADRO*1.0d-3/md_CALORIE
      real(8),parameter :: md_F_CONVERT=md_AVOGADRO*1.0d-13/md_CALORIE
      real(8),parameter :: rad2deg=57.29577951308232088d0
      real(8),parameter :: deg2rad=1.0d0/rad2deg
      real(8),parameter :: md_QQ_4PiE=md_ELEMENTARY_CHARGE **2   
     &                                 * 0.25d0 / md_PI  
     &                                 /md_VACUUM_DIELECTRIC_CONSTANT 
      real(8),parameter :: rvkbolz=1d0/md_BOLTZMANN
      end module
c----------------------------------------------------------------------
      module md_condition
      integer(4) :: mdstep = 0
      real(8) :: dt=0.0d0
      integer(4) :: degree_of_freedom=0
      real(8) :: degree_of_freedom_inverse=0.0d0
      integer(4) :: md_condition__howmany_steps=0
      end module
c----------------------------------------------------------------------
      module md_coulomb
      real(8),allocatable :: chgv(:)
      end module
c----------------------------------------------------------------------
      module md_ewald
      logical :: ewald_sterm=.false. ! default
      end module
c----------------------------------------------------------------------
      module md_file
      character(LEN=1024) :: session_name
      integer(4),parameter :: f_trj=11,f_mntr=12 
      integer(4),parameter :: f_restart_bin=16
      integer(4),parameter :: f_mdff=17
      end module
c----------------------------------------------------------------------
      module md_forces
      real(8),allocatable :: f(:,:)
      real(8),allocatable :: wk_f(:,:) 
      real(8),allocatable :: w3_f(:,:,:)
      real(8),allocatable :: chgv_table(:,:),epsilon_sqrt_table(:,:),
     &                       R_half_table(:,:)
      end module
c----------------------------------------------------------------------
      module md_monitors
      real(8) :: hamiltonian
      real(8) :: p_energy, wk_p_energy
      real(8) :: k_energy
      real(8) :: t_energy
      real(8) :: temperature
      end module
c----------------------------------------------------------------------
      module md_segment
      integer(4),allocatable :: segtop(:)
      integer(4),allocatable :: lsegtop(:),lseg_natoms(:) ! meta-data
      end module
c----------------------------------------------------------------------
      module md_periodic
      integer(4) :: nsegments=-1
      integer(4),allocatable :: seg_natoms(:)
      real(8),allocatable :: seg_cx(:)
      real(8),allocatable :: seg_cy(:)
      real(8),allocatable :: seg_cz(:)
      real(8),allocatable :: wseg_cx(:)
      real(8),allocatable :: wseg_cy(:)
      real(8),allocatable :: wseg_cz(:)
      end module
c----------------------------------------------------------------------
      module g_main
      integer(4) :: trj_start=0,trj_interval=1
      integer(4) :: restart_start=0,restart_interval=1
      integer(4) :: mntr_start=0,mntr_interval=1
      real(8) :: maxwell_temperature=0.0d0
      integer(4) :: randomseed=1235
      logical :: reset_maxwell=.false.  ! default
      end module
c----------------------------------------------------------------------
      module md_fmm_domdiv_flg
      integer(4),allocatable :: idom(:,:,:)
      integer(4),allocatable :: idcell(:,:,:)
      integer(4),allocatable :: nd2c(:)
      integer(4),allocatable :: id2c(:,:)
      integer(4),allocatable :: ixflg(:),iyflg(:),izflg(:)
      integer(4),allocatable :: ixmax(:),iymax(:),izmax(:)
      integer(4),allocatable :: ixmin(:),iymin(:),izmin(:)
      integer(4) :: lxdiv,lydiv,lzdiv
      integer(4) :: nxdiv,nydiv,nzdiv
      integer(4) :: ndcell=1
      integer(4) :: nselfatm=0
      integer(4) :: ndatm=0
      integer(4) :: nselfseg=0
      integer(4) :: max_seg
      logical :: mpi_manual_division_flg=.false.

      complex(8),allocatable :: wwm_local0(:,:,:,:,:)

      type fmm_data_t
        integer(4) :: mcell_size, nscell
        integer(4) :: lclx, lcly, lclz
        integer(4) :: nscxdiv, nscydiv, nsczdiv
        integer(4) :: nbound_xm, nbound_ym, nbound_zm
        integer(4) :: nbound_xp, nbound_yp, nbound_zp
        complex(8),allocatable :: wm_local(:,:,:,:)
        complex(8),allocatable :: wm_global(:,:,:,:)
        complex(8),allocatable :: wl_local(:,:,:,:)
        complex(8),allocatable :: wwl_local(:,:,:,:,:)
      end type fmm_data_t
      type(fmm_data_t),allocatable,target :: fmm_data(:)

      end module
c----------------------------------------------------------------------
      module mod_wk_fmmewald
      complex(8),allocatable :: wk_wl(:)
      end module
c----------------------------------------------------------------------
      module comm_base
      integer npz,npy,npx
      integer ipz,ipy,ipx
      end module
c----------------------------------------------------------------------
      module comm_d3
      integer nczdiv, ncydiv, ncxdiv
      ! integer,allocatable :: icbufp(:)[:]
      integer, POINTER :: icbufp(:) => null ()
      ! integer,allocatable :: ircbufp(:)[:]
      integer, POINTER :: ircbufp(:) => null ()
      ! integer,allocatable :: icbufm(:)[:]
      integer, POINTER :: icbufm(:) => null ()
      ! integer,allocatable :: ircbufm(:)[:]
      integer, POINTER :: ircbufm(:) => null ()
      ! integer,allocatable :: ibuffp(:)[:]
      integer, POINTER :: ibuffp(:) => null ()
      ! integer,allocatable :: irbuffp(:)[:]
      integer, POINTER :: irbuffp(:) => null ()
      ! integer,allocatable :: ibuffm(:)[:]
      integer, POINTER :: ibuffm(:) => null ()
      ! integer,allocatable :: irbuffm(:)[:]
      integer, POINTER :: irbuffm(:) => null ()
      ! real(8),allocatable :: buffp(:,:)[:]
      real(8), POINTER :: buffp(:,:) => null ()
      ! real(8),allocatable :: rbuffp(:,:)[:]
      real(8), POINTER :: rbuffp(:,:) => null ()
      ! real(8),allocatable :: buffm(:,:)[:]
      real(8), POINTER :: buffm(:,:) => null ()
      ! real(8),allocatable :: rbuffm(:,:)[:]
      real(8), POINTER :: rbuffm(:,:) => null ()
      end module
c----------------------------------------------------------------------
      module comm_bd
      integer nczdiv, ncydiv, ncxdiv
      integer max_mvatom      ! max number of atoms moving to the cell.
      integer max_mvseg       ! max number of segments moving to the cell.
      integer max_cellcbd     ! max number of cells on communication buffer.
      real(8),allocatable :: abucket(:,:,:,:,:)
      integer,allocatable :: iabucket(:,:,:,:)
      integer,allocatable :: isbucket(:,:,:,:)
      integer,allocatable :: ncseg(:,:,:)
      integer,allocatable :: ncatom(:,:,:)
      real(8),allocatable :: buffp(:,:)
      real(8),allocatable :: buffm(:,:)
      integer,allocatable :: ibuffp(:)
      integer,allocatable :: ibuffm(:)
      integer,allocatable :: isbufp(:)
      integer,allocatable :: isbufm(:)
      ! real(8),allocatable :: rbuff_p(:,:)[:]
      real(8), POINTER :: rbuff_p(:,:) => null ()
      real(8) :: rbuff_p_desc
      ! real(8),allocatable :: rbuff_m(:,:)[:]
      real(8), POINTER :: rbuff_m(:,:) => null ()
      real(8) :: rbuff_m_desc
      ! integer,allocatable :: irbuff_p(:)[:]
      integer, POINTER :: irbuff_p(:) => null ()
      integer :: irbuff_p_desc
      ! integer,allocatable :: irbuff_m(:)[:]
      integer, POINTER :: irbuff_m(:) => null ()
      integer :: irbuff_m_desc
      ! integer,allocatable :: irsbuf_p(:)[:]
      integer, POINTER :: irsbuf_p(:) => null ()
      integer :: irsbuf_p_desc
      ! integer,allocatable :: irsbuf_m(:)[:]
      integer, POINTER :: irsbuf_m(:) => null ()

      integer(4), dimension(1) :: comm_bd_img_dims

      integer :: irsbuf_m_desc

      integer,allocatable :: ncatmw(:,:,:,:)
      end module
c----------------------------------------------------------------------
      module md_multiplestep
      integer(4) :: maxMTm=1
      integer(4) :: maxMTl=1
      real(8),allocatable :: fshort(:,:),fmiddle(:,:),flong(:,:)
      real(8) :: eneshort,enemiddle,enelong
      end module
!----------------------------------------------------------------------
      module mpivar
      implicit none
      integer(4) :: myrank=0, nprocs=1, mpiout=0
      end module
!----------------------------------------------------------------------
      module ompvar
      integer(4) :: nomp=1
      end module
!----------------------------------------------------------------------
! SUBROUTINEs
!----------------------------------------------------------------------
      subroutine init_openmp
      use ompvar
      implicit none
!$    include 'omp_lib.h'
!$    nomp = omp_get_max_threads()
      return
      end
c----------------------------------------------------------------------
      subroutine fmod_set_maxsegments
      use md_fmm
      use md_periodic
      use md_fmm_domdiv_flg
      implicit none
      max_nsegments_per_cell = 
     -       max(int(nsegments/ncell**3*margin_mnpc),60)
      return
      end
c----------------------------------------------------------------------
      subroutine fmod_set_ncell()
      use trj_mpi
      use md_fmm
      use md_fmm
      use md_fmm_domdiv_flg
      implicit none
      integer(4) :: ntmp
      nimage = ncell*3
      nlevel = int(log(dble(ncell))/log(2.0d0)+0.000001d00)
      ntmp = (nmax+1)*(nmax+1)
      allocate(winput(ntmp))
      allocate(woutput(ntmp))
      allocate(wewald(ntmp))
      allocate(premm(ntmp,ntmp))
      allocate(preml(ntmp,ntmp))
      allocate(prell(ntmp,ntmp))
      allocate(shmm(ntmp,ntmp,0:1,0:1,0:1,0:nlevel))
      allocate(shml(ntmp,ntmp,-5:5,-5:5,-5:5,0:nlevel))
      allocate(shll(ntmp,ntmp,0:1,0:1,0:1,0:nlevel))
      allocate(shew(ntmp,ntmp))
      winput = dcmplx(0.0d0, 0.0d0)
      woutput = dcmplx(0.0d0, 0.0d0)
      wewald = dcmplx(0.0d0, 0.0d0)
      premm = dcmplx(0.0d0, 0.0d0)
      preml = dcmplx(0.0d0, 0.0d0)
      prell = dcmplx(0.0d0, 0.0d0)
      shmm = dcmplx(0.0d0, 0.0d0)
      shml = dcmplx(0.0d0, 0.0d0)
      shll = dcmplx(0.0d0, 0.0d0)
      shew = dcmplx(0.0d0, 0.0d0)
      end
c----------------------------------------------------------------------
      subroutine abort_with_message_a(mesg)
      implicit none
      character(LEN=*) mesg
      write(0,'(a)') mesg
      call mpistop()
      end
