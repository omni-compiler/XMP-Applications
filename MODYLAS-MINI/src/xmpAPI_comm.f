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
      subroutine init_comm_buffer()
c----------------------------------------------------------------------
      use trj_org        ! n
      use trj_mpi        ! n
      use md_forces    ! f, wk_k
      use md_monitors  ! p_energy, wk_p_energy
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use ompvar
      implicit none
      integer(4) :: ii,i0,iam
      integer(4) :: icx0,icy0,icz0
!$    include 'omp_lib.h'

      wk_p_energy = 0d0

!$omp parallel default(shared)
!$omp& private(iam,i0,ii)
!$omp& private(icx0,icy0,icz0)
!$omp do
      do iam = 0,nomp-1
      do ii=1,lxdiv*lydiv*lzdiv
      icz0=mod(ii-1,lzdiv)   +3   
      icy0=mod(ii-1,lzdiv*lydiv)
      icy0=icy0/lzdiv      +3   
      icx0=(ii-1)/(lzdiv*lydiv)+3
      do i0=tag(icz0,icy0,icx0),
     &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
        w3_f(1,i0,iam) = 0.0d0
        w3_f(2,i0,iam) = 0.0d0
        w3_f(3,i0,iam) = 0.0d0
      end do ! i0
      end do ! ii
      end do ! iam
!$omp end do
!$omp end parallel

      return
      end
c----------------------------------------------------------------------
      subroutine init_comm_bound() 
c----------------------------------------------------------------------
      use comm_base
      use comm_bd
      use trj_mpi
      use md_fmm_domdiv_flg
      use mpivar
      use xmp_api
      implicit none
      integer(4) :: itmp

      npz = nzdiv
      npy = nydiv
      npx = nxdiv
      ipz = izflg(myrank)-1
      ipy = iyflg(myrank)-1
      ipx = ixflg(myrank)-1
      nczdiv = lzdiv
      ncydiv = lydiv
      ncxdiv = lxdiv

      max_mvatom = na1cell*64
      max_mvseg = max_mvatom

      max_cellcbd = max(nczdiv*ncydiv, ncydiv*ncxdiv)
      max_cellcbd = max(max_cellcbd, nczdiv*ncxdiv)
      itmp = max(ncydiv+ncxdiv, ncxdiv+nczdiv)
      itmp = max(itmp, nczdiv+ncydiv)
      max_cellcbd = max_cellcbd + 2*itmp + 4

      allocate(abucket(6, max_mvatom, nczdiv+2, ncydiv+2, ncxdiv+2))
      allocate(  iabucket(max_mvatom, nczdiv+2, ncydiv+2, ncxdiv+2))
      allocate(  isbucket(max_mvseg,  nczdiv+2, ncydiv+2, ncxdiv+2))
      allocate(                 ncseg(nczdiv+2, ncydiv+2, ncxdiv+2))
      allocate(                ncatom(nczdiv+2, ncydiv+2, ncxdiv+2))

      allocate(buffp   (6,max_cellcbd*max_mvatom))
      buffp_lb(1) = 1
      buffp_lb(2) = 1
      buffp_ub(1) = 6
      buffp_ub(2) = max_cellcbd*max_mvatom
      call xmp_new_local_array(buffp_local_desc,8,2, 
     & buffp_lb,buffp_ub,loc(buffp))
      call xmp_new_array_section(buffp_local_sec,2)

      allocate(buffm   (6,max_cellcbd*max_mvatom))
      buffm_lb(1) = 1
      buffm_lb(2) = 1
      buffm_ub(1) = 6
      buffm_ub(2) = max_cellcbd*max_mvatom
      call xmp_new_local_array(buffm_local_desc,8,2, 
     & buffm_lb,buffm_ub,loc(buffm))
      call xmp_new_array_section(buffm_local_sec,2)

      allocate(ibuffp  (  max_cellcbd*max_mvatom))
      ibuffp_lb(1) = 1
      ibuffp_ub(1) = max_cellcbd*max_mvatom
      call xmp_new_local_array(ibuffp_local_desc,8,1,
     & ibuffp_lb,ibuffp_ub,loc(ibuffp))
      call xmp_new_array_section(ibuffp_local_sec,1)

      allocate(ibuffm  (  max_cellcbd*max_mvatom))

      allocate(isbufp  (2*max_cellcbd + 1 + max_cellcbd*max_mvseg))
      isbufp_lb(1) = 1
      isbufp_ub(1) = 2*max_cellcbd + 1 + max_cellcbd*max_mvseg
      call xmp_new_local_array(isbufp_local_desc,4,1, 
     & isbufp_lb,isbufp_ub,loc(isbufp))
      call xmp_new_array_section(isbufp_local_sec,1)

      allocate(isbufm  (2*max_cellcbd + 1 + max_cellcbd*max_mvseg))
      isbufm_lb(1) = 1
      isbufm_ub(1) = 2*max_cellcbd + 1 + max_cellcbd*max_mvseg
      call xmp_new_local_array(isbufm_local_desc,4,1,
     & isbufm_lb,isbufm_ub,loc(isbufm))
      call xmp_new_array_section(isbufm_local_sec,1)


      allocate( ncatmw(32, nczdiv+2, ncydiv+2, ncxdiv+2) )

!      !allocate(rbuff_p (6,max_cellcbd*max_mvatom)[*])
      rbuff_p_lb(1) = 1
      rbuff_p_lb(2) = 1
      rbuff_p_ub(1) = 6
      rbuff_p_ub(2) = max_cellcbd*max_mvatom
      call xmp_new_coarray(rbuff_p_desc,8,2,
     & rbuff_p_lb,rbuff_p_ub,1,img_dims)
      call xmp_coarray_bind(rbuff_p_desc,rbuff_p)
      call xmp_new_array_section(rbuff_p_sec,2)

!      !allocate(rbuff_m (6,max_cellcbd*max_mvatom)[*])
      rbuff_m_lb(1) = 1
      rbuff_m_lb(2) = 1
      rbuff_m_ub(1) = 6
      rbuff_m_ub(2) = max_cellcbd*max_mvatom
      call xmp_new_coarray(rbuff_m_desc,8,2,
     & rbuff_m_lb,rbuff_m_ub,1,img_dims)
      call xmp_coarray_bind(rbuff_m_desc,rbuff_m)
      call xmp_new_array_section(rbuff_m_sec,2)

!      !allocate(irbuff_p(  max_cellcbd*max_mvatom)[*])
      irbuff_p_lb(1) = 1
      irbuff_p_ub(1) = max_cellcbd*max_mvatom
      call xmp_new_coarray(irbuff_p_desc,4,1,
     & irbuff_p_lb,irbuff_p_ub,1,img_dims)
      call xmp_coarray_bind(irbuff_p_desc,irbuff_p)
      call xmp_new_array_section(irbuff_p_sec,1)

!      !allocate(irbuff_m(  max_cellcbd*max_mvatom)[*])
      irbuff_m_lb(1) = 1
      irbuff_m_ub(1) = max_cellcbd*max_mvatom
      call xmp_new_coarray(irbuff_m_desc,4,1,
     & irbuff_m_lb,irbuff_m_ub,1,img_dims)
      call xmp_coarray_bind(irbuff_m_desc,irbuff_m)
      call xmp_new_array_section(irbuff_m_sec,1)

!      !allocate(irsbuf_p(2*max_cellcbd + 1 + max_cellcbd*max_mvseg)[*])
       irsbuf_p_lb(1) = 1
       irsbuf_p_ub(1) = 2*max_cellcbd + 1 + max_cellcbd*max_mvseg
      call xmp_new_coarray(irsbuf_p_desc,4,1,
     & irsbuf_p_lb,irsbuf_p_ub,1,img_dims)
      call xmp_coarray_bind(irsbuf_p_desc,irsbuf_p)
      call xmp_new_array_section(irsbuf_p_sec,1)

!      !allocate(irsbuf_m(2*max_cellcbd + 1 + max_cellcbd*max_mvseg)[*])
       irsbuf_m_lb(1) = 1
       irsbuf_m_ub(1) = 2*max_cellcbd + 1 + max_cellcbd*max_mvseg
      call xmp_new_coarray(irsbuf_m_desc,4,1,
     & irsbuf_m_lb,irsbuf_m_ub,1,img_dims)
      call xmp_coarray_bind(irsbuf_m_desc,irsbuf_m)
      call xmp_new_array_section(irsbuf_m_sec,1)

      return
      end
c----------------------------------------------------------------------
      subroutine comm_bound()
c----------------------------------------------------------------------
      use trj_org
      use trj_mpi
      use md_periodic
      use md_segment
      use md_fmm
      use md_fmm_domdiv_flg
      use comm_base
      use comm_bd
      use param
      use mpivar
      use unitcell 
      use xmp_api
      implicit none
      INCLUDE 'mpif.h'
      integer nbase, nbase2
      real(8) rdcellz, rdcelly, rdcellx
      integer iczg0pr, icyg0pr, icxg0pr
      real(8) z0, y0, x0
      integer ipz_dest, ipy_dest, ipx_dest
      integer ipz_src, ipy_src, ipx_src
      integer nsa
      integer lcse(2,3,9)
      integer iczg, icyg, icxg
      integer icz, icy, icx
      integer icz0, icz1
      integer icy0, icy1
      integer icx0, icx1
      integer i,i0,k0,itmp
      integer ncc, ncs, ncsr, ics
      integer nca, ncar, ica, ncarp, ncarm
      integer nccp, ncsp, ncap
      integer nccm, ncsm, ncam
      integer ncc_p, ncs_p, nca_p
      integer ncc_m, ncs_m, nca_m
      integer ldcell
      integer loc_csbound
!
! ...  effective only when OpenMP compile.
!$    integer ncthread(2, 2, 2), itz, ity, itx, iam, nth
!$    integer omp_get_thread_num
!$    external omp_get_thread_num
!$    integer omp_get_num_threads
!$    external omp_get_num_threads
!
      integer istatus(mpi_status_size), ierr
!
      integer(4) status
      rdcellx=dble(ncell)/cellx
      rdcelly=dble(ncell)/celly
      rdcellz=dble(ncell)/cellz

! ---- 3D rank order rule. ----
!     ipx=mod(myrank,npx)
!     ipy=mod((myrank-ipx)/npx,npy)
!     ipz=mod((myrank-ipx-ipy*npx)/(npx*npy),npz)
!
!-----  bucket preparation. -----
!
! assumption ; system boundary condition is NOT applied to wseg_cx, 
!              wseg_cy, wseg_cz.
!
!
      ncseg = 0
      ncatom = 0
      ncatmw = 0

! cell index of the end of previous rank.
      iczg0pr = ncell * ipz / npz
      icyg0pr = ncell * ipy / npy
      icxg0pr = ncell * ipx / npx
!
!$omp parallel default(none)
!$omp& private(x0,y0,z0)
!$omp& private(icx,icy,icz)
!$omp& private(icxg,icyg,iczg)
!$omp& private(itx,ity,itz)
!$omp& private(ncs,nsa,nca)
!$omp& private(ics,ica)
!$omp& private(ncthread)
!$omp& private(iam)
!$omp& private(nth)
!$omp& shared(ncell)
!$omp& shared(nselfseg)
!$omp& shared(wseg_cx,wseg_cy,wseg_cz)
!$omp& shared(cellxh,cellyh,cellzh)
!$omp& shared(cellx,celly,cellz)
!$omp& shared(icxg0pr,icyg0pr,iczg0pr)
!$omp& shared(rdcellx,rdcelly,rdcellz)
!$omp& shared(ncxdiv,ncydiv,nczdiv)
!$omp& shared(isbucket)
!$omp& shared(iabucket)
!$omp& shared(abucket)
!$omp& shared(lsegtop,lseg_natoms)
!$omp& shared(ncseg)
!$omp& shared(wkxyz,wkv,m2i)
!$omp& shared(ncatmw)
!$omp& shared(myrank)
!$    iam = omp_get_thread_num()
!$    nth = omp_get_num_threads()
!$    ncthread = 0
!$    if(nth == 2) then
!$       ncthread(1,1,1) = 0
!$       ncthread(2,1,1) = 0
!$       ncthread(1,2,1) = 0
!$       ncthread(2,2,1) = 0
!$       ncthread(1,1,2) = 1
!$       ncthread(2,1,2) = 1
!$       ncthread(1,2,2) = 1
!$       ncthread(2,2,2) = 1
!$    elseif(nth == 4) then
! 1 to the 4th fractional subspace boundary.
! thread id for each subspace/domain. for FX-1.
!$       ncthread(1,1,1) = 0
!$       ncthread(2,1,1) = 0
!$       ncthread(1,2,1) = 1
!$       ncthread(2,2,1) = 1
!$       ncthread(1,1,2) = 2
!$       ncthread(2,1,2) = 2
!$       ncthread(1,2,2) = 3
!$       ncthread(2,2,2) = 3
!$    elseif(nth >= 8) then
! 1 to the 8th fractional subspace boundary.
! thread id for each subspace/domain. for K computer.
!$       ncthread(1,1,1) = 0
!$       ncthread(2,1,1) = 1
!$       ncthread(1,2,1) = 2
!$       ncthread(2,2,1) = 3
!$       ncthread(1,1,2) = 4
!$       ncthread(2,1,2) = 5
!$       ncthread(1,2,2) = 6
!$       ncthread(2,2,2) = 7
!$    endif
!
!... segment loop.
      do ics = 1, nselfseg
         x0 = wseg_cx(ics) + cellxh
         y0 = wseg_cy(ics) + cellyh
         z0 = wseg_cz(ics) + cellzh
!
! one extra cell outside process and beginning address 1 for cell index 
! is to be considered.
! segment address is negative in the negative direction system boundary.
! since intrinsic function INT truncates to the 0-direction nearest 
! whole number, addition of 1.0 before truncation is necessary.
! to avoid the addition intrinsic function FLOOR which truncate to the
! minus direction largest whole number is preferable.
! 
! global address.
         icxg = floor(x0*rdcellx) + 1
         icyg = floor(y0*rdcelly) + 1
         iczg = floor(z0*rdcellz) + 1
! global address to local address. one cell for boundary.
         icx = icxg - icxg0pr + 1
         icy = icyg - icyg0pr + 1
         icz = iczg - iczg0pr + 1
! sort cell for thread.
! ...  effective only when OpenMP compile.
!$       itx = 2
!$       if(icx <= ncxdiv/2+1) itx = 1
!$       ity = 2
!$       if(icy <= ncydiv/2+1) ity = 1
!$       itz = 2
!$       if(icz <= nczdiv/2+1) itz = 1

! ...  effective only when OpenMP compile.
!$       if(ncthread(itz,ity,itx) /= iam) cycle
!
! periodic bondary condition to coordinate wkxyz.
!
         if(icxg  <= 0) then
            do ica = 1, lseg_natoms(ics)
               nsa = lsegtop(ics) + ica - 1
               wkxyz(1,nsa) = wkxyz(1,nsa) + cellx
            end do
         endif
         if(icxg  > ncell) then
            do ica = 1, lseg_natoms(ics)
               nsa = lsegtop(ics) + ica - 1
               wkxyz(1,nsa) = wkxyz(1,nsa) - cellx
            end do
         endif
         if(icyg  <= 0) then
            do ica = 1, lseg_natoms(ics)
               nsa = lsegtop(ics) + ica - 1
               wkxyz(2,nsa) = wkxyz(2,nsa) + celly
            end do
         endif
         if(icyg  > ncell) then
            do ica = 1, lseg_natoms(ics)
               nsa = lsegtop(ics) + ica - 1
               wkxyz(2,nsa) = wkxyz(2,nsa) - celly
            end do
         endif
         if(iczg  <= 0) then
            do ica = 1, lseg_natoms(ics)
               nsa = lsegtop(ics) + ica - 1
               wkxyz(3,nsa) = wkxyz(3,nsa) + cellz
            end do
         endif
         if(iczg  > ncell) then
            do ica = 1, lseg_natoms(ics)
               nsa = lsegtop(ics) + ica - 1
               wkxyz(3,nsa) = wkxyz(3,nsa) - cellz
            end do
         endif
!
! segment data.
         ncs = ncatmw(2,icz,icy,icx) + 1
         isbucket(ncs,icz,icy,icx) = lseg_natoms(ics)
         ncatmw(2,icz,icy,icx) = ncs
! atom data.
!     In the bounding process, all of the atom belonging to the same 
!     segment must be treated in the same way.
         nca = ncatmw(1,icz,icy,icx)
         do ica = 1, lseg_natoms(ics)
            nsa = lsegtop(ics) + ica - 1
            nca = nca + 1
            abucket(1,nca,icz,icy,icx) = wkxyz(1,nsa)
            abucket(2,nca,icz,icy,icx) = wkxyz(2,nsa)
            abucket(3,nca,icz,icy,icx) = wkxyz(3,nsa)
            abucket(4,nca,icz,icy,icx) = wkv(1,nsa)
            abucket(5,nca,icz,icy,icx) = wkv(2,nsa)
            abucket(6,nca,icz,icy,icx) = wkv(3,nsa)
            iabucket(nca,icz,icy,icx) = m2i(nsa)
         end do
         ncatmw(1,icz,icy,icx) = nca
      end do
!$OMP ENDPARALLEL
      do icx = 1, ncxdiv + 2
        do icy = 1, ncydiv + 2
          do icz = 1, nczdiv + 2
            ncatom(icz,icy,icx) = ncatmw(1,icz,icy,icx)
            ncseg(icz,icy,icx) = ncatmw(2,icz,icy,icx)
          end do
        end do
      end do
!
!-----  boundary communication code starts here. ------
!
! ( bound +Z )
! nearest neighbor.
      icz = nczdiv + 2
      icy0 = 2
      icy1 = ncydiv + 1
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(1)
      icz = nczdiv + 2
      icy = 1
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(2)
      icz = nczdiv + 2
      icy = ncydiv + 2
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(3)
      icz = nczdiv + 2
      icy0 = 2
      icy1 = ncydiv + 1
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 2d-diagonal adjacent(4)
      icz = nczdiv + 2
      icy0 = 2
      icy1 = ncydiv + 1
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(5)
      icz = nczdiv + 2
      icy = 1
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(6)
      icz = nczdiv + 2
      icy = ncydiv + 2
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(7)
      icz = nczdiv + 2
      icy = 1
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(8)
      icz = nczdiv + 2
      icy = ncydiv + 2
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
      ncc = 0
      ncs = loc_csbound
      nca = 0
!
      call add_buffer(abucket,iabucket,isbucket,buffp,ibuffp,isbufp,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
!
      ipz_dest = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
      ipz_src  = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
!
      call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER,
     &             ipz_dest, myrank,
     &             ncsr, 1, MPI_INTEGER,
     &             ipz_src, ipz_src,
     &             mpi_comm_world, istatus, ierr )

      isbufp(ncs+1) = nca
!coarray      call mpi_sendrecv(isbufp, ncs + 1, MPI_INTEGER,
!coarray     &             ipz_dest, myrank,
!coarray     &             irsbuf_p, ncsr, MPI_INTEGER,
!coarray     &             ipz_src, ipz_src,
!coarray     &             mpi_comm_world, istatus, ierr )

      img_dims(1) = ipz_dest+1
      ! irsbuf_p(1:ncs+1)[ipz_dest+1] = isbufp(1:ncs+1) ! Put
      call xmp_array_section_set_triplet(isbufp_local_sec, 
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)
      call xmp_array_section_set_triplet(irsbuf_p_sec, 
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)
      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,irsbuf_p_desc,irsbuf_p_sec, 
     & isbufp_local_desc,isbufp_local_sec,status)

      !sync all
      call xmp_sync_all(status)
!!

      ncar = irsbuf_p(ncsr)
!coarray      call mpi_sendrecv(buffp, 6*nca, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_dest, myrank,
!coarray     &             rbuff_p, 6*ncar, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_src, ipz_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !rbuff_p(:,1:nca)[ipz_dest+1] = buffp(:,1:nca) ! Put
      call xmp_array_section_set_triplet(rbuff_p_sec,1,int(1,kind=8),
     & int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuff_p_sec,2,int(1,kind=8),
     & int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(buffp_local_sec, 
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(buffp_local_desc, 
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,rbuff_p_desc,rbuff_p_sec,
     & buffp_local_desc,buffp_local_sec,status)
!      sync all
      call xmp_sync_all(status)
!!

!coarray      call mpi_sendrecv(ibuffp, nca, MPI_INTEGER, 
!coarray     &             ipz_dest, myrank,
!coarray     &             irbuff_p, ncar, MPI_INTEGER, 
!coarray     &             ipz_src, ipz_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irbuff_p(1:nca)[ipz_dest+1] = ibuffp(1:nca) ! Put
      call xmp_array_section_set_triplet(irbuff_p_sec,1,int(1,kind=8),
     & int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(ibuffp_local_sec, 
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,irbuff_p_desc,irbuff_p_sec, 
     & ibuffp_local_desc,ibuffp_local_sec,status)


!      sync all
      call xmp_sync_all(status)

!!

!
! ( bound -Z )
!
! ( bound -Z )
! nearest neighbor.
      icz = 1
      icy0 = 2
      icy1 = ncydiv + 1
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(13)
      icz = 1
      icy = 1
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(14)
      icz = 1
      icy = ncydiv + 2
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(15)
      icz = 1
      icy0 = 2
      icy1 = ncydiv + 1
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 2d-diagonal adjacent(16)
      icz = 1
      icy0 = 2
      icy1 = ncydiv + 1
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(17)
      icz = 1
      icy = 1
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(18)
      icz = 1
      icy = ncydiv + 2
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(19)
      icz = 1
      icy = 1
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 3d-diagonal adjacent(20)
      icz = 1
      icy = ncydiv + 2
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
      ncc = 0
      ncs = loc_csbound
      nca = 0
!
      call add_buffer(abucket,iabucket,isbucket,buffm,ibuffm,isbufm,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
!
      ipz_dest = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
      ipz_src  = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
!
      call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER,
     &             ipz_dest, myrank,
     &             ncsr, 1, MPI_INTEGER,
     &             ipz_src, ipz_src,
     &             mpi_comm_world, istatus, ierr )

      isbufm(ncs+1) = nca
!coarray      call mpi_sendrecv(isbufm, ncs + 1, MPI_INTEGER,
!coarray     &             ipz_dest, myrank,
!coarray     &             irsbuf_m, ncsr, MPI_INTEGER,
!coarray     &             ipz_src, ipz_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irsbuf_m(1:ncs+1)[ipz_dest+1] = isbufm(1:ncs+1) ! Put

      call xmp_array_section_set_triplet(irsbuf_m_sec,1,int(1,kind=8),
     & int(ncs+1,kind=8),1,status)

      call xmp_array_section_set_triplet(isbufm_local_sec, 
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)

      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,irsbuf_m_desc,irsbuf_m_sec, 
     & isbufm_local_desc,isbufm_local_sec,status)



!      sync all
      call xmp_sync_all(status)
!!

      ncar = irsbuf_m(ncsr)
!coarray      call mpi_sendrecv(buffm, 6*nca, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_dest, myrank,
!coarray     &             rbuff_m, 6*ncar, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_src, ipz_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !rbuff_m(1:6,1:nca)[ipz_dest+1] = buffm(1:6,1:nca) ! Put

      call xmp_array_section_set_triplet(rbuff_m_sec,
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuff_m_sec,
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(buffm_local_sec, 
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(buffm_local_sec, 
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,rbuff_m_desc,rbuff_m_sec, 
     & buffm_local_desc,buffm_local_sec,status)


!      sync all
      call xmp_sync_all(status)
!!

!coarray      call mpi_sendrecv(ibuffm, nca, MPI_INTEGER, 
!coarray     &             ipz_dest, myrank,
!coarray     &             irbuff_m, ncar, MPI_INTEGER, 
!coarray     &             ipz_src, ipz_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irbuff_m(1:nca)[ipz_dest+1] = ibuffm(1:nca) ! Put

      call xmp_array_section_set_triplet(irbuff_m_sec,
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(ibuffm_local_sec, 
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,irbuff_m_desc,irbuff_m_sec, 
     & ibuffm_local_desc,ibuffm_local_sec,status)


!      sync all
      call xmp_sync_all(status)
!!

!
! merge source nearest neighbors(+Z) on receive buffer to local bucket.
      icx0 = 2
      icx1 = ncxdiv + 1
      icy0 = 2
      icy1 = ncydiv + 1
      icz = 2
      ldcell = 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
!
!  Merge part of 2d-diagonal adjacent on z-direction receive buffer
!  to cell bucket corresponding to y-directional nearest neighbor 
!  communication.
! source rank +Z 2d-diagonal adjacent(1)
      icx0 = 2
      icx1 = ncxdiv + 1
      icy = 1
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! source rank +Z 2d-diagonal adjacent(2)
      icx0 = 2
      icx1 = ncxdiv + 1
      icy = ncydiv + 2
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! source rank +Z 2d-diagonal adjacent(3)
      icx = 1
      icy0 = 2
      icy1 = ncydiv + 1
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank +Z 2d-diagonal adjacent(4)
      icx = ncxdiv+ 2
      icy0 = 2
      icy1 = ncydiv + 1
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
      ncc_p = 0
      ncs_p = loc_csbound
      nca_p = 0
!
      call add_bucket(abucket,iabucket,isbucket,
     &                rbuff_p,irbuff_p,irsbuf_p,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                ncc_p,ncs_p,nca_p)
!
! merge source rank nearest neighbors(-Z) on receive buffer to local bucket.
      icx0 = 2
      icx1 = ncxdiv + 1
      icy0 = 2
      icy1 = ncydiv + 1
      icz = nczdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
!
!  Merge part of 2d-diagonal adjacent on z-direction receive buffer
!  to cell bucket corresponding to y-directional nearest neighbor 
!  communication.
! source rank -Z 2d-diagonal adjacent(13)
      icx0 = 2
      icx1 = ncxdiv + 1
      icy = 1
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! source rank -Z 2d-diagonal adjacent(14)
      icx0 = 2
      icx1 = ncxdiv + 1
      icy = ncydiv + 2
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! source rank -Z 2d-diagonal adjacent(15)
      icx = 1
      icy0 = 2
      icy1 = ncydiv + 1
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank -Z 2d-diagonal adjacent(16)
      icx = ncxdiv + 2
      icy0 = 2
      icy1 = ncydiv + 1
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(ncydiv*ncxdiv + 2*ncxdiv + 2*ncydiv + 4) + 1
      ncc_m = 0
      ncs_m = loc_csbound
      nca_m = 0
!
      call add_bucket(abucket,iabucket,isbucket,
     &                rbuff_m,irbuff_m,irsbuf_m,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                ncc_m,ncs_m,nca_m)

!
! ( bound +Y & bound -Y )
! ( +Y buffer )
! nearest neighbor.
      icz0 = 2
      icz1 = nczdiv + 1
      icy = ncydiv + 2
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(10)
      icz0 = 2
      icz1 = nczdiv + 1
      icy = ncydiv + 2
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 2d-diagonal adjacent(12)
      icz0 = 2
      icz1 = nczdiv + 1
      icy = ncydiv + 2
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(nczdiv*ncxdiv + 2*nczdiv + 4) + 1
      nccp = 0
      ncsp = loc_csbound
      ncap = 0
!
      call add_buffer(abucket,iabucket,isbucket,buffp,ibuffp,isbufp,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                nccp,ncsp,ncap)
!
! ( -Y buffer )
! nearest neighbor.
      icz0 = 2
      icz1 = nczdiv + 1
      icy = 1
      icx0 = 2
      icx1 = ncxdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! 2d-diagonal adjacent(9)
      icz0 = 2
      icz1 = nczdiv + 1
      icy = 1
      icx = 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! 2d-diagonal adjacent(11)
      icz0 = 2
      icz1 = nczdiv + 1
      icy = 1
      icx = ncxdiv + 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(nczdiv*ncxdiv +  2*nczdiv + 4) + 1
      nccm = 0
      ncsm = loc_csbound
      ncam = 0
!
      call add_buffer(abucket,iabucket,isbucket,buffm,ibuffm,isbufm,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,nccm,ncsm,ncam)
!
! ( -Y buffer ) +Z source rank 3d-diagonal adjacent(5)
!
      ldcell = 1
      call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffm,ibuffm,isbufm,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_p,ncs_p,nca_p,nccm,ncsm,ncam,ldcell)
!
! ( +Y buffer ) +Z source rank 3d-diagonal adjacent(6)
!
      ldcell = 1
      call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffp,ibuffp,isbufp,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_p,ncs_p,nca_p,nccp,ncsp,ncap,ldcell)
!
! ( -Y buffer ) +Z source rank 3d-diagonal adjacent(7)
!
      ldcell = 1
      call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffm,ibuffm,isbufm,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_p,ncs_p,nca_p,nccm,ncsm,ncam,ldcell)
!
! ( +Y buffer ) +Z source rank 3d-diagonal adjacent(8)
!
      ldcell = 1
      call add_buffb(rbuff_p,irbuff_p,irsbuf_p,buffp,ibuffp,isbufp,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_p,ncs_p,nca_p,nccp,ncsp,ncap,ldcell)
!
! ( -Y buffer ) -Z source rank 3d-diagonal adjacent(17)
!
      ldcell = 1
      call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffm,ibuffm,isbufm,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_m,ncs_m,nca_m,nccm,ncsm,ncam,ldcell)
!
! ( +Y buffer ) -Z source rank 3d-diagonal adjacent(18)
!
      ldcell = 1
      call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffp,ibuffp,isbufp,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_m,ncs_m,nca_m,nccp,ncsp,ncap,ldcell)
!
! ( -Y buffer ) -Z source rank 3d-diagonal adjacent(19)
!
      ldcell = 1
      call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffm,ibuffm,isbufm,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_m,ncs_m,nca_m,nccm,ncsm,ncam,ldcell)
!
! ( +Y buffer ) -Z source rank 3d-diagonal adjacent(20)
!
      ldcell = 1
      call add_buffb(rbuff_m,irbuff_m,irsbuf_m,buffp,ibuffp,isbufp,
     &               max_mvatom,max_mvseg,max_cellcbd,
     &               ncc_m,ncs_m,nca_m,nccp,ncsp,ncap,ldcell)
!
!
! +Y comm.
      ipy_dest  = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
      ipy_src   = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx

!
      call mpi_sendrecv(ncsp + 1, 1, MPI_INTEGER,
     &             ipy_dest, myrank,
     &             ncsr, 1, MPI_INTEGER,
     &             ipy_src, ipy_src,
     &             mpi_comm_world, istatus, ierr )

      isbufp(ncsp+1) = ncap
!coarray      call mpi_sendrecv(isbufp, ncsp + 1, MPI_INTEGER,
!coarray     &             ipy_dest, myrank,
!coarray     &             irsbuf_p, ncsr, MPI_INTEGER,
!coarray     &             ipy_src, ipy_src,
!coarray     &             mpi_comm_world, istatus, ierr )
!      sync all ! do Not change
      call xmp_sync_all(status)
      !irsbuf_p(1:ncsp+1)[ipy_dest+1] = isbufp(1:ncsp+1) ! Put

      call xmp_array_section_set_triplet(irsbuf_p_sec,
     & 1,int(1,kind=8),int(ncsp+1,kind=8),1,status)

      call xmp_array_section_set_triplet(isbufp_local_sec, 
     & 1,int(1,kind=8),int(ncsp+1,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put_local(img_dims,irsbuf_p_desc,irsbuf_p_sec, 
     & isbufp_local_desc,isbufp_local_sec,status)



!      sync all
      call xmp_sync_all(status)
!!

      ncarp = irsbuf_p(ncsr)
!coarray      call mpi_sendrecv(buffp, 6*ncap, MPI_DOUBLE_PRECISION, 
!coarray     &             ipy_dest, myrank,
!coarray     &             rbuff_p, 6*ncarp, MPI_DOUBLE_PRECISION, 
!coarray     &             ipy_src, ipy_src,
!coarray     &             mpi_comm_world, istatus, ierr )

      !rbuff_p(1:6,1:ncap)[ipy_dest+1] = buffp(1:6,1:ncap) ! Put

      call xmp_array_section_set_triplet(rbuff_p_sec,
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuff_p_sec,
     & 2,int(1,kind=8),int(ncap,kind=8),1,status)

      call xmp_array_section_set_triplet(buffp_local_sec, 
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(buffp_local_sec, 
     & 2,int(1,kind=8),int(ncap,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put_local(img_dims,rbuff_p_desc,rbuff_p_sec, 
     & buffp_local_desc,buffp_local_sec,status)
!!

!coarray      call mpi_sendrecv(ibuffp, ncap, MPI_INTEGER, 
!coarray     &             ipy_dest, myrank,
!coarray     &             irbuff_p, ncarp, MPI_INTEGER, 
!coarray     &             ipy_src, ipy_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irbuff_p(1:ncap)[ipy_dest+1] = ibuffp(1:ncap) ! Put

      call xmp_array_section_set_triplet(irbuff_p_sec,
     & 1,int(1,kind=8),int(ncap,kind=8),1,status)

      call xmp_array_section_set_triplet(ibuffp_local_sec, 
     & 1,int(1,kind=8),int(ncap,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put_local(img_dims,irbuff_p_desc,irbuff_p_sec, 
     & ibuffp_local_desc,ibuffp_local_sec,status)


!      sync all
      call xmp_sync_all(status)
!!

!
! -Y comm.
!
      ipy_dest  = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
      ipy_src   = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
!
      call mpi_sendrecv(ncsm + 1, 1, MPI_INTEGER,
     &             ipy_dest, myrank,
     &             ncsr, 1, MPI_INTEGER,
     &             ipy_src, ipy_src,
     &             mpi_comm_world, istatus, ierr )

      isbufm(ncsm+1) = ncam
!coarray      call mpi_sendrecv(isbufm, ncsm + 1, MPI_INTEGER,
!coarray     &             ipy_dest, myrank,
!coarray     &             irsbuf_m, ncsr, MPI_INTEGER,
!coarray     &             ipy_src, ipy_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irsbuf_m(1:ncsm+1)[ipy_dest+1] = isbufm(1:ncsm+1) ! Put

      call xmp_array_section_set_triplet(irsbuf_m_sec,
     & 1,int(1,kind=8),int(ncsm+1,kind=8),1,status)

      call xmp_array_section_set_triplet(isbufm_local_sec, 
     & 1,int(1,kind=8),int(ncsm+1,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put_local(img_dims,irsbuf_m_desc,irsbuf_m_sec, 
     & isbufm_local_desc,isbufm_local_sec,status)

!      sync all
      call xmp_sync_all(status)
!!

      ncarm = irsbuf_m(ncsr)
!coarray      call mpi_sendrecv(buffm, 6*ncam, MPI_DOUBLE_PRECISION, 
!coarray     &             ipy_dest, myrank,
!coarray     &             rbuff_m, 6*ncarm, MPI_DOUBLE_PRECISION, 
!coarray     &             ipy_src, ipy_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !rbuff_m(1:6,1:ncam)[ipy_dest+1] = buffm(1:6,1:ncam) ! Put

      call xmp_array_section_set_triplet(rbuff_m_sec,
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuff_m_sec,
     & 2,int(1,kind=8),int(ncam,kind=8),1,status)

      call xmp_array_section_set_triplet(buffm_local_sec, 
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(buffm_local_sec, 
     & 2,int(1,kind=8),int(ncam,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put_local(img_dims,rbuff_m_desc,rbuff_m_sec, 
     & buffm_local_desc,buffm_local_sec,status)

!!

!coarray      call mpi_sendrecv(ibuffm, ncam, MPI_INTEGER, 
!coarray     &             ipy_dest, myrank,
!coarray     &             irbuff_m, ncarm, MPI_INTEGER, 
!coarray     &             ipy_src, ipy_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irbuff_m(1:ncam)[ipy_dest+1] = ibuffm(1:ncam) ! Put

      call xmp_array_section_set_triplet(irbuff_m_sec,
     & 1,int(1,kind=8),int(ncam,kind=8),1,status)

      call xmp_array_section_set_triplet(ibuffm_local_sec, 
     & 1,int(1,kind=8),int(ncam,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put_local(img_dims,irbuff_m_desc,irbuff_m_sec, 
     & ibuffm_local_desc,ibuffm_local_sec,status)
!      sync all
      call xmp_sync_all(status)
!!

!
! +Y receive buffer.
! merge source rank nearest neighbors(+Y) on receive buffer to local bucket.
      icx0 = 2
      icx1 = ncxdiv + 1
      icy = 2
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! source rank +Y 2d-diagonal adjacent(10)
      icx = 1
      icy = 2
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank +Y 2d-diagonal adjacent(12)
      icx = ncxdiv+ 2
      icy = 2
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank +Z 3d-diagonal adjacent(6)
      icx = 1
      icy = 2
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank +Z 3d-diagonal adjacent(8)
      icx = ncxdiv + 2
      icy = 2
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank -Z 3d-diagonal adjacent(18)
      icx = 1
      icy = 2
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank -Z 3d-diagonal adjacent(20)
      icx = ncxdiv + 2
      icy = 2
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(nczdiv*ncxdiv + 2*nczdiv + 4) + 1
      ncc_p = 0
      ncs_p = loc_csbound
      nca_p = 0
!
      call add_bucket(abucket,iabucket,isbucket,
     &                rbuff_p,irbuff_p,irsbuf_p,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                ncc_p,ncs_p,nca_p)
!
! -Y receive buffer.
! merge source rank nearest neighbors(-Y) on receive buffer to local bucket.
      icx0 = 2
      icx1 = ncxdiv + 1
      icy = ncydiv + 1
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx0
      lcse(2,3,ldcell) = icx1
! source rank -Y 2d-diagonal adjacent(9)
      icx = 1
      icy = ncydiv + 1
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank -Y 2d-diagonal adjacent(11)
      icx = ncxdiv+ 2
      icy = ncydiv + 1
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank +Z 3d-diagonal adjacent(5)
      icx = 1
      icy = ncydiv + 1
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank +Z 3d-diagonal adjacent(7)
      icx = ncxdiv + 2
      icy = ncydiv + 1
      icz = 2
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank -Z 3d-diagonal adjacent(17)
      icx = 1
      icy = ncydiv + 1
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
! source rank -Z 3d-diagonal adjacent(19)
      icx = ncxdiv + 2
      icy = ncydiv + 1
      icz = nczdiv + 1
      ldcell = ldcell + 1
      lcse(1,1,ldcell) = icz
      lcse(2,1,ldcell) = icz
      lcse(1,2,ldcell) = icy
      lcse(2,2,ldcell) = icy
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*(nczdiv*ncxdiv + 2*nczdiv + 4) + 1
      ncc_m = 0
      ncs_m = loc_csbound
      nca_m = 0
!
      call add_bucket(abucket,iabucket,isbucket,
     &                rbuff_m,irbuff_m,irsbuf_m,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                ncc_m,ncs_m,nca_m)

!
! ( bound +X )
!
! nearest neighbor.
      icz0 = 2
      icz1 = nczdiv + 1
      icy0 = 2
      icy1 = ncydiv + 1
      icx = ncxdiv + 2
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*nczdiv*ncydiv + 1
      ncc = 0
      ncs = loc_csbound
      nca = 0
!
      call add_buffer(abucket,iabucket,isbucket,buffp,ibuffp,isbufp,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
!
      ipx_dest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
      ipx_src   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
!
      call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER,
     &             ipx_dest, myrank,
     &             ncsr, 1, MPI_INTEGER,
     &             ipx_src, ipx_src,
     &             mpi_comm_world, istatus, ierr )

      isbufp(ncs+1) = nca
!coarray      call mpi_sendrecv(isbufp, ncs + 1, MPI_INTEGER,
!coarray     &             ipx_dest, myrank,
!coarray     &             irsbuf_p, ncsr, MPI_INTEGER,
!coarray     &             ipx_src, ipx_src,
!coarray     &             mpi_comm_world, istatus, ierr )
!      sync all ! do Not change
      call xmp_sync_all(status)
      !irsbuf_p(1:ncs+1)[ipx_dest+1] = isbufp(1:ncs+1) ! Put

      call xmp_array_section_set_triplet(irsbuf_p_sec,
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)

      call xmp_array_section_set_triplet(isbufp_local_sec, 
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put_local(img_dims,irsbuf_p_desc,irsbuf_p_sec, 
     & isbufp_local_desc,isbufp_local_sec,status)
!      sync all
      call xmp_sync_all(status)
!!

      ncar = irsbuf_p(ncsr)
!coarray      call mpi_sendrecv(buffp, 6*nca, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_dest, myrank,
!coarray     &             rbuff_p, 6*ncar, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_src, ipx_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !rbuff_p(1:6,1:nca)[ipx_dest+1] = buffp(1:6,1:nca) ! Put

      call xmp_array_section_set_triplet(rbuff_p_sec,
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuff_p_sec,
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(buffp_local_sec, 
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(buffp_local_sec, 
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put_local(img_dims,rbuff_p_desc,rbuff_p_sec, 
     & buffp_local_desc,buffp_local_sec,status)

!!

!coarray      call mpi_sendrecv(ibuffp, nca, MPI_INTEGER, 
!coarray     &             ipx_dest, myrank,
!coarray     &             irbuff_p, ncar, MPI_INTEGER, 
!coarray     &             ipx_src, ipx_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irbuff_p(1:nca)[ipx_dest+1] = ibuffp(1:nca) ! Put

      call xmp_array_section_set_triplet(irbuff_p_sec,
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(ibuffp_local_sec, 
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put_local(img_dims,irbuff_p_desc,irbuff_p_sec, 
     & ibuffp_local_desc,ibuffp_local_sec,status)

!      sync all
      call xmp_sync_all(status)
!!

!
! merge nearest neighbors(+X) on receive buffer to local bucket.
      icx = 2
      icy0 = 2
      icy1 = ncydiv + 1
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*nczdiv*ncydiv + 1
      ncc_p = 0
      ncs_p = loc_csbound
      nca_p = 0
!
      call add_bucket(abucket,iabucket,isbucket,
     &                rbuff_p,irbuff_p,irsbuf_p,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                ncc_p,ncs_p,nca_p)

!
! ( bound -X )
!
! nearest neighbor.
      icz0 = 2
      icz1 = nczdiv + 1
      icy0 = 2
      icy1 = ncydiv + 1
      icx = 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*nczdiv*ncydiv + 1
      ncc = 0
      ncs = loc_csbound
      nca = 0
!
      call add_buffer(abucket,iabucket,isbucket,buffm,ibuffm,isbufm,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
!
      ipx_dest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
      ipx_src   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)

!
      call mpi_sendrecv(ncs + 1, 1, MPI_INTEGER,
     &             ipx_dest, myrank,
     &             ncsr, 1, MPI_INTEGER,
     &             ipx_src, ipx_src,
     &             mpi_comm_world, istatus, ierr )

      isbufm(ncs+1) = nca
!coarray      call mpi_sendrecv(isbufm, ncs + 1, MPI_INTEGER,
!coarray     &             ipx_dest, myrank,
!coarray     &             irsbuf_m, ncsr, MPI_INTEGER,
!coarray     &             ipx_src, ipx_src,
!coarray     &             mpi_comm_world, istatus, ierr )
!      sync all
      call xmp_sync_all(status)
      !irsbuf_m(1:ncs+1)[ipx_dest+1] = isbufm(1:ncs+1) ! Put

      call xmp_array_section_set_triplet(irsbuf_m_sec,
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)

      call xmp_array_section_set_triplet(isbufm_local_sec, 
     & 1,int(1,kind=8),int(ncs+1,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put_local(img_dims,irsbuf_m_desc,irsbuf_m_sec, 
     & isbufm_local_desc,isbufm_local_sec,status)

!      sync all
      call xmp_sync_all(status)
!!

      ncar = irsbuf_m(ncsr)
!coarray      call mpi_sendrecv(buffm, 6*nca, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_dest, myrank,
!coarray     &             rbuff_m, 6*ncar, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_src, ipx_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !rbuff_m(1:6,1:nca)[ipx_dest+1] = buffm(1:6,1:nca) ! Put

      call xmp_array_section_set_triplet(rbuff_m_sec,
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuff_m_sec,
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(buffm_local_sec, 
     & 1,int(1,kind=8),int(6,kind=8),1,status)
      call xmp_array_section_set_triplet(buffm_local_sec, 
     & 2,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put_local(img_dims,rbuff_m_desc,rbuff_m_sec, 
     & buffm_local_desc,buffm_local_sec,status)

!!

!coarray      call mpi_sendrecv(ibuffm, nca, MPI_INTEGER, 
!coarray     &             ipx_dest, myrank,
!coarray     &             irbuff_m, ncar, MPI_INTEGER, 
!coarray     &             ipx_src, ipx_src,
!coarray     &             mpi_comm_world, istatus, ierr )
      !irbuff_m(1:nca)[ipx_dest+1] = ibuffm(1:nca) ! Put

      call xmp_array_section_set_triplet(irbuff_m_sec,
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      call xmp_array_section_set_triplet(ibuffm_local_sec, 
     & 1,int(1,kind=8),int(nca,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put_local(img_dims,irbuff_m_desc,irbuff_m_sec, 
     & ibuffm_local_desc,ibuffm_local_sec,status)

!      sync all
      call xmp_sync_all(status)
!!

!
! merge nearest neighbors(-X) on receive buffer to local bucket.
      icx = ncxdiv + 1
      icy0 = 2
      icy1 = ncydiv + 1
      icz0 = 2
      icz1 = nczdiv + 1
      ldcell = 1
      lcse(1,1,ldcell) = icz0
      lcse(2,1,ldcell) = icz1
      lcse(1,2,ldcell) = icy0
      lcse(2,2,ldcell) = icy1
      lcse(1,3,ldcell) = icx
      lcse(2,3,ldcell) = icx
!
      loc_csbound = 2*nczdiv*ncydiv + 1
      ncc_m = 0
      ncs_m = loc_csbound
      nca_m = 0
!
      call add_bucket(abucket,iabucket,isbucket,
     &                rbuff_m,irbuff_m,irsbuf_m,
     &                max_mvatom,max_mvseg,max_cellcbd,
     &                nczdiv,ncydiv,ncxdiv,
     &                lcse,ldcell,ncatom,ncseg,
     &                ncc_m,ncs_m,nca_m)

!
! ------- create new cell meta-data and its entity ---------------
!         create segment meta-data and its entity, 
!         ie. the number of atoms per segment.
!
! copy atom data from bucket to new meta-data structure.
! also, setup meta-data "tag" and "na_per_cell", 
! and also, setup segment meta-data "lsegtop" and "lseg_natoms".
!
      narea = na1cell * (nczdiv + 4) * (ncydiv + 4)
      naline = na1cell * (nczdiv + 4)
!
! atom data.
      nbase = narea
      do icx = 2, ncxdiv + 1
         nbase = nbase + narea
         nbase2 = nbase + naline
         do icy = 2, ncydiv + 1
            nbase2 = nbase2 + naline
            nca = nbase2 + 2*na1cell
            do icz = 2, nczdiv + 1
               do ica = 1, ncatom(icz,icy,icx)
                  nca = nca + 1
                  wkxyz(1,nca) = abucket(1,ica,icz,icy,icx)
                  wkxyz(2,nca) = abucket(2,ica,icz,icy,icx)
                  wkxyz(3,nca) = abucket(3,ica,icz,icy,icx)
                  wkv(1,nca) = abucket(4,ica,icz,icy,icx)
                  wkv(2,nca) = abucket(5,ica,icz,icy,icx)
                  wkv(3,nca) = abucket(6,ica,icz,icy,icx)
                  m2i(nca) = iabucket(ica,icz,icy,icx)
               end do
            end do
         end do
      end do
!
! segment meta-data.
      ncs = 0
      nbase = narea
      do icx = 2, ncxdiv + 1
         nbase = nbase + narea
         nbase2 = nbase + naline
         do icy = 2, ncydiv + 1
            nbase2 = nbase2 + naline
            nca = nbase2 + 2*na1cell
            do icz = 2, nczdiv + 1
                do ics = 1, ncseg(icz,icy,icx)
                  ncs = ncs + 1
                  lseg_natoms(ncs) = isbucket(ics,icz,icy,icx)
                  lsegtop(ncs) = nca + 1
                  nca = nca + lseg_natoms(ncs)
               end do
             end do
         end do
      end do
! number of segment per process.
      nselfseg = ncs
! cell meta-data.
      nbase = narea
      do icx = 2, ncxdiv + 1
         nbase = nbase + narea
         nbase2 = nbase + naline
         do icy = 2, ncydiv + 1
            nbase2 = nbase2 + naline
            nca = nbase2 + 2*na1cell
            do icz = 2, nczdiv + 1
               tag(icz + 1, icy + 1, icx + 1) = nca + 1
               nca = nca + ncatom(icz,icy,icx)
               na_per_cell(icz + 1, icy + 1, icx + 1)
     &         = ncatom(icz,icy,icx)
            end do
         end do
      end do
!
!update of i2m
!
!$omp parallel default(shared)
!$omp& private(i,i0,k0)
!$omp do
      do i=1,n
        i2m(i)=-1  !! not necesary for performance measurement
      enddo
!$omp end do
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        i=m2i(i0)
        i2m(i)=i0
      enddo
      enddo
!$omp end do
!$omp end parallel
!
      return
      end
c----------------------------------------------------------------------
      subroutine add_buffer(abucket,iabucket,isbucket,buff,ibuff,isbuf,
     &                      max_mvatom,max_mvseg,max_cell,
     &                      nczdiv,ncydiv,ncxdiv,
     &                      lcse,ldcell,ncatom,ncseg,ncc,ncs,nca)
c----------------------------------------------------------------------
!
      implicit none
      integer max_mvatom
      integer max_mvseg
      integer max_cell        ! max number of cells on communication buffer.
      integer nczdiv
      integer ncydiv
      integer ncxdiv
      real(8) abucket(6,max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
      integer iabucket(max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
      integer isbucket(max_mvseg,nczdiv+2,ncydiv+2,ncxdiv+2)
      real(8) buff(6,max_cell*max_mvatom)
      integer ibuff(max_cell*max_mvatom)
      integer isbuf(2*max_cell + 1 + max_cell*max_mvseg)
!
      integer ldcell, ldc
      integer lcse(2,3,ldcell)
      integer ncatom(nczdiv+2,ncydiv+2,ncxdiv+2)
      integer ncseg(nczdiv+2,ncydiv+2,ncxdiv+2)
      integer ncc
      integer ncs
      integer nca
      integer ics, ica
      integer icz, icy, icx
      integer icz0, icy0, icx0
      integer icz1, icy1, icx1
!
      do ldc = 1, ldcell
         icz0 = lcse(1,1,ldc)
         icz1 = lcse(2,1,ldc)
         icy0 = lcse(1,2,ldc)
         icy1 = lcse(2,2,ldc)
         icx0 = lcse(1,3,ldc)
         icx1 = lcse(2,3,ldc)
         do icx = icx0, icx1
            do icy = icy0, icy1
               do icz = icz0, icz1
                  ncc = ncc + 1
                  isbuf(ncc) = ncatom(icz,icy,icx)
                  ncc = ncc + 1
                  isbuf(ncc) = ncseg(icz,icy,icx)
                  do ics = 1, ncseg(icz,icy,icx)
                     ncs = ncs + 1
                     isbuf(ncs) = isbucket(ics,icz,icy,icx)
                  end do
                  do ica = 1, ncatom(icz,icy,icx)
                     nca = nca + 1
                     buff(1,nca) = abucket(1,ica,icz,icy,icx)
                     buff(2,nca) = abucket(2,ica,icz,icy,icx)
                     buff(3,nca) = abucket(3,ica,icz,icy,icx)
                     buff(4,nca) = abucket(4,ica,icz,icy,icx)
                     buff(5,nca) = abucket(5,ica,icz,icy,icx)
                     buff(6,nca) = abucket(6,ica,icz,icy,icx)
                     ibuff(nca)  = iabucket(ica,icz,icy,icx)
                  end do
               end do
            end do
         end do
      end do
!
      return
      end
!
c----------------------------------------------------------------------
      subroutine add_buffb(rbuff,irbuff,irsbuf,buff,ibuff,isbuf,
     &                    max_mvatom,max_mvseg,max_cell,
     &                    nccr,ncsr,ncar,ncc,ncs,nca,ldcell)
c----------------------------------------------------------------------
!
      implicit none
      integer max_mvatom
      integer max_mvseg
      integer max_cell        ! max number of cells on communication buffer.

      real(8) rbuff (6,max_cell*max_mvatom)
      integer irbuff(max_cell*max_mvatom)
      integer irsbuf(2*max_cell + 1 + max_cell*max_mvseg)
      real(8) buff  (6,max_cell*max_mvatom)
      integer ibuff (max_cell*max_mvatom)
      integer isbuf (2*max_cell + 1 + max_cell*max_mvseg)

      integer ncc, nccr
      integer ncs, ncsr
      integer nca, ncar
      integer ldcell, ldc
      integer ics, ica
!
      do ldc = 1, ldcell
         ncc = ncc + 2
         nccr = nccr + 2
         isbuf(ncc - 1) = irsbuf(nccr - 1)
         isbuf(ncc)     = irsbuf(nccr)
! segment data.
         do ics = 1, irsbuf(nccr)
            ncs = ncs + 1
            ncsr = ncsr + 1
            isbuf(ncs) = irsbuf(ncsr)
         end do
! atom data.
         do ica = 1, irsbuf(nccr - 1)
            nca = nca + 1
            ncar = ncar + 1
            buff(1,nca) = rbuff(1,ncar)
            buff(2,nca) = rbuff(2,ncar)
            buff(3,nca) = rbuff(3,ncar)
            buff(4,nca) = rbuff(4,ncar)
            buff(5,nca) = rbuff(5,ncar)
            buff(6,nca) = rbuff(6,ncar)
            ibuff(nca)  = irbuff(ncar)
         end do
      end do
!
      return
      end
!
c----------------------------------------------------------------------
      subroutine add_bucket(abucket,iabucket,isbucket,
     &                      rbuff,irbuff,irsbuf,
     &                      max_mvatom,max_mvseg,max_cell,
     &                      nczdiv,ncydiv,ncxdiv,
     &                      lcse,ldcell,ncatom,ncseg,
     &                      ncc_b,ncs_b,nca_b)
c----------------------------------------------------------------------
!
      implicit none
      integer max_mvatom
      integer max_mvseg
      integer max_cell        ! max number of cells on communication buffer.
      integer nczdiv
      integer ncydiv
      integer ncxdiv
      real(8) abucket(6,max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
      integer iabucket(max_mvatom,nczdiv+2,ncydiv+2,ncxdiv+2)
      integer isbucket(max_mvseg,nczdiv+2,ncydiv+2,ncxdiv+2)
      real(8) rbuff(6,max_cell*max_mvatom)
      integer irbuff(max_cell*max_mvatom)
      integer irsbuf(2*max_cell + 1 + max_cell*max_mvseg)

      integer ldcell, ldc
      integer lcse(2,3,ldcell)
      integer ncatom(nczdiv+2,ncydiv+2,ncxdiv+2)
      integer ncseg(nczdiv+2,ncydiv+2,ncxdiv+2)
      integer ncc_b
      integer ncs_b, ncs, ics
      integer nca_b, nca, ica
      integer icz, icy, icx
      integer icz0, icy0, icx0
      integer icz1, icy1, icx1
!
      do ldc = 1, ldcell
         icz0 = lcse(1,1,ldc)
         icz1 = lcse(2,1,ldc)
         icy0 = lcse(1,2,ldc)
         icy1 = lcse(2,2,ldc)
         icx0 = lcse(1,3,ldc)
         icx1 = lcse(2,3,ldc)
         do icx = icx0, icx1
            do icy = icy0, icy1
               do icz = icz0, icz1
! segment data.
                  ncc_b = ncc_b + 2                     ! ncc_b
                  ncs = ncseg(icz,icy,icx)
                  do ics = 1, irsbuf(ncc_b)
                     ncs_b = ncs_b + 1                  ! ncs_b
                     ncs = ncs + 1
                     isbucket(ncs,icz,icy,icx) = irsbuf(ncs_b)
                  end do
                  ncseg(icz,icy,icx) = ncs
! atom data.
                  nca = ncatom(icz,icy,icx)
                  do ica = 1, irsbuf(ncc_b - 1)
                     nca_b = nca_b + 1                  ! nca_b
                     nca = nca + 1
                     abucket(1,nca,icz,icy,icx) = rbuff(1,nca_b)
                     abucket(2,nca,icz,icy,icx) = rbuff(2,nca_b)
                     abucket(3,nca,icz,icy,icx) = rbuff(3,nca_b)
                     abucket(4,nca,icz,icy,icx) = rbuff(4,nca_b)
                     abucket(5,nca,icz,icy,icx) = rbuff(5,nca_b)
                     abucket(6,nca,icz,icy,icx) = rbuff(6,nca_b)
                     iabucket(nca,icz,icy,icx) = irbuff(nca_b)
                  end do
                  ncatom(icz,icy,icx) = nca
               end do
            end do
         end do
      end do
!
      return
      end
c----------------------------------------------------------------------
      subroutine pre_record_data
c----------------------------------------------------------------------
      use trj_mpi
      use trj_org
      use md_fmm 
      use md_fmm_domdiv_flg
      use md_segment
      use mpivar
      use xmp_api
      implicit none
      integer(4) :: i,nsum
      integer(4) :: i0,k0,i00
      include 'mpif.h'
      integer(4) :: ierr
      real(8),allocatable :: snd(:,:),rcv(:,:)
      integer(8) :: snd_local_sec
      integer(8) :: snd_local_desc
      integer(8), dimension(2) :: snd_lb, snd_ub

!coarray      integer(4),allocatable :: natmlist(:),natmdisp(:)
      integer(4),allocatable :: natmdisp(:)
!      integer(4),allocatable :: natmlist(:)[:]
      integer(4), POINTER :: natmlist(:) => null ()
      integer(8) :: natmlist_desc
      integer(8) :: natmlist_sec
      integer(8), dimension(1) :: natmlist_lb, natmlist_ub

      integer(4),allocatable :: natmlist_tmp(:)
!      integer,allocatable :: ndis(:)[:], mdis(:)[:]
      integer, POINTER :: ndis(:) => null ()
      integer, POINTER :: mdis(:) => null ()
      integer(8) :: mdis_sec
      integer(8) :: mdis_desc
      integer(8), dimension(1) :: mdis_lb, mdis_ub

      integer(8) :: ndis_desc
      integer(8), dimension(1) :: ndis_lb, ndis_ub


!      real(8),allocatable :: rcvx(:,:)[:]
      real(8), POINTER :: rcvx(:,:) => null ()
      integer(8) :: rcvx_desc
      integer(8) :: rcvx_sec
      integer(8), dimension(2) :: rcvx_lb, rcvx_ub

      integer :: me, np, ms, mm
!!
      integer(4),allocatable :: nrearrange(:)
      integer(4) :: m2i_tmp(na1cell*lxdiv*lydiv*lzdiv)
      integer(8) :: m2i_tmp_local_sec
      integer(8) :: m2i_tmp_local_desc
      integer(8),dimension(1) :: m2i_tmp_lb,m2i_tmp_ub

      integer(4) :: status


!coarray
      !me = this_image()
      me = xmp_this_image()
      !np = num_images()
      np = xmp_num_images()
!      allocate(ndis(np)[*])
!      allocate(mdis(n)[*])
!      allocate(rcvx(6,n)[*])
      ndis_lb(1) = 1
      ndis_ub(1) = np
      call xmp_new_coarray(ndis_desc,4,1,ndis_lb,ndis_ub,1, img_dims)
      call xmp_coarray_bind(ndis_desc,ndis)

      mdis_lb(1) = 1
      mdis_ub(1) = n
      call xmp_new_coarray(mdis_desc,4,1,mdis_lb,mdis_ub,1, img_dims)
      call xmp_coarray_bind(mdis_desc,mdis)

      rcvx_lb(1) = 1
      rcvx_lb(2) = 1
      rcvx_ub(1) = 6
      rcvx_ub(2) = n
      call xmp_new_coarray(rcvx_desc,8,2,rcvx_lb,rcvx_ub,1, img_dims)
      call xmp_coarray_bind(rcvx_desc,rcvx)
      call xmp_new_array_section(rcvx_sec,2)

      m2i_tmp_lb(1) = 1
      m2i_tmp_ub(1) = na1cell*lxdiv*lydiv*lzdiv
      call xmp_new_local_array(m2i_tmp_local_desc,4,1,
     & m2i_tmp_lb,m2i_tmp_ub,loc(m2i_tmp))
      call xmp_new_array_section(m2i_tmp_local_sec,1)


!!

      if(nprocs.eq.1) then
!$omp parallel do default(shared)
!$omp& private(i,i0)
        do i = 1,n
          i0=i2m(i)
          xyz(1,i) = wkxyz(1,i0)
          xyz(2,i) = wkxyz(2,i0)
          xyz(3,i) = wkxyz(3,i0)
          v(1:3,i) = wkv(1:3,i0)
        end do
      else
        allocate(nrearrange(n))
        allocate(snd(6,n))
        snd_lb(1) = 1
        snd_lb(2) = 1
        snd_ub(1) = 6
        snd_ub(2) = n
        call xmp_new_local_array(snd_local_desc,8,2,
     &   snd_lb,snd_ub,loc(snd))
        call xmp_new_array_section(snd_local_sec,2)


        allocate(rcv(6,n))
!coarray        allocate(natmlist(nprocs),natmdisp(nprocs))
        !allocate(natmlist(nprocs)[*])
        natmlist_lb(1) = 1
        natmlist_ub(1) = nprocs
        call xmp_new_coarray(natmlist_desc,4,1,
     &   natmlist_lb,natmlist_ub,1,img_dims)
        call xmp_coarray_bind(natmlist_desc,natmlist)
        call xmp_new_array_section(natmlist_sec,1)

        allocate(natmlist_tmp(nprocs))
        allocate(natmdisp(nprocs))
!!

        nselfatm=0
        do k0=1,nselfseg
        do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          nselfatm=nselfatm+1
          m2i_tmp(nselfatm)=m2i(i0)
        enddo ! i0
        enddo ! k0

        if(nselfatm.gt.na1cell*lxdiv*lydiv*lzdiv)then
          write(*,*) 'ERROR: nselfatm overflowed!', myrank,nselfatm
          call mpistop()
        endif

!coarray        call mpi_allgather(nselfatm,1,mpi_integer,
!coarray     &                     natmlist,1,mpi_integer,
!coarray     &                     mpi_comm_world,ierr)
!coarray!
!coarray        call mpi_barrier(mpi_comm_world,ierr)
!coarray!
!coarray        natmdisp(1) = 0
        do mm = 1,np
          !natmlist(me)[mm] = nselfatm ! Put
          call xmp_array_section_set_triplet(natmlist_sec,
     &     1,int(me,kind=8),int(me,kind=8),1,status)
          img_dims(1) = mm
          call xmp_coarray_put_scalar(img_dims,natmlist_desc,
     &     natmlist_sec,nselfatm,status)

!          sync all
          call xmp_sync_all(status)
        enddo
        natmdisp(1) = 1
!!
        nsum = natmlist(1)
        do i = 2,nprocs
          natmdisp(i) = natmdisp(i-1)+natmlist(i-1)
          nsum        = nsum + natmlist(i)
        end do
!
!coarray        call mpi_gatherv(m2i_tmp,nselfatm,mpi_integer,
!coarray     &       nrearrange,natmlist,natmdisp,mpi_integer,
!coarray     &                 mpiout,mpi_comm_world,ierr)
        ms = natmdisp(me)
        !mdis(ms:ms+nselfatm-1)[mpiout+1] = m2i_tmp(1:nselfatm)

        call xmp_array_section_set_triplet(mdis_sec,
     &   1,int(ms,kind=8),int(ms+nselfatm-1,kind=8),1,status)

        call xmp_array_section_set_triplet(m2i_tmp_local_sec, 
     &   1,int(1,kind=8),int(nselfatm,kind=8),1,status)

        img_dims(1) = mpiout+1
        call xmp_coarray_put_local(img_dims,mdis_desc,mdis_sec, 
     &   m2i_tmp_local_desc,m2i_tmp_local_sec,status)


!        sync all
        call xmp_sync_all(status)
        nrearrange = mdis
!!
!
!$omp parallel do default(shared)
!$omp& private(i)
        do i = 1,nprocs
          natmlist(i) = natmlist(i)*6
          natmdisp(i) = natmdisp(i)*6
        end do
!
        i00=0
        do k0=1,nselfseg
        do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          i00=i00+1
          snd(1,i00) = wkxyz(1,i0)
          snd(2,i00) = wkxyz(2,i0)
          snd(3,i00) = wkxyz(3,i0)
          snd(4,i00) = wkv(1,i0)
          snd(5,i00) = wkv(2,i0)
          snd(6,i00) = wkv(3,i0)
        end do ! i0
        end do ! k0

!coarray        call mpi_gatherv(snd,nselfatm*6,mpi_double_precision,
!coarray     &            rcv,natmlist,natmdisp,mpi_double_precision,
!coarray     &            mpiout,mpi_comm_world,ierr)
        ms = natmdisp(me)/6
        !rcvx(1:6,ms:ms+nselfatm-1)[mpiout+1] = snd(1:6,1:nselfatm)

        call xmp_array_section_set_triplet(rcvx_sec,
     &   1,int(1,kind=8),int(6,kind=8),1,status)
        call xmp_array_section_set_triplet(rcvx_sec,
     &   2,int(ms,kind=8),int(ms+nselfatm-1,kind=8),1,status)


        call xmp_array_section_set_triplet(snd_local_sec, 
     &   1,int(1,kind=8),int(6,kind=8),1,status)
        call xmp_array_section_set_triplet(snd_local_sec, 
     &   2,int(1,kind=8),int(nselfatm,kind=8),1,status)

        img_dims(1) = mpiout+1
        call xmp_coarray_put_local(img_dims,rcvx_desc,rcvx_sec, 
     &   snd_local_desc,snd_local_sec,status)

!        sync all
        call xmp_sync_all(status)
        rcv = rcvx
!!
!
        if(myrank.eq.mpiout) then
!$omp parallel do default(none)
!$omp& private(i,i0)
!$omp& shared(xyz,v,rcv,n,nrearrange)
          do i = 1,n
            i0=nrearrange(i)
            xyz(1,i0) = rcv(1,i)
            xyz(2,i0) = rcv(2,i)
            xyz(3,i0) = rcv(3,i)
              v(1,i0) = rcv(4,i)
              v(2,i0) = rcv(5,i)
              v(3,i0) = rcv(6,i)
          end do
        end if
        deallocate(snd,rcv,nrearrange)
!coarray        deallocate(natmlist,natmdisp)
        deallocate(natmdisp)
!!
      end if

      call cell_edge()

      return
      end
