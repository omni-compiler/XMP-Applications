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
      subroutine init_comm_direct_3() 
c----------------------------------------------------------------------
      use comm_base
      use comm_d3
      use trj_mpi
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use xmp_api
      implicit none
      INCLUDE 'mpif.h'

      npz = nzdiv
      npy = nydiv
      npx = nxdiv
      ipz = izflg(myrank)-1
      ipy = iyflg(myrank)-1
      ipx = ixflg(myrank)-1
      ncxdiv = lxdiv
      ncydiv = lydiv
      nczdiv = lzdiv

!     allocate(icbufp ((ncell/npy)*(ncell/npx)*2)[*])
      icbufp_lb(1) = 1
      icbufp_ub(1) = (ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(icbufp_desc,4,1,
     & icbufp_lb,icbufp_ub,1,img_dims)
      call xmp_coarray_bind(icbufp_desc,icbufp)

!     allocate(ircbufp((ncell/npy)*(ncell/npx)*2)[*])
      ircbufp_lb(1) = 1
      ircbufp_ub(1) = (ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(ircbufp_desc,4,1,
     & ircbufp_lb,ircbufp_ub,1,img_dims)
      call xmp_coarray_bind(ircbufp_desc,ircbufp)

!      allocate(icbufm ((ncell/npy)*(ncell/npx)*2)[*])
      icbufm_lb(1) = 1
      icbufm_ub(1) = (ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(icbufm_desc,4,1,
     & icbufm_lb,icbufm_ub,1,img_dims)
      call xmp_coarray_bind(icbufm_desc,icbufm)

!      allocate(ircbufm((ncell/npy)*(ncell/npx)*2)[*])
      ircbufm_lb(1) = 1
      ircbufm_ub(1) = (ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(ircbufm_desc,4,1,
     & ircbufm_lb,ircbufm_ub,1,img_dims)
      call xmp_coarray_bind(ircbufm_desc,ircbufm)

!      allocate(ibuffp (na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      ibuffp_lb(1) = 1
      ibuffp_ub(1) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(ibuffp_desc,4,1,
     & ibuffp_lb,ibuffp_ub,1,img_dims)
      call xmp_coarray_bind(ibuffp_desc,ibuffp)

!      allocate(irbuffp(na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      irbuffp_lb(1) = 1
      irbuffp_ub(1) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(irbuffp_desc,4,1,
     & irbuffp_lb,irbuffp_ub,1,img_dims)
      call xmp_coarray_bind(irbuffp_desc,irbuffp)

!      allocate(ibuffm (na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      ibuffm_lb(1) = 1
      ibuffm_ub(1) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(ibuffm_desc,4,1,
     & ibuffm_lb,ibuffm_ub,1,img_dims)
      call xmp_coarray_bind(ibuffm_desc,ibuffm)

!      allocate(irbuffm(na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      irbuffm_lb(1) = 1
      irbuffm_ub(1) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(irbuffm_desc,4,1,
     & irbuffm_lb,irbuffm_ub,1,img_dims)
      call xmp_coarray_bind(irbuffm_desc,irbuffm)

!      allocate(buffp (3,na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      buffp_lb(1) = 1
      buffp_ub(1) = 3
      buffp_lb(2) = 1
      buffp_ub(2) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(buffp_desc,8,2,
     & buffp_lb,buffp_ub,1,img_dims)
      call xmp_coarray_bind(buffp_desc,buffp)

!      allocate(rbuffp(3,na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      rbuffp_lb(1) = 1
      rbuffp_ub(1) = 3
      rbuffp_lb(2) = 1
      rbuffp_ub(2) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(rbuffp_desc,8,2,
     & rbuffp_lb,rbuffp_ub,1,img_dims)
      call xmp_coarray_bind(rbuffp_desc,rbuffp)

!      allocate(buffm (3,na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      buffm_lb(1) = 1
      buffm_ub(1) = 3
      buffm_lb(2) = 1
      buffm_ub(2) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(buffm_desc,8,2,
     & buffm_lb,buffm_ub,1,img_dims)
      call xmp_coarray_bind(buffm_desc,buffm)

!      allocate(rbuffm(3,na1cell*(ncell/npy)*(ncell/npx)*2)[*])
      rbuffm_lb(1) = 1
      rbuffm_ub(1) = 3
      rbuffm_lb(2) = 1
      rbuffm_ub(2) = na1cell*(ncell/npy)*(ncell/npx)*2
      call xmp_new_coarray(rbuffm_desc,8,2,
     & rbuffm_lb,rbuffm_ub,1,img_dims)
      call xmp_coarray_bind(rbuffm_desc,rbuffm)


      return
      end
c----------------------------------------------------------------------
       subroutine comm_direct_3()   ! ver.20120314
c----------------------------------------------------------------------
       use comm_base
       use comm_d3
       use trj_org
       use trj_mpi
       use md_forces
       use md_monitors
       use md_fmm
       use md_fmm_domdiv_flg
       use md_segment
       use md_periodic
       use unitcell
      use mpivar
      use xmp_api
      implicit none
      INCLUDE 'mpif.h'
      integer ipz_pdest, ipy_pdest, ipx_pdest
      integer ipz_psrc, ipy_psrc, ipx_psrc
      integer ipz_mdest, ipy_mdest, ipx_mdest
      integer ipz_msrc, ipy_msrc, ipx_msrc
      integer itr, nitr
      integer icz, icy, icx
      integer icz0, icz1
      integer icy0, icy1
      integer iczp0, iczp1
      integer icyp0, icyp1
      integer icxp0, icxp1
      integer iczm0, iczm1
      integer icym0, icym1
      integer icxm0, icxm1
      integer iczb, icyb, icxb
      integer iczbp0, iczbp1
      integer iczbm0, iczbm1
      integer icybp0, icybp1, icybp1st
      integer icybm0, icybm1, icybm1st
      integer icxbp0, icxbp1
#ifndef HALFDIREE
      integer icxbp1st
#endif
      integer icxbm0, icxbm1, icxbm1st
      integer ncc,ncc2
      integer nccp
      integer nccm
      integer ica, icag
      integer icasp, icarp
      integer icasm, icarm
      integer nca
      integer ncap, ncarp, ncar2p
      integer ncam, ncarm, ncar2m
      integer nbase, nbase2, nbase3
      integer ntmp
      integer istatus(mpi_status_size, 8), ierr
#ifndef SYNC_COM
      integer,dimension(8) :: irq
      integer nrq
#endif
!coarray
      integer nd
      integer(4) status
!!
      call xmp_new_array_section(rbuffm_sec,2)
      call xmp_new_array_section(buffm_sec,2)
      call xmp_new_array_section(buffp_sec,2)
      call xmp_new_array_section(rbuffp_sec,2)
      call xmp_new_array_section(irbuffm_sec,1)
      call xmp_new_array_section(ibuffm_sec,1)
      call xmp_new_array_section(irbuffp_sec,1)
      call xmp_new_array_section(ibuffp_sec,1)
      call xmp_new_array_section(ircbufm_sec,1)
      call xmp_new_array_section(icbufm_sec,1)
      call xmp_new_array_section(icbufp_sec,1)
      call xmp_new_array_section(ircbufp_sec,1)
      call xmp_new_array_section(na_per_cell_l_sec,3)
      call xmp_new_array_section(na_per_cell_r_sec,3)
      call xmp_new_array_section(wkxyz_l_sec,2)
      call xmp_new_array_section(wkxyz_r_sec,2)
      call xmp_new_array_section(m2i_l_sec,1)
      call xmp_new_array_section(m2i_r_sec,1)

c----- common parameters for coordinate communication. -----
      ipx=mod(myrank,npx)
      ipy=mod((myrank-ipx)/npx,npy)
      ipz=mod((myrank-ipx-ipy*npx)/(npx*npy),npz)

      nczdiv = (ncell - 1)/npz + 1
      ncydiv = (ncell - 1)/npy + 1
      ncxdiv = (ncell - 1)/npx + 1

      narea = na1cell * (ncell/npz + 4) * (ncell/npy + 4)
      naline = na1cell * (ncell/npz + 4)

!
!-----  coordinate communication code starts here. ------
!
!     coordinate +Z
      ipz_pdest = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx
      ipz_psrc  = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
!     coordinate -Z
      ipz_mdest = mod(ipz-1+1/npz+npz,npz)*npx*npy + ipy*npx + ipx
      ipz_msrc  = mod(ipz+1-1/npz+npz,npz)*npx*npy + ipy*npx + ipx


      nitr = (2 - 1)/nczdiv + 1

      DO itr = 1, nitr
         if (itr == 1) then
            iczp0 = 2+ nczdiv - 1
            if (nczdiv == 1) iczp0 = iczp0+1
            iczp1 = 2+ nczdiv

            nccp = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO icz = iczp0, iczp1
                     nccp = nccp + 1
                     icbufp(nccp) = na_per_cell( icz, icy, icx )
                  END DO
               END DO
            END DO

            iczm0 = 3
            iczm1 = iczm0 + 1
            if (nczdiv == 1) iczm1 = iczm0
            nccm = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO icz = iczm0, iczm1
                     nccm = nccm + 1
                     icbufm(nccm) = na_per_cell( icz, icy, icx )
                  END DO
               END DO
            END DO

#ifdef SYNC_COM
!coarray            call mpi_sendrecv(icbufp, nccp, MPI_INTEGER,
!coarray     &             ipz_pdest, myrank,
!coarray     &             ircbufp, nccp, MPI_INTEGER, ipz_psrc, ipz_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(icbufm, nccm, MPI_INTEGER,
!coarray     &             ipz_mdest, myrank,
!coarray     &             ircbufm, nccm, MPI_INTEGER, ipz_msrc, ipz_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )

      !ircbufp(1:nccp)[ipz_pdest+1] = icbufp(1:nccp) ! Put
      call xmp_array_section_set_triplet(ircbufp_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(icbufp_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,ircbufp_desc,ircbufp_sec,
     & icbufp_desc,icbufp_sec,status)

      !ircbufm(1:nccm)[ipz_mdest+1] = icbufm(1:nccm) ! Put
      call xmp_array_section_set_triplet(ircbufm_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(icbufm_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,ircbufm_desc,ircbufm_sec,
     & icbufm_desc,icbufm_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
            call mpi_irecv(ircbufp, nccp,
     &              MPI_INTEGER, ipz_psrc, ipz_psrc, mpi_comm_world,
     &              irq(1), ierr)
            call mpi_isend(icbufp, nccp,
     &              MPI_INTEGER, ipz_pdest, myrank, mpi_comm_world, 
     &              irq(2), ierr)
            call mpi_irecv(ircbufm, nccm,
     &              MPI_INTEGER, ipz_msrc, ipz_msrc, mpi_comm_world, 
     &              irq(3), ierr)
            call mpi_isend(icbufm, nccm,
     &              MPI_INTEGER, ipz_mdest, myrank, mpi_comm_world, 
     &              irq(4), ierr)
            nrq = 4
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            iczbp0 = iczp0 - nczdiv
            iczbp1 = iczp1 - nczdiv
            iczbm0 = iczm0 + nczdiv
            iczbm1 = iczm1 + nczdiv

            ncc2 = 0
            ncarp = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  nca = tag(iczbp1+1,icy,icx) - ircbufp(ncc2+1)
                  if(nczdiv > 1) then
                     nca = tag(iczbp1+1,icy,icx)
     &                    - ircbufp(ncc2+2) - ircbufp(ncc2+1)
                  END IF
                  DO iczb = iczbp0, iczbp1
                     ncc2 = ncc2 + 1
                     na_per_cell(iczb,icy,icx) = ircbufp(ncc2)
                     tag(iczb,icy,icx) = nca
                     nca = nca + na_per_cell(iczb,icy,icx)
                     ncarp = ncarp + na_per_cell(iczb,icy,icx)
                  END DO
               END DO
            END DO

            ncap = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO ica = tag(iczp0, icy, icx), tag(iczp1, icy, icx)
     &                 + na_per_cell(iczp1, icy, icx)-1
                     ncap = ncap + 1
                     buffp(1,ncap) = wkxyz(1,ica)
                     buffp(2,ncap) = wkxyz(2,ica)
                     buffp(3,ncap) = wkxyz(3,ica)
                     ibuffp(ncap) = m2i(ica)
                  END DO
               END DO
            END DO

            ncc2 = 0
            ncarm = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  nca = tag(iczbm0-1,icy,icx)
     &                 + na_per_cell(iczbm0-1,icy,icx)
                  DO iczb = iczbm0, iczbm1
                     ncc2 = ncc2 + 1
                     na_per_cell(iczb,icy,icx) = ircbufm(ncc2)
                     tag(iczb,icy,icx) = nca
                     nca = nca + na_per_cell(iczb,icy,icx)
                     ncarm = ncarm + na_per_cell(iczb,icy,icx)
                  END DO
               END DO
            END DO
            
            ncam = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO ica = tag(iczm0, icy, icx), tag(iczm1, icy, icx)
     &                     + na_per_cell(iczm1, icy, icx)-1
                     ncam = ncam + 1
                     buffm(1,ncam) = wkxyz(1,ica)
                     buffm(2,ncam) = wkxyz(2,ica)
                     buffm(3,ncam) = wkxyz(3,ica)
                     ibuffm(ncam) = m2i(ica)
                  END DO
               END DO
            END DO

#ifdef SYNC_COM
!coarray            call mpi_sendrecv(buffp, 3*ncap, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_pdest, myrank,
!coarray     &             rbuffp, 3*ncarp, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_psrc, ipz_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(ibuffp, ncap, MPI_INTEGER, 
!coarray     &             ipz_pdest, myrank,
!coarray     &             irbuffp, ncarp, MPI_INTEGER, ipz_psrc, ipz_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray
!coarray            call mpi_sendrecv(buffm, 3*ncam, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_mdest, myrank,
!coarray     &             rbuffm, 3*ncarm, MPI_DOUBLE_PRECISION,
!coarray     &             ipz_msrc, ipz_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(ibuffm, ncam, MPI_INTEGER, ipz_mdest,
!coarray     &             myrank, irbuffm, ncarm, MPI_INTEGER,
!coarray     &             ipz_msrc, ipz_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )

      !rbuffp(1:3,1:ncap)[ipz_pdest+1] = buffp(1:3,1:ncap) ! Put
      call xmp_array_section_set_triplet(rbuffp_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuffp_sec,
     & 2,int(1,kind=8),int(ncap,kind=8),1,status)
      call xmp_array_section_set_triplet(buffp_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(buffp_sec,
     & 2,int(1,kind=8),int(ncap,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,rbuffp_desc,rbuffp_sec,
     & buffp_desc,buffp_sec,status)

      !irbuffp(1:ncap)[ipz_pdest+1]    = ibuffp(1:ncap)    ! Put
      call xmp_array_section_set_triplet(irbuffp_sec,
     & 1,int(1,kind=8),int(ncap,kind=8),1,status)
      call xmp_array_section_set_triplet(ibuffp_sec,
     & 1,int(1,kind=8),int(ncap,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,irbuffp_desc,irbuffp_sec,
     & ibuffp_desc,ibuffp_sec,status)


      !rbuffm(1:3,1:ncam)[ipz_mdest+1] = buffm(1:3,1:ncam) ! Put
      call xmp_array_section_set_triplet(rbuffm_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuffm_sec,
     & 2,int(1,kind=8),int(ncam,kind=8),1,status)
      call xmp_array_section_set_triplet(buffm_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(buffm_sec,
     & 2,int(1,kind=8),int(ncam,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,rbuffm_desc,rbuffm_sec,
     & buffm_desc,buffm_sec,status)

      !irbuffm(1:ncam)[ipz_mdest+1]    = ibuffm(1:ncam)    ! Put 
      call xmp_array_section_set_triplet(irbuffm_sec,
     & 1,int(1,kind=8),int(ncam,kind=8),1,status)
      call xmp_array_section_set_triplet(ibuffm_sec,
     & 1,int(1,kind=8),int(ncam,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,irbuffm_desc,irbuffm_sec,
     & ibuffm_desc,ibuffm_sec,status)

!      sync all
      call xmp_sync_all(status)
!!

#else
            call mpi_irecv(rbuffp, 3*ncarp, 
     &             MPI_DOUBLE_PRECISION, ipz_psrc, ipz_psrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(buffp, 3*ncap, 
     &             MPI_DOUBLE_PRECISION, ipz_pdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(irbuffp, ncarp, 
     &             MPI_INTEGER, ipz_psrc, ipz_psrc,
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(ibuffp, ncap, 
     &             MPI_INTEGER, ipz_pdest, myrank,
     &             mpi_comm_world, irq(4), ierr)

            call mpi_irecv(rbuffm, 3*ncarm, 
     &             MPI_DOUBLE_PRECISION, ipz_msrc, ipz_msrc, 
     &             mpi_comm_world, irq(5), ierr)
            call mpi_isend(buffm, 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipz_mdest, myrank, 
     &             mpi_comm_world, irq(6), ierr)
            call mpi_irecv(irbuffm, ncarm, 
     &             MPI_INTEGER, ipz_msrc, ipz_msrc,
     &             mpi_comm_world, irq(7), ierr)
            call mpi_isend(ibuffm, ncam, 
     &             MPI_INTEGER, ipz_mdest, myrank,
     &             mpi_comm_world, irq(8), ierr)

            nrq = 8
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            nca = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO ica = tag(iczbp0, icy, icx), tag(iczbp1, icy, icx)
     &                      + na_per_cell(iczbp1, icy, icx)-1
                     nca = nca + 1
                     wkxyz(1,ica) = rbuffp(1,nca)
                     wkxyz(2,ica) = rbuffp(2,nca)
                     wkxyz(3,ica) = rbuffp(3,nca)
                     m2i(ica) = irbuffp(nca)
                  END DO
               END DO
            END DO

            nca = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO ica = tag(iczbm0, icy, icx), tag(iczbm1, icy, icx)
     &                 + na_per_cell(iczbm1, icy, icx)-1
                     nca = nca + 1
                     wkxyz(1,ica) = rbuffm(1,nca)
                     wkxyz(2,ica) = rbuffm(2,nca)
                     wkxyz(3,ica) = rbuffm(3,nca)
                     m2i(ica) = irbuffm(nca)
                  END DO
               END DO
            END DO

         ELSE

#ifdef SYNC_COM
!coarray            call mpi_sendrecv(ircbufp, nccp, MPI_INTEGER, 
!coarray     &              ipz_pdest, myrank,
!coarray     &              icbufp, nccp, MPI_INTEGER, ipz_psrc, ipz_psrc,
!coarray     &              mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(ircbufm, nccm, MPI_INTEGER, 
!coarray     &             ipz_mdest, myrank,
!coarray     &             icbufm, nccm, MPI_INTEGER, ipz_msrc, ipz_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )

      !icbufp(1:nccp)[ipz_pdest+1] = ircbufp(1:nccp) ! Put
      call xmp_array_section_set_triplet(icbufp_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(ircbufp_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,icbufp_desc,icbufp_sec,
     & ircbufp_desc,ircbufp_sec,status)

      !icbufm(1:nccm)[ipz_mdest+1] = ircbufm(1:nccm) ! Put
      call xmp_array_section_set_triplet(icbufm_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(ircbufm_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,icbufm_desc,icbufm_sec,
     & ircbufm_desc,ircbufm_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
            call mpi_irecv(icbufp, nccp,
     &              MPI_INTEGER, ipz_psrc, ipz_psrc, mpi_comm_world, 
     &              irq(1), ierr)
            call mpi_isend(ircbufp, nccp,
     &              MPI_INTEGER, ipz_pdest, myrank, mpi_comm_world, 
     &              irq(2), ierr)
            call mpi_irecv(icbufm, nccm,
     &              MPI_INTEGER, ipz_msrc, ipz_msrc, mpi_comm_world, 
     &              irq(3), ierr)
            call mpi_isend(ircbufm, nccm,
     &              MPI_INTEGER, ipz_mdest, myrank, mpi_comm_world, 
     &              irq(4), ierr)
            nrq = 4
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            iczbp0 = 1
            iczbp1 = 1
            ncc2 = 0
            ncar2p = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  nca = tag(iczbp1+1,icy,icx) - icbufp(ncc2+1)
                  DO iczb = iczbp0, iczbp1
                     ncc2 = ncc2 + 1
                     na_per_cell(iczbp0,icy,icx) = icbufp(ncc2)
                     tag(iczbp0,icy,icx) = nca
                     nca = nca + na_per_cell(iczbp0,icy,icx)
                     ncar2p = ncar2p + icbufp(ncc2)
                  END DO
               END DO
            END DO

            iczbm0 = 2 + nczdiv + 2
            iczbm1 = iczbm0
            ncc2 = 0
            ncar2m = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  nca = tag(iczbm0-1,icy,icx)
     &                 + na_per_cell(iczbm0-1,icy,icx)
                  ncc2 = ncc2 + 1
                  na_per_cell(iczbm0,icy,icx) = icbufm(ncc2)
                  tag(iczbm0,icy,icx) = nca
                  nca = nca + na_per_cell(iczbm0,icy,icx)
                  ncar2m = ncar2m + icbufm(ncc2)
               END DO
            END DO

#ifdef SYNC_COM
!coarray            call mpi_sendrecv(rbuffp, 3*ncarp, MPI_DOUBLE_PRECISION, 
!coarray     &           ipz_pdest, myrank,
!coarray     &           buffp, 3*ncar2p, MPI_DOUBLE_PRECISION,
!coarray     &           ipz_psrc, ipz_psrc,
!coarray     &           mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(irbuffp, ncarp, MPI_INTEGER, ipz_pdest,
!coarray     &           myrank, ibuffp, ncar2p, MPI_INTEGER,
!coarray     &           ipz_psrc, ipz_psrc,
!coarray     &           mpi_comm_world, istatus, ierr )
!coarray
!coarray            call mpi_sendrecv(rbuffm, 3*ncarm, MPI_DOUBLE_PRECISION, 
!coarray     &             ipz_mdest, myrank,
!coarray     &             buffm, 3*ncar2m, MPI_DOUBLE_PRECISION,
!coarray     &             ipz_msrc, ipz_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(irbuffm, ncarm, MPI_INTEGER,
!coarray     &             ipz_mdest, myrank,
!coarray     &             ibuffm, ncar2m, MPI_INTEGER, ipz_msrc, ipz_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
      !buffp(1:3,1:ncarp)[ipz_pdest+1] = rbuffp(1:3,1:ncarp) ! Put 
      call xmp_array_section_set_triplet(buffp_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(buffp_sec,
     & 2,int(1,kind=8),int(ncarp,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuffp_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuffp_sec,
     & 2,int(1,kind=8),int(ncarp,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,buffp_desc,buffp_sec,
     & rbuffp_desc,rbuffp_sec,status)

      !ibuffp(1:ncarp)[ipz_pdest+1]    = irbuffp(1:ncarp)    ! Put
      call xmp_array_section_set_triplet(ibuffp_sec,
     & 1,int(1,kind=8),int(ncarp,kind=8),1,status)
      call xmp_array_section_set_triplet(irbuffp_sec,
     & 1,int(1,kind=8),int(ncarp,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,ibuffp_desc,ibuffp_sec,
     & irbuffp_desc,irbuffp_sec,status)

      !buffm(1:3,1:ncarm)[ipz_mdest+1] = rbuffm(1:3,1:ncarm) ! Put
      call xmp_array_section_set_triplet(buffm_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(buffm_sec,
     & 2,int(1,kind=8),int(ncarm,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuffm_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(rbuffm_sec,
     & 2,int(1,kind=8),int(ncarm,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,buffm_desc,buffm_sec,
     & rbuffm_desc,rbuffm_sec,status)

      !ibuffm(1:ncarm)[ipz_mdest+1]    = irbuffm(1:ncarm)    ! Put
      call xmp_array_section_set_triplet(ibuffm_sec,
     & 1,int(1,kind=8),int(ncarm,kind=8),1,status)
      call xmp_array_section_set_triplet(irbuffm_sec,
     & 1,int(1,kind=8),int(ncarm,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,ibuffm_desc,ibuffm_sec,
     & irbuffm_desc,irbuffm_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
            call mpi_irecv(buffp, 3*ncar2p, 
     &             MPI_DOUBLE_PRECISION, ipz_psrc, ipz_psrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(rbuffp, 3*ncarp, 
     &             MPI_DOUBLE_PRECISION, ipz_pdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(ibuffp, ncar2p, 
     &             MPI_INTEGER, ipz_psrc, ipz_psrc,
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(irbuffp, ncarp, 
     &             MPI_INTEGER, ipz_pdest, myrank,
     &             mpi_comm_world, irq(4), ierr)

            call mpi_irecv(buffm, 3*ncar2m, 
     &             MPI_DOUBLE_PRECISION, ipz_msrc, ipz_msrc, 
     &             mpi_comm_world, irq(5), ierr)
            call mpi_isend(rbuffm, 3*ncarm, 
     &             MPI_DOUBLE_PRECISION, ipz_mdest, myrank, 
     &             mpi_comm_world, irq(6), ierr)
            call mpi_irecv(ibuffm, ncar2m, 
     &             MPI_INTEGER, ipz_msrc, ipz_msrc,
     &             mpi_comm_world, irq(7), ierr)
            call mpi_isend(irbuffm, ncarm, 
     &             MPI_INTEGER, ipz_mdest, myrank,
     &             mpi_comm_world, irq(8), ierr)

            nrq = 8
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            nca = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO ica = tag(iczbp0, icy, icx), tag(iczbp1, icy, icx)
     &                 + na_per_cell(iczbp1, icy, icx)-1
                     nca = nca + 1
                     wkxyz(1,ica) = buffp(1,nca)
                     wkxyz(2,ica) = buffp(2,nca)
                     wkxyz(3,ica) = buffp(3,nca)
                     m2i(ica) = ibuffp(nca)
                  END DO
               END DO
            END DO

            nca = 0
            DO icx = 3, 2+ncxdiv
               DO icy = 3, 2+ncydiv
                  DO ica = tag(iczbm0, icy, icx), tag(iczbm1, icy, icx)
     &                 + na_per_cell(iczbm1, icy, icx)-1
                     nca = nca + 1
                     wkxyz(1,ica) = buffm(1,nca)
                     wkxyz(2,ica) = buffm(2,nca)
                     wkxyz(3,ica) = buffm(3,nca)
                     m2i(ica) = ibuffm(nca)
                  END DO
               END DO
            END DO

         END IF
      END DO

!     coordinate +Y
      ipy_pdest  = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx
      ipy_psrc   = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
!     coordinate -Y
      ipy_mdest  = ipz*npx*npy + mod(ipy-1+1/npy+npy,npy)*npx + ipx
      ipy_msrc   = ipz*npx*npy + mod(ipy+1-1/npy+npy,npy)*npx + ipx

      icz0 = 1
      icz1 = 2 + nczdiv + 2
      nitr = (2 - 1)/ncydiv + 1

      DO icx = 3, 2+ncxdiv
         DO itr = 1, nitr
            if (itr == 1) then
               icyp0 = 2 + ncydiv - 1
               if (ncydiv == 1) icyp0 = icyp0 + 1
               icyp1 = 2 + ncydiv
               nccp = (icz1 - icz0 + 1)*(icyp1 - icyp0 + 1)
               icybp0 = icyp0 - ncydiv
               icybp1 = icyp1 - ncydiv
               icym0 = 3
               icym1 = icym0 + 1
               if (ncydiv == 1) icym1 = icym0
               nccm = (icz1 - icz0 + 1)*(icym1 - icym0 + 1)
               icybm0 = icym0 + ncydiv
               icybm1 = icym1 + ncydiv

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(na_per_cell(icz0,icyp0,icx), nccp,
!coarray     &                MPI_INTEGER, ipy_pdest, myrank,
!coarray     &                na_per_cell(icz0,icybp0,icx), nccp, MPI_INTEGER,
!coarray     &                ipy_psrc, ipy_psrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(na_per_cell(icz0,icym0,icx), nccm,
!coarray     &                MPI_INTEGER, ipy_mdest, myrank,
!coarray     &                na_per_cell(icz0,icybm0,icx), nccm, MPI_INTEGER,
!coarray     &                ipy_msrc, ipy_msrc,
!coarray     &                mpi_comm_world, istatus, ierr )
      nd = abs(icyp1 - icyp0)
!     !    na_per_cell(:, icybp0:icybp0+nd, icx)[ipy_pdest+1]
!     !. = na_per_cell(:,  icyp0:icyp0 +nd, icx) ! Put


      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icybp0,kind=8),int(icybp0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icyp0,kind=8),int(icyp0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)


      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec, 
     & na_per_cell_desc,na_per_cell_l_sec,status)



      nd = abs(icym1 - icym0)
!    !     na_per_cell(:, icybm0:icybm0+nd, icx)[ipy_mdest+1]
!    ! . = na_per_cell(:,  icym0:icym0 +nd, icx) ! Put

      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icybm0,kind=8),int(icybm0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icym0,kind=8),int(icym0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
               call mpi_irecv(na_per_cell(icz0,icybp0,icx), nccp,
     &                 MPI_INTEGER, ipy_psrc, ipy_psrc, mpi_comm_world,
     &                 irq(1), ierr)
               call mpi_isend(na_per_cell(icz0,icyp0,icx), nccp,
     &                 MPI_INTEGER, ipy_pdest, myrank, mpi_comm_world, 
     &                 irq(2), ierr)
               call mpi_irecv(na_per_cell(icz0,icybm0,icx), nccm,
     &                 MPI_INTEGER, ipy_msrc, ipy_msrc, mpi_comm_world,
     &                 irq(3), ierr)
               call mpi_isend(na_per_cell(icz0,icym0,icx), nccm,
     &                 MPI_INTEGER, ipy_mdest, myrank, mpi_comm_world,
     &                 irq(4), ierr)
               nrq = 4
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

               nbase = tag(3, 3, icx) - 3*naline
               if(ncydiv == 1) nbase = nbase + naline
               DO icyb = icybp0, icybp1
                  nbase = nbase + naline
                  nca = nbase - na_per_cell(2,icyb,icx)
     &                 - na_per_cell(1,icyb,icx)
                  DO icz = icz0, icz1
                     tag(icz,icyb,icx) = nca
                     nca = nca + na_per_cell(icz,icyb,icx)
                  END DO
               END DO

               nbase = tag(3,2+ncydiv,icx)
               DO icyb = icybm0, icybm1
                  nbase = nbase + naline
                  nca = nbase - na_per_cell(2,icyb,icx)
     &                 - na_per_cell(1,icyb,icx)
                  DO icz = icz0, icz1
                     tag(icz,icyb,icx) = nca
                     nca = nca + na_per_cell(icz,icyb,icx)
                  END DO
               END DO

               ncap = naline * 2
               if (ncydiv == 1) ncap = naline
               icasp = tag(3,icyp0,icx) - 2*na1cell
               icarp = tag(3,icybp0,icx) - 2*na1cell
               ncam = naline * 2
               if (ncydiv == 1) ncam = naline
               icasm = tag(3,icym0,icx) - 2*na1cell
               icarm = tag(3,icybm0,icx) - 2*na1cell

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(wkxyz(1,icasp),3*ncap,
!coarray     &                MPI_DOUBLE_PRECISION, ipy_pdest, myrank,
!coarray     &                wkxyz(1,icarp), 3*ncap, MPI_DOUBLE_PRECISION,
!coarray     &                ipy_psrc, ipy_psrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(m2i(icasp), ncap, MPI_INTEGER,
!coarray     &                ipy_pdest, myrank,
!coarray     &                m2i(icarp), ncap, MPI_INTEGER,
!coarray     &                ipy_psrc, ipy_psrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(wkxyz(1,icasm), 3*ncam,
!coarray     &                MPI_DOUBLE_PRECISION, ipy_mdest, myrank,
!coarray     &                wkxyz(1,icarm), 3*ncam, MPI_DOUBLE_PRECISION, 
!coarray     &                ipy_msrc, ipy_msrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(m2i(icasm), ncam, MPI_INTEGER,
!coarray     &                ipy_mdest, myrank,
!coarray     &                m2i(icarm), ncam, MPI_INTEGER, ipy_msrc, ipy_msrc,
!coarray     &                mpi_comm_world, istatus, ierr )

!     !    wkxyz(:,icarp:icarp+ncap-1)[ipy_pdest+1]
!     !. = wkxyz(:,icasp:icasp+ncap-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)


!     !    m2i(icarp:icarp+ncap-1)[ipy_pdest+1]
!     !. = m2i(icasp:icasp+ncap-1)     ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)



!      sync all
      call xmp_sync_all(status)
!     !    wkxyz(:,icarm:icarm+ncam-1)[ipy_mdest+1]
!     !. = wkxyz(:,icasm:icasm+ncam-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)



!     !    m2i(icarm:icarm+ncam-1)[ipy_mdest+1]
!     !. = m2i(icasm:icasm+ncam-1)     ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)

!      sync all
      call xmp_sync_all(status)

!!
#else
               call mpi_irecv(wkxyz(1,icarp), 3*ncap, 
     &                MPI_DOUBLE_PRECISION, ipy_psrc, ipy_psrc, 
     &                mpi_comm_world, irq(1), ierr)
               call mpi_isend(wkxyz(1,icasp), 3*ncap, 
     &                MPI_DOUBLE_PRECISION, ipy_pdest, myrank, 
     &                mpi_comm_world, irq(2), ierr)
               call mpi_irecv(m2i(icarp), ncap, 
     &                MPI_INTEGER, ipy_psrc, ipy_psrc,
     &                mpi_comm_world, irq(3), ierr)
               call mpi_isend(m2i(icasp), ncap, 
     &                MPI_INTEGER, ipy_pdest, myrank,
     &                mpi_comm_world, irq(4), ierr)

               call mpi_irecv(wkxyz(1,icarm), 3*ncam, 
     &                MPI_DOUBLE_PRECISION, ipy_msrc, ipy_msrc, 
     &                mpi_comm_world, irq(5), ierr)
               call mpi_isend(wkxyz(1,icasm), 3*ncam, 
     &                MPI_DOUBLE_PRECISION, ipy_mdest, myrank, 
     &                mpi_comm_world, irq(6), ierr)
               call mpi_irecv(m2i(icarm), ncam, 
     &                MPI_INTEGER, ipy_msrc, ipy_msrc,
     &                mpi_comm_world, irq(7), ierr)
               call mpi_isend(m2i(icasm), ncam, 
     &                MPI_INTEGER, ipy_mdest, myrank,
     &                mpi_comm_world, irq(8), ierr)
               nrq = 8
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            else
               icybp1st = icybp0
               icybp0 = 1
               icybp1 = 1
               icybm1st = icybm0
               icybm0 = 2 + ncydiv + 2
               icybm1 = icybm0

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(na_per_cell(icz0,icybp1st,icx), nccp,
!coarray     &                MPI_INTEGER, ipy_pdest, myrank,
!coarray     &                na_per_cell(icz0,icybp0,icx), nccp, MPI_INTEGER,
!coarray     &                ipy_psrc, ipy_psrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(na_per_cell(icz0,icybm1st,icx), nccm,
!coarray     &                MPI_INTEGER, ipy_mdest, myrank,
!coarray     &                na_per_cell(icz0,icybm0,icx), nccm, MPI_INTEGER,
!coarray     &                ipy_msrc, ipy_msrc,
!coarray     &                mpi_comm_world, istatus, ierr )
      nd = abs(icyp1 - icyp0)

!     !    na_per_cell(:,   icybp0:icybp0+nd,   icx)[ipy_pdest+1]
!     !. = na_per_cell(:, icybp1st:icybp1st+nd, icx) ! Put

      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icybp0,kind=8),int(icybp0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icybp1st,kind=8),int(icybp1st+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)

!      sync all
      call xmp_sync_all(status)
      nd = abs(icym1 - icym0)
!     !    na_per_cell(:,   icybm0:icybm0+nd,   icx)[ipy_pdest+1]
!     !. = na_per_cell(:, icybm1st:icybm1st+nd, icx) ! Put

      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icybm0,kind=8),int(icybm0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(1,kind=8),int(lzdiv+4,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icybm1st,kind=8),int(icybm1st+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
               call mpi_irecv(na_per_cell(icz0,icybp0,icx), nccp,
     &                 MPI_INTEGER, ipy_psrc, ipy_psrc, mpi_comm_world,
     &                 irq(1), ierr)
               call mpi_isend(na_per_cell(icz0,icybp1st,icx), nccp,
     &                 MPI_INTEGER, ipy_pdest, myrank, mpi_comm_world, 
     &                 irq(2), ierr)
               call mpi_irecv(na_per_cell(icz0,icybm0,icx), nccm,
     &                 MPI_INTEGER, ipy_msrc, ipy_msrc, mpi_comm_world,
     &                 irq(3), ierr)
               call mpi_isend(na_per_cell(icz0,icybm1st,icx), nccm,
     &                 MPI_INTEGER, ipy_mdest, myrank, mpi_comm_world, 
     &                 irq(4), ierr)
               nrq = 4
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

               nbase = tag(3, 3, icx) - 2*naline
               nca = nbase - na_per_cell(2,icybp0,icx)
     &              - na_per_cell(1,icybp0,icx)
               DO icz = icz0, icz1
                  tag(icz,icybp0,icx) = nca
                  nca = nca + na_per_cell(icz,icybp0,icx)
               END DO

               nbase = tag(3,2+ncydiv,icx) + 2*naline
               nca = nbase - na_per_cell(2,icybm0,icx)
     &              - na_per_cell(1,icybm0,icx)
               DO icz = icz0, icz1
                  tag(icz,icybm0,icx) = nca
                  nca = nca + na_per_cell(icz,icybm0,icx)
               END DO

               ncap = naline
               icasp = tag(3,icybp1st,icx) - 2*na1cell
               icarp = tag(3,icybp0,icx) - 2*na1cell
               ncam = naline
               icasm = tag(3,icybm1st,icx) - 2*na1cell
               icarm = tag(3,icybm0,icx) - 2*na1cell

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(wkxyz(1,icasp), 3*ncap,
!coarray     &                MPI_DOUBLE_PRECISION, ipy_pdest, myrank,
!coarray     &                wkxyz(1,icarp), 3*ncap, MPI_DOUBLE_PRECISION, 
!coarray     &                ipy_psrc, ipy_psrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(m2i(icasp), ncap, MPI_INTEGER,
!coarray     &                ipy_pdest, myrank,
!coarray     &                m2i(icarp), ncap, MPI_INTEGER,ipy_psrc,ipy_psrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(wkxyz(1,icasm), 3*ncam,
!coarray     &                MPI_DOUBLE_PRECISION, ipy_mdest, myrank,
!coarray     &                wkxyz(1,icarm), 3*ncam, MPI_DOUBLE_PRECISION, 
!coarray     &                ipy_msrc, ipy_msrc,
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(m2i(icasm), ncam, MPI_INTEGER,
!coarray     &                ipy_mdest, myrank,
!coarray     &                m2i(icarm), ncam, MPI_INTEGER,ipy_msrc,ipy_msrc,
!coarray     &                mpi_comm_world, istatus, ierr )


!     !    wkxyz(:,icarp:icarp+ncap-1)[ipy_pdest+1]
!     !. = wkxyz(:,icasp:icasp+ncap-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)

!     !    m2i(icarp:icarp+ncap-1)[ipy_pdest+1]
!     !. = m2i(icasp:icasp+ncap-1)     ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)


!      sync all
      call xmp_sync_all(status)

!     !    wkxyz(:,icarm:icarm+ncam-1)[ipy_mdest+1]
!     !. = wkxyz(:,icasm:icasm+ncam-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)



!     !    m2i(icarm:icarm+ncam-1)[ipy_mdest+1]
!     !. = m2i(icasm:icasm+ncam-1)     ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)

!      sync all
      call xmp_sync_all(status)


!!
#else
               call mpi_irecv(wkxyz(1,icarp), 3*ncap, 
     &                MPI_DOUBLE_PRECISION, ipy_psrc, ipy_psrc, 
     &                mpi_comm_world, irq(1), ierr)
               call mpi_isend(wkxyz(1,icasp), 3*ncap, 
     &                MPI_DOUBLE_PRECISION, ipy_pdest, myrank, 
     &                mpi_comm_world, irq(2), ierr)
               call mpi_irecv(m2i(icarp), ncap, 
     &                MPI_INTEGER, ipy_psrc, ipy_psrc,
     &                mpi_comm_world, irq(3), ierr)
               call mpi_isend(m2i(icasp), ncap, 
     &                MPI_INTEGER, ipy_pdest, myrank,
     &                mpi_comm_world, irq(4), ierr)

               call mpi_irecv(wkxyz(1,icarm), 3*ncam, 
     &                MPI_DOUBLE_PRECISION, ipy_msrc, ipy_msrc, 
     &                mpi_comm_world, irq(5), ierr)
               call mpi_isend(wkxyz(1,icasm), 3*ncam, 
     &                MPI_DOUBLE_PRECISION, ipy_mdest, myrank, 
     &                mpi_comm_world, irq(6), ierr)
               call mpi_irecv(m2i(icarm), ncam, 
     &                MPI_INTEGER, ipy_msrc, ipy_msrc,
     &                mpi_comm_world, irq(7), ierr)
               call mpi_isend(m2i(icasm), ncam, 
     &                MPI_INTEGER, ipy_mdest, myrank,
     &                mpi_comm_world, irq(8), ierr)

               nrq = 8
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            END IF
         END DO
      END DO

#ifndef HALFDIREE
!     coordinate +X
      ipx_pdest  = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
      ipx_psrc   = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
#endif
!     coordinate -X
      ipx_mdest  = ipz*npx*npy + ipy*npx + mod(ipx-1+1/npx+npx,npx)
      ipx_msrc   = ipz*npx*npy + ipy*npx + mod(ipx+1-1/npx+npx,npx)
      icz0 = 1
      icz1 = 2 + nczdiv + 2
      icy0 = 1
      icy1 = 2 + ncydiv + 2
      nitr = (2 - 1)/ncxdiv + 1

      DO itr = 1, nitr
         if (itr == 1) then

#ifndef HALFDIREE
            icxp0 = 2 + ncxdiv - 1
            if (ncxdiv == 1) icxp0 = icxp0 + 1
            icxp1 = 2 + ncxdiv
            
            nccp = (icz1 - icz0 +1)*(icy1 - icy0 +1)*(icxp1 - icxp0 +1)
            icxbp0 = icxp0 - ncxdiv
            icxbp1 = icxp1 - ncxdiv
#endif
            icxm0 = 3
            icxm1 = icxm0 + 1
            if (ncxdiv == 1) icxm1 = icxm0

            nccm = (icz1 - icz0 +1)*(icy1 - icy0 +1)*(icxm1 - icxm0 +1)
            icxbm0 = icxm0 + ncxdiv
            icxbm1 = icxm1 + ncxdiv

#ifdef SYNC_COM
#ifndef HALFDIREE
!coarray            call mpi_sendrecv(na_per_cell(icz0,icy0,icxp0), nccp,
!coarray     &             MPI_INTEGER, ipx_pdest, myrank,
!coarray     &             na_per_cell(icz0,icy0,icxbp0),
!coarray     &             nccp, MPI_INTEGER, ipx_psrc, ipx_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!     ! na_per_cell(icz0:icz1,icy0:icy1,icxbp0:icxbp0+(icxp1-icxp0))
!     !. [ipx_pdest+1]
!     !. = na_per_cell(icz0:icz1,icy0:icy1,icxp0:icxp1) ! Put
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icxbp0,kind=8),int(icxbp0+(icxp1-icxp0),kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icxp0,kind=8),int(icxp1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)


!      sync all
      call xmp_sync_all(status)

!!
#endif
!coarray            call mpi_sendrecv(na_per_cell(icz0,icy0,icxm0), nccm,
!coarray     &             MPI_INTEGER, ipx_mdest, myrank,
!coarray     &             na_per_cell(icz0,icy0,icxbm0), nccm, MPI_INTEGER,
!coarray     &             ipx_msrc, ipx_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!     ! na_per_cell(icz0:icz1,icy0:icy1,icxbm0:icxbm0+(icxm1-icxm0))
!     !. [ipx_mdest+1]
!     !. = na_per_cell(icz0:icz1,icy0:icy1,icxm0:icxm1) ! Put

      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icxbm0,kind=8),int(icxbm0+(icxm1-icxm0),kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icxm0,kind=8),int(icxm1,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
#ifndef HALFDIREE
            call mpi_irecv(na_per_cell(icz0,icy0,icxbp0), nccp,
     &              MPI_INTEGER, ipx_psrc, ipx_psrc, mpi_comm_world,
     &              irq(1), ierr)
            call mpi_isend(na_per_cell(icz0,icy0,icxp0), nccp,
     &              MPI_INTEGER, ipx_pdest, myrank, mpi_comm_world, 
     &              irq(2), ierr)
            call mpi_irecv(na_per_cell(icz0,icy0,icxbm0), nccm,
     &              MPI_INTEGER, ipx_msrc, ipx_msrc, mpi_comm_world,
     &              irq(3), ierr)
            call mpi_isend(na_per_cell(icz0,icy0,icxm0), nccm,
     &              MPI_INTEGER, ipx_mdest, myrank, mpi_comm_world, 
     &              irq(4), ierr)
            nrq = 4
#else
            call mpi_irecv(na_per_cell(icz0,icy0,icxbm0), nccm,
     &              MPI_INTEGER, ipx_msrc, ipx_msrc, mpi_comm_world,
     &              irq(1), ierr)
            call mpi_isend(na_per_cell(icz0,icy0,icxm0), nccm,
     &              MPI_INTEGER, ipx_mdest, myrank, mpi_comm_world, 
     &              irq(2), ierr)
            nrq = 2
#endif

            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifndef HALFDIREE
            nbase3 = tag(3,3,3) - 3*narea - 3*naline
            if(ncxdiv == 1) nbase3 = nbase3 + narea
            DO icxb = icxbp0, icxbp1
               nbase3 = nbase3 + narea
               nbase2 = nbase3
               DO icy = icy0, icy1
                  nbase2 = nbase2 + naline
                  nca = nbase2 - na_per_cell(2,icy,icxb)
     &                 - na_per_cell(1,icy,icxb)
                  DO icz = icz0, icz1
                     tag(icz,icy,icxb) = nca
                     nca = nca + na_per_cell(icz,icy,icxb)
                  END DO
               END DO
            END DO
#endif

            nbase3 = tag(3,3,2+ncxdiv) - 3*naline
            DO icxb = icxbm0, icxbm1
               nbase3 = nbase3 + narea
               nbase2 = nbase3
               DO icy = icy0, icy1
                  nbase2 = nbase2 + naline
                  nca = nbase2 - na_per_cell(2,icy,icxb)
     &                 - na_per_cell(1,icy,icxb)
                  DO icz = icz0, icz1
                     tag(icz,icy,icxb) = nca
                     nca = nca + na_per_cell(icz,icy,icxb)
                  END DO
               END DO
            END DO

#ifndef HALFDIREE
            ncap = narea * 2
            if (ncxdiv == 1) ncap = narea
            icasp = tag(3,icy0,icxp0) - 2*na1cell
            icarp = tag(3,icy0,icxbp0) - 2*na1cell
#endif
            ncam = narea * 2
            if (ncxdiv == 1) ncam = narea
            icasm = tag(3,icy0,icxm0) - 2*na1cell
            icarm = tag(3,icy0,icxbm0) - 2*na1cell

#ifdef SYNC_COM
#ifndef HALFDIREE
!coarray            call mpi_sendrecv(wkxyz(1,icasp), 3*ncap, 
!coarray     &             MPI_DOUBLE_PRECISION, ipx_pdest, myrank,
!coarray     &             wkxyz(1,icarp), 3*ncap, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_psrc, ipx_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(m2i(icasp), ncap, MPI_INTEGER,
!coarray     &             ipx_pdest, myrank,
!coarray     &             m2i(icarp), ncap, MPI_INTEGER, ipx_psrc, ipx_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!     !    wkxyz(:,icarp:icarp+ncap-1)[ipx_pdest+1]
!     !. = wkxyz(:,icasp:icasp+ncap-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)

!     !    m2i(icarp:icarp+ncap-1)[ipx_pdest+1]
!     !. = m2i(icasp:icasp+ncap-1)     ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#endif
!coarray            call mpi_sendrecv(wkxyz(1,icasm), 3*ncam, 
!coarray     &             MPI_DOUBLE_PRECISION, ipx_mdest, myrank,
!coarray     &             wkxyz(1,icarm), 3*ncam, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_msrc, ipx_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(m2i(icasm), ncam, MPI_INTEGER,
!coarray     &             ipx_mdest, myrank,
!coarray     &             m2i(icarm), ncam, MPI_INTEGER, ipx_msrc, ipx_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!     !    wkxyz(:,icarm:icarm+ncam-1)[ipx_mdest+1]
!     !. = wkxyz(:,icasm:icasm+ncam-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasm,kind=8),int(icasm+ncap-1,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)

!     !    m2i(icarm:icarm+ncam-1)[ipx_mdest+1]
!     !. = m2i(icasm:icasm+ncam-1)     ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)
!      sync all
      call xmp_sync_all(status)

!!
#else
#ifndef HALFDIREE
            call mpi_irecv(wkxyz(1,icarp), 3*ncap, 
     &             MPI_DOUBLE_PRECISION, ipx_psrc, ipx_psrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(wkxyz(1,icasp), 3*ncap, 
     &             MPI_DOUBLE_PRECISION, ipx_pdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(m2i(icarp), ncap, 
     &             MPI_INTEGER, ipx_psrc, ipx_psrc,
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(m2i(icasp), ncap, 
     &             MPI_INTEGER, ipx_pdest, myrank,
     &             mpi_comm_world, irq(4), ierr)

            call mpi_irecv(wkxyz(1,icarm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_msrc, ipx_msrc, 
     &             mpi_comm_world, irq(5), ierr)
            call mpi_isend(wkxyz(1,icasm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_mdest, myrank, 
     &             mpi_comm_world, irq(6), ierr)
            call mpi_irecv(m2i(icarm), ncam, 
     &             MPI_INTEGER, ipx_msrc, ipx_msrc,
     &             mpi_comm_world, irq(7), ierr)
            call mpi_isend(m2i(icasm), ncam, 
     &             MPI_INTEGER, ipx_mdest, myrank,
     &             mpi_comm_world, irq(8), ierr)

            nrq = 8
#else
            call mpi_irecv(wkxyz(1,icarm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_msrc, ipx_msrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(wkxyz(1,icasm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_mdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(m2i(icarm), ncam, 
     &             MPI_INTEGER, ipx_msrc, ipx_msrc,
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(m2i(icasm), ncam, 
     &             MPI_INTEGER, ipx_mdest, myrank,
     &             mpi_comm_world, irq(4), ierr)

            nrq = 4
#endif
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

         else
#ifndef HALFDIREE
            icxbp1st = icxbp0
            icxbp0 = 1
            icxbp1 = 1
#endif
            icxbm1st = icxbm0
            icxbm0 = 2 + ncxdiv + 2
            icxbm1 = icxbm0

#ifdef SYNC_COM
#ifndef HALFDIREE
!coarray            call mpi_sendrecv(na_per_cell(icz0,icy0,icxbp1st), nccp,
!coarray     &           MPI_INTEGER, ipx_pdest, myrank,
!coarray     &           na_per_cell(icz0,icy0,icxbp0), nccp, MPI_INTEGER,
!coarray     &           ipx_psrc, ipx_psrc,
!coarray     &           mpi_comm_world, istatus, ierr )
      nd = abs(icxp1 - icxp0)
!     ! na_per_cell(icz0:icz1,icy0:icy1,icxbp0:icxbp0+nd)[ipx_pdest+1]
!     !.= na_per_cell(icz0:icz1,icy0:icy1,icxbp1st:icxbp1st+nd)

      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icxbp0,kind=8),int(icxbp0+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icxbp1st,kind=8),int(icxbp1st+nd,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)




!      sync all
      call xmp_sync_all(status)
!!
#endif
!coarray            call mpi_sendrecv(na_per_cell(icz0,icy0,icxbm1st), nccm,
!coarray     &             MPI_INTEGER, ipx_mdest, myrank,
!coarray     &             na_per_cell(icz0,icy0,icxbm0), nccm, MPI_INTEGER,
!coarray     &             ipx_msrc, ipx_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
      nd = abs(icxm1 - icxm0)

!     ! na_per_cell(icz0:icz1,icy0:icy1,icxbm0:icxbm0+nd)[ipx_pdest+1]
!     !.= na_per_cell(icz0:icz1,icy0:icy1,icxbm1st:icxbm1st+nd)

      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_r_sec,
     & 3,int(icxbm0,kind=8),int(icxbm0+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 1,int(icz0,kind=8),int(icz1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 2,int(icy0,kind=8),int(icy1,kind=8),1,status)
      call xmp_array_section_set_triplet(na_per_cell_l_sec,
     & 3,int(icxbm1st,kind=8),int(icxbm1st+nd,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,
     & na_per_cell_desc,na_per_cell_r_sec,
     & na_per_cell_desc,na_per_cell_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
#ifndef HALFDIREE
            call mpi_irecv(na_per_cell(icz0,icy0,icxbp0), nccp,
     &             MPI_INTEGER, ipx_psrc, ipx_psrc, mpi_comm_world, 
     &             irq(1), ierr)
            call mpi_isend(na_per_cell(icz0,icy0,icxbp1st), nccp,
     &             MPI_INTEGER, ipx_pdest, myrank, mpi_comm_world, 
     &             irq(2), ierr)
            call mpi_irecv(na_per_cell(icz0,icy0,icxbm0), nccm, 
     &             MPI_INTEGER, ipx_msrc, ipx_msrc, mpi_comm_world, 
     &             irq(3), ierr)
            call mpi_isend(na_per_cell(icz0,icy0,icxbm1st), nccm,
     &             MPI_INTEGER, ipx_mdest, myrank,  mpi_comm_world, 
     &             irq(4), ierr)
            nrq = 4
#else
            call mpi_irecv(na_per_cell(icz0,icy0,icxbm0), nccm, 
     &             MPI_INTEGER, ipx_msrc, ipx_msrc, mpi_comm_world, 
     &             irq(1), ierr)
            call mpi_isend(na_per_cell(icz0,icy0,icxbm1st), nccm,
     &             MPI_INTEGER, ipx_mdest, myrank,  mpi_comm_world, 
     &             irq(2), ierr)
            nrq = 2
#endif
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

#ifndef HALFDIREE
            nbase2 = tag(3,3,3) -2*narea - 3*naline
            DO icy = icy0, icy1
               nbase2 = nbase2 + naline
               nca = nbase2 - na_per_cell(2,icy,icxbp0)
     &              - na_per_cell(1,icy,icxbp0)
               DO icz = icz0, icz1
                  tag(icz,icy,icxbp0) = nca
                  nca = nca + na_per_cell(icz,icy,icxbp0)
               END DO
            END DO
#endif
            nbase2 = tag(3,3,2+ncxdiv) + 2*narea - 3*naline
            DO icy = icy0, icy1
               nbase2 = nbase2 + naline
               nca = nbase2 - na_per_cell(2,icy,icxbm0)
     &              - na_per_cell(1,icy,icxbm0)
               DO icz = icz0, icz1
                  tag(icz,icy,icxbm0) = nca
                  nca = nca + na_per_cell(icz,icy,icxbm0)
               END DO
            END DO

#ifndef HALFDIREE
            ncap = narea
            icasp = tag(3,icy0,icxbp1st) - 2*na1cell
            icarp = tag(3,icy0,icxbp0) - 2*na1cell
#endif
            ncam = narea
            icasm = tag(3,icy0,icxbm1st) - 2*na1cell
            icarm = tag(3,icy0,icxbm0) - 2*na1cell

#ifdef SYNC_COM
#ifndef HALFDIREE
!coarray            call mpi_sendrecv(wkxyz(1,icasp), 3*ncap, 
!coarray     &             MPI_DOUBLE_PRECISION, ipx_pdest, myrank,
!coarray     &             wkxyz(1,icarp), 3*ncap, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_psrc, ipx_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(m2i(icasp), ncap, MPI_INTEGER,
!coarray     &             ipx_pdest, myrank,
!coarray     &             m2i(icarp), ncap, MPI_INTEGER, ipx_psrc, ipx_psrc,
!coarray     &             mpi_comm_world, istatus, ierr )

!     ! wkxyz(:,icarp:icarp+ncap-1)[ipx_pdest+1]
!     !. = wkxyz(:,icasp:icasp+ncap-1) ! Put

      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)



!     ! m2i(icarp:icarp+ncap-1)[ipx_pdest+1]
!     !. = m2i(icasp:icasp+ncap-1) ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarp,kind=8),int(icarp+ncap-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasp,kind=8),int(icasp+ncap-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#endif
!coarray            call mpi_sendrecv(wkxyz(1,icasm), 3*ncam, 
!coarray     &             MPI_DOUBLE_PRECISION, ipx_mdest, myrank,
!coarray     &             wkxyz(1,icarm), 3*ncam, MPI_DOUBLE_PRECISION, 
!coarray     &             ipx_msrc, ipx_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!coarray            call mpi_sendrecv(m2i(icasm), ncam, MPI_INTEGER,
!coarray     &             ipx_mdest, myrank,
!coarray     &             m2i(icarm), ncam, MPI_INTEGER, ipx_msrc, ipx_msrc,
!coarray     &             mpi_comm_world, istatus, ierr )
!     ! wkxyz(:,icarm:icarm+ncam-1)[ipx_mdest+1]
!     !. = wkxyz(:,icasm:icasm+ncam-1) ! Put
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_r_sec,
     & 2,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 1,int(1,kind=8),int(3,kind=8),1,status)
      call xmp_array_section_set_triplet(wkxyz_l_sec,
     & 2,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,wkxyz_desc,wkxyz_r_sec,
     & wkxyz_desc,wkxyz_l_sec,status)


!     ! m2i(icarm:icarm+ncam-1)[ipx_mdest+1]
!     !. = m2i(icasm:icasm+ncam-1) ! Put
      call xmp_array_section_set_triplet(m2i_r_sec,
     & 1,int(icarm,kind=8),int(icarm+ncam-1,kind=8),1,status)

      call xmp_array_section_set_triplet(m2i_l_sec,
     & 1,int(icasm,kind=8),int(icasm+ncam-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,m2i_desc,m2i_r_sec,
     & m2i_desc,m2i_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
#ifndef HALFDIREE
            call mpi_irecv(wkxyz(1,icarp), 3*ncap, 
     &             MPI_DOUBLE_PRECISION, ipx_psrc, ipx_psrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(wkxyz(1,icasp), 3*ncap, 
     &             MPI_DOUBLE_PRECISION, ipx_pdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(m2i(icarp), ncap, 
     &             MPI_INTEGER, ipx_psrc, ipx_psrc,
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(m2i(icasp), ncap, 
     &             MPI_INTEGER, ipx_pdest, myrank,
     &             mpi_comm_world, irq(4), ierr)

            call mpi_irecv(wkxyz(1,icarm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_msrc, ipx_msrc, 
     &             mpi_comm_world, irq(5), ierr)
            call mpi_isend(wkxyz(1,icasm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_mdest, myrank, 
     &             mpi_comm_world, irq(6), ierr)
            call mpi_irecv(m2i(icarm), ncam, 
     &             MPI_INTEGER, ipx_msrc, ipx_msrc,
     &             mpi_comm_world, irq(7), ierr)
            call mpi_isend(m2i(icasm), ncam, 
     &             MPI_INTEGER, ipx_mdest, myrank,
     &             mpi_comm_world, irq(8), ierr)

            nrq = 8
#else
            call mpi_irecv(wkxyz(1,icarm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_msrc, ipx_msrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(wkxyz(1,icasm), 3*ncam, 
     &             MPI_DOUBLE_PRECISION, ipx_mdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(m2i(icarm), ncam, 
     &             MPI_INTEGER, ipx_msrc, ipx_msrc,
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(m2i(icasm), ncam, 
     &             MPI_INTEGER, ipx_mdest, myrank,
     &             mpi_comm_world, irq(4), ierr)

            nrq = 4
#endif
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

         END IF
      END DO

!=== create i2m ===!
      ntmp=0
!$omp parallel default(none)
!$omp& private(ncc,icx,icy,icz,ica,icag)
!$omp& shared(ncxdiv,ncydiv,nczdiv,tag,na_per_cell)
!$omp& shared(m2i,i2m)
!$omp& reduction(+:ntmp)
!$omp do
#ifndef HALFDIREE
      do ncc=1,(ncxdiv+4)*(ncydiv+4)*(nczdiv+4)
#else
      do ncc=2*(ncydiv+4)*(nczdiv+4)+1,
     &       (ncxdiv+4)*(ncydiv+4)*(nczdiv+4)
#endif
        icz=mod(ncc-1,nczdiv+4)   +1   
        icy=mod(ncc-1,(nczdiv+4)*(ncydiv+4))
        icy=icy/(nczdiv+4)   +1   
        icx=(ncc-1)/((nczdiv+4)*(ncydiv+4))+1

        if(icx.ge.3.and.icx.le.ncxdiv+2 .and.
     &     icy.ge.3.and.icy.le.ncydiv+2 .and.  
     &     icz.ge.3.and.icz.le.nczdiv+2) cycle
        do ica=tag(icz,icy,icx),tag(icz,icy,icx)
     &        +na_per_cell(icz,icy,icx)-1
          icag=m2i(ica)
          i2m(icag)=ica
          ntmp=ntmp+1
!         if(icag .le. -1) cycle
        end do ! ica
      end do ! ncc
!$omp end do
!$omp end parallel
      ndatm=nselfatm+ntmp

      call xmp_free_array_section(icbufp_sec)
      call xmp_free_array_section(ircbufp_sec)
      call xmp_free_array_section(icbufm_sec)
      call xmp_free_array_section(ircbufm_sec)
      call xmp_free_array_section(ibuffp_sec)
      call xmp_free_array_section(irbuffp_sec)
      call xmp_free_array_section(ibuffm_sec)
      call xmp_free_array_section(irbuffm_sec)
      call xmp_free_array_section(buffp_sec)
      call xmp_free_array_section(rbuffp_sec)
      call xmp_free_array_section(buffm_sec)
      call xmp_free_array_section(rbuffm_sec)
      call xmp_free_array_section(wkxyz_l_sec)
      call xmp_free_array_section(wkxyz_r_sec)
      call xmp_free_array_section(m2i_l_sec)
      call xmp_free_array_section(m2i_r_sec)
      call xmp_free_array_section(na_per_cell_l_sec)
      call xmp_free_array_section(na_per_cell_r_sec)

      return
      end
