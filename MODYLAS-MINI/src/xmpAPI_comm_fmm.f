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
c---------------------------------------------------------------------
      subroutine comm_fmm_local_top(il0,mylm,wm,nscell,
     $                              nsczdiv, nscydiv,nscxdiv)
!! UL ver.20111004
c---------------------------------------------------------------------
      use comm_base
      use md_fmm
      use mpivar
      use xmp_api
      implicit none
      INCLUDE 'mpif.h'
      integer m1
      integer(4) nscell, nsczdiv, nscydiv, nscxdiv, mylm, il0
      complex(8) :: ccbuf(mylm*5*nscydiv*nscxdiv)
      integer(8) :: ccbuf_local_desc, ccbuf_local_sec
      integer(8), dimension(1) :: ccbuf_lb, ccbuf_ub
!coarray      complex(8) :: rccbuf(mylm*5*nscydiv*nscxdiv,2)
      complex(8) :: wm(mylm,nscell,nscell,nscell)
!coarray
      !complex(8),allocatable :: rccbuf(:,:)[:]
      complex(8), POINTER :: rccbuf(:,:) => null()
      integer(8) :: rccbuf_desc
      integer(8) :: rccbuf_l_sec, rccbuf_r_sec
      integer(8), dimension(2) :: rccbuf_lb,rccbuf_ub
      !complex(8),allocatable :: wm_tmp(:,:,:,:)[:]
      complex(8), POINTER :: wm_tmp(:,:,:,:) => null()
      integer(8) :: wm_tmp_desc
      integer(8) :: wm_tmp_l_sec, wm_tmp_r_sec
      integer(8), dimension(2) :: wm_tmp_lb,wm_tmp_ub
      !integer,allocatable :: ndis(:)[:]
      integer, POINTER :: ndis(:) => null()
      integer(8) :: ndis_desc
      integer(8) :: ndis_sec
      integer(8), dimension(1) :: ndis_lb,ndis_ub

      integer me, np, nb, nd
      integer ierrcode
!!
      integer np_supercellz,np_supercelly,np_supercellx
      integer mcell_size
      integer ipz_dest, ipy_dest, ipx_dest
      integer ipz_src,  ipy_src,  ipx_src
      integer nitr,itr
      integer icz0, icz1, iczb0, iczb1
      integer icy0, icy1, icyb0, icyb1
      integer icx0, icx1, icxb0, icxb1
      integer ncc, ncc2
      integer icx, icy, icz
      integer m
      integer ibs, ibr
      integer iczb
      integer icyb0prior, icxb0prior
      integer ierr,istatus(mpi_status_size)
      integer(4)  status

!coarray
      !allocate( rccbuf(mylm*5*nscydiv*nscxdiv,2)[*] )
      rccbuf_lb(1)=1
      rccbuf_lb(2)=1
      rccbuf_ub(1)=mylm*5*nscydiv*nscxdiv
      rccbuf_ub(2)=2
      call xmp_new_coarray(rccbuf_desc,
     & 8,2,rccbuf_lb,rccbuf_ub,1,img_dims)
      call xmp_coarray_bind(rccbuf_desc,rccbuf)
      call xmp_new_array_section(rccbuf_l_sec,2)
      call xmp_new_array_section(rccbuf_r_sec,2)


      !allocate( wm_tmp(mylm,nscell,nscell,nscell)[*] )
      wm_tmp_lb(1)=1
      wm_tmp_lb(2)=1
      wm_tmp_lb(3)=1
      wm_tmp_lb(4)=1
      wm_tmp_ub(1)=mylm
      wm_tmp_ub(2)=nscell
      wm_tmp_ub(3)=nscell
      wm_tmp_ub(4)=nscell
      call xmp_new_coarray(wm_tmp_desc,
     & 8,4,wm_tmp_lb,wm_tmp_ub,1,img_dims)
      call xmp_coarray_bind(wm_tmp_desc,wm_tmp)
      call xmp_new_array_section(wm_tmp_l_sec,4)
      call xmp_new_array_section(wm_tmp_r_sec,4)


      wm_tmp = wm
      me = xmp_this_image()
      np = xmp_num_images()
      !allocate( ndis(np)[*] )
      ndis_lb(1)=1
      ndis_ub(1)=np
      call xmp_new_coarray(ndis_desc,8,1,ndis_lb,ndis_ub,1,img_dims)
      call xmp_coarray_bind(ndis_desc,ndis)
      call xmp_new_array_section(ndis_sec,1)
!!
      ccbuf_lb(1) = 1
      ccbuf_ub(1) = mylm*5*nscydiv*nscxdiv
      call xmp_new_local_array(ccbuf_local_desc,8,1,
     & ccbuf_lb,ccbuf_ub,loc(ccbuf))
      call xmp_new_array_section(ccbuf_local_sec,1)

!=== local constant ===!
      m1 = (nmax+1)*(nmax+1)
      mcell_size = 2**il0

!=== global constant ===!
      np_supercellz = (npz * mcell_size-1) / ncell + 1
      np_supercelly = (npy * mcell_size-1) / ncell + 1
      np_supercellx = (npx * mcell_size-1) / ncell + 1

      icz0 = (ncell / mcell_size * ipz) / npz + 1
      icz1 = icz0 + nsczdiv - 1
      icy0 = (ncell / mcell_size * ipy) / npy + 1
      icy1 = icy0 + nscydiv - 1
      icx0 = (ncell / mcell_size * ipx) / npx + 1
      icx1 = icx0 + nscxdiv - 1

! ULmoment +Z
      ipz_dest = mod(ipz + np_supercellz - 1/npz + npz, npz)*npy*npx
     &         + ipy*npx + ipx
      ipz_src  = mod(ipz - np_supercellz + 1/npz + npz, npz)*npy*npx
     &         + ipy*npx + ipx
      nitr = nscell / 2
      if(nscell > npz) nitr = npz / 2

      DO itr = 1, nitr
         if (itr==1) then
            iczb0 = mod(icz0 - nsczdiv - 1 + nscell, nscell) + 1
            iczb1 = mod(icz1 - nsczdiv - 1 + nscell, nscell) + 1

            ncc = 0
            DO icx = icx0, icx1
               DO icy = icy0, icy1
                  DO icz = icz0, icz1
                     DO m = 1, m1
                        ncc = ncc + 1
                        ccbuf(ncc) = wm_tmp(m, icz, icy, icx )
                     END DO
                  END DO
               END DO
            END DO

!coarray            call mpi_sendrecv(ccbuf, ncc, MPI_DOUBLE_COMPLEX,ipz_dest,
!coarray     &           myrank, rccbuf(1,1), ncc, MPI_DOUBLE_COMPLEX, 
!coarray     &           ipz_src, ipz_src,mpi_comm_world,istatus,ierr)
      !rccbuf(1:ncc,1)[ipz_dest+1] = ccbuf(1:ncc) ! Put
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 2,int(1,kind=8),int(1,kind=8),1,status)
      call xmp_array_section_set_triplet(ccbuf_local_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,rccbuf_desc,rccbuf_r_sec, 
     & ccbuf_local_desc,ccbuf_local_sec,status)



!      sync all
      call xmp_sync_all(status)
!!

            ncc2 = 0
            DO icx = icx0, icx1
               DO icy = icy0, icy1
                  DO iczb = iczb0, iczb1
                     DO m = 1, m1
                        ncc2 = ncc2 + 1
                        wm_tmp(m,iczb,icy,icx) = rccbuf(ncc2,1)
                     END DO
                  END DO
               END DO
            END DO
         else
            ibs = mod(itr, 2) + 1
            ibr = mod(itr+1, 2) + 1
            iczb0 = mod(icz0 - nsczdiv*itr - 1 + nscell, nscell) + 1
            iczb1 = mod(iczb0 + nsczdiv - 1 - 1 + nscell,nscell) + 1

!debug
!      write(6,*) "ibs,ibr,iczb0,iczb1,ncc,mylm*5*nscydiv*nscxdiv= ",
!     & ibs,ibr,iczb0,iczb1,ncc,mylm*5*nscydiv*nscxdiv
!      call mpi_abort(mpi_comm_world,ierrcode,ierr)
!!
!coarray            call mpi_sendrecv(rccbuf(1,ibs), ncc, MPI_DOUBLE_COMPLEX,
!coarray     &           ipz_dest,
!coarray     &           myrank, rccbuf(1,ibr), ncc, MPI_DOUBLE_COMPLEX,
!coarray     &           ipz_src, ipz_src,mpi_comm_world,istatus,ierr)
      !rccbuf(1:ncc,ibr)[ipz_dest+1] = rccbuf(1:ncc,ibs) ! Put
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 2,int(ibr,kind=8),int(ibr,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_l_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_l_sec,
     & 2,int(ibs,kind=8),int(ibs,kind=8),1,status)
      img_dims(1) = ipz_dest+1
      call xmp_coarray_put(img_dims,rccbuf_desc,rccbuf_r_sec, 
     & rccbuf_desc,rccbuf_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!

            ncc2 = 0
            DO icx = icx0, icx1
               DO icy = icy0, icy1
                  DO iczb = iczb0, iczb1
                     DO m = 1, m1
                        ncc2 = ncc2 + 1
                        wm_tmp(m,iczb,icy,icx) = rccbuf(ncc2,ibr)
                     END DO
                  END DO
               END DO
            END DO
         end if
      END DO

! ULmoment -Z
      ipz_dest = mod(ipz - np_supercellz + 1/npz + npz, npz)*npy*npx
     &         + ipy*npx + ipx
      ipz_src  = mod(ipz + np_supercellz - 1/npz + npz, npz)*npy*npx
     &         + ipy*npx + ipx
      nitr = (nscell - 1) / 2
      if(nscell > npz) nitr = (npz - 1) / 2

      DO itr = 1, nitr
         if (itr==1) then
            iczb0 = mod(icz0 + nsczdiv - 1, nscell) + 1
            iczb1 = mod(icz1 + nsczdiv - 1, nscell) + 1

!coarray            call mpi_sendrecv(ccbuf, ncc, MPI_DOUBLE_COMPLEX,
!coarray     &           ipz_dest, myrank,
!coarray     &           rccbuf(1,1), ncc, MPI_DOUBLE_COMPLEX,
!coarray     &           ipz_src, ipz_src,mpi_comm_world,istatus,ierr)
      !rccbuf(1:ncc,1)[ipz_dest+1] = ccbuf(1:ncc) ! Put
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 2,int(1,kind=8),int(1,kind=8),1,status)
      call xmp_array_section_set_triplet(ccbuf_local_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      img_dims(1) = ipz_dest+1
      call xmp_coarray_put_local(img_dims,rccbuf_desc,rccbuf_r_sec, 
     & ccbuf_local_desc,ccbuf_local_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
            ncc2 = 0
            DO icx = icx0, icx1
               DO icy = icy0, icy1
                  DO iczb = iczb0, iczb1
                     DO m = 1, m1
                        ncc2 = ncc2 + 1
                        wm_tmp(m,iczb,icy,icx) = rccbuf(ncc2,1)
                     END DO
                  END DO
               END DO
            END DO
         else
            ibs = mod(itr, 2) + 1
            ibr = mod(itr+1, 2) + 1
            iczb0 = mod(icz0 + nsczdiv * itr - 1, nscell) + 1
            iczb1 = mod(iczb0 + nsczdiv - 1 - 1, nscell) + 1

!coarray            call mpi_sendrecv(rccbuf(1,ibs), ncc, MPI_DOUBLE_COMPLEX,
!coarray     &           ipz_dest,
!coarray     &           myrank, rccbuf(1,ibr), ncc, MPI_DOUBLE_COMPLEX,
!coarray     &           ipz_src, ipz_src,mpi_comm_world,istatus,ierr)
      !rccbuf(1:ncc,ibr)[ipz_dest+1] = rccbuf(1:ncc,ibs) ! Put
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_r_sec,
     & 2,int(ibr,kind=8),int(ibr,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_l_sec,
     & 1,int(1,kind=8),int(ncc,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbuf_l_sec,
     & 2,int(ibs,kind=8),int(ibs,kind=8),1,status)
      img_dims(1) = ipz_dest+1
      call xmp_coarray_put(img_dims,rccbuf_desc,rccbuf_r_sec, 
     & rccbuf_desc,rccbuf_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
            ncc2 = 0
            DO icx = icx0, icx1
               DO icy = icy0, icy1
                  DO iczb = iczb0, iczb1
                     DO m = 1, m1
                        ncc2 = ncc2 + 1
                        wm_tmp(m,iczb,icy,icx) = rccbuf(ncc2,ibr)
                     END DO
                  END DO
               END DO
            END DO
         end if
      END DO

! ULmoment +Y
      ipy_dest = ipz*npy*npx
     &     + mod(ipy + np_supercelly - 1/npy + npy, npy)*npx + ipx
      ipy_src = ipz*npy*npx
     &     + mod(ipy - np_supercelly + 1/npy + npy, npy)*npx + ipx
      nitr = nscell / 2
      if(nscell > npy) nitr = npy / 2

      DO icx = icx0, icx1
         DO itr = 1, nitr
            if (itr == 1) then
               icyb0 = mod(icy0 - nscydiv - 1 + nscell, nscell) + 1
               icyb1 = mod(icy1 - nscydiv - 1 + nscell, nscell) + 1
               ncc = m1 * nscell * (icyb1 - icyb0 + 1)
!coarray               call mpi_sendrecv(wm_tmp(1,1,icy0,icx), ncc, 
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_dest, myrank, wm_tmp(1,1,icyb0,icx), ncc, 
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_src, ipy_src,mpi_comm_world,istatus,ierr )
      nd = abs(icyb1 - icyb0)
!      ndis(me)[ipy_src+1] = icyb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icyb0,status)

!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipy_dest+1)
!         wm_tmp( :, :,   nb:nb  +nd, icx )[ipy_dest+1]
!     . = wm_tmp( :, :, icy0:icy0+nd, icx ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icy0,kind=8),int(icy0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
            else
               icyb0prior = icyb0
               icyb0 = mod(icy0 - nscydiv * itr - 1 + nscell,
     &                     nscell) + 1
!coarray               call mpi_sendrecv(wm_tmp(1,1,icyb0prior,icx), ncc,
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_dest, myrank, wm_tmp(1,1,icyb0,icx), ncc, 
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_src, ipy_src,mpi_comm_world,istatus,ierr )
!      !ndis(me)[ipy_src+1] = icyb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icyb0,status)


!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipy_dest+1)
!         wm_tmp( :, :,         nb:nb        +nd, icx )[ipy_dest+1]
!     . = wm_tmp( :, :, icyb0prior:icyb0prior+nd, icx ) ! Put

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icyb0prior,kind=8),int(icyb0prior+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
            endif
         END DO
      END DO

! ULmoment -Y
      ipy_dest = ipz*npy*npx
     &     + mod(ipy - np_supercelly + 1/npy + npy, npy)*npx + ipx
      ipy_src = ipz*npy*npx
     &     + mod(ipy + np_supercelly - 1/npy + npy, npy)*npx + ipx
      nitr = (nscell - 1) / 2
      if(nscell > npy) nitr = (npy -1) / 2

      DO icx = icx0, icx1
         DO itr = 1, nitr
            if (itr ==1) then
               icyb0 = mod(icy0 + nscydiv - 1, nscell) + 1
               icyb1 = mod(icy1 + nscydiv - 1, nscell) + 1
               ncc = m1 * nscell * (icyb1 - icyb0 + 1)
!coarray               call mpi_sendrecv(wm_tmp(1,1,icy0,icx), ncc,
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_dest, myrank, wm_tmp(1,1,icyb0,icx), ncc,
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_src, ipy_src,mpi_comm_world,istatus,ierr)
      nd = abs(icyb1 - icyb0)
!      ndis(me)[ipy_src+1] = icyb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icyb0,status)

!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipy_dest+1)
!         wm_tmp( :, :,   nb:nb  +nd, icx )[ipy_dest+1]
!     . = wm_tmp( :, :, icy0:icy0+nd, icx ) ! Put

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icy0,kind=8),int(icy0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
            else
               icyb0prior = icyb0
               icyb0 = mod(icy0 + nscydiv * itr - 1, nscell) + 1
!coarray               call mpi_sendrecv(wm_tmp(1,1,icyb0prior,icx), ncc, 
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_dest, myrank, wm_tmp(1,1,icyb0,icx), ncc, 
!coarray     &              MPI_DOUBLE_COMPLEX,
!coarray     &              ipy_src, ipy_src,mpi_comm_world,istatus,ierr)
!      ndis(me)[ipy_src+1] = icyb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icyb0,status)

!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipy_dest+1)
!         wm_tmp( :, :,         nb:nb        +nd, icx )[ipy_dest+1]
!     . = wm_tmp( :, :, icyb0prior:icyb0prior+nd, icx ) ! Put

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icyb0prior,kind=8),int(icyb0prior+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
            end if
         END DO
      END DO

! ULmoment +X
      ipx_dest = ipz*npy*npx + ipy*npx
     &     + mod(ipx + np_supercellx - 1/npx + npx, npx)
      ipx_src = ipz*npy*npx + ipy*npx
     &     + mod(ipx - np_supercellx + 1/npx + npx, npx)
      nitr = nscell / 2
      if(nscell > npx) nitr = npx / 2

      DO itr = 1, nitr
         if (itr ==1) then
            icxb0 = mod(icx0 - nscxdiv - 1 + nscell, nscell) + 1
            icxb1 = mod(icx1 - nscxdiv - 1 + nscell, nscell) + 1
            ncc = m1 * nscell * nscell * (icxb1 - icxb0 + 1)
!coarray            call mpi_sendrecv(wm_tmp(1,1,1,icx0), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_dest, myrank, wm_tmp(1,1,1,icxb0), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_src, ipx_src,mpi_comm_world,istatus,ierr )
      nd = abs(icxb1 - icxb0)
!      ndis(me)[ipx_src+1] = icxb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxb0,status)

!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipx_dest+1)
!         wm_tmp( :, :, :,   nb:nb  +nd )[ipx_dest+1]
!     . = wm_tmp( :, :, :, icx0:icx0+nd ) ! Put

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx0,kind=8),int(icx0+nd,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
         else
            icxb0prior = icxb0
            icxb0 = mod(icx0 - nscxdiv * itr - 1 + nscell, nscell) + 1
!coarray            call mpi_sendrecv(wm_tmp(1,1,1,icxb0prior), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_dest, myrank, wm_tmp(1,1,1,icxb0), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_src, ipx_src,mpi_comm_world,istatus,ierr )
!      ndis(me)[ipx_src+1] = icxb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxb0,status)


!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipx_dest+1)
!         wm_tmp( :, :, :,         nb:nb        +nd )[ipx_dest+1]
!     . = wm_tmp( :, :, :, icxb0prior:icxb0prior+nd ) ! Put

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxb0prior,kind=8),int(icxb0prior+nd,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
         end if
      END DO

! ULmoment -X
      ipx_dest = ipz*npy*npx + ipy*npx
     &     + mod(ipx - np_supercellx + 1/npx + npx, npx)
      ipx_src = ipz*npy*npx + ipy*npx
     &     + mod(ipx + np_supercellx - 1/npx + npx, npx)
      nitr = (nscell - 1) / 2
      if(nscell > npx) nitr = (npx - 1) / 2

      DO itr = 1, nitr
         if (itr ==1) then
            icxb0 = mod(icx0 + nscxdiv - 1, nscell) + 1
            icxb1 = mod(icx1 + nscxdiv - 1, nscell) + 1
            ncc = m1 * nscell * nscell * (icxb1 - icxb0 + 1)
!coarray            call mpi_sendrecv(wm_tmp(1,1,1,icx0), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_dest, myrank, wm_tmp(1,1,1,icxb0), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_src, ipx_src,mpi_comm_world,istatus,ierr )
      nd = abs(icxb1 - icxb0)
!      ndis(me)[ipx_src+1] = icxb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxb0,status)
!      sync all
      call xmp_sync_all(status)

!      nb = ndis(ipx_dest+1)
!         wm_tmp( :, :, :,   nb:nb  +nd )[ipx_dest+1]
!     . = wm_tmp( :, :, :, icx0:icx0+nd ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx0,kind=8),int(icx0+nd,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
         else
            icxb0prior = icxb0
            icxb0 = mod(icx0 + nscxdiv * itr - 1, nscell) + 1
!coarray            call mpi_sendrecv(wm_tmp(1,1,1,icxb0prior), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_dest, myrank, wm_tmp(1,1,1,icxb0), ncc, 
!coarray     &           MPI_DOUBLE_COMPLEX,
!coarray     &           ipx_src, ipx_src,mpi_comm_world,istatus,ierr )
!      ndis(me)[ipx_src+1] = icxb0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_src+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxb0,status)

!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipx_dest)
!         wm_tmp( :, :, :,         nb:nb        +nd )[ipx_dest+1]
!     . = wm_tmp( :, :, :, icxb0prior:icxb0prior+nd ) ! Put

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(nscell,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxb0prior,kind=8),int(icxb0prior+nd,kind=8),1,status)

      img_dims(1) = ipx_dest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
         end if
      END DO

!coarray
      wm = wm_tmp
!!

      return
      END
c---------------------------------------------------------------------
      subroutine comm_fmm_local_multi(ilevel, mylm, wm,
     $             lclz, lcly, lclx, nbsize, nscydiv, nscxdiv)
!! LL ver.20120211 
c---------------------------------------------------------------------
      use comm_base
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use xmp_api
      implicit none
      include 'mpif.h'
      integer(4) :: ilevel
      integer(4) :: mylm, lclz, lcly, lclx, nbsize
      integer(4) :: nscxdiv, nscydiv, nsczdiv
      complex(8) :: wm(mylm, lclz, lcly, lclx)
      integer(4) :: nbound_zm, nbound_ym, nbound_xm
      integer(4) :: nbound_zp, nbound_yp, nbound_xp

      complex(8)  ccbufp(mylm*nbsize*nscydiv*nscxdiv)
      integer(8) :: ccbufp_local_sec,ccbufp_local_desc
      integer(8), dimension(1) :: ccbufp_lb,ccbufp_ub

!coarray      complex(8) rccbufp(mylm*nbsize*nscydiv*nscxdiv, 2)
      complex(8)  ccbufm(mylm*nbsize*nscydiv*nscxdiv)
      integer(8) :: ccbufm_local_sec,ccbufm_local_desc
      integer(8), dimension(1) :: ccbufm_lb,ccbufm_ub
!coarray      complex(8) rccbufm(mylm*nbsize*nscydiv*nscxdiv, 2)
!      complex(8),allocatable :: rccbufp(:,:)[:]
      complex(8), POINTER :: rccbufp(:,:) => null()
      integer(8) :: rccbufp_desc
      integer(8) :: rccbufp_l_sec, rccbufp_r_sec
      integer(8), dimension(2) :: rccbufp_lb, rccbufp_ub
!      complex(8),allocatable :: rccbufm(:,:)[:]
      complex(8), POINTER :: rccbufm(:,:) => null()
      integer(8) :: rccbufm_desc
      integer(8) :: rccbufm_l_sec, rccbufm_r_sec
      integer(8), dimension(2) :: rccbufm_lb, rccbufm_ub
!      complex(8),allocatable :: wm_tmp(:,:,:,:)[:]
      complex(8), POINTER :: wm_tmp(:,:,:,:) => null()
      integer(8) :: wm_tmp_desc
      integer(8) :: wm_tmp_l_sec, wm_tmp_r_sec
      integer(8), dimension(2) :: wm_tmp_lb, wm_tmp_ub

!      integer,allocatable :: ndis(:)[:]
      integer, POINTER :: ndis(:) => null()
      integer(8) :: ndis_desc
      integer(8) :: ndis_sec
      integer(8), dimension(1) :: ndis_lb,ndis_ub

!      integer,allocatable :: mdis(:)[:]
      integer, POINTER :: mdis(:) => null()
      integer(8) :: mdis_desc
      integer(8) :: mdis_sec
      integer(8), dimension(1) :: mdis_lb,mdis_ub

      integer me, np, nb, nd, mb, md
!!
      integer m1
      integer np_supercell
      integer mcell_size
      integer ipz_pdest, ipy_pdest, ipx_pdest
      integer ipz_psrc, ipy_psrc, ipx_psrc
      integer ipz_mdest, ipy_mdest, ipx_mdest
      integer ipz_msrc, ipy_msrc, ipx_msrc
      integer nitr, itr
      integer nf_pprovider
      integer nf_mprovider
      integer icz, icy, icx
      integer iczp0, iczp1
      integer iczbp0, iczbp1
      integer iczm0, iczm1
      integer iczbm0, iczbm1
      integer iczb
      integer icyp0, icyp1
      integer icybp0
      integer icym0, icym1
      integer icybm0
      integer icxp0, icxp1
      integer icxbp0
      integer icxm0, icxm1
      integer icxbm0
      integer icybp0prior
      integer icxbp0prior
      integer icybm0prior
      integer icxbm0prior
      integer ncc2
      integer nccp
      integer nccm
      integer m
      integer ibs, ibr
      integer istatus(mpi_status_size, 4), ierr
      integer(4) status
#ifndef SYNC_COM
      integer,dimension(4) :: irq
      integer nrq
#endif
!coarray
!      allocate( rccbufp(mylm*nbsize*nscydiv*nscxdiv, 2)[*] )
      rccbufp_lb(1)=1
      rccbufp_lb(2)=1
      rccbufp_ub(1)=mylm*nbsize*nscydiv*nscxdiv
      rccbufp_ub(2)=2
      call xmp_new_coarray(rccbufp_desc,
     & 8,2,rccbufp_lb,rccbufp_ub,1,img_dims)
      call xmp_coarray_bind(rccbufp_desc,rccbufp)
      call xmp_new_array_section(rccbufp_l_sec,2)
      call xmp_new_array_section(rccbufp_r_sec,2)

!      allocate( rccbufm(mylm*nbsize*nscydiv*nscxdiv, 2)[*] )
      rccbufm_lb(1)=1
      rccbufm_lb(2)=1
      rccbufm_ub(1)=mylm*nbsize*nscydiv*nscxdiv
      rccbufm_ub(2)=2
      call xmp_new_coarray(rccbufm_desc,
     & 8,2,rccbufm_lb,rccbufm_ub,1,img_dims)
      call xmp_coarray_bind(rccbufm_desc,rccbufm)
      call xmp_new_array_section(rccbufm_l_sec,2)
      call xmp_new_array_section(rccbufm_r_sec,2)

!      allocate( wm_tmp(mylm, lclz, lcly, lclx)[*] )
      wm_tmp_lb(1)=1
      wm_tmp_lb(2)=1
      wm_tmp_lb(3)=1
      wm_tmp_lb(4)=1
      wm_tmp_ub(1)=mylm
      wm_tmp_ub(2)=lclz
      wm_tmp_ub(3)=lcly
      wm_tmp_ub(4)=lclx
      call xmp_new_coarray(wm_tmp_desc,
     & 8,4,wm_tmp_lb,wm_tmp_ub,1,img_dims)
      call xmp_coarray_bind(wm_tmp_desc,wm_tmp)
      call xmp_new_array_section(wm_tmp_l_sec,4)
      call xmp_new_array_section(wm_tmp_r_sec,4)

      wm_tmp = wm
!      me = this_image()
      me = xmp_this_image()
!      np = num_images()
      np = xmp_num_images()
!      allocate( ndis(np)[*] )
      ndis_lb(1)=1
      ndis_ub(1)=np
      call xmp_new_coarray(ndis_desc,4,1,ndis_lb,ndis_ub,1,img_dims)
      call xmp_coarray_bind(ndis_desc,ndis)
      call xmp_new_array_section(ndis_sec,1)
!      allocate( mdis(np)[*] )
      mdis_lb(1)=1
      mdis_ub(1)=np
      call xmp_new_coarray(mdis_desc,4,1,mdis_lb,mdis_ub,1,img_dims)
      call xmp_coarray_bind(mdis_desc,mdis)
      call xmp_new_array_section(mdis_sec,1)

      ccbufp_lb(1) = 1
      ccbufp_ub(1) = mylm*nbsize*nscydiv*nscxdiv
      call xmp_new_local_array(ccbufp_local_desc,8,1,
     & ccbufp_lb,ccbufp_ub,loc(ccbufp))
      call xmp_new_array_section(ccbufp_local_sec,1)

      ccbufm_lb(1) = 1
      ccbufm_ub(1) = mylm*nbsize*nscydiv*nscxdiv
      call xmp_new_local_array(ccbufm_local_desc,8,1,
     & ccbufm_lb,ccbufm_ub,loc(ccbufm))
      call xmp_new_array_section(ccbufm_local_sec,1)

!!

! ---- 3D rank order rule. ----
!     ipx=mod(myrank, npx)
!     ipy=mod((myrank - ipx) / npx, npy)
!     ipz=mod((myrank - ipx - ipy*npx) / (npx*npy), npz)

!=== local constant ===!
      m1 = mylm
      mcell_size = fmm_data(ilevel)%mcell_size

!=== global constant ===!
      nsczdiv = fmm_data(ilevel)%nsczdiv
      nscydiv = fmm_data(ilevel)%nscydiv
      nscxdiv = fmm_data(ilevel)%nscxdiv

      nbound_zm=fmm_data(ilevel)%nbound_zm
      nbound_zp=fmm_data(ilevel)%nbound_zp
      nbound_ym=fmm_data(ilevel)%nbound_ym
      nbound_yp=fmm_data(ilevel)%nbound_yp
      nbound_xm=fmm_data(ilevel)%nbound_xm
      nbound_xp=fmm_data(ilevel)%nbound_xp

!----- lower level moment communication starts here. -----

! LLmoment +Z
      np_supercell = (npz * mcell_size - 1) / ncell + 1
      ipz_pdest = mod(ipz + np_supercell - 1/npz + npz, npz)*npy*npx
     &  + ipy*npx + ipx
      ipz_psrc  = mod(ipz - np_supercell + 1/npz + npz, npz)*npy*npx
     &  + ipy*npx + ipx
! LLmoment -Z
      ipz_mdest = mod(ipz - np_supercell + 1/npz  + npz, npz)*npy*npx
     &  + ipy*npx + ipx
      ipz_msrc  = mod(ipz + np_supercell - 1/npz  + npz, npz)*npy*npx
     &  + ipy*npx + ipx

      nf_pprovider = 0
      nf_mprovider = 0
      if(nbound_zm == 4 .and. nsczdiv == 1) nf_pprovider = 1
      if(nbound_zp == 4 .and. nsczdiv == 1) nf_mprovider = 1
      nitr = max( (nbound_zm - 1) / nsczdiv + 1,
     &            (nbound_zp - 1) / nsczdiv + 1    )


!coarray
!      allocate( rccbufp(mylm*nbsize*nscydiv*nscxdiv, 2)[*] )
!      allocate( rccbufm(mylm*nbsize*nscydiv*nscxdiv, 2)[*] )
!!
      do itr = 1, nitr
         if (itr == 1) then        ! first iteration
            iczp0 = nsczdiv + 1
            if(nsczdiv < nbound_zm) iczp0 = nbound_zm + 1
            iczp1 = nbound_zm + nsczdiv
            iczbp0 = iczp0 - nsczdiv
            iczbp1 = iczp1 - nsczdiv

            nccp = 0
            do icx = nbound_xm + 1, nbound_xm + nscxdiv
               do icy = nbound_ym + 1, nbound_ym + nscydiv
                  do icz = iczp0, iczp1
                     do m = 1, mylm
                        nccp = nccp + 1
                        ccbufp(nccp) = wm_tmp(m, icz, icy, icx )
                     end do
                  end do
               end do
            end do

            iczm0 = nbound_zm + 1
            iczm1 = nbound_zm + nbound_zp
            if(nsczdiv < nbound_zp) iczm1 = nbound_zm + nsczdiv
            iczbm0 = iczm0 + nsczdiv
            iczbm1 = iczm1 + nsczdiv
            nccm = 0
            do icx = nbound_xm + 1, nbound_xm + nscxdiv
               do icy = nbound_ym + 1, nbound_ym + nscydiv
                  do icz = iczm0, iczm1
                     do m = 1, mylm
                        nccm = nccm + 1
                        ccbufm(nccm) = wm_tmp(m, icz, icy, icx )
                     end do
                  end do
               end do
            end do

#ifdef SYNC_COM
!coarray            call mpi_sendrecv(ccbufp, nccp, MPI_DOUBLE_COMPLEX,
!coarray     &              ipz_pdest, myrank,
!coarray     &              rccbufp(1,1), nccp, MPI_DOUBLE_COMPLEX,
!coarray     &              ipz_psrc, ipz_psrc, 
!coarray     &              mpi_comm_world, istatus, ierr)
!coarray            call mpi_sendrecv(ccbufm, nccm, MPI_DOUBLE_COMPLEX,
!coarray     &              ipz_mdest, myrank,
!coarray     &              rccbufm(1,1), nccm, MPI_DOUBLE_COMPLEX, 
!coarray     &              ipz_msrc, ipz_msrc,
!coarray     &              mpi_comm_world, istatus, ierr )
!      rccbufp(1:nccp,1)[ipz_pdest+1] = ccbufp(1:nccp) ! Put
      call xmp_array_section_set_triplet(rccbufp_r_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufp_r_sec,
     & 2,int(1,kind=8),int(1,kind=8),1,status)

      call xmp_array_section_set_triplet(ccbufp_local_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put_local(img_dims,rccbufp_desc,rccbufp_r_sec, 
     & ccbufp_local_desc,ccbufp_local_sec,status)


!      rccbufm(1:nccm,1)[ipz_mdest+1] = ccbufm(1:nccm) ! Put
      call xmp_array_section_set_triplet(rccbufm_r_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufm_r_sec,
     & 2,int(1,kind=8),int(1,kind=8),1,status)

      call xmp_array_section_set_triplet(ccbufm_local_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put_local(img_dims,rccbufm_desc,rccbufm_r_sec, 
     & ccbufm_local_desc,ccbufm_local_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
            call mpi_irecv(rccbufp(1,1), nccp,
     &              MPI_DOUBLE_COMPLEX, ipz_psrc, ipz_psrc, 
     &              mpi_comm_world, irq(1), ierr)
            call mpi_isend(ccbufp, nccp,
     &              MPI_DOUBLE_COMPLEX, ipz_pdest, myrank, 
     &              mpi_comm_world, irq(2), ierr)
            call mpi_irecv(rccbufm(1,1), nccm,
     &              MPI_DOUBLE_COMPLEX, ipz_msrc, ipz_msrc, 
     &              mpi_comm_world, irq(3), ierr)
            call mpi_isend(ccbufm, nccm,
     &              MPI_DOUBLE_COMPLEX, ipz_mdest, myrank, 
     &              mpi_comm_world, irq(4), ierr)
            nrq = 4
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            ncc2 = 0
            do icx = nbound_xm + 1, nbound_xm + nscxdiv
               do icy = nbound_ym + 1, nbound_ym + nscydiv
                  do iczb = iczbp0, iczbp1
                     do m = 1, mylm
                        ncc2 = ncc2 + 1
                        wm_tmp(m,iczb,icy,icx) = rccbufp(ncc2,1)
                     end do
                  end do
               end do
            end do

            ncc2 = 0
            do icx = nbound_xm + 1, nbound_xm + nscxdiv
               do icy = nbound_ym + 1, nbound_ym + nscydiv
                  do iczb = iczbm0, iczbm1
                     do m = 1, mylm
                        ncc2 = ncc2 + 1
                        wm_tmp(m,iczb,icy,icx) = rccbufm(ncc2,1)
                     end do
                  end do
               end do
            end do

         else          ! iteration follows

            ibs = mod(itr, 2) + 1
            ibr = mod(itr + 1, 2) + 1
            iczbp0 = nbound_zm - nsczdiv*itr  + 1
            iczbp1 = iczbp0 + nsczdiv - 1
            iczbm0 = nbound_zm + nsczdiv + nsczdiv * (itr - 1) + 1
            iczbm1 = iczbm0 + nsczdiv - 1

            if(nsczdiv /= 1 .or. (nsczdiv == 1 .and. itr < nitr)) then

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(rccbufp(1,ibs), nccp,
!coarray     &                MPI_DOUBLE_COMPLEX,
!coarray     &                ipz_pdest, myrank, 
!coarray     &                rccbufp(1,ibr), nccp, MPI_DOUBLE_COMPLEX, 
!coarray     &                ipz_psrc, ipz_psrc, 
!coarray     &                mpi_comm_world, istatus, ierr )
!coarray               call mpi_sendrecv(rccbufm(1,ibs), nccm,
!coarray     &                MPI_DOUBLE_COMPLEX,
!coarray     &                ipz_mdest, myrank, 
!coarray     &                rccbufm(1,ibr), nccm, MPI_DOUBLE_COMPLEX, 
!coarray     &                ipz_msrc, ipz_msrc, 
!coarray     &                mpi_comm_world, istatus, ierr )
!      rccbufp(1:nccp,ibr)[ipz_pdest+1] = rccbufp(1:nccp,ibs) ! Put
      call xmp_array_section_set_triplet(rccbufp_r_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufp_r_sec,
     & 2,int(ibr,kind=8),int(ibr,kind=8),1,status)

      call xmp_array_section_set_triplet(rccbufp_l_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufp_l_sec,
     & 2,int(ibs,kind=8),int(ibs,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,rccbufp_desc,rccbufp_r_sec, 
     & rccbufp_desc,rccbufp_l_sec,status)

!      rccbufm(1:nccm,ibr)[ipz_mdest+1] = rccbufm(1:nccm,ibs) ! Put
      call xmp_array_section_set_triplet(rccbufm_r_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufm_r_sec,
     & 2,int(ibr,kind=8),int(ibr,kind=8),1,status)

      call xmp_array_section_set_triplet(rccbufm_l_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufm_l_sec,
     & 2,int(ibs,kind=8),int(ibs,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,rccbufm_desc,rccbufm_r_sec,
     & rccbufm_desc,rccbufm_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
            call mpi_irecv(rccbufp(1,ibr), nccp,
     &              MPI_DOUBLE_COMPLEX, ipz_psrc, ipz_psrc, 
     &              mpi_comm_world, irq(1), ierr)
            call mpi_isend(rccbufp(1,ibs), nccp,
     &              MPI_DOUBLE_COMPLEX, ipz_pdest, myrank, 
     &              mpi_comm_world, irq(2), ierr)
            call mpi_irecv(rccbufm(1,ibr), nccm,
     &              MPI_DOUBLE_COMPLEX, ipz_msrc, ipz_msrc, 
     &              mpi_comm_world, irq(3), ierr)
            call mpi_isend(rccbufm(1,ibs), nccm,
     &              MPI_DOUBLE_COMPLEX, ipz_mdest, myrank, 
     &              mpi_comm_world, irq(4), ierr)
            nrq = 4
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

               ncc2 = 0
               do icx = nbound_xm + 1, nbound_xm + nscxdiv
                  do icy = nbound_ym + 1, nbound_ym + nscydiv
                     do iczb = iczbp0, iczbp1
                        do m = 1, mylm
                           ncc2 = ncc2 + 1
                           wm_tmp(m,iczb,icy,icx) = rccbufp(ncc2,ibr)
                        end do
                     end do
                  end do
               end do
               ncc2 = 0
               do icx = nbound_xm + 1, nbound_xm + nscxdiv
                  do icy = nbound_ym + 1, nbound_ym + nscydiv
                     do iczb = iczbm0, iczbm1
                        do m = 1, mylm
                           ncc2 = ncc2 + 1
                           wm_tmp(m,iczb,icy,icx) = rccbufm(ncc2,ibr)
                        end do
                     end do
                  end do
               end do

            else      ! = final pairing. =

               if(nf_pprovider == 1) then
#ifdef SYNC_COM
!coarray                  call mpi_send(rccbufp(1,ibs), nccp, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipz_pdest, myrank, 
!coarray     &                   mpi_comm_world, istatus, ierr )
!      rccbufp(1:nccp,ibr)[ipz_pdest+1] = rccbufp(1:nccp,ibs) ! Put
      call xmp_array_section_set_triplet(rccbufp_r_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufp_r_sec,
     & 2,int(ibr,kind=8),int(ibr,kind=8),1,status)

      call xmp_array_section_set_triplet(rccbufp_l_sec,
     & 1,int(1,kind=8),int(nccp,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufp_l_sec,
     & 2,int(ibs,kind=8),int(ibs,kind=8),1,status)
      img_dims(1) = ipz_pdest+1
      call xmp_coarray_put(img_dims,rccbufp_desc,rccbufp_r_sec, 
     & rccbufp_desc,rccbufp_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
                  call mpi_isend(rccbufp(1,ibs), nccp,
     &                   MPI_DOUBLE_COMPLEX, ipz_pdest, myrank, 
     &                   mpi_comm_world, irq(1), ierr)
#endif

               else

#ifdef SYNC_COM
!coarray                  call mpi_recv(rccbufp(1,ibr), nccp, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipz_psrc, ipz_psrc, 
!coarray     &                   mpi_comm_world, istatus, ierr )
#else
                  call mpi_irecv(rccbufp(1,ibr), nccp,
     &                   MPI_DOUBLE_COMPLEX, ipz_psrc, ipz_psrc, 
     &                   mpi_comm_world, irq(1), ierr)
#endif

               endif     ! final provider (p)

               if(nf_mprovider == 1) then
#ifdef SYNC_COM
!coarray                  call mpi_send(rccbufm(1,ibs), nccm, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipz_mdest, myrank, 
!coarray     &                   mpi_comm_world, istatus, ierr )
!      rccbufm(1:nccm,ibr)[ipz_mdest+1] = rccbufm(1:nccm,ibs) ! Put
      call xmp_array_section_set_triplet(rccbufm_r_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufm_r_sec,
     & 2,int(ibr,kind=8),int(ibr,kind=8),1,status)

      call xmp_array_section_set_triplet(rccbufm_l_sec,
     & 1,int(1,kind=8),int(nccm,kind=8),1,status)
      call xmp_array_section_set_triplet(rccbufm_l_sec,
     & 2,int(ibs,kind=8),int(ibs,kind=8),1,status)
      img_dims(1) = ipz_mdest+1
      call xmp_coarray_put(img_dims,rccbufm_desc,rccbufm_r_sec, 
     & rccbufm_desc,rccbufm_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
#else
                  call mpi_isend(rccbufm(1,ibs), nccm,
     &                   MPI_DOUBLE_COMPLEX, ipz_mdest, myrank, 
     &                   mpi_comm_world, irq(2), ierr)
#endif

               else
#ifdef SYNC_COM
!coarray                  call mpi_recv(rccbufm(1,ibr), nccm, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipz_msrc, ipz_msrc, 
!coarray     &                   mpi_comm_world, istatus, ierr )
#else
                  call mpi_irecv(rccbufm(1,ibr), nccm,
     &                   MPI_DOUBLE_COMPLEX, ipz_msrc, ipz_msrc, 
     &                   mpi_comm_world, irq(2), ierr)
#endif

               endif     ! final provider (m)
#ifndef SYNC_COM
               nrq = 2
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

               if(nf_pprovider == 0) then
                  ncc2 = 0
                  do icx = nbound_xm + 1, nbound_xm + nscxdiv
                     do icy = nbound_ym + 1, nbound_ym + nscydiv
                        do iczb = iczbp0, iczbp1
                           do m = 1, mylm
                              ncc2 = ncc2 + 1
                              wm_tmp(m,iczb,icy,icx) = rccbufp(ncc2,ibr)
                           end do
                        end do
                     end do
                  end do
               endif     ! final receiver (p)

               if(nf_mprovider == 0) then
                  ncc2 = 0
                  do icx = nbound_xm + 1, nbound_xm + nscxdiv
                     do icy = nbound_ym + 1, nbound_ym + nscydiv
                        do iczb = iczbm0, iczbm1
                           do m = 1, mylm
                              ncc2 = ncc2 + 1
                              wm_tmp(m,iczb,icy,icx) = rccbufm(ncc2,ibr)
                           end do
                        end do
                     end do
                  end do
               endif     ! final receiver (m)

            endif     ! final pairing. !ordinary root

         endif   ! iteration
      end do     ! iteration

! LLmoment +Y
      np_supercell = (npy * mcell_size - 1) / ncell + 1
      ipy_pdest = ipz*npy*npx 
     &  + mod(ipy + np_supercell - 1/npy + npy, npy)*npx + ipx
      ipy_psrc = ipz*npy*npx 
     &  + mod(ipy - np_supercell + 1/npy + npy, npy)*npx + ipx
! LLmoment -Y
      ipy_mdest = ipz*npy*npx
     &  + mod(ipy - np_supercell + 1/npy + npy, npy)*npx + ipx
      ipy_msrc = ipz*npy*npx 
     &  + mod(ipy + np_supercell - 1/npy + npy, npy)*npx + ipx

      nf_pprovider = 0
      nf_mprovider = 0
      if(nbound_ym == 4 .and. nscydiv == 1) nf_pprovider = 1
      if(nbound_yp == 4 .and. nscydiv == 1) nf_mprovider = 1
      nitr = max( (nbound_ym - 1) / nscydiv + 1,
     &            (nbound_yp - 1) / nscydiv + 1 )

      do icx = nbound_xm + 1, nbound_xm + nscxdiv
         do itr = 1, nitr
            if (itr == 1) then                ! first iteration

               icyp0 = nscydiv + 1
               if(nscydiv < nbound_ym) icyp0 = nbound_ym + 1
               icyp1 = nbound_ym + nscydiv
               icybp0 = icyp0 - nscydiv
               nccp = (nbound_zm + nsczdiv + nbound_zp)
     &               *(icyp1 - icyp0 + 1) * mylm
               icym0 = nbound_ym + 1
               icym1 = nbound_ym + nbound_yp
               if(nscydiv < nbound_yp) icym1 = nbound_ym + nscydiv
               icybm0 = icym0 + nscydiv
               nccm = (nbound_zm + nsczdiv + nbound_zp)
     &           *(icym1 - icym0 + 1) * mylm

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(wm_tmp(1,1,icyp0,icx), nccp,
!coarray     &                MPI_DOUBLE_COMPLEX, 
!coarray     &                ipy_pdest, myrank, 
!coarray     &                wm_tmp(1,1,icybp0,icx), nccp, MPI_DOUBLE_COMPLEX,
!coarray     &                ipy_psrc, ipy_psrc, 
!coarray     &                mpi_comm_world, istatus, ierr  )
      nd = abs(icyp1 - icyp0)
!      ndis(me)[ipy_psrc+1] = icybp0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_psrc+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icybp0,status)


!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipy_pdest+1)
!         wm_tmp( :, :,    nb:nb   +nd, icx )[ipy_pdest+1]
!     . = wm_tmp( :, :, icyp0:icyp0+nd, icx ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icyp0,kind=8),int(icyp0+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)


      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)



!      sync all
      call xmp_sync_all(status)
!!
!coarray               call mpi_sendrecv(wm_tmp(1,1,icym0,icx), nccm, 
!coarray     &                MPI_DOUBLE_COMPLEX, 
!coarray     &                ipy_mdest, myrank, 
!coarray     &                wm_tmp(1,1,icybm0,icx), nccm, MPI_DOUBLE_COMPLEX,
!coarray     &                ipy_msrc, ipy_msrc, 
!coarray     &                mpi_comm_world, istatus, ierr )
      md = abs(icym1 - icym0)
!      mdis(me)[ipy_msrc+1] = icybm0 ! Put
      call xmp_array_section_set_triplet(mdis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_msrc+1
      call xmp_coarray_put_scalar(img_dims,mdis_desc,
     & mdis_sec,icybm0,status)



!      sync all
      call xmp_sync_all(status)
      mb = mdis(ipy_mdest+1)
!         wm_tmp( :, :, mb:mb      +md, icx )[ipy_mdest+1]
!     . = wm_tmp( :, :, icym0:icym0+md, icx ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(mb,kind=8),int(mb+md,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icym0,kind=8),int(icym0+md,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
               call mpi_irecv(wm_tmp(1,1,icybp0,icx), nccp,
     &                MPI_DOUBLE_COMPLEX, ipy_psrc, ipy_psrc, 
     &                mpi_comm_world, irq(1), ierr)
               call mpi_isend(wm_tmp(1,1,icyp0,icx), nccp,
     &                MPI_DOUBLE_COMPLEX, ipy_pdest, myrank, 
     &                mpi_comm_world, irq(2), ierr)
               call mpi_irecv(wm_tmp(1,1,icybm0,icx), nccm,
     &                MPI_DOUBLE_COMPLEX, ipy_msrc, ipy_msrc, 
     &                mpi_comm_world, irq(3), ierr)
               call mpi_isend(wm_tmp(1,1,icym0,icx), nccm,
     &                MPI_DOUBLE_COMPLEX, ipy_mdest, myrank, 
     &                mpi_comm_world, irq(4), ierr)
               nrq = 4
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            else                  ! iteration follows

               icybp0prior = icybp0
               icybp0 = nbound_ym - nscydiv*itr + 1
               icybm0prior = icybm0
               icybm0 = nbound_ym + nscydiv + nscydiv * (itr - 1) + 1

               if(nscydiv /= 1 .or.
     &                        (nscydiv == 1 .and. itr < nitr)) then

#ifdef SYNC_COM
!coarray                  call mpi_sendrecv(wm_tmp(1,1,icybp0prior,icx), nccp,
!coarray     &                 MPI_DOUBLE_COMPLEX, 
!coarray     &                 ipy_pdest, myrank, 
!coarray     &                 wm_tmp(1,1,icybp0,icx), nccp, MPI_DOUBLE_COMPLEX,
!coarray     *                 ipy_psrc, ipy_psrc, 
!coarray     &                 mpi_comm_world, istatus, ierr )
!      ndis(me)[ipy_psrc+1] = icybp0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_psrc+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icybp0,status)
!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipy_pdest+1)
!         wm_tmp( :, :,          nb:nb         +nd, icx )[ipy_pdest+1] ! Put
!     . = wm_tmp( :, :, icybp0prior:icybp0prior+nd, icx )

      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icybp0prior,kind=8),int(icybp0prior+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
!coarray                  call mpi_sendrecv(wm_tmp(1,1,icybm0prior,icx), nccm,
!coarray     &                 MPI_DOUBLE_COMPLEX, 
!coarray     &                 ipy_mdest, myrank, 
!coarray     &                 wm_tmp(1,1,icybm0,icx), nccm, MPI_DOUBLE_COMPLEX,
!coarray     &                 ipy_msrc, ipy_msrc, 
!coarray     &                 mpi_comm_world, istatus, ierr )
!      mdis(me)[ipy_msrc+1] = icybm0 ! Put
      call xmp_array_section_set_triplet(mdis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_msrc+1
      call xmp_coarray_put_scalar(img_dims,mdis_desc,
     & mdis_sec,icybm0,status)

!      sync all
      call xmp_sync_all(status)
      mb = mdis(ipy_mdest+1)
!         wm_tmp( :, :,          mb:mb         +md, icx )[ipy_mdest+1] ! Put
!     . = wm_tmp( :, :, icybm0prior:icybm0prior+md, icx )
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(mb,kind=8),int(mb+md,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icybm0prior,kind=8),int(icybm0prior+md,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_mdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
#else
                  call mpi_irecv(wm_tmp(1,1,icybp0,icx), nccp,
     &                   MPI_DOUBLE_COMPLEX, ipy_psrc, ipy_psrc, 
     &                   mpi_comm_world, irq(1), ierr)
                  call mpi_isend(wm_tmp(1,1,icybp0prior,icx), nccp,
     &                   MPI_DOUBLE_COMPLEX, ipy_pdest, myrank, 
     &                   mpi_comm_world, irq(2), ierr)
                  call mpi_irecv(wm_tmp(1,1,icybm0,icx), nccm,
     &                   MPI_DOUBLE_COMPLEX, ipy_msrc, ipy_msrc, 
     &                   mpi_comm_world, irq(3), ierr)
                  call mpi_isend(wm_tmp(1,1,icybm0prior,icx), nccm,
     &                   MPI_DOUBLE_COMPLEX, ipy_mdest, myrank, 
     &                   mpi_comm_world, irq(4), ierr)
                  nrq = 4
                  call mpi_waitall(nrq, irq, istatus, ierr)
#endif

               else     !  = final pairing. =

                  if(nf_pprovider == 1) then
#ifdef SYNC_COM
!coarray                     call mpi_send(wm_tmp(1,1,icybp0prior,icx), nccp, 
!coarray     &                      MPI_DOUBLE_COMPLEX, 
!coarray     &                      ipy_pdest, myrank, 
!coarray     &                      mpi_comm_world, istatus, ierr )
!      ndis(me)[ipy_psrc+1] = icybp0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_psrc+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icybp0,status)
!      sync all
      call xmp_sync_all(status)
!      nb = ndis(ipy_pdest+1)
!         wm_tmp( :, :,          nb:nb         +nd, icx )[ipy_pdest+1]
!     . = wm_tmp( :, :, icybp0prior:icybp0prior+nd, icx ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(nb,kind=8),int(nb+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icybp0prior,kind=8),int(icybp0prior+nd,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      img_dims(1) = ipy_pdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
#else
                     call mpi_isend(wm_tmp(1,1,icybp0prior,icx), nccp,
     &                      MPI_DOUBLE_COMPLEX, ipy_pdest, myrank, 
     &                      mpi_comm_world, irq(1), ierr)
#endif

                  else
#ifdef SYNC_COM
!coarray                     call mpi_recv(wm_tmp(1,1,icybp0,icx), nccp, 
!coarray     &                      MPI_DOUBLE_COMPLEX, 
!coarray     &                      ipy_psrc, ipy_psrc, 
!coarray     &                      mpi_comm_world, istatus, ierr )
#else
                     call mpi_irecv(wm_tmp(1,1,icybp0,icx), nccp,
     &                      MPI_DOUBLE_COMPLEX, ipy_psrc, ipy_psrc, 
     &                      mpi_comm_world, irq(1), ierr)
#endif

                  endif         ! final provider (p)

                  if(nf_mprovider == 1) then
#ifdef SYNC_COM
!coarray                     call mpi_send(wm_tmp(1,1,icybm0prior,icx), nccm, 
!coarray     &                      MPI_DOUBLE_COMPLEX,
!coarray     &                      ipy_mdest, myrank, 
!coarray     &                      mpi_comm_world, istatus, ierr )
!      mdis(me)[ipy_msrc+1] = icybm0 ! Put
      call xmp_array_section_set_triplet(mdis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipy_msrc+1
      call xmp_coarray_put_scalar(img_dims,mdis_desc,
     & mdis_sec,icybm0,status)
!      sync all
      call xmp_sync_all(status)
      md = mdis(ipy_mdest+1)
!         wm_tmp( :, :,          mb:mb         +md, icx )[ipy_mdest+1]
!     . = wm_tmp( :, :, icybm0prior:icybm0prior+md, icx ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(mb,kind=8),int(mb+md,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(icybm0prior,kind=8),int(icybm0prior+md,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icx,kind=8),int(icx,kind=8),1,status)
!      sync all
      call xmp_sync_all(status)
#else
                     call mpi_isend(wm_tmp(1,1,icybm0prior,icx), nccm,
     &                      MPI_DOUBLE_COMPLEX, ipy_mdest, myrank, 
     &                      mpi_comm_world, irq(2), ierr)
#endif

                  else
#ifdef SYNC_COM
!coarray                     call mpi_recv(wm_tmp(1,1,icybm0,icx), nccm, 
!coarray     &                      MPI_DOUBLE_COMPLEX,
!coarray     &                      ipy_msrc, ipy_msrc, 
!coarray     &                      mpi_comm_world, istatus, ierr )
#else
                     call mpi_irecv(wm_tmp(1,1,icybm0,icx), nccm,
     &                      MPI_DOUBLE_COMPLEX, ipy_msrc, ipy_msrc, 
     &                      mpi_comm_world, irq(2), ierr)
#endif

                  endif         ! final provider (m)
#ifndef SYNC_COM
                  nrq = 2
                  call mpi_waitall(nrq, irq, istatus, ierr)
#endif

               endif         ! final pairing. !ordinary root

            endif         ! iteration
         end do           ! iteration
      end do            ! ipx

! LLmoment +X
      np_supercell = (npx* mcell_size - 1) / ncell + 1
      ipx_pdest = ipz*npy*npx + ipy*npx
     &  + mod(ipx + np_supercell - 1/npx + npx, npx)
      ipx_psrc = ipz*npy*npx + ipy*npx
     &  + mod(ipx - np_supercell + 1/npx + npx, npx)
! LLmoment -X
      ipx_mdest = ipz*npy*npx + ipy*npx
     &  + mod(ipx - np_supercell + 1/npx + npx, npx)
      ipx_msrc = ipz*npy*npx + ipy*npx
     &  + mod(ipx + np_supercell - 1/npx + npx, npx)

      nf_pprovider = 0
      nf_mprovider = 0
      if(nbound_xm == 4 .and. nscxdiv == 1) nf_pprovider = 1
      if(nbound_xp == 4 .and. nscxdiv == 1) nf_mprovider = 1

      nitr = max( (nbound_xm - 1) / nscxdiv + 1,
     &            (nbound_xp - 1) / nscxdiv + 1  )

      do itr = 1, nitr
         if (itr == 1) then    ! first iteration

            icxp0 = nscxdiv + 1
            if(nscxdiv < nbound_xm) icxp0 = nbound_xm + 1
            icxp1 = nbound_xm + nscxdiv
            icxbp0 = icxp0 - nscxdiv
            nccp = (nbound_zm + nsczdiv + nbound_zp)
     &       *(nbound_ym + nscydiv + nbound_yp)*(icxp1 - icxp0 + 1)
     &       *mylm
            icxm0 = nbound_xm + 1
            icxm1 = nbound_xm + nbound_xp
            if(nscxdiv < nbound_xp) icxm1 = nbound_xm + nscxdiv
            icxbm0 = icxm0 + nscxdiv
            nccm = (nbound_zm + nsczdiv + nbound_zp)
     &       *(nbound_ym + nscydiv + nbound_yp)*(icxm1 - icxm0 + 1)
     &       *mylm

#ifdef SYNC_COM
!coarray            call mpi_sendrecv(wm_tmp(1,1,1,icxp0), nccp,
!coarray     &             MPI_DOUBLE_COMPLEX,
!coarray     &             ipx_pdest, myrank, 
!coarray     &             wm_tmp(1,1,1,icxbp0), nccp, MPI_DOUBLE_COMPLEX,
!coarray     &             ipx_psrc, ipx_psrc, 
!coarray     &             mpi_comm_world, istatus, ierr )
!      ndis(me)[ipx_psrc+1] = icxbp0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_psrc+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxbp0,status)
!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipx_pdest+1)
!         wm_tmp( :, :, :,    nb:nb   +nd )[ipx_pdest+1]
!     . = wm_tmp( :, :, :, icxp0:icxp0+nd ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxp0,kind=8),int(icxp0+nd,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)


!      sync all
      call xmp_sync_all(status)
!!
!coarray            call mpi_sendrecv(wm_tmp(1,1,1,icxm0), nccm,
!coarray     &             MPI_DOUBLE_COMPLEX,
!coarray     &             ipx_mdest, myrank, 
!coarray     &             wm_tmp(1,1,1,icxbm0), nccm, MPI_DOUBLE_COMPLEX, 
!coarray     &             ipx_msrc, ipx_msrc, 
!coarray     &             mpi_comm_world, istatus, ierr )
!      mdis(me)[ipx_msrc+1] = icxbm0 ! Put
      call xmp_array_section_set_triplet(mdis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_msrc+1
      call xmp_coarray_put_scalar(img_dims,mdis_desc,
     & mdis_sec,icxbm0,status)
!      sync all
      call xmp_sync_all(status)
      mb = mdis(ipx_mdest+1)
!         wm_tmp( :, :, :,    mb:mb   +md )[ipx_mdest+1]
!     . = wm_tmp( :, :, :, icxm0:icxm0+md )
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(mb,kind=8),int(mb+md,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxm0,kind=8),int(icxm0+md,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
#else
            call mpi_irecv(wm_tmp(1,1,1,icxbp0), nccp,
     &             MPI_DOUBLE_COMPLEX, ipx_psrc, ipx_psrc, 
     &             mpi_comm_world, irq(1), ierr)
            call mpi_isend(wm_tmp(1,1,1,icxp0), nccp,
     &             MPI_DOUBLE_COMPLEX, ipx_pdest, myrank, 
     &             mpi_comm_world, irq(2), ierr)
            call mpi_irecv(wm_tmp(1,1,1,icxbm0), nccm,
     &             MPI_DOUBLE_COMPLEX, ipx_msrc, ipx_msrc, 
     &             mpi_comm_world, irq(3), ierr)
            call mpi_isend(wm_tmp(1,1,1,icxm0), nccm,
     &             MPI_DOUBLE_COMPLEX, ipx_mdest, myrank, 
     &             mpi_comm_world, irq(4), ierr)
            nrq = 4
            call mpi_waitall(nrq, irq, istatus, ierr)
#endif

         else                 ! iteration follows

            icxbp0prior = icxbp0
            icxbp0 = nbound_xm - nscxdiv*itr + 1
            icxbm0prior = icxbm0
            icxbm0 = nbound_xm + nscxdiv + nscxdiv * (itr - 1) + 1
!coarray
      nd = abs(icxbp0prior - icxbp0)
      md = abs(icxbm0prior - icxbm0)
!!

            if(nscxdiv /= 1 .or. (nscxdiv == 1 .and. itr < nitr)) then

#ifdef SYNC_COM
!coarray               call mpi_sendrecv(wm_tmp(1,1,1,icxbp0prior), nccp,
!coarray     &                MPI_DOUBLE_COMPLEX, 
!coarray     &                ipx_pdest, myrank, 
!coarray     &                wm_tmp(1,1,1,icxbp0), nccp, MPI_DOUBLE_COMPLEX,
!coarray     &                ipx_psrc, ipx_psrc, 
!coarray     &                mpi_comm_world, istatus, ierr )
!      ndis(me)[ipx_psrc+1] = icxbp0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_psrc+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxbp0,status)
!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipx_pdest+1)
!         wm_tmp( :, :, :,          nb:nb         +nd-1 )[ipx_pdest+1]
!     . = wm_tmp( :, :, :, icxbp0prior:icxbp0prior+nd-1 ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxbp0prior,kind=8),int(icxbp0prior+nd-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)



!      sync all
      call xmp_sync_all(status)
!!
!coarray               call mpi_sendrecv(wm_tmp(1,1,1,icxbm0prior), nccm,
!coarray     &                MPI_DOUBLE_COMPLEX, 
!coarray     &                ipx_mdest, myrank, 
!coarray     &                wm_tmp(1,1,1,icxbm0), nccm, MPI_DOUBLE_COMPLEX, 
!coarray     &                ipx_msrc, ipx_msrc, 
!coarray     &                mpi_comm_world, istatus, ierr )
!      mdis(me)[ipx_msrc+1] = icxbm0 ! Put
      call xmp_array_section_set_triplet(mdis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_msrc+1
      call xmp_coarray_put_scalar(img_dims,mdis_desc,
     & mdis_sec,icxbm0,status)
!      sync all
      call xmp_sync_all(status)
      mb = mdis(ipx_mdest+1)
!         wm_tmp( :, :, :,          mb:mb         +md-1 )[ipx_mdest+1]
!     . = wm_tmp( :, :, :, icxbm0prior:icxbm0prior+md-1 ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(mb,kind=8),int(mb+md-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxbm0prior,kind=8),int(icxbm0prior+md-1,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
#else
               call mpi_irecv(wm_tmp(1,1,1,icxbp0), nccp,
     &                MPI_DOUBLE_COMPLEX, ipx_psrc, ipx_psrc, 
     &                mpi_comm_world, irq(1), ierr)
               call mpi_isend(wm_tmp(1,1,1,icxbp0prior), nccp,
     &                MPI_DOUBLE_COMPLEX, ipx_pdest, myrank, 
     &                mpi_comm_world, irq(2), ierr)
               call mpi_irecv(wm_tmp(1,1,1,icxbm0), nccm,
     &                MPI_DOUBLE_COMPLEX, ipx_msrc, ipx_msrc, 
     &                mpi_comm_world, irq(3), ierr)
               call mpi_isend(wm_tmp(1,1,1,icxbm0prior), nccm,
     &                MPI_DOUBLE_COMPLEX, ipx_mdest, myrank, 
     &                mpi_comm_world, irq(4), ierr)
               nrq = 4
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            else     !  = final pairing. =

               if(nf_pprovider == 1) then
#ifdef SYNC_COM
!coarray                  call mpi_send(wm_tmp(1,1,1,icxbp0prior), nccp, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipx_pdest, myrank, 
!coarray     &                   mpi_comm_world, istatus, ierr )
!      ndis(me)[ipx_psrc+1] = icxbp0 ! Put
      call xmp_array_section_set_triplet(ndis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_psrc+1
      call xmp_coarray_put_scalar(img_dims,ndis_desc,
     & ndis_sec,icxbp0,status)
!      sync all
      call xmp_sync_all(status)
      nb = ndis(ipx_pdest+1)
!         wm_tmp( :, :, :,          nb:nb         +nd-1 )[ipx_pdest+1]
!     . = wm_tmp( :, :, :, icxbp0prior:icxbp0prior+nd-1 ) ! put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(nb,kind=8),int(nb+nd-1,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxbp0prior,kind=8),int(icxbp0prior+nd-1,kind=8),1,status)

      img_dims(1) = ipx_pdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)

!      sync all
      call xmp_sync_all(status)
!!
#else
                  call mpi_isend(wm_tmp(1,1,1,icxbp0prior), nccp,
     &                   MPI_DOUBLE_COMPLEX, ipx_pdest, myrank, 
     &                   mpi_comm_world, irq(1), ierr)
#endif

               else
#ifdef SYNC_COM
!coarray                  call mpi_recv(wm_tmp(1,1,1,icxbp0), nccp, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipx_psrc, ipx_psrc, 
!coarray     &                   mpi_comm_world, istatus, ierr )
#else
                  call mpi_irecv(wm_tmp(1,1,1,icxbp0), nccp,
     &                   MPI_DOUBLE_COMPLEX, ipx_psrc, ipx_psrc, 
     &                   mpi_comm_world, irq(1), ierr)
#endif

               endif                ! final provider (p)

               if(nf_mprovider == 1) then
#ifdef SYNC_COM
!coarray                  call mpi_send(wm_tmp(1,1,1,icxbm0prior), nccm, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipx_mdest, myrank, 
!coarray     &                   mpi_comm_world, istatus, ierr )
!      mdis(me)[ipx_msrc+1] = icxbm0 ! Put
      call xmp_array_section_set_triplet(mdis_sec,
     & 1,int(me,kind=8),int(me,kind=8),1,status)
      img_dims(1) = ipx_msrc+1
      call xmp_coarray_put_scalar(img_dims,mdis_desc,
     & mdis_sec,icxbm0,status)
!      sync all
      call xmp_sync_all(status)
      mb = mdis(ipx_mdest+1)
!         wm_tmp( :, :, :,          mb:mb         +md )[ipx_mdest+1]
!     . = wm_tmp( :, :, :, icxbm0prior:icxbm0prior+md ) ! Put
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_r_sec,
     & 4,int(mb,kind=8),int(mb+md,kind=8),1,status)

      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 1,int(1,kind=8),int(mylm,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 2,int(1,kind=8),int(lclz,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 3,int(1,kind=8),int(lcly,kind=8),1,status)
      call xmp_array_section_set_triplet(wm_tmp_l_sec,
     & 4,int(icxbm0prior,kind=8),int(icxbm0prior+md,kind=8),1,status)

      img_dims(1) = ipx_mdest+1
      call xmp_coarray_put(img_dims,wm_tmp_desc,wm_tmp_r_sec, 
     & wm_tmp_desc,wm_tmp_l_sec,status)
!      sync all
      call xmp_sync_all(status)
!!
#else
                  call mpi_isend(wm_tmp(1,1,1,icxbm0prior), nccm,
     &                   MPI_DOUBLE_COMPLEX, ipx_mdest, myrank, 
     &                   mpi_comm_world, irq(2), ierr)
#endif

               else
#ifdef SYNC_COM
!coarray                  call mpi_recv(wm_tmp(1,1,1,icxbm0), nccm, 
!coarray     &                   MPI_DOUBLE_COMPLEX,
!coarray     &                   ipx_msrc, ipx_msrc, 
!coarray     &                   mpi_comm_world, istatus, ierr )
#else
                  call mpi_irecv(wm_tmp(1,1,1,icxbm0), nccm,
     &                   MPI_DOUBLE_COMPLEX, ipx_msrc, ipx_msrc, 
     &                   mpi_comm_world, irq(2), ierr)
#endif

               endif                ! final provider (m)
#ifndef SYNC_COM
               nrq = 2
               call mpi_waitall(nrq, irq, istatus, ierr)
#endif

            endif                ! final pairing. !ordinary root

         endif                ! iteration
      end do                  ! iteration

!coarray
      wm = wm_tmp
!!
      return
      end
c---------------------------------------------------------------------
