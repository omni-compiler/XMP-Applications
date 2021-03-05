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
      subroutine init_fmm_domain_div()
c----------------------------------------------------------------------
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      implicit none
      integer(4) :: iwkx,iwky,iwkz
      integer(4) :: maxdiv,icnt
      integer(4) :: icx, icy, icz, icxyz1
      integer(4) :: idx, idy, idz, imx, imy, imz,i
      include 'mpif.h'
      integer(4) :: ierr
      integer(4) :: status

      maxdiv = nprocs
!ya
      IF(mpi_manual_division_flg)THEN
        continue
      ELSE
        nxdiv=1
        nydiv=1
        nzdiv=1
      ENDIF
!ya
      if(nprocs.ne.1) then
!ya
        IF(mpi_manual_division_flg)THEN
          maxdiv=1
          continue
        ELSE
          do while (mod(maxdiv,2).eq.0)
          nxdiv = nxdiv*2
          maxdiv = maxdiv/2
          if(mod(maxdiv,2).eq.0) then
            nydiv = nydiv*2
            maxdiv = maxdiv/2
          end if
          if(mod(maxdiv,2).eq.0) then
            nzdiv = nzdiv*2
            maxdiv = maxdiv/2
          end if
          end do
        ENDIF
!ya
      end if

      if(((nxdiv*nydiv*nzdiv).ne.nprocs) .or.
     &   (maxdiv.ne.1) ) then
        write(*,*) ' -error init_fmm_domain_div 1-'
        write(*,*) 'The Number of MPI-procs is incorrect. '
        write(*,*) 'nxdiv,nydiv,nzdiv,nprocs= ',nxdiv,nydiv,nzdiv,nprocs
        call mpi_abort(mpi_comm_world,ierr)
      end if

      maxdiv = max(maxdiv,nxdiv)
      maxdiv = max(maxdiv,nydiv)
      maxdiv = max(maxdiv,nzdiv)

      if(maxdiv.gt.ncell) then
        write(*,*) ' -error init_fmm_domain_div 2-'
        write(*,*) ' maxdiv.gt.ncell'
        call mpi_abort(mpi_comm_world,ierr)
      endif

      lxdiv = mod(ncell,nxdiv) 
      lydiv = mod(ncell,nydiv) 
      lzdiv = mod(ncell,nzdiv) 

      if((lxdiv.ne.0).or.(lydiv.ne.0).or.(lzdiv.ne.0)) then
        write(*,*) ' -error init_fmm_domain_div 3-'
        write(*,*) ' mod(ncell,mpidiv).ne.0'
        call mpi_abort(mpi_comm_world,ierr)
      endif

      lxdiv = ncell/nxdiv
      lydiv = ncell/nydiv
      lzdiv = ncell/nzdiv

      allocate(ixflg(0:nprocs-1))
      allocate(iyflg(0:nprocs-1))
      allocate(izflg(0:nprocs-1))
      allocate(ixmax(0:nprocs-1))
      allocate(iymax(0:nprocs-1))
      allocate(izmax(0:nprocs-1))
      allocate(ixmin(0:nprocs-1))
      allocate(iymin(0:nprocs-1))
      allocate(izmin(0:nprocs-1))
!
      allocate(idom(nxdiv,nydiv,nzdiv))    ! deallocate at calc_ia2c
      allocate(idcell(ncell,ncell,ncell))  ! deallocate at calc_ia2c
      allocate(nd2c(0:nxdiv*nydiv*nzdiv-1)) ! deallocate at calc_idcell
      allocate(id2c(lxdiv*lydiv*lzdiv,0:nxdiv*nydiv*nzdiv-1)) ! dealoc

      icnt=0
      do iwkz=1,nzdiv
        do iwky=1,nydiv
          do iwkx=1,nxdiv
            ixflg(icnt)=iwkx
            iyflg(icnt)=iwky
            izflg(icnt)=iwkz
            idom(iwkx,iwky,iwkz) = icnt
            icnt=icnt+1
          end do
        end do
      end do

      ixmax=0
      iymax=0
      izmax=0
      ixmin=ncell
      iymin=ncell
      izmin=ncell
!coarray      call mpi_barrier(mpi_comm_world,ierr)
      !sync all
      call xmp_sync_all(status)
!!

      if(myrank==0)then
        IF(mpi_manual_division_flg)THEN
          write(*,'(/,a)') '**** MPI manual division'
        ELSE
          write(*,'(/,a)') '**** MPI auto division'
        ENDIF
      endif

      nd2c = 0
      do icz = 1, ncell
        do icy = 1, ncell
          do icx = 1, ncell
            icxyz1 = (icx-1) * ncell*ncell 
     &     + (icy-1) * ncell + icz-1
            imx = mod(icx,lxdiv)             
            imy = mod(icy,lydiv)             
            imz = mod(icz,lzdiv)             
            idx = int(icx/lxdiv) 
            idy = int(icy/lydiv)
            idz = int(icz/lzdiv)
            if(imx.ne.0) idx = idx + 1
            if(imy.ne.0) idy = idy + 1
            if(imz.ne.0) idz = idz + 1
            i = idom(idx,idy,idz)
            ixmax(i) = max(ixmax(i),icx)
            iymax(i) = max(iymax(i),icy)
            izmax(i) = max(izmax(i),icz)
            ixmin(i) = min(ixmin(i),icx)
            iymin(i) = min(iymin(i),icy)
            izmin(i) = min(izmin(i),icz)
            nd2c(i) = nd2c(i)+1
            id2c(nd2c(i),i) = icxyz1
          end do
        end do
      end do

      return
      end
c----------------------------------------------------------------------
      subroutine check_cutofflength               ! YA
c----------------------------------------------------------------------
      use cutoffradius
      use md_fmm_domdiv_flg
      use unitcell
      implicit none
      real(8) :: cellxd,cellyd,cellzd
      
      cellxd=cellx/dble(lxdiv)
      cellyd=celly/dble(lydiv)
      cellzd=cellz/dble(lzdiv)
      cellxd=2d0*cellxd
      cellyd=2d0*cellyd
      cellzd=2d0*cellzd
      if(cutrad.eq.0)then
        stop'LJ-cutoff length=0, which is unlikely situation.'
      endif
      if(cellxd.lt.cutrad)then  ! Cubic cell only
        write(6,*) 'Length of LEVEL:0 cell=', cellxd*1d+10, ' Aungstrom'
        write(6,*) 'LJ-cutoff length=', cutrad*1d+10, ' Aungstrom'
        write(6,*) 'This situation is unlikely, so change ncell value.'
        stop
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine fmod_alloc_multipole
c----------------------------------------------------------------------
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use ompvar
      implicit none
!ya
      integer(4) :: m1
      integer(4) :: ieo_zst, ieo_yst, ieo_xst 
      integer(4) :: ieo_zen, ieo_yen, ieo_xen 
      integer(4) :: npz, npy, npx ! process number along each axis
      integer(4) :: ipz, ipy, ipx ! iflgx,iflgy,iflgz
      integer(4) :: mcell_size
      integer(4) :: nbound_zm, nbound_ym, nbound_xm
      integer(4) :: nbound_zp, nbound_yp, nbound_xp
      integer(4) :: nsczdiv,nscydiv,nscxdiv
      integer(4) :: lclz,lcly,lclx
      integer(4) :: nscell
      integer(4) :: tgtl
      type(fmm_data_t),pointer :: d
!ya
      m1=(nmax+1)*(nmax+1)
      ipx=ixflg(myrank)-1 ; npx=nxdiv
      ipy=iyflg(myrank)-1 ; npy=nydiv
      ipz=izflg(myrank)-1 ; npz=nzdiv

      if (nlevel < 3) then
        if (myrank == 0) write(*,*) 'ERROR: nlevel < 3'
        call mpistop()
      endif 

      allocate(fmm_data(0:nlevel))

      do tgtl = 0, nlevel
        d => fmm_data(tgtl)
        mcell_size = 2**tgtl
        nsczdiv = (ncell / mcell_size-1) / npz + 1
        nscydiv = (ncell / mcell_size-1) / npy + 1
        nscxdiv = (ncell / mcell_size-1) / npx + 1
        ieo_zst = mod((ncell / mcell_size * ipz) / npz + 1,2)
        nbound_zm = 4 ; if(ieo_zst == 0) nbound_zm = 5
        ieo_zen = mod((ncell/mcell_size * (ipz+1) -1)/npz + 1, 2)
        nbound_zp = 4 ; if(ieo_zen == 1) nbound_zp = 5
        ieo_yst = mod((ncell / mcell_size * ipy) / npy + 1,2)
        nbound_ym = 4 ; if(ieo_yst == 0) nbound_ym = 5
        ieo_yen = mod((ncell/mcell_size * (ipy+1) -1)/npy + 1, 2)
        nbound_yp = 4 ; if(ieo_yen == 1) nbound_yp = 5
        ieo_xst = mod((ncell / mcell_size * ipx) / npx + 1,2)
        nbound_xm = 4 ; if(ieo_xst == 0) nbound_xm = 5
        ieo_xen = mod((ncell/mcell_size * (ipx+1) -1)/npx + 1, 2)
        nbound_xp = 4 ; if(ieo_xen == 1) nbound_xp = 5

        nscell = (ncell - 1) / mcell_size + 1 

        lclx = nbound_xm+nscxdiv+nbound_xp
        lcly = nbound_ym+nscydiv+nbound_yp
        lclz = nbound_zm+nsczdiv+nbound_zp

        d%lclx=lclx; d%lcly=lcly; d%lclz=lclz
        d%mcell_size=mcell_size; d%nscell=nscell
        d%nscxdiv=nscxdiv; d%nscydiv=nscydiv; d%nsczdiv=nsczdiv
        d%nbound_xm=nbound_xm; d%nbound_xp=nbound_xp
        d%nbound_ym=nbound_ym; d%nbound_yp=nbound_yp
        d%nbound_zm=nbound_zm; d%nbound_zp=nbound_zp

        allocate(d%wm_local(m1,lclz,lcly,lclx))
        allocate(d%wm_global(m1,nscell,nscell,nscell))
        allocate(d%wl_local(m1,nsczdiv,nscydiv,nscxdiv))
        allocate(d%wwl_local(m1,nsczdiv,nscydiv,nscxdiv,0:nomp-1))

        if (tgtl == 0) then
          allocate(wwm_local0(m1,lclz,lcly,lclx,0:nomp-1))
        endif

      end do

      return
      end
c----------------------------------------------------------------------
      subroutine calc_idcell()
c----------------------------------------------------------------------
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      implicit none
      include 'mpif.h'
      integer(4) :: icx0, icy0, icz0
      integer(4) :: icx1, icy1, icz1
      integer(4) :: ich1,k,ih1,ntmp
      integer(4) :: ndirect2,mrcsafe,ndfmm0
      integer(4) :: ix,iy,iz
      integer(4) :: ix2,iy2,iz2
      integer(4),allocatable ::
     &  idfmm0x(:),idfmm0y(:),idfmm0z(:)

      ndirect2 = ndirect + ndcellmargin
      if (ndirect2 .gt. nimage-1) then
        mrcsafe = nimage-1
      else
        mrcsafe = ndirect2
      endif
      ndfmm0 = 0
      do ix = -mrcsafe, mrcsafe
        ix2 = ix * ix
        do iy = -mrcsafe, mrcsafe
          iy2 = iy * iy
          do iz = -mrcsafe, mrcsafe
            iz2 = iz * iz
#ifdef COMM_CUBE
            if (abs(ix) <= ndirect2 .and.
     &          abs(iy) <= ndirect2 .and.
     &          abs(iz) <= ndirect2) then
#else
            if (ix2+iy2+iz2 .le. ndirect2*ndirect2) then
#endif
              ndfmm0 = ndfmm0 + 1
            endif
          enddo
        enddo
      enddo
      allocate(idfmm0x(ndfmm0))
      allocate(idfmm0y(ndfmm0))
      allocate(idfmm0z(ndfmm0))
      idfmm0x = 0
      idfmm0y = 0
      idfmm0z = 0
      ntmp = 0
      do ix = -mrcsafe, mrcsafe
        ix2 = ix * ix
        do iy = -mrcsafe, mrcsafe
          iy2 = iy * iy
          do iz = -mrcsafe, mrcsafe
            iz2 = iz * iz
#ifdef COMM_CUBE
            if (abs(ix) <= ndirect2 .and.
     &          abs(iy) <= ndirect2 .and.
     &          abs(iz) <= ndirect2) then
#else
            if (ix2+iy2+iz2 .le. ndirect2*ndirect2) then
#endif
              ntmp = ntmp + 1
              idfmm0x(ntmp) = ix
              idfmm0y(ntmp) = iy
              idfmm0z(ntmp) = iz
            endif
          enddo
        enddo
      enddo

      idcell = 0
      ntmp   = 0
      do ih1 = 1,nd2c(myrank)
        ich1 = id2c(ih1,myrank)
        icz0 = mod(ich1,ncell)+1
        icy0 = ich1/ncell
        icx0 = icy0/ncell+1
        icy0 = mod(icy0,ncell)+1
        do k = 1,ndfmm0
          icx1 = icx0 + idfmm0x(k)
          icy1 = icy0 + idfmm0y(k)
          icz1 = icz0 + idfmm0z(k)
          if (icx1 .gt. 3*ncell .or. icx1 .le. -2*ncell .or.
     &        icy1 .gt. 3*ncell .or. icy1 .le. -2*ncell .or.
     &        icz1 .gt. 3*ncell .or. icz1 .le. -2*ncell) then
              cycle
          else
            do while (icx1 .gt. ncell)
              icx1 = icx1 - ncell
            end do
            do while (icx1 .le. 0)
              icx1 = icx1 + ncell
            end do
            do while (icy1 .gt. ncell)
              icy1 = icy1 - ncell
            end do
            do while (icy1 .le. 0)
              icy1 = icy1 + ncell
            end do
            do while (icz1 .gt. ncell)
              icz1 = icz1 - ncell
            end do
            do while (icz1 .le. 0)
              icz1 = icz1 + ncell
            end do
            if(idcell(icx1,icy1,icz1).eq.0) then
              ntmp = ntmp + 1
              idcell(icx1,icy1,icz1) = ntmp
            end if
          end if
        end do
      end do

      ndcell = ntmp

      deallocate(idfmm0x)
      deallocate(idfmm0y)
      deallocate(idfmm0z)
!ya
      deallocate(nd2c,id2c)
!ya

      return
      end
c---------------------------------------------------------------------
      subroutine fmod_alloc_metadata !init_fmm_direct_com()
c----------------------------------------------------------------------
      use trj_org
      use trj_mpi
      use md_forces
      use md_periodic
      use md_fmm
      use md_fmm_domdiv_flg
      use shakerattleroll
      use md_multiplestep
      use md_condition
      use md_segment
      use shakerattleroll
      use mpivar
      use ompvar
      implicit none
      integer(4) :: itmp
      include 'mpif.h'

!############
!  metadata
!############
      allocate(tag(lzdiv+4,lydiv+4,lxdiv+4))
      allocate(na_per_cell(lzdiv+4,lydiv+4,lxdiv+4)[*])
!############
!  segment
!############
      max_seg = max_nsegments_per_cell*lzdiv*lydiv*lxdiv
      allocate(wseg_cz(max_seg))
      allocate(wseg_cy(max_seg))
      allocate(wseg_cx(max_seg))
      allocate(ndseg_fmmn(lzdiv,lydiv,lxdiv))
      itmp=max_nsegments_per_cell*(lxdiv+4)*(lydiv+4)*(lzdiv+4)
      allocate(lsegtop(itmp))
      allocate(lseg_natoms(itmp))
!############
!  atom
!############
      na1cell=max( int((n/ncell**3)*na1cellmargin),10 )
      na5cell=na1cell*5
      nadirect=na1cell*(lxdiv+4)*(lydiv+4)*(lzdiv+4)
!Coordinate & Velocity
      allocate(wkxyz(3,nadirect)[*])
      allocate(wkv(3,nadirect))
      allocate(m2i(nadirect)[*])
!Force
      allocate(wk_f(3,nadirect))
      allocate(w3_f(3,nadirect,0:nomp-1))
!TABLE for Force_calculation
      allocate(chgv_table(na5cell,0:nomp-1))
      allocate(epsilon_sqrt_table(na5cell,0:nomp-1))
      allocate(R_half_table(na5cell,0:nomp-1))
!MT
      allocate(fshort(3,nadirect))
      allocate(fmiddle(3,nadirect))
      allocate(flong(3,nadirect))
!SHAKE
      allocate(xyzstr(3,nadirect))
      allocate(kshake(nadirect))

      return
      end
c----------------------------------------------------------------------
