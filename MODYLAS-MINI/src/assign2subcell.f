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
      subroutine update_wsegc() 
c----------------------------------------------------------------------
      use atommass
      use trj_org
      use trj_mpi
      use unitcell
      use md_fmm
      use md_fmm_domdiv_flg
      use md_periodic
      use md_segment
      use param
      implicit none
      integer(4) :: i,ipar,i0,k0
      real(8) :: totm
      real(8) :: tmpx,tmpy,tmpz

!$omp parallel do default(none)
!$omp& private(i0,k0)
!$omp& private(i,ipar,tmpx,tmpy,tmpz,totm)
!$omp& shared(nselfseg,mass)
!$omp& shared(cellx,celly,cellz,paranum)
!$omp& shared(wseg_cx,wseg_cy,wseg_cz,wkxyz)
!$omp& shared(m2i,lsegtop,lseg_natoms)
      do k0=1,nselfseg
        tmpx=0.0d0
        tmpy=0.0d0
        tmpz=0.0d0
        totm=0.0d0
        do i0 = lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
          i=m2i(i0)
          ipar=paranum(i)
          totm=totm+mass(ipar)
          tmpx=tmpx+wkxyz(1,i0)*mass(ipar)
          tmpy=tmpy+wkxyz(2,i0)*mass(ipar)
          tmpz=tmpz+wkxyz(3,i0)*mass(ipar)
        enddo ! i0
        tmpx=tmpx/totm
        tmpy=tmpy/totm
        tmpz=tmpz/totm
        wseg_cx(k0)=tmpx 
        wseg_cy(k0)=tmpy 
        wseg_cz(k0)=tmpz 
      enddo ! k0

      return
      end
c---------------------------------------------------------------------
      subroutine calc_ia2c()  !! called once in main_f90.f
c---------------------------------------------------------------------
      use trj_org
      use trj_mpi 
      use md_fmm
      use md_fmm_domdiv_flg
      use md_periodic
      use md_segment
      use mpivar
      use unitcell
      implicit none
      integer(4),allocatable :: jtok(:)
      integer(4),allocatable :: ktoj(:),ktoj2(:)
      real(8) :: x0, y0, z0
      integer(4) :: icx, icy, icz
      integer(4) :: icx0, icy0, icz0, icxyz0
      integer(4) :: ipx, ipy, ipz, ipxyz
      integer(4) :: i, k, ll
      integer(4) :: k0,i0
        real(8) :: rdcellx,rdcelly,rdcellz
        real(8) :: rxdiv,rydiv,rzdiv
      integer(4) :: k00,k000,ii,ic

      allocate(ktoj(max_seg))  ! deallocate at calc_ia2c
      allocate(ktoj2(max_seg)) ! deallocate at calc_ia2c

      allocate(jtok(nsegments))
      rdcellx=dble(ncell)/cellx
      rdcelly=dble(ncell)/celly
      rdcellz=dble(ncell)/cellz
      rxdiv  =dble(nxdiv)/cellx
      rydiv  =dble(nydiv)/celly
      rzdiv  =dble(nzdiv)/cellz

      ndseg_fmmn = 0
      ktoj =-1
      jtok =-1
      m2i = 0

      i0 = 0 ! local atom id
      ll = 0 ! local segment id
      do k = 1, nsegments
        x0 = seg_cx(k)+cellxh
        y0 = seg_cy(k)+cellyh
        z0 = seg_cz(k)+cellzh
        icx=min0( int(x0*rdcellx),ncell-1 )+1
        icy=min0( int(y0*rdcelly),ncell-1 )+1
        icz=min0( int(z0*rdcellz),ncell-1 )+1

        if(idcell(icx,icy,icz).ne.0) then
          ipx=min0( int(x0*rxdiv),nxdiv-1 )+1
          ipy=min0( int(y0*rydiv),nydiv-1 )+1
          ipz=min0( int(z0*rzdiv),nzdiv-1 )+1
          ipxyz=idom(ipx,ipy,ipz)

        if(ipxyz.eq.myrank) then
          icx0=icx-ixmin(myrank)+1 ! global -> local
          icy0=icy-iymin(myrank)+1 ! global -> local
          icz0=icz-izmin(myrank)+1 ! global -> local
          icxyz0 = (icx0-1)*(lydiv*lzdiv) + (icy0-1)*lzdiv + (icz0-1)
          ndseg_fmmn(icz0,icy0,icx0) = ndseg_fmmn(icz0,icy0,icx0) + 1
        if (ndseg_fmmn(icz0,icy0,icx0) .gt. max_nsegments_per_cell) then
            write(*,*) ndseg_fmmn(icz0,icy0,icx0),max_nsegments_per_cell
            call abort_with_message_a('Error: ndseg_fmmn overflow', 0)
        endif
          k0=max_nsegments_per_cell*icxyz0+ndseg_fmmn(icz0,icy0,icx0)
          wseg_cx(k0) = seg_cx(k)
          wseg_cy(k0) = seg_cy(k)
          wseg_cz(k0) = seg_cz(k)
          ktoj2(k0)=k
          ll=ll+1

          do i = segtop(k)+1,segtop(k)+seg_natoms(k)
            i0 = i0 + 1
            if(i0.gt.nadirect) stop'ERROR : i0 overflow in calc_ia2c'
            i2m(i)=i0
            m2i(i0)=i
          wkxyz(1,i0) = xyz(1,i)
          wkxyz(2,i0) = xyz(2,i)
          wkxyz(3,i0) = xyz(3,i)
            wkv(1,i0) = v(1,i)
            wkv(2,i0) = v(2,i)
            wkv(3,i0) = v(3,i)
          end do

        end if
        end if
      enddo
      nselfatm = i0
      nselfseg = ll

!^^^ Left-align wseg ^^^!
      k000=0
      do ii=1,lxdiv*lydiv*lzdiv
        icz0=mod(ii-1,lzdiv)     +1   
        icy0=mod(ii-1,lzdiv*lydiv)
        icy0=icy0/lzdiv          +1   
        icx0=(ii-1)/(lzdiv*lydiv)+1
        icxyz0 = (icx0-1)*(lzdiv*lydiv) + (icy0-1)*lzdiv + (icz0-1)
        k00 = max_nsegments_per_cell*icxyz0
        do ic=1,ndseg_fmmn(icz0,icy0,icx0)
          k000=k000+1
          k0=k00+ic
          seg_cx(k000)=wseg_cx(k0)
          seg_cy(k000)=wseg_cy(k0)
          seg_cz(k000)=wseg_cz(k0)
          k=ktoj2(k0)
          ktoj(k000)=k
          jtok(k)=k000
        enddo ! ic
      enddo ! ii

!^^^ reset wseg_cx ^^^!
      wseg_cx=0d0
      wseg_cy=0d0
      wseg_cz=0d0
      do k0=1,nselfseg
        wseg_cx(k0)=seg_cx(k0)
        wseg_cy(k0)=seg_cy(k0)
        wseg_cz(k0)=seg_cz(k0)
        k = ktoj(k0)
        i = segtop(k)+1 
        i0 = i2m(i)    
        lsegtop(k0)    =i0 
        lseg_natoms(k0)=seg_natoms(k)
      enddo

      deallocate(ktoj,ktoj2)
      deallocate(jtok)
      deallocate(idcell)
      deallocate(idom)
      return
      end
c----------------------------------------------------------------------
      subroutine atom2cell
c----------------------------------------------------------------------
      use md_fmm
      use md_fmm_domdiv_flg
      use md_segment
      use md_periodic
      use trj_org
      use trj_mpi
      use mpivar
      implicit none
      real(8),allocatable :: wkx2(:), wky2(:), wkz2(:), wkv2(:,:)
      integer(4),allocatable :: m2i2(:)
      integer(4) :: i, ii,ic, i0, i00, iz,iy,ix
      integer(4) :: j0x,j0y,j0z
      integer(4) :: icx0,icy0,icz0
      integer(4) :: jx,jy,jz
      include 'mpif.h' 
!     integer(4) :: ntmp
      integer(4) :: k0

      allocate(wkx2(nadirect))
      allocate(wky2(nadirect))
      allocate(wkz2(nadirect))
      allocate(wkv2(3,nadirect))
      allocate(m2i2(nadirect))

!### initialize ###    
!$omp parallel default(shared)
!$omp& private(ix,iy,iz,i,i0)
!$omp do
      do ix=1,lxdiv+4
        do iy=1,lydiv+4
          do iz=1,lzdiv+4
        tag(iz,iy,ix)=0
        na_per_cell(iz,iy,ix)=0
          enddo
        enddo
      enddo
!$omp end do nowait
!$omp do
      do i=1,n
        i2m(i) =-1       ! initialize
      enddo
!$omp end do
!$omp do
      do i0=1,nadirect
        wkx2(i0)=0d0
        wky2(i0)=0d0
        wkz2(i0)=0d0
        wkv2(1:3,i0)=0d0
        m2i2(i0)=m2i(i0) ! store m2i
        m2i(i0)=-1       ! initialize
      enddo
!$omp end do
!$omp end parallel

!### self range ###    
      k0=0
!		ii=0
      do j0x=1,lxdiv+4
        jx=ixmin(myrank)-3+j0x
        do j0y=1,lydiv+4
          jy=iymin(myrank)-3+j0y
          do j0z=1,lzdiv+4
            jz=izmin(myrank)-3+j0z
        IF(j0x.ge.3.and.j0x.le.lxdiv+2 .and.
     &     j0y.ge.3.and.j0y.le.lydiv+2 .and.
     &     j0z.ge.3.and.j0z.le.lzdiv+2)THEN
          if(j0z==3)THEN
          i00=(j0x-1)*na1cell*(lydiv+4)*(lzdiv+4)
     &       +(j0y-1)*na1cell          *(lzdiv+4)
     &       +(j0z-1)*na1cell
          endif
        ELSE
          cycle ! outside of myrank
        ENDIF

        icx0 = j0x-2
        icy0 = j0y-2
        icz0 = j0z-2
      IF(ndseg_fmmn(icz0,icy0,icx0)==0)then
        write(*,*) 'WARNNING: ndseg_fmmn=0',
     &  icx0+ixmin(myrank)-1,icy0+iymin(myrank)-1,icz0+izmin(myrank)-1
      ENDIF
        tag(j0z,j0y,j0x)=i00+1
        do ic = 1, ndseg_fmmn(icz0,icy0,icx0)
          k0=k0+1 ! local
          do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
            i00=i00+1 ! local
            m2i(i00)=m2i2(i0)
            wkx2(i00)    =wkxyz(1,i0)
            wky2(i00)    =wkxyz(2,i0)
            wkz2(i00)    =wkxyz(3,i0)
            wkv2(1:3,i00)=wkv(1:3,i0)
          enddo ! i0
          lsegtop(k0)=i00-lseg_natoms(k0)+1 ! renew lsegtop
        enddo ! ic

      IF(ndseg_fmmn(icz0,icy0,icx0)==0)then
        na_per_cell(j0z,j0y,j0x)=0
      ELSE
        na_per_cell(j0z,j0y,j0x)=i00-tag(j0z,j0y,j0x)+1
        if(na_per_cell(j0z,j0y,j0x).gt.na1cell)then
          write(*,*) 'ERROR: na_per_cell over flowed!'
          call mpistop()
        endif
      ENDIF

          enddo ! jz
        enddo ! jy
      enddo ! jx

!$omp parallel default(shared)
!$omp& private(ii,i0,i,icx0,icy0,icz0)
!$omp do
      do ii=1,lxdiv*lydiv*lzdiv
      icz0=mod(ii-1,lzdiv)     +3   
      icy0=mod(ii-1,lzdiv*lydiv)
      icy0=icy0/lzdiv          +3   
      icx0=(ii-1)/(lzdiv*lydiv)+3
      do i0=tag(icz0,icy0,icx0),
     &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
        i=m2i(i0)
        i2m(i)=i0
        wkxyz(1,i0)=wkx2(i0)
        wkxyz(2,i0)=wky2(i0)
        wkxyz(3,i0)=wkz2(i0)
        wkv(1:3,i0)=wkv2(1:3,i0)
      enddo ! i0
      enddo ! ii
!$omp end do
!$omp end parallel

      deallocate(wkx2,wky2,wkz2,wkv2,m2i2)

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_pbc()
c-----------------------------------------------------------------------
      use trj_mpi
      use md_fmm
      use md_fmm_domdiv_flg
      use md_periodic
      use md_segment
      use param
      use mpivar
      use unitcell
      implicit none
      integer(4)::j0
      integer(4)::jxb,jyb,iz,jzb
      integer(4)::ncellx,ncelly,ncellz
 
      ncellx=ncell
      ncelly=ncell
      ncellz=ncell

!$omp parallel default(shared)
!$omp& private(jxb,jyb,jzb,iz,j0)
      do jxb=1,ncellx/nxdiv+4
         do jyb=1,ncelly/nydiv+4
            do iz =3,ncellz/nzdiv+2
               jzb=iz

!++++++++++++  Boundary phase start
!
!------------  Z-direction
               if(jzb==3 .and. izmin(myrank)==1) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(2,jyb,jxb)
     &                  + na_per_cell(2,jyb,jxb)-1
                     wkxyz(3,j0)=wkxyz(3,j0)-cellz
                  enddo
!$omp end do
               endif
               if(lzdiv .eq. 1) then
               if(jzb==3 .and. izmin(myrank)==2) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(1,jyb,jxb)
     &                  + na_per_cell(1,jyb,jxb)-1
                     wkxyz(3,j0)=wkxyz(3,j0)-cellz
                  enddo
!$omp end do
               endif
               endif
               if(jzb==ncellz/nzdiv+1 .and. izmax(myrank)==ncellz) then
!$omp do
                  do j0=tag(ncellz/nzdiv+3,jyb,jxb),
     &                  tag(ncellz/nzdiv+4,jyb,jxb)
     &                  + na_per_cell(ncellz/nzdiv+4,jyb,jxb)-1
                     wkxyz(3,j0)=wkxyz(3,j0)+cellz
                  enddo
!$omp end do
               endif
               if(lzdiv .eq. 1) then
               if(izmax(myrank)==ncellz) then
!$omp do
                  do j0=tag(4,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(3,j0)=wkxyz(3,j0)+cellz
                  enddo
!$omp end do
               endif
               if(izmax(myrank)==ncellz-1) then
!$omp do
                  do j0=tag(ncellz/nzdiv+4,jyb,jxb),
     &                  tag(ncellz/nzdiv+4,jyb,jxb)
     &                  + na_per_cell(ncellz/nzdiv+4,jyb,jxb)-1
                     wkxyz(3,j0)=wkxyz(3,j0)+cellz
                  enddo
!$omp end do
               endif
               endif

!------------  Y-direction
               if(jzb==3 .and. jyb<=2 .and. iymin(myrank)==1) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)-celly
                  enddo
!$omp end do
               else if(jyb<=2 .and. iymin(myrank)==1) then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)-celly
                  enddo
!$omp end do
               endif
               if(lydiv .eq. 1) then
               if(jzb==3 .and. jyb<=1 .and. iymin(myrank)==2) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)-celly
                  enddo
!$omp end do
               else if(jyb<=1 .and. iymin(myrank)==2) then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)-celly
                  enddo
!$omp end do
               endif
               endif
               if(jzb==3.and.jyb>=ncelly/nydiv+3.and.
     &            iymax(myrank)==ncelly) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)+celly
                  enddo
!$omp end do
               elseif(jyb>=ncelly/nydiv+3.and.iymax(myrank)==ncelly)then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)+celly
                  enddo
!$omp end do
               endif
               if(lydiv .eq. 1) then
               if(jzb==3.and.jyb>=ncelly/nydiv+4.and.
     &            iymax(myrank)==ncelly-1) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)+celly
                  enddo
!$omp end do
               else if(jyb>=ncelly/nydiv+4.and.
     &                 iymax(myrank)==ncelly-1) then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(2,j0)=wkxyz(2,j0)+celly
                  enddo
!$omp end do
               endif
               endif

!------------  X-direction
               if(jzb==3 .and. jxb<=2 .and. ixmin(myrank)==1) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)-cellx
                  enddo
!$omp end do
               else if(jxb<=2 .and. ixmin(myrank)==1) then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)-cellx
                  enddo
!$omp end do
               endif
               if(lxdiv .eq. 1) then
               if(jzb==3 .and. jxb<=1 .and. ixmin(myrank)==2) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)-cellx
                  enddo
!$omp end do
               else if(jxb<=1 .and. ixmin(myrank)==2) then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)-cellx
                  enddo
!$omp end do
               endif
               endif
               if(jzb==3.and.jxb>=ncellx/nxdiv+3.and.
     &            ixmax(myrank)==ncellx) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)+cellx
                  enddo
!$omp end do
               elseif(jxb>=ncellx/nxdiv+3.and.ixmax(myrank)==ncellx)then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)+cellx
                  enddo
!$omp end do
               endif
               if(lxdiv .eq. 1) then
               if(jzb==3.and.jxb>=ncellx/nxdiv+4.and.
     &            ixmax(myrank)==ncellx-1) then
!$omp do
                  do j0=tag(1,jyb,jxb),
     &                  tag(5,jyb,jxb)
     &                  + na_per_cell(5,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)+cellx
                  enddo
!$omp end do
               else if(jxb>=ncellx/nxdiv+4.and.
     &            ixmax(myrank)==ncellx-1) then
!$omp do
                  do j0=tag(jzb+2,jyb,jxb),
     &                  tag(jzb+2,jyb,jxb)
     &                  + na_per_cell(jzb+2,jyb,jxb)-1
                     wkxyz(1,j0)=wkxyz(1,j0)+cellx
                  enddo
!$omp end do
               endif
               endif
!
!++++++++++++  Boundary phase ended
!
           enddo !iz
         enddo !jyb
      enddo !jxb
!$omp end parallel

      return
      end
