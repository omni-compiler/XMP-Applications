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
      subroutine calc_mm
c---------------------------------------------------------------------
      use trj_org
      use comm_base
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use unitcell
#include "timing.h90"
      implicit none
      integer(4) :: il, il0
      type(fmm_data_t),pointer :: d
      integer nbsize

      TIME_START(TM_P2M)
      call calc_fmm
      TIME_STOP(TM_P2M)

      DO il = 1, nlevel

        il0 = il-1
        d => fmm_data(il0)

        nbsize = max(d%nbound_zm, d%nbound_zp)

        TIME_START(TM_COMM_FMM)
        if(il0.lt.lgflg+1)then
          call comm_fmm_local_multi(il0, (nmax+1)*(nmax+1),
     &           d%wm_local, d%lclz, d%lcly,d %lclx,
     &           nbsize, d%nscydiv, d%nscxdiv)
        else
          call comm_fmm_local_top(il0, (nmax+1)*(nmax+1),
     &           d%wm_global, d%nscell, 
     &           d%nsczdiv, d%nscydiv, d%nscxdiv)
        endif
        TIME_STOP(TM_COMM_FMM)

        TIME_START(TM_M2M)
        call merge_fmm(il0, il)
        TIME_STOP(TM_M2M)

      ENDDO ! il

      return
      end
c---------------------------------------------------------------------
      subroutine calc_fmm 
c---------------------------------------------------------------------
      use trj_org
      use trj_mpi
      use md_coulomb
      use md_fmm
      use md_segment
      use md_periodic
      use md_fmm_domdiv_flg
      use param
      use comm_base
      use mpivar
      use ompvar
      use unitcell
      implicit none
      real(8) :: rnc, r_xbox, r_ybox, r_zbox
      real(8) :: x0, y0, z0, q0, rx, ry, rz
      integer(4) :: icx, icy, icz
      integer(4) :: icx0, icy0, icz0
      integer(4) :: i,ii
      integer(4) :: j,k, m1
      real(8) :: xta, yta, zta, rad, the, csthe, phi
      real(8) :: f, g
      real(8) :: algndr
      complex(8) :: zphi
      integer(4) :: izst,iyst,ixst
      integer(4) :: i0,ipar,iam=0
      integer(4) :: icx00,icy00,icz00 
      integer(4) :: mcell_size
      complex(8),pointer :: wl_local0(:,:,:,:)
      complex(8),pointer :: wm_local0(:,:,:,:)
      integer(4) :: nbound_zm, nbound_ym, nbound_xm
      integer(4) :: nbound_zp, nbound_yp, nbound_xp

!$    include 'omp_lib.h'

      wl_local0 => fmm_data(0)%wl_local
      wm_local0 => fmm_data(0)%wm_local

      rnc = 1.0d0 / dble(ncell)
      r_xbox = 1.0d0 / cellx
      r_ybox = 1.0d0 / celly
      r_zbox = 1.0d0 / cellz

      mcell_size=fmm_data(0)%mcell_size
      nbound_zm=fmm_data(0)%nbound_zm
      nbound_zp=fmm_data(0)%nbound_zp
      nbound_ym=fmm_data(0)%nbound_ym
      nbound_yp=fmm_data(0)%nbound_yp
      nbound_xm=fmm_data(0)%nbound_xm
      nbound_xp=fmm_data(0)%nbound_xp

      izst = izmin(myrank)
      iyst = iymin(myrank)
      ixst = ixmin(myrank)

      iam = 0

!$omp parallel default(none)
!$omp& private(iam) 
!$omp& private(ii,icx,icy,icz) 
!$omp& private(icx0,icy0,icz0,icx00,icy00,icz00)
!$omp& private(i,x0,y0,z0)
!$omp& private(q0,rx,ry,rz,xta,yta,zta,rad,the,csthe,phi,j,f,k)
!$omp& private(zphi,g,m1)
!$omp& private(i0,ipar)
!$omp& shared(nmax)
!$omp& shared(wkxyz,chgv,pre)
!$omp& shared(rnc,cellx,celly,cellz,paranum)
!$omp& shared(izst,ixst,iyst)
!$omp& shared(tag,na_per_cell,m2i)
!$omp& shared(lxdiv,lydiv,lzdiv)
!$omp& shared(nbound_xm,nbound_ym,nbound_zm)
!$omp& shared(wm_local0,wwm_local0,nomp)
!$    iam=omp_get_thread_num()
CC!$omp do
      do ii=1,lxdiv*lydiv*lzdiv
         icz0=mod(ii-1,lzdiv)     +3   
         icy0=mod(ii-1,lzdiv*lydiv)
         icy0=icy0/lzdiv          +3   
         icx0=(ii-1)/(lzdiv*lydiv)+3
          icx=icx0+ixst-3 ! local -> global
          icy=icy0+iyst-3 ! local -> global
          icz=icz0+izst-3 ! local -> global
        icx00=icx0-2+nbound_xm ! local -> local FMM
        icy00=icy0-2+nbound_ym ! local -> local FMM
        icz00=icz0-2+nbound_zm ! local -> local FMM
!       wm_local0(:,icz00,icy00,icx00)=0d0
        wwm_local0(:,icz00,icy00,icx00,iam)=0d0
!$omp do
            do i0=tag(icz0,icy0,icx0),
     &            tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
              x0 = wkxyz(1,i0)
              y0 = wkxyz(2,i0)
              z0 = wkxyz(3,i0)
              i = m2i(i0) 
              ipar = paranum(i)
              q0 = chgv(ipar)
              rx = x0 - ((icx-0.5d0) * rnc - 0.5d0) * cellx
              ry = y0 - ((icy-0.5d0) * rnc - 0.5d0) * celly
              rz = z0 - ((icz-0.5d0) * rnc - 0.5d0) * cellz
              xta=rx
              yta=ry
              zta=rz
              call cart2angle(xta,yta,zta,rad,the,csthe,phi)
**** calculate m-matrices
              do j=0,nmax
                f=q0*rad**j
                do k=-j,j
                  zphi=dcmplx(0.d0,k*phi)
                  g=f*pre(j,iabs(k))*algndr(j,iabs(k),csthe)
                  m1 = j*j+j+1+k
!                 wm_local0(m1,icz00,icy00,icx00)
!    $           =wm_local0(m1,icz00,icy00,icx00)+g*exp(-zphi)
                  wwm_local0(m1,icz00,icy00,icx00,iam)
     $           =wwm_local0(m1,icz00,icy00,icx00,iam)+g*exp(-zphi)
                enddo ! k
              enddo ! j
            enddo ! i0
!$omp end do
      enddo ! ii
CC!$omp end do
!$omp do
      do ii=1,lxdiv*lydiv*lzdiv
        icz0=mod(ii-1,lzdiv)     +3
        icy0=mod(ii-1,lzdiv*lydiv)
        icy0=icy0/lzdiv          +3
        icx0=(ii-1)/(lzdiv*lydiv)+3
         icx=icx0+ixst-3 ! local -> global
         icy=icy0+iyst-3 ! local -> global
         icz=icz0+izst-3 ! local -> global
       icx00=icx0-2+nbound_xm ! local -> local FMM
       icy00=icy0-2+nbound_ym ! local -> local FMM
       icz00=icz0-2+nbound_zm ! local -> local FMM
       wm_local0(:,icz00,icy00,icx00)=0d0
      do iam=0,nomp-1
        wm_local0(:,icz00,icy00,icx00)=
     &                   wm_local0(:,icz00,icy00,icx00)
     &                 +wwm_local0(:,icz00,icy00,icx00,iam)
      enddo ! iam
      enddo ! ii
!$omp end do
!$omp end parallel
      return
      end
c---------------------------------------------------------------------
      subroutine merge_fmm(nl1, nl2)
c---------------------------------------------------------------------
      use comm_base
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use ompvar
      use unitcell
      implicit none
      integer(4) :: nl1, nl2
      integer(4) :: lclzpre, lclypre, lclxpre
      integer(4) :: lclzpost, lclypost, lclxpost
      complex(8),pointer :: wm_local_pre(:,:,:,:) 
      complex(8),pointer :: wm_local_post(:,:,:,:) 
      type(fmm_data_t),pointer :: d1, d2

!NOTE: nl2 -> level n
!      nl1 -> level n-1

      d1 => fmm_data(nl1)
      d2 => fmm_data(nl2)

      if(nl1.lt.lgflg) then
        lclzpre =d1%lclz; lclypre =d1%lcly; lclxpre =d1%lclx
        wm_local_pre => d1%wm_local
        lclzpost=d2%lclz; lclypost=d2%lcly; lclxpost=d2%lclx
        wm_local_post => d2%wm_local
      else if(nl1.eq.lgflg) then
        lclzpre =d1%lclz; lclypre =d1%lcly; lclxpre =d1%lclx
        wm_local_pre => d1%wm_local
        lclzpost=d2%nscell; lclypost=d2%nscell; lclxpost=d2%nscell
        wm_local_post => d2%wm_global
      else
        lclzpre =d1%nscell; lclypre =d1%nscell; lclxpre =d1%nscell
        wm_local_pre => d1%wm_global
        lclzpost=d2%nscell; lclypost=d2%nscell; lclxpost=d2%nscell
        wm_local_post => d2%wm_global
      endif

      call M2M(nl1,nl2,(nmax+1)*(nmax+1),lclzpre, lclypre, lclxpre,
     $         lclzpost,lclypost,lclxpost,
     $         wm_local_pre,wm_local_post,nomp)

      return
      end
c---------------------------------------------------------------------
      subroutine M2M(nl1,nl2,mylm,lclzpre, lclypre, lclxpre,
     $                       lclzpost,lclypost,lclxpost,
     $                       wm_local_preM2M,wm_local_postM2M,nomp)
c---------------------------------------------------------------------
      use comm_base
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use unitcell
      implicit none
      complex(8) wm_local_preM2M(mylm, lclzpre, lclypre, lclxpre)
      complex(8) wm_local_postM2M(mylm,lclzpost,lclypost,lclxpost)
      complex(8) postM2M_omp(mylm,lclzpost,lclypost,lclxpost,0:nomp-1)
      integer(4) :: nl1, nl2
      integer(4) :: ixp, iyp, izp
      integer(4) :: ii,j, m1, m2, m3 
      integer(4) :: icx0,icy0,icz0,icx00,icy00,icz00
      integer(4) :: mcell_size
      integer(4) :: nbound_xm00,nbound_ym00,nbound_zm00
      integer(4) :: nsczdiv,nscydiv,nscxdiv
      integer(4) :: lclzpost,lclypost,lclxpost,lclzpre ,lclypre ,lclxpre
      integer(4) :: mylm,izst2,iyst2,ixst2,izst1,iyst1,ixst1
      integer(4) :: nbound_zm, nbound_ym, nbound_xm
      integer(4) :: iam=0,nomp
      type(fmm_data_t),pointer :: d1, d2
!$    include 'omp_lib.h'

!NOTE: nl2 -> level n
!      nl1 -> level n-1
!
      d1 => fmm_data(nl1)
      d2 => fmm_data(nl2)

!### parameters for level n ###
      nbound_xm00 = d2%nbound_xm
      nbound_ym00 = d2%nbound_ym
      nbound_zm00 = d2%nbound_zm
      nsczdiv = d2%nsczdiv
      nscydiv = d2%nscydiv
      nscxdiv = d2%nscxdiv
      izst2 = (izmin(myrank)-1) / d2%mcell_size+1
      iyst2 = (iymin(myrank)-1) / d2%mcell_size+1
      ixst2 = (ixmin(myrank)-1) / d2%mcell_size+1

!### parameters for level n-1 ###
      nbound_xm = d1%nbound_xm
      nbound_ym = d1%nbound_ym
      nbound_zm = d1%nbound_zm
      izst1 = (izmin(myrank)-1) / d1%mcell_size+1
      iyst1 = (iymin(myrank)-1) / d1%mcell_size+1
      ixst1 = (ixmin(myrank)-1) / d1%mcell_size+1

c process >= number of supercell
      mcell_size = d1%mcell_size
      if(ncell / mcell_size / npz .le.1) nbound_zm=4
      if(ncell / mcell_size / npy .le.1) nbound_ym=4
      if(ncell / mcell_size / npx .le.1) nbound_xm=4

      iam = 0

!$omp parallel default(shared)
!$omp& private(iam)
!$omp& private(ii,m1,m2,m3,j)
!$omp& private(icx0,icy0,icz0,icx00,icy00,icz00)
!$omp& private(ixp,iyp,izp)
!$omp& shared(m_wk,nmerge_fmm,shmm)
!$omp& shared(nsczdiv,nscydiv,nscxdiv)
!$omp& shared(nbound_xm  ,nbound_ym  ,nbound_zm  )
!$omp& shared(nbound_xm00,nbound_ym00,nbound_zm00)
!$omp& shared(wm_local_postM2M)
!$omp& shared(postM2M_omp)
!$    iam=omp_get_thread_num()
      do ii=1,nscxdiv*nscydiv*nsczdiv
        icz0=mod(ii-1,nsczdiv)       +1
        icy0=mod(ii-1,nsczdiv*nscydiv)
        icy0=icy0/nsczdiv            +1 
        icx0=(ii-1)/(nsczdiv*nscydiv)+1

            if(nl1.lt.lgflg)then
              icx00=icx0+nbound_xm00 ! local (nl) -> local FMM (nl) 
              icy00=icy0+nbound_ym00 ! local (nl) -> local FMM (nl) 
              icz00=icz0+nbound_zm00 ! local (nl) -> local FMM (nl) 
              izp = icz0*2-1+nbound_zm ! local (nl) -> local FMM (nl-1)
              iyp = icy0*2-1+nbound_ym ! local (nl) -> local FMM (nl-1)
              ixp = icx0*2-1+nbound_xm ! local (nl) -> local FMM (nl-1)
            else if(nl1.eq.lgflg)then
              icx00=icx0+ixst2-1 ! local (nl) -> global (nl) 
              icy00=icy0+iyst2-1 ! local (nl) -> global (nl) 
              icz00=icz0+izst2-1 ! local (nl) -> global (nl) 
              izp = icz0*2-1+nbound_zm ! local (nl) -> local FMM (nl-1)
              iyp = icy0*2-1+nbound_ym ! local (nl) -> local FMM (nl-1)
              ixp = icx0*2-1+nbound_xm ! local (nl) -> local FMM (nl-1)
            else if(nl1.gt.lgflg)then
              icx00=icx0+ixst2-1 ! local (nl) -> global (nl) 
              icy00=icy0+iyst2-1 ! local (nl) -> global (nl) 
              icz00=icz0+izst2-1 ! local (nl) -> global (nl) 
              izp = icz00*2-1 ! local (nl) -> global (nl-1)
              iyp = icy00*2-1 ! local (nl) -> global (nl-1)
              ixp = icx00*2-1 ! local (nl) -> global (nl-1)
            endif
            postM2M_omp(:,icz00,icy00,icx00,iam)=0d0
!$omp do
            do j = 1,nmerge_fmm
              m1 = m_wk(1,j)
              m2 = m_wk(2,j)
              m3 = m_wk(3,j)
              postM2M_omp(m1,icz00,icy00,icx00,iam)
     $       =postM2M_omp(m1,icz00,icy00,icx00,iam)
     $     +wm_local_preM2M(m3,izp  ,iyp  ,ixp  )*shmm(m2,m1,0,0,0,nl1)
     $     +wm_local_preM2M(m3,izp  ,iyp  ,ixp+1)*shmm(m2,m1,1,0,0,nl1)
     $     +wm_local_preM2M(m3,izp  ,iyp+1,ixp  )*shmm(m2,m1,0,1,0,nl1)
     $     +wm_local_preM2M(m3,izp  ,iyp+1,ixp+1)*shmm(m2,m1,1,1,0,nl1)
     $     +wm_local_preM2M(m3,izp+1,iyp  ,ixp  )*shmm(m2,m1,0,0,1,nl1)
     $     +wm_local_preM2M(m3,izp+1,iyp  ,ixp+1)*shmm(m2,m1,1,0,1,nl1)
     $     +wm_local_preM2M(m3,izp+1,iyp+1,ixp  )*shmm(m2,m1,0,1,1,nl1)
     $     +wm_local_preM2M(m3,izp+1,iyp+1,ixp+1)*shmm(m2,m1,1,1,1,nl1)
            enddo
!$omp end do
      enddo ! ii
!$omp do
      do ii=1,nscxdiv*nscydiv*nsczdiv
        icz0=mod(ii-1,nsczdiv)       +1
        icy0=mod(ii-1,nsczdiv*nscydiv)
        icy0=icy0/nsczdiv            +1 
        icx0=(ii-1)/(nsczdiv*nscydiv)+1
        if(nl1.lt.lgflg)then
          icx00=icx0+nbound_xm00 ! local (nl) -> local FMM (nl) 
          icy00=icy0+nbound_ym00 ! local (nl) -> local FMM (nl) 
          icz00=icz0+nbound_zm00 ! local (nl) -> local FMM (nl) 
            izp=icz0*2-1+nbound_zm ! local (nl) -> local FMM (nl-1)
            iyp=icy0*2-1+nbound_ym ! local (nl) -> local FMM (nl-1)
            ixp=icx0*2-1+nbound_xm ! local (nl) -> local FMM (nl-1)
        else if(nl1.eq.lgflg)then
          icx00=icx0+ixst2-1 ! local (nl) -> global (nl) 
          icy00=icy0+iyst2-1 ! local (nl) -> global (nl) 
          icz00=icz0+izst2-1 ! local (nl) -> global (nl) 
            izp=icz0*2-1+nbound_zm ! local (nl) -> local FMM (nl-1)
            iyp=icy0*2-1+nbound_ym ! local (nl) -> local FMM (nl-1)
            ixp=icx0*2-1+nbound_xm ! local (nl) -> local FMM (nl-1)
        else if(nl1.gt.lgflg)then
          icx00=icx0+ixst2-1 ! local (nl) -> global (nl) 
          icy00=icy0+iyst2-1 ! local (nl) -> global (nl) 
          icz00=icz0+izst2-1 ! local (nl) -> global (nl) 
            izp=icz00*2-1 ! local (nl) -> global (nl-1)
            iyp=icy00*2-1 ! local (nl) -> global (nl-1)
            ixp=icx00*2-1 ! local (nl) -> global (nl-1)
        endif
        wm_local_postM2M(:,icz00,icy00,icx00)=0d0
        do iam=0,nomp-1
          wm_local_postM2M(:,icz00,icy00,icx00)=
     &              wm_local_postM2M(:,icz00,icy00,icx00)
     &                  +postM2M_omp(:,icz00,icy00,icx00,iam)
        enddo
      enddo
!$omp end do
!$omp end parallel

      return
      end
c---------------------------------------------------------------------
      subroutine energy_fmm()
c---------------------------------------------------------------------
      use comm_base
      use trj_org
      use trj_mpi
      use md_coulomb
      use md_fmm
      use md_forces
      use md_segment
      use md_periodic
      use md_monitors
      use md_fmm_domdiv_flg
      use md_const
      use mpivar
      use ompvar
      use param
      use unitcell
#include "timing.h90"
      implicit none
      integer(4) :: izst0,iyst0,ixst0
      integer(4) :: izsts,iysts,ixsts,izens,iyens,ixens
      integer(4) :: ixyzsize,iysize,izsize
      real(8) :: rnnnc
      real(8) :: fxfmm, fyfmm, fzfmm
      integer(4) :: nl
      integer(4) :: i,ii,i0,j0,k0,ipar
      integer(4) :: icx, icy, icz
      integer(4) :: icx0, icy0, icz0
      integer(4) :: icx00, icy00, icz00
      integer(4) :: iam=0
      integer(4) :: m1,m3
      real(8) :: rad, the, csthe, phi, f0p,f0f,algndr
      real(8) :: pot1, qta, xta, yta, zta, g, g0
      real(8) :: xcnt,ycnt,zcnt
      complex(8) :: zphi
      complex(8) :: fi1,fi2,fi3,ilx,ily,ilz,opx,opy,opz
      real*8 f1,f2
      integer(4) :: mcell_size
      complex(8),pointer :: wl_local0(:,:,:,:)
      type(fmm_data_t),pointer :: d1, d2
!$    include 'omp_lib.h'

      wl_local0 => fmm_data(0)%wl_local

      m3=(nmax+1)*(nmax+1) 
      pot1 = 0.0d0

      do nl = 0, nlevel
        fmm_data(nl)%wl_local = 0d0
      end do

******** multipole to local translation
**** contribution of rank 0
      do j0=0,nmax
        do k0=-j0,j0
          m1 = j0*j0+j0+1+k0
          winput(m1) 
     &      = fmm_data(nlevel)%wm_global(m1,1,1,1)/cellx**j0
        enddo
      enddo

      call fmmewald(nmax,winput,woutput,wewald,shew)

      do j0=0,nmax
        do k0=-j0,j0
          m1 = j0*j0+j0+1+k0
          fmm_data(nlevel)%wl_local(m1,1,1,1)
     &      = woutput(m1)/cellx**(j0+1)
        enddo
      enddo

      TIME_START(TM_M2L_L2L)
      DO nl = nlevel, 0, -1

        d1 => fmm_data(nl)
        if (nl > 0) then
          d2 => fmm_data(nl-1)
        else
          d2 => fmm_data(nl)
        end if

!
! OpenMP share variables, depending on nl
!
        mcell_size = d1%mcell_size !! MODIF
        izst0 =(izmin(myrank)-1)/mcell_size+1
        iyst0 =(iymin(myrank)-1)/mcell_size+1
        ixst0 =(ixmin(myrank)-1)/mcell_size+1
c
c     interaction energy by cell mulipole method
c
        if(nl.eq.0) then
          izsts = izmin(myrank)
          iysts = iymin(myrank)
          ixsts = ixmin(myrank)
          izens = izmax(myrank)
          iyens = iymax(myrank)
          ixens = ixmax(myrank)
        else 
          izsts = izmin(myrank)/mcell_size
          iysts = iymin(myrank)/mcell_size
          ixsts = ixmin(myrank)/mcell_size
          if(mod(izmin(myrank),mcell_size).ne.0) izsts = izsts +1
          if(mod(iymin(myrank),mcell_size).ne.0) iysts = iysts +1
          if(mod(ixmin(myrank),mcell_size).ne.0) ixsts = ixsts +1
          izens = izmax(myrank)/mcell_size
          iyens = iymax(myrank)/mcell_size
          ixens = ixmax(myrank)/mcell_size
          if(mod(izmax(myrank),mcell_size).ne.0) izens = izens +1
          if(mod(iymax(myrank),mcell_size).ne.0) iyens = iyens +1
          if(mod(ixmax(myrank),mcell_size).ne.0) ixens = ixens +1
        endif
        ixyzsize= (ixens-ixsts+1)*(iyens-iysts+1)*(izens-izsts+1)
        iysize  = (iyens-iysts+1)
        izsize  = (izens-izsts+1)

        if(nl-1.lt.lgflg)then
          call energy_fmm_2(nl, m3,
     &           d1%wm_local, d1%wl_local, d2%wl_local,
     &           d1%nbound_zm, d1%nbound_ym, d1%nbound_xm,
     &           d1%nbound_zp, d1%nbound_yp, d1%nbound_xp,
     &           d1%nsczdiv, d1%nscydiv, d1%nscxdiv,
     &           d2%nsczdiv, d2%nscydiv, d2%nscxdiv,
     &           d1%nscell, ixyzsize, izsize, iysize,
     &           izsts,iysts,ixsts,izst0,iyst0,ixst0,
     &           d1%wwl_local,nomp)
        else if(nl-1.ge.lgflg)then
          call energy_fmm_3(nl, m3,
     &           d1%wl_local, d1%wm_global, d2%wl_local,
     &           d1%nsczdiv, d1%nscydiv, d1%nscxdiv,
     &           d2%nsczdiv, d2%nscydiv, d2%nscxdiv,
     &           d1%nscell, ixyzsize,izsize,iysize,
     &           izsts,iysts,ixsts,izst0,iyst0,ixst0,
     &           d1%wwl_local,nomp)
        endif

      ENDDO ! nl
      TIME_STOP(TM_M2L_L2L)

      TIME_START(TM_L2P)
      rnnnc  = 1.0d0 / dble(ncell)
!$omp parallel default(none)
!$omp& private(ii,iam,i,i0,ipar)
!$omp& private(j0,f0p,f0f,k0)
!$omp& private(icx0,icy0,icz0,icx,icy,icz)
!$omp& private(icx00,icy00,icz00)
!$omp& private(qta,xta,yta,zta,rad,the,csthe,phi)
!$omp& private(zphi,g,g0,m1,fxfmm,fyfmm,fzfmm,fi1,fi2,fi3)
!$omp& private(f1,f2,ilx,ily,ilz,opx,opy,opz)
!$omp& private(xcnt,ycnt,zcnt)
!$omp& shared(nmax,chgv,paranum,wl_local0)
!$omp& shared(wkxyz,rnnnc,cellx,celly,cellz)
!$omp& shared(pre,w3_f)
!$omp& shared(tag,na_per_cell,m2i)
!$omp& shared(lxdiv,lydiv,lzdiv,ixmin,iymin,izmin,myrank)
!$omp& reduction(+:pot1)
!$    iam = omp_get_thread_num()
CC!$omp do
      do ii=1,lxdiv*lydiv*lzdiv
        icz0=mod(ii-1,lzdiv)     +3   
        icy0=mod(ii-1,lzdiv*lydiv)
        icy0=icy0/lzdiv          +3   
        icx0=(ii-1)/(lzdiv*lydiv)+3
        icx = icx0-3 + ixmin(myrank) ! local -> global
        icy = icy0-3 + iymin(myrank) ! local -> global
        icz = icz0-3 + izmin(myrank) ! local -> global
        icx00=icx0-2                 ! local -> local FMM
        icy00=icy0-2                 ! local -> local FMM
        icz00=icz0-2                 ! local -> local FMM
        xcnt=((icx-0.5d0)*rnnnc-0.5d0)*cellx
        ycnt=((icy-0.5d0)*rnnnc-0.5d0)*celly
        zcnt=((icz-0.5d0)*rnnnc-0.5d0)*cellz
!$omp do
        do i0=tag(icz0,icy0,icx0),
     &        tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
        i=m2i(i0)
        ipar= paranum(i)
        qta = md_QQ_4PiE*chgv(ipar)
        xta = wkxyz(1,i0) - xcnt 
        yta = wkxyz(2,i0) - ycnt 
        zta = wkxyz(3,i0) - zcnt 
        call cart2angle(xta,yta,zta,rad,the,csthe,phi)
        fxfmm=0.d0
        fyfmm=0.d0
        fzfmm=0.d0
**** calculate multipole interaction
        do j0=0,nmax
        f0p=rad**j0
        f0f=qta*rad**(j0-2)

        do k0=-j0,j0
          m1 = j0*j0+j0+1+k0
          zphi=dcmplx(0.d0,k0*phi)
!!!       g=f0p*pre(j0,iabs(k0))*algndr(j0,iabs(k0),csthe)
          g0=pre(j0,iabs(k0))*algndr(j0,iabs(k0),csthe)
          g=f0p*g0
          pot1=pot1+0.5d0*qta
     &        *dble(wl_local0(m1,icz00,icy00,icx00)*g*exp(zphi))
        if(j0==0) cycle

******** calculate multipole force
        if(iabs(k0-1).le.j0) then
          zphi=dcmplx(0.d0,(k0-1)*phi)
          fi1 =-pre(j0,iabs(k0-1))*algndr(j0,iabs(k0-1),csthe)*exp(zphi)
          fi1 =fi1*(-1)**((-1+iabs(k0-1)-iabs(k0))/2)
        else
          fi1 =dcmplx(0.d0,0.d0)
        endif
        zphi=dcmplx(0.d0,k0*phi)
!!!     fi2 =pre(j0,iabs(k0))     *algndr(j0,iabs(k0  ),csthe)*exp(zphi)
        fi2 =g0*exp(zphi)  !! recycle of g0 
        if(iabs(k0+1).le.j0) then
          zphi=dcmplx(0.d0,(k0+1)*phi)
          fi3 =-pre(j0,iabs(k0+1))*algndr(j0,iabs(k0+1),csthe)*exp(zphi)
          fi3 =fi3*(-1)**((+1+iabs(k0+1)-iabs(k0))/2)
        else
          fi3 =dcmplx(0.d0,0.d0)
        endif
        f1 =dsqrt(dble(j0*(j0+1)-k0*(k0+1)))
        f2 =dsqrt(dble(j0*(j0+1)-k0*(k0-1)))
        ilx=dcmplx(0.d0,1.d0)*0.5d0*(f1*fi3+f2*fi1)
        ily=                  0.5d0*(f1*fi3-f2*fi1)
        ilz=dcmplx(0.d0,1.d0)*k0*fi2
        opx=yta*ilz-zta*ily
        opy=zta*ilx-xta*ilz
        opz=xta*ily-yta*ilx
        fxfmm=fxfmm
     $    +f0f*dble(wl_local0(m1,icz00,icy00,icx00)*(-j0*fi2*xta+opx))
        fyfmm=fyfmm
     $    +f0f*dble(wl_local0(m1,icz00,icy00,icx00)*(-j0*fi2*yta+opy))
        fzfmm=fzfmm
     $    +f0f*dble(wl_local0(m1,icz00,icy00,icx00)*(-j0*fi2*zta+opz))

        enddo ! k0
        enddo ! j0
        w3_f(1,i0,iam)=w3_f(1,i0,iam)+fxfmm
        w3_f(2,i0,iam)=w3_f(2,i0,iam)+fyfmm
        w3_f(3,i0,iam)=w3_f(3,i0,iam)+fzfmm
        enddo ! i0
!$omp end do
      enddo ! ii
CC!$omp end do
!$omp end parallel

      wk_p_energy = wk_p_energy + pot1

      TIME_STOP(TM_L2P)

      return
      end
c---------------------------------------------------------------------
      subroutine energy_fmm_2(nl,m3,wm_localx,wl_localx,
     &                        wl_localxx,
     &                        nbound_zm, nbound_ym, nbound_xm,
     &                        nbound_zp, nbound_yp, nbound_xp,
     &                        nsczdiv,nscydiv,nscxdiv,
     &                        nsczdiv_2,nscydiv_2,nscxdiv_2,nscell,
     &                        ixyzsize,izsize,iysize,
     &                        izsts,iysts,ixsts,izst0,iyst0,ixst0,
     &                        wwl_localx,nomp)
c---------------------------------------------------------------------
      use comm_base
      use trj_org
      use trj_mpi
      use md_fmm
      use md_forces
      use md_segment
      use md_periodic
      use md_monitors
      use md_fmm_domdiv_flg
      use md_const
      use mpivar
      use param
      use unitcell
      implicit none
      integer(4) :: izst,iyst,ixst,izen,iyen,ixen
      integer(4) :: izst0,iyst0,ixst0
      integer(4) :: izsts,iysts,ixsts
      integer(4) :: ixyzsize,iysize,izsize
      integer(4) :: iall, nl
      integer(4) :: icx0, icy0, icz0
      integer(4) :: icx00, icy00, icz00
      integer(4) :: icx1, icy1, icz1
      integer(4) :: icx, icy, icz 
      integer(4) :: j,nomp
      integer(4) :: k, iam
      integer(4) :: jxp, jyp, jzp, ic, jc, kc
      integer(4) :: l, m, m1, m2
      integer(4) :: m3,nsczdiv,nscydiv,nscxdiv
      integer(4) :: nsczdiv_2,nscydiv_2,nscxdiv_2
      integer(4) :: nbound_zm, nbound_ym, nbound_xm
      integer(4) :: nbound_zp, nbound_yp, nbound_xp
      complex(8) :: wm_localx(m3,
     &              nbound_zm+nsczdiv+nbound_zp,
     &              nbound_ym+nscydiv+nbound_yp,
     &              nbound_xm+nscxdiv+nbound_xp)
      complex(8) :: wl_localxx(m3,nsczdiv_2,nscydiv_2,nscxdiv_2) 
      complex(8) :: wl_localx(m3,nsczdiv,nscydiv,nscxdiv) 
      complex(8) :: wwl_localx(m3,nsczdiv,nscydiv,nscxdiv,0:nomp-1) 
      integer(4) :: nscell
      integer(4) :: icxg0,icyg0,iczg0,load,ieo
!$    include 'omp_lib.h'

      wwl_localx=0.0d0 

      iam=0 

      if(nl.eq.nlevel) goto 1000

! global address of starting cell.
      icxg0 = (nscell * ipx) / npx + 1
      icyg0 = (nscell * ipy) / npy + 1
      iczg0 = (nscell * ipz) / npz + 1

!$omp parallel default(none)
!$omp& private(iam)
!$omp& private(load)
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1)
!$omp& private(ic,jc,kc,ieo)
!$omp& private(m1,m2)
!$omp& shared(nl,shml)
!$omp& shared(wwl_localx,wm_localx)
!$omp& shared(nload,lddir)
!$omp& shared(nscell)
!$omp& shared(icxg0,icyg0,iczg0)
!$omp& shared(nscxdiv,nscydiv,nsczdiv)
!$omp& shared(nbound_xm,nbound_ym,nbound_zm)
!$omp& shared(nmax)
!$omp& shared(nchunk)
!$    iam=omp_get_thread_num()
!$omp do schedule(static,nchunk)
      DO load = 1, nload
        ic = lddir(1,load)
        jc = lddir(2,load)
        kc = lddir(3,load)

        do icx0 = 1, nscxdiv
          icx1 = icx0 + ic

          if( icxg0 + icx1 - 1 >   3*nscell .or.
     &        icxg0 + icx1 - 1 <= -2*nscell        ) CYCLE

          ieo = mod(icxg0 + icx0, 2)        ! =0 when current cell is odd.
          if(ic < 0) ieo = mod(ieo + 1, 2)  ! reverse when negative direction.
          if(iabs(ic) + ieo > 5) CYCLE      ! odd-even truncation.

          icx1 = icx1 + nbound_xm

        do icy0 = 1, nscydiv
          icy1 = icy0 + jc

          if( icyg0 + icy1 - 1 >   3*nscell .or.
     &        icyg0 + icy1 - 1 <= -2*nscell        ) CYCLE

          ieo = mod(icyg0 + icy0, 2)        ! =0 when current cell is odd.
          if(jc < 0) ieo = mod(ieo + 1, 2)  ! reverse when negative direction.
          if(iabs(jc) + ieo > 5) CYCLE      ! odd-even truncation.

          icy1 = icy1 + nbound_ym

        do icz0 = 1, nsczdiv
          icz1 = icz0 + kc

          if( iczg0 + icz1 - 1 >   3*nscell .or.
     &        iczg0 + icz1 - 1 <= -2*nscell        ) CYCLE

          ieo = mod(iczg0 + icz0, 2)        ! =0 when current cell is odd.
          if(kc < 0) ieo = mod(ieo + 1, 2)  ! reverse when negative direction.
          if(iabs(kc) + ieo > 5) CYCLE      ! odd-even truncation.

          icz1 = icz1 + nbound_zm

**** multipole to local translation
          do m1=1,(nmax+1)*(nmax+1)
            do m2=1,(nmax+1)*(nmax+1)
              wwl_localx(m1,icz0,icy0,icx0,iam)
     $      = wwl_localx(m1,icz0,icy0,icx0,iam)
     $      + wm_localx(m2,icz1,icy1,icx1)*shml(m2,m1,kc,jc,ic,nl)
            enddo
          enddo

        end do ! icz0

        end do ! icy0

        end do ! icx0

      ENDDO ! load
!$omp end do
!$omp end parallel

1000  continue

!$omp parallel default(none)
!$omp& private(iall,icx0,icy0,icz0)
!$omp& private(icx,icy,icz)
!$omp& private(m1,iam)
!$omp& shared(nmax,nomp)
!$omp& shared(izst0,iyst0,ixst0)
!$omp& shared(ixyzsize,iysize,izsize)
!$omp& shared(izsts,iysts,ixsts)
!$omp& shared(wl_localx,wwl_localx,nl)
CC!$omp do
      do iall=0,ixyzsize-1
        icz0 = mod(iall,izsize)+izsts
        icy0 = iall/izsize
        icx0 = icy0/iysize+ixsts
        icy0 = mod(icy0,iysize)+iysts
        icx  =icx0-ixst0+1 ! global -> local FMM M2L L
        icy  =icy0-iyst0+1 ! global -> local FMM M2L L
        icz  =icz0-izst0+1 ! global -> local FMM M2L L
!$omp do
        do m1=1,(nmax+1)*(nmax+1) 
          do iam=0,nomp-1
            wl_localx(m1,icz,icy,icx)=wl_localx(m1,icz,icy,icx)
     &        + wwl_localx(m1,icz,icy,icx,iam) 
          enddo
        enddo 
!$omp end do
      enddo 
CC!$omp end do
!$omp end parallel

      IF (nl /= 0) THEN
!$omp parallel default(none)
!$omp& private(iall,icz0,icy0,icx0)
!$omp& private(izst,iyst,ixst,izen,iyen,ixen,icx00,icy00,icz00)
!$omp& private(jxp,jyp,jzp,ic,jc,kc,icx,icy,icz)
!$omp& private(j,k,l,m,m1,m2)
!$omp& shared(ixyzsize,izsize,iysize,nmax)
!$omp& shared(izsts,iysts,ixsts,izst0,iyst0,ixst0)
!$omp& shared(izmin,iymin,ixmin,izmax,iymax,ixmax,myrank,nl)
!$omp& shared(wl_localxx,wl_localx,shll)
!$omp do
      do iall=0,ixyzsize-1
        icz0 = mod(iall,izsize)+izsts
        icy0 = iall/izsize
        icx0 = icy0/iysize+ixsts
        icy0 = mod(icy0,iysize)+iysts
        izst=(izmin(myrank)-1)/2**(nl-1)+1
        iyst=(iymin(myrank)-1)/2**(nl-1)+1
        ixst=(ixmin(myrank)-1)/2**(nl-1)+1
        izen=(izmax(myrank)-1)/2**(nl-1)+1
        iyen=(iymax(myrank)-1)/2**(nl-1)+1
        ixen=(ixmax(myrank)-1)/2**(nl-1)+1
        icx00=icx0-ixst0+1 ! global -> local FMM L2L
        icy00=icy0-iyst0+1 ! global -> local FMM L2L
        icz00=icz0-izst0+1 ! global -> local FMM L2L

        do jxp = 0, 1
        do jyp = 0, 1
        do jzp = 0, 1
          ic = icx0*2 - 1 + jxp
          jc = icy0*2 - 1 + jyp
          kc = icz0*2 - 1 + jzp

          if(ic.ge.ixst.and.ic.le.ixen.and.
     &       jc.ge.iyst.and.jc.le.iyen.and.
     &       kc.ge.izst.and.kc.le.izen
     &      )then
            icx=ic-ixst+1 ! global (nl-1) -> local FMM L2L (nl-1)
            icy=jc-iyst+1 ! global (nl-1) -> local FMM L2L (nl-1)
            icz=kc-izst+1 ! global (nl-1) -> local FMM L2L (nl-1)
**** local to local translation
            do j=0,nmax ; do k=-j,j ; do l=j,nmax ; do m=-l,l
              if(l-j.lt.iabs(m-k)) cycle
              m1 = j*j+j+1+k
              m2 = l*l+l+1+m
              wl_localxx(m1,icz,icy,icx)=wl_localxx(m1,icz,icy,icx)
     &         +wl_localx(m2,icz00,icy00,icx00)
     &         *shll(m2,m1,jxp,jyp,jzp,nl)
            enddo ; enddo ; enddo ; enddo
          endif

        enddo ! jzp = 0, 1
        enddo ! jyp = 0, 1
        enddo ! jxp = 0, 1
      enddo ! iall=0,ixyzsize-1
!$omp end do
!$omp end parallel
      ELSE  ! nl == 0
!$omp parallel default(none)
!$omp& private(iall,icz0,icy0,icx0,icx,icy,icz,m1)
!$omp& shared(ixyzsize,izsize,iysize)
!$omp& shared(izsts,iysts,ixsts,izst0,iyst0,ixst0)
!$omp& shared(nmax,wl_localxx,wl_localx)
!$omp do
      do iall=0,ixyzsize-1
        icz0 = mod(iall,izsize)+izsts
        icy0 = iall/izsize
        icx0 = icy0/iysize+ixsts
        icy0 = mod(icy0,iysize)+iysts
        icx  =icx0-ixst0+1 ! global -> local FMM M2L L
        icy  =icy0-iyst0+1 ! global -> local FMM M2L L
        icz  =icz0-izst0+1 ! global -> local FMM M2L L
        do m1=1,(nmax+1)*(nmax+1)
          wl_localxx(m1,icz,icy,icx)=wl_localx(m1,icz,icy,icx)
        enddo ! m1
      enddo ! iall
!$omp end do
!$omp end parallel
      ENDIF ! if (nl /= 0)
  
      return
      end 
c---------------------------------------------------------------------
      subroutine energy_fmm_3(nl,m3,wl_localx,wm_globalx,
     &                        wl_localxx,nsczdiv,nscydiv,nscxdiv,
     &                        nsczdiv_2,nscydiv_2,nscxdiv_2,nscell,
     &                        ixyzsize,izsize,iysize,
     &                        izsts,iysts,ixsts,izst0,iyst0,ixst0,
     &                        wwl_localx,nomp)
c---------------------------------------------------------------------
      use trj_org
      use trj_mpi
      use md_fmm
      use md_fmm_domdiv_flg
      use comm_base
      use md_forces
      use md_segment
      use md_periodic
      use md_monitors
      use md_const
      use mpivar
      use param
      use unitcell
      implicit none
      integer(4) :: izst,iyst,ixst,izen,iyen,ixen
      integer(4) :: izst0,iyst0,ixst0
      integer(4) :: izsts,iysts,ixsts 
      integer(4) :: ixyzsize,iysize,izsize
      integer(4) :: iall, nl
      integer(4) :: icx0, icy0, icz0
      integer(4) :: icx00, icy00, icz00
      integer(4) :: icx1, icy1, icz1
      integer(4) :: icx, icy, icz 
      integer(4) :: j
      integer(4) :: k
      integer(4) :: jxp, jyp, jzp, ic, jc, kc
      integer(4) :: l, m, m1, m2
      integer(4) :: nomp,iam
      integer(4) :: m3,nsczdiv,nscydiv,nscxdiv
      integer(4) :: nsczdiv_2,nscydiv_2,nscxdiv_2
      complex(8) :: wm_globalx(m3,nscell,nscell,nscell)
      complex(8) :: wl_localxx(m3,nsczdiv_2,nscydiv_2,nscxdiv_2) 
      complex(8) :: wl_localx(m3,nsczdiv,nscydiv,nscxdiv) 
      complex(8) :: wwl_localx(m3,nsczdiv,nscydiv,nscxdiv,0:nomp-1) 
      integer(4) :: nscell
      integer(4) :: icxg0,icyg0,iczg0,load,ieo
!$    include 'omp_lib.h'

      wwl_localx=0.0d0 
      iam=0 

! global address of starting cell.
      icxg0 = (nscell * ipx) / npx + 1
      icyg0 = (nscell * ipy) / npy + 1
      iczg0 = (nscell * ipz) / npz + 1

      if(nl.eq.nlevel) goto 1000

!$omp parallel default(none)
!$omp& private(iam)
!$omp& private(load)
!$omp& private(icx0,icy0,icz0,icx1,icy1,icz1)
!$omp& private(ic,jc,kc,ieo)
!$omp& private(m1,m2)
!$omp& shared(nl,shml)
!$omp& shared(wwl_localx,wm_globalx)
!$omp& shared(nload,lddir)
!$omp& shared(nscell)
!$omp& shared(icxg0,icyg0,iczg0)
!$omp& shared(nscxdiv,nscydiv,nsczdiv)
!$omp& shared(nmax)
!$omp& shared(nchunk)
!$     iam=omp_get_thread_num() 
!$omp do schedule(static,nchunk)
      DO load = 1, nload
        ic = lddir(1,load)
        jc = lddir(2,load)
        kc = lddir(3,load)

        do icx0 = 1, nscxdiv
          icx1 = icxg0 + icx0 - 1 + ic

          if( icx1 >  3*nscell .or. icx1 <= -2*nscell ) CYCLE

          ieo = mod(icxg0 + icx0, 2)        ! =0 when current cell is odd.
          if(ic < 0) ieo = mod(ieo + 1, 2)  ! reverse when negative direction.
          if(iabs(ic) + ieo > 5) CYCLE      ! odd-even truncation.

          icx1 = mod( icx1 + 6*nscell - 1, nscell) + 1

        do icy0 = 1, nscydiv
          icy1 = icyg0 + icy0 - 1 + jc

          if( icy1 >  3*nscell .or. icy1 <= -2*nscell ) CYCLE

          ieo = mod(icyg0 + icy0, 2)        ! =0 when current cell is odd.
          if(jc < 0) ieo = mod(ieo + 1, 2)  ! reverse when negative direction.
          if(iabs(jc) + ieo > 5) CYCLE      ! odd-even truncation.

          icy1 = mod( icy1 + 6*nscell - 1, nscell) + 1

        do icz0 = 1, nsczdiv
          icz1 = iczg0 + icz0 - 1 + kc

          if( icz1 >  3*nscell .or. icz1 <= -2*nscell ) CYCLE

          ieo = mod(iczg0 + icz0, 2)        ! =0 when current cell is odd.
          if(kc < 0) ieo = mod(ieo + 1, 2)  ! reverse when negative direction.
          if(iabs(kc) + ieo > 5) CYCLE      ! odd-even truncation.

          icz1 = mod( icz1 + 6*nscell - 1, nscell) + 1

**** multipole to local translation
          do m1=1,(nmax+1)*(nmax+1)
            do m2=1,(nmax+1)*(nmax+1)
              wwl_localx(m1,icz0,icy0,icx0,iam)
     $      = wwl_localx(m1,icz0,icy0,icx0,iam)
     $      + wm_globalx(m2,icz1,icy1,icx1)*shml(m2,m1,kc,jc,ic,nl)
            enddo
          enddo

        end do ! icz0

        end do ! icy0

        end do ! icx0

      enddo ! load
!$omp end do
!$omp end parallel

1000    continue

      DO iall=0,ixyzsize-1

        icz0 = mod(iall,izsize)+izsts
        icy0 = iall/izsize
        icx0 = icy0/iysize+ixsts
        icy0 = mod(icy0,iysize)+iysts

        ixst=(ixmin(myrank)-1)/2**nl+1
        iyst=(iymin(myrank)-1)/2**nl+1
        izst=(izmin(myrank)-1)/2**nl+1
        icx=icx0-ixst+1
        icy=icy0-iyst+1
        icz=icz0-izst+1
!$omp parallel default(none) 
!$omp& private(m1,iam) 
!$omp& shared(nmax,icz,icy,icx,nomp)
!$omp& shared(wl_localx,wwl_localx)
!$omp do 
        do m1=1,(nmax+1)*(nmax+1) 
          do iam=0,nomp-1
            wl_localx(m1,icz,icy,icx)=wl_localx(m1,icz,icy,icx) 
     &           +wwl_localx(m1,icz,icy,icx,iam) 
          enddo 
        enddo 
!$omp end do 
!$omp end parallel 

      if (nl /= 0) then

        ixst=(ixmin(myrank)-1)/2**(nl-1)+1
        iyst=(iymin(myrank)-1)/2**(nl-1)+1
        izst=(izmin(myrank)-1)/2**(nl-1)+1
        ixen=(ixmax(myrank)-1)/2**(nl-1)+1
        iyen=(iymax(myrank)-1)/2**(nl-1)+1
        izen=(izmax(myrank)-1)/2**(nl-1)+1
        icx00=icx0-ixst0+1 ! global -> local FMM L2L
        icy00=icy0-iyst0+1 ! global -> local FMM L2L
        icz00=icz0-izst0+1 ! global -> local FMM L2L

        do jxp = 0, 1
        do jyp = 0, 1
        do jzp = 0, 1
          ic = icx0*2 - 1 + jxp
          jc = icy0*2 - 1 + jyp
          kc = icz0*2 - 1 + jzp
          if(ic.ge.ixst.and.ic.le.ixen.and.
     &       jc.ge.iyst.and.jc.le.iyen.and.
     &       kc.ge.izst.and.kc.le.izen
     &      )then
            icx=ic-ixst+1 ! global (nl-1) -> local FMM L2L (nl-1)
            icy=jc-iyst+1 ! global (nl-1) -> local FMM L2L (nl-1)
            icz=kc-izst+1 ! global (nl-1) -> local FMM L2L (nl-1)
**** local to local translation
            do j=0,nmax ; do k=-j,j ; do l=j,nmax ; do m=-l,l
              if(l-j.lt.iabs(m-k)) cycle
              m1 = j*j+j+1+k
              m2 = l*l+l+1+m
              wl_localxx(m1,icz,icy,icx)=wl_localxx(m1,icz,icy,icx)
     &         +wl_localx(m2,icz00,icy00,icx00)
     &         *shll(m2,m1,jxp,jyp,jzp,nl)
            enddo ; enddo ; enddo ; enddo
          endif

        enddo ! jzp = 0, 1
        enddo ! jyp = 0, 1
        enddo ! jxp = 0, 1
      endif ! if (nl /= 0)

      ENDDO ! iall=0,ixyzsize-1

      return
      end 
c----------------------------------------------------------------------
      subroutine init_fmm()
c---------------------------------------------------------------------
      use trj_mpi ! ya
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use ompvar
      use unitcell
      implicit none
      integer(4) :: i, j, k
      integer(4) :: l, m, m1, m2
      integer(4) :: ix, iy, iz
      integer(4) :: il
      integer(4) :: ncref,ncdir
      integer(4) :: nshell,ishl
      integer(4) :: ithr
      integer(4) :: irx, iry, irz

      integer(4) :: iht
      integer(4) :: ith,nscan
      integer(4) :: ist,ien
      integer(4) :: isign,mst
      integer(4),allocatable,dimension(:) :: ixst,ixen
      integer(4),allocatable,dimension(:) :: iyst,iyen
      integer(4),allocatable,dimension(:) :: izst,izen
      real(8) :: factorial
      real(8) :: xtb,ytb,ztb,rad,the,csthe,phi,rnnnc
      real(8) :: exp, algndr
      complex(8) :: zphi

!*** For any number of threads at energy_fmm.              ************
!... Generation of relative cell address list.                      ...
!
! maximum referencing cell-range for M2L translation..
      ncref = 5
! cell range of direct force calculation.
      ncdir = 2
! relative cell address list.
      allocate( lddir(3, (2*ncref+1)**3 - (2*ncdir+1)**3) )
! openmp chunk size definition. to avoid difference of OMP implementation.
!$    nchunk = ((2*ncref+1)**3 - (2*ncdir+1)**3 - 1)/nomp + 1

!... A multi-level cubic shells ordering with cyclic assignment is  ...
!... employed to balance computational load of M2L translation among...
!... threads is employed. Cyclic assignment causes frequent         ...
!... replacement of j-cell moment.                                  ...
! working array in this subroutine.
      allocate( ixst( (6+2)*ncref ) )
      allocate( ixen( (6+2)*ncref ) )
      allocate( iyst( (6+2)*ncref ) )
      allocate( iyen( (6+2)*ncref ) )
      allocate( izst( (6+2)*ncref ) )
      allocate( izen( (6+2)*ncref ) )

! Definition of multi-layer cubic shell. Ordered from center.
      nshell = 0
      do i = 1, ncref
        nshell = nshell + 1
        izst(nshell) = -i + 1 ; izen(nshell) = i
        iyst(nshell) = -i     ; iyen(nshell) = -i
        ixst(nshell) = -i     ; ixen(nshell) = i - 1
        nshell = nshell + 1
        izst(nshell) = -i     ; izen(nshell) = i - 1
        iyst(nshell) = i      ; iyen(nshell) = i
        ixst(nshell) = -i + 1 ; ixen(nshell) = i
        nshell = nshell + 1
        izst(nshell) = -i     ; izen(nshell) = i - 1
        iyst(nshell) = -i + 1 ; iyen(nshell) = i
        ixst(nshell) = -i     ; ixen(nshell) = -i
        nshell = nshell + 1
        izst(nshell) = -i + 1 ; izen(nshell) = i
        iyst(nshell) = -i     ; iyen(nshell) = i - 1
        ixst(nshell) = i      ; ixen(nshell) = i
        nshell = nshell + 1
        izst(nshell) = -i     ; izen(nshell) = -i
        iyst(nshell) = -i     ; iyen(nshell) = i - 1
        ixst(nshell) = -i + 1 ; ixen(nshell) = i
        nshell = nshell + 1
        izst(nshell) = i      ; izen(nshell) = i
        iyst(nshell) = -i + 1 ; iyen(nshell) = i
        ixst(nshell) = -i     ; ixen(nshell) = i - 1
        nshell = nshell + 1
        izst(nshell) = -i     ; izen(nshell) = -i
        iyst(nshell) = -i     ; iyen(nshell) = -i
        ixst(nshell) = -i     ; ixen(nshell) = -i
        nshell = nshell + 1
        izst(nshell) = i      ; izen(nshell) = i
        iyst(nshell) = i      ; iyen(nshell) = i
        ixst(nshell) = i      ; ixen(nshell) = i
      end do

! Relative cell address list ordered cyclic on multi-layer cubic shell.
! for scan order including turning cyclic. turning might be wrong. use no turning.
!      isign = -1
      mst = -1
      isign = 1
!      mst = 1
      nload = 0

      do ithr = 1, nomp
! randomize thread selection during cyclic scan. hard to be effective.
!        ith = mod(ithr*3-1, nomp) + 1
        ith = ithr
        mst = mst * isign
        if(mst>0) then
          ist = 1
          ien = nshell
        else
          ist = nshell
          ien = 1
        endif
        nscan = 0

      do ishl = ist, ien, mst
        do irx = ixst(ishl), ixen(ishl)
          do iry = iyst(ishl), iyen(ishl)
            do irz = izst(ishl), izen(ishl)
              iht = 0
              if( irx < -2 .OR. irx > 2 ) iht = 1
              if( iry < -2 .OR. iry > 2 ) iht = 1
              if( irz < -2 .OR. irz > 2 ) iht = 1
              if( iht == 1 ) then
                nscan = nscan + 1
                if(mod(nscan - 1, nomp) + 1 == ith) then
                  nload = nload + 1
                  lddir(1,nload) = irx
                  lddir(2,nload) = iry
                  lddir(3,nload) = irz
                endif
              endif
            end do  ! irz
          end do  ! iry
        end do  ! irx

      end do  ! ishl

      end do  ! ithr
!
      deallocate( ixst )
      deallocate( ixen )
      deallocate( iyst )
      deallocate( iyen )
      deallocate( izst )
      deallocate( izen )

!*** End of relative cell address list generation.  *****************
Cfj
c
      mdg=2*nmax   !! moved from fmodules.f

      allocate(fa(0:mdg,-mdg:mdg))
      do j=0,mdg
        do k=-mdg, mdg
          fa(j,k) = (-1)**j/dsqrt(factorial(j-k)*factorial(j+k))
        enddo
      enddo
c
      allocate(pre(0:mdg,0:mdg))
******** prepare the pre-factor
      do j=0,mdg
        do k=0,j
          pre(j,k)=dsqrt(factorial(j-k)/factorial(j+k))
        enddo
      enddo
c
      do j=0,nmax
        do k=-j,j
          do l=0,j
            do m=-l,l
              if ((j-l) < abs(k-m)) cycle
              m1 = j*j+j+1+k
              m2 = l*l+l+1+m
              premm(m1,m2) = fa(l,m)*fa(j-l,k-m)/fa(j,k)
     &                   *dcmplx(0.d0,1.d0)**(iabs(k)-iabs(m)-iabs(k-m))
            enddo
          enddo
        enddo
      enddo
c
      do j=0,nmax
        do k=-j,j
          do l=0,nmax
            do m=-l,l
              if(j+l.lt.iabs(m-k)) cycle
              m1 = j*j+j+1+k
              m2 = l*l+l+1+m
              preml(m1,m2) = fa(j,k)*fa(l,m)/fa(j+l,m-k)*(-1)**l
     &                   *dcmplx(0.d0,1.d0)**(iabs(k-m)-iabs(k)-iabs(m))
     &                   *pre(j+l,iabs(m-k))
            enddo
          enddo
        enddo
      enddo
c
      do j=0,nmax
        do k=-j,j
          do l=j,nmax
            do m=-l,l
              if(l-j.lt.iabs(m-k)) cycle
              m1 = j*j+j+1+k
              m2 = l*l+l+1+m
              prell(m1,m2) = fa(l-j,m-k)*fa(j,k)/fa(l,m)/(-1)**(l+j)
     &               *pre(l-j,iabs(m-k))
     &               *dcmplx(0.d0,1.d0)**(iabs(m)-iabs(k)-iabs(m-k))
            enddo
          enddo
        enddo
      enddo
c
      do il=0,nlevel
        rnnnc = 2**il / dble(ncell)
        do ix=0,1
        do iy=0,1
        do iz=0,1
          xtb = (ix-0.5d0)*rnnnc*cellx
          ytb = (iy-0.5d0)*rnnnc*celly
          ztb = (iz-0.5d0)*rnnnc*cellz
          call cart2angle(xtb,ytb,ztb,rad,the,csthe,phi)
          do j=0,nmax
            do k=-j,j
              m1 = j*j+j+1+k
              do l=0,j
                do m=-l,l
                  if ((j-l) < abs(k-m)) cycle
                  m2 = l*l+l+1+m
                  zphi=dcmplx(0.d0,m*phi)
                  shmm(m2,m1,ix,iy,iz,il) = rad**l*pre(l,iabs(m))
     &                 *algndr(l,iabs(m),csthe)*exp(-zphi)*premm(m1,m2)
                enddo
              enddo
            enddo
          enddo
        enddo
        enddo
        enddo
c
        do ix=-5,5
        do iy=-5,5
        do iz=-5,5
          xtb = ix * rnnnc * cellx
          ytb = iy * rnnnc * celly
          ztb = iz * rnnnc * cellz
          call cart2angle(xtb,ytb,ztb,rad,the,csthe,phi)
          if(rad.ne.0d0)then
          do j=0,nmax
            do k=-j,j
              m1 = j*j+j+1+k
              do l=0,nmax
                do m=-l,l
                  if(j+l.lt.iabs(m-k)) cycle
                  zphi=dcmplx(0.d0,(m-k)*phi)
                  m2 = l*l+l+1+m
C                  shml(m2,m1,ix,iy,iz,il) = 1.0d0/rad**(j+l+1)
                  shml(m2,m1,iz,iy,ix,il) = 1.0d0/rad**(j+l+1)
     &               *algndr(j+l,iabs(m-k),csthe)*exp(zphi)*preml(m1,m2)
                enddo
              enddo
            enddo
          enddo
          endif
        enddo
        enddo
        enddo
c
        do ix=0,1
        do iy=0,1
        do iz=0,1
          xtb = -(ix-0.5d0)*rnnnc*0.5d0*cellx
          ytb = -(iy-0.5d0)*rnnnc*0.5d0*celly
          ztb = -(iz-0.5d0)*rnnnc*0.5d0*cellz
          call cart2angle(xtb,ytb,ztb,rad,the,csthe,phi)
          do j=0,nmax
          do k=-j,j
          do l=j,nmax
            do m=-l,l
              if(l-j.lt.iabs(m-k)) cycle
              m1 = j*j+j+1+k
              m2 = l*l+l+1+m
              zphi=dcmplx(0.d0,(m-k)*phi)
              shll(m2,m1,ix,iy,iz,il) = rad**(l-j)
     &               *algndr(l-j,iabs(m-k),csthe)*exp(zphi)*prell(m1,m2)
            enddo
          enddo
          enddo
          enddo
        enddo
        enddo
        enddo
      enddo
      deallocate(premm)
      deallocate(preml)
      deallocate(prell)
      call init_fmmewald(mdg,nmax,pre)
      call init_merge_fmm()
c
      do j=0,nmax
        do k=-j,j
          m1 = j*j+j+1+k
          do l=0,nmax-j
            do m=-l,l
              if(j+l.lt.iabs(m-k)) cycle
              m2 = l*l+l+1+m
              shew(m2,m1)=fa(l,m)*fa(j,k)/fa(j+l,m-k)
     &           *dcmplx(0.d0,1.d0)**(iabs(k-m)-iabs(k)-iabs(m))*(-1)**l
            enddo
          enddo
        enddo
      enddo
C
      return
      end
c---------------------------------------------------------------------
      subroutine init_merge_fmm()
c---------------------------------------------------------------------
      use md_fmm
      implicit none
      integer(4) :: j,k,m,l, m1, m20, m30,ntmp1

      ntmp1 = 0
      do j=0,nmax
        do k=-j,j
          m1 = j*j+j+1+k
          do l=0,j
            do m=-l,l
              if ((j-l) < abs(k-m)) cycle
              ntmp1 = ntmp1 + 1
            enddo
          enddo
        enddo
      enddo
      allocate(m_wk(3,ntmp1))
      ntmp1 = 0
      do j=0,nmax
        do k=-j,j
          m1 = j*j+j+1+k
          do l=0,j
            m20 = l*l+l+1
            m30 = (j-l)*(j-l)+(j-l)+1+k
            do m=-l,l
              if ((j-l) < abs(k-m)) cycle
              ntmp1 = ntmp1 + 1
              m_wk(1,ntmp1) = m1
              m_wk(2,ntmp1) = m20+m
              m_wk(3,ntmp1) = m30-m
            enddo
          enddo
        enddo
      enddo
      nmerge_fmm = ntmp1

      return
      end
c----------------------------------------------------------------------
      subroutine init_fmmewald(mdg,nmax,pre)
c----------------------------------------------------------------------
      use mod_wk_fmmewald
      implicit none
      integer(4) :: nmax,nhmax,nrmax
      real(8) :: wkappa,nhcut,nrcut,wkappa2
      real(8) :: pi,nhcut2,nrcut2,vol,nh2,nr2
      integer(4) :: mdg
      integer(4) :: i,j,k,m,n,m1
      real(8) :: pre(0:mdg,0:mdg)
      real(8) :: xta,yta,zta,rad,the,csthe,phi
      real(8) :: fvh,lngamma,incmpgamm,algndr
      complex(8) :: fi,zphi

******** set Ewald parameter
      nhmax=20
      nrmax=20
      wkappa=2.d0
      nhcut=nhmax
      nrcut=nrmax
      vol=1.d0
******** set Ewald parameter
      wkappa2=wkappa*wkappa
      nhcut2=nhcut*nhcut
      nrcut2=nrcut*nrcut
******** set constants
      pi=dacos(-1.d0)

      allocate(wk_wl((nmax+1)*(nmax+1)))
      wk_wl = dcmplx(0.d0,0.d0)

      do i=-nhmax,nhmax
      do j=-nhmax,nhmax
      do k=-nhmax,nhmax
        nh2=i*i+j*j+k*k
        if(nh2.gt.nhcut2)goto 1100
        if(nh2.eq.0.d0)goto 1100
        xta=dble(i)
        yta=dble(j)
        zta=dble(k)
        call cart2angle(xta,yta,zta,rad,the,csthe,phi)
        do n=2,nmax,2
          fi=dcmplx(0.d0,1.d0)**n*pi**(n-0.5d0)
          fvh=dexp(-pi**2*nh2/wkappa2-lngamma(n+0.5d0))*rad**(n-2)/vol
          do m=-n,n,2
            m1 = n*n+n+1+m
            zphi=dcmplx(0.d0,m*phi)
            wk_wl(m1) = wk_wl(m1)
     &              +fi*fvh*pre(n,iabs(m))
     &              *algndr(n,iabs(m),csthe)*exp(zphi)
          enddo
        enddo
1100  enddo
      enddo
      enddo

      do i=-nrmax,nrmax
      do j=-nrmax,nrmax
      do k=-nrmax,nrmax
        nr2=i*i+j*j+k*k
        if(nr2.gt.nrcut2)goto 1200
        if(iabs(i).le.2.and.iabs(j).le.2.and.iabs(k).le.2)then
          if(nr2.eq.0.d0)goto 1200
          xta=dble(i)
          yta=dble(j)
          zta=dble(k)
          call cart2angle(xta,yta,zta,rad,the,csthe,phi)
          if(rad.ne.0d0)then
          do n=2,nmax,2
            fvh=incmpgamm(n+0.5d0,wkappa2*rad**2)/rad**(n+1)
            do m=-n,n,2
              m1 = n*n+n+1+m
              zphi=dcmplx(0.d0,m*phi)
              wk_wl(m1) = wk_wl(m1)
     &              -fvh*pre(n,iabs(m))
     &              *algndr(n,iabs(m),csthe)*exp(zphi)
            enddo
          enddo
          endif
        else
          xta=dble(i)
          yta=dble(j)
          zta=dble(k)
          call cart2angle(xta,yta,zta,rad,the,csthe,phi)
          if(rad.ne.0d0)then
          do n=2,nmax,2
            fvh=(1d0-incmpgamm(n+0.5d0,wkappa2*rad**2))/rad**(n+1)
            do m=-n,n,2
              m1 = n*n+n+1+m
              zphi=dcmplx(0.d0,m*phi)
              wk_wl(m1) = wk_wl(m1)
     &              +fvh*pre(n,iabs(m))
     &              *algndr(n,iabs(m),csthe)*exp(zphi)
            enddo
          enddo
          endif
        endif
1200  enddo
      enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine fmmewald(nmax,wm,we,wl,shew)
c----------------------------------------------------------------------
      use mod_wk_fmmewald
      implicit none
      integer(4) :: j,k,m,n,m1,m2,m3
      integer(4) :: nmax
      complex(8) :: wm((nmax+1)*(nmax+1))
      complex(8) :: wl((nmax+1)*(nmax+1))
      complex(8) :: shew((nmax+1)*(nmax+1),(nmax+1)*(nmax+1))
      complex(8) :: we((nmax+1)*(nmax+1))

      we=dcmplx(0.d0,0.d0)
      do k = 1,(nmax+1)*(nmax+1)
          wl(k) = wk_wl(k)
      end do
      do j=0,nmax
        do k=-j,j
          m1 = j*j+j+1+k
          do n=0,nmax-j
            do m=-n,n
              if(j+n.lt.iabs(m-k)) cycle
              m2 = n*n+n+1+m
              m3 = (j+n)*(j+n)+(j+n)+1+(m-k)
              we(m1)=we(m1)+wm(m2)*shew(m2,m1)*wl(m3)
            enddo
          enddo
        enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine calc_system_dipole
c----------------------------------------------------------------------
      use md_fmm
      use md_fmm_domdiv_flg
      use mpivar
      use unitcell
      implicit none
      integer(4) :: j0,k0,m1
      complex(8) :: M1m1    ! M_1^(-1)
      complex(8) :: M10     ! M_1^0
      complex(8) :: M1p1    ! M_1^(+1)
      complex(8), parameter :: iu = (0.0d0,1.0d0)

!### calc M_1^(-1) ###!
      j0=+1
      k0=-1
      m1=j0*j0+j0+1+k0
      M1m1=fmm_data(nlevel)%wm_global(m1,1,1,1)

!### calc M_1^( 0) ###!
      j0=+1
      k0= 0
      m1=j0*j0+j0+1+k0
      M10=fmm_data(nlevel)%wm_global(m1,1,1,1)

!### calc M_1^(+1) ###!
      j0=+1
      k0=+1
      m1=j0*j0+j0+1+k0
      M1p1=fmm_data(nlevel)%wm_global(m1,1,1,1)

!### calc system dipole ###!
      sysdpl(1)=-1d0/sqrt(2d0)*(M1p1+M1m1)
      sysdpl(2)=+1d0/sqrt(2d0)*(M1p1-M1m1)/iu
      sysdpl(3)=M10

      return
      end
c----------------------------------------------------------------------
      subroutine remove_ewaldsurfaceterm
c----------------------------------------------------------------------
      use trj_mpi
      use md_coulomb
      use md_const
      use md_fmm
      use md_fmm_domdiv_flg
      use md_forces
      use md_monitors
      use md_segment
      use param
      use unitcell
      use mpivar
      implicit none
      real(8) :: DD, coef0, coefi
      real(8) :: ewald_st
      integer(4) :: k0,i0,ipar
!$    include 'omp_lib.h'
      integer(4) :: iam=0

!### calc. system-dipole ###
      sysdpl=sysdpl*md_ELEMENTARY_CHARGE
      DD=sysdpl(1)*sysdpl(1)+sysdpl(2)*sysdpl(2) +sysdpl(3)*sysdpl(3)

!### calc. force to be removed ###
      coef0=md_ELEMENTARY_CHARGE/
     &  (3d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)

!$omp parallel default(none)
!$omp& private(iam,k0,i0,ipar,coefi)
!$omp& shared(nselfseg,lsegtop,lseg_natoms)
!$omp& shared(chgv,paranum,m2i,coef0)
!$omp& shared(w3_f,sysdpl)
!$    iam=omp_get_thread_num()
!$omp do
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        ipar=paranum(m2i(i0))
        coefi=chgv(ipar)*coef0
        w3_f(1,i0,iam)=w3_f(1,i0,iam)+coefi*sysdpl(1)
        w3_f(2,i0,iam)=w3_f(2,i0,iam)+coefi*sysdpl(2)
        w3_f(3,i0,iam)=w3_f(3,i0,iam)+coefi*sysdpl(3)
      enddo
      enddo
!$omp end do
!$omp end parallel

      if(myrank==0)then
!### calc. poteintail to be removed ###
      ewald_st=DD/(6d0*cellvol*md_VACUUM_DIELECTRIC_CONSTANT)
      wk_p_energy = wk_p_energy - ewald_st
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine cart2angle(xta,yta,zta,rad,the,csthe,phi)
c----------------------------------------------------------------------
      implicit none
      real*8 xta,yta,zta,rad,the,csthe,phi
      real*8 pi
******** set constants
      pi=dacos(-1.d0)
******** calclate angles
      rad=dsqrt(xta*xta+yta*yta+zta*zta)
      if(rad.eq.0.d0) then
        the=0.d0
        phi=0.d0
        csthe=1.d0
      elseif(zta.ne.0.d0) then
        the=dacos(zta/rad)
        csthe=zta/rad
        if(xta.ne.0.d0) then
          phi=datan(yta/xta)
          if(xta.le.0.d0)then
            if(yta.ge.0.d0) phi=phi+pi
            if(yta.lt.0.d0) phi=phi-pi
          endif
        else
          if(yta.gt.0.d0) phi=0.5d0*pi
          if(yta.lt.0.d0) phi=-0.5d0*pi
        endif
      else
        the=0.5d0*pi
        csthe=0.d0
        if(xta.ne.0.d0) then
          phi=datan(yta/xta)
          if(xta.le.0.d0)then
            if(yta.ge.0.d0) phi=phi+pi
            if(yta.lt.0.d0) phi=phi-pi
          endif
        else
          if(yta.gt.0.d0) phi=0.5d0*pi
          if(yta.lt.0.d0) phi=-0.5d0*pi
        endif
      endif
      return
      end
**********************************************************************
      function algndr(n,m,x)
**********************************************************************
      implicit none
      integer(4)::m,n
      real(8)::x,algndr
      real(8)::apnn2,apnn1,apnn0,p
      integer(4)::i,j,k

      apnn2=1.d0
      if(m.gt.0) then
        p=dsqrt(1.d0-x*x)
        do i=1,m
          apnn2=-apnn2*(2.d0*i-1.d0)*p
        enddo
      endif
      if(n.eq.m) then
        algndr=apnn2
      else
        apnn1=x*(2*m+1)*apnn2
        if(n.eq.m+1) then
          algndr=apnn1
        else
          do k=m+2,n
            p=(x*(2*k-1)*apnn1-(k+m-1)*apnn2)/(k-m)
            apnn2=apnn1
            apnn1=p
          enddo
          algndr=p
        endif
      endif
      return
      end
**********************************************************************
      function factorial(n)
**********************************************************************
      implicit none
      real*8 factorial
      integer n,i
      factorial = 1.d0
      do i=1,n
         factorial = factorial * i
      enddo
      return
      end
**********************************************************************
      function incmpgamm(a,x)
**********************************************************************
      real*8 a,incmpgamm,x
      real*8 ren,kyu
      if(x.lt.a+1.d0)then
        incmpgamm=kyu(a,x)
      else
        incmpgamm=1.d0-ren(a,x)
      endif
      return
      end
**********************************************************************
      function ren(a,x)
**********************************************************************
      integer itmax
      real*8 a,ren,gln,x,eps,fpmin
      parameter (itmax=100,eps=3.d-7,fpmin=1.d-30)
      integer i
      real*8 an,b,c,d,del,h,lngamma
      gln=lngamma(a)
      b=x+1.d0-a
      c=1.d0/fpmin
      d=1.d0/b
      h=d
      do 11 i=1,itmax
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0).lt.eps)goto 1
11    continue
1     ren=exp(-x+a*log(x)-gln)*h
      return
      end
**********************************************************************
      function kyu(a,x)
**********************************************************************
      integer itmax
      real*8 a,kyu,gln,x,eps
      parameter (itmax=100,eps=3.d-7)
      integer n
      real*8 ap,del,sum,lngamma
      gln=lngamma(a)
      if(x.le.0.d0)then
        kyu=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,itmax
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)goto 1
11    continue
1     kyu=sum*exp(-x+a*log(x)-gln)
      return
      end
**********************************************************************
      function lngamma(xx)
**********************************************************************
      implicit none
      real*8 lngamma,xx,ser,tmp

      tmp  = xx + 4.5d0
      tmp  = tmp-(xx - 0.5d0)*dlog(tmp)
      ser = 1.000000000190015d0
     $     + (76.18009172947146d0   / xx)
     $     - (86.50532032941677d0   / (xx + 1.0d0))
     $     + (24.01409824083091d0   / (xx + 2.0d0))
     $     - (1.231739572450155d0   / (xx + 3.0d0))
     $     + (0.1208650973866179d-2 / (xx + 4.0d0))
     $     - (0.5395239384953d-5    / (xx + 5.0d0))

      lngamma=dlog(2.5066282746310005d0 * ser) - tmp
      return
      end
**********************************************************************
