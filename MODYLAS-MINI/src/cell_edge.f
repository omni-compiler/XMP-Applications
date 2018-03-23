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
      subroutine cell_edge
c----------------------------------------------------------------------
      use atommass
      use trj_org
      use unitcell
      use md_periodic
      use md_segment
      use param
      implicit none
      integer(4) :: i,isegstart,isegend,k,ipar
      real(8) :: totm,tmpx,tmpy,tmpz,segrx,segry,segrz
      real(8) :: segsx,segsy,segsz
      real(8) :: r_cellx,r_celly,r_cellz

      r_cellx = 1.0d0/cellx
      r_celly = 1.0d0/celly
      r_cellz = 1.0d0/cellz

!$omp parallel do default(none)
!$omp& private(k,i,isegstart,isegend,ipar)
!$omp& private(tmpx,tmpy,tmpz,totm)
!$omp& private(segsx,segsy,segsz,segrx,segry,segrz)
!$omp& shared(nsegments,segtop,seg_natoms)
!$omp& shared(cellx,celly,cellz)
!$omp& shared(r_cellx, r_celly, r_cellz)
!$omp& shared(paranum,mass,xyz,seg_cx,seg_cy,seg_cz)
      do k = 1, nsegments
        isegstart=segtop(k)+1
        isegend  =isegstart+seg_natoms(k)-1
        tmpx=0.0d0;tmpy=0.0d0;tmpz=0.0d0;totm=0d0
        do i = isegstart, isegend
          ipar=paranum(i)
          totm=totm+mass(ipar)
          tmpx=tmpx+xyz(1,i)*mass(ipar)
          tmpy=tmpy+xyz(2,i)*mass(ipar)
          tmpz=tmpz+xyz(3,i)*mass(ipar)
        enddo ! i
        tmpx=tmpx/totm
        tmpy=tmpy/totm
        tmpz=tmpz/totm
!       ### normarize 
        segsx = r_cellx*tmpx
        segsy = r_celly*tmpy
        segsz = r_cellz*tmpz
!       ### check periodic boundary condition
        segsx=segsx-floor(segsx+0.5d0)  ! -1,0,+1
        segsy=segsy-floor(segsy+0.5d0)  ! -1,0,+1
        segsz=segsz-floor(segsz+0.5d0)  ! -1,0,+1
!       ### recover unit
        segrx=cellx*segsx
        segry=celly*segsy
        segrz=cellz*segsz
!       ### segment
        seg_cx(k)=segrx  
        seg_cy(k)=segry
        seg_cz(k)=segrz
!       ### atom
        do i = isegstart, isegend
          xyz(1,i)=xyz(1,i)+segrx-tmpx
          xyz(2,i)=xyz(2,i)+segry-tmpy
          xyz(3,i)=xyz(3,i)+segrz-tmpz
        enddo ! i
      enddo ! k

      return
      end
