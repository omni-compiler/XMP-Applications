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
      subroutine energy_direct
c---------------------------------------------------------------------
      use trj_org
      use trj_mpi
      use cutoffradius
      use md_charmm_lj
      use md_fmm
      use md_forces
      use md_coulomb
      use md_segment
      use md_periodic
      use md_monitors
      use md_fmm_domdiv_flg
      use md_charmm_lj
      use md_const
      use param
      use mpivar
      use ompvar
      implicit none

      real(8) :: R, Rr6, Rr12, Cij, coef, eps
      real(8) :: xi, yi, zi
      real(8) :: rx, ry, rz, r2, r2_r
      real(8) :: rc, rc_r, rc2_r
      real(8) :: tlx, tly, tlz, tcx, tcy, tcz, tmp
      real(8) :: stlcx, stlcy, stlcz
      real(8) :: chgv_i0
      real(8) :: epsilon_sqrt_i0,R_half_i0
      integer(4) :: i,j,i0,j0
      integer(4) :: jxb,jyb,jzb,ix,iy,iz,ic
      integer(4) :: iam
      real(8) :: Ulj12, Ulj6, Ucoulomb
      real(8) :: sUlj12, sUlj6
      real(8) :: sUcoulomb
      include 'mpif.h' 
!$    include 'omp_lib.h'

      sUlj12 = 0d0
      sUlj6  = 0d0 
      sUcoulomb = 0d0
      iam = 0

!$omp parallel
!$omp&  private(jxb,jyb,jzb,ix,iy,iz)
!$omp&  private(ic,j,iam)
!$omp&  private(i0,j0,i,epsilon_sqrt_i0,R_half_i0,xi,yi,zi)
!$omp&  private(stlcx,stlcy,stlcz,eps,R,rx,ry,rz,r2,rc2_r)
!$omp&  private(r2_r,Rr6,Rr12,coef,tlx,tly,tlz,Ulj6,Ulj12)
!$omp&  private(chgv_i0)
!$omp&  private(rc,rc_r,Cij,tmp,tcx,tcy,tcz)
!$omp&  private(Ucoulomb)
!$omp&  reduction(+:sUlj6,sUlj12,sUcoulomb)
!$    iam = omp_get_thread_num()
      do jxb=1,ncell/nxdiv+4
         do jyb=1,ncell/nydiv+4
            do iz =3,ncell/nzdiv+2
               jzb=iz

!++++++++++++  Make table for [eps,R,chgv] start
               ic=1
               do j0=tag(jzb-2,jyb,jxb),
     &               tag(jzb+2,jyb,jxb)+na_per_cell(jzb+2,jyb,jxb)-1
                  j=m2i(j0)
                  epsilon_sqrt_table(ic,iam)=epsilon_sqrt(paranum(j))
                  R_half_table(ic,iam)=R_half(paranum(j))
                  chgv_table(ic,iam)=md_QQ_4PiE*chgv(paranum(j))
                  ic=ic+1
               enddo
!++++++++++++  Make table for [eps,R,chgv] ended

!++++++++++++  Main calculation start
               do ix=3,ncell/nxdiv+2
                  do iy=3,ncell/nydiv+2
                     if(abs(jxb-ix)<=2 .and. abs(jyb-iy)<=2) then
!
!+++++++++ Lennard-Jones and Coulomb force calculation start ++++++++++
!
                     if(iy/=jyb .or. ix/=jxb) then
!  Not include own cell
!$omp do
                        do i0=tag(iz,iy,ix),
     &                        tag(iz,iy,ix)+na_per_cell(iz,iy,ix)-1
                           i=m2i(i0)
                           epsilon_sqrt_i0=epsilon_sqrt(paranum(i))
                           R_half_i0=R_half(paranum(i))
                           chgv_i0=chgv(paranum(i))
                           xi=wkxyz(1,i0)
                           yi=wkxyz(2,i0)
                           zi=wkxyz(3,i0)

                           ic=1
                           stlcx=0.d0
                           stlcy=0.d0
                           stlcz=0.d0
!ocl nounroll
                           do j0=tag(jzb-2,jyb,jxb),
     &                           tag(jzb+2,jyb,jxb)
     &                           + na_per_cell(jzb+2,jyb,jxb)-1
                              rx=xi-wkxyz(1,j0)
                              ry=yi-wkxyz(2,j0)
                              rz=zi-wkxyz(3,j0)
                              r2=rx*rx+ry*ry+rz*rz
                              r2_r=1.d0/r2
!------ Lennard-Jones part start
                              ! ^^^ spherical cut-off ^^^
                              if(r2<=cutrad2) then
                                 eps=epsilon_sqrt_i0
     &                               *epsilon_sqrt_table(ic,iam)
                              else
                                 eps=0d0
                              endif !cut-off
                              R=R_half_i0+R_half_table(ic,iam)
                              Rr6=R * R * r2_r
                              Rr6=Rr6 * Rr6 * Rr6
                              Rr12=Rr6 * Rr6
                              coef=12.d0 * eps * r2_r * (Rr12-Rr6)
                              tlx=coef*rx
                              tly=coef*ry
                              tlz=coef*rz
                              ! ^^^ potential ^^^
                              Ulj12=     eps*Rr12
                              Ulj6 =-2d0*eps*Rr6
                              sUlj12=sUlj12+Ulj12
                              sUlj6 =sUlj6 +Ulj6
                              ! ^^^ force ^^^
                              stlcx=stlcx+tlx
                              stlcy=stlcy+tly
                              stlcz=stlcz+tlz
!------ Coulomb part start
                              rc =sqrt(r2)
                              rc_r=1.d0/rc
                              rc2_r=rc_r*rc_r
                              Cij=chgv_i0*chgv_table(ic,iam)
                              Cij=Cij*rc_r
                              tmp=Cij*rc2_r
                              tcx=tmp*rx
                              tcy=tmp*ry
                              tcz=tmp*rz
                              Ucoulomb=Cij
                              sUcoulomb=sUcoulomb+Ucoulomb
                              stlcx=stlcx+tcx
                              stlcy=stlcy+tcy
                              stlcz=stlcz+tcz
                              ic=ic+1
                           enddo !j0
                           w3_f(1,i0,0)=w3_f(1,i0,0)+stlcx
                           w3_f(2,i0,0)=w3_f(2,i0,0)+stlcy
                           w3_f(3,i0,0)=w3_f(3,i0,0)+stlcz
                        enddo !i0
!$omp end do

                     else  !own_cell==other_cell

!  include own cell
!$omp do
                       do i0=tag(iz,iy,ix),
     &                       tag(iz,iy,ix)+na_per_cell(iz,iy,ix)-1
                           i=m2i(i0)
                           epsilon_sqrt_i0=epsilon_sqrt(paranum(i))
                           R_half_i0=R_half(paranum(i))
                           chgv_i0=chgv(paranum(i))
                           xi=wkxyz(1,i0)
                           yi=wkxyz(2,i0)
                           zi=wkxyz(3,i0)

                           ic=1
                           stlcx=0.d0
                           stlcy=0.d0
                           stlcz=0.d0
!ocl nounroll
                           do j0=tag(jzb-2,jyb,jxb),
     &                           tag(jzb+2,jyb,jxb)
     &                           + na_per_cell(jzb+2,jyb,jxb)-1
                              rx=xi-wkxyz(1,j0)
                              ry=yi-wkxyz(2,j0)
                              rz=zi-wkxyz(3,j0)
                              if(i0/=j0) then
                              r2=rx*rx+ry*ry+rz*rz
                              r2_r=1.d0/r2
!------ Lennard-Jones part start
                              ! ^^^ spherical cut-off ^^^
                              if(r2<=cutrad2) then
                                 eps=epsilon_sqrt_i0
     &                               *epsilon_sqrt_table(ic,iam)
                              else
                                 eps=0d0
                              endif !cut-off
                              R=R_half_i0+R_half_table(ic,iam)
                              Rr6=R * R * r2_r
                              Rr6=Rr6 * Rr6 * Rr6
                              Rr12=Rr6 * Rr6
                              coef=12.d0 * eps * r2_r * (Rr12-Rr6)
                              tlx=coef*rx
                              tly=coef*ry
                              tlz=coef*rz
                              ! ^^^ potential ^^^
                              Ulj12=     eps*Rr12
                              Ulj6 =-2d0*eps*Rr6
                              sUlj12=sUlj12+Ulj12
                              sUlj6 =sUlj6 +Ulj6
                              ! ^^^ force ^^^
                              stlcx=stlcx+tlx
                              stlcy=stlcy+tly
                              stlcz=stlcz+tlz
!------ Coulomb part start
                              rc =sqrt(r2)
                              rc_r=1.d0/rc
                              rc2_r=rc_r*rc_r
                              Cij=chgv_i0*chgv_table(ic,iam)
                              Cij=Cij*rc_r
                              tmp=Cij*rc2_r
                              tcx=tmp*rx
                              tcy=tmp*ry
                              tcz=tmp*rz
                              Ucoulomb=Cij
                              sUcoulomb=sUcoulomb+Ucoulomb
                              stlcx=stlcx+tcx
                              stlcy=stlcy+tcy
                              stlcz=stlcz+tcz
                              endif !own_cell/=other_cell
                              ic=ic+1
                           enddo !j0
                           w3_f(1,i0,0)=w3_f(1,i0,0)+stlcx
                           w3_f(2,i0,0)=w3_f(2,i0,0)+stlcy
                           w3_f(3,i0,0)=w3_f(3,i0,0)+stlcz
                       enddo !i0
!$omp end do
                    endif  !own_cell/=other_cell
!
!+++++++++ Lennard-Jones and Coulomb force calculation ended ++++++++++
                    endif !jxb-ix
                 enddo !iy
              enddo !ix
!++++++++++++  Main calculation ended
           enddo !iz
         enddo !jyb
      enddo !jxb
!$omp end parallel
 
      sUlj12=sUlj12*0.5d0
      sUlj6 =sUlj6 *0.5d0
      sUcoulomb=sUcoulomb*0.5d0

      wk_p_energy = wk_p_energy + (sUlj12+sUlj6) + sUcoulomb

      return
      end
c----------------------------------------------------------------------
      subroutine remove_void123lj()
c----------------------------------------------------------------------
      use trj_mpi
      use md_fmm
      use md_fmm_domdiv_flg
      use md_charmm_lj
      use md_void
      use param
      use md_monitors
      use md_const
      use md_forces
      use md_segment
      use mpivar
      implicit none
      integer(4) :: ii,i,i0,j,j0,ipar,jpar
      integer(4) :: isp
      integer(4) :: iam
      integer(4) :: icx0,icy0,icz0
      real(8) :: epsi,Ri,eij,Rij
      real(8) :: xli,yli,zli,rx,ry,rz,r2
      real(8) :: r2_r, Rr6, Rr12, rUlj6,rUlj12, coef_lj
      real(8) :: tlx,tly,tlz
      real(8) :: void_scale=0.5d0 !! double-counted (void)
!$    include 'omp_lib.h'

      if(nvoid==0)return

      rUlj6 =0d0
      rUlj12=0d0

      iam = 0
!$omp parallel default(none)
!$omp& private(ii,iam,i0,j0,isp)
!$omp& private(icx0,icy0,icz0)
!$omp& private(i,ipar,epsi,Ri)
!$omp& private(j,jpar,eij,Rij)
!$omp& private(xli,yli,zli,rx,ry,rz)
!$omp& private(r2,r2_r,Rr6,Rr12,coef_lj)
!$omp& private(tlx,tly,tlz)
!$omp& shared(void_n,void_atom2)
!$omp& shared(m2i,i2m,paranum,epsilon_sqrt,R_half)
!$omp& shared(wkxyz,w3_f)
!$omp& shared(tag,na_per_cell)
!$omp& shared(lxdiv,lydiv,lzdiv)
!$omp& reduction(+:rUlj12,rUlj6)
!$    iam = omp_get_thread_num()
      do ii=1,lxdiv*lydiv*lzdiv
      icz0=mod(ii-1,lzdiv)   +3   
      icy0=mod(ii-1,lzdiv*lydiv)
      icy0=icy0/lzdiv      +3   
      icx0=(ii-1)/(lzdiv*lydiv)+3
!$omp do
      do i0=tag(icz0,icy0,icx0),
     &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
        i   = m2i(i0)
        ipar= paranum(i)
        epsi= epsilon_sqrt(ipar)
        Ri  = R_half(ipar)
        xli = wkxyz(1,i0)
        yli = wkxyz(2,i0)
        zli = wkxyz(3,i0)
        do isp=1,void_n(ipar)
          j    = void_atom2(ipar,isp)+(i-ipar)
          j0   = i2m(j)
          jpar = paranum(j)
          eij = epsi*epsilon_sqrt(jpar) 
          Rij = Ri +R_half(jpar)
          rx = xli - wkxyz(1,j0)
          ry = yli - wkxyz(2,j0)
          rz = zli - wkxyz(3,j0)
          !^^^ LJ ^^^
          r2   = rx**2+ry**2+rz**2
          r2_r = 1d0/r2
          Rr6  = Rij*Rij*r2_r
          Rr6  = Rr6*Rr6*Rr6
          Rr12 = Rr6*Rr6
          rUlj12 = rUlj12 + eij*Rr12
          rUlj6  = rUlj6  + eij*(-2d0*Rr6)
          coef_lj = 12d0*eij*r2_r*(Rr12-Rr6)
          tlx=coef_lj*rx
          tly=coef_lj*ry
          tlz=coef_lj*rz
          w3_f(1,i0,iam)=w3_f(1,i0,iam)-tlx
          w3_f(2,i0,iam)=w3_f(2,i0,iam)-tly
          w3_f(3,i0,iam)=w3_f(3,i0,iam)-tlz
        enddo ! isp
      enddo ! i0
!$omp end do
      enddo ! ii
!$omp end parallel

      rUlj6 =void_scale*rUlj6
      rUlj12=void_scale*rUlj12

      wk_p_energy = wk_p_energy - (rUlj6+rUlj12)

      return
      end
c----------------------------------------------------------------------
      subroutine remove_void123cl()
c----------------------------------------------------------------------
      use trj_mpi
      use md_fmm
      use md_fmm_domdiv_flg
      use md_void
      use param
      use md_monitors
      use md_const
      use md_coulomb
      use md_forces
      use md_segment
      use mpivar
      implicit none
      integer(4) :: ii,i,i0,j,j0,ipar,jpar
      integer(4) :: isp
      integer(4) :: iam
      integer(4) :: icx0,icy0,icz0
      real(8) :: Ci,Cij
      real(8) :: xci,yci,zci,rcx,rcy,rcz
      real(8) :: rc_r, rc2_r, rUcl, coef_cl
      real(8) :: tcx,tcy,tcz,tmp_mdQQ4PiE
      real(8) :: void_scale=0.5d0 !! double-counted (void)
!$    include 'omp_lib.h'

      if(nvoid==0)return

      tmp_mdQQ4PiE=md_QQ_4PiE
      rUcl  =0d0
      iam = 0

!$omp parallel default(none)
!$omp& private(ii,iam,i0,j0,isp)
!$omp& private(icx0,icy0,icz0)
!$omp& private(i,ipar,Ci)
!$omp& private(j,jpar,Cij)
!$omp& private(xci,yci,zci,rcx,rcy,rcz)
!$omp& private(rc_r,rc2_r,coef_cl)
!$omp& private(tcx,tcy,tcz)
!$omp& shared(ndatm,void_n,void_atom2)
!$omp& shared(m2i,i2m,paranum,chgv,tmp_mdQQ4PiE)
!$omp& shared(wkxyz,w3_f)
!$omp& shared(lxdiv,lydiv,lzdiv)
!$omp& shared(tag,na_per_cell)
!$omp& reduction(+:rUcl)
!$    iam = omp_get_thread_num()
      do ii=1,lxdiv*lydiv*lzdiv
      icz0=mod(ii-1,lzdiv)   +3   
      icy0=mod(ii-1,lzdiv*lydiv)
      icy0=icy0/lzdiv      +3   
      icx0=(ii-1)/(lzdiv*lydiv)+3
!$omp do
      do i0=tag(icz0,icy0,icx0),
     &      tag(icz0,icy0,icx0)+na_per_cell(icz0,icy0,icx0)-1
        i   = m2i(i0)
        ipar= paranum(i)
        Ci  = tmp_mdQQ4PiE*chgv(ipar) 
        xci = wkxyz(1,i0)
        yci = wkxyz(2,i0)
        zci = wkxyz(3,i0)
        do isp=1,void_n(ipar)
          j    = void_atom2(ipar,isp)+(i-ipar)
          j0   = i2m(j)
          jpar = paranum(j)
          Cij = Ci *chgv(jpar) 
          rcx = xci - wkxyz(1,j0)
          rcy = yci - wkxyz(2,j0)
          rcz = zci - wkxyz(3,j0)
          !^^^ Coulomb ^^^
          rc2_r = 1d0/(rcx**2+rcy**2+rcz**2)
          rc_r  = sqrt(rc2_r)
          rUcl = rUcl + Cij*rc_r
          coef_cl = Cij*rc_r*rc2_r
          tcx=coef_cl*rcx
          tcy=coef_cl*rcy
          tcz=coef_cl*rcz
          w3_f(1,i0  ,iam)=w3_f(1,i0  ,iam)-tcx
          w3_f(2,i0  ,iam)=w3_f(2,i0  ,iam)-tcy
          w3_f(3,i0  ,iam)=w3_f(3,i0  ,iam)-tcz
        enddo ! isp
      enddo ! i0
!$omp end do
      enddo ! ii
!$omp end parallel

      rUcl  =void_scale*rUcl

      wk_p_energy = wk_p_energy - (rUcl)

      return
      end
