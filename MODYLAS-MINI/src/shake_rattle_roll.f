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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     pshake, shake/roll
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shake_roll(dtin)
      use pshake
      use shakerattleroll
      use trj_mpi
      implicit none
      real(8) :: shktol
      integer(4) :: iterchk, iam
      integer(4) :: i_ps,l0,l1,iconstraints
      real(8) :: dtin,rdtin,rdtbin,rdt2in
!$    include 'omp_lib.h'

      iam=0

      rdtin =1d0/dtin
      rdtbin=rdtin+rdtin
      rdt2in=rdtin*rdtin
!
      iterchk=0
      shktol = shake_tolerance

!$omp parallel default(shared)
!$omp& private(iam,l0,l1,i_ps,iconstraints)
!$    iam=omp_get_thread_num()
!$omp do
      do l0=1,l0max
        i_ps = type_psL(l0)
          call pshake_roll_main(l0, rdtin, rdtbin,
     &                                 iterchk,shktol,i_ps, iam)
      end do
!$omp end do
!$omp end parallel
      if(iterchk.ne.0) then
        write(*,*) 'SHAKE/ROLL iteration no convergence.'
        call mpistop()
      end if
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     prattle, rattle/roll
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rattle_roll(dtin)
      use shakerattleroll
      use pshake
      implicit none
      real(8) :: shktol
      integer(4) :: iterchk, iam
      integer(4) :: i_ps,l0,l1,iconstraints
      real(8) :: dtin,rdtin,rdtbin
!$    include 'omp_lib.h'

      iam = 0
!
      rdtin =1d0/dtin
      rdtbin=rdtin+rdtin
!
      iterchk=0
      shktol = shake_tolerance

!$omp parallel default(shared)
!$omp& private(iam,l0,l1,i_ps,iconstraints)
!$    iam = omp_get_thread_num()
!$omp do
      do l0=1,l0max
        i_ps = type_psL(l0)
          call prattle_roll_main(l0,dtin,rdtbin,iterchk,shktol,i_ps,iam)
      end do
!$omp end do
!$omp end parallel
      if(iterchk.ne.0) then
        write(*,*) 'RATTLE/ROLL iteration no convergence.'
        call mpistop()
      end if
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     pshake/roll main
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pshake_roll_main(l0, rdtin, rdtbin,
     &                             iterchk,shktol,i_ps, iam)
      use pshake
      use shakerattleroll
      use trj_mpi
      implicit none
      integer(4) :: i,ida,idb, j,l0,l1,iconstraints
      integer(4) :: icycle
      real(8) :: rdtin,rdtbin
      real(8) :: shktol,spro
      real(8) :: rij2
      real(8) :: gmi,gmj
      real(8) ::deltax(n_atom_max),deltay(n_atom_max),deltaz(n_atom_max)
      integer(4) :: iflgloop
      integer(4) :: iterchk, iam
      integer(4) :: i_ps, n_atom
      real(8) ::dx0_ps(n_cnst_max),dy0_ps(n_cnst_max),dz0_ps(n_cnst_max)
      real(8) ::dx1_ps(n_cnst_max),dy1_ps(n_cnst_max),dz1_ps(n_cnst_max)
      real(8) :: sabun_ps(n_cnst_max)
      real(8) :: tmp_ps(n_cnst_max)
      real(8) :: v_mat(3,n_cnst_max,n_atom_max)
      real(8) :: v_mat_old(3,n_cnst_max,n_atom_max)
      integer(4) :: atom_ps(n_cnst_max), l, k
      integer(4) :: couple_wk(2,n_cnst_max)
      real(8) :: mass_wk(n_atom_max)
      real(8) :: ma

      iconstraints=ibseL(l0)

      n_atom = n_atom_ps(i_ps)
      atom_ps(1) = shkijL(1,1,l0)

      v_mat_old(:,:,:) = 0.0d0

      do l1=1, iconstraints
        couple_wk(1,l1) = couple_ps(1,l1,i_ps)
        couple_wk(2,l1) = couple_ps(2,l1,i_ps)
      enddo

      do i = 2, n_atom
        atom_ps(i) = atom_ps(1) + i - 1
      enddo

      do i=1, n_atom
        mass_wk(i) = mass_ps(i,i_ps)
      enddo

      do l1= 1, iconstraints !n_pshake = iconstraints
        ida = shkijL(1,l1,l0)
        idb = shkijL(2,l1,l0)
        gmi = rmassL(1,l1,l0)
        gmj = rmassL(2,l1,l0)
        dx0_ps(l1) = xyzstr(1,idb) - xyzstr(1,ida)
        dy0_ps(l1) = xyzstr(2,idb) - xyzstr(2,ida)
        dz0_ps(l1) = xyzstr(3,idb) - xyzstr(3,ida)
        v_mat_old(1,l1,couple_wk(1,l1)) = dx0_ps(l1)*gmi
        v_mat_old(2,l1,couple_wk(1,l1)) = dy0_ps(l1)*gmi
        v_mat_old(3,l1,couple_wk(1,l1)) = dz0_ps(l1)*gmi
        v_mat_old(1,l1,couple_wk(2,l1)) =-dx0_ps(l1)*gmj
        v_mat_old(2,l1,couple_wk(2,l1)) =-dy0_ps(l1)*gmj
        v_mat_old(3,l1,couple_wk(2,l1)) =-dz0_ps(l1)*gmj
      enddo

      do j=1, n_atom
      v_mat(:,:,j) = 0.0d0
      do i=1, iconstraints
      do l=1, iconstraints
      do k=1, 3
        v_mat(k,i,j) = v_mat(k,i,j) +
     &                   v_mat_old(k,l,j) * a_0_sym(l,i,i_ps)
      enddo
      enddo
      enddo
      enddo
    
      do icycle=1,maxshakecycle
        iflgloop = 0

        do l1 = 1, iconstraints
          ida = shkijL(1,l1,l0)
          idb = shkijL(2,l1,l0)
          dx1_ps(l1) = wkxyz(1,idb) - wkxyz(1,ida)
          dy1_ps(l1) = wkxyz(2,idb) - wkxyz(2,ida)
          dz1_ps(l1) = wkxyz(3,idb) - wkxyz(3,ida)
          rij2 = dx1_ps(l1)*dx1_ps(l1) 
     &         + dy1_ps(l1)*dy1_ps(l1) 
     &         + dz1_ps(l1)*dz1_ps(l1)
          sabun_ps(l1) = dij2L(l1,l0) - rij2
          !judge of a constraint for pshake
          if (dabs(sabun_ps(l1))/dij2L(l1,l0) .le. shktol) then !(relative val)**2
            iflgloop=iflgloop+1
          end if
        enddo

        if(iflgloop.eq.iconstraints) then
          goto 123
        endif

        do l1 = 1, iconstraints
          ida = couple_wk(1,l1)
          idb = couple_wk(2,l1)
          spro = dx1_ps(l1)*(v_mat(1,l1,ida)-v_mat(1,l1,idb))+
     &          dy1_ps(l1)*(v_mat(2,l1,ida)-v_mat(2,l1,idb))+
     &          dz1_ps(l1)*(v_mat(3,l1,ida)-v_mat(3,l1,idb))
          tmp_ps(l1) = -sabun_ps(l1) * 0.5d0 / spro
        end do

        do i = 1, n_atom
          deltax(i) = 0.0d0
          deltay(i) = 0.0d0
          deltaz(i) = 0.0d0
          do l1=1, iconstraints
            deltax(i) = deltax(i) + v_mat(1,l1,i) * tmp_ps(l1)
            deltay(i) = deltay(i) + v_mat(2,l1,i) * tmp_ps(l1)
            deltaz(i) = deltaz(i) + v_mat(3,l1,i) * tmp_ps(l1)
          enddo
        enddo 
!ocl norecurrence(wkxyz,wkv)
        do i = 1, n_atom
          ida = atom_ps(i)
          ma = mass_wk(i)
          wkxyz(1,ida) = wkxyz(1,ida) + deltax(i)
          wkxyz(2,ida) = wkxyz(2,ida) + deltay(i)
          wkxyz(3,ida) = wkxyz(3,ida) + deltaz(i)
          deltax(i) = deltax(i) * rdtin
          deltay(i) = deltay(i) * rdtin
          deltaz(i) = deltaz(i) * rdtin
          wkv(1,ida) = wkv(1,ida) + deltax(i)
          wkv(2,ida) = wkv(2,ida) + deltay(i)
          wkv(3,ida) = wkv(3,ida) + deltaz(i)
          deltax(i) = deltax(i) * rdtbin * ma
          deltay(i) = deltay(i) * rdtbin * ma
          deltaz(i) = deltaz(i) * rdtbin * ma
        enddo

      end do
123   continue

      if (icycle .ge. maxshakecycle) then
        iterchk = iterchk + 1
      end if
      
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     prattle/roll main
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prattle_roll_main(
     &             l0,dtin,rdtbin,iterchk,shktol,i_ps,iam)
      use pshake
      use shakerattleroll
      use trj_mpi
      implicit none
      real(8) :: dtin, rdtbin
      integer(4) :: i,ida,idb, j,l0,l1,iconstraints
      integer(4) :: icycle, iam
      real(8) :: dvx,dvy,dvz
      real(8) :: shktol,spro
      real(8) :: gmi,gmj
      real(8) :: tmp, tmp2
      real(8) ::deltax(n_atom_max),deltay(n_atom_max),deltaz(n_atom_max)
      real(8) :: v11,v22,v33
      integer(4) :: iflgloop
      integer(4) :: iterchk
      integer(4) :: i_ps, n_atom
      real(8) ::dx1_ps(n_cnst_max),dy1_ps(n_cnst_max),dz1_ps(n_cnst_max)
      real(8) :: tmp2_ps(n_cnst_max)
      real(8) :: v_mat(3,n_cnst_max,n_atom_max)
      real(8) :: v_mat_old(3,n_cnst_max,n_atom_max)
      integer(4) :: atom_ps(n_cnst_max), l, k
      integer(4) :: couple_wk(2,n_cnst_max)
      real(8) :: ma
      real(8) :: mass_wk(n_atom_max)

      iconstraints=ibseL(l0)

      n_atom = n_atom_ps(i_ps)
      atom_ps(1) = shkijL(1,1,l0)

      v_mat_old(:,:,:) = 0.0d0

      do l1=1, iconstraints
        couple_wk(1,l1) = couple_ps(1,l1,i_ps)
        couple_wk(2,l1) = couple_ps(2,l1,i_ps)
      enddo

      do l1 = 2, n_atom
        atom_ps(l1) = atom_ps(1) + l1 - 1
      enddo

      do i=1, n_atom
        mass_wk(i) = mass_ps(i,i_ps)
      enddo

      do l1=1, iconstraints !n_pshake = iconstraints
        ida = shkijL(1,l1,l0)
        idb = shkijL(2,l1,l0)
        gmi = rmassL(1,l1,l0)
        gmj = rmassL(2,l1,l0)
        dx1_ps(l1) = wkxyz(1,idb) - wkxyz(1,ida)
        dy1_ps(l1) = wkxyz(2,idb) - wkxyz(2,ida)
        dz1_ps(l1) = wkxyz(3,idb) - wkxyz(3,ida)
        v_mat_old(1,l1,couple_wk(1,l1)) = dx1_ps(l1)*gmi
        v_mat_old(2,l1,couple_wk(1,l1)) = dy1_ps(l1)*gmi
        v_mat_old(3,l1,couple_wk(1,l1)) = dz1_ps(l1)*gmi
        v_mat_old(1,l1,couple_wk(2,l1)) =-dx1_ps(l1)*gmj
        v_mat_old(2,l1,couple_wk(2,l1)) =-dy1_ps(l1)*gmj
        v_mat_old(3,l1,couple_wk(2,l1)) =-dz1_ps(l1)*gmj
      enddo

      do j=1, n_atom
      v_mat(:,:,j) = 0.0d0
      do i=1, iconstraints
      do l=1, iconstraints
      do k=1, 3
        v_mat(k,i,j) = v_mat(k,i,j)
     &               + v_mat_old(k,l,j) * a_0_sym(l,i,i_ps)
      enddo
      enddo
      enddo
      enddo

!
!       apply the rattle algorithm to correct the velocities
!      
      do icycle=1,maxshakecycle
        iflgloop=0

        do l1 = 1, iconstraints
          ida = shkijL(1,l1,l0)
          idb = shkijL(2,l1,l0)
          dvx = wkv(1,idb) - wkv(1,ida)
          dvy = wkv(2,idb) - wkv(2,ida)
          dvz = wkv(3,idb) - wkv(3,ida)
          spro = dx1_ps(l1)*dvx + dy1_ps(l1)*dvy + dz1_ps(l1)*dvz 
          tmp2 = spro / dij2L(l1,l0)
          tmp2_ps(l1) = spro * rdij2_ps(l1,i_ps) ! unit: 1/s
          tmp = tmp2*dtin              ! unit: none
          if (abs(tmp) .le. shktol) then ! absolute value
            iflgloop=iflgloop+1
          endif
        enddo

        if(iflgloop.eq.iconstraints) then
          goto 123
        endif

        do i = 1, n_atom
          deltax(i) = 0.0d0
          deltay(i) = 0.0d0
          deltaz(i) = 0.0d0
          do l1=1, iconstraints
            deltax(i) = deltax(i) + v_mat(1,l1,i) * tmp2_ps(l1)
            deltay(i) = deltay(i) + v_mat(2,l1,i) * tmp2_ps(l1)
            deltaz(i) = deltaz(i) + v_mat(3,l1,i) * tmp2_ps(l1)
          enddo
        enddo 
!ocl norecurrence(wkv)
        do i = 1, n_atom
          ida = atom_ps(i)
          ma = mass_wk(i)
          wkv(1,ida) =  wkv(1,ida) + deltax(i)
          wkv(2,ida) =  wkv(2,ida) + deltay(i)
          wkv(3,ida) =  wkv(3,ida) + deltaz(i)
          deltax(i) = deltax(i) * rdtbin * ma
          deltay(i) = deltay(i) * rdtbin * ma
          deltaz(i) = deltaz(i) * rdtbin * ma
          v11 = wkxyz(1,ida) * deltax(i)
          v22 = wkxyz(2,ida) * deltay(i)
          v33 = wkxyz(3,ida) * deltaz(i)
          !v21 = -wkxyz(2,ia) * deltax(i)
          !v31 = -wkxyz(3,ia) * deltax(i)
          !v32 = -wkxyz(3,ia) * deltay(i)
        enddo

      end do
123   continue
      if(icycle.eq.maxshakecycle)then
        iterchk = iterchk + 1
      endif

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     pshake initialize1    2012/10/23
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pshake_initialize1
      use atommass
      use mpivar
      use pshake
      use pshake_init
      use param
      use shakerattleroll
      use trj_mpi
      use trj_org
      implicit none
!      integer(4):: n_type_temp_max=100
      integer(4):: num, GroupNumber, cnst_temp, atom_temp
      integer(4):: ia, ib, ia_max, ia_min
      integer(4):: i, j, k
      integer(4),allocatable:: n_cnst_temp(:), n_atom_temp(:)
      integer(4),allocatable:: couple_ps_temp(:,:,:)
      real(8),allocatable:: r_init_temp(:,:,:), mass_ps_temp(:,:)
      integer(4):: count_type_ps, i_atom, i_cnst, i_type
      integer(4):: nl_cnst_max, nl_atom_max
      logical,allocatable:: ps_flag(:)
      include 'mpif.h'

      n_type_ps = 0
      if(totnconst==0) return  !whole this subroutine

      allocate(gn2ips(rngrp_ps))
      gn2ips(:) = 0

      allocate(ps_flag(rngrp_ps))
      ps_flag(:) = .true.

      allocate(n_cnst_temp(rngrp_ps))
      allocate(n_atom_temp(rngrp_ps))
      allocate(couple_ps_temp(2,n_cnst_max,rngrp_ps))
      allocate(r_init_temp(3,n_atom_max,rngrp_ps))
      allocate(mass_ps_temp(n_atom_max,rngrp_ps))
      n_cnst_temp(:) = 0
      n_atom_temp(:) = 0
      couple_ps_temp(:,:,:) = 0
      r_init_temp(:,:,:) = 0.0d0
      mass_ps_temp(:,:) = 0.0d0

      nl_cnst_max = 0
      nl_atom_max = 0
      count_type_ps = 0
      do i=1, n
        num = paranum(i)
        GroupNumber = ShakeGroupLeader(num)
        if(GroupNumber == 0)cycle
        
        if(ps_flag(GroupNumber))then
          ps_flag(GroupNumber) = .false.
          j = count_type_ps + 1

          !get the number of constraints in the shake group
          cnst_temp = nconstraints(GroupNumber)
          !not rigid
!kojima
        if(cnst_temp .gt. n_cnst_max)then 
          if(myrank==0)then
            write(*,*)'STOP: cnst_temp .gt. n_cnst_max at GroupNumber=',
     &                  GroupNumber
            write(*,*)'n_cnst_max=', n_cnst_max
            write(*,*)'cnst_temp =', cnst_temp
            write(*,*)
          endif
          stop ! for update_shake_local
        endif
!kojima
          !! less than 3 constraints
          !if(cnst_temp < 3)then
          !  goto 100
          !endif

          !get the number of atoms in the shake group
          ia_max = 0
          ia_min = 100000000
          do i_cnst= 1, cnst_temp
            ia=atom1S(GroupNumber,i_cnst)
            ib=atom2S(GroupNumber,i_cnst)
            if(ia .lt. ia_min)then
              ia_min = ia
            endif
            if(ib .gt. ia_max)then
              ia_max = ib
            endif
          enddo
          atom_temp = ia_max - ia_min + 1
!kojima
        if(atom_temp .gt. n_atom_max)then
          if(myrank==0)then
            write(*,*)'WARN: atom_temp .gt. n_atom_max at GroupNumber=',
     &                  GroupNumber
            write(*,*)'n_atom_max=', n_atom_max
            write(*,*)'atom_temp =', atom_temp
            write(*,*)
          endif
          goto 100
          !stop
        endif
!kojima
          !!not rigid
          !if(atom_temp == 3 .and. cnst_temp /= 3)then
          !  goto 100
          !else if(atom_temp == 4 .and. cnst_temp /= 6)then
          !  goto 100
          !endif

          !temporary assignment for P-SHAKE
          !get the numbers of atom & constraint
          n_cnst_temp(j) = cnst_temp
          n_atom_temp(j) = atom_temp
          if(atom_temp > nl_atom_max)then
            nl_atom_max = atom_temp
          endif
          if(cnst_temp > nl_cnst_max)then
            nl_cnst_max = cnst_temp
          endif

          !get the constraint couples in the SHAKE group
          do i_cnst=1, n_cnst_temp(j)
            couple_ps_temp(1,i_cnst,j) = 
     &             atom1S(GroupNumber,i_cnst) - ia_min + 1
            couple_ps_temp(2,i_cnst,j) = 
     &             atom2S(GroupNumber,i_cnst) - ia_min + 1
          enddo

          !get mass of atoms in the SHAKE group
          do i_atom=1, n_atom_temp(j) !check
            ia = i + i_atom - 1
!           write(*,*)ia, paranum(ia)
            mass_ps_temp(i_atom,j) = mass(paranum(ia))
            r_init_temp(1,i_atom,j) = xyz(1,ia)
            r_init_temp(2,i_atom,j) = xyz(2,ia)
            r_init_temp(3,i_atom,j) = xyz(3,ia)
          enddo

          count_type_ps = count_type_ps + 1
          gn2ips(GroupNumber) = count_type_ps
        endif

100     continue
        
      enddo

      n_type_ps = count_type_ps

      if(n_type_ps==0)then
        if(myrank==0)then 
          write(*,'(/,a)') '**** PSHAKE has NOT been applied.'
        endif
        return
      end if

      if(myrank==0)then
        write(*,'(/,a, i0, a, i0)')
     &   '**** PSHAKE applied: ', n_type_ps, '/', rngrp_ps
      endif

      !assign the P-SHAKE informations

      allocate(n_cnst_ps(n_type_ps))
      allocate(n_atom_ps(n_type_ps))
      allocate(couple_ps(2,nl_cnst_max,n_type_ps))
      allocate(r_init_ps(3,nl_atom_max,n_type_ps))
      allocate(mass_ps(nl_atom_max,n_type_ps))
      n_cnst_ps(:)=0
      n_atom_ps(:)=0
      couple_ps(:,:,:)=0
      r_init_ps(:,:,:)=0.0d0
      mass_ps(:,:)=0.0d0

      allocate(a_0(nl_cnst_max,nl_cnst_max,n_type_ps))
      allocate(a_0_sym(nl_cnst_max,nl_cnst_max,n_type_ps))
      allocate(rdij2_ps(nl_cnst_max,n_type_ps))
      a_0(:,:,:)=0.0d0
      a_0_sym(:,:,:)=0.0d0
      rdij2_ps(:,:)=0.0d0

      do i_type=1, n_type_ps
        n_cnst_ps(i_type) = n_cnst_temp(i_type)
        n_atom_ps(i_type) = n_atom_temp(i_type)
        do i_cnst=1, n_cnst_temp(i_type)
          couple_ps(1,i_cnst,i_type)
     $  = couple_ps_temp(1,i_cnst,i_type)
          couple_ps(2,i_cnst,i_type)
     $  = couple_ps_temp(2,i_cnst,i_type)
        enddo
        do i_atom=1, n_atom_temp(i_type)
          do k=1, 3
            r_init_ps(k,i_atom,i_type)
     &  =   r_init_temp(k,i_atom,i_type)
          enddo
        enddo
        do i_atom=1, n_atom_temp(i_type)
          mass_ps(i_atom,i_type) = mass_ps_temp(i_atom,i_type)
        enddo
      enddo

      deallocate(ps_flag)
      deallocate(n_cnst_temp)
      deallocate(n_atom_temp)
      deallocate(couple_ps_temp)
      deallocate(r_init_temp)
      deallocate(mass_ps_temp)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     finish p-shake 1
!     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pshake_finish1
      use pshake_init
      use pshake
      implicit none

      if(n_type_ps==0) return

      deallocate(n_cnst_ps)
      deallocate(r_init_ps)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     finish p-shake 2
!     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pshake_finish2
      use pshake
      implicit none

      if(n_type_ps==0) return

      deallocate(couple_ps)
      deallocate(n_atom_ps)
      deallocate(a_0)
      deallocate(a_0_sym)
      deallocate(type_psL)
      deallocate(rdij2_ps)
      deallocate(gn2ips)
      deallocate(mass_ps)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     initialize p-shake 2
!     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pshake_initialize2
      use shakerattleroll
      use pshake
      use pshake_init
      implicit none
      integer(4),parameter:: n_piter_max=10
      integer(4):: k, i, j, l, n_piter, i_cnst, i_atom
      integer(4):: ia, ib, ja, jb, i_count, j_count
      real(8):: dr(3,n_cnst_max), dv(3,n_cnst_max,n_cnst_max)
      real(8):: dsdr(3,n_cnst_max,n_atom_max)
      real(8):: r_mat(n_cnst_max-1,n_cnst_max-1), r_clm(n_cnst_max-1)
      real(8):: r_mat_all(n_cnst_max,n_cnst_max)
      real(8):: r_mat_inv(n_cnst_max-1,n_cnst_max-1)
      real(8):: v_pshake(3,n_cnst_max,n_atom_max)
      real(8):: v_pshake_new(3,n_cnst_max,n_atom_max)
      real(8):: v_mat_old(3,n_cnst_max,n_atom_max)
      real(8):: v_mat(3,n_cnst_max,n_atom_max)
      real(8):: a_mat(n_cnst_max,n_cnst_max)
      real(8):: a_0_new(n_cnst_max,n_cnst_max)
      real(8):: a_0_wk(n_cnst_max,n_cnst_max)
      real(8):: gmi, gmj, r_mass_ps(n_atom_max)
      integer(4):: i_ps
      integer(4):: n_cnst_wk, n_atom_wk
      integer(4):: couple_wk(2,n_cnst_max)
      real(8):: r_init_wk(3,n_atom_max)
      real(8):: mass_wk(n_atom_max), len_cnst_wk(n_cnst_max)
      real(8):: temp
      real(8):: det
      !for check number of constraits, N, in all constraint group
      !for version of N = 3 in all group

      do i_ps=1, n_type_ps  !whole subroutine

      !initialize
      couple_wk(:,:) = 0
      r_init_wk(:,:) = 0.0d0
      mass_wk(:) = 0.0d0
      len_cnst_wk(:) = 0.0d0

      n_cnst_wk = n_cnst_ps(i_ps)
      n_atom_wk = n_atom_ps(i_ps)
      !write(*,*)n_cnst_wk, n_atom_wk
      do i_cnst=1, n_cnst_wk
        couple_wk(1,i_cnst)
     $= couple_ps(1,i_cnst,i_ps)
        couple_wk(2,i_cnst)
     $= couple_ps(2,i_cnst,i_ps)
!       write(*,*)'couple:',couple_wk(1,i_cnst), couple_wk(2,i_cnst)
      enddo
      do i_atom=1, n_atom_wk
        do k=1, 3
          r_init_wk(k,i_atom)
     &=   r_init_ps(k,i_atom,i_ps)
        enddo
        !write(*,*)'r:',(r_init_wk(i,i_atom),i=1,3)
      enddo
      do i_atom=1, n_atom_wk
        mass_wk(i_atom) = mass_ps(i_atom,i_ps)
        !write(*,*)'mass:',mass_wk(i_atom)
      enddo

      !check
      do i=1, n_atom_wk
        !mass_pshake(i) = mass_pshake(i) * md_ATOMIC_MASS_UNIT
        r_mass_ps(i) = 1.0d0 / mass_wk(i)
      enddo

      if(n_cnst_wk > 1)then !multi constraints in shake group

        !check
        !dt2 = dt*dt
       
        !calc eq. 6
        do i=1, n_cnst_wk
          ia = couple_wk(1,i)
          ib = couple_wk(2,i)
          !write(*,*)ia, ib
          do k=1, 3
            dr(k,i) = r_init_wk(k,ia) - r_init_wk(k,ib)
          enddo
          !temp = dr(1,i)*dr(1,i) + dr(2,i)*dr(2,i) + dr(3,i)*dr(3,i)
          !temp = dsqrt(temp) 
          !write(*,*)temp, dr(1,i), dr(2,i), dr(3,i)
        enddo
        !write(*,*)
       
        dsdr(:,:,:) = 0.0d0
        do i=1, n_cnst_wk
          ia = couple_wk(1,i)
          ib = couple_wk(2,i)
          do k=1, 3
            dsdr(k,i,ia) = 2.0d0*dr(k,i)
            dsdr(k,i,ib) = -1.0d0*dsdr(k,i,ia)
          enddo
        enddo
       
        do j=1, n_atom_wk
        do i=1, n_cnst_wk
        do k=1, 3
          v_pshake(k,i,j) = dsdr(k,i,j)*r_mass_ps(j)
          !v_pshake(k,i,j) = dsdr(k,i,j)*dt2*r_mass_pshake(j)
          !!There is (dt)^2 at every matrix element.
        enddo
        !write(*,*)v_pshake(1,i,j), v_pshake(2,i,j), v_pshake(3,i,j)
        enddo
        enddo
       
        !initialize A0
        a_0_wk(:,:) = 0.0d0
        do i=1, n_cnst_wk
          a_0_wk(i,i) = 1.0d0
        enddo
        
        do n_piter=1, n_piter_max
       
          !calc eq. 7
          dv(:,:,:) = 0.0d0
          do j=1, n_cnst_wk  !Yoko
            ja = couple_wk(1,j)
            jb = couple_wk(2,j)
            do i=1, n_cnst_wk  !tate
              do k=1, 3
                dv(k,i,j) = v_pshake(k,i,ja) - v_pshake(k,i,jb)
              enddo
            !write(*,*)i, j,dv(1,i,j), dv(2,i,j), dv(3,i,j)
            enddo
          enddo
          !write(*,*)
       
          !calc r matrix (Nc x Nc)
          do j=1, n_cnst_wk  !tate
          do i=1, n_cnst_wk  !Yoko
          r_mat_all(i,j) =
     &     dv(1,i,j)*dr(1,j) + dv(2,i,j)*dr(2,j) + dv(3,i,j)*dr(3,j)
          enddo
          !write(*,*)i, r_mat_all(i,1), r_mat_all(i,2)
          enddo
       
          a_mat(:,:) = 0.0d0
          do i=1, n_cnst_wk
            a_mat(i,i) = 1.0d0
          enddo
       
          !iteration for A0
          do k=1, n_cnst_wk
       
            !calc r matrix (Nc-1 x Nc-1)
            j_count = 1
            do j=1, n_cnst_wk !Yoko
              if(j == k)cycle
              i_count = 1
              do i=1, n_cnst_wk !tate
                if(i == k)cycle
                r_mat(i_count,j_count) = r_mat_all(i,j)
                !write(*,*)'a', i_count, j_count, r_mat(i_count,j_count)
                i_count = i_count + 1
              enddo
              !write(*,*)r_mat(i_count,1), r_mat(i_count,2)
              j_count = j_count + 1
            enddo
            !write(*,*)
       
            !calc r column (Nc-1 x 1) in eq. 29
            j_count = 1
            do j=1, n_cnst_wk
              if(j == k)cycle
              r_clm(j_count) = -r_mat_all(k,j)
              !write(*,*)'b', j_count, r_clm(j_count)
              j_count = j_count + 1
            enddo
       
            !!!! matrix inversion (MiniApp version)
            if (n_cnst_wk == 2) then
              r_mat_inv(1,1) = 1.0d0 / r_mat(1,1)
            else if (n_cnst_wk == 3) then
              det = r_mat(1,1) * r_mat(2,2) - r_mat(1,2) * r_mat(2,1)
              det = 1.0d0 / det
              r_mat_inv(1,1) = r_mat(2,2) * det
              r_mat_inv(2,1) = -r_mat(1,2) * det
              r_mat_inv(1,2) = -r_mat(2,1) * det
              r_mat_inv(2,2) = r_mat(1,1) * det
            else
              write(*,*) 'STOP: n_cnst_wk > 3'
              call mpistop
            end if
       
            !calculate the elements in A
            i_count=1
            do i=1, n_cnst_wk ! tate for A
              if(i == k)then
                a_mat(k,i) = 1.0d0
                cycle
              endif
       
              do j=1, n_cnst_wk-1
                a_mat(k,i) = a_mat(k,i) + r_mat_inv(j,i_count)*r_clm(j)
                !write(*,*)i, j, r_mat_inv(j,i_count), a_mat(k,i)
              enddo
       
              i_count = i_count + 1
            enddo
       
          enddo
       
          !"V <- VA"
          v_pshake_new(:,:,:) = 0.0d0
          a_0_new(:,:) = 0.0d0
          do j=1, n_atom_wk  !tate
          do i=1, n_cnst_wk !Yoko
          do l=1, n_cnst_wk
            do k=1, 3
              v_pshake_new(k,i,j) =
     &        v_pshake_new(k,i,j) + v_pshake(k,l,j)*a_mat(i,l)
            enddo
          enddo
          enddo
          enddo
       
          do j=1, n_atom_wk
          do i=1, n_cnst_wk
          do k=1, 3
            v_pshake(k,i,j) = v_pshake_new(k,i,j)
          enddo
          enddo
          enddo
          
          !"(A_0) <- (A_0)A"
          do j=1, n_cnst_wk
          do i=1, n_cnst_wk
          do l=1, n_cnst_wk
            a_0_new(i,j) = a_0_new(i,j) + a_0_wk(l,j)*a_mat(i,l)
          enddo
          enddo
          enddo
       
          do j=1, n_cnst_wk
          do i=1, n_cnst_wk
            a_0_wk(i,j) = a_0_new(i,j)
          enddo
          enddo
          
       
        enddo
       
        !set transposed matrix of a_0 for cost cut
        do j=1, n_cnst_wk
        do i=1, n_cnst_wk
          a_0(i,j,i_ps) = a_0_wk(i,j)
          a_0_sym(j,i,i_ps) = a_0_wk(i,j)
        enddo
        enddo

      else ! 1 constraints in shake group

        do j=1, n_cnst_wk
        do i=1, n_cnst_wk
          a_0(i,j,i_ps) = 1.0d0
          a_0_sym(j,i,i_ps) = 1.0d0
        enddo
        enddo

      endif

      v_mat_old(:,:,:)=0.0d0
      do i_cnst=1, n_cnst_wk
        ia = couple_wk(1,i_cnst)
        ib = couple_wk(2,i_cnst)
        do k=1, 3
          dr(k,i_cnst) = r_init_wk(k,ia) - r_init_wk(k,ib)
        enddo
        gmi= r_mass_ps(ia)
        gmj= r_mass_ps(ib)
        v_mat_old(1,i_cnst,ia) = dr(1,i_cnst)*gmi
        v_mat_old(2,i_cnst,ia) = dr(2,i_cnst)*gmi
        v_mat_old(3,i_cnst,ia) = dr(3,i_cnst)*gmi
        v_mat_old(1,i_cnst,ib) =-dr(1,i_cnst)*gmj
        v_mat_old(2,i_cnst,ib) =-dr(2,i_cnst)*gmj
        v_mat_old(3,i_cnst,ib) =-dr(3,i_cnst)*gmj
      enddo

      v_mat(:,:,:) = 0.0d0
      do i=1, n_cnst_wk
      do j=1, n_atom_wk
      do l=1, n_cnst_wk
      do k=1, 3
        v_mat(k,i,j) = v_mat(k,i,j) +
     &                   v_mat_old(k,l,j) * a_0_sym(l,i,i_ps)
      enddo
      enddo
      enddo
      enddo

      do i_cnst=1, n_cnst_wk
        ia = couple_wk(1,i_cnst)
        ib = couple_wk(2,i_cnst)
        temp = 
     &   (v_mat(1,i_cnst,ia)-v_mat(1,i_cnst,ib))*dr(1,i_cnst) +
     &   (v_mat(2,i_cnst,ia)-v_mat(2,i_cnst,ib))*dr(2,i_cnst) +
     &   (v_mat(3,i_cnst,ia)-v_mat(3,i_cnst,ib))*dr(3,i_cnst)
        rdij2_ps(i_cnst,i_ps) = 1.0d0/temp
      enddo

      enddo
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_shake_local()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use shakerattleroll
      use md_fmm
      use md_fmm_domdiv_flg
      !pshake
      use pshake
      implicit none

!     maxibn=5
!     maxibn=10
      totnconstL=max(
     &   int(totnconst/(ncell**3)*(lxdiv*lydiv*lzdiv)*2.5d0),30)
      allocate(ibseL(totnconstL))
      allocate(shkijL(2,n_cnst_max,totnconstL))
      allocate(rmassL(2,n_cnst_max,totnconstL))
      allocate(dij2L(n_cnst_max,totnconstL))
      allocate(type_psL(totnconstL))

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine update_shake_local()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use atommass
      use md_fmm
      use md_fmm_domdiv_flg
      use md_segment
      use mpivar
      use param
      use pshake
      use shakerattleroll
      use trj_mpi
      implicit none
      include 'mpif.h'
      integer(4) :: k0,i0,l0,l1
      integer(4) :: i,ibn
      integer(4) :: ia,ib,ida,idb
      integer(4) :: num, GroupNumber, dn1, dn2

!=== copy constants to local arrays ===!
      l0=0
      kshake=0
      do k0=1,nselfseg
      do i0=lsegtop(k0),lsegtop(k0)+lseg_natoms(k0)-1
        if(kshake(i0).ne.0) cycle
        i = m2i(i0) ! G (global)
        num=paranum(i)
        if(ShakeGroupLeader(num)==0) cycle
        GroupNumber=ShakeGroupLeader(num)
        l0=l0+1     
        if(l0.ge.totnconstL)then
           write(*,*) l0,totnconstL
           stop'ERROR: l0.ge.totnconstL'
        endif
        ibseL(l0)=nconstraints(GroupNumber)
        type_psL(l0)=gn2ips(GroupNumber) !pshake

        l1=0
        do ibn = 1, nconstraints(GroupNumber)
          l1=l1+1 ; if(l1.ge.n_cnst_max) stop'ERROR: l1.ge.n_cnst_max'
!         l1=l1+1 ; if(l1.ge.maxibn) stop'ERROR: l1.ge.maxibn'
          dn1=abs(atom1S(GroupNumber,ibn)-num)
          dn2=abs(atom2S(GroupNumber,ibn)-num)
          ia=i + dn1
          ib=i + dn2
          ida = i2m(ia)  ! L (local)
          idb = i2m(ib)
          kshake(ida)=+1
          kshake(idb)=+1
            shkijL(1,l1,l0) = ida
            shkijL(2,l1,l0) = idb
            dij2L(l1,l0) = slength(GroupNumber,ibn)**2
            rmassL(1,l1,l0) = r_mass(paranum(ia)) ! G
            rmassL(2,l1,l0) = r_mass(paranum(ib)) ! G
        enddo ! ibn
      enddo ! i0
      enddo ! ii
      l0max=l0

      return
      end
