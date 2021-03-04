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
      subroutine initialize_application
c----------------------------------------------------------------------
      implicit none
      call check_parallel_condition
c     initialize MD objects
      call init_g_main
      call init_md_check
      call init_md_condition
      call init_md_velocity
      call init_fmm_domain_div  
      call check_cutofflength              
      call calc_idcell
      call fmod_set_maxsegments
      call fmod_alloc_metadata  
      call fmod_alloc_multipole 
      call init_fmm
      call init_comm_direct_3() 
      call init_comm_bound() 
      call init_shake_local() 
      call pshake_initialize1
      call pshake_initialize2
      call pshake_finish1

      return
      end
c----------------------------------------------------------------------
      subroutine open_trj
c----------------------------------------------------------------------
      use md_file
      implicit none
      integer(4) :: io
      open(f_trj, file=trim(session_name)// '.mdtrj.bin', iostat=io,   
     &         status='replace', access='sequential',form='unformatted')
      if (io /= 0) then
        call abort_with_message_a('Cannot create mdtrj file.')
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine record_current_trajectory
c----------------------------------------------------------------------
      use trj_org
      use unitcell
      use trj_org
      use nhc
      use md_condition
      use md_file
      use md_periodic
      implicit none

      write(f_trj) mdstep, mdstep*dt
      !  Write coordinates and velocities of atoms.
      write(f_trj) n
      write(f_trj) xyz(1:3,:), v(1:3,:)
      !  Write positions and velocities of thermostats.
      write(f_trj) nnhc
      write(f_trj) rss, vss
      !  Write positions and velocities of barostats.
      write(f_trj) nnhc
      write(f_trj) rssb,vssb
      !  Write cell parameters (length and angles).
      write(f_trj) cellx,celly,cellz,alpha,beta,gamma,vboxg
      call flush(f_trj)

      return
      end
c-----------------------------------------------------------------------
      subroutine open_mntr
c----------------------------------------------------------------------
      use md_file
      implicit none
      integer(4) :: io
      open(f_mntr, file=trim(session_name)// '.mdmntr', iostat=io,  
     &           status='replace', access='sequential',form='formatted')
      if (io /= 0) then
        call abort_with_message_a('Cannot create mdmntr file.')
      endif
      write(f_mntr,'(a,a,a)') 
     & '## ', trim(session_name), '.mdmntr' //
     & ' -- monitor variables output from MD calculation by modylas'
      write(f_mntr,'(a)') '#'
      write(f_mntr,'(a)') '# datas below are formated as:'
      write(f_mntr,'(13a)') 
     &  '#',' step     '          , 
     &      ' time               ', 
     &      ' Hamiltonian        ', 
     &      ' potential-E        ', 
     &      ' kinetic-E          ', 
     &      ' total energy       ', 
     &      ' temperature        '
      write(f_mntr,'(13a)') 
     &  '#', '          '          , 
     &       '              [sec] ', 
     &       '           [J/cell] ', 
     &       '           [J/cell] ', 
     &       '           [J/cell] ', 
     &       '           [J/cell] ', 
     &       '                [K] '
      write(f_mntr,'(a)') '#'
      call flush(f_mntr)
      return
      end
c-----------------------------------------------------------------------
      subroutine record_current_monitors
c----------------------------------------------------------------------
      use unitcell
      use md_condition
      use md_file
      use md_monitors
      implicit none
      
      write(f_mntr,'(i10,6es20.12)') 
     &   mdstep, mdstep*dt, hamiltonian,    
     &   p_energy, k_energy, 
     &   t_energy, temperature
      call flush(f_mntr)

      return
      end
c-----------------------------------------------------------------------
      subroutine open_restart
c----------------------------------------------------------------------
      use md_file
      implicit none
      integer(4) :: io
      open(f_restart_bin, file=trim(session_name)// '.restart.bin', 
     &         iostat=io, status='replace', 
     &         access='sequential',form='unformatted')
      if (io /= 0) then
        call abort_with_message_a('Cannot create restart.bin file.')
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine record_restart_binary
c----------------------------------------------------------------------
      use trj_org
      use unitcell
      use trj_org
      use nhc
      use md_condition
      use md_file
      use md_periodic
      implicit none

      !  Write coordinates and velocities of atoms.
      rewind(f_restart_bin)
      write(f_restart_bin) n
      write(f_restart_bin) xyz(1:3,:), v(1:3,:)
      !  Write positions and velocities of thermostats.
      write(f_restart_bin) nnhc
      write(f_restart_bin) rss, vss
      !  Write positions and velocities of barostats.
      write(f_restart_bin) nnhc
      write(f_restart_bin) rssb, vssb
      !  Write cell parameters (length and angles).
      write(f_restart_bin) cellx,celly,cellz, alpha,beta,gamma, vboxg
      call flush(f_restart_bin)

      return
      end
c-----------------------------------------------------------------------
      subroutine record_current_state
c----------------------------------------------------------------------
      use md_condition
      use g_main
      use md_condition
      use mpivar
      implicit none

!### record mdrun & mdmntr ###
      if (mod((mdstep - mntr_start), mntr_interval) == 0) then
      IF (myrank.eq.mpiout) THEN
        call record_current_monitors
      ENDIF
      endif

      if (mod((mdstep-trj_start),trj_interval)==0.or.
     &    mod((mdstep-restart_start),restart_interval)==0) then
          call pre_record_data
      endif

!### record mdtrj.bin ###
      if (mod((mdstep-trj_start),trj_interval)==0) then
      IF (myrank.eq.mpiout) THEN
          call record_current_trajectory
      ENDIF
      endif

!### record restart.bin ###
      if (mod((mdstep-restart_start),restart_interval)==0)then
      IF (myrank.eq.mpiout) THEN
          call record_restart_binary
      ENDIF
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine check_parallel_condition
c-----------------------------------------------------------------------
      use mpivar
      implicit none
      include 'mpif.h'
      integer(4) :: ierr

!### checm nprocs ###!
      if(mod(nprocs,2).ne.0)then
        if(myrank==0)then
        write(*,*) 'ERROR: nprocs is not equal to 2 powers.'
        endif
!coarray        call mpi_barrier(mpi_comm_world, ierr)
        sync all
!!
        call mpiend
        stop 1
      endif
      if(nprocs.lt.8)then
        if(myrank==0)then
        write(*,*) 'ERROR: nprocs less than 8 is not supported.'
        endif
!coarray        call mpi_barrier(mpi_comm_world, ierr)
        sync all
!!
        call mpiend
        stop 1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine init_g_main
c----------------------------------------------------------------------
      use mpivar
      implicit none
      if (myrank.eq.mpiout) then
c       open *.mdtrj file to which trajectory of MD calculation
c       will be written
        call open_trj
c      
c       open *.mdmntr file to which monitor variables (Hamiltonian etc.)
c       will be written
        call open_mntr
c      
c       open *.restart.bin file to which force log will be written
        call open_restart
        endif
      return
      end
c-----------------------------------------------------------------------
      subroutine generate_velocity
c----------------------------------------------------------------------
      use atommass
      use trj_org
      use trj_org
      use md_condition
      use md_const
      use g_main
      use param
      implicit none
      real(8),allocatable :: randbuf(:)
      real(8) :: w1,w2,M,vG__x,vG__y,vG__z,coef,T
      real(8) :: rand0, rand1, genrand_real2
      integer(4) :: i,nrand
c     generate maxwell distribution of velocity of atom,
c     in which posibility variable sqrt(m / kB / T) * v(x, y or z)
c     obeys to normal distribution with
c      expectation = 0
c      standard deviation = 1
c
c     after random generation of velocity, they are adjusted for
c     total momentum of the sytem to be 0 and temperature to be
c     md_generic__maxwell temperature
c
c    handling special case: maxwell_temperature == 0
      if (maxwell_temperature < 1.0d-50) then
        v = 0.0d0
        return
      endif
c     generate uniform random numbers
      if (mod(n,2) == 0) then
        nrand = 3*n
      else
        nrand = 3*n+1
      endif
      allocate(randbuf(nrand))
      call init_genrand(randomseed)
c     conver uniform random numbers to normal random numbers
      do i=1,nrand-1,2
        rand0 = genrand_real2()
        rand1 = genrand_real2()
        w1=sqrt(-2.0d0*log(rand0))*cos(2.0d0*md_PI*rand1)
        w2=sqrt(-2.0d0*log(rand1))*cos(2.0d0*md_PI*rand0)
        randbuf(i+0) = w1
        randbuf(i+1) = w2
      enddo
c     add velocities
      do i=1,n
        coef = sqrt(md_BOLTZMANN * maxwell_temperature * 
     &               r_mass(paranum(i)))
        v(1,i) = randbuf(i*3-2) * coef
        v(2,i) = randbuf(i*3-1) * coef
        v(3,i) = randbuf(i*3-0) * coef
      enddo
      deallocate(randbuf)
c     velocity of center of mass
      M = 0.0d0
      vG__x = 0.0d0;  vG__y = 0.0d0;  vG__z = 0.0d0
      do i=1,n
        if (mass(paranum(i)) .lt. 1.0e+10) then
          M = M + mass(paranum(i))
          vG__x = vG__x + v(1,i) * mass(paranum(i))
          vG__y = vG__y + v(2,i) * mass(paranum(i))
          vG__z = vG__z + v(3,i) * mass(paranum(i))
        endif
      enddo
      vG__x = vG__x / M
      vG__y = vG__y / M
      vG__z = vG__z / M
c     subtracting momentum of the system
      do i=1,n
        if (mass(paranum(i)) .lt. 1.0e+10) then
          v(1,i) = v(1,i) - vG__x
          v(2,i) = v(2,i) - vG__y
          v(3,i) = v(3,i) - vG__z
        endif
      enddo
c     velocity scaling
      T = 0.0d0
      do i=1,n
        if (mass(paranum(i)) .lt. 1.0e+10) then
          T = T + mass(paranum(i))
     &    *(v(1,i)*v(1,i)+v(2,i)*v(2,i)+v(3,i)*v(3,i))
        endif
      enddo
      T = T * degree_of_freedom_inverse  * rvkbolz
      v = v * sqrt(maxwell_temperature / T)
      return
      end
c-----------------------------------------------------------------------
      subroutine remove_system_momentum
c----------------------------------------------------------------------
      use atommass
      use trj_org
      use trj_org
      use param
      implicit none
      real(8) :: M, vG__x, vG__y, vG__z
      integer(4) :: i
c     velocity of center of mass
      M = 0.0d0
      vG__x = 0.0d0
      vG__y = 0.0d0
      vG__z = 0.0d0
      do i=1, n
        M = M + mass(paranum(i))
        vG__x = vG__x + v(1,i) * mass(paranum(i))
        vG__y = vG__y + v(2,i) * mass(paranum(i))
        vG__z = vG__z + v(3,i) * mass(paranum(i))
      enddo
      vG__x = vG__x / M
      vG__y = vG__y / M
      vG__z = vG__z / M
c     subtracting momentum of the system
      do i=1, n
        v(1,i) = v(1,i) - vG__x
        v(2,i) = v(2,i) - vG__y
        v(3,i) = v(3,i) - vG__z
      enddo
      totalmass=M
      return
      end
c-----------------------------------------------------------------------
      subroutine init_md_velocity
c----------------------------------------------------------------------
      use md_segment
      use g_main
      implicit none
c     generate velocity of atom obeying to Maxwell distribution
      if (maxwell_temperature > 0.0d0 .or. reset_maxwell) then
        call generate_velocity
c             No need to exchange by jrearrange
      else
        call remove_system_momentum
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine init_md_check
c-----------------------------------------------------------------------
      use trj_org
      use cutoffradius
      use md_condition
      use md_file
      use md_periodic
      use md_segment
      use g_main
      use mpivar
      use unitcell
      implicit none
      if(myrank.eq.mpiout) then
c
c     check if the last trajectory will be saved
      if (mod(md_condition__howmany_steps,restart_interval) /= 0) then
        write(6,*) 'WARN: the last trajectory will not be saved'
      endif
c
c     check periodic boundary

c     check validity of cut-off length
      if(cutrad > cellxh) then
        write(6,*) 
     &       'cut-off length for force is greater than a half of cell'
        call mpistop
      endif
      if(cutrad > cellyh) then
        write(6,*) 
     &       'cut-off length for force is greater than a half of cell'
        call mpistop
      endif
      if(cutrad > cellzh) then
        write(6,*) 
     &       'cut-off length for force is greater than a half of cell'
        call mpistop
      endif
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine init_md_condition
c----------------------------------------------------------------------
      use trj_org
      use md_condition
      use shakerattleroll
      implicit none
c     if using SHAKE, degree_of_freedom -= howmany_constraints
c     but this operation ought to be done by md_shake
      degree_of_freedom = n * 3 - 3
      if(totnconst .ne. 0) then
        degree_of_freedom = degree_of_freedom - totnconst
      endif
      degree_of_freedom_inverse = 1.0d0 / degree_of_freedom
      return
      end
c-----------------------------------------------------------------------
      subroutine cleanup
c----------------------------------------------------------------------
      use md_file
      close(f_trj)
      close(f_mntr)
      close(f_restart_bin)
      return
      end
c-----------------------------------------------------------------------
