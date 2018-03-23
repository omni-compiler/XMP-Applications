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
subroutine closing
    use mpivar, only: myrank
#include "timing.h90"
    implicit none

#ifdef PROF_MAPROF
    if (myrank == 0) write(*,'(/,a)') '**** timings'
    call maprof_print_time_mpi(TM_MAIN_LOOP,     'Main Loop:     ')
    if (myrank == 0) write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_FMM,           'FMM:           ')
    call maprof_print_time_mpi(TM_ENERGY_DIRECT, 'Energy Direct: ')
    call maprof_print_time_mpi(TM_SHAKE,         'Shake:         ')
    call maprof_print_time_mpi(TM_RATTLE,        'Rattle         ')
    call maprof_print_time_mpi(TM_COMM_DIRECT,   'Comm Direct:   ')
    call maprof_print_time_mpi(TM_MIGRATION,     'Migration:     ')
    call maprof_print_time_mpi(TM_ENE_REDUCTION, 'Ene Reduction: ')
    call maprof_print_time_mpi(TM_OUTPUT,        'Output:        ')
    if (myrank == 0) write(*, '(a)') '---------------'
    call maprof_print_time_mpi(TM_P2M,           'Particle2M:    ')
    call maprof_print_time_mpi(TM_M2M,           'M2M:           ')
    call maprof_print_time_mpi(TM_M2L_L2L,       'M2L+L2L:       ')
    call maprof_print_time_mpi(TM_L2P,           'L2Particle:    ')
    call maprof_print_time_mpi(TM_COMM_FMM,      'Comm FMM       ')

    call maprof_output()
#endif

end subroutine closing
