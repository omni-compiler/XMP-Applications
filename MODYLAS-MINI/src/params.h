/*
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
*/
/*
 *  { "key woard", type, "default value" },
 */
{ "trj_start", CONFIG_INT, "0" },
{ "trj_interval", CONFIG_INT, "1" },
{ "restart_start", CONFIG_INT, "0" },
{ "restart_interval", CONFIG_INT, "1" },
{ "mntr_start", CONFIG_INT, "0" },
{ "mntr_interval", CONFIG_INT, "1" },

{ "randomseed", CONFIG_INT, "1235" },
{ "dt", CONFIG_DOUBLE, NULL },
{ "step", CONFIG_INT, NULL },
{ "maxwell_velocities", CONFIG_BOOL, "no" },
{ "temperature", CONFIG_DOUBLE, "300.0" },

{ "nstep_skip_middle", CONFIG_INT, "1" },
{ "nstep_skip_long", CONFIG_INT, "1" },

{ "manual_division", CONFIG_BOOL, "no" },
{ "nxdiv", CONFIG_INT, "1" },
{ "nydiv", CONFIG_INT, "1" },
{ "nzdiv", CONFIG_INT, "1" },

{ "shake_max_iteration", CONFIG_INT, "100" },
{ "shake_tolerance", CONFIG_DOUBLE, "1.0e-8" },

{ "cutoff", CONFIG_DOUBLE, NULL },

{ "ndirect", CONFIG_INT, "2" },
{ "nmax", CONFIG_INT, "4" },
{ "ULswitch", CONFIG_INT, "1" },

{ "ewald_surface_term", CONFIG_BOOL, "no" },

{ "ncell", CONFIG_INT, NULL },
