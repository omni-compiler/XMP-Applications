!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran reset field kernel.
!>  @author Wayne Gaudin
!>  @details Copies all of the final end of step filed data to the begining of
!>  step data, ready for the next timestep.

MODULE reset_field_kernel_module

CONTAINS

  SUBROUTINE reset_field_kernel(x_min,x_max,y_min,y_max, &
    density0,                                            &
    density1,                                            &
    energy0,                                             &
    energy1,                                             &
    xvel0,                                               &
    xvel1,                                               &
    yvel0,                                               &
    yvel1)

    IMPLICIT NONE
    use xmp_defs_module

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1

    INTEGER :: j,k
!$xmp nodes p(*,*)
!$xmp template t(x_min-2:x_max+3,y_min-2:y_max+3)
!$xmp distribute t(gblock(blocksize_x), gblock(blocksize_y)) onto p
!$xmp align (i,j) with t(i,j) :: density0, density1, energy0, energy1, xvel0, xvel1, yvel0, yvel1
!$xmp shadow (2:2, 2:2) :: density0, density1, energy0, energy1
!$xmp shadow (2:3, 2:3) :: xvel0, xvel1, yvel0, yvel1
!$xmp save_desc :: p, t, density0, density1, energy0, energy1, xvel0, xvel1, yvel0, yvel1
!!$OMP PARALLEL

!$xmp loop (j,k) on t(j,k)
!$OMP PARALLEL DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        density0(j,k)=density1(j,k)
        energy0(j,k)=energy1(j,k)
      ENDDO
    ENDDO

!$xmp loop (j,k) on t(j,k) expand(0:1,0:1)
!$OMP PARALLEL DO
    DO k=y_min,y_max+1
      DO j=x_min,x_max+1
        xvel0(j,k)=xvel1(j,k)
        yvel0(j,k)=yvel1(j,k)
      ENDDO
    ENDDO
!!$OMP END PARALLEL

  END SUBROUTINE reset_field_kernel

END MODULE reset_field_kernel_module
