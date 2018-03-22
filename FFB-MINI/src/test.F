      program main
      include "mpif.h"
      implicit REAL*4 (A-H,O-Z)
      REAL*4  SU
      REAL*4  C_SU
      REAL*4  C(10)
      INTEGER*4 I(10)
      SU=1.0
      CALL sub2(C_SU,C,I,10)
      !write(*,*) C_SU
      end program main
C
      subroutine sub2(c1,c2,i,n)
      include "mpif.h"
      implicit  REAL*4(A-H,O-Z)
      REAL*8 A(n)
      REAL*4 B(n)
      DIMENSION I(n),C2(n)
      call DDCOMX(1,1,I,I,n,I,I,N,B,B,B,N,1,II,B,B,n)
      end subroutine sub2
