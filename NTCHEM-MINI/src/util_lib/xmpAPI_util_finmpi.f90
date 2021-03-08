      SUBROUTINE Util_FinMPI
!
      USE MPI_Module, ONLY : MPIMain
      USE XMP_API
!
!     o Terminate MPI
!
      IMPLICIT NONE
!
!coarray
!     INCLUDE "mpif.h"
!
!     INTEGER :: IErr
!
!     CALL MPI_FINALIZE(IErr)
!
      CALL xmp_api_finalize
      IF (MPIMain) THEN
         WRITE(*, '("MPI has been terminated")')
      END IF
!
      END SUBROUTINE
