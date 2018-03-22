subroutine ffb_main
!use mpi 
include "mpif.h"
use makemesh
implicit none
integer :: icnt

icnt = test()
end subroutine  ffb_main
