module x 
private
!include "mpif.h"
#define MAPROF_ROOT 0
#define MAPROF_AVE  1
#define MAPROF_MIN  2
#define MAPROF_MAX  3
#define MAPROF_SD   4
  implicit none

integer ::i
   interface

    subroutine maprof_time_start(id) !!bind(c)
!      use iso_c_binding
!      integer(kind=4) :: id 
      integer(kind=4) :: id
    end subroutine maprof_time_start

  end interface

contains

function c_string(f_str) result(c_str)
!  use iso_c_binding
  character(len=*), intent(in) :: f_str
!!character(len=1, len=1) :: c_str(len_trim(f_str)+1)
  character(len=1) :: c_str(len(f_str)+1)
  integer :: i, n

!!n = len_trim(f_str)
  n = len(f_str)
  do i = 1, n
    c_str(i) = f_str(i:i)
  end do
  c_str(n+1) = char(0)!c_null_char

end function c_string
end module x
