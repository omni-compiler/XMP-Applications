!
! Copyright (C) 2014 RIKEN AICS
!
module ConfigRead
  implicit none

  private
  public ConfigRead_parse
  public ConfigRead_get_int, ConfigRead_get_double
  public ConfigRead_get_bool
  public ConfigRead_get_max_str_len
  public ConfigRead_set_max_str_len

  integer :: max_str_len = 256

contains

function to_c_string(f_str) result(c_str)
  character(len=*), intent(in) :: f_str
  character(len=1) :: c_str(len(f_str)+1)
  integer :: i, n

  n = len(f_str)
  do i = 1, n
    c_str(i) = f_str(i:i)
  end do
  c_str(n+1) = char(0)
end function to_c_string


integer function ConfigRead_get_max_str_len()
  ConfigRead_get_max_str_len = max_str_len
end function ConfigRead_get_max_str_len


subroutine ConfigRead_set_max_str_len(len)
  integer, intent(in) :: len
  max_str_len = len
end subroutine ConfigRead_set_max_str_len


subroutine ConfigRead_parse(file)
  character(*), intent(in) :: file
  interface
    subroutine parse(file)
      character(len=1), intent(in) :: file(*)
    end subroutine parse
  end interface
  call parse(to_c_string(file))
end subroutine ConfigRead_parse


integer function ConfigRead_get_int(key)
  character(*), intent(in) :: key
  interface
    integer(kind=4) function get_int(key)
      character(len=1), intent(in) :: key(*)
    end function get_int
  end interface
  ConfigRead_get_int = get_int(to_c_string(key))
end function ConfigRead_get_int


real(8) function ConfigRead_get_double(key)
  character(*), intent(in) :: key
  interface
    real(kind=8) function get_double(key)
      character(len=1), intent(in) :: key(*)
    end function get_double
  end interface
  ConfigRead_get_double = get_double(to_c_string(key))
end function ConfigRead_get_double


logical function ConfigRead_get_bool(key)
  character(*), intent(in) :: key
  interface
    integer(kind=4) function get_bool(key)
      character(len=1), intent(in) :: key(*)
    end function get_bool
  end interface
  if (get_bool(to_c_string(key)) /= 0) then
    ConfigRead_get_bool = .true.
  else
    ConfigRead_get_bool = .false.
  end if
end function ConfigRead_get_bool


end module ConfigRead
