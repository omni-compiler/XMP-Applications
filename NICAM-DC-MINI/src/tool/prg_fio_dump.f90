!-------------------------------------------------------------------------------
!
!+  Program FIO dump
!
!-------------------------------------------------------------------------------
program prg_fio_dump
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !      header/data veiwer for new format data
  !      ( packaged NICAM data format : PaNDa )
  !
  !++ Current Corresponding Author: H. Yashiro
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.90      11-09-01  H.Yashiro : [NEW]
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_fio, only : &
    FIO_HLONG,         &
    FIO_LITTLE_ENDIAN, &
    FIO_BIG_ENDIAN,    &
    FIO_DUMP_HEADER,   &
    FIO_DUMP_ALL,      &
    FIO_DUMP_ALL_MORE
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  character(LEN=FIO_HLONG) :: fname     = ""
  integer                  :: mode      = FIO_DUMP_HEADER
  integer                  :: endian    = FIO_BIG_ENDIAN
  logical                  :: filelok   = .false.
  logical                  :: modelok   = .false.
  logical                  :: endianlok = .false.

  character(LEN=FIO_HLONG) :: argstr
  integer :: n, narg

  integer :: command_argument_count

  integer :: fid ! return from C program
  !=============================================================================

  narg = command_argument_count()

  if ( narg == 0 ) then
     write(*,*) 'Usage : fio_dump [option] [file]'
     write(*,*) '  -h show header only'
     write(*,*) '  -d dump all data   '
     write(*,*) '  -e dump all data (60 digit mode)'
     write(*,*) '  -b force dump with big-endian'
     write(*,*) '  -l force dump with little-endian'
     stop
  endif

  do n = 1, narg
     call get_command_argument(n,argstr)

     if ( argstr(1:1) == '-' ) then
        select case(argstr(2:2)) 
        case('h')
           if(.not. modelok) mode = FIO_DUMP_HEADER
           modelok = .true.
        case('d')
           if(.not. modelok) mode = FIO_DUMP_ALL
           modelok = .true.
        case('e') ! [add] 20120621 H.Yashiro
           if(.not. modelok) mode = FIO_DUMP_ALL_MORE
           modelok = .true.
        case('b')
           if(.not. endianlok) endian = FIO_BIG_ENDIAN
           endianlok = .true.
        case('l')
           if(.not. endianlok) endian = FIO_LITTLE_ENDIAN
           endianlok = .true.
        endselect
     else
        if(.not. filelok) fname = trim(argstr)
        filelok = .true.
     endif
  enddo

  call fio_syscheck()

  call fio_register_file(fid,trim(fname))

  call fio_dump_finfo(fid,endian,mode)

end program prg_fio_dump
!-------------------------------------------------------------------------------
