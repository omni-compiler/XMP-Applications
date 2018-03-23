module mod_perf
  use iso_c_binding
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  private
  public perf_start, perf_stop
  public perf_add_fp_ops, perf_add_ld_ops, perf_add_st_ops
  public perf_add_ld_min_ops, perf_add_st_min_ops
  public perf_set_ops
  public perf_set_fp_ops, perf_set_ld_ops, perf_set_st_ops
  public perf_set_ld_min_ops, perf_set_st_min_ops
  public perf_get_time, perf_get_flops
  public perf_get_throughput, perf_get_effective_throughput
  public perf_print, perf_print_time
#ifdef USE_MPI
  public perf_print_time_mpi, perf_print_time_mpi_full
#endif

  ! PERF_ValueType
  enum, bind(c)
    enumerator :: PERF_ROOT, PERF_AVE, PERF_MIN, PERF_MAX, PERF_SD
  end enum

  interface

    subroutine perf_start(id) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
    end subroutine perf_start

    subroutine perf_stop(id) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
    end subroutine perf_stop

    subroutine perf_add_fp_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine perf_add_fp_ops

    subroutine perf_add_ld_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine perf_add_ld_ops

    subroutine perf_add_st_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine perf_add_st_ops

    subroutine perf_add_ld_min_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine perf_add_ld_min_ops

    subroutine perf_add_st_min_ops(id, ops) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      real(c_double), value :: ops
    end subroutine perf_add_st_min_ops

    real(c_double) function perf_get_time(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function perf_get_time

    real(c_double) function perf_get_flops(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function perf_get_flops

    real(c_double) function perf_get_throughput(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function perf_get_throughput

    real(c_double) function perf_get_effective_throughput(id, type) bind(c)
      use iso_c_binding
      integer(c_int), value :: id 
      integer(c_int), value :: type
    end function perf_get_effective_throughput

  end interface

  interface perf_set_fp_ops
    module procedure perf_set_fp_ops_i4, perf_set_fp_ops_i8, &
         perf_set_fp_ops_r8
  end interface perf_set_fp_ops

  interface perf_set_ld_ops
    module procedure perf_set_ld_ops_i4, perf_set_ld_ops_i8, &
         perf_set_ld_ops_r8
  end interface perf_set_ld_ops

  interface perf_set_st_ops
    module procedure perf_set_st_ops_i4, perf_set_st_ops_i8, &
         perf_set_st_ops_r8
  end interface perf_set_st_ops

  interface perf_set_ld_min_ops
    module procedure perf_set_ld_min_ops_i4, &
         perf_set_ld_min_ops_i8, &
         perf_set_ld_min_ops_r8
  end interface perf_set_ld_min_ops

  interface perf_set_st_min_ops
    module procedure perf_set_st_min_ops_i4, &
         perf_set_st_min_ops_i8, &
         perf_set_st_min_ops_r8
  end interface perf_set_st_min_ops

  interface perf_set_ops
    module procedure perf_set_ops_i4, perf_set_ops_i8, &
         perf_set_ops_r8
  end interface perf_set_ops

contains


subroutine perf_set_fp_ops_r8(id, fp, cnt)
  integer, intent(in) :: id
  real(8), value :: fp
  integer, optional, intent(in) :: cnt

  if (present(cnt)) fp = fp * cnt;
  call perf_add_fp_ops(id, fp)

end subroutine perf_set_fp_ops_r8


subroutine perf_set_fp_ops_i8(id, fp, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: fp
  integer, optional, intent(in) :: cnt

  call perf_set_fp_ops_r8(id, dble(fp), cnt)
  
end subroutine perf_set_fp_ops_i8


subroutine perf_set_fp_ops_i4(id, fp, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: fp
  integer, optional, intent(in) :: cnt

  call perf_set_fp_ops_r8(id, dble(fp), cnt)

end subroutine perf_set_fp_ops_i4


subroutine perf_set_ld_ops_r8(id, ld, cnt)
  integer, intent(in) :: id
  real(8), value :: ld
  integer, optional, intent(in) :: cnt

  if (present(cnt)) ld = ld * cnt
  call perf_add_ld_ops(id, ld)

end subroutine perf_set_ld_ops_r8


subroutine perf_set_ld_ops_i8(id, ld, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: ld
  integer, optional, intent(in) :: cnt

  call perf_set_ld_ops_r8(id, dble(ld), cnt)

end subroutine perf_set_ld_ops_i8


subroutine perf_set_ld_ops_i4(id, ld, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: ld
  integer, optional, intent(in) :: cnt

  call perf_set_ld_ops_r8(id, dble(ld), cnt)

end subroutine perf_set_ld_ops_i4


subroutine perf_set_st_ops_r8(id, st, cnt)
  integer, intent(in) :: id
  real(8), value :: st
  integer, optional, intent(in) :: cnt

  if (present(cnt)) st = st * cnt
  call perf_add_st_ops(id, st)

end subroutine perf_set_st_ops_r8


subroutine perf_set_st_ops_i8(id, st, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: st
  integer, optional, intent(in) :: cnt

  call perf_set_st_ops_r8(id, dble(st), cnt)
  
end subroutine perf_set_st_ops_i8

subroutine perf_set_st_ops_i4(id, st, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: st
  integer, optional, intent(in) :: cnt

  call perf_set_st_ops_r8(id, dble(st), cnt)

end subroutine perf_set_st_ops_i4


subroutine perf_set_ld_min_ops_r8(id, ld_min, cnt)
  integer, intent(in) :: id
  real(8), value :: ld_min
  integer, optional, intent(in) :: cnt

  if (present(cnt)) ld_min = ld_min * cnt
  call perf_add_ld_min_ops(id, ld_min)

end subroutine perf_set_ld_min_ops_r8


subroutine perf_set_ld_min_ops_i8(id, ld_min, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: ld_min
  integer, optional, intent(in) :: cnt

  call perf_set_ld_min_ops_r8(id, dble(ld_min), cnt)

end subroutine perf_set_ld_min_ops_i8


subroutine perf_set_ld_min_ops_i4(id, ld_min, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: ld_min
  integer, optional, intent(in) :: cnt

  call perf_set_ld_min_ops_r8(id, dble(ld_min), cnt)

end subroutine perf_set_ld_min_ops_i4


subroutine perf_set_st_min_ops_r8(id, st_min, cnt)
  integer, intent(in) :: id
  real(8), value :: st_min
  integer, optional, intent(in) :: cnt

  if (present(cnt)) st_min = st_min * cnt
  call perf_add_st_min_ops(id, st_min)

end subroutine perf_set_st_min_ops_r8


subroutine perf_set_st_min_ops_i8(id, st_min, cnt)
  integer, intent(in) :: id
  integer(8), intent(in) :: st_min
  integer, optional, intent(in) :: cnt

  call perf_set_st_min_ops_r8(id, dble(st_min), cnt)

end subroutine perf_set_st_min_ops_i8


subroutine perf_set_st_min_ops_i4(id, st_min, cnt)
  integer, intent(in) :: id
  integer, intent(in) :: st_min
  integer, optional, intent(in) :: cnt

  call perf_set_st_min_ops_r8(id, dble(st_min), cnt)

end subroutine perf_set_st_min_ops_i4


subroutine perf_set_ops_r8(sec, fp, ld, st, ld_min, st_min, cnt)
  integer, intent(in) :: sec
  real(8), intent(in) :: fp, ld, st
  real(8), intent(in) :: ld_min, st_min
  integer, optional, intent(in) :: cnt

  call perf_set_fp_ops_r8(sec, fp, cnt)
  call perf_set_ld_ops_r8(sec, ld, cnt)
  call perf_set_st_ops_r8(sec, st, cnt)
  call perf_set_ld_min_ops_r8(sec, ld_min, cnt)
  call perf_set_st_min_ops_r8(sec, st_min, cnt)

end subroutine perf_set_ops_r8


subroutine perf_set_ops_i8(sec, fp, ld, st, ld_min, st_min, cnt)
  integer, intent(in) :: sec
  integer(8), intent(in) :: fp, ld, st
  integer(8), intent(in) :: ld_min, st_min
  integer, optional, intent(in) :: cnt

  call perf_set_ops_r8(sec, dble(fp), dble(ld), dble(st), &
       dble(ld_min), dble(st_min), cnt)
  
end subroutine perf_set_ops_i8


subroutine perf_set_ops_i4(sec, fp, ld, st, ld_min, st_min, cnt)
  integer, intent(in) :: sec
  integer, intent(in) :: fp, ld, st
  integer, intent(in) :: ld_min, st_min
  integer, optional, intent(in) :: cnt

  call perf_set_ops_r8(sec, dble(fp), dble(ld), dble(st), &
                      dble(ld_min), dble(st_min), cnt)

end subroutine perf_set_ops_i4


subroutine perf_print(id, name)
  integer, intent(in) :: id
  character(*), intent(in) :: name

  integer :: myrank = 0
  real(8) :: time, flops, tput, etput
#ifdef USE_MPI
  integer :: ierr
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
#endif

  time  = perf_get_time(id, PERF_ROOT)
  flops = perf_get_flops(id, PERF_ROOT)
  tput  = perf_get_throughput(id, PERF_ROOT)
  etput = perf_get_effective_throughput(id, PERF_ROOT)

  if (myrank == 0) then

  if (time >  0.0) then
    write(*, '(a, " Time:                 ", f20.3, " (S)")') name, time
  else
    write(*, '(a, " Time:                 ", a20)') name, "UNKNOWN"
  end if
  if (time > 0.0 .and. flops > 0.0) then
    write(*, '(a, " FLOPS:                ", f20.3, " (GFLOPS)")') name, flops
  else
    write(*, '(a, " FLOPS:                ", a20)') name, "UNKNOWN"
  end if
  if (time > 0.0 .and. tput > 0.0) then
    write(*, '(a, " Throughput:           ", f20.3, " (GB/S)")') name, tput
  else
    write(*, '(a, " Throughput:           ", a20)') name, "UNKNOWN"
  end if
  if (time > 0.0 .and. etput > 0.0) then
    write(*, '(a, " Effective Throughput: ", f20.3, " (GB/S)")') name, etput
  else
    write(*, '(a, " Effective Throughput: ", a20)') name, "UNKNOWN"
  end if

  end if

end subroutine perf_print


subroutine perf_print_time(id, name)
  integer :: id
  character(*) :: name
  real(8) :: time

  time = perf_get_time(id, PERF_ROOT)
  write(*, '(a, f20.3, " [sec]")') name, time

end subroutine perf_print_time

#ifdef USE_MPI
subroutine perf_print_time_mpi(id, name)
  integer :: id
  character(*) :: name
  integer :: myrank, ierr
  real(8) :: time_root, time_min, time_max

  time_root = perf_get_time(id, PERF_ROOT)
  time_min  = perf_get_time(id, PERF_MIN)
  time_max  = perf_get_time(id, PERF_MAX)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  if (myrank == 0) then
    write(*, '(a, "(rank0)", f14.5, " / (min)", f14.5, " / (max)", f14.5)')  &
      name, time_root, time_min, time_max
  end if

end subroutine perf_print_time_mpi

subroutine perf_print_time_mpi_full(id, name)
  integer :: id
  character(*) :: name
  integer :: myrank, ierr
  real(8) :: time_root, time_min, time_max, time_ave, time_sd

  time_root = perf_get_time(id, PERF_ROOT)
  time_min  = perf_get_time(id, PERF_MIN)
  time_max  = perf_get_time(id, PERF_MAX)
  time_ave  = perf_get_time(id, PERF_AVE)
  time_sd   = perf_get_time(id, PERF_SD)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  if (myrank == 0) then
    write(*, '(a)') name
    write(*, '("    rank0: ", f14.5, " [sec]")') time_root
    write(*, '("    min:   ", f14.5, " [sec]")') time_min
    write(*, '("    max:   ", f14.5, " [sec]")') time_max
    write(*, '("    ave:   ", f14.5, " [sec]")') time_ave
    write(*, '("    sd:    ", f14.5, " [sec]")') time_sd
  end if

end subroutine perf_print_time_mpi_full

#endif

end module mod_perf
