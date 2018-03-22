module mod_maprof

  interface maprof_set_ops
    module procedure maprof_set_ops_i4
  end interface maprof_set_ops

  interface

    subroutine maprof_output()
    end subroutine maprof_output

  end interface

contains

subroutine maprof_set_ops_i4
end subroutine maprof_set_ops_i4


end module mod_maprof
