module statistics
    integer, private, parameter :: DP = kind(1.d0)
    integer, save :: rpn2_total_calls = 0
    real(kind=DP) :: rpn2_time_start = 0.d0, rpn2_time_end = 0.d0
    real(kind=DP) :: rpn2_total_time = 0.d0
contains
    subroutine print_stats()
        write(*,*) "Total number of rpn2 calls:", rpn2_total_calls
        write(*,*) "Total time spent in rpn2():", rpn2_total_time
    end subroutine

    subroutine rpn2_start_timer()
        rpn2_total_calls = rpn2_total_calls + 1
        call cpu_time(rpn2_time_start)
    end subroutine

    subroutine rpn2_stop_timer()
        call cpu_time(rpn2_time_end)
        rpn2_total_time = rpn2_total_time + rpn2_time_end - rpn2_time_start
    end subroutine
end module statistics
