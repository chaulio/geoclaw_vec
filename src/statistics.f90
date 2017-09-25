module statistics
!$  use omp_lib    

    integer, private, parameter :: DP = kind(1.d0)
    integer, allocatable       :: rpn2_total_calls(:)
    real(kind=DP), allocatable :: rpn2_total_time(:)
    integer, allocatable       :: rpt2_total_calls(:)
    real(kind=DP), allocatable :: rpt2_total_time(:)
contains

#   if !defined(_OPENMP)
        ! if openmp is not enabled, we need these functions to exist!
       
        function omp_get_thread_num() result(i_thread)
            integer :: i_thread
            i_thread = 0
        end function

        function omp_get_max_threads() result(i_threads)
            integer :: i_threads
            i_threads = 1
        end function
#   endif

    ! Creates arrays - each thread will have its own timers and counters
    subroutine init_stats()
        integer     :: num_threads
        
        num_threads = omp_get_max_threads()

        allocate(rpn2_total_calls(num_threads))
        allocate(rpn2_total_time(num_threads))
        allocate(rpt2_total_calls(num_threads))
        allocate(rpt2_total_time(num_threads))
        
        rpn2_total_calls = 0
        rpn2_total_time = 0.0
        rpt2_total_calls = 0
        rpt2_total_time = 0.0
    end subroutine
    
    function get_thread_num() result(thread_num)
        integer :: thread_num
        thread_num = omp_get_thread_num()
    end function

    function get_wtime() result(wtime)
        double precision :: wtime

#       if defined(__OPENMP)
            wtime = omp_get_wtime()
#       else
            integer(kind = selected_int_kind(16))   :: counts, count_rate

            call system_clock(counts, count_rate)

            wtime = dble(counts) / dble(count_rate)
#       endif
    end function

    subroutine print_stats()
        integer          :: i, rpn2_calls, rpt2_calls
        double precision :: rpn2_time, rpt2_time
        
        ! sum up all counters/timers
        rpn2_calls = 0
        rpt2_calls = 0
        rpn2_time = 0.0
        rpt2_time = 0.0
        do i=1,omp_get_max_threads()
            rpn2_calls = rpn2_calls + rpn2_total_calls(i)
            rpt2_calls = rpt2_calls + rpt2_total_calls(i)
            rpn2_time = rpn2_time + rpn2_total_time(i)
            rpt2_time = rpt2_time + rpt2_total_time(i)
        end do
        
        ! get the average time
        rpn2_time = rpn2_time / omp_get_max_threads()
        rpt2_time = rpt2_time / omp_get_max_threads()
        
        write(*,*) "rpn2 (calls/time):", rpn2_calls, rpn2_time, "(avg. thread time)"
        write(*,*) "rpt2 (calls/time):", rpt2_calls, rpt2_time, "(avg. thread time)"
    end subroutine

    subroutine rpn2_start_timer()
        integer     :: idx
        idx = get_thread_num() + 1
        
        rpn2_total_calls(idx) = rpn2_total_calls(idx) + 1
        rpn2_total_time(idx) = rpn2_total_time(idx) - get_wtime()
    end subroutine

    subroutine rpn2_stop_timer()
        integer     :: idx
        idx = get_thread_num() + 1

        rpn2_total_time(idx) = rpn2_total_time(idx) + get_wtime()
    end subroutine
    
    subroutine rpt2_start_timer()
        integer     :: idx
        idx = get_thread_num() + 1

        rpt2_total_calls(idx) = rpt2_total_calls(idx) + 1
        rpt2_total_time(idx) = rpt2_total_time(idx) - get_wtime()
    end subroutine

    subroutine rpt2_stop_timer()
        integer     :: idx
        idx = get_thread_num() + 1

        rpt2_total_time(idx) = rpt2_total_time(idx) + get_wtime()
    end subroutine
end module statistics
