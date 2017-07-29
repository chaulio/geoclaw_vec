module papi_module
    implicit none
#ifndef _PAPI_MODULE_
#define _PAPI_MODULE_
    include 'f90papi.h'
#endif
    private
    ! Floating point types    
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.d0)

    ! FLOPS counting related variables
    integer,   parameter            :: numevents = 1
    integer,   dimension(numevents) :: events = (/ PAPI_DP_OPS /)
    integer*8, dimension(numevents) :: values ! To store event results
    integer*8                       :: flpops ! Absolute # of FP ops
    integer*8                       :: maxflpops = 0, minflpops = huge(1_8)
    integer                         :: eventset = PAPI_NULL
    real(kind=DP)                   :: mflops
    real(kind=DP)                   :: avg_mflops = 0.d0
    
    ! Timing
    integer*8                       :: clock_start, clock_end
    real(kind=DP)                   :: total_time = 0.d0, acc_time = 0.d0
    real(kind=DP)                   :: avg_rims = 0.d0 ! Avg. Riemann solves/second

    ! Other
    integer                         :: maxrs = 0, minrs = huge(1_4) ! Min/Max 1d domain length
    logical                         :: riemannstats = .false.
    
    ! PAPI - general
    integer                         :: ret = PAPI_VER_CURRENT
    integer                         :: calls = 0, zerotimes = 0, zeroflops = 0

    save 
    public :: papi_init, papi_start, papi_stop, papi_summary
contains
    subroutine papi_init()
        implicit none

        ! Initialize PAPI library
        call PAPIF_library_init(ret)
        if (ret .lt. 0) then
            print *, "FATAL: An error occured while initializing!", ret
            call exit(1)
        end if

        call PAPIF_create_eventset(eventset, ret)
        ! Add events to be measured
        call PAPIF_add_event(eventset, events, ret)
    end subroutine

    subroutine papi_start()
        calls = calls + 1
        ! ======= TIMING =======
        !call PAPIF_get_virt_usec(clock_start)
        call PAPIF_get_real_nsec(clock_start)
        ! ======= TIMING =======
        
        call PAPIF_start(eventset, ret)
    end subroutine

    subroutine papi_stop(riemann_solves)
        integer, intent(in), optional :: riemann_solves

        call PAPIF_stop(eventset, values, ret)
        
        ! ======= TIMING =======
        !call PAPIF_get_virt_usec(clock_end)
        call PAPIF_get_real_nsec(clock_end)
        ! ======= TIMING =======
        ! Count number of zero times (=> insufficient timer resultion)
        if (clock_end - clock_start .eq. 0) zerotimes = zerotimes + 1
        
        ! Total time in microseconds
        total_time = real(clock_end - clock_start, kind=DP)
        flpops = values(1)
        if (flpops .lt. 1) then
            zeroflops = zeroflops + 1
            mflops = 0.d0
        else
            ! Note: if ns timing is used this is actually GFLOPS not MFLOPS
            mflops = real(flpops, kind=DP)/(total_time)
        endif
        
        ! Total time in seconds
        total_time = total_time * 1d-9
        
        ! Riemann solves per second
        if (present(riemann_solves)) then
            riemannstats = .true.
            if (riemann_solves .eq. 0) write(*,*) "Warning: Domain length is 0!"
            avg_rims = (avg_rims*(calls-1) &
                + real(riemann_solves,kind=DP)/total_time) &
                / real(calls, kind=DP)
            maxflpops = max(maxflpops, flpops/riemann_solves)
            minflpops = min(minflpops, flpops/riemann_solves)
            minrs = min(minrs, riemann_solves)
            maxrs = max(maxrs, riemann_solves)
        endif
        ! Total time in seconds
        acc_time = acc_time + total_time

        ! Determine FLOP-related values
        avg_mflops = (avg_mflops*(calls-1)+mflops) / real(calls, kind=DP)
    end subroutine

    subroutine papi_summary()
        !write(*,"(a25,i10)") "Number of FLOPs:", flpops
        !write(*,'(a,i5,f14.8)') "=> TIMER TIME [s]:", N, total_time
        write(*,'(a,i14)')  "=> # OF CALLS:", calls
        write(*,'(a,f14.8)')"=> ACCUM. TIME [s]:", acc_time
        write(*,'(a,f14.8)')"=> PAPI [G|M]FLOPS (AVG):", avg_mflops
        write(*,'(a,2i8)')  "=> # ZERO TIMES/FP OPS:", zerotimes, zeroflops
        if (riemannstats) then
            write(*,*)          "======= RPN2 STATS ========"
            write(*,'(a,f14.3)')"=> RIEMANN SOLVES/SECOND (AVG):", avg_rims
            write(*,'(a,2i8)')  "=> MIN/MAX FP OPS PER RIEMANN SOLVE (AVG):", &
                minflpops, maxflpops 
            write(*,'(a,2i8)')  "=> MIN/MAX RIEMANN SOLVES PER CALL:", minrs, maxrs
        endif
        !write(42,'(i5,f14.8)') mflops
    end subroutine
end module
