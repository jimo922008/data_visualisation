module timing_module

    implicit none
    private
    public :: start_timer, stop_timer, elapsed_time

    integer :: start_clock, end_clock, count_rate
    real(8) :: elapsed_seconds

contains

    subroutine start_timer()
        call system_clock(count_rate=count_rate)
        call system_clock(start_clock)
    end subroutine start_timer

    subroutine stop_timer()
        call system_clock(end_clock)
        elapsed_seconds = real(end_clock - start_clock, 8) / real(count_rate, 8)
    end subroutine stop_timer

    function elapsed_time() result(time)
        real(8) :: time
        time = elapsed_seconds
    end function elapsed_time

end module timing_module