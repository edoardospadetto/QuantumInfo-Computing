

module debugqic
contains
    subroutine breakifn(d_string, cndtion, debug)
        implicit none
        character(len = *) :: d_string
        logical :: cndtion, debug
        ! if debug is enabled continue otherwise noop
        if (debug) then

            !if error occurs stop the program
            if (.not. cndtion) then
                write (*,*) d_string 
                write(*,*) "Execution terminated"
                STOP
            end if
        end if
    end subroutine breakifn

end module debugqic
