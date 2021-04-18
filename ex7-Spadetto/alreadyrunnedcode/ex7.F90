include "debugqic.F90"
include "matrixqic.F90"
include "printer.F90"
include "quantumop.F90"
include "tdschrodinger.F90"


!###############################################
!used for debbugging
module test_mod
!###############################################
    contains
!Computes eigenstate of a moving harmonic oscillator.
!not the solution of the problem.
        subroutine realeigf(x1,x2,m,deltat,N,vel,idx,w,val)
            use printer
            implicit none
            double precision:: vale,val
            integer, intent(IN):: N,idx
            integer ::ii,jj
            double complex, dimension(N) :: result
            double precision , intent(IN):: w,x1, x2, deltat,vel, m
            double precision ::Aim, Areal, Cim, Creal, Bim, Breal,eigzero,x,t

            t = deltat*(idx)
            do ii = 1, N

                x = x1 +(x2-x1)*ii/dble(N)

                if((x-t*vel)**2 < 2) then
                    eigzero = exp(-0.5*w*(x-t*vel)**2)
                else
                    eigzero =dble(0.0)
                end if
                Areal = cos(m*x*vel - 0.5*m*t*(vel**2))
                !print*, t
                Breal = cos(w*t)
                Aim=sin(m*x*vel - 0.5*m*t*(vel**2))
                Bim=sin(w*t)
                Creal = Areal*Breal - Aim*Bim
                Cim = Aim*Breal + Areal*Bim

                result(ii)%re = eigzero*Creal
                result(ii)%im = eigzero*Cim

            end do
                call ozgrid(x1-vel*t,x2-vel*t,result,idx, "tth")
                vale = 0.D0
                do jj = 1, size(result)
                    vale = vale + (x1+jj*((x2-x1)/N))*(result(jj)%re**2+result(jj)%im**2)

                end do
                vale = vale / (norm2(result%re)**2+norm2(result%im**2))
                write(*,*) vale, vel*t
                call wddo("data7.txt",vale*(w/vel), t)

        end subroutine


end module test_mod

!########################################
!MAIN program
! #####################################
program sch_td
     use printer
     use debugqic
     use quantumop
     use time_dependent_schrodinger
     use test_mod

     implicit none

     integer :: L,ii, timesteps,jj
     double precision :: x1, x2, h , m , w, dt, velocity, now, value, val
     double precision :: eps
     double complex, dimension(:), allocatable :: ezero, evolved
     h=1.
     m=1.
     w=10.
     x1 = -5.D0
     x2 = 5.D0
     L = 600
     velocity = 0.5
     dt = 1e-4
     timesteps = int(1/(velocity*dt))
     write(*,*) "######### TIME EVOLUTION ###########"
     write(*,*) "steps : ", timesteps
     write(*,*) "T :" , dble(1)/velocity
     write(*,*) "dt :", dt
     eps = (x2-x1)/(2*dble(L)+1)
    allocate(ezero(2*L+1))
    allocate(evolved(2*L+1))
    call setenv(x1,x2,L)
    ezero = init(x1,x2,L,w,h,w)
    call ozgrid(x1,x2,ezero,0,"obt")
    evolved = ezero

    !evolve the function, and compute average, and std deviation.
    do ii = 1, timesteps
        now = dble(ii)*dt
        value = 0.0
        val = 0.0
        evolved = update(evolved,velocity,x1,x2,L,w,now-dt,dt,m,h)
        !average
        do jj = 1 , size(evolved)
            value = value +(eps*dble(jj-1))*(evolved(jj)*CONJG(evolved(jj)))
        end do
        !std deviation
        do jj = 1 , size(evolved)
            val = val +((eps*dble(jj-1)-value))**2*(evolved(jj)*CONJG(evolved(jj)))
        end do
            value = value+x1
        call wdddo("info2.txt", dble(ii)*dt, (value), sqrt(val))
        if(mod(ii,int(timesteps/100)).eq.0) then
                call ozgrid(x1,x2,evolved,ii,"obt")
            !print*, ii
        end if
    end do
    print*, "######## DONE! ########"


end program
