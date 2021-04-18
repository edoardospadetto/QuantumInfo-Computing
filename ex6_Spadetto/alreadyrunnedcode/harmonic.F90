include "modules/debugqic.F90"
include "modules/matrixqic.F90"

module quantumop
    use debugqic
        !variables of the simulation
        double precision :: xmin,xmax,e
        integer:: dim



    contains
        ! set variables of the simulation
        subroutine setenv(x1,x2,L)
            implicit none
            double precision :: x1,x2
            integer :: L
            xmin = x1
            xmax = x2
            dim  = L
            e = (x2-x1)/real(2*L+1)
        end subroutine setenv
        !###############################################
        ! compute kinetic energy therm
        function K_1d(cnst) result(K)
        !###############################################
            implicit none
            integer :: N,ii
            double  complex :: cnst
            double  complex :: K(2*dim+1,2*dim+1)
            K= (.0,.0)

            do ii = 2, 2*dim
                    K(ii,ii) =(-2.0,.0)
                    K(ii+1,ii)=(1.0,.0)
                    K(ii-1,ii)=(1.0,.0)
            end do
            K(1,1) =(-2.0,.0)
            K(2*dim+1,2*dim+1) =(-2.0,.0)
            K(2,1)=(1.0,.0)
            K(2*dim,2*dim+1)=(1.0,.0)

            K=cnst*K*(1/e**2)
        end function K_1d
        !###############################################
        !compute the potential energy therm
        function x2_1d(cnst) result(x2)
        !###############################################
            implicit none
            integer :: ii
            double complex :: x2(2*dim+1,2*dim+1), cnst
            double complex :: tmp
            x2 = (0.0,0.0)
            do ii = 1, 2*dim+1
                tmp%im=0.0
                tmp%re = e*(ii-1-dim)
                x2(ii,ii) = tmp**2
            end do
            x2=cnst*x2
        end function x2_1d
    
        !###############################################
        !compute and output eigenvectors real part and probabilyty density function
        subroutine eigvecout(x1,x2,hamiltonian,which)
        !###############################################
            use matrixqic
            implicit none
            double complex, dimension(:,:) :: hamiltonian
            character*20 :: numb
            integer :: which, L,ii
            double precision :: x1, x2,eps

            L = size(hamiltonian, 1)
            eps= (x2-x1)/real(L)
            
            write(numb,"(I3,'_',I4,'_',E8.2)") which,L,x2
            numb = trim(adjustl(numb))
            !compute prob
            hamiltonian(:, which) = hamiltonian(:, which)/sqrt(norm2(hamiltonian(:, which)%re)**2&
                                                        +norm2(hamiltonian(:, which)%im)**2)
            !output data
            do ii = 1,L
                call wddo(trim("./results/real")//trim(numb)//".txt" , x1 + (x2-x1)/real(L)*ii, (hamiltonian(ii,which)%re/sqrt(e)))
                call wddo(trim("./results/prob")//trim(numb)//".txt" , x1 + (x2-x1)/real(L)*ii, &
                        (realpart(hamiltonian(ii,which))**2+imagpart(hamiltonian(ii,which))**2)/e )
            end do
        end subroutine eigvecout
    
        !###############################################
        ! Compute and output stats about eigvals
        subroutine deltaeigvalhs(x1,x2,eigv,cnst)
        !###############################################
            use matrixqic
            implicit none
            double precision, dimension(:) :: eigv
            integer :: L,ii
            double precision :: x1, x2,eps, dsum, cnst
            character*20 :: numb, nb

            L = size(eigv, 1)
            eps= (x2-x1)/real(L)

            write(numb,"(I4,'_',E8.2)") L,abs(x1)
            numb = trim(adjustl(numb))
            write(nb,"('_',E8.2)") abs(x1)
            nb = trim(adjustl(nb))
            dsum =0.
            !compute s300
            do ii = 1, 300
                dsum= dsum + (eigv(ii)-cnst*(ii-0.5))**2
            end do
                !output s300
                call widdo("./results/score"//trim(nb)//".txt", L, eps, dsum/(300.0)**2)
                dsum =0
                !compute stot (not used in the report)
                do ii = 1, size(eigv)

                    dsum= dsum + (eigv(ii)-cnst*(ii-0.5))**2
                    call widdo("./results/deltaeig"//trim(numb)//".txt", ii, eigv(ii), (eigv(ii)-cnst*(ii-0.5))**2)

                end do
                !compute and output s20
                dsum =0
                do ii = 1, 20

                    dsum= dsum + (eigv(ii)-cnst*(ii-0.5))**2

                end do
                call widdo("./results/score20_"//trim(nb)//".txt", L, eps, dsum/(10.0**2))
        end subroutine deltaeigvalhs

end module quantumop


!###############################################
!used for debbugging
module test_mod
!###############################################
    contains
        subroutine test0(op, vec, eigvals, eps, cnst)
            implicit none
            double complex, dimension(:,:) :: op
            double complex, dimension(:,:):: vec
            double precision, dimension(:):: eigvals
            double precision:: eps, cnst

            integer*2 :: ii
            double complex, dimension(:), allocatable:: cgvec

            allocate (cgvec(size(eigvals)))
            print*,eps


            do ii = 1, 100
                cgvec=conjg(vec(:,ii))
                print*, eigvals(ii)
                print*, cnst*(ii-0.5)
                print*, dot_product(cgvec, matmul(op,vec(:,ii)) / dot_product(cgvec,vec(:,ii)))
            end do


            deallocate (cgvec)

        end subroutine test0


end module test_mod




!########################################################################################
!The program
program harmonic
use test_mod
use  debugqic
use matrixqic
use quantumop
implicit none
    integer*4 :: L,ii,jj
    double complex, dimension(:,:), allocatable :: hamiltonian,ham
    double precision, dimension(:), allocatable :: eigf
    double precision :: x1 =-10.0, x2 =10.0, myw, hm , ko, h=1. , m=1. , w = 20., eps, cnst0,dd(5)
    double complex :: ek = 0., harm=0.
    real :: start, finish
    dd = (/3.0,4.0,5.0,7.0,10.0/)


    do ii = 1,10  !iterate over L
        do jj = 1,5 !loop domain
            x1 =- dd(jj)
            x2 = dd(jj)
            L = 70*ii+100

            eps= (x2-x1)/(2*L+1)
            call CPU_TIME(start)
            allocate(hamiltonian(2*L+1,2*L+1))
            allocate(eigf(2*L+1))
            allocate(ham(2*L+1,2*L+1))

            hm=h*0.5*m !set constants
            ko=m*w*w*0.5
            myw=w*h

            harm%re = ko
            ek%re = -hm

            call setenv(x1,x2,L)
            hamiltonian = +K_1d(ek)+x2_1d(harm)  !here we get the hamiltonian
            ham = hamiltonian ! copy hamiltonia for debug

            !compute eigenvalues and eigenvectors
            call eigz(hamiltonian,2*L+1,eigf)


            call CPU_TIME(finish)
            !print eigenvectors, stats of eigenvalues and time
            call eigvecout(x1,x2,hamiltonian,1)
            call eigvecout(x1,x2,hamiltonian,2)
            call eigvecout(x1,x2,hamiltonian,3)
            call eigvecout(x1,x2,hamiltonian,4)
            call deltaeigvalhs(x1,x2,eigf,myw)
            call wfio("./results/time.txt", 2*L+1, finish-start)

            deallocate(hamiltonian)
            deallocate(eigf)
            deallocate(ham)
        end do
    end do

end program
