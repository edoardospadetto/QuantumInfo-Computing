include "modules/debugqic.F90"
include "modules/matrixqic.F90"
include "modules/printer.F90"
include "modules/qstates.F90"
include "modules/isingop.F90"



!########################################
!MAIN program
! #####################################
program sch_td
     use printer
     use debugqic
     use quantumstatesqic
     use matrixqic
     use isingop
     implicit none



     double complex, dimension(:,:), allocatable :: test_system, test_systemsep, test_systemns
     double complex, dimension(:,:,:), allocatable :: test3
     real*8, dimension(:), allocatable :: energies
     double precision :: lam
     integer ::ii,jj,kk, teststate, testbody,ll
     integer, dimension(:), allocatable :: tokeep
     real :: start, start1, start2 , finish, finish1, finish2
     character(len = 1000):: fn, GG
     integer :: bodies(4)
     

         bodies = (/3,4,5,8/)
         do kk = 1, 4
             testbody = bodies(kk) ! Number of bodies
             print*, "running for number of bodies = ", testbody
             allocate(test_system(2**testbody,2**testbody))
             allocate(energies(2**testbody))


             do ii = 0, 100 ! 100 steps till lambda = 3
                 lam = dble(ii)*(3.0/100.00)
                 print*, lam
                 test_system = fullham(lam,testbody) ! Get the hamiltonian

                call eigz(test_system, size(test_system, dim = 1), energies) !Matrixqic module
                do jj = 1,5
                    write(fn,"(I5)") jj
                    write(GG,"(I2)") testbody
                    call wddo("./results"//trim(adjustl(GG))//"/data"//trim(adjustl(fn))//".txt",lam, energies(jj)/(testbody))
                end do
                
                !### ########Used to get info of eigenstates ####################
                !energies = norm2(test_system(:,1)%re)**2 + norm2(test_system(:,1)%im)**2 ! To get info of eigenstates
                !do jj = 1, size(energies)
                !    write(fn,"(I5.5)") ii+100
                !    call wddo("./results/data"//trim(adjustl(fn))//".txt",dble(jj), energies(jj))
                !end do
            end do
            print*, "      -   DONE!"
            deallocate(test_system)
            deallocate(energies)
          end do 



end program
