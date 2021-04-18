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



     double complex, dimension(:,:), allocatable :: test_system, test_systemsep, test_systemns, Prf, Prb, rdHam, opn,opn1
     double complex, dimension(:,:,:), allocatable :: test3
     double complex :: paulix(2,2)
     real*8, dimension(:), allocatable :: energies
     double precision :: lam
     integer ::ii,jj,kk, teststate, testbody,ll,ss, N
     integer, dimension(:), allocatable :: tokeep
     real :: start, start1, start2 , finish, finish1, finish2
     character(len = 1000):: fn, GG
     integer :: bodies(3), iterations(3)
            
            bodies = (/2,3,4/)
            iterations = (/5,10,20/)
            
            do ss = 1, 3
                do ll = 1, 3

                    N = bodies(ll)
                    testbody = 2*N !2N
                    lam = 0.5
                     allocate(test_system(2**testbody,2**testbody))

                     allocate(Prf(2**testbody,2**N))
                     allocate(Prb(2**N,2**testbody))
                     allocate(rdHam(2**N,2**N))
                     allocate(opn(2**N,2**N))
                     allocate(opn1(2**N,2**N))
                     allocate(energies(2**testbody))

                       
                        
                    do kk = 1, 100
                        lam = dble(kk)*(3.0/100.00)

                        paulix = 0.0
                        paulix(1,1)%re=0.0
                        paulix(1,2)%re=1.0
                        paulix(2,1)%re=1.0
                        paulix(2,2)%re=0.0
                        
                        
                        opn = getbigmat(paulix,testbody,testbody)
                        opn1 = getbigmat(paulix,1,testbody)

                        print*, lam

                        !Compute the Hamiltonian for N particles
                        test_system = fullham(lam,testbody) ! Get the hamiltonian
                        !Diagonalize 
                        print*, "iterations = " ,ss 
                        print*, "Nop = " , N 
                        do jj = 1, iterations(ss)
                            call eigz(test_system, size(test_system, dim = 1), energies) !Matrixqic module          
                            
                            Prf = test_system(:,1:2**N) !Got projector
                            call invz(test_system, size(test_system, dim = 1) )
                            Prb = test_system(1:2**N,:)      !Got "inverse" Projector

                            test_system = cmplx( dble(0.0),dble(0.0) )
                            do ii = 1, 2**N 
                                rdHam(ii,ii) = energies(ii) ! Got reduced hamiltonian
                            end do 

                            

                            opn=  matmul(Prb, matmul(opn,Prf))
                            opn1= matmul(Prb, matmul(opn1,Prf))
                            print*, "here!"
                            opn= getbigmat(opn,2,2) 
                            opn1= getbigmat(opn1,1,2) 
                            test_system = getbigmat(rdHam,1,2) + getbigmat(rdHam, 2,2) + matmul(opn,opn1)
                            
                        end do 
                        write(GG,"(I5)") iterations(ss)
                        write(fn,"(I5)") N
                        call wddo("./results/data"//trim(adjustl(fn))//"_"//trim(adjustl(GG))//".txt",lam, & 
                                     energies(1)/(testbody*2**(iterations(ss)-1)) )
                        
                        
                    
                    end do 

                  
                    print*, "      -   DONE!", N 
                    deallocate(test_system)
                    deallocate(energies)
                    deallocate(Prf)
                    deallocate(Prb)
                    deallocate(rdHam)
                    deallocate(opn)
                    deallocate(opn1) 
                end do 
            end do 
      

          



end program
