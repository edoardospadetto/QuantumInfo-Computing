include "modules/debugqic.F90"
include "modules/matrixqic.F90"
include "modules/printer.F90"
include "modules/qstates.F90"
include "modules/isingop.F90"



!########################################
!MAIN program
! #####################################

module exercise10 
    contains 
    subroutine realspacerenormalization
     use printer
     use debugqic
     use quantumstatesqic
     use matrixqic
     use isingop


     double complex, dimension(:,:), allocatable :: test_system, test_systemsep, test_systemns, Prf, Prb, rdHam, opn,opn1
     double complex, dimension(:,:,:), allocatable :: test3
     double complex :: paulix(2,2)
     real*8, dimension(:), allocatable :: energies
     double precision :: lam
     integer ::ii,jj,kk, teststate, testbody,ll, N
     real :: start, start1, start2 , finish, finish1, finish2
     character(len = 1000):: fn, GG
     integer :: bodies(4)
     
            N = 3
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
                do jj = 1, 10

                    call eigz(test_system, size(test_system, dim = 1), energies) !Matrixqic module          
                    Prf = test_system(:,1:2**N) !Got projector
                    call invz(test_system, size(test_system, dim = 1) )
                    Prb = test_system(1:2**N,:)      !Got "inverse" Projector

                    test_system = cmplx( dble(0.0),dble(0.0) )
                    
                    do ii = 1, 2**N 
                        rdHam(ii,ii) = energies(ii) ! Got reduced hamiltonian
                    end do 

                    print*, size(prf, dim = 1), size(prf, dim = 2), size(prb, dim = 1), size(prb, dim = 2), size(opn, dim = 1)

                    opn=  matmul(Prb, matmul(opn,Prf))
                    opn1= matmul(Prb, matmul(opn1,Prf))
                    print*, "here!"
                    opn= getbigmat(opn,2,2) 
                    opn1= getbigmat(opn1,1,2) 
                    test_system = getbigmat(rdHam,1,2) + getbigmat(rdHam, 2,2) + matmul(opn,opn1)
                    
                end do 

                do jj = 1,5
                    write(fn,"(I5)") jj
                    call wddo("./results/data"//trim(adjustl(fn))//".txt",lam, energies(jj)/(testbody*2**5))
                end do


                
            end do 

            print*, "      -   DONE!"
            deallocate(test_system)
            deallocate(energies)

    end subroutine

    subroutine reduceddensitymatrixrg
     use printer
     use debugqic
     use quantumstatesqic
     use matrixqic
     use isingop



     double complex, dimension(:,:), allocatable :: test_system,Prf, Prb, rdHam, opn1, opn4, bigdns
     double complex, dimension(:,:), allocatable :: mvham
     type(dnsmat) :: dnsbig,temp
     double complex ,dimension(:,:), allocatable :: dnstiny
     double complex :: paulix(2,2), pauliz(2,2)
     real*8, dimension(:), allocatable :: energies, eigvals
     double precision :: lam
     integer ::ii,jj,kk, teststate, testbody,ll, N
     integer, dimension(:), allocatable :: tokeep
     real :: start, start1, start2 , finish, finish1, finish2
     character(len = 1000):: fn, GG
     integer :: bodies(4)

            N = 1
            testbody = (N+1)*2 !2N
            lam = 0.5
            allocate(test_system(2**testbody,2**testbody))  !Htot = H1+H2+H3+H4+H12+H23+H34
            allocate(Prf(2**(N+1),2**N))                 !Projector 
            allocate(Prb(2**N,2**(N+1)))                 !Projector inverse
            allocate(tokeep(N+1))
            allocate(dnstiny(2**N,2**N))
            allocate(mvHam(2**(N+1), 2**(N+1)))          !mvham = H1 +H2 +H12
            allocate(rdHam(2**N,2**N))                   !H'1 = P-1 (mvham)  P after projection
            allocate(opn1(2**N,2**N))                   ! interaction  projectorinv x 1L(2**n) x sigmax(2) x projector before
            allocate(opn4(2**(N),2**(N)))                    ! interaction  projectorinv x 1L(2**n) x sigmax(2) x projector before
            allocate(energies(2**testbody))
            allocate(eigvals(2**N))

            do ii = 1, N+1
                tokeep(ii) = ii
            end do
               
                
            do kk = 1, 100
                 lam = dble(kk)*(3.0/100.00)

                 !set pauli matrix
                 paulix = 0.0
                 paulix(1,1)%re=0.0
                 paulix(1,2)%re=1.0
                 paulix(2,1)%re=1.0
                 paulix(2,2)%re=0.0

                pauliz = 0.0
                pauliz(1,1)%re=1.0
                pauliz(1,2)%re=0.0
                pauliz(2,1)%re=0.0
                pauliz(2,2)%re=-1.0
    
                
                !opn first round 

                

                print*, lam

            !Compute the Hamiltonian for 2*N+1 particles
                test_system = fullham(lam,testbody) ! Get the hamiltonian full for all 4 parts of the system for now 2*(N+1) particles
                   !First get the new hamiltonian for half of the total number of particles
                    mvham = fullham(lam,N+1)
            

                do jj = 1, 10
                  
                    !Diagonalize 
                    call eigz(test_system, size(test_system, dim = 1), energies) !Matrixqic module  
                    !Get density matrix
                    dnsbig = BDM(   (RESHAPE(test_system(:,1),(/2**(2*N+2), 1/))), (/dble(1.0)/)) 
                    dnsbig%nobodies = 2*N+2
                    dnsbig%nostatesxb = 2
                    !Get reduceddensitymatrix 
                    temp = marginalizeDMMB2(dnsbig, tokeep)
                    dnstiny = temp%mat !dimension should be 2**(N+1)
                    print*, "here0"
                    !Diagonalize density matrix
                    call eigz(dnstiny, size(dnstiny, dim = 1), eigvals)  
                    print*, "here0"
                    ! Now i have the projector  
                    !Get forward and backward projector
                    Prf = dnstiny(:,1:2**N) !Got projector
                    print*, "here21"
                    call invz(dnstiny, size(dnstiny, dim = 1) )
                     print*, "here22"
                    Prb = dnstiny(1:2**N,:)      !Got "inverse" Projector

                    dnstiny = cmplx( dble(0.0),dble(0.0) )

                    ! Condensate the first N+1 particles on 2**N dimension
                            !     print*, size(prf, dim = 1), size(prf, dim = 2), size(prb, dim = 1), size(prb, dim = 2), size(opn, dim = 1)
                    rdHam = matmul(Prb, matmul(mvham,Prf))
                    print*, "here23"
                    print*, size(opn1,dim=1), size(getbigmatdiff(2**(N),paulix,1), dim = 1) 
                    print*, "qui?"
                    opn1 = matmul(Prb, matmul(getbigmatdiff(2**(N),paulix,1),Prf)) !interaction term for the next iteration
                    print*,"this"
                    opn4 = matmul(Prb, matmul(getbigmatdiff(1,paulix,2**(N)),Prf))
                    print*, "Here"
                    !Add a particle 
                    mvham = getbigmatdiff(1,rdHam,2) + getbigmatdiff(2**N, pauliz, 1) &
                            + matmul( getbigmatdiff(1,opn1,2) , getbigmatdiff(2**N,paulix,1) ) 
                    !Full Hamiltonian

                    test_system = getbigmatdiff(1,rdHam,2**(N+2))+getbigmatdiff(2**N, pauliz, 2**(N+1))+ & 
                                getbigmatdiff(2**(N+2),rdHam,1)+getbigmatdiff(2**(N+1), pauliz, 2**(N)) &  
                                 +matmul( getbigmatdiff(1,opn1,2**(N+2)) , getbigmatdiff(2**N,paulix,2**(N+1)) ) + & ! Interaction H12 
                                 matmul( getbigmatdiff(2**N,paulix,2**(N+1)) , getbigmatdiff(2**(N+1),paulix,2**(N)) ) + & ! interaction H23
                                 matmul( getbigmatdiff(2**(N+1),paulix,2**(N)), getbigmatdiff(2**(N+2),opn1,1)) !interaction H34 

                
                end do 

                do jj = 1,5
                    write(fn,"(I5)") jj
                    call wddo("./results/datardm"//trim(adjustl(fn))//".txt",lam, energies(jj)/(testbody*2**5))
                end do


                
            end do 

            print*, "      -   DONE!"
            deallocate(test_system)
            deallocate(energies)


    end subroutine 

end module 

program sch_td
   
     use exercise10
     implicit none

     call reduceddensitymatrixrg




            
             




        
          



end program
