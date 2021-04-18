include "modules/debugqic.F90"
include "modules/matrixqic.F90"
include "modules/printer.F90"
include "modules/qstates.F90"



!########################################
!MAIN program
! #####################################
program sch_td
     use printer
     use debugqic
     use quantumstatesqic
     implicit none


     type(compositewf) :: sys,syscheck
     type(dnsmat) :: density_matrix, mat2, mat2check
     double complex, dimension(:,:), allocatable :: test_system, test_systemsep, test_systemns
     double precision :: tessst
     integer ::ii,jj,kk, teststate, testbody,ll
     integer, dimension(:), allocatable :: tokeep
     real :: start, start1, start2 , finish, finish1, finish2

     !1) Check timing of wavefunction generation
     write(*,*) "######## Test Initialization of Systems ################"
     teststate = 2
     do ii = 2,25
         allocate(test_systemsep(teststate,ii))
         allocate(test_systemns(teststate**ii,1))
         write(*,*)"    -allocated"
         call cpu_time(start1)
         test_systemsep = newsystemsep(teststate, ii)
         call cpu_time(finish1)
         write(*,*)"    -donesep"
         call cpu_time(start2)
         test_systemns = newsystemnonseparable(teststate, ii)
         call cpu_time(finish2)
         call wdddo("./results/generation2states.txt", dble(ii) , dble(finish1-start1), dble(finish2-start2))
         deallocate(test_systemsep)
         deallocate(test_systemns)
     end do
     write(*,*) "######## Done! ################"
     !2) Test Reliability
     !2.1) Pure states
     write(*,*) "######## Test Reliability: Pure states ################"
     teststate = 3
     testbody = 3
     allocate(tokeep(1))
     tokeep(1) = 1
     allocate(test_system(teststate,testbody))
     test_system = newsystemsep(size(test_system,dim = 1) , size(test_system,dim = 2))
     write(*,*) "      -test state builded"
     sys =  MergeWF(test_system)
     write(*,*) "      -Got total wavefunction"
     syscheck = MergeWF(test_system(:,tokeep ))
     write(*,*) "      -Got total wavefunction for test"
     mat2check =  BDMmb(syscheck)
     write(*,*) "      -Got density matrix for test"
     density_matrix = BDMmb(sys)
     write(*,*) "      -Got density matrix"
     mat2 = marginalizeDMMB2(density_matrix , tokeep)
     write(*,*) "      -Reduced"
     write(*,*) "------- Control density Matrix --------"
     call pzm(mat2check%mat)
     write(*,*) "------- Obtained density Matrix --------"
     call pzm(mat2%mat)
     deallocate (test_system)
     deallocate(tokeep)
     write(*,*) "######## Done! ################"

     !3) Test timing of reduction
     write(*,*) "######## Test Reduction Timing ################"

     !testbody = 6
     do jj = 2,5
     do ii = 2,4
     do kk = 1,ii-1
             testbody = ii
             teststate = jj
             allocate(test_system(teststate,testbody))
             allocate(tokeep(kk))
             do ll = 1,size(tokeep)
                 tokeep(ll) = ll
             end do
             sys =  MergeWF(test_system)
             density_matrix = BDMmb(sys)
             call cpu_time(start)
             mat2 = marginalizeDMMB2(density_matrix , tokeep)
             call cpu_time(finish)
              call wddddo("./results/red_bodies2.txt", dble(testbody), dble(teststate),dble(kk) , dble(finish-start))
             deallocate(tokeep)
             deallocate(test_system)
    end do
    end do
    end do



end program
!first 2 states 13 bodies
!second 3 states 9 bodies
