include "./matrixqic.f90"
include "./mat_help.f90"
program test
    use mmatrixqic
    use zero
    implicit none

    type(matrixqic):: test_matrix, adj_test_matrix, tmp
    double complex, dimension(:,:), allocatable :: tmpmatrix
    integer :: N, pow, ierr
    character*40 :: stats_file, file1
    real:: start, finish, tracet, adjt
    print*, ""
    print*, "--- Excercise 2 --- Edoardo Spadetto"
    print*,""

    N=3

    allocate(tmpmatrix(N,N))
    tmpmatrix =  gen_rnd_matrc(N,N)
    !Here we init the empty type matrixqic with dimension N
    call Init(N,test_matrix)
    !and we see that it is empty
    print*, "Empty Matrix"
    call printmatrixc(test_matrix%el)
    print*, ""

    !we set the elements of the matrix equal to the tmpmatrix
    test_matrix%el = tmpmatrix
    print*, "Matrix"
    call printmatrixc(test_matrix%el)
    print*, ""

    !TRACE
    !Now the trace is not computed
    print*, "Empty Trace"
    print*, test_matrix%trace
    print*, ""
    !We can use the overloaded operator to compute, store and get the value
    print*, "Trace operator"
    print*, .tr.test_matrix
    print*, ""
    !or we can compute it and store it in the trace field
    call ctrace(test_matrix)
    !and then explicitely ask for it
    print*, "Trace explicit call"
    print*, test_matrix%trace
    print*, ""


    !ADJOINT
    !We can use the overloaded operator to get the adjoint
    tmp =  .adj.test_matrix
    !print it
    print*, "Adjoint operator"
    call printmatrixc(tmp%el)
    print*, ""
    !Or explicitely set a matrixqic to the adjoint of the matrix
    call cgadj(test_matrix, adj_test_matrix)
    print*, "Adjoint explicit call"
    call printmatrixc(adj_test_matrix%el)
    print*, ""

    !PRINT ON FILE
    file1 = "outex2.txt"
    call ofile(adj_test_matrix,file1)

    !deallocate used memory
    deallocate(tmpmatrix)
    deallocate(test_matrix%el)
    deallocate(adj_test_matrix%el)

    !RUN Performance tests on the computation time
    print*, "Performances"

    N=3

    stats_file = "ex2stats.txt"

    !open(unit = 1, file=stats_file, status = "old")

    do pow = 1, 8
        allocate(tmpmatrix(N**pow,N**pow),stat=ierr)
        if(ierr .eq. 0 ) then
            tmpmatrix =  gen_rnd_matrc(N**pow,N**pow)
            call Init(N**pow,test_matrix)
            test_matrix%el = tmpmatrix
            CALL CPU_TIME(start)
            call ctrace(test_matrix)
            CALL CPU_TIME(finish)
            tracet = finish-start
            CALL CPU_TIME(start)
            call cgadj(test_matrix, adj_test_matrix)
            CALL CPU_TIME(finish)
            adjt= finish-start
            write(*,*) N**(pow*2), tracet, adjt
            !write(1,*) N**(pow*2), tracet, adjt


            deallocate(tmpmatrix)
            deallocate(test_matrix%el)
            deallocate(adj_test_matrix%el)
        else
            print*, ierr
        end if

    end do
    !close(1)



end program test
