
include "matmulqic.f90"

program ex3

    use debugqic
    use matmulqic

    IMPLICIT NONE
    !3 matrices
    REAL, DIMENSION(:,:), ALLOCATABLE :: imata, imatb, omat,wrongmat
    INTEGER :: arows, bcols, freedim, dim

    !cpu time
    REAL :: start, finish

    dim = read_dim("dim.txt")
    !dimensions of the matrices that we have to multiply
    arows=dim
    bcols=dim
    freedim=dim

    ALLOCATE(imata(arows,freedim))
    ALLOCATE(imatb(freedim,bcols))
    ALLOCATE(omat(arows,bcols))
    ALLOCATE(wrongmat(7,8))

    !generate matrix a
    imata=GEN_RND_MATR(arows,freedim)
    !generate matrix b
    imatb=GEN_RND_MATR(freedim,bcols)

    call set_debug(.TRUE.)

    CALL CPU_TIME(start)
    omat = MATMUL(imata, imatb)
    CALL CPU_TIME(finish)
    write(*,*) "matmul  ", finish -start
    call write_float_out("./out/matmul.txt", dim, finish-start)



    CALL CPU_TIME(start)
    CALL MATMUL1(imata, imatb, omat)
    CALL CPU_TIME(finish)
    write(*,*) "matmul1 ", finish -start
    call write_float_out("./out/matmul1.txt", dim, finish-start)



    CALL CPU_TIME(start)
    CALL MATMUL2(imata, imatb, omat)
    CALL CPU_TIME(finish)
    write(*,*) "matmul2 ", finish -start
    call write_float_out("./out/matmul2.txt", dim, finish-start)


    CALL CPU_TIME(start)
    CALL MATMUL3(imata, imatb, omat)
    CALL CPU_TIME(finish)
    write(*,*) "matmul3 ", finish -start
    call write_float_out("./out/matmul3.txt", dim, finish-start)


    CALL CPU_TIME(start)
    CALL MATMUL4(imata, imatb, omat)
    CALL CPU_TIME(finish)
    write(*,*) "matmul4 ", finish -start
    call write_float_out("./out/matmul4.txt", dim, finish-start)








end program
