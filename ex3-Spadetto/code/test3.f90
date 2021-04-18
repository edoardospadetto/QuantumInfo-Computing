
include "matmulqic.f90"

program ex3

    use debugqic
    use matmulqic

    IMPLICIT NONE
    !3 matrices
    REAL, DIMENSION(:,:), ALLOCATABLE :: imata, imatb, omat,wrongmat
    INTEGER :: arows, bcols, freedim

    !cpu time
    REAL :: start, finish


    !dimensions of the matrices that we have to multiply
    arows=302
    bcols=222
    freedim=412

    ALLOCATE(imata(arows,freedim))
    ALLOCATE(imatb(freedim,bcols))
    ALLOCATE(omat(arows,bcols))
    ALLOCATE(wrongmat(7,8))

    !generate matrix a
    imata=GEN_RND_MATR(arows,freedim)
    !generate matrix b
    imatb=GEN_RND_MATR(freedim,bcols)

    !CALL PRINTMATRIX(imata)
    !CALL PRINTMATRIX(imatb)
    PRINT*," "
    WRITE (*,*) "------------------------------------"
    call set_debug(.TRUE.)

    CALL CPU_TIME(start)
    omat = MATMUL(imata, imatb)
    CALL CPU_TIME(finish)
    PRINT*," "
    WRITE (*,*) "Fortran MATMUL Elapsed Time:  ",  finish-start
    !CALL PRINTMATRIX(omat)
    WRITE (*,*) "------------------------------------"
    PRINT*," "


    CALL CPU_TIME(start)
    CALL MATMUL1(imata, imatb, omat)
    CALL CPU_TIME(finish)
    PRINT*," "
    WRITE (*,*) "Custom MATMUL1 Elapsed Time:  ",  finish-start
    !CALL PRINTMATRIX(omat)
    WRITE (*,*) "------------------------------------"
    PRINT*," "


    CALL CPU_TIME(start)
    CALL MATMUL2(imata, imatb, omat)
    CALL CPU_TIME(finish)
    PRINT*," "
    WRITE (*,*) "Custom MATMUL2 Elapsed Time:  ",  finish-start
    !CALL PRINTMATRIX(omat)
    WRITE (*,*) "------------------------------------"
    PRINT*," "

    CALL CPU_TIME(start)
    CALL MATMUL3(imata, imatb, omat)
    CALL CPU_TIME(finish)
    PRINT*," "
    WRITE (*,*) "Custom MATMUL3 Elapsed Time:  ",  finish-start
    !CALL PRINTMATRIX(omat)
    WRITE (*,*) "------------------------------------"
    PRINT*," "

    CALL CPU_TIME(start)
    CALL MATMUL4(imata, imatb, omat)
    CALL CPU_TIME(finish)
    PRINT*," "
    WRITE (*,*) "Custom MATMUL4 Elapsed Time:  ",  finish-start
    !CALL PRINTMATRIX(omat)
    WRITE (*,*) "------------------------------------"
    PRINT*," "


    print*, "If i insert a bad condition in the debug subroutine?"
    CALL MATMUL1(imata, wrongmat, omat)




end program
