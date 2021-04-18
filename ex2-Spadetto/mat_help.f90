!Module to perform matrix multiplication, 2 different functions with different for
!loops order, then other two functions one to generate random matrices, one to print matrices.
!COMPILE : gfortran -g ex1.f90 -o ex1_c -O1/2/3 to change optimization
!For high dimensions fortran standard function is always better

MODULE zero
  CONTAINS
    !generates random matrix
    FUNCTION gen_rnd_matr(sizerow,sizecol) RESULT(gen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: sizerow,  sizecol
      REAL, DIMENSION(sizerow,sizecol) :: gen
      call RANDOM_NUMBER(gen)
    END FUNCTION gen_rnd_matr

    FUNCTION gen_rnd_matrc(sizerow,sizecol) RESULT(gen)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: sizerow,  sizecol
      INTEGER ::  ii, jj
      DOUBLE precision :: a, b
      DOUBLE COMPLEX, DIMENSION(sizerow,sizecol) :: gen
      
      do ii= 1, sizerow
          do jj = 1, sizecol
              call RANDOM_NUMBER(a)
              call RANDOM_NUMBER(b)
              gen(ii,jj) = cmplx(a,b,kind(0D0))
          end do
      end do



  END FUNCTION gen_rnd_matrc

    !print a matrix on terminal
    SUBROUTINE printmatrix(mat)
      IMPLICIT NONE
      INTEGER :: ii,jj
      REAL, DIMENSION(:,:), INTENT(IN) :: mat
      INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
      sizem=SHAPE(mat)

      PRINT*, ""
      DO ii=1,sizem(1)
        DO jj=1,sizem(2)
          WRITE(*, fmt="(f0.2, tr2)", advance="no") mat(ii,jj)
        END DO
        PRINT*, ""
      END DO


    END SUBROUTINE printmatrix

    !print a complexmatrix on terminal
    SUBROUTINE printmatrixc(mat)
      IMPLICIT NONE
      INTEGER :: ii,jj
      DOUBLE COMPLEX, DIMENSION(:,:), INTENT(IN) :: mat
      INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
      sizem=SHAPE(mat)

      PRINT*, ""
      DO ii=1,sizem(1)
           write(*, "(*('('sf6.2xspf6.2x'i)':x))")  mat(ii,:)
        PRINT*, ""
      END DO


  END SUBROUTINE printmatrixc

    !perform matrix multiplication
    SUBROUTINE matmul1(aa,bb,result)
      IMPLICIT NONE
      INTEGER :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))

      IF (sizea(2) .NE. sizeb(1)) THEN

        PRINT *, "Invalid matrix shapes"

      ELSE

        DO ii=1,sizea(1)
          DO jj=1,sizeb(2)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb(kk,jj)
            END DO
          END DO
        END DO


      END IF
    END SUBROUTINE matmul1

    !perform matrix multiplication with different order
    SUBROUTINE matmul2(aa,bb,result)
      IMPLICIT NONE
      INTEGER :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))

      IF (sizea(2) .NE. sizeb(1)) THEN

        PRINT *, "Invalid matrix shapes"

      ELSE

        DO jj=1,sizeb(2)
          DO ii=1,sizea(1)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb(kk,jj)
            END DO
          END DO
        END DO


      END IF
    END SUBROUTINE matmul2

END MODULE zero
