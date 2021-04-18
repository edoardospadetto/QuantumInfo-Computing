#include "debugqic.F90"

module matrixqic
    use debugqic

    contains
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !generate random complex double matrix hermitian
    FUNCTION rghcm(tsize) RESULT(gen)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: tsize
      INTEGER ::  ii, jj
      DOUBLE precision :: a, b
      DOUBLE COMPLEX, DIMENSION(tsize,tsize) :: gen

      do ii= 1, tsize
          do jj = 1, ii
              call RANDOM_NUMBER(a)
              call RANDOM_NUMBER(b)
              gen(ii,jj) = cmplx(a,b,kind(0D0))
              if (ii .NE. jj) then
                  gen(jj,ii) = cmplx(a,-b,kind(0D0))
              end if
          end do
      end do
  END FUNCTION rghcm

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !generate random complex double matrix
  FUNCTION rgcm(sizea,sizeb) RESULT(gen)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sizea,sizeb
    INTEGER ::  ii, jj
    DOUBLE precision :: a, b
    DOUBLE COMPLEX, DIMENSION(sizea,sizeb) :: gen

    do ii= 1, sizea
        do jj = 1, sizeb
            call RANDOM_NUMBER(a)
            call RANDOM_NUMBER(b)
            gen(ii,jj) = cmplx(a,b,kind(0D0))
        end do
    end do
END FUNCTION rgcm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!generate random double diagonal matrix
FUNCTION rgddm(sizem) RESULT(gen)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sizem
  INTEGER ::  ii, jj
  DOUBLE precision :: a
  DOUBLE precision, DIMENSION(sizem,sizem) :: gen

  do ii= 1, sizem
      do jj = 1, sizem
          gen(ii,jj) =0.0
      end do
      call RANDOM_NUMBER(a)
      gen(ii,ii) = a
  end do
END FUNCTION rgddm


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !generates random real matrix
    FUNCTION rgrm(sizerow,sizecol) RESULT(gen)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: sizerow,  sizecol
      REAL, DIMENSION(sizerow,sizecol) :: gen
      call RANDOM_NUMBER(gen)
  END FUNCTION rgrm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !print a matrix on terminal
    SUBROUTINE prm(mat)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE
      !inout vars
      REAL, DIMENSION(:,:), INTENT(IN) :: mat
      !tmp vars
      INTEGER*2 :: ii,jj
      INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem

      sizem=SHAPE(mat)

      PRINT*, ""
      DO ii=1,sizem(1)
        DO jj=1,sizem(2)
          WRITE(*, fmt="(f0.2, tr2)", advance="no") mat(ii,jj)
        END DO
        PRINT*, ""
      END DO

  END SUBROUTINE prm

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !read dimension from file.
    FUNCTION read_dim(path) RESULT(dim)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    character(len = *) :: path
    integer :: dim
    logical :: file_exist = .FALSE.
    INQUIRE(FILE=path, EXIST=file_exist)

    call breakifn(path // "does not exists",file_exist, .TRUE.)

    open(unit = 3, file = path , status = "old" )

    read(3, fmt= "(i8)") dim

    close(3)

    END FUNCTION

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!print a complexmatrix on terminal
SUBROUTINE pcm(mat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   IMPLICIT NONE
   INTEGER :: ii,jj
   DOUBLE COMPLEX, DIMENSION(:,:), INTENT(IN) :: mat
   INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
   sizem=SHAPE(mat)

   PRINT*, ""
   DO ii=1,sizem(1)
        write(*, "(*(' 'sf6.2xspf6.2x'i ':x))")  mat(ii,:)
     PRINT*, ""
   END DO


END SUBROUTINE pcm


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as INT   REAL
    SUBROUTINE wfio(outfile, myint, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        real , intent(IN) :: myreal
        integer , intent(IN):: myint
        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "new"

        inquire(FILE=outfile, EXIST=file_exist)

        if (file_exist) then
            stat = "old"
        end if

        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(I5,E22.3)") myint,myreal

        close(2)


    END SUBROUTINE


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as real   REAL
    SUBROUTINE wffo(outfile, myr, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        real , intent(IN) :: myreal,myr

        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "old"

        inquire(FILE=outfile, EXIST=file_exist)

        if (.not. file_exist) then
            stat = "new"
        else if (file_exist) then
            stat = "old"
        end if



        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(E22.3,E22.3)") myr,myreal

        close(2)


    END SUBROUTINE

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !add line to file outfile as REAL
    SUBROUTINE wdo(outfile, myreal)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        double precision , intent(IN) :: myreal
        character(len = *), intent(IN):: outfile
        logical::  file_exist = .FALSE.
        character(len = 3) :: stat = "new"

        inquire(FILE=outfile, EXIST=file_exist)

        if (file_exist) then
            stat = "old"
        else
            stat="new"
        end if

        open(unit = 2, file=outfile, status = stat, access = 'append')

        write(2, "(E22.3)") myreal

        close(2)


    END SUBROUTINE










end module
