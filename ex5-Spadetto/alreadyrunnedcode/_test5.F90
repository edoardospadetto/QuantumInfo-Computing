# 1 "test5.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "test5.F90"


# 1 "modules/matmulqic.F90" 1


# 1 "modules/debugqic.F90" 1



    module debugqic
    contains
        subroutine breakifn(d_string, cndtion, debug)
            implicit none
            character(len = *), intent(IN) :: d_string

            logical, intent(IN), optional :: debug
            logical, intent(IN) :: cndtion
            logical :: de
            ! if debug is enabled continue otherwise noop
            ! set debug to its right value
            de=.false.
            if( present(debug)) then
                de=debug
            end if

            if (de) then
                !if error occurs stop the program and print
                if (.not. cndtion) then
                    write (*,*) d_string
                    write(*,*) "Execution terminated"
                    STOP

                end if
            end if
        end subroutine breakifn

    end module debugqic


# 3 "modules/matmulqic.F90" 2

MODULE matmulqic
    use debugqic
    !variable to enable debug
    LOGICAL :: debug_var = .FALSE.

  CONTAINS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !enable debug for the function declared in this module
    SUBROUTINE set_debug(mybool)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LOGICAL, INTENT(IN) :: mybool
        debug_var = mybool
    END SUBROUTINE set_debug






    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication
    SUBROUTINE matmul1(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IMPLICIT NONE
      !inout vars
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result
      !temp vars
      INTEGER*2 :: ii,jj,kk
      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)
      !allocate memory to the result of multiplication
      ALLOCATE (result(sizea(1),sizeb(2)))

      !check equality between rows and columns using the debug subroutine
      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)

      !compute the multiplication using the loop
        DO ii=1,sizea(1)
          DO jj=1,sizeb(2)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb(kk,jj)
            END DO
          END DO
        END DO
    END SUBROUTINE matmul1



    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication with different order
    SUBROUTINE matmul2(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !everything is equal to matmul 1 except the order of the nested for loops
      IMPLICIT NONE
      INTEGER*2 :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))

      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)

        DO jj=1,sizeb(2)
          DO ii=1,sizea(1)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb(kk,jj)
            END DO
          END DO
        END DO


    END SUBROUTINE matmul2

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication with different order
    SUBROUTINE matmul3(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ! I introduced another variable, aa1 that's the transpose of aa. I use this in the loop to perform 2 loops over columns and one over rows.

      IMPLICIT NONE
      INTEGER*2 :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result
      REAL, DIMENSION(:,:), ALLOCATABLE :: aa1

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))
      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)

        ALLOCATE(aa1(sizea(2), sizea(1)))
        aa1=TRANSPOSE(aa)
        DO jj=1,sizeb(2)
          DO ii=1,sizea(1)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa1(kk,ii)*bb(kk,jj)
            END DO
          END DO
        END DO


    END SUBROUTINE matmul3

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !perform matrix multiplication with different order
    SUBROUTINE matmul4(aa,bb,result)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !as matmul3 but 2 loops on rows and one on columns transposing bb

      IMPLICIT NONE
      INTEGER*2 :: ii,jj,kk
      REAL, DIMENSION(:,:), INTENT(IN) :: aa,bb
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: result
      REAL, DIMENSION(:,:), ALLOCATABLE :: bb1

      INTEGER, DIMENSION(SIZE(SHAPE(aa))) :: sizea
      INTEGER, DIMENSION(SIZE(SHAPE(bb))) :: sizeb

      sizea=SHAPE(aa)
      sizeb=SHAPE(bb)

      ALLOCATE (result(sizea(1),sizeb(2)))

      call breakifn("Invalid Matrix Shapes", (sizea(2) .EQ. sizeb(1)), debug_var)


        ALLOCATE(bb1(sizeb(2), sizeb(1)))
        bb1=TRANSPOSE(bb)
        DO ii=1,sizea(1)
          DO jj=1,sizeb(2)
            result(ii,jj)=0
            DO kk=1,sizea(2)
              result(ii,jj) = result(ii,jj)+ aa(ii,kk)*bb1(jj,kk)
            END DO
          END DO
        END DO

    END SUBROUTINE matmul4

END MODULE matmulqic
# 3 "test5.F90" 2

# 1 "modules/matrixqic.F90" 1


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
# 4 "test5.F90" 2

module tmphelp
    use debugqic
    implicit none
    contains
    !###########################################################
    !computes mean of 1d array
    function mean(a) result(themean)
    !###########################################################
        implicit none

        real*8, dimension(:) :: a
        real*8 :: themean
        themean = dble(sum(a))/dble(size(a))
    end function

    !###########################################################
    !computes MOVINGAVERAGE of 1d array
    function movavg(a, period) result(movavgarray)
    !###########################################################
        implicit none
        real*8, dimension(:) :: a
        real*8, dimension(size(a)):: movavgarray
        real*8:: tmpavg
        integer :: ii, period
        call breakifn("2*Period > size", size(a) .ge. 2*period, .true.)
        do ii = period+1,size(a)-period
            tmpavg = dble(size(a(ii-period: ii+period)))/sum(a(ii-period: ii+period))
            movavgarray(ii) = a(ii)*tmpavg
        end do
        do ii = 1, period
            tmpavg = dble(size(a(1: 2*period)))/sum(a(1: 2*period))
            movavgarray(ii)=a(ii)*tmpavg
        end do
        do ii = size(a)-period+1, size(a)
            tmpavg = dble(size(a(size(a)-2*period: size(a))))/sum(a(size(a)-2*period: size(a)))
            movavgarray(ii)=a(ii)*tmpavg
        end do

    end function





    !###########################################################
    !computes eigenvalues of diagonal matrix ;)
    !and sort them in ascending order
    function deigs(a) result(eigs)
    !###########################################################
        implicit none
        real*8, dimension(:,:) :: a
        real*8 :: eigs(size(a, dim=1)), tmp
        integer :: ii , flag, dim
        !compute eigenvalues
        dim = size(a, dim=1)
        call breakifn("Square Matrix is needed " , size(a, dim=1) .eq. size(a, dim=2))
        do ii = 1, dim
                eigs(ii) = a(ii,ii)
        end do
        !sort
        flag =1
        do while (flag .eq. 1)
            flag = 0
            do ii = 1, dim-1
                if (eigs(ii) .ge. eigs(ii+1)) then
                    tmp = eigs(ii)
                    eigs(ii) = eigs(ii+1)
                    eigs(ii+1)= tmp
                    flag=1
                    !print*, ii
                end if
            end do
        end do
    end function

    !###########################################################
    !computes 2d array of numerical density distribution
    !first the array of frequency is computed from minv to maxv given a
    !certain number of bins nob, then the probability density is found by
    ! P_i approx n_i / N and p(x) = d(P(x))/dx apprx deltaP /delta x
    ! since P(x) is the cumulative distribution deltaP = P_i approx = n_i / N
    ! then we can extimate p(x) = n_i/(N*deltax)
    function density(array, nobs, minv, maxv) result(data)
    !array of observations, number of bins, where to start histogram, where to end histogram
    !###########################################################

        implicit none

        real*8 , dimension(:) :: array
        integer ii, nobs, where,cnter
        real*8 maxv, minv, data(nobs,2),delta

        delta = (maxv-minv)/dble(nobs)
        cnter =0
        do ii = 1,size(array)
            where = floor((array(ii)/delta) -minv +1)
            if ((where .ge. 1) .and.(where .le. nobs)) then
                cnter = cnter +1
                data(where,1) = data(where,1)+1
            end if
        end do
        do ii = 1,nobs
            data(ii,1) = data(ii,1)/(delta*cnter)
            data(ii,2) = minv+delta*(ii-1+0.5)
        end do



    end function




end module

program ex5
    use debugqic
    use matmulqic
    use matrixqic
    use tmphelp

    IMPLICIT NONE

    INTEGER          N, nummat, nob
    INTEGER*2        ii, jj, dolop , kk
    PARAMETER        (nob = 100)
    PARAMETER        (nummat = 20)
    PARAMETER        ( N = 1000 )
    INTEGER          LDA, LDVL, LDVR
    PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
    INTEGER          LWMAX, nws(1)
    PARAMETER        ( LWMAX = 100000 )
    INTEGER          INFO, LWORK, periods(4)
    DOUBLE PRECISION RWORK( 2*N-1 ), Aa(N,N)
    COMPLEX*16       A( LDA, N ), VL( LDVL, N ), VR( LDVR, N )
    double precision   W( N ), avgdelta, deltaeig(N-1)
    COMPLEX*16       WORK( LWMAX )
    double precision, dimension(nummat,N-1):: alldelta
    double precision, dimension(nob,2):: dns
    double precision :: movavgar((N-1)*nummat,4), rs(N-2), rsm, lims(4)
    character(len=12):: outfiles(2)
    character(len=1):: numb(4)
    numb =  (/'1','2','3','4'/)
    outfiles(1) = 'densityd.txt'
    outfiles(2) = 'densityh.txt'
    periods= (/5,25,50,75/)
    lims=(/3.0,3.0,3.0,3.0/)

do dolop = 1,2
rsm = 0

print*, dolop
! for each of the l random matrix
    do kk = 1, nummat
        !print*, kk
        if (dolop .eq. 2) then
            !#############compute eigenvalues################
             A = rghcm(N)
             LWORK = -1
             call zheev ( 'N','U',N,A,LDA,W,WORK,LWORK,RWORK,INFO)
             LWORK = MIN(LWMAX, INT(WORK(1)))
             call zheev ( 'N','U',N,A,LDA,W,WORK,LWORK,RWORK,INFO)
        else

            Aa=rgddm(N)
            W=deigs(Aa)

        end if
        !#######compute the normalized spacings#############
        avgdelta =0
        deltaeig =0
        do ii = 1, N-1
            deltaeig(ii)=W(ii+1)-W(ii) !compute spacing
        end do

        !moving average normalized spacings
        ! compute frequencies
        do ii = 1,4
            movavgar((kk-1)*(N-1)+1:kk*(N-1),ii) = movavg(deltaeig,periods(ii))
        end do

        !average r coefficient
        do jj = 1,N-2
            rs(jj) = min(deltaeig(jj),deltaeig(jj+1)) / max(deltaeig(jj),deltaeig(jj+1))
        end do
        rsm = rsm+ mean(rs)/nummat



        avgdelta = mean(deltaeig) !normalize spacings
        deltaeig = deltaeig/avgdelta


        do ii = 1, N-1
            !save avg spacing in a matrix
            alldelta(kk,ii) = deltaeig(ii) !compute spacing
        end do


    end do
    !new shape for the array containing all the s_i.
    nws(1)=(N-1)*nummat

    ! if DIAGONAL change the limits of the distributio
    if (dolop .eq. 1) then
        dns = density(reshape(alldelta, nws),nob,dble(0.0),dble(5.0))
    else
    !if hermitian other limits
        dns = density(reshape(alldelta, nws),nob,dble(0.0),dble(0.3))
    end if


    !stamp everything to file one for diagonal one for hermitian
    print*, "out to file"

    do kk = 1, nob
            call wffo("./results/" // outfiles(dolop),real(dns(kk,2)),real(dns(kk,1)))
    end do

    do ii = 1,4
        dns =  density(movavgar(:,ii), nob, dble(0.0), lims(ii))
        do kk = 1, nob
            call wffo("./results/prd" //numb(ii) // numb(dolop) // ".txt",real(dns(kk,2)),real(dns(kk,1)))
        end do
    end do

    !Print average r coefficient on screen
    print*, rsm

end do
end program
