module mmatrixqic
    type matrixqic
        integer :: dim
        double complex, dimension(:,:), allocatable :: el
        double complex :: trace
        double complex :: det
    end type matrixqic


    interface operator(.adj.)
        module procedure adjoint
    end interface

    interface operator(.tr.)
        module procedure trace
    end interface


    contains
        !Compute the determinant of a matrixqic
        subroutine Cdet(mat)
            implicit none

            type(matrixqic), intent(INOUT) :: mat

            !declarations
            integer :: info, ii, sgn
            double complex :: determinant
            integer, dimension(mat%dim):: ipiv
            double complex, dimension(mat%dim,mat%dim):: tmpmat

            tmpmat= mat%el
            ipiv =0


            call zgetrf(mat%dim, mat%dim, tmpmat, mat%dim, ipiv, info)

            !DO ii=1,mat%dim
            !     write(*, "(*('('sf6.2xspf6.2x'i)':x))")  tmpmat(ii,:)
            !  PRINT*, ""
            !END DO
            !print*, ipiv


            if (info .ne. 0) then
                write(*,*) "Error broken determinant computation"
                write(*,*) "Determinant Not Computed"


            else
                determinant =1
                sgn =1
                do ii=1, mat%dim
                   determinant = determinant*tmpmat(ii,ii)
                   if(ipiv(ii) .ne. ii) then
                    sgn = -sgn
                   end if
                end do
                determinant = sgn*determinant
                mat%det = determinant
            end if
            return
        end subroutine Cdet

        subroutine Gdet(matr, mydet)
            implicit none
            type(matrixqic), intent(IN) :: matr
            double complex, intent(OUT) :: mydet
            mydet = matr%det

            return
        end subroutine Gdet

        subroutine Gtrace(matr, mytrace)
            implicit none
            type(matrixqic), intent(IN) :: matr
            double complex, intent(OUT) :: mytrace
            mytrace = matr%trace
            return
        end subroutine Gtrace

        subroutine Init(N, result)
            implicit none

            integer, intent(IN) :: N
            type (matrixqic), intent(OUT) :: result
            !Print*, "Constructor of empty matrixqic.."

            allocate(result%el(N,N))
            result%dim = N
            result%el =0
            result%det =0
            result%trace =0
        end subroutine Init

        !Constructor of matrixqic
        subroutine Build(mat, result)
            implicit none
            !declarations

            double complex, dimension(:,:), intent(IN) :: mat

            type(matrixqic), intent(OUT) :: result

            integer :: ii, jj
            integer, dimension(size(shape(mat))) :: thissize
            integer :: N

            !exec-
            !Print*, "Run Construnctor of matrixqic"
            thissize = shape(mat)
            if(size(thissize) .ne. 2) then
                write(*,*) "Error, mat must be of dim = 2\n"
            else
                if  (thissize(1) .ne. thissize(2))  then
                    write(*,*) "Error, input be a square mat\n"
                else
                    N= thissize(1)
                    result%dim = N
                    result%trace = 0
                    allocate (result%el(N,N))
                    do ii=1,N
                        do jj=1,N
                            result%el(ii,jj) = mat(ii,jj)
                        end do
                        result%trace = result%trace + result%el(ii,ii)
                    end do
                    !write(*,*)"    - determinant"
                    call cdet(result)
                    !write(*,*) "                 done!",result%det

                end if
            end if
            Print*,""
            return
        end subroutine Build

        subroutine build_empty(mat, result)
            implicit none
            !declarations

            double complex, dimension(:,:), intent(IN) :: mat

            type(matrixqic), intent(OUT) :: result

            integer :: ii, jj
            integer, dimension(size(shape(mat))) :: thissize
            integer :: N

            !exec-
            !Print*, "Run Construnctor of matrixqic, no det, no trace"
            thissize = shape(mat)
            if(size(thissize) .ne. 2) then
                write(*,*) "Error, mat must be of dim = 2\n"
            else
                if (thissize(1) .ne. thissize(2))  then
                    write(*,*) "Error, input be a square mat\n"
                else
                    N=thissize(1)
                    result%dim = N
                    result%trace = 0
                    allocate (result%el(N,N))
                    result%el = mat

                end if
            end if
            !Print*,""
            return
        end subroutine build_empty

        !Compute and get inverse matrix
        subroutine CGinv(mat, invmat)
            implicit none
            type(matrixqic), intent(IN) ::mat
            type(matrixqic), intent(OUT) ::invmat
            double complex,dimension(mat%dim,mat%dim)::tmpmati
            double complex, dimension(mat%dim) :: work
            integer, dimension(mat%dim) :: ipiv

            integer:: info,ii,error



            !Print*, "Computing Inverse..."


            info =0
            ipiv=0
            work=0
            tmpmati=mat%el




            !print*,"    - Compute Inverse Matrix... "
            !print*,"            - PLU decomposition "
            call zgetrf(mat%dim,mat%dim, tmpmati,mat%dim, ipiv, info)
            !print*,"                    done! "
            !print*,"            - Inversion "
            call zgetri(mat%dim,tmpmati,mat%dim,ipiv,work,mat%dim,info )
            !print*,"                    done! "

            call Build(tmpmati,invmat)
            !Print*,""




            return
        end subroutine CGinv

        !output matrixqic in a file
        subroutine Ofile(mat,outfile)
            implicit none
            type(matrixqic), intent(IN) :: mat
            character*40, intent(INOUT):: outfile

            character*3 :: stat = "new"
            logical :: file_exist = .FALSE.
            integer :: ii

            outfile = trim(outfile)
            inquire(FILE=outfile, EXIST=file_exist)

            if (file_exist) then
                stat = "old"
            end if

            open(unit = 2, file=outfile, status = stat)

            write(2,*) "dimension = "
            write(2,*) mat%dim
            write(2,*) "trace = "
            write(2, "(*(''sf6.2xspf6.2x'i ':x))") mat%trace
            write(2,*) "determinant = "
            write(2, "(*(''sf6.2xspf6.2x'i ':x))") mat%det
            write(2,*) "full matrix = "


            do ii=1,mat%dim
                 write(2, "(*(''sf6.2xspf6.2x'i ':x))")  mat%el(ii,:)
            end do
            close(2)
            return
        end subroutine Ofile

        subroutine cgadj(mat, adjmat)
            implicit none
            type(matrixqic), intent(IN) :: mat
            type(matrixqic), intent(OUT) :: adjmat

            integer ii , jj
            double complex, dimension(mat%dim, mat%dim):: tmpmat
            double precision:: im, re

            !Print*, "Computing Adjoint..."

            tmpmat = transpose(mat%el)

            tmpmat = conjg(tmpmat)
            !do ii=1,mat%dim
            !    do jj=1,mat%dim
            !        im =  -aimag(tmpmat(ii,jj))
            !        re =  real(tmpmat(ii,jj))
            !        tmpmat(ii,jj) = cmplx(re, im, kind (0D0))
            !    end do
            !end do
            !call build(mat%dim, tmpmat,adjmat)
            call build_empty(tmpmat,adjmat)
            !Print*,""
            return
        end subroutine cgadj

        subroutine ctrace(mat)
            implicit none
            type(matrixqic) , intent(INOUT) :: mat
            integer:: ii
            mat%trace=0
            do ii =1, mat%dim
                mat%trace= mat%trace + mat%el(ii,ii)
            end do
        end subroutine ctrace

        function Adjoint(mat)
            implicit none
            type(matrixqic), intent(IN)::mat
            type(matrixqic) :: adjoint

            call cgadj(mat,adjoint)

        end function Adjoint

        function Trace(mat)
            implicit none
            type(matrixqic), intent(IN)::mat
            type(matrixqic) :: tmp
            double complex :: trace
            tmp = mat
            call ctrace(tmp)
            trace = tmp%trace

        end function Trace



end module mmatrixqic
