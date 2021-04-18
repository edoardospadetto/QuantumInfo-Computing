
#include "modules/matmulqic.F90"
#include "modules/matrixqic.F90"

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
