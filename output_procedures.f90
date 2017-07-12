module output
    use mydata, only: kind
    implicit none

contains

subroutine matrix_to_file (A, dir, name, e_format)                   !--- writes down matrix A to the file "log/[dir(8)][name(8)]" in current directory
    real(kind), intent(in), dimension(:,:) :: A                 !--- with format "([e_format(15)],\)" of each element
    character(*), intent(in) :: dir, name, e_format

    integer :: i, j

    open(unit=2,file=dir//name)
    do i = 1,size(A,1)
        do j = 1,size(A,2)
            write(2,"("//e_format//",\)") A(i,j)
        enddo
        write(2,"(/,\)")
    enddo
    close(2)
end subroutine matrix_to_file



function histogram(A, nbin, fl, x_min, range)        !--- splits array A to "nbin"-bins histogram
    use mydata, only: kind                           !--- takes "min element" and "range" values from arguments (if fl = true), 
    implicit none                                    !--- or calculate it

    real(kind), intent(in), dimension(:) :: A
    integer, intent(in) :: nbin
    logical, intent(in) :: fl
    real(kind), intent(in) :: x_min, range
    real(kind), dimension(nbin,2) :: histogram

    real(kind) :: xmin, bin
    integer :: i, j

    histogram = 0.D0

    if ( fl ) then
        xmin = x_min
        bin = range / real(nbin)
    else
        xmin = minval(A)
        bin = (maxval(A) - xmin) / real(nbin)
    endif

    histogram(:,1) = (/ ( xmin+(i-0.5D0)*bin, i = 1, nbin) /)

    do i = 1,size(A)
        j = ceiling( (A(i) - xmin) / bin )
        if (j .EQ. 0) j = 1
        histogram(j,2) = histogram(j,2) + 1.D0
    enddo

endfunction histogram


subroutine short_3d_diagnostic(A, B, nbin, dir, name)         !--- split array to [nbin] parts and write down to the file (3 columns: x,y,f(x,y))
    real(kind), intent(in), dimension(:) :: A, B
    integer, intent(in) :: nbin
    character(*), intent(in) :: dir
    character(*), intent(in) :: name

    integer, dimension(1:nbin,1:nbin) :: T
    integer i, j, k, ii

    T = 0
    do i = 1,size(A)
        j = ceiling((A(i)+1)/2.*nbin)
        k = ceiling(B(i)*nbin)
        T(j,k) = T(j,k) + 1
    enddo
    open(unit=2,file='log/'//trim(dir)//trim(name))
    do i = 1,nbin
        do ii = 1,nbin
            write(2,"(f4.2,f4.2,i)") 2./nbin*i - 1., 1./nbin*ii, T(i,ii)
        end do      
    end do
    close(2)
end subroutine short_3d_diagnostic

function ww_diagnostic(E, F, S, A, nbin, norm, xmin, range)
    real(kind), dimension(nbin,4) :: ww_diagnostic
    real(kind), dimension(:), intent(in):: E, F, S, A    !-- Eps, WW: Full, Symmetric, Asymmetric
    integer, intent(in) :: nbin
    real(kind), intent(in) :: norm, xmin, range
    
    real(kind) :: bin
    integer :: i, j

    bin = range / real(nbin)
    ww_diagnostic(:,1) = (/ ( xmin+(i-0.5D0)*bin, i = 1, nbin) /)
    ww_diagnostic(:,2:4) = 0.D0
    do i = 1,size(E)
        j = ceiling( (E(i) - xmin) / bin )
        if (j .EQ. 0) j = 1
        ww_diagnostic(j,2) = ww_diagnostic(j,2) + F(i)
        ww_diagnostic(j,3) = ww_diagnostic(j,3) + S(i)
        ww_diagnostic(j,4) = ww_diagnostic(j,4) + A(i)
    enddo 
    ww_diagnostic(:,2:4) = ww_diagnostic(:,2:4) * norm
endfunction

function ww_3d_diagnostic(E, C, F, S, A, nbin, norm)
    real(kind), dimension(nbin*nbin,5) :: ww_3d_diagnostic
    real(kind), dimension(:), intent(in):: E, C, F, S, A    !-- {Eps, Cos Theta}, WW: Full, Symmetric, Asymmetric
    integer, intent(in) :: nbin
    real(kind), intent(in) :: norm
    
    integer :: i, ii, j, k

    do i = 1, nbin
        do ii = 1, nbin
            ww_3d_diagnostic((i-1)*nbin+ii,1) = (i - 0.5D0) / nbin     !-- epsilon {0, 1}
            ww_3d_diagnostic((i-1)*nbin+ii,2) = (ii - 0.5D0) * 2.D0 / nbin - 1.D0   !-- cos theta {-1, 1}
        enddo
    enddo
    ww_3d_diagnostic(:,3:5) = 0.D0
    do i = 1,size(E)
        j = ceiling( E(i)*nbin )    !-- epsilon
        k = ceiling( (C(i)+1.D0)/2.D0*nbin )   !-- cos theta
        if (j .EQ. 0) j = 1
        if (k .EQ. 0) k = 1
        ww_3d_diagnostic((j-1)*nbin+k,3) = ww_3d_diagnostic((j-1)*nbin+k,3) + F(i)
        ww_3d_diagnostic((j-1)*nbin+k,4) = ww_3d_diagnostic((j-1)*nbin+k,4) + S(i)
        ww_3d_diagnostic((j-1)*nbin+k,5) = ww_3d_diagnostic((j-1)*nbin+k,5) + A(i)
    enddo 
    ww_3d_diagnostic(:,3:5) = ww_3d_diagnostic(:,3:5) * norm
endfunction

!-- read from input stream one of list numbers ------------
!----------------------------------------------------------
integer function read_num(stream, Num)
    integer, dimension(:), intent(in) :: Num
    integer, intent(in) :: stream

    logical :: fl
    integer :: i

    fl = .FALSE.
    read(stream,*) read_num
    do i = 1, size(Num)
        if ( read_num .EQ. Num(i) ) fl = .TRUE.
    enddo
    if (.NOT. fl) stop 'wrong symbol found'
endfunction
!----------------------------------------------------------

!-- it makes smth like '2-35_1-8_0-3' ---------------------
!-- from e = 2.35, er = 1.8, tet = 0.3 --------------------
    character(12) function make_name(e, er, tet)
        real(kind), intent(in) :: e, er, tet

        write( make_name, "(I1,A1,I2,A1,I1,A1,I1,A1,I1,A1,I1)" )  floor(e),'-',int(mod(e*100,100.D0)),'_',floor(er),'-',int(mod(er*10,10.D0)), '_',floor(tet),'-',int(mod(tet*10,10.D0))
    endfunction
!----------------------------------------------------------
endmodule


