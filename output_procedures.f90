module output
    use cons, only: kind
    use types
    implicit none

contains

!-- read from input stream one of numbers from list -------
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
end function
!----------------------------------------------------------

!-- count length of file (lines) --------------------------
integer function count_lines(stream)
    integer, intent(in) :: stream
    
    character(50) :: serv_c
    
    count_lines = 0
    do while (.NOT.eof(stream))
        count_lines = count_lines + 1
        read(stream,*) serv_c
    enddo
end function
!----------------------------------------------------------



subroutine matrix_to_file (A, dir, name, e_format)   !--- writes down matrix A to the file "log/[dir(8)][name(8)]" in current directory
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

end function histogram



subroutine recoup_print(C_ls, j)
    type(recoup), dimension(-1:1,0:1), intent(in) :: C_ls
    integer, intent(in) :: j
    
    integer :: i, ii
    
    write(5,"(/,'recouple coefficients',/,A5,2(A10),/)") 'L \ S', '0', '1'
    print   "(/,'recouple coefficients',/,A5,2(A10),/)", 'L \ S', '0', '1'
    do i = -1, 1
        write(5,"(I5,\)") j+i
        print   "(I5,\)", j+i
        do ii = 0, 1
            if ( C_ls(i,ii).not_null ) then
                write(5,"(F10.2,\)") C_ls(i,ii).c
                print   "(F10.2,\)", C_ls(i,ii).c
            else
                write(5,"(L10,\)") C_ls(i,ii).not_null
                print   "(L10,\)", C_ls(i,ii).not_null
            endif
        enddo
        write(5,"('')")
        print   "('')"
    enddo
end subroutine

subroutine clgd_print(C_lml, j)
    integer, intent(in) :: j
    type(ar_clgd), dimension(-1:1, -j-1:j+1), intent(in) :: C_lml
      
    integer :: i, ii, iii
        
    write(5,"(/,'Clebsh-Gordan coefficients',/,A5,\)") 'L \Ml'
    print   "(/,'Clebsh-Gordan coefficients',/,A5,\)", 'L \Ml'
    do i = -j-1, j+1
        write(5,"(I10,\)") i
        print   "(I10,\)", i
    enddo
    write(5,"('')")
    print   "('')"
    do i = -1, 1
        write(5,"(I5,\)") j+i
        print   "(I5,\)", j+i
        do ii = -j-1, j+1
            if ( C_lml(i,ii).n .NE. 0 ) then
                write(5,"(A7,I2,A1,\)") '[', C_lml(i,ii).n, ']'
                print   "(A7,I2,A1,\)", '[', C_lml(i,ii).n, ']'
            else
                write(5,"(L10,\)") C_lml(i,ii).n
                print   "(L10,\)", C_lml(i,ii).n
            endif
        enddo
        write(5,"('')")
        print   "('')"
    enddo
    
    write(5,"(A,/,4(A5),A8)") 'not-null only:', 'L', 'Ml', 'ml1', 'ml2', 'coef'
    print   "(A,/,4(A5),A8)", 'not-null only:', 'L', 'Ml', 'ml1', 'ml2', 'coef'
    do i = -1, 1
        do ii = -j-1, j+1
            if ( C_lml(i,ii).n .NE. 0 ) then
                do iii = 1, size(C_lml(i,ii).ar)
                    write(5, "(4(I5),f10.5)") C_lml(i,ii).ar(iii).l, C_lml(i,ii).ar(iii).ml, C_lml(i,ii).ar(iii).m1, C_lml(i,ii).ar(iii).m2, C_lml(i,ii).ar(iii).c  
                    print    "(4(I5),f10.5)", C_lml(i,ii).ar(iii).l, C_lml(i,ii).ar(iii).ml, C_lml(i,ii).ar(iii).m1, C_lml(i,ii).ar(iii).m2, C_lml(i,ii).ar(iii).c  
               enddo
           endif
       enddo
    enddo
    
end subroutine

end module


