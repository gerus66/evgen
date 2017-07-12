!-- histogram and output all about {Evs} 
!-- requied: log/ - folder for log information in current directory
!-- used units: 11 - data_output()
!--             12 - diagn_output()
!-- ATTENTION! don't forget to use 'set_name_out' before to output any files

module output_ev
    use cons, only: kind
    use vars, only: Evs
    use types, only: correlation
    implicit none

!== private ===============================================
    character(12), private :: name_out   !-- basical name for all output files (template is '..-._.-._.-.')
    integer, private :: nbin   !-- number of bins for all histograms
    logical, private :: is_ww
    real(kind), private :: norm
    character(1), private :: d
                   
    contains
    
    subroutine set_nbin(n)
        integer, intent(in) :: n
        
        nbin = n
    end subroutine
    
    subroutine set_name_out(e1, e2, e3)
        real(kind), intent(in) :: e1, e2, e3
        
        write( name_out, "(I1,'-',I2,'_',I1,'-',I1,'_',I1,'-',I1)" ) &
             floor(abs(e1)), int(mod(abs(e1)*100,100.D0)), floor(abs(e2)), int(mod(abs(e2)*10,10.D0)), floor(abs(e3)), int(mod(abs(e3)*10,10.D0))
    end subroutine
    
    subroutine set_is_ww(not_gen)
        logical, intent(in) :: not_gen
        
        is_ww = not_gen
        if ( not_gen ) then
            norm = size(Evs) / sum(Evs.ww)
            d = 'm'
        else
            norm = 1.D0
            d = 'g'
        endif    
    end subroutine
    
    subroutine data_output()
        integer :: i
        
        open(unit=11,file='log/'//name_out//'.data')
        do i = 1,size(Evs)
            write(11,"(10(1x,e12.5))"), Evs(i).eng, Evs(i).p1.x, Evs(i).p1.y, Evs(i).p1.z, &
                                                    Evs(i).p2.x, Evs(i).p2.y, Evs(i).p2.z, &
                                                    Evs(i).core.x, Evs(i).core.y, Evs(i).core.z 
        enddo
        close(11)
    end subroutine
    
    subroutine diagn_output()
        integer :: i
        
        call to_file(hist(Evs.y1.eps, 0.D0, 1.D0), '_eps_y1')
        call to_file(hist(Evs.y1.cos_teta, -1.D0, 2.D0), '_cos_y1')
        call to_file(hist(Evs.y2.eps, 0.D0, 1.D0), '_eps_y2')
        call to_file(hist(Evs.y2.cos_teta, -1.D0, 2.D0), '_cos_y2')
        call to_file(hist(Evs.t.eps, 0.D0, 1.D0), '_eps_t')
        call to_file(hist(Evs.t.cos_teta, -1.D0, 2.D0), '_cos_t')   
        call to_file_3d(hist_3d(Evs.y1), '_3d_y1')
        call to_file_3d(hist_3d(Evs.y2), '_3d_y2')
        call to_file_3d(hist_3d(Evs.t), '_3d_t')
        print "(/,'diagnostic has been done',/)" 
    end subroutine
   
    function hist(X, xmin, xrange)   
        real(kind), dimension(4,nbin) :: hist 
        real(kind), dimension(:), intent(in) :: X
        real(kind), intent(in) :: xmin, xrange
    
            real(kind) :: bin
            integer :: i, j

            bin = xrange / real(nbin)

            hist(1,:) = (/ ( xmin+(i-0.5D0)*bin, i = 1, nbin) /)   !-- set every point in the center of x-bin
            hist(2:4,:) = 0.D0
            do i = 1,size(X)
                j = ceiling( (X(i) - xmin) / bin )
                if (j .EQ. 0) j = 1
                if ( is_ww ) then
                    hist(2,j) = hist(2,j) + Evs(i).ww     !-- sum up data in every x-bin
                    hist(3,j) = hist(3,j) + Evs(i).ww_sym
                    hist(4,j) = hist(4,j) + Evs(i).ww_asym
                else
                    hist(2,j) = hist(2,j) + 1.D0    !-- add 1 event to x-bin of histogram
                endif
            enddo
            if ( is_ww ) hist(2:4,:) = hist(2:4,:) * norm   !-- normalize ww(total) to nbins-hystogram   
        end function
   
        function hist_3d(X)
            integer, dimension(nbin,nbin) :: hist_3d
            type(correlation), dimension(:), intent(in) :: X
    
            real(kind), dimension(nbin,nbin) :: R
            integer :: i, j, k

            R = 0.D0
            do i = 1,size(X)
                j = ceiling( (X(i).eps)*nbin )     !-- epsilon
                k = ceiling( (X(i).cos_teta+1.D0)/2.D0*nbin )   !-- cos theta
                if (j .EQ. 0) j = 1
                if (k .EQ. 0) k = 1
                if ( is_ww ) then
                    R(k,j) = R(k,j) + Evs(i).ww     !-- sum up data in xy-bin of plot (only WW-total)
                else
                    R(k,j) = R(k,j) + 1.D0    !-- add 1 event to xy-bin of histogram
                endif
            enddo 
            if ( is_ww ) R = R * norm   !-- normalize plot(WW-total) to histogram
            hist_3d = nint(R)
        end function 

        subroutine to_file(A, path)  
            real(kind), intent(in), dimension(4,nbin) :: A 
            character(*), intent(in) :: path

            integer :: i

            open(unit=12,file='log/'//name_out//d//path//'.csv')
            do i = 1, nbin
                write(12,"(E13.6,' ',\)") A(1,i)
                write(12,"(3(I,' '),\)") nint(A(2:4,i))
                write(12,"(/,\)")
            enddo
            close(12)
        end subroutine
        
        subroutine to_file_3d(A, path)  
            integer, intent(in), dimension(nbin,nbin) :: A
            character(*), intent(in) :: path

            integer :: i

            open(unit=12,file='log/'//name_out//d//path//'.csv')
            do i = 1, nbin
                write(12,"(I5,' ',\)") A(:,i)
                write(12,"(/,\)")
            enddo
            close(12)
        end subroutine
       
end module
