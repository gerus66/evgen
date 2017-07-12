!-- pp-interaction
!-- requied: pp.txt - config. file in current directory
!--          pp/ - folder with pp-data files (for extrapolation) in current directory
!--          log/ - folder for log information in current directory
!-- used units: 21 - pp_interact()
!--            22 - read_pp()
!-- used pp-data files (separate for S=0 and S=1) should have the same step for energy
!-- ATTENTION! don't forget to use 'dealloc_pp'

module pp
    use cons, only: kind, prec, pi
    implicit none

!== private ===============================================
    real(kind), dimension(:,:), allocatable, private :: Pdata !-- (:,1:3) for S=0, (:,4:6) for S=1 
    real(kind), private :: step = 0.D0  !-- step for energy in pp-data file
        
    contains 
    
    real(kind) function get_Wpp(is_sym, en)   !-- Wpp(xi) = a*xi + b
        logical, intent(in) :: is_sym
        real(kind), intent(in) :: en
        
        integer :: i, i_st
        
        i = int( en / step ) + 1
        if ( is_sym ) then
            i_st = 1
        else
            i_st = 4
        endif
        get_Wpp = Pdata(i, i_st+1) * en + Pdata(i, i_st+2)
    end function 
    
    subroutine pp_dealloc()  
        deallocate (Pdata)
    end subroutine
     
    subroutine pp_interact(e)
        real(kind), intent(in) :: e  !-- reading pp-data up to this energy
        
        integer :: i
        character(50) :: pp_0, pp_1  !-- names of pp-data file S=0 and S=1 at pp/
        
        open(unit=21,file='pp.txt')
        read(21,*) pp_0; read(21,*) pp_1
        close(21)
        
        call read_pp(pp_0, 1)
        call read_pp(pp_1, 4)
       
        call normalize(1)
        call normalize(4)
        
        call prepare(1)
        call prepare(4)
        
        open(unit=21,file='log/diagn_pp.log') 
        write(21,'(5(a13, 1x))') '#', 'W0pp(E) rebld', 'W1pp(E) rebld'
        do i=1, size(Pdata,1)
            write(21,'(5(e13.6, 1x))') step*i, get_Wpp(.TRUE., step*i), get_Wpp(.FALSE., step*i)
        enddo
        close(21)
 
        print    "(/,20('-'),/,'pp-interaction is ready to use',/,20('-'),/)"
        
        contains
            
        subroutine read_pp(pp_name, i_st)
            character(*), intent(in) :: pp_name
            integer, intent(in) :: i_st
            
            integer :: i, n
            real(kind) :: serv_k
            
            open(unit=22, file='pp/'//pp_name)
            read(22,*) n    !-- length of file
            read(22,*) serv_k
            if ( step.LT.prec ) then
                step = serv_k !-- global for pp.mod
            else
                if ( dabs(step-serv_k).GT.prec ) stop 'we need equal energy-step for 0 and 1'
            endif
            i = int(e/step)+1   !-- requied length of file for 'e'(e comes from parent)
            if ( n.LT.i ) stop "length of pp-data is not enough"  
            if ( .NOT.allocated(Pdata) ) allocate(Pdata(i,6))
            do i=1, size(Pdata,1)
                read(22,*) Pdata(i,i_st)
            enddo
            close(22)
        end subroutine
        
        subroutine normalize(i_st)     !-- normilize to  pi / 8
            integer, intent(in) :: i_st
            
            integer :: i
            real(kind) :: serv_k
            
            serv_k = 0.D0
            do i=1, size(Pdata,1)-1
                serv_k = serv_k + Pdata(i,i_st) * sqrt(1.D0 - step*i/e)
            enddo
            Pdata(:,i_st) = Pdata(:,i_st) * pi / ( serv_k * step * 8.D0 )
        end subroutine
        
        subroutine prepare(i_st)
            integer, intent(in) :: i_st
            
            real(kind) :: serv_k
            integer :: i
            
            serv_k = sqrt( e**3 / step ) 
            Pdata(1,i_st) = Pdata(1,i_st) * serv_k
            Pdata(1,i_st+1) = Pdata(1,i_st) / step
            Pdata(1,i_st+2) = 0.D0
            do i=2, size(Pdata,1)            !-- / V
                Pdata(i,i_st) = Pdata(i,i_st) * serv_k / sqrt( real(i) )
                Pdata(i,i_st+1) = (Pdata(i,i_st) - Pdata(i-1,i_st)) / step
                Pdata(i,i_st+2) = Pdata(i-1,i_st)*i - Pdata(i,i_st)*(i-1)
            enddo 
        end subroutine
        
    end subroutine
    
end module
