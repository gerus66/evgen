!-- fit two square arrays
!-- requied: fit.txt - config. file in current directory
!--          fit/ - folder with array for fitting in current directory
!--          log/ - folder for log information in current directory
!-- used units: 42 - fit_read()
!--             43 - log/fit.log
!--             44 - log/fitplot.csv

module fit
    use cons, only: kind, sq2, prec, pi
    use output
    implicit none

!== private ===============================================
!-- parameters --------------------------------------------
    real(kind), parameter, private :: alph_k = 1.358D0   !-- alpha = 0.05 (for Kolmogorov test)
    integer, parameter, private :: l_stat = 20  !-- lower limit for statistic (for Kolmogorov test)
!-- types -------------------------------------------------   
    type grid 
        real(kind) :: cur, st, fin, step   !-- current, start, finish, step
    end type    
    interface set_null
        module procedure set_null_grid
    end interface
    interface operator(*)   !-- multiply 'grid' to real in form 'k*grid'
        module procedure mult
    end interface
!-- variables ---------------------------------------------    
    type(grid), private :: a, phi_a, b, phi_b     !-- grid parameters  ( a + b + c = 1, phi_a|phi_b <= 2pi, phi_c = 0 )
    integer, dimension(:,:), allocatable, private :: T, Y   !-- arrays for compare
!-- interfaces --------------------------------------------
    interface xi_sq
        module procedure xi   !-- one-dimension array: real, expected
        module procedure xi_3d   !-- two-dimension array: real, expected
    end interface
!==========================================================     
    
    contains

!== facility for 'grid' type interfaces ===================      
    type(grid) function set_null_grid(g)
        type(grid), intent(in) :: g
        
        set_null_grid = grid(0.D0, 0.D0, 0.D0, 0.D0)
    end function
    
    type(grid) function mult(k, g)
        real(kind), intent(in) :: k
        type(grid), intent(in) :: g
               
        mult = grid(g.cur*k, g.st*k, g.fin*k, g.step*k)
    end function
!==========================================================    
 
!== general functions =====================================

!-- xi-squared for 1-dimension array ----------------------
    real(kind) function xi(A1, A2)
        integer, dimension(:), intent(in) :: A1, A2  !-- A1(real) and A2(expected) should have equal size
        
        integer :: i
        real(kind), dimension(size(A1)) :: X1, X2
        
        xi = 0.D0
        X1 = dfloat(A1) / dfloat(sum(A1)) 
        X2 = dfloat(A2) / dfloat(sum(A2))
        do i = 1, size(A1)
            xi = xi + (X1(i) - X2(i))**2 / X2(i)
        enddo
        xi = xi * dfloat(sum(A1)) / dfloat(size(A1))
    end function

!-- xi-squared for 2-dimension array ---------------------- 
    real(kind) function xi_3d(A1, A2)
        integer, dimension(:,:), intent(in) :: A1, A2   !-- A1(real) and A2(expected) should have equal size
        
        integer :: i , ii
        real(kind), dimension(size(A1,1),size(A1,2)) :: X1, X2

        xi_3d = 0.D0
        X1 = dfloat(A1) / dfloat(sum(A1)) 
        X2 = dfloat(A2) / dfloat(sum(A2))
        do i = 1, size(X1,1)
            do ii = 1, size(X1,2)
                if ( X2(i,ii).NE.0 ) xi_3d = xi_3d + (X1(i,ii) - X2(i,ii))**2 / X2(i,ii)
            enddo
        enddo 
        xi_3d = xi_3d * dfloat(sum(A1)) / dfloat(size(A1,1)) / dfloat(size(A1,2))
    end function     





!========================================================== 
    
!== switch-on     
!-- Don't forget to use 'fit_dealloc()' in the end of program
    subroutine fit_read()
        character(50) :: path_t, path_y
        integer :: i, n
        
        open(unit=42,file='fit.txt')
        read(42,*) path_t   !-- path to 1st array for fitting
        read(42,*) path_y   !-- path to 2nd array for fitting
        a = read_grid(42)   !-- grid parameters
        phi_a = ( pi / 180.D0 ) * read_grid(42)   !-- grad to rad
      !  phi_a = phi_a * ( pi / 180.D0 )  
        b = read_grid(42)   !-- grid parameters
        phi_b = ( pi / 180.D0 ) * read_grid(42)   !-- grad to rad
     !   phi_b = phi_b * ( pi / 180.D0 )  !-- grad to rad
        close(42)
        open(unit=42,file='fit/'//path_t)  !-- file with 1st array for fitting
        n = count_lines(42)
        allocate(T(n,n)) 
        rewind(42)
        do i = 1, n   !-- read array
            read(42,*) T(:,i)
        enddo
        close(42)
        
        open(unit=42,file='fit/'//path_y)  !-- file with 2nd array for fitting
        n = count_lines(42)
        allocate(Y(n,n)) 
        rewind(42)
        do i = 1, n   !-- read array
            read(42,*) Y(:,i)
        enddo
        close(42)
        open(unit=43,file='log/fit.log')
        open(unit=44,file='log/fitplot.csv')
                
        contains
        
        type(grid) function read_grid(stream)
            integer, intent(in) :: stream
       !     type(grid), intent(out) :: g
            
            read(stream,*) read_grid.st
            read(stream,*) read_grid.fin
            read(stream,*) read_grid.step
            read_grid.cur = read_grid.st            
        end function       
    end subroutine
    
    subroutine fit_dealloc()
        deallocate(T,Y)
        close(43)
        close(44)
    end subroutine

!----------------------------------------------------------
    
    subroutine get_ab(norma, normb, ca, cb, is_end)
        real(kind), intent(in) :: norma, normb
        complex(kind), intent(out) :: ca, cb
        logical, intent(out) :: is_end
        
        ca = cmplx(sqrt(a.cur/norma)*cos(phi_a.cur), sqrt(a.cur/norma)*sin(phi_a.cur))
        cb = cmplx( sqrt(abs(1.D0-a.cur)/normb), 0.D0)
        write(44,"(2(F5.3,' '),\)") a.cur, phi_a.cur
        write(43,"('a: ',I3,'% ',I3,'grad',/,'b: ',I3,'% ',I3,'grad')") &
               nint(a.cur*100), nint(phi_a.cur/pi*180.D0), nint((1.D0-a.cur)*100), 0
        print    "('a: ',I3,'% ',I3,'grad',/,'b: ',I3,'% ',I3,'grad')", &
               nint(a.cur*100), nint(phi_a.cur/pi*180.D0), nint((1.D0-a.cur)*100), 0
        is_end = .FALSE.
        phi_a.cur = phi_a.cur + phi_a.step
        if ( (phi_a.cur-phi_a.fin).GT.prec ) then
            phi_a.cur = phi_a.st
            a.cur = a.cur + a.step 
            if ( (a.cur-a.fin).GT.prec ) then
                is_end = .TRUE.
                a.cur = a.st
            endif
        endif              
    end subroutine
    
    subroutine get_abc(norma, normb, normc, ca, cb, cc, is_end)
        real(kind), intent(in) :: norma, normb, normc
        complex(kind), intent(out) :: ca, cb, cc
        logical, intent(out) :: is_end
        
        logical :: is_valid
        
        ca = cmplx(sqrt(a.cur/norma)*cos(phi_a.cur), sqrt(a.cur/norma)*sin(phi_a.cur))
        cb = cmplx(sqrt(b.cur/normb)*cos(phi_b.cur), sqrt(b.cur/normb)*sin(phi_b.cur))
        cc = cmplx( sqrt(abs(1.D0-a.cur-b.cur)/normc), 0.D0)
        write(44,"(4(F5.3,' '),\)") a.cur, phi_a.cur, b.cur, phi_b.cur
        write(43,"('a: ',I3,'% ',I3,'grad',/,'b: ',I3,'% ',I3,'grad',/,'c: ',I3,'% ',I3,'grad')") &
             nint(a.cur*100), nint(phi_a.cur/pi*180.D0), nint(b.cur*100), nint(phi_b.cur/pi*180.D0), nint((1.D0-a.cur-b.cur)*100), 0
        print    "('a: ',I3,'% ',I3,'grad',/,'b: ',I3,'% ',I3,'grad',/,'c: ',I3,'% ',I3,'grad')", &
             nint(a.cur*100), nint(phi_a.cur/pi*180.D0), nint(b.cur*100), nint(phi_b.cur/pi*180.D0), nint((1.D0-a.cur-b.cur)*100), 0
        is_end = .FALSE.
        is_valid = .FALSE.
        do while ( .NOT.is_valid )
            phi_b.cur = phi_b.cur + phi_b.step
            if ( (phi_b.cur-phi_b.fin).GT.prec ) then
                phi_b.cur = phi_b.st
                b.cur = b.cur + b.step 
                if ( (b.cur-b.fin).GT.prec ) then
                    b.cur = b.st
                    phi_a.cur = phi_a.cur + phi_a.step
                    if ( (phi_a.cur-phi_a.fin).GT.prec ) then
                        phi_a.cur = phi_a.st
                        a.cur = a.cur + a.step 
                        if ( (a.cur-a.fin).GT.prec ) then
                            is_end = .TRUE.
                            a.cur = a.st
                        endif
                    endif  
                endif
            endif  
            if ( (1.D0-a.cur-b.cur).GT.(-prec) ) is_valid = .TRUE.
        enddo            
    end subroutine
    
!----------------------------------------------------------
    
    subroutine frame_kolm_test(G2)
        integer, dimension(:,:), intent(in) :: G2
        
        call kolm_test_3d(T,G2)
        call kolm_test_3d(Y,G2)
    end subroutine
     
    subroutine kolm_test_3d(T1, T2)
        integer, dimension(:,:), intent(in) :: T1, T2    !-- T1 and T2 should have equal sizes
        
        logical :: is_equal
        real(kind) :: d, sum_fail
        integer :: i , n_fail

        sum_fail = 0.D0
        n_fail = 0
        do i = 1, size(T1,1)
            write(43,"('cos theta = ',F5.3,\)") dfloat(i) * 2.D0 / size(T1,1)
            print    "('cos theta = ',F5.3,\)", dfloat(i) * 2.D0 / size(T1,1)
            if ( sum(T1(i,:)).LT.l_stat .OR. sum(T2(i,:)).LT.l_stat ) then
                write(43,"(' low statistic: ',2I)") sum(T1(i,:)), sum(T2(i,:))
                print    "(' low statistic: ',2I)", sum(T1(i,:)), sum(T2(i,:))
                n_fail = n_fail + 1
            else
                call kolm_test(T1(i,:), T2(i,:), is_equal, d)
                if ( is_equal ) then
                    write(43,"(L5)") is_equal
                    print    "(L5)", is_equal
                else
                    write(43,"(L5,F10.5)") is_equal, d
                    print    "(L5,F10.5)", is_equal, d 
                    n_fail = n_fail + 1
                    sum_fail = sum_fail + d
                endif
            endif
        enddo
        write(43,"('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5)") n_fail, size(T1,1), sum_fail
        print    "('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5)", n_fail, size(T1,1), sum_fail   
        do i = 1, size(T1,2)
            write(43,"('eps = ',F5.3,\)") dfloat(i) / size(T1,2)
            print    "('eps = ',F5.3,\)", dfloat(i) / size(T1,2)
            if ( sum(T1(:,i)).LT.l_stat .OR. sum(T2(:,i)).LT.l_stat ) then
                write(43,"(' low statistic: ',2I)") sum(T1(:,i)), sum(T2(:,i))
                print    "(' low statistic: ',2I)", sum(T1(:,i)), sum(T2(:,i))
                n_fail = n_fail + 1
            else
                call kolm_test(T1(:,i), T2(:,i), is_equal, d)
                if ( is_equal ) then
                    write(43,"(L5)") is_equal
                    print    "(L5)", is_equal
                else
                    write(43,"(L5,F10.5)") is_equal, d
                    print    "(L5,F10.5)", is_equal, d
                    n_fail = n_fail + 1
                    sum_fail = sum_fail + d
                endif
            endif
        enddo
        write(43,"('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5,/,50('-'))") n_fail, 2*size(T1,1), sum_fail
        print    "('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5,/,50('-'))", n_fail, 2*size(T1,1), sum_fail
      !  i_log = i_log + 1
        write(44,"(I5,' ',2(F5.3,' '),I5,' ',F10.5)") a.cur-a.step, b.cur-b.step, n_fail, sum_fail
    end subroutine 
        
    subroutine kolm_test(A1, A2, is_equal, alph)
        integer, dimension(:), intent(in) :: A1, A2  !-- A1 and A2 should have equal size
        logical, intent(out) :: is_equal
        real(kind), intent(out) :: alph
        
        integer :: i
        real(kind), dimension(size(A1)) :: X1, X2
        
        is_equal = .FALSE.
        alph = 0.D0
        X1 = dfloat(A1) / dfloat(sum(A1)) 
        X2 = dfloat(A2) / dfloat(sum(A2))
        do i = 2, size(A1)
            X1(i) = X1(i-1) + X1(i)
            X2(i) = X2(i-1) + X2(i)
        enddo
 !       alph = maxval(dabs(X1 - X2)) * sqrt( dfloat(sum(A1)) * dfloat(sum(A2)) / dfloat(sum(A1)+sum(A2)) )
        print *,'sum A = ', sum(A1)
        alph = maxval(dabs(X1 - X2)) * sqrt( dfloat(sum(A1)) )
        if ( alph.LT.alph_k ) then
            is_equal = .TRUE.
            alph = 0.D0
        else  
            alph = alph - alph_k
        endif   
    end subroutine

!== Xi^2 test =============================================
    subroutine xi_test(Tt, Yt)
        integer, dimension(:,:), intent(in) :: Tt, Yt
        
        real(kind) :: serv_k
        
        serv_k = xi_sq(T, Tt)
        write(43,"('T: xi-sq ',F6.3)") serv_k
        print    "('T: xi-sq ',F6.3)", serv_k
        write(44,"(F6.3,' ',\)") serv_k
        serv_k = xi_sq(Y, Yt)
        write(43,"('Y: xi-sq ',F6.3,/)") serv_k
        print    "('Y: xi-sq ',F6.3,/)", serv_k
        write(44,"(F6.3)") serv_k
    end subroutine
     
    subroutine xi_test_3d(T1, T2)
        integer, dimension(:,:), intent(in) :: T1, T2    !-- T1 and T2 should have equal sizes
        
      !  logical :: is_equal
        real(kind) :: d !, sum_fail
        integer :: i !, n_fail

    !    sum_fail = 0.D0
    !    n_fail = 0
        do i = 1, size(T1,1)
            write(43,"('cos theta = ',F5.3,\)") dfloat(i) * 2.D0 / size(T1,1)
            print    "('cos theta = ',F5.3,\)", dfloat(i) * 2.D0 / size(T1,1)
            if ( sum(T1(i,:)).LT.l_stat .OR. sum(T2(i,:)).LT.l_stat ) then
                write(43,"(' low statistic: ',2I)") sum(T1(i,:)), sum(T2(i,:))
                print    "(' low statistic: ',2I)", sum(T1(i,:)), sum(T2(i,:))
                
          !      n_fail = n_fail + 1
            else
   !!             call xi_test(T1(i,:), T2(i,:), d)
              !  if ( is_equal ) then
               !     write(43,"(L5)") is_equal
             !       print    "(L5)", is_equal
             !   else
                    write(43,"(F10.5)") d
                    print    "(F10.5)", d 
               !     n_fail = n_fail + 1
               !     sum_fail = sum_fail + d
              !  endif
            endif
        enddo
      !  write(43,"('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5)") n_fail, size(T1,1), sum_fail
     !   print    "('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5)", n_fail, size(T1,1), sum_fail   
        do i = 1, size(T1,2)
            write(43,"('eps = ',F5.3,\)") dfloat(i) / size(T1,2)
            print    "('eps = ',F5.3,\)", dfloat(i) / size(T1,2)
            if ( sum(T1(:,i)).LT.l_stat .OR. sum(T2(:,i)).LT.l_stat ) then
                write(43,"(' low statistic: ',2I)") sum(T1(:,i)), sum(T2(:,i))
                print    "(' low statistic: ',2I)", sum(T1(:,i)), sum(T2(:,i))
          !      n_fail = n_fail + 1
            else
       !!         call xi_test(T1(:,i), T2(:,i), d)
            !    if ( is_equal ) then
            !        write(43,"(L5)") is_equal
            !        print    "(L5)", is_equal
            !    else
                    write(43,"(F10.5)") d
                    print    "(F10.5)", d
             !       n_fail = n_fail + 1
            !        sum_fail = sum_fail + d
            !    endif
            endif
        enddo
     !   write(43,"('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5,/,50('-'))") n_fail, 2*size(T1,1), sum_fail
     !   print    "('number of fails: ',I3,'( from ',I3,' )',/,'summary fail: ',F10.5,/,50('-'))", n_fail, 2*size(T1,1), sum_fail
  !      i_log = i_log + 1
    !    write(44,"(I5,' ',2(F5.3,' '),I5,' ',F10.5)") i_log, a.cur-a.step, b.cur-b.step, n_fail, sum_fail
    end subroutine     

end module
