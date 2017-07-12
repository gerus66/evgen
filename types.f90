module types
    use cons, only: kind, proton_mass, neutron_mass, sq3
    implicit none

!==========================================================
!-- particle ----------------------------------------------
    type part         
        real(kind) :: m
        integer :: z, a
    endtype
!----------------------------------------------------------
!-- two particles -----------------------------------------
    type twopart
        real(kind) :: m, m_m1, m_m2, mx2e   !-- m2 -> m1
        integer :: z1, z2
    endtype
!----------------------------------------------------------
!-- system of particles -----------------------------------
    type sopart
        type(twopart) :: x, y
    endtype
!----------------------------------------------------------
!-- vector in decart coordinates --------------------------
    type vector
        real(kind) :: x, y, z
    endtype
!----------------------------------------------------------
!-- vector in spherical coordinates -----------------------
    type ort
        real(kind) :: cos_teta, phi
    endtype
!----------------------------------------------------------
!-- correlation parameters of particle system -------------
    type correlation
        real(kind) :: cos_teta, eps
        type(ort) :: kx, ky
    endtype
!----------------------------------------------------------
!-- event -------------------------------------------------
    type event
        real(kind) :: eng, ww, ww_sym, ww_asym, tet
        type(correlation) :: y1, y2, t  
        type(vector) :: p1, p2, core 
    endtype
!----------------------------------------------------------
!-- Clebsh-Gordan coefficients ----------------------------
    type clgd
        real(kind) :: c
        integer :: l1, m1, l2, m2, L, Ml
    endtype
!----------------------------------------------------------
!-- array of Cl-Gd coefs with defined L, Ml ---------------
    type ar_clgd
        integer :: n, l, ml  !-- size of ar
        type(clgd), dimension(:), allocatable :: ar 
    endtype
!----------------------------------------------------------
!-- recouple coefficients {j1,l1,j2,l2} -> {J,L,S} --------
    type recoup
        logical :: not_null
        real(kind) :: c
        integer :: j1x2, l1, j2x2, l2, j, l, s
 !       type(subar) :: sub
    endtype
!----------------------------------------------------------
!-- resonance ---------------------------------------------
    type res
        real(kind) :: eng, gam, tet, rch
        integer :: l, jx2
    endtype
!----------------------------------------------------------
!-- configuration -----------------------------------------
    type conf
        real(kind) :: dE, g3
        complex(kind) :: coef, a11, a22, a12, a21
        type(res) :: r1, r2
        type(ar_clgd), dimension(:,:), allocatable :: C_lml   !-- Clebsh-Gordan coefficients
   !-- structure of C_lml:   [-1:1, -j-j:j+1]
   !-- L \ Ml  -J-1     -J   ....   J+1
   !--   J-1   False   False  ..   False
   !--    J    False    [x]   ..   False
   !--   J+1     [x]    [x]   ..    [x]
        type(recoup), dimension(-1:1,0:1) :: C_ls !-- recouple coefficients 
   !-- structure of C_ls:
   !--  L \ S     0      1
   !--   J-1   False     x
   !--    J       x      x
   !--   J+1   False     x
    end type
!----------------------------------------------------------

!==========================================================
!----------------------------------------------------------
    interface set_null
        module procedure set_null_part
        module procedure set_null_twopart
        module procedure set_null_sopart
        module procedure set_null_res
        module procedure set_null_vec
        module procedure set_null_ort
        module procedure set_null_corr
        module procedure set_null_event
        module procedure set_null_clgd
        module procedure set_null_ar_clgd
        module procedure set_null_recoup
    endinterface
!----------------------------------------------------------
    interface abs
        module procedure modul
    endinterface
!----------------------------------------------------------
    interface operator(-)
        module procedure minp
        module procedure revp
    endinterface
!----------------------------------------------------------
    interface operator(+)
        module procedure sump
    endinterface
!----------------------------------------------------------
    interface operator(*)
        module procedure multp
    endinterface
!----------------------------------------------------------
    interface make_all
        module procedure make_all_clgd
        module procedure make_all_recoup
    endinterface
!----------------------------------------------------------

    contains

!== SET_NULL block ========================================
!-- create null particle ----------------------------------
    type(part) function set_null_part(p)
        type(part), intent(in) :: p

        set_null_part = part(0.D0, 0, 0)
    endfunction
!----------------------------------------------------------
!-- create null tow particles -----------------------------
   type(twopart) function set_null_twopart(tp)
        type(twopart), intent(in) :: tp

        set_null_twopart = twopart(0.D0, 0.D0, 0.D0, 0.D0, 0, 0)
    endfunction
!----------------------------------------------------------
!-- create null system of particles -----------------------
    type(sopart) function set_null_sopart(sp)
        type(sopart), intent(in) :: sp

        type(twopart) :: tp

        set_null_sopart = sopart(set_null(tp), set_null(tp))
    endfunction
!----------------------------------------------------------
!-- create null resonance ---------------------------------
    type(res) function set_null_res(r)
        type(res), intent(in) :: r

        set_null_res = res(0.D0, 0.D0, 0.D0, 0.D0, 0, 0)
    endfunction
!----------------------------------------------------------
!-- create null vector ------------------------------------
    type(vector) function set_null_vec(v)
        type(vector), intent(in) :: v

        set_null_vec = vector(0.D0, 0.D0, 0.D0)
    endfunction
!----------------------------------------------------------
!-- create null spherical vector --------------------------
    type(ort) function set_null_ort(o)
        type(ort), intent(in) :: o

        set_null_ort = ort(0.D0, 0.D0)
    endfunction
!----------------------------------------------------------
!-- create null correlation -------------------------------
    type(correlation) function set_null_corr(corr)
        type(correlation), intent(in) :: corr

        type(ort) :: o

        set_null_corr = correlation(0.D0, 0.D0, set_null(o), set_null(o))
    endfunction
!----------------------------------------------------------
!-- create null event -------------------------------------
    type(event) function set_null_event(ev)
        type(event), intent(in) :: ev

        type(correlation):: corr
        type(vector) :: v

        set_null_event = event(0.D0, 0.D0, 0.D0, 0.D0, 0.D0, set_null(corr), set_null(corr), set_null(corr), set_null(v), set_null(v), set_null(v))
    endfunction
!----------------------------------------------------------
!-- create null Clebsh-Gordan coefficient -----------------
    type(clgd) function set_null_clgd(cg)
        type(clgd), intent(in) :: cg

        set_null_clgd = clgd(0.D0, 0, 0, 0, 0, 0, 0)
    endfunction
!----------------------------------------------------------
!-- create null array of Clebsh-Gordan coeffs -------------
    type(ar_clgd) function set_null_ar_clgd(ar_cg)
        type(ar_clgd), intent(in) :: ar_cg

        set_null_ar_clgd.n = 0
        set_null_ar_clgd.l = 0
        set_null_ar_clgd.ml = 0
        allocate(set_null_ar_clgd.ar(set_null_ar_clgd.n))
    end function
!----------------------------------------------------------
!-- create null recouple coefficient ----------------------
    type(recoup) function set_null_recoup(rc)
        type(recoup), intent(in) :: rc

        set_null_recoup = recoup(.FALSE.,0.D0, 0, 0, 0, 0, 0, 0, 0)
    endfunction
!----------------------------------------------------------

!==========================================================

!-- set mass by Z and A numbers ---------------------------
    type(part) function setmass(p)
        type(part), intent(in) :: p
        
        setmass = p
        if ( p.a .LT. p.z ) then
            print *, 'mass cant be negative'
        else
            setmass.m = p.z*proton_mass + (p.a-p.z)*neutron_mass 
        endif
    endfunction
!----------------------------------------------------------
!-- set X-part from {p1, p2} and Y-part from {p1+p2, p3} --
    type(sopart) function set_system (p1, p2, p3, e)
        type(part), intent(in) :: p1, p2, p3
        real(kind), intent(in) :: e

        real(kind) :: m

        m = p1.m*p2.m/(p1.m+p2.m)
        set_system.x = twopart(m, m/p1.m, m/p2.m, 2.D0*e*m, p1.z, p2.z)
        m = (p1.m+p2.m)*p3.m/(p1.m+p2.m+p3.m)
        set_system.y = twopart(m, m/(p1.m+p2.m), m/p3.m, 2.D0*e*m, p1.z+p2.z, p3.z)
    endfunction
!----------------------------------------------------------
!-- increase vector P in A times --------------------------
    type(vector) function multp(a, p)
        real(kind), intent(in) :: a
        type(vector), intent(in) :: p

        multp = vector(a*p.x, a*p.y, a*p.z)
    endfunction
!----------------------------------------------------------
!-- vector P1 + vector P2 ---------------------------------
    type(vector) function sump(p1, p2)
        type(vector), intent(in) :: p1, p2

        sump = vector(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
    endfunction
!----------------------------------------------------------
!-- vector P1 - vector P2 ---------------------------------
    type(vector) function minp(p1, p2)
        type(vector), intent(in) :: p1, p2

        minp = vector(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
    endfunction
!----------------------------------------------------------
!--  - vector P -------------------------------------------
    type(vector) function revp(p)
        type(vector), intent(in) :: p

        revp = vector(-p.x, -p.y, -p.z)
    endfunction
!----------------------------------------------------------
!-- | P | -------------------------------------------------
    elemental real(kind) function modul(p)
        type(vector), intent(in) :: p

        modul = sqrt(p.x*p.x + p.y*p.y + p.z*p.z)
    endfunction
!----------------------------------------------------------
!-- (P1, P2) ----------------------------------------------
    real(kind) function dot(p1, p2)
        type(vector), intent(in) :: p1, p2

        dot = p1.x*p2.x + p1.y*p2.y + p1.z*p2.z
    endfunction
!----------------------------------------------------------
!-- make all Cl-Gd coefficients for the (l1,l2)->J,L-------
    subroutine make_all_clgd(l1, l2, j, l, Clbs)
        integer, intent(in) :: l1, l2, j, l
        type(ar_clgd), dimension(-j-1:j+1), intent(out) :: Clbs

        integer :: i, m1, m2, ml
        type(clgd) :: cl
        real(kind) :: clb
        
        Clbs = set_null(Clbs(-j-1))
        do ml = -j-1, j+1 
            if ( abs(ml).LE.l ) then
                Clbs(ml).l = l; Clbs(ml).ml = ml; Clbs(ml).n = 0
                m1 = -l1
                do while (m1 .LE. l1)   
                    m2 = ml - m1
                    if ( abs(m2) .LE. l2 ) then  
                        cl = clgd(clb(dfloat(l1), dfloat(m1), dfloat(l2), dfloat(m2), dfloat(L), dfloat(Ml)), l1, m1, l2, m2, L, Ml) 
                        if ( dabs(cl.c) .GT. 1.D-15 ) Clbs(ml).n = Clbs(ml).n + 1
                    endif
                    m1 = m1 + 1
                enddo  
                if ( Clbs(ml).n.NE.0 ) then
                    deallocate(Clbs(ml).ar)
                    allocate(Clbs(ml).ar(Clbs(ml).n)) 
                    m1 = -l1
                    i = 1
                    do while (m1 .LE. l1)   
                        m2 = ml - m1
                        if ( abs(m2) .LE. l2 ) then  
                            cl = clgd(clb(dfloat(l1), dfloat(m1), dfloat(l2), dfloat(m2), dfloat(L), dfloat(Ml)), l1, m1, l2, m2, L, Ml) 
                            if ( dabs(cl.c) .GT. 1.D-15 ) then
                                Clbs(ml).ar(i) = cl
                                i = i + 1
                            endif
                        endif
                        m1 = m1 + 1
                    enddo 
                endif                 
            endif
        enddo                   
    end subroutine
!----------------------------------------------------------
!-- make all recouple coefs of (j1x2,l1,j2x2,l2)->J,S -----
    subroutine make_all_recoup(j1x2, l1, j2x2, l2, j, s, Rcs)
        integer, intent(in) :: j1x2, l1, j2x2, l2, j, s
        type(recoup), dimension(-1:1), intent(out) :: Rcs

        integer :: i
        real(kind) :: s9j, h, hj1, hj2, j1, j2
        type(recoup) :: rc
        
        Rcs = set_null(rc)
        j1 = dfloat(j1x2) / 2.D0; j2 = dfloat(j2x2) / 2.D0
        hj1 = h(j1); hj2 = h(j2)
        do i = -1, 1    !-- l = j + i
            if (s.EQ.0 .AND. i.EQ.0 .OR. j.GE.-i) then
                Rcs(i) = recoup( .TRUE., h(dfloat(j+i)) * h(dfloat(s)) * hj1 * hj2 * s9j(dfloat(l1), dfloat(l2), dfloat(j+i), 0.5D0, 0.5D0, dfloat(s), j1, j2, dfloat(j)), j1x2, l1, j2x2, l2, j, j+i, s )
                if ( dabs(Rcs(i).c) .LT. 1.D-15 ) Rcs(i) = set_null(rc)
            endif          
        enddo
    end subroutine
!----------------------------------------------------------
end module
