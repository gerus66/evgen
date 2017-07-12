module types
    use mydata, only: kind, proton_mass, neutron_mass, sq3
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
!-- start and finish indexes for subarray -----------------
    type subar
        integer :: i_st, i_fin
    endtype
!----------------------------------------------------------
!-- recouple coefficients {j1,l1,j2,l2} -> {J,L,S} --------
    type recoup
        real(kind) :: c
        integer :: j1x2, l1, j2x2, l2, j, l, s
        type(subar) :: sub
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
        type(res) :: r1, r2
        type(clgd), dimension(:), allocatable :: Clbs
    endtype
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
        module procedure set_null_subar
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
!-- create null subarray ----------------------------------
    type(subar) function set_null_subar(sa)
        type(subar), intent(in) :: sa

        set_null_subar = subar(0, 0)
    endfunction
!----------------------------------------------------------
!-- create null recouple coefficient ----------------------
    type(recoup) function set_null_recoup(rc)
        type(recoup), intent(in) :: rc

        type(subar) :: sb

        set_null_recoup = recoup(0.D0, 0, 0, 0, 0, 0, 0, 0, set_null(sb))
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
!-- make all Cl-Gd coefficients for the (l1,l2)->J --------
    subroutine make_all_clgd(l1, l2, j, Clbs, count_null)
        integer, intent(in) :: l1, l2, j
        type(clgd), dimension(:), intent(out) :: Clbs
        integer, intent(out) :: count_null

        integer :: i, m1, m2, l, ml
        type(clgd) :: cl
        real(kind) :: clb

        i = 1
        if (j .EQ. 0) then
            l = 0
        else
            l = j - 1
        endif
        do while (l .LE. j+1)
            ml = -l
            do while (ml .LE. l)
                m1 = -l1
                do while (m1 .LE. l1)   
                    m2 = ml - m1
                    if ( abs(m2) .LE. l2 ) then  
                        clbs(i) = clgd(clb(dfloat(l1), dfloat(m1), dfloat(l2), dfloat(m2), dfloat(L), dfloat(Ml)), l1, m1, l2, m2, L, Ml) 
                        if ( dabs(clbs(i).c) .LT. 1.D-15 ) then
                            write(5,"(4(I5),A12)") l, ml, m1, m2, 'equal null'
                            print   "(4(I5),A12)", l, ml, m1, m2, 'equal null'
                        else
                            i = i + 1
                        endif 
                    else
                        write(5,"(4(I5),A12)") l, ml, m1, m2, 'not exist'
                        print   "(4(I5),A12)", l, ml, m1, m2, 'not exist'
                    endif
                    m1 = m1 + 1
                enddo
                ml = ml + 1
            enddo
            l = l + 1
        enddo
        count_null = size(clbs)-i+1
        do i=size(clbs)-count_null+1, size(clbs)
            clbs(i) = set_null(cl)
        enddo
    endsubroutine
!----------------------------------------------------------
!-- make all recouple coefs of (j1x2,l1,j2x2,l2)->J,S -----
    subroutine make_all_recoup(j1x2, l1, j2x2, l2, j, s, clbs, rcs, count_null)
        integer, intent(in) :: j1x2, l1, j2x2, l2, j, s
        type(clgd), dimension(:), intent(in) :: clbs
        type(recoup), dimension(:), intent(out) :: rcs
        integer, intent(out) :: count_null

        integer :: i, l
        real(kind) :: s9j, h, hj1, hj2, j1, j2
        type(recoup) :: rc

        j1 = dfloat(j1x2) / 2.D0; j2 = dfloat(j2x2) / 2.D0
        hj1 = h(j1); hj2 = h(j2)
        i = 1
        if (j .EQ. 0) then
            l = 0
        else
            l = j - s
        endif
        do while (l .LE. j+s)
            rcs(i) = recoup( h(dfloat(l)) * h(dfloat(s)) * hj1 * hj2 * s9j(dfloat(l1), dfloat(l2), dfloat(l), 0.5D0, 0.5D0, dfloat(s), j1, j2, dfloat(j)), &
                              j1x2, l1, j2x2, l2, j, l, s, make_sub(clbs, l) )
            if ( dabs(rcs(i).c) .LT. 1.D-15 ) then
                print "(2(I5),A12)", l, s, 'equal null'
            else
                i = i + 1
            endif 
            l = l + 1
        enddo
        count_null = size(rcs)-i+1
        do i=size(rcs)-count_null+1, size(rcs)
            rcs(i) = set_null(rc)
        enddo
    endsubroutine
!----------------------------------------------------------
!-- make subarray of Cl-Gd coefficients with the L --------
    type(subar) function make_sub(clbs, l)
        type(clgd), dimension(:), intent(in) :: clbs
        integer, intent(in) :: l

        integer :: i

        make_sub = set_null(make_sub)
        do i = 1, size(clbs)
            if ( clbs(i).l .EQ. l ) then
                if ( make_sub.i_st .EQ. 0 ) then
                    make_sub.i_st = i
                    make_sub.i_fin = i
                else
                    make_sub.i_fin = i
                endif
            endif
        enddo
    endfunction
!----------------------------------------------------------
endmodule