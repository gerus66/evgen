module amplitudes
    use mydata, only: kind, fm_mev, alpha_const, pi, dvapi
    use types
    use other
    implicit none

    contains
!==========================================================
!-- width of single resonance -----------------------------
    real(kind) function gamm(r, tpart, e)
        type(res), intent(in) :: r
        type(twopart), intent(in) :: tpart
        real(kind), intent(in) ::  e 

        gamm = fm_mev*fm_mev * r.tet / (tpart.m * r.rch*r.rch) * pen(r.rch, r.l, tpart, e)   !--- fm_mev - constant from module mydata.f90
    endfunction
!----------------------------------------------------------
!-- penetrability -----------------------------------------
    real(kind) function pen(rrch, rl, tpart, e)
        real(kind), intent(in) :: rrch
        integer, intent(in) :: rl
        type(twopart), intent(in) :: tpart
        real(kind), intent(in) ::  e 
        
        real(kind) :: ekk, rho
        integer :: ifail, l
        real(kind), dimension(1:30) :: FC, GC, FCP, GCP 
  
        ekk = dsqrt( 2.D0 * tpart.m * e )
        rho = rrch * ekk / fm_mev        !--- fm_mev - constant from module mydata.f90
        l = rl + 1 
        if ( (tpart.z1 .EQ. 0) .OR. ( (tpart.z2) .EQ. 0) ) then
            call coulfg(rho, 0.D0, dfloat(rl), dfloat(rl), FC, GC, FCP, GCP, 1, 1, ifail)
            pen = 1.D0 / (FC(l)**2 + GC(l)**2) / rho
        else
            call coulfg(rho, dfloat(tpart.z1)*dfloat(tpart.z2)*tpart.m*alpha_const/ekk, dfloat(rl), dfloat(rl), FC, GC, FCP, GCP, 1, 0, ifail)       !--- alpha_const - constant from module mydata.f90
            pen = rho / (FC(l)**2 + GC(l)**2)
        endif    	
    endfunction
!----------------------------------------------------------
!-- amplitude of "r"-resonance of "p"-two-particles -------
    complex(kind) function Arp(r, tpart, e)
        type(res), intent(in) :: r
        type(twopart), intent(in) :: tpart
        real(kind), intent(in) :: e
    
        real(kind) :: g

        g = gamm(r, tpart, e)
        Arp = sqrt(g) / cmplx( r.eng-e, -g/2 )
    endfunction
!----------------------------------------------------------
!-- [l1 x l2] -> L,Ml -------------------------------------
    complex(kind) function lxl(kx, ky, clbs)
        type(ort), intent(in) :: kx, ky
        type(clgd), dimension(:), intent(in) :: clbs

        integer :: i

        lxl = cmplx(0.D0, 0.D0)
        do i = 1, size(clbs)
            lxl = lxl + clbs(i).c * Ylm(clbs(i).l1, clbs(i).m1, kx) * Ylm(clbs(i).l2, clbs(i).m2, ky)
        enddo
    endfunction
!----------------------------------------------------------
!-- Amplitude^2 of p-p interaction ------------------------
!    real(kind) function Wpp(is_sym, tpart, e, r_pp)
!        logical, intent(in) :: is_sym
!        type(twopart), intent(in) :: tpart
!        real(kind), intent(in) ::  e, r_pp
!
!        if (is_sym) then
!            Wpp = 4.D0 / pen(r_pp, 0, tpart, e)
!        else
!            Wpp = pen(r_pp, 1, tpart, e) * fm_mev / r_pp / sqrt( 2.D0 * e * tpart.m )
!        endif
!    endfunction
!----------------------------------------------------------
!-- Sin^2(delta_0) ----------------------------------------
!    real(kind) function sind0(tpart, e)
!        type(twopart), intent(in) :: tpart
!        real(kind), intent(in) :: e

!        sind0 = dvapi * tpart.z1 * tpart.z2 * alpha_const * tpart.m * (1.D0 + 8.166*e) / &                      !-- tg(delta0)
!                ( exp( pi * tpart.z1 * tpart.z2 * alpha_const * sqrt( 2.D0 * tpart.m / e ) ) - 1.D0 ) / &
!                ( 46.05D0*e*e + 113.9*e + 25.21 )
!        sind0 = sind0 * sind0         !-- tg^2(delta0)
!        sind0 = sind0 / ( sind0 + 1.D0 )         !-- sin^2(delta0)
!    endfunction
!----------------------------------------------------------
!-- Amplitude of potential scattering ---------------------
    complex(kind) function Aps(r_ps, r, tpart, e)
        type(res), intent(in) :: r
        type(twopart), intent(in) :: tpart
        real(kind), intent(in) ::  e, r_ps

        real(kind) :: ekk, rho, pen, gf
        integer :: ifail, l
        real(kind), dimension(1:30) :: FC, GC, FCP, GCP

        ekk = dsqrt( 2.D0 * tpart.m * e )
        rho = r_ps * ekk / fm_mev        !--- fm_mev - constant from module mydata.f90
        l = r.l + 1 
        if ( (tpart.z1 .EQ. 0) .OR. ( (tpart.z2) .EQ. 0) ) then
            call coulfg(rho, 0.D0, dfloat(r.l), dfloat(r.l), FC, GC, FCP, GCP, 1, 1, ifail)
            pen = 1.D0 / (FC(l)**2 + GC(l)**2) / rho
        else
            call coulfg(rho, dfloat(tpart.z1)*dfloat(tpart.z2)*tpart.m*alpha_const/ekk, dfloat(r.l), dfloat(r.l), FC, GC, FCP, GCP, 1, 0, ifail)       !--- alpha_const - constant from module mydata.f90
            pen = rho / (FC(l)**2 + GC(l)**2)
        endif  
        gf = GC(l) / FC(l)
        Aps = - 2.D0 * r_ps / fm_mev * sqrt( tpart.m/(r.tet*pen) ) / ( 1.D0 + gf*gf ) * cmplx(gf,-1.D0)
    endfunction
!----------------------------------------------------------
endmodule