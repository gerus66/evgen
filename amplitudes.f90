module amplitudes
    use cons, only: kind, fm_mev, alpha_const, pi, dvapi
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
        if ( r.eng.LE.0.D0 ) then
            Arp = sqrt(g) / cmplx( -e, r.eng )
        else
            Arp = sqrt(g) / cmplx( r.eng-e, -g/2 )
        endif
    end function
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
    end function
!----------------------------------------------------------

end module
