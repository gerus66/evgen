!-- potential scattering 
!-- requied: ps.txt - config. file in current directory
!--          coulfg4.for - subroutine 'coulfg'
!-- used units: 31 - ps_read()

module ps
    use cons, only: kind, fm_mev, alpha_const
    use types, only: res, twopart
    implicit none

!== private ===============================================
    real(kind), private :: rps = 0.D0  !-- radius
        
    contains
    
    subroutine ps_read()
        open(unit=31,file='ps.txt')
        read(31,*) rps
        close(31)
    end subroutine
        
    complex(kind) function Aps(r, tpart, e)
        type(res), intent(in) :: r
        type(twopart), intent(in) :: tpart
        real(kind), intent(in) ::  e

        real(kind) :: ekk, rho, pen, gf
        integer :: ifail, l
        real(kind), dimension(1:30) :: FC, GC, FCP, GCP

        ekk = dsqrt( 2.D0 * tpart.m * e )
        rho = rps * ekk / fm_mev 
        l = r.l + 1 
        if ( (tpart.z1 .EQ. 0) .OR. ( (tpart.z2) .EQ. 0) ) then
            call coulfg(rho, 0.D0, dfloat(r.l), dfloat(r.l), FC, GC, FCP, GCP, 1, 1, ifail)
            pen = 1.D0 / (FC(l)**2 + GC(l)**2) / rho
        else
            call coulfg(rho, dfloat(tpart.z1)*dfloat(tpart.z2)*tpart.m*alpha_const/ekk, dfloat(r.l), dfloat(r.l), &
                        FC, GC, FCP, GCP, 1, 0, ifail)      
            pen = rho / (FC(l)**2 + GC(l)**2)
        endif  
        gf = GC(l) / FC(l)
        Aps = - 2.D0 * rps / fm_mev * sqrt( tpart.m/(r.tet*pen) ) / ( 1.D0 + gf*gf ) * cmplx(gf,-1.D0)
    end function

end module
