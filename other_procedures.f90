module other
    use cons, only: kind, pi, sqpi, sq3, sqdvapi
    use types
    implicit none

    contains
!==========================================================
!-- X : m2->m1, Y : m3->X ---------------------------------
    type(correlation) function CmToJacobi(p1, p2, p3, sp, ex2, prec) 
        real(kind), intent(in) :: ex2, prec    
        type(vector), intent(in) :: p1, p3, p2   
        type(sopart), intent(in) :: sp
                                          
        type(vector) :: px, py
        real(kind) :: abspx, abspy

        CmToJacobi = set_null(CmToJacobi)        
        px = sp.x.m_m1 * p1 - sp.x.m_m2 * p2 
        py = sp.y.m_m1 * p1 + sp.y.m_m1 * p2 - sp.y.m_m2 * p3 
        abspx = abs(px)
        abspy = abs(py) 
        if ( ( abspx*abspx / sp.x.m + abspy*abspy / sp.y.m - ex2 ) .GT. prec ) then
            print *, 'Y px, py : delta energy > ', prec/2.D0, ( abspx*abspx / sp.x.m / 2.D0 + abspy*abspy / sp.y.m / 2.D0 - ex2/2.D0 )
        endif
        CmToJacobi.eps = abspx*abspx / sp.x.mx2e
        CmToJacobi.cos_teta = dot(px,py) / abspx / abspy
        CmToJacobi.kx.cos_teta = px.z / abspx
        CmToJacobi.kx.phi = atan2(px.y, px.x) + pi
        CmToJacobi.ky.cos_teta = py.z / abspy
        CmToJacobi.ky.phi = atan2(py.y, py.x) + pi
    endfunction
!----------------------------------------------------------
!-- Spherical Ylm(cos_teta, phi) for l = {0, .., 2} -------
    complex(kind) function Ylm(l, m, k)    
        integer, intent(in) :: l, m
        type(ort), intent(in) :: k

        select case (l)
            case (0)
                Ylm = cmplx( 0.5D0 / sqpi, 0.D0 )
            case (1)
                select case (m)
                    case (-1)
                        Ylm = sq3 / 2.D0 / sqdvapi * sqrt(1.D0-k.cos_teta*k.cos_teta) * cmplx(cos(k.phi), -sin(k.phi))
                    case (0)
                        Ylm = cmplx( sq3 / 2.D0 / sqpi * k.cos_teta, 0.D0 )
                    case (1)
                        Ylm = -sq3 / 2.D0 / sqdvapi * sqrt(1.D0-k.cos_teta*k.cos_teta) * cmplx(cos(k.phi), sin(k.phi))
                    case default
                        print *, 'Ylm not defined'
                endselect
            case (2)
                select case (m)
                    case (-2)
                        Ylm = sqrt(15.D0) / 4.D0 / sqdvapi * (1.D0-k.cos_teta*k.cos_teta) * cmplx(cos(2.D0*k.phi), -sin(2.D0*k.phi))
                    case (-1)
                        Ylm = sqrt(15.D0) / 2.D0 / sqdvapi * sqrt(1.D0-k.cos_teta*k.cos_teta) * k.cos_teta * cmplx(cos(k.phi), -sin(k.phi))
                    case (0)
                        Ylm = cmplx( sqrt(5.D0 / pi) / 4.D0 * ( 3.D0*k.cos_teta*k.cos_teta - 1.D0 ), 0.D0 )
                    case (1)
                        Ylm = -sqrt(15.D0) / 2.D0 / sqdvapi * sqrt(1.D0-k.cos_teta*k.cos_teta) * k.cos_teta * cmplx(cos(k.phi), sin(k.phi))
                    case (2)
                        Ylm = -sqrt(15.D0) / 4.D0 / sqdvapi * (1.D0-k.cos_teta*k.cos_teta) * cmplx(cos(2.D0*k.phi), sin(2.D0*k.phi))
                    case default
                        print *, 'Ylm not defined'
                endselect
            case default
                print *, 'Ylm not defined'
            endselect
    endfunction
!----------------------------------------------------------
endmodule
