module generator
    use cons, only: kind, pi, dvapi, sq2, nping, prec, eps_min
    use types
    use amplitudes
    use other
    use output
    use output_ev
    use pp
    use ps
    use vars
    implicit none
       
!== private ===============================================
    type(sopart), private :: t, y1, y2   !-- calculated at "Events_gen" from - global core, part1, part2
    real(kind), private :: ex2       !-- 2*Etotal
!==========================================================
 
contains 

!----------------------------------------------------------
!-- array of events generation ----------------------------
!----------------------------------------------------------
!-- changed global: Evs                                  --
!-- changed private: t, y1, y2, ex2              --
    subroutine Events_gen()
   !     character(*), intent(in) :: name

        type(event) :: serv_ev   !-- for any needs
        real(kind) :: wwmax, &     !-- max of probability (majorating block)
                     g1, g2, &  !-- width of resonances
                     time, &    !-- time diagnostic
                     serv_k, &   !-- for any needs
                     v3    !-- coefficient v3 / 2  ( 2 comes from symmetric A )
        integer :: i, ii, &   !-- counter
                   serv_i  !-- for any needs
        logical :: serv_l   !-- for any needs
      
        y1 = set_system(core, part1, part2, e)  !-- set particles system "Y1"
        y2 = set_system(part2, core, part1, e)  !-- set particles system "Y2"
        t = set_system(part1, part2, core, e)   !-- set particles system "T"
        ex2 = e * 2.D0   !-- usefull constant
        
        write(5,"(/,A14,F10.3)") "E total, MeV: ", e
        print   "(/,A14,F10.3)", "E total, MeV: ", e
        do i = 1, size(Confs)
            g1 = gamm(Confs(i).r1, y1.x, Confs(i).r1.eng)  
            g2 = gamm(Confs(i).r2, y2.x, Confs(i).r2.eng)
            write(5, "(/,A2,I2,4(A10),2(/,A4,F10.3,ES10.2,2(F10.3)))") '# ', i, &
                'Er, MeV', 'Gr, MeV', 'Thet^2', 'Rch, fm', '', Confs(i).r1.eng, g1, Confs(i).r1.tet, &
                 Confs(i).r1.rch, '', Confs(i).r2.eng, g2, Confs(i).r2.tet, Confs(i).r2.rch   
            print    "(/,A2,I2,4(A10),2(/,A4,F10.3,ES10.2,2(F10.3)))", '# ', i, &
                'Er, MeV', 'Gr, MeV', 'Thet^2', 'Rch, fm', '', Confs(i).r1.eng, g1, Confs(i).r1.tet, &
                 Confs(i).r1.rch, '', Confs(i).r2.eng, g2, Confs(i).r2.tet, Confs(i).r2.rch       
            Confs(i).dE = e - Confs(i).r1.eng - Confs(i).r2.eng  !-- usefull constant   
            if ( Confs(i).r1.eng.LE.0.D0 .OR. Confs(i).r2.eng.LE.0.D0 ) then
                v3 = sqrt( e )
            else
                v3 = sqrt( (4.D0*Confs(i).dE*Confs(i).dE + g1*g1 + g2*g2) * e / pi ) / 4.D0      !-- coefficient v3 / sq2  ( sq2 comes from symmetric A 
            endif
            Confs(i).C_ls.c = Confs(i).C_ls.c * v3   !-- multiple recoup.coef to v3 for the same config.           
        enddo
!----------------------------------------------------------
! majorating (and diagnostic also)
!----------------------------------------------------------
        print "(/,20('-'),/,'majorating',/,20('-'),/)"
        time = dsecnds(0.D0)
        Evs = set_null(Evs(1))
        wwmax = 0.D0
        call random_seed()
        do i = 1, size(Evs)
            Evs(i) = event_gen()
            if ( Evs(i).ww.GT.wwmax ) then
                wwmax = Evs(i).ww
                print "(I,ES)", i, wwmax
            endif
        enddo
        write(5,"('max was found: ',ES,/,'it takes ',F10.4,' sekonds')") wwmax, dsecnds(time)
        print   "('max was found: ',ES,/,'it takes ',F10.4,' sekonds')", wwmax, dsecnds(time)
        serv_k = sum(Evs.ww) * 16.D0 * pi**2 / size(Evs)   !-- G total
        write(5,"(/,'Gtotal, MeV (',I8,' events):',ES)") size(Evs), serv_k
        print   "(/,'Gtotal, MeV (',I8,' events):',ES)", size(Evs), serv_k 
        if ( size(Confs).EQ.1 ) then
            serv_k = serv_k * ( 4.D0*Confs(1).dE*Confs(1).dE + 2.D0*serv_k*serv_k  ) / &
                                  ( 4.D0*Confs(1).dE*Confs(1).dE + g1*g1 + g2*g2 )
            write(5,"(/,'Gtotal, MeV ( 2nd step ):',ES)") serv_k
            print   "(/,'Gtotal, MeV ( 2nd step ):',ES)", serv_k             
        endif
        if ( do_out ) then
            call set_name_out(e, Confs(1).r1.eng, Confs(1).r1.tet)
            call set_is_ww(.TRUE.)
            call diagn_output()
        endif
!----------------------------------------------------------

        !-- generation ------------------------------------
        if (do_gen) then
            time = dsecnds(0.D0)
            print    "(/,A,I,A)", 'generation of ', size(Evs), ' events...'
            write(5, "(/,A,I,A)") 'generation of ', size(Evs), ' events...'
            i = 1
            serv_i = 0   !-- counter of all generation (for integrating)
            do while (i .LE. size(Evs))
                if ( mod(i, nping) .EQ. 0 ) then
                    write (5, "(I,A6)") i, 'ping!'
                endif
                serv_ev = Evs(i)   !-- try to use existed event from majorating
                serv_l = .FALSE.
                do while (.NOT. serv_l)
                    call random_number(serv_k)
                    if ( serv_k * wwmax .LT. serv_ev.ww ) then
                        Evs(i) = serv_ev
                        i = i + 1
                        serv_l = .TRUE.
                    else
                        serv_ev = event_gen()
                    endif
                    serv_i = serv_i + 1
                enddo
            enddo
            write (5, "(I,A)") serv_i, ' events was generated at all'
            serv_k = wwmax * 16.D0 * pi**2 * size(Evs) / serv_i
            print    "(/,A14,I8,A8,ES,/)", 'Gtotal, MeV (', serv_i, 'events) ', serv_k
            write(5, "(/,A,ES,/)") 'Gtotal, MeV ', serv_k
            if ( size(Confs).EQ.1 ) then
                serv_k = serv_k * ( 4.D0*Confs(1).dE*Confs(1).dE + 2.D0*serv_k*serv_k  ) / &
                                  ( 4.D0*Confs(1).dE*Confs(1).dE + g1*g1 + g2*g2 )
                write(5, "(/,A30,ES,/)") 'Gtotal, MeV ( 2nd step ) ', serv_k
                print    "(/,A30,ES,/)", 'Gtotal, MeV ( 2nd step ) ', serv_k               
            endif
            print "(A,F10.4,A)",  'it takes ', dsecnds(time), ' sekonds'
            call set_is_ww(.FALSE.)
            call diagn_output()
            call data_output()
        endif  
        do i = 1, size(Confs)
            g1 = gamm(Confs(i).r1, y1.x, Confs(i).r1.eng)  
            g2 = gamm(Confs(i).r2, y2.x, Confs(i).r2.eng)      
            Confs(i).dE = e - Confs(i).r1.eng - Confs(i).r2.eng  !-- usefull constant   
            if ( Confs(i).r1.eng.LE.0.D0 .OR. Confs(i).r2.eng.LE.0.D0 ) then
                v3 = sqrt( e )
            else
                v3 = sqrt( (4.D0*Confs(i).dE*Confs(i).dE + g1*g1 + g2*g2) * e / pi ) / 4.D0      !-- coefficient v3 / sq2  ( sq2 comes from symmetric A 
            endif
            Confs(i).C_ls.c = Confs(i).C_ls.c / v3   !-- divide recoup.coef to v3 for the same config.           
        enddo
    end subroutine
!----------------------------------------------------------
!-- one event generation ----------------------------------
    type(event) function event_gen()
        real(kind), dimension(6) :: rand ! = &
      !  (/0.249572452257228, 0.994163934376058, 0.640602593985953, 0.590190373438494, 0.280130785801967, 0.253168121687737/)    !-- 5 random parameters and 6th for Monte-Carlo
        logical :: fl
        real(kind) :: t_kx_sin_teta, t_ky_sin_teta, modp
        type(vector) :: px, py, v_rel

        event_gen = set_null(event_gen)
        event_gen.eng = e
        fl = .FALSE.
        do while (.NOT. fl)
            call random_number(rand)
            if ( ( rand(6)*0.5D0 .LT. sqrt(rand(1)*(1.D0-rand(1))) ) .AND. ( rand(1) .GT. eps_min ) ) then
                !-- system "T" ----------------------------
                event_gen.t.eps = rand(1)
                event_gen.t.kx.cos_teta = rand(2)*2.D0 - 1.D0
                event_gen.t.kx.phi = rand(3) * dvapi
                event_gen.t.ky.cos_teta = rand(4)*2.D0 - 1.D0
                event_gen.t.ky.phi = rand(5) * dvapi
                t_kx_sin_teta = sqrt(1.D0-event_gen.t.kx.cos_teta*event_gen.t.kx.cos_teta)
                t_ky_sin_teta = sqrt(1.D0-event_gen.t.ky.cos_teta*event_gen.t.ky.cos_teta)
                event_gen.t.cos_teta = t_kx_sin_teta * t_ky_sin_teta * cos(event_gen.t.kx.phi - event_gen.t.ky.phi) + &
                                       event_gen.t.kx.cos_teta * event_gen.t.ky.cos_teta
                modp = sqrt(t.x.mx2e * event_gen.t.eps)
                px = vector(modp*t_kx_sin_teta*cos(event_gen.t.kx.phi), modp*t_kx_sin_teta*sin(event_gen.t.kx.phi), modp*event_gen.t.kx.cos_teta)
                modp = sqrt(t.y.mx2e * (1.D0-event_gen.t.eps))
                py = vector(modp*t_ky_sin_teta*cos(event_gen.t.ky.phi), modp*t_ky_sin_teta*sin(event_gen.t.ky.phi), modp*event_gen.t.ky.cos_teta)
                if ( abs( abs(px)**2 / t.x.m + abs(py)**2 / t.y.m - ex2 ) .GT. prec ) then
                    print *, 'T px, py : delta energy > ', prec/2.D0, ( abs(px)**2 / 2.D0 / t.x.m + abs(py)**2 / 2.D0 / t.y.m - e )
                endif
                !-- "Center of Mass" system ---------------
                event_gen.core = -py
                event_gen.p1 = t.x.m_m2 * py + px
                event_gen.p2 = t.x.m_m1 * py - px
                if ( abs( abs(event_gen.core)**2 / core.m + abs(event_gen.p1)**2 / part1.m + abs(event_gen.p2)**2 / part2.m - ex2 ) .GT. prec ) then
                    print *, 'CM p1, p2, core : delta energy > ', prec/2.D0, ( abs(event_gen.core)**2 / 2.D0 / core.m + abs(event_gen.p1)**2 / 2.D0 / part1.m + abs(event_gen.p2)**2 / 2.D0 / part2.m - e )
                endif
                if ( abs(event_gen.p1 + event_gen.p2 + event_gen.core) .GT. prec*100 ) then
                    print *, 'CM not in null : delta |Pr| > ', prec, abs(event_gen.p1 + event_gen.p2 + event_gen.core)
                endif
                !-- "Y1" system ---------------------------
                event_gen.y1 = CmToJacobi(event_gen.core, event_gen.p1, event_gen.p2, y1, ex2, prec)
                !-- "Y2" system ---------------------------
                event_gen.y2 = CmToJacobi(event_gen.p2, event_gen.core, event_gen.p1, y2, ex2, prec)
                if ( ( event_gen.y1.eps .GT. eps_min ) .AND. ( event_gen.y2.eps .GT. eps_min ) ) then
                    call ww(event_gen)
                    v_rel = (1.D0 / part1.m) * event_gen.p1 - (1.D0 / core.m) * event_gen.core
                    event_gen.tet = 2.D0 * atan( sqrt(v_rel.x*v_rel.x + v_rel.y*v_rel.y) * v_beam )
                    fl = .TRUE.
                endif
            endif
        enddo
    end function
!----------------------------------------------------------
!-- probability calculating -------------------------------
    subroutine ww(ev)
        type(event), intent(inout) :: ev

        complex(kind) :: znam
        real(kind) :: ey1x, ey2x
        integer :: i

        ey1x = e * ev.y1.eps
        ey2x = e * ev.y2.eps
        do i = 1, size(Confs)
            znam = cmplx(Confs(i).dE, (gamm(Confs(i).r1, y1.x, ey1x) + gamm(Confs(i).r2, y2.x, ey2x))/2.D0)   
            if ( do_ps ) then
                Confs(i).a11 = ( Arp(Confs(i).r1, y1.x, ey1x) + Aps(Confs(i).r1, y1.x, ey1x) )* sqrt(gamm(Confs(i).r2, y1.y, e-ey1x)) / znam  !-- 1st resonance of 1st two-particles
                Confs(i).a22 = ( Arp(Confs(i).r2, y2.x, ey2x) + Aps(Confs(i).r2, y2.x, ey2x) ) * sqrt(gamm(Confs(i).r1, y2.y, e-ey2x)) / znam  !-- 2st resonance of 2st two-particles
                znam = cmplx(Confs(i).dE, (gamm(Confs(i).r2, y1.x, ey1x) + gamm(Confs(i).r1, y2.x, ey2x))/2.D0)
                Confs(i).a12 = ( Arp(Confs(i).r1, y2.x, ey2x) + Aps(Confs(i).r1, y2.x, ey2x) ) * sqrt(gamm(Confs(i).r2, y2.y, e-ey2x)) / znam  !-- 1st resonance of 2st two-particles
                Confs(i).a21 = ( Arp(Confs(i).r2, y1.x, ey1x) + Aps(Confs(i).r2, y1.x, ey1x) ) * sqrt(gamm(Confs(i).r1, y1.y, e-ey1x)) / znam !-- 2st resonance of 1st two-particles
            else
                Confs(i).a11 = Arp(Confs(i).r1, y1.x, ey1x) * sqrt(gamm(Confs(i).r2, y1.y, e-ey1x)) / znam     !-- 1st resonance of 1st two-particles                
                Confs(i).a22 = Arp(Confs(i).r2, y2.x, ey2x) * sqrt(gamm(Confs(i).r1, y2.y, e-ey2x)) / znam     !-- 2st resonance of 2st two-particles
                znam = cmplx(Confs(i).dE, (gamm(Confs(i).r2, y1.x, ey1x) + gamm(Confs(i).r1, y2.x, ey2x))/2.D0)
                Confs(i).a12 = Arp(Confs(i).r1, y2.x, ey2x) * sqrt(gamm(Confs(i).r2, y2.y, e-ey2x)) / znam    !-- 1st resonance of 2st two-particles
                Confs(i).a21 = Arp(Confs(i).r2, y1.x, ey1x) * sqrt(gamm(Confs(i).r1, y1.y, e-ey1x)) / znam    !-- 2st resonance of 1st two-particles
            endif
        enddo
        ev.ww_sym = half_ww(.TRUE.)
        ev.ww_asym = half_ww(.FALSE.) 
        if (do_pp) then
            ey1x = ev.t.eps *e 
            ev.ww_sym = get_Wpp(.TRUE., ey1x) * ev.ww_sym
            ev.ww_asym = get_Wpp(.FALSE., ey1x) * ev.ww_asym
        endif
        ev.ww = ev.ww_sym + ev.ww_asym
        contains

        real(kind) function half_ww(is_sym)
            logical, intent(in) :: is_sym

            integer :: i, ml, ii, s
            complex(kind) :: ww_g
            real(kind) :: ww_ml
            
            if ( is_sym ) s = 0
            if ( .NOT.is_sym ) s = 1
            half_ww = 0.D0
            do i = -1, 1   !-- L = J + i
                ww_ml = 0.D0
                do ml = -j-1, j+1   
                    ww_g = cmplx(0.D0, 0.D0)
                    do ii = 1, size(Confs)
                        if ( Confs(ii).C_ls(i, s).not_null .AND. Confs(ii).C_lml(i, ml).n.NE.0 ) then
                            ww_g = ww_g + Confs(ii).coef * Confs(ii).C_ls(i, s).c * a_sym(ii, is_sym, Confs(ii).C_lml(i, ml).ar)  
                        endif
                    enddo
                    ww_ml = ww_ml + cdabs(ww_g)**2    
                enddo
                half_ww = half_ww + ww_ml / (2*(j+i)+1)
            enddo 
        end function

        complex(kind) function a_sym(i, is_sym, Sub_clbs)
            integer, intent(in) :: i   !-- № of config.
            logical, intent(in) :: is_sym
            type(clgd), dimension(:), intent(in) :: Sub_clbs

            a_sym = lxl(ev.y1.kx, ev.y1.ky, Sub_clbs) * Confs(i).a11 + (-1.D0)**(Sub_clbs(1).l1)*lxl(ev.y2.ky, ev.y2.kx, Sub_clbs) * Confs(i).a22
            if (is_sym) then
                a_sym = a_sym + lxl(ort(-ev.y2.kx.cos_teta, ev.y2.kx.phi+pi), ev.y2.ky, Sub_clbs) * Confs(i).a12 + &
                                (-1.D0)**(Sub_clbs(1).l1)*lxl(ev.y1.ky, ort(-ev.y1.kx.cos_teta, ev.y1.kx.phi+pi), Sub_clbs) * Confs(i).a21 
            else
                a_sym = a_sym - lxl(ort(-ev.y2.kx.cos_teta, ev.y2.kx.phi+pi), ev.y2.ky, Sub_clbs) * Confs(i).a12 - &
                                (-1.D0)**(Sub_clbs(1).l1)*lxl(ev.y1.ky, ort(-ev.y1.kx.cos_teta, ev.y1.kx.phi+pi), Sub_clbs) * Confs(i).a21 
            endif
        end function

    end subroutine
!----------------------------------------------------------
end module
