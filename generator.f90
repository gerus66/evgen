module generator
    use mydata, only: kind, pi, dvapi, sq2
    use types
    use amplitudes
    use other
    use output
    implicit none
!== global ================================================
    type(event), dimension(:), allocatable :: Evs   !-- all events
    type(conf), dimension(:), allocatable :: Confs  !-- all configuration
  !  type(clgd), dimension(:), allocatable :: Clbs   !-- Clebsh-Gordan coefficients (all notnull)
    type(recoup), dimension(:), allocatable :: Cjl0, Cjl1  !-- recouple coefficients (all notnull): S=0, S=1
    real(kind), dimension(:,:), allocatable :: Pp0, Pp1 !-- proton-proton interaction (f(xi), a, b: f(x)=ax+b)
    type(part) :: core, part1, part2  !-- physical particles
 !   type(res), dimension(:), allocatable :: R1, R2   !-- resonances
    real(kind) :: e, eps_min, prec, v_beam, step_pp, r_ps   !-- total energy; min epsilon for "coulfg" dont produce trash; precision of check; 1/2Vbeam; p-p energy step in file; potential scatering radius
    integer :: nping, j    !-- number events for ping; total J
    logical :: do_grid, do_gen, do_pp, do_ps
    character(50) :: pp_0, pp_1  !-- names of p-p data file for extrapolation
    
!== private ===============================================
    type(sopart), private :: t, y1, y2   !-- calculated at "Events_gen" from - global core, part1, part2
    real(kind), private :: v3, dE, ex2       !-- calculated at "Events_gen" from - global e, r1, r2
!==========================================================
 
contains 
!-- some preparation for pp-interaction -------------------
    subroutine pp_read()
        integer :: i, serv_i
        real(kind) :: serv_k

        open(unit=1, file=pp_0)
        read(1,*) serv_i
        read(1,*) step_pp
        i = int(e/step_pp)+1
        if ( serv_i .LT. i ) stop 'pp-files have not enough length'    
        allocate(Pp0(i,3), Pp1(i,3))  !-- read from files up to E_total
        do i=1, size(Pp0,1)
            read(1,*) Pp0(i,1)
        enddo
        close(1)
        open(unit=1, file=pp_1)
        read(1,*) serv_i; 
        if ( serv_i .LT. size(Pp1,1) ) stop 'pp-files have not enough length' 
        read(1,*) serv_k
        if ( abs(serv_k-step_pp) .GT. prec ) stop 'we need equal energy-step for 0 and 1'
        do i=1, size(Pp1,1)
            read(1,*) Pp1(i,1)
        enddo
        close(1)
        serv_k = 0.D0
        do i=1, size(Pp0,1)-1
            serv_k = serv_k + Pp0(i,1) * sqrt(1.D0 - step_pp*i/e)
        enddo
        Pp0(:,1) = Pp0(:,1) * pi / ( serv_k * step_pp * 8.D0 )  !-- normilize to pi / 8
        serv_k = 0.D0
        do i=1, size(Pp1,1)-1
            serv_k = serv_k + Pp1(i,1) * sqrt(1.D0 - step_pp*i/e)
        enddo
        Pp1(:,1) = Pp1(:,1) * pi / ( serv_k * step_pp * 8.D0 )  !-- normilize to pi / 8

        serv_k = sqrt( e**3 / step_pp ) 
        Pp0(1,1) = Pp0(1,1) * serv_k
        Pp0(1,2) = Pp0(1,1) / step_pp
        Pp0(1,3) = 0.D0
        Pp1(1,1) = Pp1(1,1) * serv_k
        Pp1(1,2) = Pp1(1,1) / step_pp
        Pp1(1,3) = 0.D0
        do i=2, size(Pp0,1)            !-- / V
            Pp0(i,1) = Pp0(i,1) * serv_k / sqrt( real(i) )
            Pp0(i,2) = (Pp0(i,1) - Pp0(i-1,1)) / step_pp
            Pp0(i,3) = Pp0(i-1,1)*i - Pp0(i,1)*(i-1)
            Pp1(i,1) = Pp1(i,1) * serv_k / sqrt( real(i) )
            Pp1(i,2) = (Pp1(i,1) - Pp1(i-1,1)) / step_pp
            Pp1(i,3) = Pp1(i-1,1)*i - Pp1(i,1)*(i-1)
        enddo       
    endsubroutine
!----------------------------------------------------------
!-- array of events generation ----------------------------
!----------------------------------------------------------
!-- changed global: Evs                                  --
!-- changed private: t, y1, y2, v3, dE, ex2              --
    subroutine Events_gen(name)
        character(*), intent(in) :: name

        type(event) :: serv_ev   !-- for any needs
        real(kind) :: wwmax, &     !-- max of probability (majorating block)
                     g1, g2, &  !-- width of resonances
                     time, &    !-- time diagnostic
                     serv_k   !-- for any needs
        integer :: i, &   !-- counter
                   serv_i  !-- for any needs
        logical :: serv_l   !-- for any needs
      
        y1 = set_system(core, part1, part2, e)  !-- set particles system "Y1"
        y2 = set_system(part2, core, part1, e)  !-- set particles system "Y2"
        t = set_system(part1, part2, core, e)   !-- set particles system "T"

        if ( do_pp ) then
            open(unit=4,file='log/'//make_name(e,Confs(1).r1.eng,Confs(1).r1.tet)//'pp_diagn.log')   
            serv_k = 0.D0
            i = 1
            write(4,'(4(a13, 1x))') '', 'W0pp(E)', 'W1pp(E)', ''
            do while ( serv_k .LE. e )
                serv_i = int(serv_k / step_pp) + 1
                write(4,'(4(e13.6, 1x))') serv_k, Pp0(i,1), Pp1(i,1), Pp0(serv_i,2) * serv_k + Pp0(serv_i,3)
                serv_k = serv_k + step_pp
                i = i + 1
            enddo
            close(4)
            print "(/,A,/)", 'log/'//make_name(e,Confs(1).r1.eng,Confs(1).r1.tet)//'_pp_diagn.log created'
        endif

        if ( do_ps ) then
            open(unit=6,file='log/psr_diagn.log')
            serv_k = 0.1D0
            !write(6,'(4(a13, 1x))') '', '|A(E)|^2', '|Aps(E)|^2', '|A+Aps(E)|^2'
            write(6,'(4(a13, 1x))') '', 'Re(A(E))', 'Re(Aps(E))', 'Re(A+Aps(E))'
            do while ( serv_k .LE. 5.D0 )
                !write(6,'(4(e13.6, 1x))') serv_k, cdabs(Arp(r1, y1.x, serv_k))**2, cdabs(Aps(r_ps, r1, y1.x, serv_k))**2, cdabs(Arp(r1, y1.x, serv_k)+Aps(r_ps, r1, y1.x, serv_k))**2
       !         write(6,'(4(e13.6, 1x))') serv_k, dreal(Arp(r1, y1.x, serv_k)), dreal(Aps(r_ps, r1, y1.x, serv_k)), dreal(Arp(r1, y1.x, serv_k)+Aps(r_ps, r1, y1.x, serv_k))
                serv_k = serv_k + 0.01D0
            enddo
            close(6)
            print "(/,A,/)", 'log/psr_diagn.log created'

            open(unit=6,file='log/psi_diagn.log')
            serv_k = 0.1D0
            write(6,'(4(a13, 1x))') '', 'Im(A(E))', 'Im(Aps(E))', 'Im(A+Aps(E))'
            do while ( serv_k .LE. 5.D0 )
                !write(6,'(4(e13.6, 1x))') serv_k, cdabs(Arp(r1, y1.x, serv_k))**2, cdabs(Aps(r_ps, r1, y1.x, serv_k))**2, cdabs(Arp(r1, y1.x, serv_k)+Aps(r_ps, r1, y1.x, serv_k))**2
       !         write(6,'(4(e13.6, 1x))') serv_k, dimag(Arp(r1, y1.x, serv_k)), dimag(Aps(r_ps, r1, y1.x, serv_k)), dimag(Arp(r1, y1.x, serv_k)+Aps(r_ps, r1, y1.x, serv_k))
                serv_k = serv_k + 0.01D0
            enddo
            close(6)
            print "(/,A,/)", 'log/psi_diagn.log created'

        !    open(unit=4,file='log/ps2_diagn.log')
       !     serv_k = 0.1D0
        !    write(4,'(4(a13, 1x))') '', '|A(E)|^2', '|Aps(E)|^2', '|A+Aps(E)|^2'
       !     do while ( serv_k .LE. 5.D0 )
        !        write(4,'(4(e13.6, 1x))') serv_k, cdabs(Arp(r2, y1.x, serv_k))**2, cdabs(Aps(r_ps, r2, y1.x, serv_k))**2, cdabs(Arp(r2, y1.x, serv_k)+Aps(r_ps, r2, y1.x, serv_k))**2
       !         serv_k = serv_k + 0.01D0
        !    enddo
        !    close(4)
       !     print "(/,A,/)", 'log/ps2_diagn.log created'
        endif

        g1 = gamm(Confs(1).r1, y1.x, Confs(1).r1.eng)  
        g2 = gamm(Confs(1).r2, y2.x, Confs(1).r2.eng)   
        print    "(/,A11,F10.3)", 'Etotal, MeV', e
        write(5, "(/,A11,F10.3)") 'Etotal, MeV', e
        print    "(4(A10),2(/,F10.3,ES10.2,2(F10.3)))", 'Er, MeV', 'Gr, MeV', 'Thet^2', 'Rch, fm', Confs(1).r1.eng, g1, Confs(1).r1.tet, Confs(1).r1.rch, Confs(1).r2.eng, g2, Confs(1).r2.tet, Confs(1).r2.rch
        write(5, "(4(A10),2(/,F10.3,ES10.2,2(F10.3)))") 'Er, MeV', 'Gr, MeV', 'Thet^2', 'Rch, fm', Confs(1).r1.eng, g1, Confs(1).r1.tet, Confs(1).r1.rch, Confs(1).r2.eng, g2, Confs(1).r2.tet, Confs(1).r2.rch
        dE = e - Confs(1).r1.eng - Confs(1).r2.eng  !-- usefull constant   
        ex2 = e * 2.D0   !-- usefull constant
        v3 = ( 4.D0*dE*dE + g1*g1 + g2*g2 ) * e / 16.D0 / pi       !-- coefficient v3 / 2  ( 2 comes from symmetric A**2 )
        !-- majorating ------------------------------------
        time = dsecnds(0.D0)
        Evs = set_null(serv_ev)
        i = 1
        wwmax = 0.D0
        call random_seed()
        do while (i .LE. size(Evs))
            Evs(i) = event_gen()
            if ( Evs(i).ww .GT. wwmax ) then
                wwmax = Evs(i).ww
                write(5,"(I,ES)") i, wwmax
                print   "(I,ES)", i, wwmax
            endif
            i = i + 1
        enddo
        open(unit=2,file='log/'//make_name(e,Confs(1).r1.eng,Confs(1).r1.tet)//'diagn.data') 
        write(2,"(9(1x,A12))") 'Ww', 'Ww_sym', 'Ww_asym', 'Y1_eps', 'Y1_cos', 'Y2_eps', 'Y2_cos', 'T_eps', 'T_cos'
        do i = 1,size(Evs)
            write(2,"(9(1x,e12.5))"), Evs(i).ww, Evs(i).ww_sym, Evs(i).ww_asym, &
                                      Evs(i).y1.eps, Evs(i).y1.cos_teta, Evs(i).y2.eps, Evs(i).y2.cos_teta, Evs(i).t.eps, Evs(i).t.cos_teta 
        enddo
        close(2)
        print "(/,A)", '-diagn.data created'
        print    "(/,A,ES)", 'max was found: ', wwmax !, 'generation of ', size(Evs), ' events...'
        write(5, "(/,A,ES)") 'max was found: ', wwmax !, 'generation of ', size(Evs), ' events...'
        print "(A,F10.4,A)",  'it takes ', dsecnds(time), ' sekonds'

        serv_k = sum(Evs.ww) * 16.D0 * pi**2 / size(Evs)
        print    "(/,A14,I8,A8,ES,/)", 'Gtotal, MeV (', size(Evs), 'events) ', serv_k

        !-- generation ------------------------------------
        if (do_gen) then
            time = dsecnds(0.D0)
            print    "(/,A,I,A)", 'generation of ', size(Evs), ' events...'
            write(5, "(/,A,I,A)") 'generation of ', size(Evs), ' events...'
            i = 1
            serv_i = 0   !-- counter of all generation (for integrating)
            do while (i .LE. size(Evs))
                if ( mod(i, nping) .EQ. 0 ) then
                    write (5, "(I,A6)") i, 'ping'
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
            print "(A,F10.4,A)",  'it takes ', dsecnds(time), ' sekonds'
        endif  
    endsubroutine
!----------------------------------------------------------
!-- one event generation ----------------------------------
    type(event) function event_gen()
        real(kind), dimension(6) :: rand != &
    !    (/0.249572452257228, 0.994163934376058, 0.640602593985953, 0.590190373438494, 0.280130785801967, 0.253168121687737/)  !-- 5 random parameters and 6th for Monte-Carlo
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
    endfunction
!----------------------------------------------------------
!-- probability calculating -------------------------------
    subroutine ww(ev)
        type(event), intent(inout) :: ev

        complex(kind) :: a11, a22, a12, a21, znam
        real(kind) :: ey1x, ey2x
        integer :: i

        ey1x = e * ev.y1.eps
        ey2x = e * ev.y2.eps
        znam = cmplx(dE, (gamm(Confs(1).r1, y1.x, ey1x) + gamm(Confs(1).r2, y2.x, ey2x))/2.D0)   
        if ( do_ps ) then
            a11 = ( Arp(Confs(1).r1, y1.x, ey1x) + Aps(r_ps, Confs(1).r1, y1.x, ey1x) )* sqrt(gamm(Confs(1).r2, y1.y, e-ey1x)) / znam     !-- 1st resonance of 1st two-particles
            a22 = ( Arp(Confs(1).r2, y2.x, ey2x) + Aps(r_ps, Confs(1).r2, y2.x, ey2x) ) * sqrt(gamm(Confs(1).r1, y2.y, e-ey2x)) / znam     !-- 2st resonance of 2st two-particles
            znam = cmplx(dE, (gamm(Confs(1).r2, y1.x, ey1x) + gamm(Confs(1).r1, y2.x, ey2x))/2.D0)
            a12 = ( Arp(Confs(1).r1, y2.x, ey2x) + Aps(r_ps, Confs(1).r1, y2.x, ey2x) ) * sqrt(gamm(Confs(1).r2, y2.y, e-ey2x)) / znam    !-- 1st resonance of 2st two-particles
            a21 = ( Arp(Confs(1).r2, y1.x, ey1x) + Aps(r_ps, Confs(1).r2, y1.x, ey1x) ) * sqrt(gamm(Confs(1).r1, y1.y, e-ey1x)) / znam    !-- 2st resonance of 1st two-particles
        else
            a11 = Arp(Confs(1).r1, y1.x, ey1x) * sqrt(gamm(Confs(1).r2, y1.y, e-ey1x)) / znam     !-- 1st resonance of 1st two-particles
            a22 = Arp(Confs(1).r2, y2.x, ey2x) * sqrt(gamm(Confs(1).r1, y2.y, e-ey2x)) / znam     !-- 2st resonance of 2st two-particles
            znam = cmplx(dE, (gamm(Confs(1).r2, y1.x, ey1x) + gamm(Confs(1).r1, y2.x, ey2x))/2.D0)
            a12 = Arp(Confs(1).r1, y2.x, ey2x) * sqrt(gamm(Confs(1).r2, y2.y, e-ey2x)) / znam    !-- 1st resonance of 2st two-particles
            a21 = Arp(Confs(1).r2, y1.x, ey1x) * sqrt(gamm(Confs(1).r1, y1.y, e-ey1x)) / znam    !-- 2st resonance of 1st two-particles
        endif
        ev.ww_sym = v3 * half_ww(Cjl0, .TRUE.)
    !    print *, 'ww_sym : ', ev.ww_sym
    !    pause
        ev.ww_asym = v3 * half_ww(Cjl1, .FALSE.) 
    !    print *, 'ww_asym : ', ev.ww_asym
    !    pause
        if (do_pp) then
            ey1x = ev.t.eps *e
            i = int(ey1x / step_pp) + 1
            ev.ww_sym = ( Pp0(i,2) * ey1x + Pp0(i,3) ) * ev.ww_sym
            ev.ww_asym = ( Pp1(i,2) * ey1x + Pp1(i,3) ) * ev.ww_asym
        endif
        ev.ww = ev.ww_sym + ev.ww_asym
   !     print *, 'ww : ', ev.ww
   !     pause
        
        contains

        real(kind) function half_ww(Cjl, is_sym)
            type(recoup), dimension(:), intent(in) :: Cjl
            logical, intent(in) :: is_sym

            integer :: i, ist, ifin
            real(kind) :: ww_l

            half_ww = 0.D0
            do i = 1, size(Cjl)
                ww_l = 0.D0
                ist = Cjl(i).sub.i_st
                ifin = ist
                do while (ist .LE. Cjl(i).sub.i_fin)
                    do while ( ( ifin .LE. Cjl(i).sub.i_fin ) .AND. ( Confs(1).Clbs(ifin).ml .EQ. Confs(1).Clbs(ist).ml ) )
                        ifin = ifin + 1
                    enddo
                    ww_l = ww_l + cdabs( Cjl(i).c * a_sym(is_sym, Confs(1).Clbs(ist:ifin-1)) )**2
                 !   print *, 'recoup = ', Cjl(i).c * v3
                !    print *, 'a_sym = ', a_sym(is_sym, Confs(1).Clbs(ist:ifin-1))
                !    print *, 'sum = ', cdabs( v3 * Cjl(i).c * a_sym(is_sym, Confs(1).Clbs(ist:ifin-1)) )**2
                !    pause
                    ist = ifin
                enddo
                half_ww = half_ww + ww_l / (2.D0*Cjl(i).l + 1)
            enddo  
        end function

        complex(kind) function a_sym(is_sym, Sub_clbs)
            logical, intent(in) :: is_sym
            type(clgd), dimension(:), intent(in) :: Sub_clbs

            a_sym = lxl(ev.y1.kx, ev.y1.ky, Sub_clbs) * a11 + (-1.D0)**(Sub_clbs(1).l1)*lxl(ev.y2.ky, ev.y2.kx, Sub_clbs) * a22
            if (is_sym) then
                a_sym = a_sym + lxl(ort(-ev.y2.kx.cos_teta, ev.y2.kx.phi+pi), ev.y2.ky, Sub_clbs) * a12 + &
                                (-1.D0)**(Sub_clbs(1).l1)*lxl(ev.y1.ky, ort(-ev.y1.kx.cos_teta, ev.y1.kx.phi+pi), Sub_clbs) * a21 
            else
                a_sym = a_sym - lxl(ort(-ev.y2.kx.cos_teta, ev.y2.kx.phi+pi), ev.y2.ky, Sub_clbs) * a12 - &
                                (-1.D0)**(Sub_clbs(1).l1)*lxl(ev.y1.ky, ort(-ev.y1.kx.cos_teta, ev.y1.kx.phi+pi), Sub_clbs) * a21 
            endif
        end function

    end subroutine
!----------------------------------------------------------
end module
