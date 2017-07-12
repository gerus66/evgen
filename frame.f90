!=== units (1) - input: config.txt, grid.txt, ps.txt, pp.txt etc
!===       (2) - output diagnostic: log/{%%}diagn.data
!===       (3) - output: log/{%%}.data
!===       (4) - diagnostic: log/pp_diagn.log, log/ps1_diagn.log, etc
!===       (5) - log: log/gen.log
!===       (6) - parallel input: {%path_to_element_file%}.txt
program evgen
    use mydata, only: kind, pi
    use types
    use output
    use generator
    implicit none

    type(clgd), dimension(:), allocatable :: Clbs_all   !-- all Clebsh-Gordan coefficients
    type(recoup), dimension(:), allocatable :: Cjl_all  !-- all recouple coefficients
    
    real(kind) :: e_nuc, &  !-- input: energy of beam MeV/nuc
                  et_st, er_st, tet_st, et_step, er_step, tet_step, et_fin, er_fin, tet_fin, &  !-- for generation of grid of [Events]
                  serv_k   !-- for any needs
    integer :: i, ii, i_et, i_er, i_tet, &    !-- counters
               serv_i, &    !-- for any needs
               nbin   !-- number of bins for diagnostic
    character(50) :: serv_c   !-- for any needs
    logical :: serv_l, serv_l2    !-- for any needs
    type(event) :: serv_ev

!============ open log-file ===============================
    open(unit=5,file='log/gen.log')

!============ read config-file ============================
    open(unit=1,file='config.txt')
    read(1,*) serv_c; read(1,*) serv_i   !-- number of events
    allocate( Evs(serv_i) )
    write(5,"(/,A7,I/)") 'events ', serv_i
    print   "(/,A7,I/)", 'events ', serv_i
    read(1,*) nping   !-- period for progress check 
    read(1,*) nbin; read(1,*) serv_c  !-- number of bins
    do_grid = read_num( 1, (/0, 1/) )  !-- one case / grid
    do_gen = read_num( 1, (/0, 1/) )    !-- diagnostics / diagnostics&generation
    do_pp = read_num( 1, (/0, 1/) )    !-- no App / include App (proton-proton)
    do_ps = read_num( 1, (/0, 1/) )    !-- no Aps / include Aps (potential scattering)
    write(5,"(4(A15,I3,/))") 'grid of cases', do_grid, 'generation', do_gen, 'proton-proton', do_pp, 'potential scat.', do_ps
    print   "(4(A15,I3,/))", 'grid of cases', do_grid, 'generation', do_gen, 'proton-proton', do_pp, 'potential scat.', do_ps
    read(1,*) serv_c; read(1,*) serv_i   !-- number of configurations
    allocate( Confs(serv_i) )
    do i = 1, size(Confs)
        read(1,*) serv_c
        read(1,*) serv_k
        open(unit=6,file='elements/'//serv_c) 
        if ( i .EQ. 1 ) then
            read(6,*) serv_c; read(6,*) core.z
            read(6,*) part1.z
            read(6,*) part2.z 
            read(6,*) core.a; core = setmass(core)
            read(6,*) part1.a; part1 = setmass(part1)
            read(6,*) part2.a; part2 = setmass(part2)
            read(6,*) serv_c; read(6,*) j  !-- J total
            read(6,*) e                     !-- E total
        else
            do ii = 1, 10
                read(6,*) serv_c
            enddo
        endif
        write(5,"(A6,/,A,/)") '------', serv_c
        print   "(A6,/,A,/)", '------', serv_c
        read(6,*) serv_c; read(6,*) Confs(i).r1.l   !-- 1st resonance
        read(6,*) Confs(i).r1.jx2
        read(6,*) Confs(i).r1.rch
        read(6,*) Confs(i).r1.eng
        read(6,*) Confs(i).r1.tet
        read(6,*) serv_c; read(6,*) Confs(i).r2.l   !-- 2nd resonance
        read(6,*) Confs(i).r2.jx2
        read(6,*) Confs(i).r2.rch
        read(6,*) Confs(i).r2.eng
        read(6,*) Confs(i).r2.tet
        close(6)  
        write(5,"(5(A5),/,A1,I2,A2,I5,I3,A2,I5,A2,I3,/)") 'j1', 'l1', 'j2', 'l2', 'J', '{', Confs(i).r1.jx2, '/2', Confs(i).r1.l, Confs(i).r2.jx2, '/2', Confs(i).r2.l, '}', j  !-- all quantum information
        print   "(5(A5),/,A1,I2,A2,I5,I3,A2,I5,A2,I3,/)", 'j1', 'l1', 'j2', 'l2', 'J', '{', Confs(i).r1.jx2, '/2', Confs(i).r1.l, Confs(i).r2.jx2, '/2', Confs(i).r2.l, '}', j
!--- Clebsh-Gordan coeficients ----------------------------
        write(5,"(A,/,4(A5),A8)") 'Clebsh-Gordan coefficients', 'L', 'Ml', 'ml1', 'ml2', 'coef'
        print   "(A,/,4(A5),A8)", 'Clebsh-Gordan coefficients', 'L', 'Ml', 'ml1', 'ml2', 'coef'
        allocate( Clbs_all( (6*j+4)*(2*Confs(i).r1.l+1) ) )    !-- size is a bit more then needed, but sufficient for sure
        call make_all(Confs(i).r1.l, Confs(i).r2.l, j, Clbs_all, serv_i)    !-- Clebsh-Gordan coefs for the [l1xl2]->J configuration: all; number of null ones 
        allocate ( Confs(i).Clbs(size(Clbs_all) - serv_i) )
        Confs(i).Clbs = Clbs_all(1:size(Clbs_all) - serv_i)   !-- Clebsh-Gordan coefs for the [l1xl2]->J configuration: notnull only
        deallocate (Clbs_all)
        do ii = 1, size(Confs(i).Clbs)
            print    "(4(I5),f10.5)", Confs(i).Clbs(ii).l, Confs(i).Clbs(ii).ml, Confs(i).Clbs(ii).m1, Confs(i).Clbs(ii).m2, Confs(i).Clbs(ii).c  !-- show all Cl-Gd coefs
            write(5, "(4(I5),f10.5)") Confs(i).Clbs(ii).l, Confs(i).Clbs(ii).ml, Confs(i).Clbs(ii).m1, Confs(i).Clbs(ii).m2, Confs(i).Clbs(ii).c  !-- and write it to log
        enddo

    enddo
    read(1,*) serv_c; read(1,*) prec  !-- precision of check that everything going allright
    read(1,*) eps_min    !-- min epsilon for coulfg don't produce trash
    read(1,*) e_nuc      !-- energy of beam
    close(1)

!============ read other files ============================
    if ( do_grid ) then
        open(unit=1,file='grid.txt')
        read(1,*) et_step 
        read(1,*) et_fin
        read(1,*) er_step 
        read(1,*) er_fin
        read(1,*) tet_step
        read(1,*) tet_fin
        close(1)
    endif
    if ( do_ps ) then
        open(unit=1,file='ps.txt')
        read(1,*) r_ps
        close(1)
    endif
    if ( do_pp ) then
        open(unit=1,file='pp.txt')
        read(1,*) pp_0
        read(1,*) pp_1
        close(1)
    endif

!==========================================================
!============ some preparations ===========================
    serv_k = (core.a + part1.a + part2.a) * e_nuc / (core.m + part1.m + part2.m)
    v_beam = sqrt( 1.D0 + 1.D0 / (serv_k*serv_k + 2.D0*serv_k) ) / 2.D0  !-- for [event].tet calculation; really it's 1/2*v_beam; GLOBAL from generator.f90

    print    "(/,A,/,2(A5),A8,A15)", 'recouple coefficients', 'L', 'S', 'coef', 'Cl-Gd coefs'
    write(5, "(/,A,/,2(A5),A8,A15)") 'recouple coefficients', 'L', 'S', 'coef', 'Cl-Gd coefs'
    allocate(Cjl_all(1))
    call make_all(Confs(1).r1.jx2, Confs(1).r1.l, Confs(1).r2.jx2, Confs(1).r2.l, j, 0, Confs(1).Clbs, Cjl_all, serv_i)  !-- recouple coefs for the [j1,l1,j2,l2]->S=0,J configuration: all; number of null ones
    allocate (Cjl0(size(Cjl_all)-serv_i))
    Cjl0 = Cjl_all(1:size(Cjl_all)-serv_i)  !-- recouple coefs for the [j1,l1,j2,l2]->S=0,J configuration: notnull only
    deallocate(Cjl_all)
    do i = 1, size(Cjl0) 
        print    "(2(I5),f10.5,I5,A4,I2)", Cjl0(i).l, 0, Cjl0(i).c, Cjl0(i).sub.i_st, ' to ', Cjl0(i).sub.i_fin  !-- show all recouple coefs for S=0
        write(5, "(2(I5),f10.5,I5,A4,I2)") Cjl0(i).l, 0, Cjl0(i).c, Cjl0(i).sub.i_st, ' to ', Cjl0(i).sub.i_fin  !-- and write it to log
    enddo
    allocate(Cjl_all(3))
    call make_all(Confs(1).r1.jx2, Confs(1).r1.l, Confs(1).r2.jx2, Confs(1).r2.l, j, 1, Confs(1).Clbs, Cjl_all, serv_i)  !-- recouple coefs for the [j1,l1,j2,l2]->S=1,J configuration: all; number of null ones
    allocate (Cjl1(size(Cjl_all)-serv_i))
    Cjl1 = Cjl_all(1:size(Cjl_all)-serv_i)  !-- recouple coefs for the [j1,l1,j2,l2]->S=1,J configuration: notnull only
    deallocate(Cjl_all)
    do i = 1, size(Cjl1)
        print    "(2(I5),f10.5,I5,A4,I2)", Cjl1(i).l, 1, Cjl1(i).c, Cjl1(i).sub.i_st, ' to ', Cjl1(i).sub.i_fin  !-- show all recouple coefs for S=1
        write(5, "(2(I5),f10.5,I5,A4,I2)") Cjl1(i).l, 1, Cjl1(i).c, Cjl1(i).sub.i_st, ' to ', Cjl1(i).sub.i_fin  !-- and write it to log
    enddo

    print *, 'local end'
!==========================================================   
!============ general block ===============================
    select case ( do_grid )   !-- one case / grid of cases
        case (.False.)     !-- one case
            if ( do_pp ) call pp_read
            call Events_gen(make_name(e,Confs(1).r1.eng,Confs(1).r1.tet))
            if (do_gen) then
                call Events_output()
            else
                call diagn_output()
            endif
            if ( do_pp ) deallocate (Pp0, Pp1)
        case (.True.)     !-- grid of cases
            print *, 'sorry, we still cant calculate a grid of cases'
 !           et_st = e; er_st = r1.eng; tet_st = r1.tet
 !           print "(a17,a10,a10,/,3(a7,3(f10.4),/))", 'start', 'final', 'step', 'E_total', et_st, et_fin, et_step, 'E_res', er_st, er_fin, er_step, 'theta^2', tet_st, tet_fin, tet_step
 !           do while (e .LE. et_fin)
 !               if ( do_pp ) call pp_read
 !               r1.eng = er_st; r2.eng = er_st
 !               do while (r1.eng .LE. er_fin)
 !                   r1.tet = tet_st; r2.tet = tet_st
 !                   do while (r1.tet .LE. tet_fin) 
 !                       call Events_gen(make_name(e,r1.eng,r1.tet))
 !                       if (do_gen) then
 !                           call Events_output()
 !                       else
 !                           call diagn_output()
 !                       endif
 !                       r1.tet = r1.tet + tet_step; r2.tet = r2.tet + tet_step
 !                   enddo
 !                   r1.eng = r1.eng + er_step; r2.eng = r2.eng + er_step
 !               enddo
 !               e = e + et_step
 !               if ( do_pp ) deallocate (Pp0, Pp1)
 !           enddo
    endselect

!==========================================================
    deallocate (Cjl0, Cjl1)
    deallocate (Evs)
    do i = 1, size(Confs)
        deallocate (Confs(i).Clbs)
    enddo
    deallocate (Confs)
    print *, 'the end'
    close(5)
!==========================================================
    contains

    subroutine Events_output()
        write (5, "(A20,2(A10),A20)") 'dx/l', 'dy/l', 'dz/l', 'average, %'
        serv_k = 100 / sum(abs(Evs.p1))
        write (5, "(A10,3(f10.2))") 'p1 ', sum(Evs.p1.x)*serv_k, sum(Evs.p1.y)*serv_k, sum(Evs.p1.z)*serv_k
        serv_k = 100 / sum(abs(Evs.p2))
        write (5, "(A10,3(f10.2))") 'p2 ', sum(Evs.p2.x)*serv_k, sum(Evs.p2.y)*serv_k, sum(Evs.p2.z)*serv_k
        serv_k = 100 / sum(abs(Evs.core))
        write (5, "(A10,3(f10.2))") 'core ', sum(Evs.core.x)*serv_k, sum(Evs.core.y)*serv_k, sum(Evs.core.z)*serv_k
        write( serv_c, "(I1,A1,I2,A1,I1,A1,I1,A1,I1,A1,I1)" )  floor(e),'-',int(mod(e*100,100.D0)),'_',floor(Confs(1).r1.eng),'-',int(mod(Confs(1).r1.eng*10,10.D0)), '_',floor(Confs(1).r1.tet),'-',int(mod(Confs(1).r1.tet*10,10.D0))
        call matrix_to_file(histogram(Evs.t.cos_teta, nbin, .TRUE., -1.D0, 2.D0), 'log/'//serv_c, '_cos_t.log', 'e13.6, 1x')
        call matrix_to_file(histogram(Evs.t.eps, nbin, .TRUE., 0.D0, 1.D0), 'log/'//serv_c, '_eps_t.log', 'e13.6, 1x')
        call short_3d_diagnostic(Evs.t.cos_teta, Evs.t.eps, nbin, serv_c, '_3d_t.log')
        call matrix_to_file(histogram(Evs.y1.cos_teta, nbin, .TRUE., -1.D0, 2.D0), 'log/'//serv_c, '_cos_y1.log', 'e13.6, 1x')
        call matrix_to_file(histogram(Evs.y1.eps, nbin, .TRUE., 0.D0, 1.D0), 'log/'//serv_c, '_eps_y1.log', 'e13.6, 1x')
        call short_3d_diagnostic(Evs.y1.cos_teta, Evs.y1.eps, nbin, serv_c, '_3d_y1.log')
        call matrix_to_file(histogram(Evs.y2.cos_teta, nbin, .TRUE., -1.D0, 2.D0), 'log/'//serv_c, '_cos_y2.log', 'e13.6, 1x')
        call matrix_to_file(histogram(Evs.y2.eps, nbin, .TRUE., 0.D0, 1.D0), 'log/'//serv_c, '_eps_y2.log', 'e13.6, 1x')
        call short_3d_diagnostic(Evs.y2.cos_teta, Evs.y2.eps, nbin, serv_c, '_3d_y2.log')
        call matrix_to_file(histogram(Evs.tet, nbin, .TRUE., 0.D0, 0.1D0), 'log/'//serv_c, 'tett.log', 'e13.6, 1x')
        open(unit=3,file='log/'//serv_c//'.data')
        do i = 1,size(Evs)
            write(3,"(10(1x,e12.5))"), Evs(i).eng, Evs(i).p1.x, Evs(i).p1.y, Evs(i).p1.z, Evs(i).p2.x, Evs(i).p2.y, Evs(i).p2.z, Evs(i).core.x, Evs(i).core.y, Evs(i).core.z 
        enddo
        close(3)
        print *, '====================================='
        write (5, "(A)") '====================================='
    endsubroutine


    subroutine diagn_output()
        serv_k = size(Evs) / sum(Evs.ww)
        serv_c = make_name(e, Confs(1).r1.eng, Confs(1).r1.tet)
        call matrix_to_file(ww_diagnostic(Evs.y1.eps, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, serv_k, 0.D0, 1.D0), 'log/'//serv_c(1:12), 'diagn_eps_y1.log', 'e13.6, 1x')
        call matrix_to_file(ww_diagnostic(Evs.y1.cos_teta, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, serv_k, -1.D0, 2.D0), 'log/'//serv_c(1:12), 'diagn_cos_y1.log', 'e13.6, 1x')
        call matrix_to_file(ww_diagnostic(Evs.y2.eps, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, serv_k, 0.D0, 1.D0), 'log/'//serv_c(1:12), 'diagn_eps_y2.log', 'e13.6, 1x')
        call matrix_to_file(ww_diagnostic(Evs.y2.cos_teta, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, serv_k, -1.D0, 2.D0), 'log/'//serv_c(1:12), 'diagn_cos_y2.log', 'e13.6, 1x')
        call matrix_to_file(ww_diagnostic(Evs.t.eps, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, serv_k, 0.D0, 1.D0), 'log/'//serv_c(1:12), 'diagn_eps_t.log', 'e13.6, 1x')
        call matrix_to_file(ww_diagnostic(Evs.t.cos_teta, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, serv_k, -1.D0, 2.D0), 'log/'//serv_c(1:12), 'diagn_cos_t.log', 'e13.6, 1x')
        call matrix_to_file(ww_3d_diagnostic(Evs.y1.eps, Evs.y1.cos_teta, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, 1.D0), 'log/'//serv_c(1:12), 'diagn_3d_y1.log', 'e13.6, 1x')
        call matrix_to_file(ww_3d_diagnostic(Evs.y2.eps, Evs.y2.cos_teta, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, 1.D0), 'log/'//serv_c(1:12), 'diagn_3d_y2.log', 'e13.6, 1x')
        call matrix_to_file(ww_3d_diagnostic(Evs.t.eps, Evs.t.cos_teta, Evs.ww, Evs.ww_sym, Evs.ww_asym, nbin, 1.D0), 'log/'//serv_c(1:12), 'diagn_3d_t.log', 'e13.6, 1x')
    endsubroutine
endprogram
