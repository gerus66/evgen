!=== units (1) - input: config.txt, grid.txt, ps.txt, pp.txt etc
!===       (2) - output diagnostic: log/{%%}diagn.data
!===       (3) - output: log/{%%}.data
!===       (4) - diagnostic: log/pp_diagn.log, log/ps1_diagn.log, etc
!===       (5) - log: log/gen.log
!===       (6) - parallel input: {%path_to_element_file%}.txt
program evgen
    use cons, only: kind, pi
    use types
    use output
    use output_ev
    use generator
    use pp
    use ps
    use fit
    implicit none
  
    real(kind) :: et_st, er_st, tet_st, et_step, er_step, tet_step, et_fin, er_fin, tet_fin, &  !-- for generation 
                                                                                                !-- of grid of [Events]
                  serv_k, serv_kk   !-- for any needs
    integer :: i, ii, iii, i_et, i_er, i_tet, &    !-- counters
               serv_i    !-- for any needs
    character(50) :: serv_c   !-- for any needs
    logical :: serv_l, serv_l2    !-- for any needs
    type(event) :: serv_ev

!============ open log-file ===============================
    open(unit=5,file='log/gen.log')

!============ read config files ===========================
    open(unit=1,file='config.txt')
    read(1,*) serv_c; read(1,*) serv_i   !-- number of events
    allocate( Evs(serv_i) )
    write(5,"(/,A7,I/)") 'events ', serv_i
    print   "(/,A7,I/)", 'events ', serv_i
    read(1,*) serv_i; call set_nbin(serv_i); read(1,*) serv_c  !-- number of bins
    read(1,*) v_beam; read(1,*) serv_c      !-- energy of beam
    do_out = read_num( 1, (/0, 1/) )  !-- no output / usual output
    do_grid = read_num( 1, (/0, 1/) )  !-- one case / grid
    do_gen = read_num( 1, (/0, 1/) )    !-- diagnostics / diagnostics&generation
    do_pp = read_num( 1, (/0, 1/) )    !-- no App / include App (proton-proton)
    do_ps = read_num( 1, (/0, 1/) )    !-- no Aps / include Aps (potential scattering)
    do_fit = read_num( 1, (/0, 1/) )    !-- usual work / fitting
    write(5,"(5(A15,I3,/))") 'grid', do_grid, 'gen', do_gen, 'p-p', do_pp, 'ps', do_ps, 'fit', do_fit
    print   "(5(A15,I3,/))", 'grid', do_grid, 'gen', do_gen, 'p-p', do_pp, 'ps', do_ps, 'fit', do_fit
    read(1,*) serv_c; read(1,*) serv_i   !-- number of configurations
    allocate( Confs(serv_i) )
    do i = 1, size(Confs)
        read(1,*) serv_c        !-- filename in directory elements/
        read(1,*) serv_k      !-- [Re]
        read(1,*) serv_kk      !-- [Im]
        Confs(i).coef = cmplx(serv_k, serv_kk)  !-- coefficient for configuration
        open(unit=6,file='elements/'//serv_c) 
        if ( i .EQ. 1 ) then     !-- read element's general data only for 1st case 
            read(6,*) serv_c; read(6,*) core.z   !-- Z
            read(6,*) part1.z
            read(6,*) part2.z 
            read(6,*) core.a; core = setmass(core)  !-- A, M
            read(6,*) part1.a; part1 = setmass(part1)
            read(6,*) part2.a; part2 = setmass(part2)
            read(6,*) serv_c; read(6,*) j  !-- J total
            read(6,*) e                     !-- E total
        else
            do ii = 1, 10
                read(6,*) serv_c
            enddo
        endif
        write(5,"(/,20('-'),/,'#', I2 ,' | coef = (', F4.2 ,', ', F4.2 ,')',/,20('-'))") i, Confs(i).coef
        print   "(/,20('-'),/,'#', I2 ,' | coef = (', F4.2 ,', ', F4.2 ,')',/,20('-'))", i, Confs(i).coef
        read(6,*) serv_c; read(6,*) Confs(i).r1.l   !-- 1st resonance  l
        read(6,*) Confs(i).r1.jx2      !-- j * 2
        read(6,*) Confs(i).r1.rch      !-- radius of chanel
        read(6,*) Confs(i).r1.eng      !-- E res
        read(6,*) Confs(i).r1.tet      !-- theta^2
        read(6,*) serv_c; read(6,*) Confs(i).r2.l   !-- 2nd resonance
        read(6,*) Confs(i).r2.jx2
        read(6,*) Confs(i).r2.rch
        read(6,*) Confs(i).r2.eng
        read(6,*) Confs(i).r2.tet
        read(6,*) serv_c; read(6,*) Confs(i).g3
        close(6) 
        write(5,"(5(A5),/,'{',I2,'/2',I5,I3,'/2',I5,'}',I3,/)") 'j1', 'l1', 'j2', 'l2', 'J', &
             Confs(i).r1.jx2, Confs(i).r1.l, Confs(i).r2.jx2, Confs(i).r2.l, j 
        print   "(5(A5),/,'{',I2,'/2',I5,I3,'/2',I5,'}',I3,/)", 'j1', 'l1', 'j2', 'l2', 'J', &
             Confs(i).r1.jx2, Confs(i).r1.l, Confs(i).r2.jx2, Confs(i).r2.l, j 
!--- recouple coeficients ---------------------------------
        call make_all(Confs(i).r1.jx2,Confs(i).r1.l,Confs(i).r2.jx2,Confs(i).r2.l,j,0,Confs(i).C_ls(:,0)) !-- S=0  
        call make_all(Confs(i).r1.jx2,Confs(i).r1.l,Confs(i).r2.jx2,Confs(i).r2.l,j,1,Confs(i).C_ls(:,1)) !-- S=1
        call recoup_print(Confs(i).C_ls, j)
!--- Clebsh-Gordan coeficients ----------------------------     
        allocate( Confs(i).C_lml(-1:1, -j-1:j+1) )    
        do ii = -1, 1    !-- calculate Cl-Gl for only notnull by recouple-coefs L
            if ( Confs(i).C_ls(ii,0).not_null ) then
                call make_all(Confs(i).r1.l, Confs(i).r2.l, j, Confs(i).C_ls(ii,0).l, Confs(i).C_lml(ii,:))
            else
                if ( Confs(i).C_ls(ii,1).not_null ) then
                    call make_all(Confs(i).r1.l, Confs(i).r2.l, j, Confs(i).C_ls(ii,1).l, Confs(i).C_lml(ii,:))
                else
                    Confs(i).C_lml(ii,:) = set_null(Confs(i).C_lml(-1,0))
                endif
            endif 
        enddo
        call clgd_print(Confs(i).C_lml, j)
!----------------------------------------------------------        
    enddo
    
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
    if ( do_ps ) call ps_read()
!==========================================================
!============ some preparations ===========================
    serv_k = (core.a + part1.a + part2.a) * v_beam / (core.m + part1.m + part2.m)
    v_beam = sqrt( 1.D0 + 1.D0 / (serv_k*serv_k + 2.D0*serv_k) ) / 2.D0  !-- for [event].tet calculation; really it's 1/2*v_beam
!==========================================================   
!============ general block ===============================
    select case ( do_grid )   !-- one case / grid of cases
    case (.False.)     !-- one case
        if ( do_pp ) call pp_interact(e)
        if ( .NOT.do_fit ) then
            call Events_gen()
        else
            call fit_read()
            select case ( size(Confs) )
            case ( 1 )
                Confs(1).coef = cmplx( 1.D0, 0.D0 )
                call Events_gen()
            case ( 2 )
                serv_l = .FALSE.
                do while ( .NOT.serv_l )
                    call get_ab(Confs(1).g3, Confs(2).g3, Confs(1).coef, Confs(2).coef, serv_l)                    
                    write(43,"(2('(',F8.3,' ',F8.3,')'))") Confs(1).coef, Confs(2).coef
                    print    "(2('(',F8.3,' ',F8.3,')'))", Confs(1).coef, Confs(2).coef
                    call Events_gen()
                    call set_is_ww(.TRUE.)
                    call xi_test(hist_3d(Evs.t), hist_3d(Evs.y2))
                  !  if ( get_s().EQ.'t' ) call frame_kolm_test(hist_3d(Evs.t))
                 !   if ( get_s().EQ.'y' ) call frame_kolm_test(hist_3d(Evs.y2))
                enddo
            case ( 3 )
                serv_l = .FALSE.
                i = 1
                do while ( .NOT.serv_l .AND. i.LE.100 )
                    call get_abc(Confs(1).g3, Confs(2).g3, Confs(3).g3, Confs(1).coef, Confs(2).coef, Confs(3).coef, serv_l)
                    write(43,"(3('(',F8.3,' ',F8.3,')'))") Confs(1).coef, Confs(2).coef, Confs(3).coef
                    print    "(3('(',F8.3,' ',F8.3,')'))", Confs(1).coef, Confs(2).coef, Confs(3).coef
                    print *, i
                    call Events_gen()
                    call xi_test(hist_3d(Evs.t), hist_3d(Evs.y2))
                    i = i + 1
               !     if ( get_s().EQ.'t' ) call frame_xi_test(hist_3d(Evs.t))
               !     if ( get_s().EQ.'y' ) call frame_xi_test(hist_3d(Evs.y2))
                enddo
            case default
                print *, 'sorry, we cant fit it'
            endselect
            call fit_dealloc()
        endif
        if ( do_pp ) call pp_dealloc()
    case (.True.)     !-- grid of cases
        if  ( size(Confs) .GT. 1 ) then
            print *, 'sorry, cant calculate a grid of cases with multiple configurations'
        else
            et_st = e; er_st = Confs(1).r1.eng; tet_st = Confs(1).r1.tet
            write (5, "(/,20('-'),/,A20,/,20('-'),/,A17,A10,A10,/,3(A7,3(F10.4),/))") 'grid of parameters:', &
                       'start', 'final', 'step', 'E_total', et_st, et_fin, et_step, &
                       'E_res', er_st, er_fin, er_step, 'theta^2', tet_st, tet_fin, tet_step
            print     "(/,20('-'),/,A20,/,20('-'),/,A17,A10,A10,/,3(A7,3(F10.4),/))", 'grid of parameters:',&
                      'start', 'final', 'step', 'E_total', et_st, et_fin, et_step, &
                      'E_res', er_st, er_fin, er_step, 'theta^2', tet_st, tet_fin, tet_step
            do while (e .LE. et_fin)
                if ( do_pp ) call pp_interact(e)    !-- load pp-data up to this energy
                Confs(1).r1.eng = er_st; Confs(1).r2.eng = er_st
                do while (Confs(1).r1.eng .LE. er_fin)
                    Confs(1).r1.tet = tet_st; Confs(1).r2.tet = tet_st
                    do while (Confs(1).r1.tet .LE. tet_fin) 
                        call Events_gen()
                        Confs(1).r1.tet = Confs(1).r1.tet + tet_step; Confs(1).r2.tet = Confs(1).r2.tet + tet_step
                    enddo
                    Confs(1).r1.eng = Confs(1).r1.eng + er_step; Confs(1).r2.eng = Confs(1).r2.eng + er_step
                enddo
                e = e + et_step
                if ( do_pp ) call pp_dealloc()
            enddo
        endif
    endselect

!==========================================================
    deallocate (Evs)
    do i = 1, size(Confs)
        do ii = -1, 1
            do iii = -j-1, j+1
                deallocate (Confs(i).C_lml(ii,iii).ar)
            enddo
        enddo
        deallocate(Confs(i).C_lml)
    enddo
    deallocate (Confs)
    close(5)

end program
