!-- all global variables for "evgen" project

module vars
    use cons, only: kind
    use types
    implicit none
    
!== global ================================================
    type(event), dimension(:), allocatable :: Evs   !-- all events
    type(conf), dimension(:), allocatable :: Confs  !-- all configurations
    type(part) :: core, part1, part2  !-- physical particles
    real(kind) :: e, &   !-- E total
                  v_beam   !-- 1/2*Vbeam
    integer :: j    !-- J total
    logical :: do_out, do_grid, do_gen, do_pp, do_ps, do_fit   
    
end module
