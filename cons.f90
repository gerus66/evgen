!-- usefull constants

module cons
    implicit none
    
!== global ================================================
    integer, parameter :: kind = 8   !-- precision of real
    real(kind), parameter :: prec = 2.D-14   !-- precision for every checks

    real(kind), parameter :: eps_min = 0.02D0   !-- min epsilon for "coulfg4" don't produce trash
    integer, parameter :: nping = 100000  !-- period for progress check (evgen will say "ping!" to you)

    real(kind), parameter :: pi = 3.141592653589793238462643383279502884D0, &  ! Pi
                             dvapi = 6.283185307179586476925286766559006D0, &  ! 2*Pi
                             sqpi = 1.7724538509055160272981674833411452D0, &  ! Sqrt(Pi)
                             sqdvapi = 2.506628274631000502415765284811045D0   ! Sqrt(2*Pi)

    real(kind), parameter :: neutron_mass = 937.56563D0, &
                             proton_mass = 938.27647D0, &
                             alpha_const = 0.007281886D0, &  !-- 1/137,..
                             fm_mev = 197.327D0   
    
    real(kind), parameter :: sq3 = 1.7320508075688772935274463415058724D0, &    ! Sqrt(3)
                             sq2 = 1.4142135623730950488016887242096981D0       ! Sqrt(2)
    
    REAL(KIND), PARAMETER :: PI2 = 9.869604401089358618834490999876151D0      ! Pi**2
    REAL(KIND), PARAMETER :: PI3 = 31.006276680299820175476315067101395D0     ! Pi**3 
    REAL(KIND), PARAMETER :: TRIPI = 9.424777960769379715387930149838509D0    ! 3*Pi
    REAL(KIND), PARAMETER :: CHEPI = 12.566370614359172953850573533118012D0   ! 4*Pi
    REAL(KIND), PARAMETER :: PIDVA = 1.5707963267948966192313216916397514D0   ! Pi/2
    REAL(KIND), PARAMETER :: PICHE = 0.7853981633974483096156608458198757D0   ! Pi/4
    REAL(KIND), PARAMETER :: RT2DPI = 0.79788456080286535587989211986876373D0 ! Sqrt(2/Pi)
    REAL(KIND), PARAMETER :: CHISLOE = 2.7182818284590452353602874713526625D0   ! E
    REAL(KIND), PARAMETER :: CHISLOE2 = 7.389056098930650227230427460575008D0   ! E**2
    REAL(KIND), PARAMETER :: CHISLOE3 = 20.085536923187667740928529654581718D0  ! E**3
    REAL(KIND), PARAMETER :: SQM2 = 0.7071067811865475244008443621048490D0   ! 1/Sqrt(2)
    REAL(KIND), PARAMETER :: SQM3 = 0.5773502691896257645091487805019575D0   ! 1/Sqrt(3)
    REAL(KIND), PARAMETER :: SQ5 = 2.2360679774997896964091736687312762D0    ! Sqrt(5)
    REAL(KIND), PARAMETER :: SQM5 = 0.4472135954999579392818347337462553D0   ! 1/Sqrt(5)
 
    REAL(KIND), PARAMETER :: AH2 = 1.0D0 / 20.748D0     !  AMH  = 1.D0 / 20.7575925 D0
   
end module cons
