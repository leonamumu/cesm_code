module Hydrology1Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  Hydrology1Mod
!
! !DESCRIPTION:
! Calculation of
! (1) water storage of intercepted precipitation
! (2) direct throughfall and canopy drainage of precipitation
! (3) the fraction of foliage covered by water and the fraction
!     of foliage that is dry and transpiring.
! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
!
! !PUBLIC TYPES:
   use shr_kind_mod,       only: r8=>shr_kind_r8
   use clm_varctl,         only: iulog

   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: Hydrology1
   public :: CanopyFallFrac_EASS

! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hydrology1
!
! !INTERFACE:
   subroutine Hydrology1(lbc, ubc, lbp, ubp, num_nolakec, filter_nolakec, &
                         num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Calculation of
! (1) water storage of intercepted precipitation
! (2) direct throughfall and canopy drainage of precipitation
! (3) the fraction of foliage covered by water and the fraction
!     of foliage that is dry and transpiring.
! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
! Note:  The evaporation loss is taken off after the calculation of leaf
! temperature in the subroutine clm\_leaftem.f90, not in this subroutine.
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd   , only : clm_a2l
    use clm_varcon   , only : tfrz, istice, istwet, istsoil, istice_mec, isturb, &
                              icol_roof, icol_sunwall, icol_shadewall, denh2o
    use clm_varcon   , only : istcrop
    use FracWetMod   , only : FracWet
    use clm_time_manager , only : get_step_size
    use subgridAveMod, only : p2c
    use SNICARMod    , only : snw_rds_min

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: lbc, ubc                     ! column bounds
    integer, intent(in) :: num_nolakec                  ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)    ! column filter for non-lake points
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/15/02, Peter Thornton: Migrated to new data structures. Required
! adding a PFT loop.
! 4/26/05, Peter Thornton: Made the canopy interception factor fpi max=0.25
!   the default behavior
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arrays
!
    integer , pointer :: cgridcell(:)      ! columns's gridcell
    integer , pointer :: clandunit(:)      ! columns's landunit
    integer , pointer :: pgridcell(:)      ! pft's gridcell
    integer , pointer :: plandunit(:)      ! pft's landunit
    integer , pointer :: pcolumn(:)        ! pft's column
    integer , pointer :: npfts(:)          ! number of pfts in column
    integer , pointer :: pfti(:)           ! column's beginning pft index
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: ctype(:)          ! column type
    real(r8), pointer :: forc_rain(:)      ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)      ! snow rate [mm/s]
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    logical , pointer :: do_capsnow(:)     ! true => do snow capping
    real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)
    real(r8), pointer :: dewmx(:)          ! Maximum allowed dew [mm]
    integer , pointer :: frac_veg_nosno(:) ! fraction of veg not covered by snow (0/1 now) [-]
    real(r8), pointer :: elai(:)           ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)           ! one-sided stem area index with burying by snow
    real(r8), pointer :: h2ocan_loss(:)    ! canopy water mass balance term (column)
    real(r8), pointer :: irrig_rate(:)     ! current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]
    !--------------------------
    !added by Jing Chen, 17 May 2012
    real(r8), pointer :: clumping_EASS(:)  ! clumping index in EASS[-]
    real(r8), pointer :: elai_over_EASS(:) ! elai of overstory[-]
    real(r8), pointer :: elai_under_EASS(:)! elai of understroy[-]
    real(r8), pointer :: t_grnd_EASS(:)    ! ground temperature (Kelvin)
!
! local pointers to original implicit inout arrays
!
    integer , pointer :: snl(:)            ! number of snow layers
    real(r8), pointer :: snowdp(:)         ! snow height (m)
    real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
    real(r8), pointer :: h2ocan(:)         ! total canopy water (mm H2O)
    real(r8), pointer :: qflx_irrig(:)          ! irrigation amount (mm/s)
    integer , pointer :: n_irrig_steps_left(:)  ! number of time steps for which we still need to irrigate today
    !---------------------
    ! added by Jing Chen Oct 20 2012
    integer , pointer :: snl_EASS(:)       ! number of snow layers
    real(r8), pointer :: snowdp_EASS(:)    ! snow height (m)
    real(r8), pointer :: h2osno_EASS(:)    ! snow water (mm H2O)
    real(r8), pointer :: h2ocan_EASS(:)    ! total canopy water (mm H2O)
    !-------------------
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: qflx_prec_intr(:)     ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd(:)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snwcp_liq(:)     ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:)     ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snow_grnd_pft(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd_col(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_rain_grnd(:)     ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: fwet(:)               ! fraction of canopy that is wet (0 to 1)
    real(r8), pointer :: fdry(:)               ! fraction of foliage that is green and dry [-] (new)
    real(r8), pointer :: zi(:,:)               ! interface level below a "z" level (m)
    real(r8), pointer :: dz(:,:)               ! layer depth (m)
    real(r8), pointer :: z(:,:)                ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water (kg/m2)
    real(r8), pointer :: frac_iceold(:,:)      ! fraction of ice relative to the tot water
    real(r8), pointer :: snw_rds(:,:)          ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: mss_bcpho(:,:)        ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcphi(:,:)        ! mass of hydrophilic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bctot(:,:)        ! total mass of BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bc_col(:)         ! total column mass of BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bc_top(:)         ! total top-layer mass of BC (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)        ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)        ! mass of hydrophilic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_octot(:,:)        ! total mass of OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_oc_col(:)         ! total column mass of OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_oc_top(:)         ! total top-layer mass of OC (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)         ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)         ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)         ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)         ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dsttot(:,:)       ! total mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_col(:)        ! total column mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_top(:)        ! total top-layer mass of dust in snow (col,lyr) [kg]

    !--------------------------
    !added by Jing Chen, 17 May 2012
    real(r8), pointer :: densnow_EASS(:)          ! snow denisity on the ground [kg m-3]
    real(r8), pointer :: frac_rain_over_EASS(:)   ! fraction of rainfall on the overstory(0~1)
    real(r8), pointer :: frac_rain_under_EASS(:)  ! fraction of rainfall on the understory(0~1)
    real(r8), pointer :: frac_snow_grnd_EASS(:)   ! fraction of snow on the ground(0~1)
    real(r8), pointer :: frac_snow_over_EASS(:)   ! fraction of snow on the overstory[m s-1]
    real(r8), pointer :: frac_snow_under_EASS(:)  ! fraction of snow on the understory[m s-1]
    real(r8), pointer :: mass_rain_over_EASS(:)   ! rainfall mass on the overstory[kg m-2]
    real(r8), pointer :: mass_rain_under_EASS(:)  ! rainfall mass on the understory[kg m-2]
    real(r8), pointer :: mass_snow_grnd_EASS(:)   ! snow mass on the ground[kg m-2]
    real(r8), pointer :: mass_snow_over_EASS(:)   ! snow mass on the overstory[kg m-2]
    real(r8), pointer :: mass_snow_under_EASS(:)  ! snow mass on the understory[kg m-2]
    real(r8), pointer :: qflx_snow_grnd_pft_EASS(:)    !snow rate on the grond at pft-level  [m s-1]
    real(r8), pointer :: qflx_snow_grnd_col_EASS(:)    !snow rate on the grond at column-level  [m s-1]
    real(r8), pointer :: sum_rain_over_EASS(:)    ! area of overstory covered by rain [m2 m-2]
    real(r8), pointer :: sum_rain_under_EASS(:)   ! area of understory covered by rain [m2 m-2]
    real(r8), pointer :: sum_snow_over_EASS(:)    ! area of overstory covered by snow [m2 m-2]
    real(r8), pointer :: sum_snow_under_EASS(:)   ! area of understory covered by snow [m2 m-2]
    real(r8), pointer :: qflx_prec_intr_EASS(:)   ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd_EASS(:)   ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snwcp_liq_EASS(:)   ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice_EASS(:)   ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_rain_grnd_EASS(:)   ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: zi_EASS(:,:)             ! interface level below a "z" level (m)
    real(r8), pointer :: dz_EASS(:,:)             ! layer depth (m)
    real(r8), pointer :: z_EASS(:,:)              ! layer thickness (m)
    real(r8), pointer :: t_soisno_EASS(:,:)       ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice_EASS(:,:)     ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq_EASS(:,:)     ! liquid water (kg/m2)
    real(r8), pointer :: frac_iceold_EASS(:,:)    ! fraction of ice relative to the tot water

! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: f                            ! filter index
    integer  :: pi                           ! pft index
    integer  :: p                            ! pft index
    integer  :: c                            ! column index
    integer  :: l                            ! landunit index
    integer  :: g                            ! gridcell index
    integer  :: newnode                      ! flag when new snow node is set, (1=yes, 0=no)
    real(r8) :: dtime                        ! land model time step (sec)
    real(r8) :: h2ocanmx                     ! maximum allowed water on canopy [mm]
    real(r8) :: fpi                          ! coefficient of interception
    real(r8) :: xrun                         ! excess water that exceeds the leaf capacity [mm/s]
    real(r8) :: dz_snowf                     ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(lbp:ubp)            ! frac of precipitation that is snow
    real(r8) :: fracrain(lbp:ubp)            ! frac of precipitation that is rain
    real(r8) :: qflx_candrip(lbp:ubp)        ! rate of canopy runoff and snow falling off canopy [mm/s]
    real(r8) :: qflx_through_rain(lbp:ubp)   ! direct rain throughfall [mm/s]
    real(r8) :: qflx_through_snow(lbp:ubp)   ! direct snow throughfall [mm/s]
    real(r8) :: qflx_prec_grnd_snow(lbp:ubp) ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain(lbp:ubp) ! rain precipitation incident on ground [mm/s]
    !--------------------------
    !added by Jing Chen, 17 May 2012
    real(r8) :: dmass_rain_over_EASS        ! mass difference of rain fall between two time steps, overstory[kg m-2]
    real(r8) :: dmass_rain_under_EASS       ! mass difference of rain fall between two time steps, understory[kg m-2]
    real(r8) :: dmass_snow_over_EASS        ! mass difference of snow between two time steps, overstory[kg m-2]
    real(r8) :: dmass_snow_under_EASS       ! mass difference of snow between two time steps, understory[kg m-2]
    real(r8) :: mass_rain_over_EASS_bef     ! rainfall mass on the overstory at step (p-1)[kg m-2]
    real(r8) :: mass_rain_under_EASS_bef    ! rainfall mass on the understory at step (p-1)[kg m-2] 
    real(r8) :: mass_snow_over_EASS_bef     ! snow mass on the overstory at step (p-1)[kg m-2]
    real(r8) :: mass_snow_under_EASS_bef    ! snow mass on the understory at step (p-1)[kg m-2]
    real(r8) :: rain_over_EASS(lbp:ubp)     !rain rate , overstory [m s-1](=forc_rain(:)/1000)
    real(r8) :: rain_under_EASS(lbp:ubp)    !rain rate at p,understory [m s-1]
    real(r8) :: snow_over_EASS(lbp:ubp)     !snow rate ,overstory [m s-1]
    real(r8) :: snow_under_EASS(lbp:ubp)    !snow rate ,understory  [m s-1](=forc_snow(:)/1000)
    real(r8) :: dz_snowf_EASS               ! layer thickness rate change due to precipitation [mm/s] 
    real(r8) :: qflx_prec_grnd_snow_EASS(lbp:ubp) ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain_EASS(lbp:ubp) ! rain precipitation incident on ground [mm/s]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    pgridcell          => clm3%g%l%c%p%gridcell
    forc_rain          => clm_a2l%forc_rain
    forc_snow          => clm_a2l%forc_snow

    ! Assign local pointers to derived type members (landunit-level)

    ltype              => clm3%g%l%itype

    ! Assign local pointers to derived type members (column-level)

    cgridcell          => clm3%g%l%c%gridcell
    clandunit          => clm3%g%l%c%landunit
    ctype              => clm3%g%l%c%itype
    pfti               => clm3%g%l%c%pfti
    npfts              => clm3%g%l%c%npfts
    do_capsnow         => clm3%g%l%c%cps%do_capsnow
    forc_t             => clm3%g%l%c%ces%forc_t
    t_grnd             => clm3%g%l%c%ces%t_grnd
    snl                => clm3%g%l%c%cps%snl
    snowdp             => clm3%g%l%c%cps%snowdp
    h2osno             => clm3%g%l%c%cws%h2osno
    zi                 => clm3%g%l%c%cps%zi
    dz                 => clm3%g%l%c%cps%dz
    z                  => clm3%g%l%c%cps%z
    frac_iceold        => clm3%g%l%c%cps%frac_iceold
    t_soisno           => clm3%g%l%c%ces%t_soisno
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    qflx_snow_grnd_col => clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd
    h2ocan_loss        => clm3%g%l%c%cwf%h2ocan_loss
    snw_rds            => clm3%g%l%c%cps%snw_rds
    mss_bcpho          => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi          => clm3%g%l%c%cps%mss_bcphi
    mss_bctot          => clm3%g%l%c%cps%mss_bctot
    mss_bc_col         => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top         => clm3%g%l%c%cps%mss_bc_top
    mss_ocpho          => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi          => clm3%g%l%c%cps%mss_ocphi
    mss_octot          => clm3%g%l%c%cps%mss_octot
    mss_oc_col         => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top         => clm3%g%l%c%cps%mss_oc_top
    mss_dst1           => clm3%g%l%c%cps%mss_dst1
    mss_dst2           => clm3%g%l%c%cps%mss_dst2
    mss_dst3           => clm3%g%l%c%cps%mss_dst3
    mss_dst4           => clm3%g%l%c%cps%mss_dst4
    mss_dsttot         => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col        => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top        => clm3%g%l%c%cps%mss_dst_top
    !--------------------------------------
    !added by Jing Chen, 17 May 2012
    densnow_EASS          => clm3%g%l%c%cwf%densnow_EASS
    snowdp_EASS           => clm3%g%l%c%cps%snowdp_EASS
    frac_snow_grnd_EASS   => clm3%g%l%c%cps%frac_snow_grnd_EASS
    mass_snow_grnd_EASS   => clm3%g%l%c%cws%mass_snow_grnd_EASS
    qflx_snow_grnd_col_EASS => clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd_EASS
    snl_EASS              => clm3%g%l%c%cps%snl_EASS
    h2osno_EASS           => clm3%g%l%c%cws%h2osno_EASS
    zi_EASS               => clm3%g%l%c%cps%zi_EASS
    dz_EASS               => clm3%g%l%c%cps%dz_EASS
    z_EASS                => clm3%g%l%c%cps%z_EASS
    frac_iceold_EASS      => clm3%g%l%c%cps%frac_iceold_EASS
    t_soisno_EASS         => clm3%g%l%c%ces%t_soisno_EASS
    h2osoi_ice_EASS       => clm3%g%l%c%cws%h2osoi_ice_EASS
    h2osoi_liq_EASS       => clm3%g%l%c%cws%h2osoi_liq_EASS
    t_grnd_EASS           => clm3%g%l%c%ces%t_grnd_EASS

    !--------------------------------------

    ! Assign local pointers to derived type members (pft-level)

    plandunit          => clm3%g%l%c%p%landunit
    pcolumn            => clm3%g%l%c%p%column
    dewmx              => clm3%g%l%c%p%pps%dewmx
    frac_veg_nosno     => clm3%g%l%c%p%pps%frac_veg_nosno
    elai               => clm3%g%l%c%p%pps%elai
    esai               => clm3%g%l%c%p%pps%esai
    h2ocan             => clm3%g%l%c%p%pws%h2ocan
    qflx_prec_intr     => clm3%g%l%c%p%pwf%qflx_prec_intr
    qflx_prec_grnd     => clm3%g%l%c%p%pwf%qflx_prec_grnd
    qflx_snwcp_liq     => clm3%g%l%c%p%pwf%qflx_snwcp_liq
    qflx_snwcp_ice     => clm3%g%l%c%p%pwf%qflx_snwcp_ice
    qflx_snow_grnd_pft => clm3%g%l%c%p%pwf%qflx_snow_grnd
    qflx_rain_grnd     => clm3%g%l%c%p%pwf%qflx_rain_grnd
    fwet               => clm3%g%l%c%p%pps%fwet
    fdry               => clm3%g%l%c%p%pps%fdry
    irrig_rate         => clm3%g%l%c%cps%irrig_rate
    n_irrig_steps_left => clm3%g%l%c%cps%n_irrig_steps_left
    qflx_irrig         => clm3%g%l%c%cwf%qflx_irrig
    !-----------------------------------
    !added by Jing Chen, 17 May 2012
    frac_rain_over_EASS   => clm3%g%l%c%p%pps%frac_rain_over_EASS
    frac_rain_under_EASS  => clm3%g%l%c%p%pps%frac_rain_under_EASS
    frac_snow_over_EASS   => clm3%g%l%c%p%pps%frac_snow_over_EASS
    frac_snow_under_EASS  => clm3%g%l%c%p%pps%frac_snow_under_EASS
    mass_rain_over_EASS   => clm3%g%l%c%p%pws%mass_rain_over_EASS
    mass_rain_under_EASS  => clm3%g%l%c%p%pws%mass_rain_under_EASS
    mass_snow_over_EASS   => clm3%g%l%c%p%pws%mass_snow_over_EASS
    mass_snow_under_EASS  => clm3%g%l%c%p%pws%mass_snow_under_EASS
    qflx_snow_grnd_pft_EASS=> clm3%g%l%c%p%pwf%qflx_snow_grnd_EASS
    sum_rain_over_EASS    => clm3%g%l%c%p%pwf%sum_rain_over_EASS
    sum_rain_under_EASS   => clm3%g%l%c%p%pwf%sum_rain_under_EASS
    sum_snow_over_EASS    => clm3%g%l%c%p%pwf%sum_snow_over_EASS
    sum_snow_under_EASS   => clm3%g%l%c%p%pwf%sum_snow_under_EASS
    clumping_EASS         => clm3%g%l%c%p%pps%clumping_EASS
    elai_over_EASS        => clm3%g%l%c%p%pps%elai_over_EASS
    elai_under_EASS       => clm3%g%l%c%p%pps%elai_under_EASS
    h2ocan_EASS           => clm3%g%l%c%p%pws%h2ocan_EASS
    qflx_prec_intr_EASS   => clm3%g%l%c%p%pwf%qflx_prec_intr_EASS
    qflx_prec_grnd_EASS   => clm3%g%l%c%p%pwf%qflx_prec_grnd_EASS
    qflx_snwcp_liq_EASS   => clm3%g%l%c%p%pwf%qflx_snwcp_liq_EASS
    qflx_snwcp_ice_EASS   => clm3%g%l%c%p%pwf%qflx_snwcp_ice_EASS
    qflx_rain_grnd_EASS   => clm3%g%l%c%p%pwf%qflx_rain_grnd_EASS
    !-------------------------------------

    ! Compute time step

    dtime = get_step_size()

    ! Start pft loop

    do f = 1, num_nolakep
       p = filter_nolakep(f)
       g = pgridcell(p)
       l = plandunit(p)
       c = pcolumn(p)

       ! Canopy interception and precipitation onto ground surface
       ! Add precipitation to leaf water

       if (ltype(l)==istsoil .or. ltype(l)==istwet .or. ltype(l)==isturb .or. &
           ltype(l)==istcrop) then

          qflx_candrip(p) = 0._r8      ! rate of canopy runoff
          qflx_through_snow(p) = 0._r8 ! rain precipitation direct through canopy
          qflx_through_rain(p) = 0._r8 ! snow precipitation direct through canopy
          qflx_prec_intr(p) = 0._r8    ! total intercepted precipitation
          fracsnow(p) = 0._r8          ! fraction of input precip that is snow
          fracrain(p) = 0._r8          ! fraction of input precip that is rain
          !---------------
          ! added by Jing Chen Oct 20 2012
          qflx_prec_intr_EASS(p) = 0._r8    ! total intercepted precipitation
          bifall = 0.001_r8

          rain_over_EASS(p)      = 0._r8
          mass_rain_over_EASS(p) = 0._r8
          frac_rain_over_EASS(p) = 0._r8
          sum_rain_over_EASS(p)  = 0._r8
          dmass_rain_over_EASS   = 0._r8

          rain_under_EASS(p)     = 0._r8
          mass_rain_under_EASS(p)= 0._r8
          frac_rain_under_EASS(p)= 0._r8
          sum_rain_under_EASS(p) = 0._r8
          dmass_rain_under_EASS  = 0._r8

          snow_over_EASS(p)      = 0.0_r8
          mass_snow_over_EASS(p) = 0.0_r8
          frac_snow_over_EASS(p) = 0.0_r8
          sum_snow_over_EASS(p)  = 0.0_r8
          dmass_snow_over_EASS   = 0._r8

          snow_under_EASS(p)     = 0.0_r8
          mass_snow_under_EASS(p)= 0.0_r8
          frac_snow_under_EASS(p)= 0.0_r8
          sum_snow_under_EASS(p) = 0.0_r8
          dmass_snow_under_EASS  = 0._r8

          !-----------------

          if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then
             ! canopy covered by percepitaion
             if (frac_veg_nosno(p) == 1 .and. (forc_rain(g) + forc_snow(g)) > 0._r8) then   

                ! determine fraction of input precipitation that is snow and rain

                fracsnow(p) = forc_snow(g)/(forc_snow(g) + forc_rain(g))
                fracrain(p) = forc_rain(g)/(forc_snow(g) + forc_rain(g))

                ! The leaf water capacities for solid and liquid are different,
                ! generally double for snow, but these are of somewhat less
                ! significance for the water budget because of lower evap. rate at
                ! lower temperature.  Hence, it is reasonable to assume that
                ! vegetation storage of solid water is the same as liquid water.
                h2ocanmx = dewmx(p) * (elai(p) + esai(p))     ![mm]
                
                ! Coefficient of interception
                ! set fraction of potential interception to max 0.25
                fpi = 0.25_r8*(1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))
                
                ! Direct throughfall[mm/s]
                qflx_through_snow(p) = forc_snow(g) * (1._r8-fpi)
                qflx_through_rain(p) = forc_rain(g) * (1._r8-fpi)
                
                ! Intercepted precipitation [mm/s]
                qflx_prec_intr(p) = (forc_snow(g) + forc_rain(g)) * fpi
                
                ! Water storage of intercepted precipitation and dew [mm H2O]
                h2ocan(p) = max(0._r8, h2ocan(p) + dtime*qflx_prec_intr(p))
 
                ! Initialize rate of canopy runoff and snow falling off canopy[mm/s]
                qflx_candrip(p) = 0._r8
                
                ! Excess water that exceeds the leaf capacity[mm/s]
                xrun = (h2ocan(p) - h2ocanmx)/dtime
                
                ! Test on maximum dew on leaf
                ! Note if xrun > 0 then h2ocan must be at least h2ocanmx
                if (xrun > 0._r8) then
                   qflx_candrip(p) = xrun
                   h2ocan(p) = h2ocanmx
                end if

                !-------------------------
                ! added by Jing Chen, 17 May 2012
                if (ltype(l)==istsoil .or. ltype(l)==istcrop ) then

                   ! fraction of rainfall and snow on canopy including overstory and understroy
                   ! used to calculate evaporation from liqid water and snow of leaves respectly 
                   ! in EASS
                   rain_over_EASS(p)=forc_rain(g)/1000.0_r8    ! rain rate , unit change[mm s-1]->[m s-1]
                   snow_over_EASS(p)=forc_snow(g)/1000.0_r8    ! snow rate , unit change[mm s-1]->[m s-1]

                   !(1) rainfall on the canopy
                   if (rain_over_EASS(p) .gt. 1.0e-10_r8 ) then
                        !(1.1) rainfall on the overstroy
                        mass_rain_over_EASS_bef=mass_rain_over_EASS(p)

                        call CanopyFallFrac_EASS(dtime, clumping_EASS(p), elai_over_EASS(p) , &         ! input var.
                                              rain_over_EASS(p), denh2o, mass_rain_over_EASS_bef, &  ! input var.
                                              mass_rain_over_EASS(p), frac_rain_over_EASS(p), &      ! output var.
                                              sum_rain_over_EASS(p))                                 ! output var.

                        !(1.2) rainfall on the understroy
                        mass_rain_under_EASS_bef=mass_rain_under_EASS(p)
                        dmass_rain_over_EASS=mass_rain_over_EASS(p)-mass_rain_under_EASS_bef   ![kg m-2]
                        dmass_rain_over_EASS=max(0.0_r8, dmass_rain_over_EASS)
                        rain_under_EASS(p)=rain_over_EASS(p)-dmass_rain_over_EASS/denh2o/dtime  ![m s-1](=kg m-2 kg-1 m3 s-1)

                        call CanopyFallFrac_EASS(dtime, clumping_EASS(p),elai_under_EASS(p), &           ! input var.
                                              rain_under_EASS(p), denh2o, mass_rain_under_EASS_bef, & ! input var. 
                                              mass_rain_under_EASS(p), frac_rain_under_EASS(p),  &    ! output var.
                                              sum_rain_under_EASS(p))                                 ! output var.

                        dmass_rain_under_EASS=mass_rain_under_EASS(p)-mass_rain_under_EASS_bef
                        dmass_rain_under_EASS=max(0.0_r8, dmass_rain_under_EASS)

                   end if

                   !(2)snow on the canopy
                   if ( snow_over_EASS(p) .gt. 1.0e-10_r8 .and. frac_veg_nosno(p) == 1 ) then

                        ! snow density on the overstroy based on temperature
                        if ( clm_a2l%forc_t(g) .gt. (tfrz + 2.0_r8) ) then
                            bifall=50.0_r8 + 1.7_r8*(17.0_r8)**1.5_r8
                        else if (clm_a2l%forc_t(g) .gt. (tfrz - 15.0_r8) ) then
                            bifall=50.0_r8 + 1.7_r8*(clm_a2l%forc_t(g) - tfrz + 15.0_r8)**1.5_r8
                        else
                            bifall=50.0_r8
                        end if

                        !(2.1) snow on the overstroy
                        mass_snow_over_EASS_bef=mass_snow_over_EASS(p)
                        call CanopyFallFrac_EASS(dtime, clumping_EASS(p),elai_over_EASS(p), &          ! input var.
                                                 snow_over_EASS(p), bifall, mass_snow_over_EASS_bef, & ! input var.
                                                 mass_snow_over_EASS(p), frac_snow_over_EASS(p), &     ! output var.
                                                 sum_snow_over_EASS(p))      ! output var.

                        dmass_snow_over_EASS=mass_snow_over_EASS(p)-mass_snow_over_EASS_bef
                        dmass_snow_over_EASS=max(0.0_r8, dmass_snow_over_EASS)

                        !(2.2) snow on the understroy
                        mass_snow_under_EASS_bef=mass_snow_under_EASS(p)
                        snow_under_EASS(p)=snow_over_EASS(p)-dmass_snow_over_EASS/bifall/dtime  ![m s-1]
                        call CanopyFallFrac_EASS(dtime, clumping_EASS(p),elai_under_EASS(p), &             ! input var.
                                              snow_under_EASS(p), bifall, mass_snow_under_EASS_bef,  &  ! input var.
                                              mass_snow_under_EASS(p), frac_snow_under_EASS(p),&        ! output var.
                                              sum_snow_under_EASS(p))      ! output var.

                        dmass_snow_under_EASS=mass_snow_under_EASS(p)-mass_snow_under_EASS_bef
                        dmass_snow_under_EASS=max(0.0_r8, dmass_snow_under_EASS)

                   end if

                   !---------------------------
                   ! Intercepted precipitation [mm/s]
                   qflx_prec_intr_EASS(p) = ( ( mass_rain_over_EASS(p) + mass_rain_under_EASS(p) ) / denh2o &
                                       +   ( mass_snow_over_EASS(p) + mass_snow_under_EASS(p) ) / bifall ) * 1000._r8 / dtime

                   ! Water storage of intercepted precipitation and dew [mm H2O]
                   h2ocan_EASS(p) = max(0._r8, h2ocan_EASS(p) + dtime*qflx_prec_intr_EASS(p))
                else
                   qflx_prec_intr_EASS(p) =qflx_prec_intr(p)
                   h2ocan_EASS(p) = h2ocan(p)
                end if
             end if
             !---------------------------
          end if

       else if (ltype(l)==istice .or. ltype(l)==istice_mec) then

          h2ocan(p)            = 0._r8
          qflx_candrip(p)      = 0._r8
          qflx_through_snow(p) = 0._r8
          qflx_through_rain(p) = 0._r8
          qflx_prec_intr(p)    = 0._r8
          fracsnow(p)          = 0._r8
          fracrain(p)          = 0._r8
          !---------------
          ! added by Jing Chen Oct 20 2012
          h2ocan_EASS(p)       = 0._r8
          qflx_prec_intr_EASS(p) = 0._r8
          !----------------

       end if

       ! Precipitation onto ground (kg/(m2 s))
       ! PET, 1/18/2005: Added new terms for mass balance correction
       ! due to dynamic pft weight shifting (column-level h2ocan_loss)
       ! Because the fractionation between rain and snow is indeterminate if
       ! rain + snow = 0, I am adding this very small flux only to the rain
       ! components.

       if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then
          if (frac_veg_nosno(p) == 0) then            ! Bare ground
             qflx_prec_grnd_snow(p) = forc_snow(g)                     ![mm/s]
             qflx_prec_grnd_rain(p) = forc_rain(g) + h2ocan_loss(c)    ![mm/s]
             !-----------------------
             ! added by Jing Chen Oct 20 2012
             qflx_prec_grnd_snow_EASS(p) = forc_snow(g)                     ![mm/s]
             qflx_prec_grnd_rain_EASS(p) = forc_rain(g) + h2ocan_loss(c)    ![mm/s]
             !-----------------------
          else
             qflx_prec_grnd_snow(p) = qflx_through_snow(p) + (qflx_candrip(p) * fracsnow(p))
             qflx_prec_grnd_rain(p) = qflx_through_rain(p) + (qflx_candrip(p) * fracrain(p)) + h2ocan_loss(c)
             !-----------------------
             ! added by Jing Chen Oct 20 2012
             ! snow on the ground under the canopy[mm s-1]
             qflx_prec_grnd_snow_EASS(p) = (snow_under_EASS(p)-dmass_snow_under_EASS/bifall/dtime) * 1000._r8   ![mm/s]
             qflx_prec_grnd_rain_EASS(p) = (rain_under_EASS(p)-dmass_rain_under_EASS/denh2o/dtime) * 1000._r8 + h2ocan_loss(c)    ![mm/s]
             !-----------------------
          end if
       ! Urban sunwall and shadewall have no intercepted precipitation
       else
          qflx_prec_grnd_snow(p) = 0._r8
          qflx_prec_grnd_rain(p) = 0._r8
          !------------------
          ! added by Jing Chen Oct 20 2012
          qflx_prec_grnd_snow_EASS(p) = 0._r8
          qflx_prec_grnd_rain_EASS(p) = 0._r8
          !---------------
       end if

       ! Determine whether we're irrigating here; set qflx_irrig appropriately
       if (n_irrig_steps_left(c) > 0) then
          qflx_irrig(c)         = irrig_rate(c)
          n_irrig_steps_left(c) = n_irrig_steps_left(c) - 1
       else
          qflx_irrig(c) = 0._r8
       end if

       ! Add irrigation water directly onto ground (bypassing canopy interception)
       ! Note that it's still possible that (some of) this irrigation water will runoff (as runoff is computed later)
       qflx_prec_grnd_rain(p) = qflx_prec_grnd_rain(p) + qflx_irrig(c)
       !-----------------
       ! added by Jing Chen Oct 20 2012
       qflx_prec_grnd_rain_EASS(p) = qflx_prec_grnd_rain_EASS(p) + qflx_irrig(c)
       !------------------

       ! Done irrigation

       qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)
       !-----------------
       ! added by Jing Chen Oct 20 2012
       qflx_prec_grnd_EASS(p) = qflx_prec_grnd_snow_EASS(p) + qflx_prec_grnd_rain_EASS(p)
       !------------------

       if (do_capsnow(c)) then
          qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
          qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)
          qflx_snow_grnd_pft(p) = 0._r8
          qflx_rain_grnd(p) = 0._r8
          !-------------
          ! added by Jing Chen Oct 20 2012
          qflx_snwcp_liq_EASS(p) = qflx_prec_grnd_rain_EASS(p)
          qflx_snwcp_ice_EASS(p) = qflx_prec_grnd_snow_EASS(p)
          qflx_snow_grnd_pft_EASS(p) = 0._r8
          qflx_rain_grnd_EASS(p) = 0._r8
          !-----------------
       else
          qflx_snwcp_liq(p) = 0._r8
          qflx_snwcp_ice(p) = 0._r8
          qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
          qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
          !-------------
          ! added by Jing Chen Oct 20 2012
          qflx_snwcp_liq_EASS(p) = 0._r8
          qflx_snwcp_ice_EASS(p) = 0._r8
          qflx_snow_grnd_pft_EASS(p) = qflx_prec_grnd_snow_EASS(p)
          qflx_rain_grnd_EASS(p)     = qflx_prec_grnd_rain_EASS(p)
          !-----------------
       end if
    end do ! (end pft loop)

    ! Determine the fraction of foliage covered by water and the
    ! fraction of foliage that is dry and transpiring.

    call FracWet(num_nolakep, filter_nolakep)

    ! Update column level state variables for snow.

    call p2c(num_nolakec, filter_nolakec, qflx_snow_grnd_pft, qflx_snow_grnd_col)
    call p2c(num_nolakec, filter_nolakec, qflx_snow_grnd_pft_EASS, qflx_snow_grnd_col_EASS)

    ! Determine snow height and snow water

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       l = clandunit(c)
       g = cgridcell(c)

       ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
       ! U.S.Department of Agriculture Forest Service, Project F,
       ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

       if (do_capsnow(c)) then
          dz_snowf = 0._r8
          !-------------
          ! added by Jing Chen 0ct 20 2012
          dz_snowf_EASS = 0._r8
          !---------------
       else
          if (forc_t(c) > tfrz + 2._r8) then
             bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
          else if (forc_t(c) > tfrz - 15._r8) then
             bifall=50._r8 + 1.7_r8*(forc_t(c) - tfrz + 15._r8)**1.5_r8
          else
             bifall=50._r8
          end if

          dz_snowf = qflx_snow_grnd_col(c)/bifall
          snowdp(c) = snowdp(c) + dz_snowf*dtime
          h2osno(c) = h2osno(c) + qflx_snow_grnd_col(c)*dtime  ! snow water equivalent (mm)
          !------------------------
          ! added by Jing Chen 0ct 20 2012
          dz_snowf_EASS = qflx_snow_grnd_col_EASS(c)/bifall
          snowdp_EASS(c) = snowdp_EASS(c) + dz_snowf_EASS*dtime
          h2osno_EASS(c) = h2osno_EASS(c) + qflx_snow_grnd_col_EASS(c)*dtime  ! snow water equivalent (mm)
          !------------------------

       end if

       if (ltype(l)==istwet ) then
          if ( t_grnd(c)>tfrz ) then
              h2osno(c)=0._r8
              snowdp(c)=0._r8
          end if
          !------------------------
          ! added by Jing Chen 0ct 20 2012
          if ( t_grnd_EASS(c)>tfrz ) then
              h2osno_EASS(c)=0._r8
              snowdp_EASS(c)=0._r8
          end if
          !-------------------------
       end if

       ! When the snow accumulation exceeds 10 mm, initialize snow layer
       ! Currently, the water temperature for the precipitation is simply set
       ! as the surface air temperature

       newnode = 0    ! flag for when snow node will be initialized
       if (snl(c) == 0 .and. qflx_snow_grnd_col(c) > 0.0_r8 .and. snowdp(c) >= 0.01_r8) then
          newnode         = 1
          snl(c)          = -1
          dz(c,0)         = snowdp(c)                       ! meter
          z(c,0)          = -0.5_r8*dz(c,0)
          zi(c,-1)        = -dz(c,0)
          t_soisno(c,0)   = min(tfrz, forc_t(c))      ! K
          h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
          h2osoi_liq(c,0) = 0._r8                   ! kg/m2
          frac_iceold(c,0) = 1._r8

          ! intitialize SNICAR variables for fresh snow:
          snw_rds(c,0)    = snw_rds_min

          mss_bcpho(c,:)  = 0._r8
          mss_bcphi(c,:)  = 0._r8
          mss_bctot(c,:)  = 0._r8
          mss_bc_col(c)   = 0._r8
          mss_bc_top(c)   = 0._r8

          mss_ocpho(c,:)  = 0._r8
          mss_ocphi(c,:)  = 0._r8
          mss_octot(c,:)  = 0._r8
          mss_oc_col(c)   = 0._r8
          mss_oc_top(c)   = 0._r8

          mss_dst1(c,:)   = 0._r8
          mss_dst2(c,:)   = 0._r8
          mss_dst3(c,:)   = 0._r8
          mss_dst4(c,:)   = 0._r8
          mss_dsttot(c,:) = 0._r8
          mss_dst_col(c)  = 0._r8
          mss_dst_top(c)  = 0._r8
       end if

       ! The change of ice partial density of surface node due to precipitation.
       ! Only ice part of snowfall is added here, the liquid part will be added
       ! later.

       if (snl(c) < 0 .and. newnode == 0) then
          h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+dtime*qflx_snow_grnd_col(c)
          dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
       end if

       !--------------------------
       !added by Jing Chen, 17 May 2012
       ! 1) snow depth on the ground
       newnode = 0    ! flag for when snow node will be initialized
       if (snl_EASS(c) == 0 .and. qflx_snow_grnd_col_EASS(c) > 0.0_r8 .and. snowdp_EASS(c) >= 0.01_r8) then
          newnode              = 1
          snl_EASS(c)          = -1
          dz_EASS(c,0)         = snowdp_EASS(c)            ! meter
          z_EASS(c,0)          = -0.5_r8*dz_EASS(c,0)
          zi_EASS(c,-1)        = -dz_EASS(c,0)
          t_soisno_EASS(c,0)   = min(tfrz, forc_t(c))      ! K
          h2osoi_ice_EASS(c,0) = h2osno_EASS(c)            ! kg/m2
          h2osoi_liq_EASS(c,0) = 0._r8                     ! kg/m2
          frac_iceold_EASS(c,0)= 1._r8
       end if

       if (snl_EASS(c) < 0 .and. newnode == 0) then
          h2osoi_ice_EASS(c,snl_EASS(c)+1) = h2osoi_ice_EASS(c,snl_EASS(c)+1)+dtime*qflx_snow_grnd_col_EASS(c)
          dz_EASS(c,snl_EASS(c)+1) = dz_EASS(c,snl_EASS(c)+1)+dz_snowf_EASS*dtime
       end if

       
        ! 2) snow mass on the ground
        mass_snow_grnd_EASS(c)=mass_snow_grnd_EASS(c)+snowdp_EASS(c)*bifall
        mass_snow_grnd_EASS(c)=max( 0.0_r8, mass_snow_grnd_EASS(c) )

        ! 3) snow fraction on the ground
        ! update the snow density on the ground
        if ( ( dz_EASS(c,snl_EASS(c)) .ge. 1.0e-10_r8 ) .or. ( dz_EASS(c,snl_EASS(c)) .le. -1.0e-10_r8 ) ) then
                densnow_EASS(c)=( densnow_EASS(c)*snowdp_EASS(c) + bifall*dz_EASS(c,snl_EASS(c)) ) &
                               / (snowdp_EASS(c)+dz_EASS(c,snl_EASS(c)) )
        else
                densnow_EASS(c)=( densnow_EASS(c)-250.0_r8) * exp(-0.001_r8*dtime/3600)+250.0_r8
        end if

        ! the fraction of snow on the ground
        frac_snow_grnd_EASS(c)= mass_snow_grnd_EASS(c) / ( 0.05_r8 * densnow_EASS(c) )
        frac_snow_grnd_EASS(c)=min( frac_snow_grnd_EASS(c), 1.0_r8 )
       !--------------------------

    end do

  end subroutine Hydrology1
  
  !----------------------
  ! BOP
  !!INROTINE: CanopyFallFrac_EASS()

  !!INTERFACE:
  subroutine CanopyFallFrac_EASS(dti, clumping, LAI, fall, denfall, fallcan_bef, fallcan, frac_fallcan, sum_fallcan )
  !!USES:

  !!ARGUMENTS:
    implicit none
    real(r8), intent(in) :: dti           
    real(r8), intent(in) :: LAI         ! canopy lai in EASS, including overstroy and understroy
    real(r8), intent(in) :: clumping    ! clumping index
    real(r8), intent(in) :: fall        ! rainfall or snow [m s-1]  
    real(r8), intent(in) :: denfall     ! density of rainfall or snow [kg m-3] 
    real(r8), intent(in) :: fallcan_bef
    real(r8), intent(out) :: fallcan    ! rainfall or snow mass on the canopy[kg m-2]   
    real(r8), intent(out) :: frac_fallcan  ! fraction of rainfall or snow on the canopy [0~1] 
    real(r8), intent(out) :: sum_fallcan   ! area of canopy covered by fall [m2 m-2]   
 
  !! LOCAL VARIABLES:

  !! OTHER LOCAL VARIABLES:
    real(r8) :: chi              ![-]
    real(r8) :: fallcan_max      ! max. value of rainfall or snow on canopy[kg m-2](???why the unit is kg m-2)
    real(r8) :: sum_fallcan_max  ! max. area of canopy covered by snow [m2 m-2]	  

   !-------------------------

    chi=1.0_r8-exp(-1.0_r8*LAI*clumping)

    fallcan_max=0.1_r8*LAI
    fallcan=fallcan_bef + fall * chi * denfall * dti

    fallcan=max(0.0_r8, fallcan)
    fallcan=min(fallcan,fallcan_max)

    frac_fallcan=min(1.0_r8, fallcan/fallcan_max)

    sum_fallcan_max=0.01_r8*LAI
    sum_fallcan=frac_fallcan*sum_fallcan_max

  end subroutine CanopyFallFrac_EASS

end module Hydrology1Mod
