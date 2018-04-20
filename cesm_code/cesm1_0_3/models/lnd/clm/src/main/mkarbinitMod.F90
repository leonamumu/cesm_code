module mkarbinitMod
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: mkarbinitMod
!
! !DESCRIPTION:
!
!
!---------------------------------------------------------------------------

! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_varctl   , only : iulog
    use shr_sys_mod  , only : shr_sys_flush
    use spmdMod      , only : masterproc

    implicit none

    SAVE
    private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

    public mkarbinit   ! Make arbitrary initial conditions
    public perturbIC   ! Perturb the initial conditions by pertlim

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkarbinit
!
! !INTERFACE:
  subroutine mkarbinit()
!
! !DESCRIPTION:
! Initializes the following time varying variables:
! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow       : snowdp, snl, dz, z, zi
! temperature: t_soisno, t_veg, t_grnd
!
! !USES:
    use shr_const_mod, only : SHR_CONST_TKFRZ
    use clmtype
    use clm_varpar   , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon   , only : bdsno, istice, istwet, istsoil, isturb, &
                              denice, denh2o, spval, sb, icol_road_perv, &
                              icol_road_imperv, icol_roof, icol_sunwall, &
                              icol_shadewall
    use clm_varcon   , only : istcrop
    use clm_varcon   , only : istice_mec, h2osno_max
    use clm_varctl   , only : iulog, pertlim
    use spmdMod      , only : masterproc
    use decompMod    , only : get_proc_bounds
    use SNICARMod    , only : snw_rds_min
    use clm_atmlnd   , only : clm_a2l
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 3/07/08 Keith Oleson: initialize h2osoi_vol for all soil layers to 0.3
! 3/18/08 David Lawrence, initialize deep layers
! 03/28/08 Mark Flanner, initialize snow aerosols and grain size
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)        ! column index associated with each pft
    integer , pointer :: ctype(:)          ! column type
    integer , pointer :: clandunit(:)      ! landunit index associated with each column
    integer , pointer :: ltype(:)          ! landunit type
    logical , pointer :: lakpoi(:)         ! true => landunit is a lake point
    integer , pointer :: plandunit(:)      ! landunit index associated with each pft
    logical , pointer :: urbpoi(:)         ! true => landunit is an urban point
    logical , pointer :: ifspecial(:)      ! true => landunit is not vegetated
    real(r8), pointer :: dz(:,:)           ! layer thickness depth (m)
    real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: bsw2(:,:)         ! Clapp and Hornberger "b" for CN code
    real(r8), pointer :: psisat(:,:)       ! soil water potential at saturation for CN code (MPa)
    real(r8), pointer :: vwcsat(:,:)       ! volumetric water content at saturation for CN code (m3/m3)
    real(r8), pointer :: zi(:,:)           ! interface level below a "z" level (m)
    real(r8), pointer :: wa(:)             ! water in the unconfined aquifer (mm)
    real(r8), pointer :: wt(:)             ! total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: zwt(:)            ! water table depth (m)
    !--------------------------
    !added by Jing Chen, 18 May 2012
    integer , pointer :: ivt(:)            ! pft vegetation type
    real(r8), pointer :: albsnow_EASS(:,:) ! albedo of new snow
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature at the column level [Kelvin]
    integer , pointer :: pgridcell(:)      ! gridcell index associated with each pft
    real(r8), pointer :: forc_pco2(:)      ! partial pressure co2 [Pa]
!-------------------------
!
! local pointers to implicit out arguments
!
    integer , pointer :: snl(:)             ! number of snow layers
    real(r8), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_lake(:,:)        ! lake temperature (Kelvin)  (1:nlevlak)
    real(r8), pointer :: t_grnd(:)          ! ground temperature (Kelvin)
    real(r8), pointer :: t_veg(:)           ! vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m(:)         ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_u(:)       ! Urban 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_r(:)       ! Rural 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: h2osoi_vol(:,:)    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: h2ocan_col(:)      ! canopy water (mm H2O) (column-level)
    real(r8), pointer :: h2ocan_pft(:)      ! canopy water (mm H2O) (pft-level)
    real(r8), pointer :: h2osno(:)          ! snow water (mm H2O)
    real(r8), pointer :: snowdp(:)          ! snow height (m)
    real(r8), pointer :: qflx_irrig(:)      ! irrigation flux (mm H2O/s)
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: soilpsi(:,:)       ! soil water potential in each soil layer (MPa)
    real(r8), pointer :: snw_rds(:,:)       ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: snw_rds_top(:)     ! snow grain size, top (col) [microns]
    real(r8), pointer :: sno_liq_top(:)     ! liquid water fraction (mass) in top snow layer (col) [frc]
    real(r8), pointer :: mss_bcpho(:,:)     ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcphi(:,:)     ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bctot(:,:)     ! total mass of BC (pho+phi) (col,lyr) [kg]
    real(r8), pointer :: mss_bc_col(:)      ! total mass of BC in snow column (col) [kg]
    real(r8), pointer :: mss_bc_top(:)      ! total mass of BC in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_bcphi(:,:) ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_bcpho(:,:) ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_ocpho(:,:)     ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)     ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_octot(:,:)     ! total mass of OC (pho+phi) (col,lyr) [kg]
    real(r8), pointer :: mss_oc_col(:)      ! total mass of OC in snow column (col) [kg]
    real(r8), pointer :: mss_oc_top(:)      ! total mass of OC in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_ocphi(:,:) ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocpho(:,:) ! mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_dst1(:,:)      ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)      ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)      ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)      ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dsttot(:,:)    ! total mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_col(:)     ! total mass of dust in snow column (col) [kg]
    real(r8), pointer :: mss_dst_top(:)     ! total mass of dust in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_dst1(:,:)  ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst2(:,:)  ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst3(:,:)  ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst4(:,:)  ! mass concentration of dust species 4 (col,lyr) [kg/kg]
    real(r8), pointer :: irrig_rate(:)         ! current irrigation rate [mm/s]
    integer,  pointer :: n_irrig_steps_left(:) ! number of time steps for which we still need to irrigate today (if 0, ignore irrig_rate)
 
    !--------------------------
    !added by Jing Chen, 17 May 2012
    real(r8), pointer :: mass_rain_over_EASS(:)   ! rainfall mass on the overstory[kg m-2]
    real(r8), pointer :: mass_rain_under_EASS(:)  ! rainfall mass on the understory[kg m-2]
    real(r8), pointer :: mass_snow_over_EASS(:)   ! snow mass on the overstory[kg m-2]
    real(r8), pointer :: mass_snow_under_EASS(:)  ! snow mass on the understory[kg m-2]
    real(r8), pointer :: frac_rain_over_EASS(:)   ! fraction of rainfall on the overstory(0~1)
    real(r8), pointer :: frac_rain_under_EASS(:)  ! fraction of rainfall on the understory(0~1)
    real(r8), pointer :: frac_snow_grnd_EASS(:)   ! fraction of snow on the ground(0~1)
    real(r8), pointer :: frac_snow_over_EASS(:)   ! fraction of snow on the overstory[m s-1]
    real(r8), pointer :: frac_snow_under_EASS(:)  ! fraction of snow on the understory(0~1)
    real(r8), pointer :: sum_rain_over_EASS(:)    ! area of overstory covered by rain [m2 m-2]
    real(r8), pointer :: sum_rain_under_EASS(:)   ! area of understory covered by rain [m2 m-2]
    real(r8), pointer :: sum_snow_over_EASS(:)    ! area of overstory covered by snow [m2 m-2]
    real(r8), pointer :: sum_snow_under_EASS(:)   ! area of understory covered by snow [m2 m-2]
    real(r8), pointer :: sum_snow_EASS(:)         ! accumulation of snow mass base on forc_snow[m s-1]
    real(r8), pointer :: densnow_EASS(:)          ! snow denisity on the ground [kg m-3]
    real(r8), pointer :: mass_snow_grnd_EASS(:)   ! snow mass on the ground[kg m-2]
    real(r8), pointer :: snowdp_EASS(:)           ! snow height [m]
    real(r8), pointer :: ponddp_EASS(:)           ! water pond depth[m]
    real(r8), pointer :: albedo_snow_EASS(:)      ! albedo of new snow based on pft
    real(r8), pointer :: Tleaf_over_sunlit_EASS(:)   ! leaf temperature of sunlit leaves, overstory [K]
    real(r8), pointer :: Tleaf_under_sunlit_EASS(:)  ! leaf temperature of sunlit leaves, understory [K]
    real(r8), pointer :: Tleaf_over_shaded_EASS(:)   ! leaf temperature of shaded leaves, overstory [K]
    real(r8), pointer :: Tleaf_under_shaded_EASS(:)  ! leaf temperature of shaded leaves, understory [K]
    real(r8), pointer :: eflx_sh_over_EASS(:)        ! sensible heat flux, overstory [w m-2]
    real(r8), pointer :: ci_over_shaded_EASS(:)      ! intercellular co2 concentration of shaded leaves,overstory [pa]
    real(r8), pointer :: ci_over_sunlit_EASS(:)      ! intercellular co2 concentration of sunlit leaves,overstory [pa]
    real(r8), pointer :: ci_under_shaded_EASS(:)     ! intercellular co2 concentration of shaded leaves,understory [pa]
    real(r8), pointer :: ci_under_sunlit_EASS(:)     ! intercellular co2 concentration of sunlit leaves,understory [pa]
    real(r8), pointer :: cs_over_shaded_EASS(:)      ! co2 concentration on the surfaces of shaded leaves, overstory[pa]
    real(r8), pointer :: cs_over_sunlit_EASS(:)      ! co2 concentration on the surfaces of sunlit leaves, overstory[pa]
    real(r8), pointer :: cs_under_shaded_EASS(:)     ! co2 concentration on the surfaces of shaded leaves, understory[pa]
    real(r8), pointer :: cs_under_sunlit_EASS(:)     ! co2 concentration on the surfaces of sunlit leaves, understory[pa]
    real(r8), pointer :: psn_over_shaded_EASS(:)     ! net photosynthesis rate for shaded leaves, overstroy[W m-2]
    real(r8), pointer :: psn_over_sunlit_EASS(:)     ! net photosynthesis rate for sunlit leaves, overstroy[W m-2]
    real(r8), pointer :: psn_under_shaded_EASS(:)    ! net photosynthesis rate for shaded leaves, understroy[W m-2]
    real(r8), pointer :: psn_under_sunlit_EASS(:)    ! net photosynthesis rate for sunlit leaves, understroy[W m-2]
    real(r8), pointer :: rs_over_shaded_EASS(:)      ! stomatal resistant of shaded leaves, overstroy[s m-1]
    real(r8), pointer :: rs_over_sunlit_EASS(:)      ! stomatal resistant of sunlit leaves, overstroy[s m-1]
    real(r8), pointer :: rs_under_shaded_EASS(:)     ! stomatal resistant of shaded leaves, understroy[s m-1]
    real(r8), pointer :: rs_under_sunlit_EASS(:)     ! stomatal resistant of sunlit leaves, understroy[s m-1]
    real(r8), pointer :: h2osoi_ice_EASS(:,:)        ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq_EASS(:,:)        ! liquid water (kg/m2)
    integer , pointer :: snl_EASS(:)             ! number of snow layers
    real(r8), pointer :: h2ocan_col_EASS(:)      ! canopy water (mm H2O) (column-level)
    real(r8), pointer :: h2ocan_pft_EASS(:)      ! canopy water (mm H2O) (pft-level)
    real(r8), pointer :: h2osno_EASS(:)          ! snow water (mm H2O)
    real(r8), pointer :: h2osoi_vol_EASS(:,:)    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: t_soisno_EASS(:,:)      ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_grnd_EASS(:)          ! ground temperature (Kelvin)
    real(r8), pointer :: wa_EASS(:)              ! water in the unconfined aquifer (mm)
    real(r8), pointer :: wt_EASS(:)              ! total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: zwt_EASS(:)             ! water table depth (m)

    !-------------------------
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: j,l,c,p,g    ! indices
    integer :: nlevs        ! number of levels
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8):: vwc,psi      ! for calculating soilpsi
!-----------------------------------------------------------------------

    if ( masterproc )then
        write(iulog,*) 'Setting initial data to non-spun up values'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype
    lakpoi     => clm3%g%l%lakpoi
    ifspecial  => clm3%g%l%ifspecial
    urbpoi     => clm3%g%l%urbpoi

    ! Assign local pointers to derived subtypes components (column-level)

    ctype            => clm3%g%l%c%itype
    clandunit        => clm3%g%l%c%landunit
    snl              => clm3%g%l%c%cps%snl
    dz               => clm3%g%l%c%cps%dz
    watsat           => clm3%g%l%c%cps%watsat
    bsw2             => clm3%g%l%c%cps%bsw2
    vwcsat           => clm3%g%l%c%cps%vwcsat
    psisat           => clm3%g%l%c%cps%psisat
    soilpsi          => clm3%g%l%c%cps%soilpsi
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol       => clm3%g%l%c%cws%h2osoi_vol
    h2ocan_col       => clm3%g%l%c%cws%pws_a%h2ocan
    qflx_irrig       => clm3%g%l%c%cwf%qflx_irrig
    snowdp           => clm3%g%l%c%cps%snowdp
    h2osno           => clm3%g%l%c%cws%h2osno
    t_soisno         => clm3%g%l%c%ces%t_soisno
    t_lake           => clm3%g%l%c%ces%t_lake
    t_grnd           => clm3%g%l%c%ces%t_grnd
    zi               => clm3%g%l%c%cps%zi
    wa               => clm3%g%l%c%cws%wa
    wt               => clm3%g%l%c%cws%wt
    zwt              => clm3%g%l%c%cws%zwt
    snw_rds          => clm3%g%l%c%cps%snw_rds
    snw_rds_top      => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top      => clm3%g%l%c%cps%sno_liq_top
    mss_bcpho        => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi        => clm3%g%l%c%cps%mss_bcphi
    mss_bctot        => clm3%g%l%c%cps%mss_bctot
    mss_bc_col       => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top       => clm3%g%l%c%cps%mss_bc_top
    mss_cnc_bcphi    => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho    => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_ocpho        => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi        => clm3%g%l%c%cps%mss_ocphi
    mss_octot        => clm3%g%l%c%cps%mss_octot
    mss_oc_col       => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top       => clm3%g%l%c%cps%mss_oc_top
    mss_cnc_ocphi    => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho    => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_dst1         => clm3%g%l%c%cps%mss_dst1
    mss_dst2         => clm3%g%l%c%cps%mss_dst2
    mss_dst3         => clm3%g%l%c%cps%mss_dst3
    mss_dst4         => clm3%g%l%c%cps%mss_dst4
    mss_dsttot       => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col      => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top      => clm3%g%l%c%cps%mss_dst_top
    mss_cnc_dst1     => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2     => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3     => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4     => clm3%g%l%c%cps%mss_cnc_dst4
    n_irrig_steps_left => clm3%g%l%c%cps%n_irrig_steps_left
    irrig_rate       => clm3%g%l%c%cps%irrig_rate
    !-----------------------
    ! added by Jing Chen, 17 May 2012
    densnow_EASS          => clm3%g%l%c%cwf%densnow_EASS
    mass_snow_grnd_EASS   => clm3%g%l%c%cws%mass_snow_grnd_EASS
    snowdp_EASS           => clm3%g%l%c%cps%snowdp_EASS
    ponddp_EASS           => clm3%g%l%c%cps%ponddp_EASS
    albedo_snow_EASS      => pftcon%albedo_snow_EASS
    albsnow_EASS          => clm3%g%l%c%cps%albsnow_EASS
    h2osoi_ice_EASS       => clm3%g%l%c%cws%h2osoi_ice_EASS
    h2osoi_liq_EASS       => clm3%g%l%c%cws%h2osoi_liq_EASS
    h2osoi_vol_EASS       => clm3%g%l%c%cws%h2osoi_vol_EASS
    h2osno_EASS           => clm3%g%l%c%cws%h2osno_EASS
    snl_EASS              => clm3%g%l%c%cps%snl_EASS
    h2ocan_col_EASS       => clm3%g%l%c%cws%pws_a%h2ocan_EASS
    t_soisno_EASS         => clm3%g%l%c%ces%t_soisno_EASS
    t_grnd_EASS           => clm3%g%l%c%ces%t_grnd_EASS
    wa_EASS               => clm3%g%l%c%cws%wa_EASS
    wt_EASS               => clm3%g%l%c%cws%wt_EASS
    zwt_EASS              => clm3%g%l%c%cws%zwt_EASS

    ! Assign local pointers to derived type members (gridcell-level)
    forc_pco2      => clm_a2l%forc_pco2
    forc_t         => clm_a2l%forc_t
    !-------------------------

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn        => clm3%g%l%c%p%column
    h2ocan_pft     => clm3%g%l%c%p%pws%h2ocan
    t_veg          => clm3%g%l%c%p%pes%t_veg
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    t_ref2m_u      => clm3%g%l%c%p%pes%t_ref2m_u
    t_ref2m_r      => clm3%g%l%c%p%pes%t_ref2m_r
    plandunit      => clm3%g%l%c%p%landunit
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out  
    !-----------------------
    ! added by Jing Chen, 17 May 2012
    mass_rain_over_EASS   => clm3%g%l%c%p%pws%mass_rain_over_EASS
    mass_rain_under_EASS  => clm3%g%l%c%p%pws%mass_rain_under_EASS
    mass_snow_over_EASS   => clm3%g%l%c%p%pws%mass_snow_over_EASS
    mass_snow_under_EASS  => clm3%g%l%c%p%pws%mass_snow_under_EASS
    frac_rain_over_EASS   => clm3%g%l%c%p%pps%frac_rain_over_EASS
    frac_rain_under_EASS  => clm3%g%l%c%p%pps%frac_rain_under_EASS
    frac_snow_over_EASS   => clm3%g%l%c%p%pps%frac_snow_over_EASS
    frac_snow_under_EASS  => clm3%g%l%c%p%pps%frac_snow_under_EASS
    sum_rain_over_EASS    => clm3%g%l%c%p%pwf%sum_rain_over_EASS
    sum_rain_under_EASS   => clm3%g%l%c%p%pwf%sum_rain_under_EASS
    sum_snow_over_EASS    => clm3%g%l%c%p%pwf%sum_snow_over_EASS
    sum_snow_under_EASS   => clm3%g%l%c%p%pwf%sum_snow_under_EASS
    sum_snow_EASS         => clm3%g%l%c%p%pwf%sum_snow_EASS
    albedo_snow_EASS      => pftcon%albedo_snow_EASS
    ivt                   => clm3%g%l%c%p%itype
    Tleaf_over_sunlit_EASS    => clm3%g%l%c%p%pes%Tleaf_over_sunlit_EASS
    Tleaf_under_sunlit_EASS   => clm3%g%l%c%p%pes%Tleaf_under_sunlit_EASS
    Tleaf_over_shaded_EASS    => clm3%g%l%c%p%pes%Tleaf_over_shaded_EASS
    Tleaf_under_shaded_EASS   => clm3%g%l%c%p%pes%Tleaf_under_shaded_EASS
    eflx_sh_over_EASS         => clm3%g%l%c%p%pef%eflx_sh_over_EASS
    rs_over_shaded_EASS       => clm3%g%l%c%p%pps%rs_over_shaded_EASS
    rs_over_sunlit_EASS       => clm3%g%l%c%p%pps%rs_over_sunlit_EASS
    rs_under_shaded_EASS      => clm3%g%l%c%p%pps%rs_under_shaded_EASS
    rs_under_sunlit_EASS      => clm3%g%l%c%p%pps%rs_under_sunlit_EASS
    ci_over_shaded_EASS       => clm3%g%l%c%p%pps%ci_over_shaded_EASS
    ci_over_sunlit_EASS       => clm3%g%l%c%p%pps%ci_over_sunlit_EASS
    ci_under_shaded_EASS      => clm3%g%l%c%p%pps%ci_under_shaded_EASS
    ci_under_sunlit_EASS      => clm3%g%l%c%p%pps%ci_under_sunlit_EASS
    cs_over_shaded_EASS       => clm3%g%l%c%p%pps%cs_over_shaded_EASS
    cs_over_sunlit_EASS       => clm3%g%l%c%p%pps%cs_over_sunlit_EASS
    cs_under_shaded_EASS      => clm3%g%l%c%p%pps%cs_under_shaded_EASS
    cs_under_sunlit_EASS      => clm3%g%l%c%p%pps%cs_under_sunlit_EASS
    pgridcell                 => clm3%g%l%c%p%gridcell
    h2ocan_pft_EASS           => clm3%g%l%c%p%pws%h2ocan_EASS
    !-----------------------

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! NOTE: h2ocan, h2osno, and snowdp has valid values everywhere
    ! canopy water (pft level)

    do p = begp, endp
       h2ocan_pft(p) = 0._r8

       !-----------------------
       ! added by Jing Chen, 17 May 2012
       h2ocan_pft_EASS(p)      = 0._r8

       mass_rain_over_EASS(p)  = 0.0_r8
       mass_rain_under_EASS(p) = 0.0_r8
       mass_snow_over_EASS(p)  = 0.0_r8  
       mass_snow_under_EASS(p) = 0.0_r8

       frac_snow_over_EASS(p)  = 0.0_r8
       frac_snow_under_EASS(p) = 0.0_r8
       frac_rain_over_EASS(p)  = 0.0_r8
       frac_rain_under_EASS(p) = 0.0_r8

       sum_snow_over_EASS(p)   = 0.0_r8
       sum_snow_under_EASS(p)  = 0.0_r8
       sum_rain_over_EASS(p)   = 0.0_r8
       sum_rain_under_EASS(p)  = 0.0_r8

       sum_snow_EASS(p)        = 0.0_r8
       
       !-----------------------
       ! added for canopy water mass balance under dynamic pft weights
       !clm3%g%l%c%p%pps%tlai(p) = 0._r8
       !clm3%g%l%c%p%pps%tsai(p) = 0._r8
       !clm3%g%l%c%p%pps%elai(p) = 0._r8
       !clm3%g%l%c%p%pps%esai(p) = 0._r8
       !clm3%g%l%c%p%pps%htop(p) = 0._r8
       !clm3%g%l%c%p%pps%hbot(p) = 0._r8
       !clm3%g%l%c%p%pps%frac_veg_nosno_alb(p) = 0._r8
    end do

    do c = begc,endc

       ! canopy water (column level)

       h2ocan_col(c) = 0._r8
       !-------------
       !added by Jing Chen Oct 19 2012
       h2ocan_col_EASS(c) = 0._r8
       !-------------

       ! snow water

       l = clandunit(c)

       ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
       ! This gives more realistic values of qflx_glcice sooner in the simulation
       !  for columns with net ablation, at the cost of delaying ice formation
       !  in columns with net accumulation.
       if (ltype(l)==istice) then
          h2osno(c) = h2osno_max
       elseif (ltype(l)==istice_mec) then
          h2osno(c) = 0.5_r8 * h2osno_max   ! 50 cm if h2osno_max = 1 m
       else
          h2osno(c) = 0._r8
       endif
       !---------------
       ! added by Jing Chen Oct 17 2012
       h2osno_EASS(c)=h2osno(c)
       !---------------
       
       ! snow depth

       snowdp(c)  = h2osno(c) / bdsno

       !-----------------
       !added by Jing Chen, 17 May 2012
       mass_snow_grnd_EASS(c)=h2osno(c)
       densnow_EASS(c)=100.0_r8 * 2.0_r8    ![kg m-3]
       snowdp_EASS(c) =snowdp(c)
       ponddp_EASS(c) = 0.0_r8
       !-----------------

       ! Initialize Irrigation to zero
       if (ltype(l)==istsoil) then
          n_irrig_steps_left(c) = 0
          irrig_rate(c)         = 0.0_r8
       end if

    end do

    ! Set snow layer number, depth and thickiness

    call snowdp2lev(begc, endc)

    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

    ! NOTE: THESE MEMORY COPIES ARE INEFFICIENT -- SINCE nlev LOOP IS NESTED FIRST!!!!
    do c = begc,endc

       t_soisno(c,-nlevsno+1:nlevgrnd) = spval
       t_lake(c,1:nlevlak) = spval

       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          t_soisno(c,-nlevsno+1:0) = spval
          if (snl(c) < 0) then    !snow layer temperatures
             do j = snl(c)+1, 0
                t_soisno(c,j) = 250._r8
             enddo
          endif
          if (ltype(l)==istice .or. ltype(l)==istice_mec) then
             do j = 1, nlevgrnd
                t_soisno(c,j) = 250._r8
             end do
          else if (ltype(l) == istwet) then
             do j = 1, nlevgrnd
                t_soisno(c,j) = 277._r8
             end do
          else if (ltype(l) == isturb) then
#if (defined VANCOUVER)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
               ! Set road top layer to initial air temperature and interpolate other
               ! layers down to 20C in bottom layer
               do j = 1, nlevurb
                  t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevurb-1)) 
               end do
             ! Set wall and roof layers to initial air temperature
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
               do j = 1, nlevurb
                  t_soisno(c,j) = 297.56
               end do
             else
               do j = 1, nlevurb
                  t_soisno(c,j) = 283._r8
               end do
             end if
#elif (defined MEXICOCITY)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
               ! Set road top layer to initial air temperature and interpolate other
               ! layers down to 22C in bottom layer
               do j = 1, nlevurb
                  t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevurb-1)) 
               end do
             ! Set wall and roof layers to initial air temperature
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
               do j = 1, nlevurb
                  t_soisno(c,j) = 289.46
               end do
             else
               do j = 1, nlevurb
                  t_soisno(c,j) = 283._r8
               end do
             end if
#else
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
               do j = 1, nlevurb
                 t_soisno(c,j) = 274._r8
               end do
             ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
             ! shock from large heating/air conditioning flux
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                      .or. ctype(c) == icol_roof) then
               do j = 1, nlevurb
                 t_soisno(c,j) = 292._r8
               end do
             end if
#endif
          else
             do j = 1, nlevgrnd
                t_soisno(c,j) = 274._r8
             end do
          endif
          t_grnd(c) = t_soisno(c,snl(c)+1)
       else                     !lake
          t_lake(c,1:nlevlak) = 277._r8
          t_grnd(c) = t_lake(c,1)
       endif
       !-------------------
       ! added by Jing Chen, Oct 18 2012
       do j = 1, nlevgrnd 
          t_soisno_EASS(c,j)=t_soisno(c,j)
       end do
       snl_EASS(c)   = snl(c)
       t_grnd_EASS(c)= t_grnd(c)
       !-------------------
    end do

    call perturbIC( clm3%g%l )

    do p = begp, endp
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)

       !-----------------------
       !added by Jing Chen, 17 May 2012
       albsnow_EASS(c,1)=albedo_snow_EASS(ivt(p))   !parameter1[12]
       albsnow_EASS(c,2)=albedo_snow_EASS(ivt(p))
       !--------------------------------------------

       ! Initialize Irrigation to zero
       if (ltype(l)==istsoil) then
          qflx_irrig(c)      = 0.0_r8
       end if

#if (defined VANCOUVER)
       t_veg(p) = 297.56
       t_ref2m(p) = 297.56
       if (urbpoi(l)) then
         t_ref2m_u(p) = 297.56
       else
         t_ref2m_u(p) = spval
       end if
       if (ifspecial(l)) then
         t_ref2m_r(p) = spval
       else
         t_ref2m_r(p) = 297.56
       end if 
#elif (defined MEXICOCITY)
       t_veg(p) = 289.46
       t_ref2m(p) = 289.46
       if (urbpoi(l)) then
         t_ref2m_u(p) = 289.46
       else
         t_ref2m_u(p) = spval
       end if
       if (ifspecial(l)) then
         t_ref2m_r(p) = spval
       else
         t_ref2m_r(p) = 289.46
       end if 
#else
       t_veg(p) = 283._r8
       t_ref2m(p) = 283._r8
       if (urbpoi(l)) then
         t_ref2m_u(p) = 283._r8
       else
         t_ref2m_u(p) = spval
       end if
       if (ifspecial(l)) then
         t_ref2m_r(p) = spval
       else
         t_ref2m_r(p) = 283._r8
       end if
#endif


       !-----------------------
       !added by Jing Chen, 20 May 2012

       Tleaf_over_sunlit_EASS(p) = t_veg(p)!clm_a2l%forc_t(g)-0.5_r8
       Tleaf_under_sunlit_EASS(p)= t_veg(p)!clm_a2l%forc_t(g)-0.5_r8
       Tleaf_over_shaded_EASS(p) = t_veg(p)!clm_a2l%forc_t(g)-0.5_r8
       Tleaf_under_shaded_EASS(p)= t_veg(p)!clm_a2l%forc_t(g)-0.5_r8

       ! sensible flux of overstory, used in AirResist_EASS, updated in 6.2
       eflx_sh_over_EASS(p)=0.0_r8

       ! stomatal resistant of leaves
       rs_over_sunlit_EASS(p)=200.0_r8
       rs_over_shaded_EASS(p)=200.0_r8
       rs_under_sunlit_EASS(p)=300.0_r8
       rs_under_shaded_EASS(p)=300.0_r8

       ! intercellular leaf co2 concentration
   !#    ci_over_sunlit_EASS(p)=0.7_r8 * forc_pco2(g)
   !#    ci_over_shaded_EASS(p)=0.7_r8 * forc_pco2(g)
   !#    ci_under_sunlit_EASS(p)=0.7_r8 * forc_pco2(g)
   !#    ci_under_shaded_EASS(p)=0.7_r8 * forc_pco2(g)

       !co2 concentration on the surfaces of leaves
  !#     cs_over_sunlit_EASS(p)=forc_pco2(g)
  !#     cs_over_shaded_EASS(p)=forc_pco2(g)
  !#     cs_under_sunlit_EASS(p)=forc_pco2(g)
  !#     cs_under_shaded_EASS(p)=forc_pco2(g)


       !-----------------------
       eflx_lwrad_out(p) = sb * (t_grnd(c))**4

    end do

    ! Set snow/soil ice and liquid mass

    ! volumetric water is set first and liquid content and ice lens are obtained
    ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
    ! and urban pervious road (other urban columns have zero soil water)

    h2osoi_vol(begc:endc,         1:) = spval
    h2osoi_liq(begc:endc,-nlevsno+1:) = spval
    h2osoi_ice(begc:endc,-nlevsno+1:) = spval
   !-----------------
   ! added by Jing Chen Oct 17 2012
    h2osoi_vol_EASS(begc:endc,         1:) = spval
    h2osoi_liq_EASS(begc:endc,-nlevsno+1:) = spval
    h2osoi_ice_EASS(begc:endc,-nlevsno+1:) = spval
   !----------------

    wa(begc:endc)  = 5000._r8
    wt(begc:endc)  = 5000._r8
    zwt(begc:endc) = 0._r8

    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_road_perv) then
                wa(c)  = 4800._r8
                wt(c)  = wa(c)
                zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
             else
                wa(c)  = spval
                wt(c)  = spval
                zwt(c) = spval
             end if
          else
             wa(c)  = 4800._r8
             wt(c)  = wa(c)
             zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
          end if
       end if
       !------------
       !added by Jing Chen Oct 19 2012
       wa_EASS(c) = wa(c)
       wt_EASS(c) = wt(c)
       zwt_EASS(c) = zwt(c)
       !-----------
    end do

    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake

          ! volumetric water
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   h2osoi_vol(c,j) = 0.0_r8
                else
                   h2osoi_vol(c,j) = 0.3_r8
                endif
             end do
          else if (ltype(l) == isturb) then 
             nlevs = nlevurb
             do j = 1, nlevs
                if (ctype(c) == icol_road_perv .and. j <= nlevsoi) then
                   h2osoi_vol(c,j) = 0.3_r8
                else
                   h2osoi_vol(c,j) = 0.0_r8
                end if
             end do
          else if (ltype(l) == istwet) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   h2osoi_vol(c,j) = 0.0_r8
                else
                   h2osoi_vol(c,j) = 1.0_r8
                endif
             end do
          else if (ltype(l) == istice .or. ltype(l) == istice_mec) then
             nlevs = nlevgrnd 
             do j = 1, nlevs
                h2osoi_vol(c,j) = 1.0_r8
             end do
          endif
          do j = 1, nlevs
             h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
        
             ! soil layers
             if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
                h2osoi_liq(c,j) = 0._r8
             else
                h2osoi_ice(c,j) = 0._r8
                h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
             endif
             !------------
             ! added by Jing Chen Oct 17 2012
             h2osoi_vol_EASS(c,j) = h2osoi_vol(c,j)
             h2osoi_ice_EASS(c,j) = h2osoi_ice(c,j)
             h2osoi_liq_EASS(c,j) = h2osoi_liq(c,j)
             !---------------
          end do

#if (defined CN) || (defined CASA)
          ! soil water potential (added 10/21/03, PET)
          ! required for CN and CASA code
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (h2osoi_liq(c,j) > 0._r8) then
                   vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
                   psi = psisat(c,j) * (vwc/vwcsat(c,j))**bsw2(c,j)
                   soilpsi(c,j) = max(psi, -15.0_r8)
                   soilpsi(c,j) = min(soilpsi(c,j),0.0_r8)
                end if
             end do
          end if
#endif
       end if

    end do

    ! Set snow

    do j = -nlevsno+1, 0
       do c = begc,endc
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake
             if (j > snl(c)) then
                h2osoi_ice(c,j) = dz(c,j)*250._r8
                h2osoi_liq(c,j) = 0._r8
                !------------
                ! added by Jing Chen Oct 17 2012
                h2osoi_ice_EASS(c,j) = h2osoi_ice(c,j)
                h2osoi_liq_EASS(c,j) = h2osoi_liq(c,j)
                !---------------
             end if
          end if
       end do
    end do


    ! initialize SNICAR fields:
    do c = begc,endc
       mss_bctot(c,:) = 0._r8
       mss_bcpho(c,:) = 0._r8
       mss_bcphi(c,:) = 0._r8
       mss_cnc_bcphi(c,:)=0._r8
       mss_cnc_bcpho(c,:)=0._r8

       mss_octot(c,:) = 0._r8
       mss_ocpho(c,:) = 0._r8
       mss_ocphi(c,:) = 0._r8
       mss_cnc_ocphi(c,:)=0._r8
       mss_cnc_ocpho(c,:)=0._r8
       
       mss_dst1(c,:) = 0._r8
       mss_dst2(c,:) = 0._r8
       mss_dst3(c,:) = 0._r8
       mss_dst4(c,:) = 0._r8
       mss_dsttot(c,:) = 0._r8
       mss_cnc_dst1(c,:)=0._r8
       mss_cnc_dst2(c,:)=0._r8
       mss_cnc_dst3(c,:)=0._r8
       mss_cnc_dst4(c,:)=0._r8
       
       if (snl(c) < 0) then
          snw_rds(c,snl(c)+1:0)        = snw_rds_min
          snw_rds(c,-nlevsno+1:snl(c)) = 0._r8
          snw_rds_top(c)               = snw_rds_min
          sno_liq_top(c) = h2osoi_liq(c,snl(c)+1) / (h2osoi_liq(c,snl(c)+1)+h2osoi_ice(c,snl(c)+1))
       elseif (h2osno(c) > 0._r8) then
          snw_rds(c,0)             = snw_rds_min
          snw_rds(c,-nlevsno+1:-1) = 0._r8
          snw_rds_top(c)           = spval
          sno_liq_top(c)           = spval
       else
          snw_rds(c,:)   = 0._r8
          snw_rds_top(c) = spval
          sno_liq_top(c) = spval
       endif
    enddo

  end subroutine mkarbinit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: perturbIC
!
! !INTERFACE:
  subroutine perturbIC( landunit )
!
! !DESCRIPTION:
!   Perturbs initial conditions by the amount in the namelist variable pertlim.
!
! !USES:
    use clm_varpar   , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varctl   , only : pertlim
    use clm_varcon   , only : isturb
    use decompMod    , only : get_proc_bounds
    use clmtype      , only : landunit_type
    implicit none
!
! !ARGUMENTS:
    type(landunit_type), intent(INOUT) :: landunit
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
! !LOCAL VARIABLES:
!
    integer :: j,l,c        ! indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    real(r8):: pertval      ! for calculating temperature perturbation
    integer :: nlevs        ! number of levels
    integer , pointer :: clandunit(:)   ! landunit index associated with each column
    integer , pointer :: ltype(:)       ! landunit type
    logical , pointer :: lakpoi(:)      ! true => landunit is a lake point
    integer , pointer :: snl(:)         ! number of snow layers
    real(r8), pointer :: t_soisno(:,:)  ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_lake(:,:)    ! lake temperature (Kelvin)  (1:nlevlak)
    real(r8), pointer :: t_grnd(:)      ! ground temperature (Kelvin)
    !---------------
    ! added by Jing Chen Oct 17 2012
    real(r8), pointer :: t_soisno_EASS(:,:)  ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_grnd_EASS(:)      ! ground temperature (Kelvin)
    !--------------------

!EOP
!-----------------------------------------------------------------------

    if ( pertlim /= 0.0_r8 )then

       if ( masterproc ) write(iulog,*) 'Applying perturbation to initial soil temperature'

       clandunit  => landunit%c%landunit
       lakpoi     => landunit%lakpoi
       ltype      => landunit%itype
       t_soisno   => landunit%c%ces%t_soisno
       t_lake     => landunit%c%ces%t_lake
       t_grnd     => landunit%c%ces%t_grnd
       snl        => landunit%c%cps%snl
       !-------------
       ! added by Jing Chen Oct 17 2012
       t_soisno_EASS   => landunit%c%ces%t_soisno_EASS
       t_grnd_EASS     => landunit%c%ces%t_grnd_EASS
       !-------------------

       ! Determine subgrid bounds on this processor

       call get_proc_bounds( begc=begc, endc=endc )

       do c = begc,endc
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake
             if ( ltype(l) == isturb) then
                nlevs = nlevurb
             else
                nlevs = nlevgrnd
             end if
             
             ! Randomly perturb soil temperature
             do j = 1, nlevs
                call random_number (pertval)
                pertval       = 2._r8*pertlim*(0.5_r8 - pertval)
                t_soisno(c,j) = t_soisno(c,j)*(1._r8 + pertval)
                !-----------
                ! added by Jing Chen Oct 19 2012
                t_soisno_EASS(c,j) = t_soisno_EASS(c,j)*(1._r8 + pertval)
                !-----------------
             end do
          else                       !lake
             ! Randomly perturb lake temperature
             do j = 1, nlevlak
                call random_number (pertval)
                pertval     = 2._r8*pertlim*(0.5_r8 - pertval)
                t_lake(c,j) = t_lake(c,j)*(1._r8 + pertval)
             end do
          endif
          ! Randomly perturb surface ground temp
          call random_number (pertval)
          pertval   = 2._r8*pertlim*(0.5_r8 - pertval)
          t_grnd(c) = t_grnd(c)*(1._r8 + pertval)
          !-------------
          ! added by Jing Chen Oct 19 2012
          t_grnd_EASS(c) = t_grnd_EASS(c)*(1._r8 + pertval)
          !--------------
       end do
    end if
    !-----------------------------------------------------------------------

  end subroutine perturbIC

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: snowdp2lev
!
! !INTERFACE:
  subroutine snowdp2lev(lbc, ubc)
!
! !DESCRIPTION:
! Create snow layers and interfaces given snow depth.
! Note that cps%zi(0) is set in routine iniTimeConst.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use clm_varpar  , only : nlevsno
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: lbc, ubc                    ! column bounds
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:)  ! landunit index associated with each column
  real(r8), pointer :: snowdp(:)     ! snow height (m)
  logical , pointer :: lakpoi(:)     ! true => landunit is a lake point
!
! local pointers to implicit out arguments
!
  integer , pointer :: snl(:)        ! number of snow layers
  real(r8), pointer :: z(:,:)        ! layer depth  (m) over snow only
  real(r8), pointer :: dz(:,:)       ! layer thickness depth (m) over snow only
  real(r8), pointer :: zi(:,:)       ! interface depth (m) over snow only
  !----------------------
  ! added by Jing Chen Oct 19 2012
  integer , pointer :: snl_EASS(:)        ! number of snow layers
  real(r8), pointer :: z_EASS(:,:)        ! layer depth  (m) over snow only
  real(r8), pointer :: dz_EASS(:,:)       ! layer thickness depth (m) over snow only
  real(r8), pointer :: zi_EASS(:,:)       ! interface depth (m) over snow only
  !--------------------
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: c,l,j      !indices
!-----------------------------------------------------------------------
 
  ! Assign local pointers to derived subtypes components (landunit-level)

  lakpoi => clm3%g%l%lakpoi

  ! Assign local pointers to derived type members (column-level)

  clandunit => clm3%g%l%c%landunit
  snowdp    => clm3%g%l%c%cps%snowdp
  snl       => clm3%g%l%c%cps%snl
  zi        => clm3%g%l%c%cps%zi
  dz        => clm3%g%l%c%cps%dz
  z         => clm3%g%l%c%cps%z
  !------------------
  ! added by Jing Chen Oct 19 2012
  snl_EASS  => clm3%g%l%c%cps%snl_EASS
  zi_EASS   => clm3%g%l%c%cps%zi_EASS
  dz_EASS   => clm3%g%l%c%cps%dz_EASS
  z_EASS    => clm3%g%l%c%cps%z_EASS
  !--------------------

  ! Initialize snow levels and interfaces (lake and non-lake points)

  do c = lbc, ubc
     dz(c,-nlevsno+1: 0) = 1.e36_r8
     z (c,-nlevsno+1: 0) = 1.e36_r8
     zi(c,-nlevsno  :-1) = 1.e36_r8
  end do

  ! Determine snow levels and interfaces for non-lake points

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        if (snowdp(c) < 0.01_r8) then
           snl(c) = 0
           dz(c,-nlevsno+1:0) = 0._r8
           z (c,-nlevsno+1:0) = 0._r8
           zi(c,-nlevsno+0:0) = 0._r8
        else
           if ((snowdp(c) >= 0.01_r8) .and. (snowdp(c) <= 0.03_r8)) then
              snl(c) = -1
              dz(c,0)  = snowdp(c)
           else if ((snowdp(c) > 0.03_r8) .and. (snowdp(c) <= 0.04_r8)) then
              snl(c) = -2
              dz(c,-1) = snowdp(c)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.04_r8) .and. (snowdp(c) <= 0.07_r8)) then
              snl(c) = -2
              dz(c,-1) = 0.02_r8
              dz(c, 0) = snowdp(c) - dz(c,-1)
           else if ((snowdp(c) > 0.07_r8) .and. (snowdp(c) <= 0.12_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = (snowdp(c) - 0.02_r8)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.12_r8) .and. (snowdp(c) <= 0.18_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = 0.05_r8
              dz(c, 0) = snowdp(c) - dz(c,-2) - dz(c,-1)
           else if ((snowdp(c) > 0.18_r8) .and. (snowdp(c) <= 0.29_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = (snowdp(c) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.29_r8) .and. (snowdp(c) <= 0.41_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = 0.11_r8
              dz(c, 0) = snowdp(c) - dz(c,-3) - dz(c,-2) - dz(c,-1)
           else if ((snowdp(c) > 0.41_r8) .and. (snowdp(c) <= 0.64_r8)) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = (snowdp(c) - dz(c,-4) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if (snowdp(c) > 0.64_r8) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = 0.23_r8
              dz(c, 0)=snowdp(c)-dz(c,-4)-dz(c,-3)-dz(c,-2)-dz(c,-1)
           endif
        end if
     end if
  end do

  ! The following loop is currently not vectorized

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        do j = 0, snl(c)+1, -1
           z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
           zi(c,j-1) = zi(c,j) - dz(c,j)
           !-----------------
           ! added by Jing Chen Oct 19 2012
           dz_EASS(c,j) = dz(c,j)
           z_EASS(c,j)  = z(c,j)
           zi_EASS(c,j) = zi(c,j)
           zi_EASS(c,j-1)=zi(c,j-1)
           !------------------ 
        end do
     end if
  end do

  ! Determine snow levels and interfaces for lake points

  do c = lbc,ubc
     l = clandunit(c)
     if (lakpoi(l)) then
        snl(c) = 0
        dz(c,-nlevsno+1:0) = 0._r8
        z (c,-nlevsno+1:0) = 0._r8
        zi(c,-nlevsno+0:0) = 0._r8
        !-----------------
        ! added by Jing Chen Oct 19 2012
        snl_EASS(c) = 0
        dz_EASS(c,-nlevsno+1:0) = dz(c,-nlevsno+1:0)
        z_EASS(c,-nlevsno+1:0)  = z(c,-nlevsno+1:0)
        zi_EASS(c,-nlevsno+1:0) = zi(c,-nlevsno+1:0)
        !------------------
     end if
  end do

  end subroutine snowdp2lev


!-----------------------------------------------------------------------

end module mkarbinitMod
