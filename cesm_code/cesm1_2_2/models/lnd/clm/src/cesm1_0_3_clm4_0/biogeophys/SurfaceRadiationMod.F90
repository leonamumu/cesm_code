module SurfaceRadiationMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceRadiationMod
!
! !DESCRIPTION:
! Calculate solar fluxes absorbed by vegetation and ground surface
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varctl  , only: iulog
   use clm_varpar  , only: numrad

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface
  !------------------------------------------
  !added by Jing Chen, 17 May 2012
  public :: SunlitLAI_EASS   ! LAI of sunlit leaves based on clumping index in EASS model

!! PUBLIC DATA MEMBERS:
  real(r8), public :: albedo_new_snow_EASS(numrad)=(/0.94_r8,0.8_r8/)      !albedo of snow by waveband(1=vis, 2=nir)
  real(r8), public :: albsnow_const_EASS(numrad)=(/0.70_r8,0.42_r8/)      !albedo of snow by waveband(1=vis, 2=nir)
  !-------------------------------------------

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 11/26/03, Peter Thornton: Added new routine for improved treatment of
!    sunlit/shaded canopy radiation.
! 4/26/05, Peter Thornton: Adopted the sun/shade algorithm as the default,
!    removed the old SurfaceRadiation(), and renamed SurfaceRadiationSunShade()
!    as SurfaceRadiation().
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRadiation
!
! !INTERFACE:
   subroutine SurfaceRadiation(lbp, ubp, num_nourbanp, filter_nourbanp)
!
! !DESCRIPTION: 
! Solar fluxes absorbed by vegetation and ground surface
! Note possible problem when land is on different grid than atmosphere.
! Land may have sun above the horizon (coszen > 0) but atmosphere may
! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
! because all fluxes (absorbed, reflected, transmitted) are multiplied
! by the incoming flux and all will equal zero.
! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
! land may have sun below horizon. This is okay because fabd, fabi,
! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
! the radiation is reflected. NDVI should equal zero in this case.
! However, the way the code is currently implemented this is only true
! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
!
! !USES:
     use clmtype
     use clm_atmlnd      , only : clm_a2l
     use clm_varpar      , only : numrad
     use clm_varcon      , only : spval, istsoil, degpsec, isecspday
     use clm_varcon      , only : istice_mec
     use clm_varcon      , only : istcrop
     use clm_time_manager, only : get_curr_date, get_step_size
     use clm_varpar      , only : nlevsno
     use SNICARMod       , only : DO_SNO_OC
     use abortutils      , only : endrun
!---------------------------------------------------------------------
     !added by Jing Chen, 26 March 2012
     use pftvarcon,   only : nc4_grass, ncorn
!#     use clm_varcon,  only : istcrop, istsoil
!--------------------------------------------------------------------

!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: lbp, ubp                   ! pft upper and lower bounds
     integer, intent(in) :: num_nourbanp               ! number of pfts in non-urban points in pft filter
     integer, intent(in) :: filter_nourbanp(ubp-lbp+1) ! pft filter for non-urban points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/18/02, Peter Thornton: Migrated to new data structures. Added a pft loop.
! 6/05/03, Peter Thornton: Modified sunlit/shaded canopy treatment. Original code
! had all radiation being absorbed in the sunlit canopy, and now the sunlit and shaded
! canopies are each given the appropriate fluxes.  There was also an inconsistency in
! the original code, where parsun was not being scaled by leaf area, and so represented
! the entire canopy flux.  This goes into Stomata (in CanopyFluxes) where it is assumed
! to be a flux per unit leaf area. In addition, the fpsn flux coming out of Stomata was
! being scaled back up to the canopy by multiplying by lai, but the input radiation flux was
! for the entire canopy to begin with.  Corrected this inconsistency in this version, so that
! the parsun and parsha fluxes going into canopy fluxes are per unit lai in the sunlit and
! shaded canopies.
! 6/9/03, Peter Thornton: Moved coszen from g%gps to c%cps to avoid problem
! with OpenMP threading over columns, where different columns hit the radiation
! time step at different times during execution.
! 6/10/03, Peter Thornton: Added constraint on negative tot_aid, instead of
! exiting with error. Appears to be happening only at roundoff level.
! 6/11/03, Peter Thornton: Moved calculation of ext inside if (coszen),
! and added check on laisun = 0 and laisha = 0 in calculation of sun_aperlai
! and sha_aperlai.
! 11/26/03, Peter Thornton: During migration to new vector code, created 
!   this as a new routine to handle sunlit/shaded canopy calculations.
! 03/28/08, Mark Flanner: Incorporated SNICAR, including absorbed solar radiation
!   in each snow layer and top soil layer, and optional radiative forcing calculation
! 26 March 2012: Jing Chen; EASS algorithm
! Added clumping index in lai calculation, and divid the values of lai into overstory(including sunlit and shaded leaves) 
! and overstory(including sunlit and shaded leaves), based on EASS_region_1028 
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
     integer , pointer :: ivt(:)           ! pft vegetation type
     integer , pointer :: pcolumn(:)       ! pft's column index
     integer , pointer :: pgridcell(:)     ! pft's gridcell index
     real(r8), pointer :: pwtgcell(:)      ! pft's weight relative to corresponding gridcell
     real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
     real(r8), pointer :: esai(:)          ! one-sided stem area index with burying by snow
     real(r8), pointer :: londeg(:)        ! longitude (degrees)
     real(r8), pointer :: latdeg(:)        ! latitude (degrees)
     real(r8), pointer :: slasun(:)        ! specific leaf area for sunlit canopy, projected area basis (m^2/gC)
     real(r8), pointer :: slasha(:)        ! specific leaf area for shaded canopy, projected area basis (m^2/gC)
     real(r8), pointer :: gdir(:)	   ! leaf projection in solar direction (0 to 1)
     real(r8), pointer :: omega(:,:)       ! fraction of intercepted radiation that is scattered (0 to 1)
     real(r8), pointer :: coszen(:)	   ! cosine of solar zenith angle
     real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (W/m**2)
     real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation (W/m**2)
     real(r8), pointer :: fabd(:,:)        ! flux absorbed by veg per unit direct flux
     real(r8), pointer :: fabi(:,:)        ! flux absorbed by veg per unit diffuse flux
     real(r8), pointer :: ftdd(:,:)        ! down direct flux below veg per unit dir flx
     real(r8), pointer :: ftid(:,:)        ! down diffuse flux below veg per unit dir flx
     real(r8), pointer :: ftii(:,:)        ! down diffuse flux below veg per unit dif flx
     real(r8), pointer :: albgrd(:,:)      ! ground albedo (direct)
     real(r8), pointer :: albgri(:,:)      ! ground albedo (diffuse)
     real(r8), pointer :: albd(:,:)        ! surface albedo (direct)
     real(r8), pointer :: albi(:,:)        ! surface albedo (diffuse)
     real(r8), pointer :: slatop(:)        ! specific leaf area at top of canopy, projected area basis [m^2/gC]
     real(r8), pointer :: dsladlai(:)      ! dSLA/dLAI, projected area basis [m^2/gC]
!-------------------------------------------------------
     !added by Jing Chen, 26 March 2012
     real(r8), pointer :: clumping_EASS(:)   ! clumping index in EASS(input)
     real(r8), pointer :: elai_over_EASS(:)  ! elai of overstory in EASS(input)
     real(r8), pointer :: elai_under_EASS(:) ! elai of understory in EASS(input)
     real(r8), pointer :: esai_over_EASS(:)  ! esai of overstory canopy, EASS(input)
     real(r8), pointer :: esai_under_EASS(:) ! esai of understory canopy, EASS(input)
     real(r8), pointer :: elai_EASS(:)       ! one-sided leaf area index with burying by snow, corrected by clumping index, in EASS
     real(r8), pointer :: lai_over_max_EASS(:)  ! ecophys const - max lai of overstory in EASS
     real(r8), pointer :: lai_under_max_EASS(:) ! ecophys const - max lai of understory in EASS
     real(r8), pointer :: clumping_pft_EASS(:)  ! ecophys const - clumping index
    real(r8), pointer :: albsod(:,:)               ! soil albedo (direct)
    real(r8), pointer :: forc_solar(:)             ! incident solar radiation [W m-2]
    real(r8), pointer :: forc_snow(:)              ! snow rate [mm s-1]
    real(r8), pointer :: albedo_nir_EASS(:)        ! albedo of canopy without snow covered based on pft, nir band 
    real(r8), pointer :: albedo_vis_EASS(:)        ! albedo of canopy without snow covered based on pft, vis band 
    real(r8), pointer :: sum_snow_over_EASS(:)     ! area of overstory covered by snow[m2 m-2]
    real(r8), pointer :: sum_snow_under_EASS(:)    ! area of the understory covered by snow[m2 s-2]
    real(r8), pointer :: frac_snow_grnd_EASS(:)    ! fraction of snow on the ground(0~1)
    real(r8), pointer :: albedo_snow_EASS(:)       ! albedo of new snow based on pft
    real(r8), pointer :: rhol(:,:)                 ! leaf reflectance: 1=vis, 2=nir

   ! local pointers to original implicit inout arguments
    real(r8), pointer :: sum_snow_EASS(:)        ! accumulation of snow mass base on forc_snow[m s-1]        
    real(r8), pointer :: snowdp_EASS(:)          ! snow depth on the ground, and will be updated for new snow fall down[m]
    real(r8), pointer :: albsnow_EASS(:,:)       ! albedo of new snow

!-------------------------------------------------------
   !
   ! local pointers to original implicit out arguments
   !
    real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
    real(r8), pointer :: laisun(:)        ! sunlit leaf area
    real(r8), pointer :: laisha(:)        ! shaded leaf area
    real(r8), pointer :: sabg(:)          ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: sabv(:)          ! solar radiation absorbed by vegetation (W/m**2)
    real(r8), pointer :: fsa(:)           ! solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsa_r(:)         ! rural solar radiation absorbed (total) (W/m**2)
    integer , pointer :: ityplun(:)       ! landunit type
    integer , pointer :: plandunit(:)     ! index into landunit level quantities
    real(r8), pointer :: parsun(:)        ! average absorbed PAR for sunlit leaves (W/m**2)
    real(r8), pointer :: parsha(:)        ! average absorbed PAR for shaded leaves (W/m**2)
    real(r8), pointer :: fsr(:)           ! solar radiation reflected (W/m**2)
    real(r8), pointer :: fsds_vis_d(:)    ! incident direct beam vis solar radiation (W/m**2)
    real(r8), pointer :: fsds_nir_d(:)    ! incident direct beam nir solar radiation (W/m**2)
    real(r8), pointer :: fsds_vis_i(:)    ! incident diffuse vis solar radiation (W/m**2)
    real(r8), pointer :: fsds_nir_i(:)    ! incident diffuse nir solar radiation (W/m**2)
    real(r8), pointer :: fsr_vis_d(:)     ! reflected direct beam vis solar radiation (W/m**2)
    real(r8), pointer :: fsr_nir_d(:)     ! reflected direct beam nir solar radiation (W/m**2)
    real(r8), pointer :: fsr_vis_i(:)     ! reflected diffuse vis solar radiation (W/m**2)
    real(r8), pointer :: fsr_nir_i(:)     ! reflected diffuse nir solar radiation (W/m**2)
    real(r8), pointer :: fsds_vis_d_ln(:) ! incident direct beam vis solar rad at local noon (W/m**2)
    real(r8), pointer :: fsds_nir_d_ln(:) ! incident direct beam nir solar rad at local noon (W/m**2)
    real(r8), pointer :: fsr_vis_d_ln(:)  ! reflected direct beam vis solar rad at local noon (W/m**2)
    real(r8), pointer :: fsr_nir_d_ln(:)  ! reflected direct beam nir solar rad at local noon (W/m**2)
    real(r8), pointer :: eff_kid(:,:)     ! effective extinction coefficient for indirect from direct
    real(r8), pointer :: eff_kii(:,:)     ! effective extinction coefficient for indirect from indirect
    real(r8), pointer :: sun_faid(:,:)    ! fraction sun canopy absorbed indirect from direct
    real(r8), pointer :: sun_faii(:,:)    ! fraction sun canopy absorbed indirect from indirect
    real(r8), pointer :: sha_faid(:,:)    ! fraction shade canopy absorbed indirect from direct
    real(r8), pointer :: sha_faii(:,:)    ! fraction shade canopy absorbed indirect from indirect
    real(r8), pointer :: sun_add(:,:)     ! sun canopy absorbed direct from direct (W/m**2)
    real(r8), pointer :: tot_aid(:,:)     ! total canopy absorbed indirect from direct (W/m**2)
    real(r8), pointer :: sun_aid(:,:)     ! sun canopy absorbed indirect from direct (W/m**2)
    real(r8), pointer :: sun_aii(:,:)     ! sun canopy absorbed indirect from indirect (W/m**2)
    real(r8), pointer :: sha_aid(:,:)     ! shade canopy absorbed indirect from direct (W/m**2)
    real(r8), pointer :: sha_aii(:,:)     ! shade canopy absorbed indirect from indirect (W/m**2)
    real(r8), pointer :: sun_atot(:,:)    ! sun canopy total absorbed (W/m**2)
    real(r8), pointer :: sha_atot(:,:)    ! shade canopy total absorbed (W/m**2)
    real(r8), pointer :: sun_alf(:,:)     ! sun canopy total absorbed by leaves (W/m**2)
    real(r8), pointer :: sha_alf(:,:)     ! shade canopy total absored by leaves (W/m**2)
    real(r8), pointer :: sun_aperlai(:,:) ! sun canopy total absorbed per unit LAI (W/m**2)
    real(r8), pointer :: sha_aperlai(:,:) ! shade canopy total absorbed per unit LAI (W/m**2)
    real(r8), pointer :: flx_absdv(:,:)   ! direct flux absorption factor (col,lyr): VIS [frc]
    real(r8), pointer :: flx_absdn(:,:)   ! direct flux absorption factor (col,lyr): NIR [frc]
    real(r8), pointer :: flx_absiv(:,:)   ! diffuse flux absorption factor (col,lyr): VIS [frc]
    real(r8), pointer :: flx_absin(:,:)   ! diffuse flux absorption factor (col,lyr): NIR [frc]
    integer , pointer :: snl(:)           ! negative number of snow layers [nbr]
    real(r8), pointer :: albgrd_pur(:,:)    ! pure snow ground albedo (direct)
    real(r8), pointer :: albgri_pur(:,:)    ! pure snow ground albedo (diffuse)
    real(r8), pointer :: albgrd_bc(:,:)     ! ground albedo without BC (direct) (col,bnd)
    real(r8), pointer :: albgri_bc(:,:)     ! ground albedo without BC (diffuse) (col,bnd)
    real(r8), pointer :: albgrd_oc(:,:)     ! ground albedo without OC (direct) (col,bnd)
    real(r8), pointer :: albgri_oc(:,:)     ! ground albedo without OC (diffuse) (col,bnd)
    real(r8), pointer :: albgrd_dst(:,:)    ! ground albedo without dust (direct) (col,bnd)
    real(r8), pointer :: albgri_dst(:,:)    ! ground albedo without dust (diffuse) (col,bnd)
    real(r8), pointer :: albsnd_hst(:,:)    ! snow albedo, direct, for history files (col,bnd) [frc]
    real(r8), pointer :: albsni_hst(:,:)    ! snow ground albedo, diffuse, for history files (col,bnd
    real(r8), pointer :: sabg_lyr(:,:)      ! absorbed radiative flux (pft,lyr) [W/m2]
    real(r8), pointer :: sfc_frc_aer(:)     ! surface forcing of snow with all aerosols (pft) [W/m2]
    real(r8), pointer :: sfc_frc_bc(:)      ! surface forcing of snow with BC (pft) [W/m2]
    real(r8), pointer :: sfc_frc_oc(:)      ! surface forcing of snow with OC (pft) [W/m2]
    real(r8), pointer :: sfc_frc_dst(:)     ! surface forcing of snow with dust (pft) [W/m2]
    real(r8), pointer :: sfc_frc_aer_sno(:) ! surface forcing of snow with all aerosols, averaged only when snow is present (pft) [W/m2]
    real(r8), pointer :: sfc_frc_bc_sno(:)  ! surface forcing of snow with BC, averaged only when snow is present (pft) [W/m2]
    real(r8), pointer :: sfc_frc_oc_sno(:)  ! surface forcing of snow with OC, averaged only when snow is present (pft) [W/m2]
    real(r8), pointer :: sfc_frc_dst_sno(:) ! surface forcing of snow with dust, averaged only when snow is present (pft) [W/m2]
    real(r8), pointer :: frac_sno(:)      ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: fsr_sno_vd(:)    ! reflected visible, direct radiation from snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsr_sno_nd(:)    ! reflected near-IR, direct radiation from snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsr_sno_vi(:)    ! reflected visible, diffuse radiation from snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsr_sno_ni(:)    ! reflected near-IR, diffuse radiation from snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsds_sno_vd(:)   ! incident visible, direct radiation on snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsds_sno_nd(:)   ! incident near-IR, direct radiation on snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsds_sno_vi(:)   ! incident visible, diffuse radiation on snow (for history files) (pft) [W/m2]
    real(r8), pointer :: fsds_sno_ni(:)   ! incident near-IR, diffuse radiation on snow (for history files) (pft) [W/m2]
    real(r8), pointer :: snowdp(:)        ! snow height (m)
    !--------------------------
    !added by Jing Chen, 18 May 2012
    real(r8), pointer :: elai_over_shaded_EASS(:) ! elai of shaded leaf in overstory canopy, EASS(output)
    real(r8), pointer :: elai_over_sunlit_EASS(:) ! elai of sunlit leaf in overstory canopy, EASS(output)
    real(r8), pointer :: elai_under_shaded_EASS(:)! elai of shaded leaf in understory canopy, EASS(output)
    real(r8), pointer :: elai_under_sunlit_EASS(:)! elai of sunlit leaf in understory canopy, EASS(output)
    real(r8), pointer :: lt_over_shaded_EASS(:)   ! (elai+esai) of shaded leaf in overstory, EASS(output)
    real(r8), pointer :: lt_over_sunlit_EASS(:)   ! (elai+esai) of sunlit leaf in overstory, EASS(output)
    real(r8), pointer :: lt_under_shaded_EASS(:)  ! (elai+esai) of shaded leaf in understory, EASS(output)
    real(r8), pointer :: lt_under_sunlit_EASS(:)  ! (elai+esai) of sunlit leaf in understory, EASS(output
    real(r8), pointer :: alb_can_EASS(:,:)        ! albedo of canopy without snow covered (1=vis,2=nir)
    real(r8), pointer :: albsno_over_EASS(:,:)    ! albedo of the overstroy covered by snow
    real(r8), pointer :: albsno_under_EASS(:,:)   ! albedo of the understory covered by snow
    real(r8), pointer :: alb_grnd_EASS(:,:)       ! albedo of ground(1=vis,2=nir)
    real(r8), pointer :: albsno_grnd_EASS(:,:)    ! albedo of the ground covered by snow
    real(r8), pointer :: dalb_over_EASS(:)        ! diff. between 1.0 and albedo of overstory
    real(r8), pointer :: dalb_under_EASS(:)       ! diff. between 1.0 and albedo of understory
    real(r8), pointer :: dalb_grnd_EASS(:)        ! diff. between 1.0 and albedo of ground
    real(r8), pointer :: swnet_EASS(:)                 ! net solar (shortwave) radiation [W m-2]
    !--------------------------
    !
    !
    ! !OTHER LOCAL VARIABLES:
    !EOP
    !
    integer , parameter :: nband = numrad    ! number of solar radiation waveband classes
    real(r8), parameter :: mpe = 1.e-06_r8   ! prevents overflow for division by zero
    integer  :: fp                  ! non-urban filter pft index
    integer  :: p                   ! pft index
    integer  :: c                   ! column index
    integer  :: l                   ! landunit index
    integer  :: g                   ! grid cell index
    integer  :: ib                  ! waveband number (1=vis, 2=nir)
    real(r8) :: absrad              ! absorbed solar radiation (W/m**2)
    real(r8) :: rnir                ! reflected solar radiation [nir] (W/m**2)
    real(r8) :: rvis                ! reflected solar radiation [vis] (W/m**2)
    real(r8) :: laifra              ! leaf area fraction of canopy
    real(r8) :: trd(lbp:ubp,numrad) ! transmitted solar radiation: direct (W/m**2)
    real(r8) :: tri(lbp:ubp,numrad) ! transmitted solar radiation: diffuse (W/m**2)
    real(r8) :: cad(lbp:ubp,numrad) ! direct beam absorbed by canopy (W/m**2)
    real(r8) :: cai(lbp:ubp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
    real(r8) :: vai(lbp:ubp)        ! total leaf area index + stem area index, one sided
    real(r8) :: ext                 ! optical depth direct beam per unit LAI+SAI
    real(r8) :: t1, t2              ! temporary variables
    real(r8) :: cosz
    integer  :: local_secp1         ! seconds into current date in local time
    real(r8) :: dtime               ! land model time step (sec)
    integer  :: year,month,day,secs !  calendar info for current time step
    integer  :: i                   ! layer index [idx]
    real(r8) :: sabg_snl_sum        ! temporary, absorbed energy in all active snow layers [W/m2]
    real(r8) :: absrad_pur          ! temp: absorbed solar radiation by pure snow [W/m2]
    real(r8) :: absrad_bc           ! temp: absorbed solar radiation without BC [W/m2]
    real(r8) :: absrad_oc           ! temp: absorbed solar radiation without OC [W/m2]
    real(r8) :: absrad_dst          ! temp: absorbed solar radiation without dust [W/m2]
    real(r8) :: sabg_pur(lbp:ubp)   ! solar radiation absorbed by ground with pure snow [W/m2]
    real(r8) :: sabg_bc(lbp:ubp)    ! solar radiation absorbed by ground without BC [W/m2]
    real(r8) :: sabg_oc(lbp:ubp)    ! solar radiation absorbed by ground without OC [W/m2]
    real(r8) :: sabg_dst(lbp:ubp)   ! solar radiation absorbed by ground without dust [W/m2]
    !-------------------------------------------------------
    !added by Jing Chen, 26 March 2012
    real(r8) :: lt_inter               !elai_over_EASS(:) + esai_over_EASS(:) + elai_under_EASS(:) + esai_under_EASS(:)
    real(r8) :: lt_EASS                
    real(r8) :: elai_can_inter         !elai_over_EASS(:) + elai_under_EASS(:)
    real(r8) :: elai_over_under_EASS           
    real(r8) :: lt_over_EASS           !(elai+esai) of overstory in EASS
    real(r8) :: lt_under_EASS          !(elai+esai) of understory in EASS
    character(len=20) :: phase         ! "LAI" or "energy" ,used in subroutine SunlitLAI_EASS
    real(r8) :: snow_over_EASS(lbp:ubp)!snow rate ,overstory [m s-1]
!------------------------------------------------------------------------------

    ! Assign local pointers to multi-level derived type members (gridcell level)

    londeg        => clm3%g%londeg
    latdeg        => clm3%g%latdeg
    forc_solad    => clm_a2l%forc_solad
    forc_solai    => clm_a2l%forc_solai
    !-------------------------------
    !added by Jing Chen, 18 May 2012
    forc_solar    => clm_a2l%forc_solar
    forc_snow     => clm_a2l%forc_snow
    !-------------------------------

    ! Assign lccal pointers to multi-level derived type members (landunit level)

    ityplun       => clm3%g%l%itype

    ! Assign local pointers to multi-level derived type members (column level)

    albgrd        => clm3%g%l%c%cps%albgrd
    albgri        => clm3%g%l%c%cps%albgri
    coszen        => clm3%g%l%c%cps%coszen
    !------------------------------
    !added by Jing Chen, 18 May 2012
    albsod              => clm3%g%l%c%cps%albsod
    frac_snow_grnd_EASS => clm3%g%l%c%cps%frac_snow_grnd_EASS
    snowdp_EASS         => clm3%g%l%c%cps%snowdp_EASS
    albsnow_EASS        => clm3%g%l%c%cps%albsnow_EASS
    alb_can_EASS        => clm3%g%l%c%cps%alb_can_EASS
    albsno_over_EASS    => clm3%g%l%c%cps%albsno_over_EASS
    albsno_under_EASS   => clm3%g%l%c%cps%albsno_under_EASS
    alb_grnd_EASS       => clm3%g%l%c%cps%alb_grnd_EASS
    albsno_grnd_EASS    => clm3%g%l%c%cps%albsno_grnd_EASS
    !-----------------------------

    ! Assign local pointers to derived type members (pft-level)

    plandunit     => clm3%g%l%c%p%landunit
    ivt           => clm3%g%l%c%p%itype
    pcolumn       => clm3%g%l%c%p%column
    pgridcell     => clm3%g%l%c%p%gridcell
    pwtgcell      => clm3%g%l%c%p%wtgcell
    elai          => clm3%g%l%c%p%pps%elai
    esai          => clm3%g%l%c%p%pps%esai
    slasun        => clm3%g%l%c%p%pps%slasun
    slasha        => clm3%g%l%c%p%pps%slasha
    gdir          => clm3%g%l%c%p%pps%gdir
    omega         => clm3%g%l%c%p%pps%omega
    laisun        => clm3%g%l%c%p%pps%laisun
    laisha        => clm3%g%l%c%p%pps%laisha
    fabd          => clm3%g%l%c%p%pps%fabd
    fabi          => clm3%g%l%c%p%pps%fabi
    ftdd          => clm3%g%l%c%p%pps%ftdd
    ftid          => clm3%g%l%c%p%pps%ftid
    ftii          => clm3%g%l%c%p%pps%ftii
    albd          => clm3%g%l%c%p%pps%albd
    albi          => clm3%g%l%c%p%pps%albi
    fsun          => clm3%g%l%c%p%pps%fsun
    sabg          => clm3%g%l%c%p%pef%sabg
    sabv          => clm3%g%l%c%p%pef%sabv
    snowdp        => clm3%g%l%c%cps%snowdp
    fsa           => clm3%g%l%c%p%pef%fsa
    fsa_r         => clm3%g%l%c%p%pef%fsa_r
    fsr           => clm3%g%l%c%p%pef%fsr
    parsun        => clm3%g%l%c%p%pef%parsun
    parsha        => clm3%g%l%c%p%pef%parsha
     fsds_vis_d    => clm3%g%l%c%p%pef%fsds_vis_d
     fsds_nir_d    => clm3%g%l%c%p%pef%fsds_nir_d
     fsds_vis_i    => clm3%g%l%c%p%pef%fsds_vis_i
     fsds_nir_i    => clm3%g%l%c%p%pef%fsds_nir_i
     fsr_vis_d     => clm3%g%l%c%p%pef%fsr_vis_d
     fsr_nir_d     => clm3%g%l%c%p%pef%fsr_nir_d
     fsr_vis_i     => clm3%g%l%c%p%pef%fsr_vis_i
     fsr_nir_i     => clm3%g%l%c%p%pef%fsr_nir_i
     fsds_vis_d_ln => clm3%g%l%c%p%pef%fsds_vis_d_ln
     fsds_nir_d_ln => clm3%g%l%c%p%pef%fsds_nir_d_ln
     fsr_vis_d_ln  => clm3%g%l%c%p%pef%fsr_vis_d_ln
     fsr_nir_d_ln  => clm3%g%l%c%p%pef%fsr_nir_d_ln
     eff_kid       => clm3%g%l%c%p%pps%eff_kid
     eff_kii       => clm3%g%l%c%p%pps%eff_kii
     sun_faid      => clm3%g%l%c%p%pps%sun_faid
     sun_faii      => clm3%g%l%c%p%pps%sun_faii
     sha_faid      => clm3%g%l%c%p%pps%sha_faid
     sha_faii      => clm3%g%l%c%p%pps%sha_faii
     sun_add       => clm3%g%l%c%p%pef%sun_add
     tot_aid       => clm3%g%l%c%p%pef%tot_aid
     sun_aid       => clm3%g%l%c%p%pef%sun_aid
     sun_aii       => clm3%g%l%c%p%pef%sun_aii
     sha_aid       => clm3%g%l%c%p%pef%sha_aid
     sha_aii       => clm3%g%l%c%p%pef%sha_aii
     sun_atot      => clm3%g%l%c%p%pef%sun_atot
     sha_atot      => clm3%g%l%c%p%pef%sha_atot
     sun_alf       => clm3%g%l%c%p%pef%sun_alf
     sha_alf       => clm3%g%l%c%p%pef%sha_alf
     sun_aperlai   => clm3%g%l%c%p%pef%sun_aperlai
     sha_aperlai   => clm3%g%l%c%p%pef%sha_aperlai
!----------------------------------------------------------------
     !added by Jing Chen, 26 March 2012     
     clumping_EASS         => clm3%g%l%c%p%pps%clumping_EASS
     elai_over_EASS        => clm3%g%l%c%p%pps%elai_over_EASS
     elai_over_shaded_EASS => clm3%g%l%c%p%pps%elai_over_shaded_EASS
     elai_over_sunlit_EASS => clm3%g%l%c%p%pps%elai_over_sunlit_EASS
     elai_under_EASS       => clm3%g%l%c%p%pps%elai_under_EASS
     elai_under_shaded_EASS=> clm3%g%l%c%p%pps%elai_under_shaded_EASS
     elai_under_sunlit_EASS=> clm3%g%l%c%p%pps%elai_under_sunlit_EASS
     esai_over_EASS        => clm3%g%l%c%p%pps%esai_over_EASS
     esai_under_EASS       => clm3%g%l%c%p%pps%esai_under_EASS
     lt_over_shaded_EASS   => clm3%g%l%c%p%pps%lt_over_shaded_EASS
     lt_over_sunlit_EASS   => clm3%g%l%c%p%pps%lt_over_sunlit_EASS
     lt_under_shaded_EASS  => clm3%g%l%c%p%pps%lt_under_shaded_EASS
     lt_under_sunlit_EASS  => clm3%g%l%c%p%pps%lt_under_sunlit_EASS
     albedo_nir_EASS       => pftcon%albedo_nir_EASS
     albedo_vis_EASS       => pftcon%albedo_vis_EASS
     sum_snow_over_EASS    => clm3%g%l%c%p%pwf%sum_snow_over_EASS
     sum_snow_under_EASS   => clm3%g%l%c%p%pwf%sum_snow_under_EASS
     sum_snow_EASS         => clm3%g%l%c%p%pwf%sum_snow_EASS
     dalb_over_EASS        => clm3%g%l%c%p%pef%dalb_over_EASS
     dalb_under_EASS       => clm3%g%l%c%p%pef%dalb_under_EASS
     dalb_grnd_EASS        => clm3%g%l%c%p%pef%dalb_grnd_EASS
     swnet_EASS            => clm3%g%l%c%p%pef%swnet_EASS
     elai_EASS             => clm3%g%l%c%p%pps%elai_EASS
!----------------------------------------------------------------
     
     ! Assign local pointers to derived type members (ecophysiological)

     slatop           => pftcon%slatop
     dsladlai         => pftcon%dsladlai
     albedo_snow_EASS => pftcon%albedo_snow_EASS
     frac_sno         => clm3%g%l%c%cps%frac_sno
     flx_absdv        => clm3%g%l%c%cps%flx_absdv
     flx_absdn        => clm3%g%l%c%cps%flx_absdn
     flx_absiv        => clm3%g%l%c%cps%flx_absiv
     flx_absin        => clm3%g%l%c%cps%flx_absin
     sabg_lyr         => clm3%g%l%c%p%pef%sabg_lyr
     snl              => clm3%g%l%c%cps%snl
     sfc_frc_aer      => clm3%g%l%c%p%pef%sfc_frc_aer
     sfc_frc_aer_sno  => clm3%g%l%c%p%pef%sfc_frc_aer_sno
     albgrd_pur       => clm3%g%l%c%cps%albgrd_pur
     albgri_pur       => clm3%g%l%c%cps%albgri_pur
     sfc_frc_bc       => clm3%g%l%c%p%pef%sfc_frc_bc
     sfc_frc_bc_sno   => clm3%g%l%c%p%pef%sfc_frc_bc_sno
     albgrd_bc        => clm3%g%l%c%cps%albgrd_bc
     albgri_bc        => clm3%g%l%c%cps%albgri_bc
     sfc_frc_oc       => clm3%g%l%c%p%pef%sfc_frc_oc
     sfc_frc_oc_sno   => clm3%g%l%c%p%pef%sfc_frc_oc_sno
     albgrd_oc        => clm3%g%l%c%cps%albgrd_oc
     albgri_oc        => clm3%g%l%c%cps%albgri_oc
     sfc_frc_dst      => clm3%g%l%c%p%pef%sfc_frc_dst
     sfc_frc_dst_sno  => clm3%g%l%c%p%pef%sfc_frc_dst_sno
     albgrd_dst       => clm3%g%l%c%cps%albgrd_dst
     albgri_dst       => clm3%g%l%c%cps%albgri_dst
     albsnd_hst       => clm3%g%l%c%cps%albsnd_hst
     albsni_hst       => clm3%g%l%c%cps%albsni_hst
     fsr_sno_vd       => clm3%g%l%c%p%pef%fsr_sno_vd
     fsr_sno_nd       => clm3%g%l%c%p%pef%fsr_sno_nd
     fsr_sno_vi       => clm3%g%l%c%p%pef%fsr_sno_vi
     fsr_sno_ni       => clm3%g%l%c%p%pef%fsr_sno_ni
     fsds_sno_vd      => clm3%g%l%c%p%pef%fsds_sno_vd
     fsds_sno_nd      => clm3%g%l%c%p%pef%fsds_sno_nd
     fsds_sno_vi      => clm3%g%l%c%p%pef%fsds_sno_vi
     fsds_sno_ni      => clm3%g%l%c%p%pef%fsds_sno_ni
     !---------------------------------
     ! added by Jing Chen Aug 30 2012
     rhol      => pftcon%rhol
     lai_over_max_EASS    => pftcon%lai_over_max_EASS
     lai_under_max_EASS   => pftcon%lai_under_max_EASS
     clumping_pft_EASS    => pftcon%clumping_pft_EASS
     !-------------------------------------

     ! Determine seconds off current time step
     
     dtime = get_step_size()
     call get_curr_date (year, month, day, secs)

     ! Determine fluxes

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        ! was redundant b/c filter already included wt>0; 
        ! not redundant anymore with chg in filter definition
        l = plandunit(p)
        !Note: Some glacier_mec pfts may have zero weight
        if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           sabg(p)       = 0._r8
           sabv(p)       = 0._r8
           fsa(p)        = 0._r8
           !---------------
           !added by Jing Chen, 26 May,2012
           dalb_over_EASS(p)  = 0.0_r8
           dalb_under_EASS(p) = 0.0_r8
           dalb_grnd_EASS(p)  = 0.0_r8
           !--------------
           l = plandunit(p)
           if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
             fsa_r(p)    = 0._r8
           end if
           sabg_lyr(p,:) = 0._r8
           sabg_pur(p)   = 0._r8
           sabg_bc(p)    = 0._r8
           sabg_oc(p)    = 0._r8
           sabg_dst(p)   = 0._r8
        end if
     end do 

     ! Loop over pfts to calculate fsun, etc
     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        l = plandunit(p)
        if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           c = pcolumn(p)
           g = pgridcell(p)
        
           vai(p) = elai(p) + esai(p)

           if (coszen(c) > 0._r8 .and. elai(p) > 0._r8 .and. gdir(p) > 0._r8) then
              cosz = max(0.001_r8, coszen(c))

              ext = gdir(p)/cosz              !CLM   for Sen.

              t1 = min(ext*elai(p), 40.0_r8)  !CLM
              t2 = exp(-t1)
              fsun(p) = (1._r8-t2)/t1
              
    
              ! new control on low lai, to avoid numerical problems in
              ! calculation of slasun, slasha
              ! PET: 2/29/04
              
              if (elai(p) > 0.01_r8) then
                 laisun(p) = elai(p)*fsun(p)
                 laisha(p) = elai(p)*(1._r8-fsun(p))

                 ! calculate the average specific leaf area for sunlit and shaded
                 ! canopies, when effective LAI > 0
                 slasun(p) = (t2*dsladlai(ivt(p))*ext*elai(p) + &
                              t2*dsladlai(ivt(p)) + &
                              t2*slatop(ivt(p))*ext - &
                              dsladlai(ivt(p)) - &
                              slatop(ivt(p))*ext) / &
                              (ext*(t2-1._r8))
                 slasha(p) = ((slatop(ivt(p)) + &
                             (dsladlai(ivt(p)) * elai(p)/2.0_r8)) * elai(p) - &
                             laisun(p)*slasun(p)) / laisha(p)
              else
                 ! special case for low elai
                 fsun(p) = 1._r8
                 laisun(p) = elai(p)
                 laisha(p) = 0._r8
                 slasun(p) = slatop(ivt(p))
                 slasha(p) = 0._r8
              end if
           else
              fsun(p)   = 0._r8
              laisun(p) = 0._r8
              laisha(p) = elai(p)
              slasun(p) = 0._r8
              slasha(p) = 0._r8
           end if

           !-------------------------------------------------------
           !added by Jing Chen, 26 March 2012

#if (defined CN)
           ! recalculate LAI
           if ( elai(p) == 0._r8 ) then
               elai_EASS(p)=0.01_r8
           else
               elai_EASS(p)=elai(p)
           end if

           ! divide the elai and esai into overstory and understory

           elai_over_EASS(p)=elai_EASS(p)
           if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
                if ( ivt(p)==nc4_grass .or. ivt(p)==ncorn ) then   !SEE module pftvarcon
                       elai_under_EASS(p)=0.1_r8
                else
                       elai_under_EASS(p)=1.18_r8*exp(-0.99_r8*elai_over_EASS(p))
                end if
                if (elai_under_EASS(p) > elai_over_EASS(p)) elai_under_EASS(p)=0.01_r8
           else
                elai_under_EASS(p)=0._r8
           end if

           esai_over_EASS(p)=esai(p)
           !#  esai_over_EASS(p)= lai_over_max_EASS(ivt(p))*0.2_r8   !!algorithm in EASS
           esai_under_EASS(p)= lai_under_max_EASS(ivt(p))*0.2_r8     !algroithm in EASS

#endif

           ! (LAI+SAI)and LAI of overstory canopy(sunlit and shaded leaves)
           phase='LAI'
           lt_over_EASS = elai_over_EASS(p) + esai_over_EASS(p)
           if ( coszen(c)>0.0_r8 ) then
                    call SunlitLAI_EASS(lt_over_EASS,clumping_pft_EASS(ivt(p)), coszen(c), phase, &
                                        lt_over_sunlit_EASS(p))     ! LAI_o_sunlit
                    call SunlitLAI_EASS(elai_over_EASS(p), clumping_pft_EASS(ivt(p)), coszen(c), phase,&
                                        elai_over_sunlit_EASS(p))   ! LAIo_sunlit	
           else
                   lt_over_sunlit_EASS(p)   = 0._r8
                   elai_over_sunlit_EASS(p) = 0._r8
           end if
           lt_over_shaded_EASS(p) = lt_over_EASS - lt_over_sunlit_EASS(p)             ! to calculate temperature of overstory canopy
           elai_over_shaded_EASS(p) = max ( 0._r8, elai_over_EASS(p) - elai_over_sunlit_EASS(p) )  ! to calculate GPP of overstory canopy

           ! (LAI+SAI)and LAI of understory canopy(sunlit and shaded leaves)
           lt_under_EASS = elai_under_EASS(p) + esai_under_EASS(p)
           lt_inter = elai_over_EASS(p) + esai_over_EASS(p) + elai_under_EASS(p) + esai_under_EASS(p)
           elai_can_inter = elai_over_EASS(p) + elai_under_EASS(p)

           if (coszen(c)>0.0_r8) then
                   call SunlitLAI_EASS(lt_inter, clumping_pft_EASS(ivt(p)), coszen(c),  phase, lt_EASS)    
                   lt_under_sunlit_EASS(p)=lt_EASS - lt_over_sunlit_EASS(p)       !LAI_u_sunlit

                   call SunlitLAI_EASS(elai_can_inter, clumping_pft_EASS(ivt(p)), coszen(c), phase, elai_over_under_EASS)  !elai_over_under_EASS(p)
                   elai_under_sunlit_EASS(p) = elai_over_under_EASS-elai_over_sunlit_EASS(p) !LAIu_sunli
           else
                   lt_under_sunlit_EASS(p)= 0._r8
                   elai_under_sunlit_EASS(p)= 0.0_r8
           end if
           lt_under_shaded_EASS(p)= lt_under_EASS - lt_under_sunlit_EASS(p)                    ! to calculate temperature of understory canopy
           elai_under_shaded_EASS(p)=max(0._r8, elai_under_EASS(p)-lt_under_sunlit_EASS(p))    ! to calculate GPP of understory canopy

          !-------------------------------------------------------

        end if
     end do
        
     ! Loop over nband wavebands
     do ib = 1, nband
        do fp = 1,num_nourbanp
           p = filter_nourbanp(fp)
           l = plandunit(p)
           if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
              c = pcolumn(p)
              g = pgridcell(p)
              
              ! Absorbed by canopy

              cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
              cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
              sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
              fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)

              l = plandunit(p)
              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
              end if
              
              ! Transmitted = solar fluxes incident on ground
              
              trd(p,ib) = forc_solad(g,ib)*ftdd(p,ib)
              tri(p,ib) = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)
     
              ! Solar radiation absorbed by ground surface
 
              absrad  = trd(p,ib)*(1._r8-albgrd(c,ib)) + tri(p,ib)*(1._r8-albgri(c,ib))
              sabg(p) = sabg(p) + absrad
              fsa(p)  = fsa(p)  + absrad
              !----------------
              ! added by Jing Chen
              swnet_EASS(p) = fsa(p)
              !-----------------
              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + absrad
              end if

#if (defined SNICAR_FRC)
              ! Solar radiation absorbed by ground surface without BC
              absrad_bc = trd(p,ib)*(1._r8-albgrd_bc(c,ib)) + tri(p,ib)*(1._r8-albgri_bc(c,ib))
              sabg_bc(p) = sabg_bc(p) + absrad_bc

              ! Solar radiation absorbed by ground surface without OC
              absrad_oc = trd(p,ib)*(1._r8-albgrd_oc(c,ib)) + tri(p,ib)*(1._r8-albgri_oc(c,ib))
              sabg_oc(p) = sabg_oc(p) + absrad_oc

              ! Solar radiation absorbed by ground surface without dust
              absrad_dst = trd(p,ib)*(1._r8-albgrd_dst(c,ib)) + tri(p,ib)*(1._r8-albgri_dst(c,ib))
              sabg_dst(p) = sabg_dst(p) + absrad_dst

              ! Solar radiation absorbed by ground surface without any aerosols
              absrad_pur = trd(p,ib)*(1._r8-albgrd_pur(c,ib)) + tri(p,ib)*(1._r8-albgri_pur(c,ib))
              sabg_pur(p) = sabg_pur(p) + absrad_pur
#endif

              ! New sunlit.shaded canopy algorithm
              
              if (coszen(c) > 0._r8 .and. elai(p) > 0._r8 .and. gdir(p) > 0._r8 ) then
                 
                 ! 1. calculate flux of direct beam radiation absorbed in the 
                 ! sunlit canopy as direct (sun_add), and the flux of direct
                 ! beam radiation absorbed in the total canopy as indirect
                 
                 sun_add(p,ib) = forc_solad(g,ib) * (1._r8-ftdd(p,ib)) * (1._r8-omega(p,ib))
                 tot_aid(p,ib) = (forc_solad(g,ib) * fabd(p,ib)) - sun_add(p,ib)
                 
                 ! the following constraint set to catch round-off level errors
                 ! that can cause negative tot_aid
                 
                 tot_aid(p,ib) = max(tot_aid(p,ib), 0._r8)
                 
                 ! 2. calculate the effective extinction coefficients for indirect
                 ! transmission originating from direct and indirect streams,
                 ! using ftid and ftii
                 
                 !eff_kid(p,ib) = -(log(ftid(p,ib)))/vai(p)
                 !eff_kii(p,ib) = -(log(ftii(p,ib)))/vai(p)
                 
                 ! 3. calculate the fraction of indirect radiation being absorbed 
                 ! in the sunlit and shaded canopy fraction. Some of this indirect originates in
                 ! the direct beam and some originates in the indirect beam.

                 !sun_faid(p,ib) = 1.-exp(-eff_kid(p,ib) * vaisun(p))
                 !sun_faii(p,ib) = 1.-exp(-eff_kii(p,ib) * vaisun(p))
                 sun_faid(p,ib) = fsun(p)
                 sun_faii(p,ib) = fsun(p)
                 sha_faid(p,ib) = 1._r8-sun_faid(p,ib)
                 sha_faii(p,ib) = 1._r8-sun_faii(p,ib)

                 ! 4. calculate the total indirect flux absorbed by the sunlit
                 ! and shaded canopy based on these fractions and the fabd and
                 ! fabi from surface albedo calculations

                 sun_aid(p,ib) = tot_aid(p,ib) * sun_faid(p,ib)
                 sun_aii(p,ib) = forc_solai(g,ib)*fabi(p,ib)*sun_faii(p,ib)
                 sha_aid(p,ib) = tot_aid(p,ib) * sha_faid(p,ib)
                 sha_aii(p,ib) = forc_solai(g,ib)*fabi(p,ib)*sha_faii(p,ib)
                 
                 ! 5. calculate the total flux absorbed in the sunlit and shaded
                 ! canopy as the sum of these terms
                 
                 sun_atot(p,ib) = sun_add(p,ib) + sun_aid(p,ib) + sun_aii(p,ib)
                 sha_atot(p,ib) = sha_aid(p,ib) + sha_aii(p,ib)
                 
                 ! 6. calculate the total flux absorbed by leaves in the sunlit
                 ! and shaded canopies
                 
                 laifra = elai(p)/vai(p)
                 sun_alf(p,ib) = sun_atot(p,ib) * laifra
                 sha_alf(p,ib) = sha_atot(p,ib) * laifra
                 
                 ! 7. calculate the fluxes per unit lai in the sunlit and shaded
                 ! canopies
                 
                 if (laisun(p) > 0._r8) then
                    sun_aperlai(p,ib) = sun_alf(p,ib)/laisun(p)
                 else
                    sun_aperlai(p,ib) = 0._r8
                 endif
                 if (laisha(p) > 0._r8) then
                    sha_aperlai(p,ib) = sha_alf(p,ib)/laisha(p)
                 else
                    sha_aperlai(p,ib) = 0._r8
                 endif

              else   ! coszen = 0 or elai = 0
                 
                 sun_add(p,ib)     = 0._r8
                 tot_aid(p,ib)     = 0._r8
                 eff_kid(p,ib)     = 0._r8
                 eff_kii(p,ib)     = 0._r8
                 sun_faid(p,ib)    = 0._r8
                 sun_faii(p,ib)    = 0._r8
                 sha_faid(p,ib)    = 0._r8
                 sha_faii(p,ib)    = 0._r8
                 sun_aid(p,ib)     = 0._r8
                 sun_aii(p,ib)     = 0._r8
                 sha_aid(p,ib)     = 0._r8
                 sha_aii(p,ib)     = 0._r8
                 sun_atot(p,ib)    = 0._r8
                 sha_atot(p,ib)    = 0._r8
                 sun_alf(p,ib)     = 0._r8
                 sha_alf(p,ib)     = 0._r8
                 sun_aperlai(p,ib) = 0._r8
                 sha_aperlai(p,ib) = 0._r8
                 
              end if
              !-------------------------------
              !added by Jing Chen, 18 May 2012

              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then

                 snow_over_EASS(p)=forc_snow(g)/1000.0_r8    ! snow rate , unit change[mm s-1]->[m s-1]

                 !albedo of snow 
                 if (snow_over_EASS(p) .gt. 1.0e-10_r8) sum_snow_EASS(p)=sum_snow_EASS(p)+snow_over_EASS(p)

                 if ((sum_snow_EASS(p) .le. 0._r8 ) .and. (snowdp_EASS(c) .le. 0._r8)) then
                        ! no snowing
                        albsnow_EASS(c,ib)=0._r8
                 elseif ( (sum_snow_EASS(p) .gt. 0.01_r8) .and. (snowdp_EASS(c) .ge. 0.0001_r8) ) then
                        ! new snowing
                        albsnow_EASS(c,ib) = albedo_new_snow_EASS(ib) 
                 else
                        ! no new snowing, the snow albedo is decrease with snow older
                        albsnow_EASS(c,ib) = (albedo_snow_EASS(ivt(p)) - albsnow_const_EASS(ib)) * &
                                          exp(-0.005_r8 * dtime / 3600.0_r8) + albsnow_const_EASS(ib)
                 end if 
  
                 !albedo of canopy without snow covered
                 if ( (coszen(c) .gt. 0._r8) .and. (forc_solar(g) .gt. 0._r8) ) then
                     if (ib==1)  alb_can_EASS(c,ib)=albedo_vis_EASS(ivt(p))    
                     if (ib==2)  alb_can_EASS(c,ib)=albedo_nir_EASS(ivt(p))    
                 else
                     alb_can_EASS(c,ib)=0.0_r8
                 end if

                 !albedo of the canopy covered by snow
                 albsno_over_EASS(c,ib)=alb_can_EASS(c,ib)*(1.0_r8-sum_snow_over_EASS(p)/elai_over_EASS(p)/2.0_r8) &
                                    +albsnow_EASS(c,ib)*sum_snow_over_EASS(p)/elai_over_EASS(p)/2.0_r8

                 albsno_under_EASS(c,ib)=alb_can_EASS(c,ib)*(1.0_r8-sum_snow_under_EASS(p)/elai_under_EASS(p)/2.0_r8) &
                                     +albsnow_EASS(c,ib)*sum_snow_under_EASS(p)/elai_under_EASS(p)/2.0_r8

                 ! albedo of ground without snow covered
                 !in subroutine SoilAlbedo(module SurfaceAlbedoMod)
                 !albsod(c,ib)=albsoi(c,ib)   ! albsod, direct-beam soil albedo; albsoi,diffuse soil albedo
   
                 if ( (coszen(c) .gt. 0._r8) .and. (forc_solar(g) .gt. 0._r8) ) then
                       alb_grnd_EASS(c,ib)=albsod(c,ib)   ! 2=nir
     
                       !albedo of ground with snow covered
                       albsno_grnd_EASS(c,ib)=alb_grnd_EASS(c,ib)*(1.0_r8-frac_snow_grnd_EASS(c)) &
                                           +albsnow_EASS(c,ib)*frac_snow_grnd_EASS(c)
                 else
                       albsno_grnd_EASS(c,ib)=0._r8
                 end if

                 ! albedo of overstory , understory and ground and with the difference between 1 
                 dalb_over_EASS(p)  = dalb_over_EASS(p) + 0.5_r8 * ( 1.0_r8 - albsno_over_EASS(c,ib) )
                 dalb_under_EASS(p) = dalb_under_EASS(p) + 0.5_r8 * ( 1.0_r8 - albsno_under_EASS(c,ib) )
                 dalb_grnd_EASS(p)  = dalb_grnd_EASS(p) + 0.5_r8 * ( 1.0_r8 - albsno_grnd_EASS(c,ib) )

              end if
              !-----------------------------
           
           end if   !(pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec)
        end do ! end of pft loop
     end do ! end nbands loop   

     !   compute absorbed flux in each snow layer and top soil layer,
     !   based on flux factors computed in the radiative transfer portion of SNICAR.
     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
           l = plandunit(p)
           if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           c = pcolumn(p)
           sabg_snl_sum = 0._r8

           ! CASE1: No snow layers: all energy is absorbed in top soil layer
           if (snl(c) == 0) then
              sabg_lyr(p,:) = 0._r8
              sabg_lyr(p,1) = sabg(p)
              sabg_snl_sum  = sabg_lyr(p,1)
   
           ! CASE 2: Snow layers present: absorbed radiation is scaled according to 
           ! flux factors computed by SNICAR
           else
              do i = -nlevsno+1,1,1
                 sabg_lyr(p,i) = flx_absdv(c,i)*trd(p,1) + flx_absdn(c,i)*trd(p,2) + &
                                 flx_absiv(c,i)*tri(p,1) + flx_absin(c,i)*tri(p,2)
                 ! summed radiation in active snow layers:
                 if (i >= snl(c)+1) then
                    sabg_snl_sum = sabg_snl_sum + sabg_lyr(p,i)
                 endif
              enddo
   
              ! Error handling: The situation below can occur when solar radiation is 
              ! NOT computed every timestep.
              ! When the number of snow layers has changed in between computations of the 
              ! absorbed solar energy in each layer, we must redistribute the absorbed energy
              ! to avoid physically unrealistic conditions. The assumptions made below are 
              ! somewhat arbitrary, but this situation does not arise very frequently. 
              ! This error handling is implemented to accomodate any value of the
              ! radiation frequency.
              if (abs(sabg_snl_sum-sabg(p)) > 0.00001_r8) then
                 if (snl(c) == 0) then
                    sabg_lyr(p,-4:0) = 0._r8
                    sabg_lyr(p,1) = sabg(p)
                 elseif (snl(c) == -1) then
                    sabg_lyr(p,-4:-1) = 0._r8
                    sabg_lyr(p,0) = sabg(p)*0.6_r8
                    sabg_lyr(p,1) = sabg(p)*0.4_r8
                 else
                    sabg_lyr(p,:) = 0._r8
                    sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                    sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                 endif
              endif

              ! If shallow snow depth, all solar radiation absorbed in top or top two snow layers
              ! to prevent unrealistic timestep soil warming 
              if (snowdp(c) < 0.10_r8) then
                 if (snl(c) == 0) then
                    sabg_lyr(p,-4:0) = 0._r8
                    sabg_lyr(p,1) = sabg(p)
                 elseif (snl(c) == -1) then
                    sabg_lyr(p,-4:-1) = 0._r8
                    sabg_lyr(p,0) = sabg(p)
                    sabg_lyr(p,1) = 0._r8
                 else
                    sabg_lyr(p,:) = 0._r8
                    sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                    sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                 endif
              endif

           endif

           ! This situation should not happen:
           if (abs(sum(sabg_lyr(p,:))-sabg(p)) > 0.00001_r8) then
              write(iulog,*) "SNICAR ERROR: Absorbed ground radiation not equal to summed snow layer radiation. pft = ",   &
                             p," Col= ", c, " Diff= ",sum(sabg_lyr(p,:))-sabg(p), " sabg(p)= ", sabg(p), " sabg_sum(p)= ", &
                             sum(sabg_lyr(p,:)), " snl(c)= ", snl(c)
              write(iulog,*) "flx_absdv1= ", trd(p,1)*(1.-albgrd(c,1)), "flx_absdv2= ", sum(flx_absdv(c,:))*trd(p,1)
              write(iulog,*) "flx_absiv1= ", tri(p,1)*(1.-albgri(c,1))," flx_absiv2= ", sum(flx_absiv(c,:))*tri(p,1)
              write(iulog,*) "flx_absdn1= ", trd(p,2)*(1.-albgrd(c,2))," flx_absdn2= ", sum(flx_absdn(c,:))*trd(p,2)
              write(iulog,*) "flx_absin1= ", tri(p,2)*(1.-albgri(c,2))," flx_absin2= ", sum(flx_absin(c,:))*tri(p,2)
   
              write(iulog,*) "albgrd_nir= ", albgrd(c,2)
              write(iulog,*) "coszen= ", coszen(c)
              call endrun()
           endif

 
#if (defined SNICAR_FRC)

           ! BC aerosol forcing (pft-level):
           sfc_frc_bc(p) = sabg(p) - sabg_bc(p)
   
           ! OC aerosol forcing (pft-level):
           if (DO_SNO_OC) then
              sfc_frc_oc(p) = sabg(p) - sabg_oc(p)
           else
              sfc_frc_oc(p) = 0._r8
           endif
   
           ! dust aerosol forcing (pft-level):
           sfc_frc_dst(p) = sabg(p) - sabg_dst(p)
   
           ! all-aerosol forcing (pft-level):
           sfc_frc_aer(p) = sabg(p) - sabg_pur(p)        
           
           ! forcings averaged only over snow:
           if (frac_sno(c) > 0._r8) then
              sfc_frc_bc_sno(p)  = sfc_frc_bc(p)/frac_sno(c)
              sfc_frc_oc_sno(p)  = sfc_frc_oc(p)/frac_sno(c)
              sfc_frc_dst_sno(p) = sfc_frc_dst(p)/frac_sno(c)
              sfc_frc_aer_sno(p) = sfc_frc_aer(p)/frac_sno(c)
           else
              sfc_frc_bc_sno(p)  = spval
              sfc_frc_oc_sno(p)  = spval
              sfc_frc_dst_sno(p) = spval
              sfc_frc_aer_sno(p) = spval
           endif

#endif
        endif
        
     enddo


     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        l = plandunit(p)
        if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           g = pgridcell(p)
        
           ! Final step of new sunlit/shaded canopy algorithm
           ! 8. calculate the total and per-unit-lai fluxes for PAR in the
           ! sunlit and shaded canopy leaf fractions
           
           parsun(p) = sun_aperlai(p,1)
           parsha(p) = sha_aperlai(p,1)
           
           ! The following code is duplicated from SurfaceRadiation
           ! NDVI and reflected solar radiation
           
           rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
           rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
           fsr(p) = rvis + rnir
           
           fsds_vis_d(p) = forc_solad(g,1)
           fsds_nir_d(p) = forc_solad(g,2)
           fsds_vis_i(p) = forc_solai(g,1)
           fsds_nir_i(p) = forc_solai(g,2)
           fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
           fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
           fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
           fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)
           
           local_secp1 = secs + nint((londeg(g)/degpsec)/dtime)*dtime
           local_secp1 = mod(local_secp1,isecspday)
           if (local_secp1 == isecspday/2) then
              fsds_vis_d_ln(p) = forc_solad(g,1)
              fsds_nir_d_ln(p) = forc_solad(g,2)
              fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
              fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
           else
              fsds_vis_d_ln(p) = spval
              fsds_nir_d_ln(p) = spval
              fsr_vis_d_ln(p) = spval
              fsr_nir_d_ln(p) = spval
           end if

           ! diagnostic variables (downwelling and absorbed radiation partitioning) for history files
           ! (OPTIONAL)
           c = pcolumn(p)
           if (snl(c) < 0) then
              fsds_sno_vd(p) = forc_solad(g,1)
              fsds_sno_nd(p) = forc_solad(g,2)
              fsds_sno_vi(p) = forc_solai(g,1)
              fsds_sno_ni(p) = forc_solai(g,2)

              fsr_sno_vd(p) = fsds_vis_d(p)*albsnd_hst(c,1)
              fsr_sno_nd(p) = fsds_nir_d(p)*albsnd_hst(c,2)
              fsr_sno_vi(p) = fsds_vis_i(p)*albsni_hst(c,1)
              fsr_sno_ni(p) = fsds_nir_i(p)*albsni_hst(c,2)
           else
              fsds_sno_vd(p) = spval
              fsds_sno_nd(p) = spval
              fsds_sno_vi(p) = spval
              fsds_sno_ni(p) = spval

              fsr_sno_vd(p) = spval
              fsr_sno_nd(p) = spval
              fsr_sno_vi(p) = spval
              fsr_sno_ni(p) = spval
           endif

        end if
     end do 

   end subroutine SurfaceRadiation

   !-------------------------------------
   !! BOP
   !! IROUTINE: SunlitLAI_EASS
   !! INTERFACE:   
   subroutine SunlitLAI_EASS(LAI,clumping, cosz,  phase, canleaf)

  !! DESCRIPTION:
  !  calculate the LAI of sunlit leaves based on clumping index in EASS model

  !! USE:

  !! ARGUMENT:
        implicit none
	real(r8), intent(in) :: LAI         ! LAI
	real(r8), intent(in) :: clumping    !clumping index in EASS(input)
	real(r8), intent(in) :: cosz        ! cosine of solar zenith angle						
        character(len=*), intent(in):: phase    ! "LAI" in subroutine LeafLAI_EASS or "energy" in subroutine EnergyBlance_EASS
	real(r8), intent(out):: canleaf         ! output var as canlai_EASS or fsunfleck_EASS, depending on phase

  !! CALLED FROM:
  ! subroutine SurfaceRadiation in this module
  ! subroutine EnergyBlance_EASS in EnergyBlance_EASSMod

  !! REVERTION HISTORY:
  ! 5/17/12, Jing Chen: calculate sunlit leaves index based on clumping index. the c++ version is EASS_region_1028.

  !! LOCAL VARIABLES:
  ! local pointers to implicite in varialbes
  ! local pointers to implicite inout arguments
  ! local pointers to implicite out arguments

  !! OTHER LOCAL VARIALBES:
	real(r8):: leaf_pro                  ! leaf projection
	real(r8):: t1, t2                    ! temporary var.
	real(r8):: fsun_EASS                 ! sunlit fraction of canopy in EASS 
	real(r8):: fsunfleck_EASS            ! tractional area of sunflecks on the horizontal plane below the leaf area index
        real(r8):: lt_over_sunlit_EASS       ! elai+esai of sunlit leaves, overstory
        real(r8):: elai_over_sunlit_EASS     ! elai of sunlit leaves, overstory
        real(r8):: lt_EASS                   ! elai_over+esai_over+elai_under+esai_under
        real(r8):: elai_over_under_EASS      ! elai_over+elai_under
	real(r8):: canlai_EASS               ! lai of canopy

    !--------------------------------------------------------------			

        leaf_pro =0.5_r8* clumping
        t1=min(leaf_pro /cosz * LAI, 40.0_r8)
        t2=1.0_r8-exp(-t1)

        canlai_EASS=2.0_r8*cosz*t2
        fsunfleck_EASS=exp(-t1)

        ! to judge the output value depending on phses
        if (phase=='LAI') then
                canleaf=canlai_EASS
        else
                canleaf=fsunfleck_EASS
        end if

   end subroutine SunlitLAI_EASS

end module SurfaceRadiationMod
