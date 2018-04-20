module LeafTemper_EASSMod
!--------------------------------------------------
!BOP
!! USE:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use clm_varctl         , only : iulog

!! PUBLIC TYPES:
   implicit none
   save

!! PUBLIC MEMBER FUNCTIONS:
   public :: LeafTemper_EASS         ! calculate air resistnace and leaf boundry resistance

!! RVESION HISTORY

!EOP
!--------------------------------------------------
contains
!-------------------------------------------------
!! BOP
!! IROUTINE: AirResist_EASS
!! INTERFACE
   subroutine LeafTemper_EASS(p, c, g, es_EASS, espar_EASS, cpair_EASS, frac_cov,dayl_factor, phase)
          
!!DESCRIPTION:

!! USES:
    use clmtype
    use PenmanEqu_EASSMod    , only : LeafFlux_EASS
    use PhotoSynth_EASSMod   , only : PhotoSynth_EASS
    use clm_varcon           , only : tfrz, rgas      ! freezing T of fresh water[K] 
    use clm_atmlnd           , only : clm_a2l

!! ARGUMENTS:
    implicit none
!#    integer , intent(in) :: fn           ! size of pft filter
!#    integer , intent(in) :: filterp(fn)  ! pft filter
!#    integer , intent(in) :: lbp, ubp     ! pft bounds
    real(r8), intent(in) :: frac_cov       ! fraction of fall coverd on leaf, to calculate the old value of traspiration flux
    character(len=*), intent(in) :: phase  ! 'over_sunlit', 'over_shaded','under_sunlit' or 'under_shaded'
    integer , intent(in) :: p              ! pft index
    integer , intent(in) :: c              ! column index
    integer , intent(in) :: g              ! gridcell index
    real(r8), intent(in) :: cpair_EASS     ! specific heat of dry air in EASS[J kg-1 k-1],
                                           ! depanding on saturated water vapor specific humidity
    real(r8), intent(in) :: es_EASS        ! saturated water vapor[pa]
    real(r8), intent(in) :: espar_EASS     ! unsaturated water vapor[pa]
    real(r8), intent(in) :: dayl_factor    ! scalar (0-1) for daylength

!! CALLED FROM:

!! LOCAL VARIABLES:
! local pointer to implicit in variables
!#    integer , pointer :: pcolumn(:)        ! pft's column index
!#    integer , pointer :: pgridcell(:)      ! pft's gridcell index
    real(r8), pointer :: forc_pco2(:)      ! partial pressure co2 [Pa]
    real(r8), pointer :: forc_t(:)         ! temperature [K]
    real(r8), pointer :: forc_rho(:)       ! surface air density[kg m-3]    p17
    real(r8), pointer :: coszen(:)         ! cosine of solar zenith angle
    real(r8), pointer :: clumping_EASS(:)  ! clumping index in EASS
    real(r8), pointer :: elai_EASS(:)      ! one-sided leaf area index with burying by snow, 
                                           ! corrected by clumping index, in EASS
    real(r8), pointer :: Tleaf_EASS(:)     ! leaf temperature  [K]
    real(r8), pointer :: par_EASS(:)       ! average absorbed PAR by leaves
    real(r8), pointer :: ra_EASS(:)        ! air resitance[s m-1]
    real(r8), pointer :: rb_EASS(:)        ! leaf boundary resistance[s m-1]
    real(r8), pointer :: frac_rain_can(:)  ! fraction of rainfall on the canopy[0~1]
    real(r8), pointer :: frac_snow_can(:)  ! fraction of snow on the canopy(0~1)
    real(r8), pointer :: radsb_EASS(:)     ! total solar radiation aborbed by leaves[W m-2]
    real(r8), pointer :: Gheat_EASS(:)     ! heat conductance[m s-1]
    real(r8), pointer :: Gevap_EASS(:)     ! intercepted water conductance[m s-1]
    real(r8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
    real(r8), pointer :: c3psn(:)          ! photosynthetic pathway: 0. = c4, 1. = c3
    integer , pointer :: ivt(:)            ! pft vegetation type

! local pointer to implicit inout variables 
    real(r8), pointer :: Gleaf_EASS(:)     ! leaf conductance [m s-1]
    real(r8), pointer :: ci_EASS(:)        ! intercellular co2 concentration of  leaves [pa]

! local pointer to implicit out variables  
    real(r8), pointer :: psn_EASS(:)       ! foliage photosynthesis [umol co2 m-2 s-1]
    real(r8), pointer :: cs_EASS(:)        ! co2 concentration on the surfaces of leaves[pa]
    real(r8), pointer :: rs_EASS(:)        ! stomatal resistant leaves[s m-1]
    real(r8), pointer :: Tleaf_new_EASS(:) ! updated leaf temperature of leaves[K]
    real(r8), pointer :: vcmx_EASS(:)      ! maximum rate of carboxylation(for debug)[umol CO2 m-2 s-1](for debug)
    real(r8), pointer :: jmx_EASS(:)       ! light-saturated rate of electron transport 
                                           ! in the photosynthetic carbon reduction cycle in leaf cells(for debug)[umol m-2 s-1]
    real(r8), pointer :: dTleaf_EASS(:)    ! leaf temperature change(for debug)[K]
    !------------

!! OTHER LOCAL VAIRABLES:
!#    integer :: f                  ! filter index
!#    integer :: p                  ! pft index
!#    integer :: c                  ! column index
!#    integer :: g                  ! gridcell index
    real(r8) :: qflx_trans_EASS     ! transpiration water vapor flux, based on un-update leaf temperature[w m-2]
    real(r8) :: de_EASS                    ! unsaturated atmopsheric pressure [kpa]
    real(r8) :: slope_EASS                 ! slope of vapor pressure curve [kpa K-1]
    real(r8) :: p_star
    real(r8) :: cf                         ! s m2/umol -> s/m
    real(r8), parameter :: psy=0.066_r8    ! psychrometer constant[kPa K-1]
                                           ! at a standard pressure of 101.kPa, the value is about 66 Pa/K at 0,
                                           ! and increase to 67 Pa/K at 20C
    real(r8), parameter :: bp=2000._r8        ! minimum leaf conductance (umol/m2/s)
    real(r8), parameter :: rsmax0 = 2.e4_r8   ! maximum stomatal resistance [s/m]

    ! Assign local pointers to pft constants
    c3psn            => pftcon%c3psn

    ! Assign local pointers to derived type members(gridcell-level)
    forc_rho        => clm_a2l%forc_rho
    forc_pco2       => clm_a2l%forc_pco2
    forc_pbot       => clm_a2l%forc_pbot

    ! Assign local pointers to derived type members(column-level)
    forc_t          => clm3%g%l%c%ces%forc_t
    coszen          => clm3%g%l%c%cps%coszen

    ! Assign local pointers to derived type members(pft-level)
 !# pcolumn         => clm3%g%l%c%p%column
 !# pgridcell       => clm3%g%l%c%p%gridcell
    clumping_EASS   => clm3%g%l%c%p%pps%clumping_EASS
    elai_EASS       => clm3%g%l%c%p%pps%elai_EASS
    ivt             => clm3%g%l%c%p%itype

    if (phase=='over_sunlit') then
          Tleaf_EASS     => clm3%g%l%c%p%pes%Tleaf_over_sunlit_EASS
          par_EASS       => clm3%g%l%c%p%pef%par_over_sunlit_EASS
          ra_EASS        => clm3%g%l%c%p%pps%ra_over_EASS
          rb_EASS        => clm3%g%l%c%p%pps%rb_over_EASS
          frac_rain_can  => clm3%g%l%c%p%pps%frac_rain_over_EASS
          frac_snow_can  => clm3%g%l%c%p%pps%frac_snow_over_EASS
          radsb_EASS     => clm3%g%l%c%p%pef%radsb_over_sunlit_EASS
          Gheat_EASS     => clm3%g%l%c%p%pps%Gheat_over_EASS
          Gevap_EASS     => clm3%g%l%c%p%pps%Gevap_over_EASS
          Gleaf_EASS     => clm3%g%l%c%p%pps%Gleaf_over_sunlit_EASS
          ci_EASS        => clm3%g%l%c%p%pps%ci_over_sunlit_EASS
          psn_EASS       => clm3%g%l%c%p%pcf%psn_over_sunlit_EASS
          cs_EASS        => clm3%g%l%c%p%pps%cs_over_sunlit_EASS
          rs_EASS        => clm3%g%l%c%p%pps%rs_over_sunlit_EASS
          Tleaf_new_EASS => clm3%g%l%c%p%pes%Tleaf_new_over_sunlit_EASS
          vcmx_EASS      => clm3%g%l%c%p%pps%vcmx_over_sunlit_EASS
          jmx_EASS       => clm3%g%l%c%p%pps%jmx_over_sunlit_EASS
          dTleaf_EASS    => clm3%g%l%c%p%pps%dTleaf_over_sunlit_EASS
     else if (phase=='over_shaded') then
          Tleaf_EASS     => clm3%g%l%c%p%pes%Tleaf_over_shaded_EASS
          par_EASS       => clm3%g%l%c%p%pef%par_over_shaded_EASS
          ra_EASS        => clm3%g%l%c%p%pps%ra_over_EASS
          rb_EASS        => clm3%g%l%c%p%pps%rb_over_EASS
          frac_rain_can  => clm3%g%l%c%p%pps%frac_rain_over_EASS
          frac_snow_can  => clm3%g%l%c%p%pps%frac_snow_over_EASS
          radsb_EASS     => clm3%g%l%c%p%pef%radsb_over_shaded_EASS
          Gheat_EASS     => clm3%g%l%c%p%pps%Gheat_over_EASS
          Gevap_EASS     => clm3%g%l%c%p%pps%Gevap_over_EASS
          Gleaf_EASS     => clm3%g%l%c%p%pps%Gleaf_over_shaded_EASS
          ci_EASS        => clm3%g%l%c%p%pps%ci_over_shaded_EASS
          psn_EASS       => clm3%g%l%c%p%pcf%psn_over_shaded_EASS
          cs_EASS        => clm3%g%l%c%p%pps%cs_over_shaded_EASS
          rs_EASS        => clm3%g%l%c%p%pps%rs_over_shaded_EASS
          Tleaf_new_EASS => clm3%g%l%c%p%pes%Tleaf_new_over_shaded_EASS
          vcmx_EASS      => clm3%g%l%c%p%pps%vcmx_over_shaded_EASS
          jmx_EASS       => clm3%g%l%c%p%pps%jmx_over_shaded_EASS
          dTleaf_EASS    => clm3%g%l%c%p%pps%dTleaf_over_shaded_EASS
     else if (phase=='under_sunlit') then
          Tleaf_EASS     => clm3%g%l%c%p%pes%Tleaf_under_sunlit_EASS
          par_EASS       => clm3%g%l%c%p%pef%par_under_sunlit_EASS
          ra_EASS        => clm3%g%l%c%p%pps%ra_under_EASS
          rb_EASS        => clm3%g%l%c%p%pps%rb_under_EASS
          frac_rain_can  => clm3%g%l%c%p%pps%frac_rain_under_EASS
          frac_snow_can  => clm3%g%l%c%p%pps%frac_snow_under_EASS
          radsb_EASS     => clm3%g%l%c%p%pef%radsb_under_sunlit_EASS
          Gheat_EASS     => clm3%g%l%c%p%pps%Gheat_under_EASS
          Gevap_EASS     => clm3%g%l%c%p%pps%Gevap_under_EASS
          Gleaf_EASS     => clm3%g%l%c%p%pps%Gleaf_under_sunlit_EASS
          ci_EASS        => clm3%g%l%c%p%pps%ci_under_sunlit_EASS
          psn_EASS       => clm3%g%l%c%p%pcf%psn_under_sunlit_EASS
          cs_EASS        => clm3%g%l%c%p%pps%cs_under_sunlit_EASS
          rs_EASS        => clm3%g%l%c%p%pps%rs_under_sunlit_EASS
          Tleaf_new_EASS => clm3%g%l%c%p%pes%Tleaf_new_under_sunlit_EASS
          vcmx_EASS      => clm3%g%l%c%p%pps%vcmx_under_sunlit_EASS
          jmx_EASS       => clm3%g%l%c%p%pps%jmx_under_sunlit_EASS
          dTleaf_EASS    => clm3%g%l%c%p%pps%dTleaf_under_sunlit_EASS
      else if (phase=='under_shaded') then
          Tleaf_EASS     => clm3%g%l%c%p%pes%Tleaf_under_shaded_EASS
          par_EASS       => clm3%g%l%c%p%pef%par_under_shaded_EASS
          ra_EASS        => clm3%g%l%c%p%pps%ra_under_EASS
          rb_EASS        => clm3%g%l%c%p%pps%rb_under_EASS
          frac_rain_can  => clm3%g%l%c%p%pps%frac_rain_under_EASS
          frac_snow_can  => clm3%g%l%c%p%pps%frac_snow_under_EASS
          radsb_EASS     => clm3%g%l%c%p%pef%radsb_under_shaded_EASS 
          Gheat_EASS     => clm3%g%l%c%p%pps%Gheat_under_EASS
          Gevap_EASS     => clm3%g%l%c%p%pps%Gevap_under_EASS
          Gleaf_EASS     => clm3%g%l%c%p%pps%Gleaf_under_shaded_EASS
          ci_EASS        => clm3%g%l%c%p%pps%ci_under_shaded_EASS
          psn_EASS       => clm3%g%l%c%p%pcf%psn_under_shaded_EASS
          cs_EASS        => clm3%g%l%c%p%pps%cs_under_shaded_EASS
          rs_EASS        => clm3%g%l%c%p%pps%rs_under_shaded_EASS
          Tleaf_new_EASS => clm3%g%l%c%p%pes%Tleaf_new_under_shaded_EASS
          vcmx_EASS      => clm3%g%l%c%p%pps%vcmx_under_shaded_EASS
          jmx_EASS       => clm3%g%l%c%p%pps%jmx_under_shaded_EASS
          dTleaf_EASS    => clm3%g%l%c%p%pps%dTleaf_under_shaded_EASS
      end if

        !-------------------------------------------------

        cf=forc_pbot(g)/(rgas*0.001_r8*Tleaf_EASS(p))*1.0e6_r8  ![umol pa J-1]=[umol m-3]

        if ( (par_EASS(p) .gt. 0._r8) .and. (coszen(c) .gt. 0.0_r8) ) then 
        ! Changed by Jing Chen, the old one is 'if ( coszen(c) .gt. 0.0_r8 ) then '

               ! to calculate qflx_trans_EASS
               call PhotoSynth_EASS(p, g, c, Tleaf_EASS(p), par_EASS(p), &           ! input var.
                                    rb_EASS(p),dayl_factor, phase, ci_EASS(p),  &    ! input var.
                                    psn_EASS(p),cs_EASS(p), rs_EASS(p), &            ! output var.
                                    vcmx_EASS(p), jmx_EASS(p) )                      ! var. for debug
         else
               ci_EASS(p) =0.7_r8 * forc_pco2(g) * c3psn(ivt(p)) + 0.4_r8 * forc_pco2(g) * (1._r8-c3psn(ivt(p))) 
               psn_EASS(p)=0.0_r8 
               cs_EASS(p) =forc_pco2(g)
               rs_EASS(p) =min(rsmax0, 1._r8/bp*cf) 
               vcmx_EASS(p)=0.0_r8
               jmx_EASS(p)=0.0_r8
         end if

         ! difference between saturated and partical water vapor pressure at the reference height(de_EASS) 
         de_EASS=es_EASS-espar_EASS          ! [pa]
         de_EASS=de_EASS/1000.0_r8           ! unit change: [pa]->[kpa] 

         slope_EASS= 4098.0_r8 * (es_EASS/1000.0_r8 )/((forc_t(c)-tfrz+237.3_r8)**2.0_r8)  !slope of vapor pressure curve[kpa K-1] 

         Gleaf_EASS(p)=1.0_r8 /(ra_EASS(p)+rb_EASS(p)+rs_EASS(p))    ! [m s-1]
         p_star=(Gleaf_EASS(p)+Gevap_EASS(p)*(frac_rain_can(p)+frac_snow_can(p)))/psy  ! [m s-1] /[kPa K-1]

         dTleaf_EASS(p)=(radsb_EASS(p)-de_EASS*forc_rho(g)*cpair_EASS*p_star)/(forc_rho(g)*cpair_EASS*(Gheat_EASS(p)+slope_EASS*p_star)) 
        !      =   [W m-2] - Kpa   * kg m-3 *J kg-1 K-1*m s-1 kPa-1 K/(kg m-3 *J kg-1 K-1*m s-1)
        !      =   (J s-1 m-2 -   m-2 *J * s-1  /( m-2 *J K-1*s-1)
        !      = J s-1 m-2/( m-2 *J K-1*s-1)=K         
         Tleaf_new_EASS(p)=forc_t(c)+dTleaf_EASS(p)

    !#   Tleaf_new_EASS(p)=Tleaf_EASS(p)+dTleaf_EASS(p)
         Tleaf_new_EASS(p)=max(forc_t(c)-3.0_r8, Tleaf_new_EASS(p))
         Tleaf_new_EASS(p)=min(forc_t(c)+5.0_r8, Tleaf_new_EASS(p))

!#      end do  

   end subroutine LeafTemper_EASS 

end module LeafTemper_EASSMod
