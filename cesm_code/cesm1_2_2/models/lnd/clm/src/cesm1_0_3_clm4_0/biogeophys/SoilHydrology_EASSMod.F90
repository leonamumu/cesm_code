module SoilHydrology_EASSMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SoilHydrology_EASSMod
!
! !DESCRIPTION:
! Calculate soil hydrology
!
! ! USES:
   use clm_varctl,    only : iulog 
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRunoff_EASS  ! Calculate surface runoff
  public :: Infiltration_EASS   ! Calculate infiltration into surface soil layer
  public :: SoilWater_EASS      ! Calculate soil hydrology
  public :: Drainage_EASS       ! Calculate subsurface drainage
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 04/25/07 Keith Oleson: CLM3.5 hydrology
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRunoff_EASS
!
! !INTERFACE:
  subroutine SurfaceRunoff_EASS (lbc, ubc, lbp, ubp, num_hydrologyc, filter_hydrologyc, &
                            num_urbanc, filter_urbanc, vol_liq_EASS, icefrac_EASS)
!
! !DESCRIPTION:
! Calculate surface runoff
!
! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clmtype
    use clm_varcon      , only : denice, denh2o, wimp, pondmx_urban, &
                                 icol_roof, icol_sunwall, icol_shadewall, &
                                 icol_road_imperv, icol_road_perv
    use clm_varpar      , only : nlevsoi, maxpatch_pft
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                     ! column bounds
    integer , intent(in)  :: lbp, ubp                     ! pft bounds   
    integer , intent(in)  :: num_hydrologyc               ! number of column soil points in column filter
    integer , intent(in)  :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
    integer , intent(in)  :: num_urbanc                   ! number of column urban points in column filter
    integer , intent(in)  :: filter_urbanc(ubc-lbc+1)     ! column filter for urban points
    real(r8), intent(out) :: vol_liq_EASS(lbc:ubc,1:nlevsoi)   ! partial volume of liquid water in layer
    real(r8), intent(out) :: icefrac_EASS(lbc:ubc,1:nlevsoi)   ! fraction of ice in layer (-)
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/26/02, Peter Thornton: Migrated to new data structures.
! 4/26/05, David Lawrence: Made surface runoff for dry soils a function
!   of rooting fraction in top three soil layers.
! 04/25/07  Keith Oleson: Completely new routine for CLM3.5 hydrology
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: cgridcell(:)           ! gridcell index for each column
    integer , pointer :: ctype(:)               ! column type index
    real(r8), pointer :: qflx_top_soil_EASS(:)  ! net water input into soil from top (mm/s)
    real(r8), pointer :: watsat(:,:)            ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: hkdepth(:)             ! decay factor (m)
    real(r8), pointer :: zwt_EASS(:)            ! water table depth (m)
!#urban    real(r8), pointer :: fcov(:)           !fractional impermeable area
!#urban    real(r8), pointer :: fsat(:)           !fractional area with water table at surface
    real(r8), pointer :: dz_EASS(:,:)           ! layer depth (m)
    real(r8), pointer :: h2osoi_ice_EASS(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq_EASS(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: wtfact(:)              ! maximum saturated fraction for a gridcell
    real(r8), pointer :: hksat(:,:)             ! hydraulic conductivity at saturation (mm H2O /s)
    real(r8), pointer :: bsw(:,:)               ! Clapp and Hornberger "b"
    real(r8), pointer :: sucsat(:,:)            ! minimum soil suction (mm)
    integer , pointer :: snl_EASS(:)            ! minus number of snow layers
    real(r8), pointer :: qflx_evap_grnd_EASS(:) ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: zi_EASS(:,:)           ! interface level below a "z" level (m)
    real(r8), pointer :: qflx_surf(:)           ! surface runoff (mm H2O /s)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_surf_EASS(:)      ! surface runoff (mm H2O /s)
    real(r8), pointer :: eff_porosity_EASS(:,:) ! effective porosity = porosity - vol_ice
    real(r8), pointer :: fracice_EASS(:,:)      ! fractional impermeability (-)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,j,fc,g                   ! indices
    real(r8) :: dtime                      ! land model time step (sec)
    real(r8) :: xs(lbc:ubc)                ! excess soil water above urban ponding limit
    real(r8) :: vol_ice(lbc:ubc,1:nlevsoi) ! partial volume of ice lens in layer
    real(r8) :: fff(lbc:ubc)               ! decay factor (m-1)
    real(r8) :: s1                         ! variable to calculate qinmax
    real(r8) :: su                         ! variable to calculate qinmax
    real(r8) :: v                          ! variable to calculate qinmax
    real(r8) :: qinmax                     ! maximum infiltration capacity (mm/s)
    real(r8) :: fcov(lbc:ubc)                    ! fractional impermeable area
    real(r8) :: fsat(lbc:ubc)                    ! fractional area with water table at surface

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    ctype              => clm3%g%l%c%itype
    qflx_top_soil_EASS => clm3%g%l%c%cwf%qflx_top_soil_EASS
    qflx_surf_EASS     => clm3%g%l%c%cwf%qflx_surf_EASS
    watsat             => clm3%g%l%c%cps%watsat
    hkdepth            => clm3%g%l%c%cps%hkdepth
    dz_EASS            => clm3%g%l%c%cps%dz_EASS
    h2osoi_ice_EASS    => clm3%g%l%c%cws%h2osoi_ice_EASS
    h2osoi_liq_EASS    => clm3%g%l%c%cws%h2osoi_liq_EASS
!#urban    fcov          => clm3%g%l%c%cws%fcov
!#urban    fsat          => clm3%g%l%c%cws%fsat
    eff_porosity_EASS  => clm3%g%l%c%cps%eff_porosity_EASS
    wtfact             => clm3%g%l%c%cps%wtfact
    zwt_EASS           => clm3%g%l%c%cws%zwt_EASS
    fracice_EASS       => clm3%g%l%c%cps%fracice_EASS
    hksat              => clm3%g%l%c%cps%hksat
    bsw                => clm3%g%l%c%cps%bsw
    sucsat             => clm3%g%l%c%cps%sucsat
    snl_EASS           => clm3%g%l%c%cps%snl_EASS
    qflx_evap_grnd_EASS=> clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd_EASS
    zi_EASS            => clm3%g%l%c%cps%zi_EASS
    qflx_surf          => clm3%g%l%c%cwf%qflx_surf

    ! Get time step

    dtime = get_step_size()

    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Porosity of soil, partial volume of ice and liquid, fraction of ice in each layer,
          ! fractional impermeability
   
          vol_ice(c,j)      = min(watsat(c,j), h2osoi_ice_EASS(c,j)/(dz_EASS(c,j)*denice))
          eff_porosity_EASS(c,j) = max(0.01_r8,watsat(c,j)-vol_ice(c,j))
          vol_liq_EASS(c,j) = min(eff_porosity_EASS(c,j), h2osoi_liq_EASS(c,j)/(dz_EASS(c,j)*denh2o))
          icefrac_EASS(c,j) = min(1._r8,h2osoi_ice_EASS(c,j)/(h2osoi_ice_EASS(c,j)+h2osoi_liq_EASS(c,j)))
          fracice_EASS(c,j) = max(0._r8,exp(-3._r8*(1._r8-icefrac_EASS(c,j)))- exp(-3._r8))/(1.0_r8-exp(-3._r8))

       end do
    end do

    ! Saturated fraction

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       fff(c) = 0.5_r8
       fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt_EASS(c))
       fcov(c) = (1._r8 - fracice_EASS(c,1)) * fsat(c) + fracice_EASS(c,1)
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       ! Maximum infiltration capacity
       s1        = max(0.01_r8,vol_liq_EASS(c,1)/max(wimp,eff_porosity_EASS(c,1)))
       su        = max(0._r8,(s1-fcov(c)) / (max(0.01_r8,1._r8-fcov(c))))
       v         = -bsw(c,1)*sucsat(c,1)/(0.5_r8*dz_EASS(c,1)*1000._r8)
       qinmax    = (1._r8+v*(su-1._r8))*hksat(c,1)

       ! Surface runoff
       qflx_surf_EASS(c) =  fcov(c) * qflx_top_soil_EASS(c) + &
                       (1._r8-fcov(c)) * max(0._r8, qflx_top_soil_EASS(c)-qinmax)
    end do

    ! Determine water in excess of ponding limit for urban roof and impervious road.
    ! Excess goes to surface runoff. No surface runoff for sunwall and shadewall.

    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
           qflx_surf_EASS(c) = qflx_surf(c)
!#       if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then
!#
!#          ! If there are snow layers then all qflx_top_soil goes to surface runoff
!#          if (snl_EASS(c) < 0) then
!#             qflx_surf(c) = max(0._r8,qflx_top_soil(c))
!#          else
!#             xs(c) = max(0._r8, &
!#                         h2osoi_liq(c,1)/dtime + qflx_top_soil(c) - qflx_evap_grnd(c) - &
!#                         pondmx_urban/dtime)
!#             if (xs(c) > 0.) then
!#                h2osoi_liq(c,1) = pondmx_urban
!#             else
!#                h2osoi_liq(c,1) = max(0._r8,h2osoi_liq(c,1)+ &
!#                                     (qflx_top_soil(c)-qflx_evap_grnd(c))*dtime)
!#             end if
!#             qflx_surf(c) = xs(c)
!#          end if
!#       else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall) then
!#         qflx_surf(c) = 0._r8
!#       end if
    end do

  end subroutine SurfaceRunoff_EASS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Infiltration_EASS
!
! !INTERFACE:
  subroutine Infiltration_EASS(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                          num_urbanc, filter_urbanc)
!
! !DESCRIPTION:
! Calculate infiltration into surface soil layer (minus the evaporation)
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : icol_roof, icol_road_imperv, icol_sunwall, icol_shadewall, &
                             icol_road_perv
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                     ! column bounds
    integer, intent(in) :: num_hydrologyc               ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
    integer, intent(in) :: num_urbanc                   ! number of column urban points in column filter
    integer, intent(in) :: filter_urbanc(ubc-lbc+1)     ! column filter for urban points
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/27/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: ctype(:)               ! column type index
    integer , pointer :: snl_EASS(:)            ! minus number of snow layers
    real(r8), pointer :: qflx_top_soil_EASS(:)  ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_surf_EASS(:)      ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_evap_grnd_EASS(:) ! ground surface evaporation rate (mm H2O/s) [+]
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_infl_EASS(:)      !infiltration (mm H2O /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, fc    !indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    ctype               => clm3%g%l%c%itype
    snl_EASS            => clm3%g%l%c%cps%snl_EASS
    qflx_top_soil_EASS  => clm3%g%l%c%cwf%qflx_top_soil_EASS
    qflx_surf_EASS      => clm3%g%l%c%cwf%qflx_surf_EASS
    qflx_infl_EASS      => clm3%g%l%c%cwf%qflx_infl_EASS
    qflx_evap_grnd_EASS => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd_EASS

    ! Infiltration into surface soil layer (minus the evaporation)

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (snl_EASS(c) >= 0) then
          qflx_infl_EASS(c) = qflx_top_soil_EASS(c) - qflx_surf_EASS(c) - qflx_evap_grnd_EASS(c)
       else
          qflx_infl_EASS(c) = qflx_top_soil_EASS(c) - qflx_surf_EASS(c)
       end if
    end do

    ! No infiltration for impervious urban surfaces

    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       if (ctype(c) /= icol_road_perv) then
          qflx_infl_EASS(c) = 0._r8
       end if
    end do
    
  end subroutine Infiltration_EASS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilWater_EASS
!
! !INTERFACE:
  subroutine SoilWater_EASS(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                       num_urbanc, filter_urbanc, &
                       vol_liq_EASS, dwat_EASS, hk_EASS, dhkdw_EASS)
!
! !DESCRIPTION:
! Soil hydrology
! Soil moisture is predicted from a 10-layer model (as with soil
! temperature), in which the vertical soil moisture transport is governed
! by infiltration, runoff, gradient diffusion, gravity, and root
! extraction through canopy transpiration.  The net water applied to the
! surface layer is the snowmelt plus precipitation plus the throughfall
! of canopy dew minus surface runoff and evaporation.
! CLM3.5 uses a zero-flow bottom boundary condition.
!
! The vertical water flow in an unsaturated porous media is described by
! Darcy's law, and the hydraulic conductivity and the soil negative
! potential vary with soil water content and soil texture based on the work
! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
! integrated over the layer thickness, in which the time rate of change in
! water mass must equal the net flow across the bounding interface, plus the
! rate of internal source or sink. The terms of water flow across the layer
! interfaces are linearly expanded by using first-order Taylor expansion.
! The equations result in a tridiagonal system equation.
!
! Note: length units here are all millimeter
! (in temperature subroutine uses same soil layer
! structure required but lengths are m)
!
! Richards equation:
!
! d wat      d     d wat d psi
! ----- = - -- [ k(----- ----- - 1) ] + S
!   dt      dz       dz  d wat
!
! where: wat = volume of water per volume of soil (mm**3/mm**3)
! psi = soil matrix potential (mm)
! dt  = time step (s)
! z   = depth (mm)
! dz  = thickness (mm)
! qin = inflow at top (mm h2o /s)
! qout= outflow at bottom (mm h2o /s)
! s   = source/sink flux (mm h2o /s)
! k   = hydraulic conductivity (mm h2o /s)
!
!                       d qin                  d qin
! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!                       d wat(j-1)             d wat(j)
!                ==================|=================
!                                  < qin
!
!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
!
!                                  > qout
!                ==================|=================
!                        d qout               d qout
! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                        d wat(j)             d wat(j+1)
!
!
! Solution: linearize k and psi about d wat and use tridiagonal
! system of equations to solve for d wat,
! where for layer j
!
!
! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
! !USES:
    use shr_kind_mod  , only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon    , only : wimp, icol_roof, icol_road_imperv, denice, e_ice ! add by Shaobo Sun
    use clm_varpar    , only : nlevsoi, max_pft_per_col
    use clm_varctl    , only : iulog
    use shr_const_mod , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use TridiagonalMod, only : Tridiagonal
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                     ! column bounds
    integer , intent(in)  :: num_hydrologyc               ! number of column soil points in column filter
    integer , intent(in)  :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
    integer , intent(in)  :: num_urbanc                   ! number of column urban points in column filter
    integer , intent(in)  :: filter_urbanc(ubc-lbc+1)     ! column filter for urban points
    real(r8), intent(in)  :: vol_liq_EASS(lbc:ubc,1:nlevsoi)   ! soil water per unit volume [mm/mm]
    real(r8), intent(out) :: dwat_EASS(lbc:ubc,1:nlevsoi)      ! change of soil water [m3/m3]
    real(r8), intent(out) :: hk_EASS(lbc:ubc,1:nlevsoi)        ! hydraulic conductivity [mm h2o/s]
    real(r8), intent(out) :: dhkdw_EASS(lbc:ubc,1:nlevsoi)     ! d(hk)/d(vol_liq)
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/27/02, Peter Thornton: Migrated to new data structures. Includes
! treatment of multiple PFTs on a single soil column.
! 04/25/07 Keith Oleson: CLM3.5 hydrology
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: ctype(:)                  ! column type index
    integer , pointer :: npfts(:)                  ! column's number of pfts - ADD
    real(r8), pointer :: pwtcol(:)                 ! weight relative to column for each pft
    real(r8), pointer :: pwtgcell(:)               ! weight relative to gridcell for each pft
    real(r8), pointer :: z_EASS(:,:)               ! layer depth (m)
    real(r8), pointer :: dz_EASS(:,:)              ! layer thickness (m)
    real(r8), pointer :: smpmin(:)                 ! restriction for min of soil potential (mm)
    real(r8), pointer :: qflx_infl_EASS(:)         ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg_pft_EASS(:) ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg_col_EASS(:) ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: eff_porosity_EASS(:,:)    ! effective porosity = porosity - vol_ice
    real(r8), pointer :: watsat(:,:)               ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: hksat(:,:)                ! hydraulic conductivity at saturation (mm H2O /s)
    real(r8), pointer :: bsw(:,:)                  ! Clapp and Hornberger "b"
    real(r8), pointer :: sucsat(:,:)               ! minimum soil suction (mm)
    real(r8), pointer :: t_soisno_EASS(:,:)        ! soil temperature (Kelvin)
    real(r8), pointer :: rootr_pft(:,:)            ! effective fraction of roots in each soil layer
    integer , pointer :: pfti(:)                   ! beginning pft index for each column
    real(r8), pointer :: fracice_EASS(:,:)         ! fractional impermeability (-)
    real(r8), pointer :: h2osoi_vol_EASS(:,:)      ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: qcharge_EASS(:)           ! aquifer recharge rate (mm/s)
    real(r8), pointer :: hkdepth(:)                ! decay factor (m)
    real(r8), pointer :: zwt_EASS(:)               ! water table depth (m)
    real(r8), pointer :: zi_EASS(:,:)              ! interface level below a "z" level (m)
    real(r8), pointer :: h2osoi_ice_EASS(:,:)      ! ice water (kg/m2), Added by Shaobo Sun
!
! local pointers to original implicit inout arguments
!
    real(r8), pointer :: h2osoi_liq_EASS(:,:)      ! liquid water (kg/m2)
!
! local pointer s to original implicit out arguments
!
    real(r8), pointer :: rootr_col_EASS(:,:)       ! effective fraction of roots in each soil layer
!#    real(r8), pointer :: smp_l_EASS(:,:)           ! soil matrix potential [mm]
!#    real(r8), pointer :: hk_l_EASS(:,:)            ! hydraulic conductivity (mm/s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: p,c,fc,j                  ! do loop indices
    integer  :: jtop(lbc:ubc)             ! top level at each column
    real(r8) :: dtime                     ! land model time step (sec)
    real(r8) :: amx(lbc:ubc,1:nlevsoi+1)  ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(lbc:ubc,1:nlevsoi+1)  ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(lbc:ubc,1:nlevsoi+1)  ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(lbc:ubc,1:nlevsoi+1)  ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(lbc:ubc,1:nlevsoi+1)  ! layer depth [mm]
    real(r8) :: dzmm(lbc:ubc,1:nlevsoi+1) ! layer thickness [mm]
    real(r8) :: den                       ! used in calculating qin, qout
    real(r8) :: dqidw0(lbc:ubc,1:nlevsoi+1) ! d(qin)/d(vol_liq(i-1))
    real(r8) :: dqidw1(lbc:ubc,1:nlevsoi+1) ! d(qin)/d(vol_liq(i))
    real(r8) :: dqodw1(lbc:ubc,1:nlevsoi+1) ! d(qout)/d(vol_liq(i))
    real(r8) :: dqodw2(lbc:ubc,1:nlevsoi+1) ! d(qout)/d(vol_liq(i+1))
    real(r8) :: dsmpdw(lbc:ubc,1:nlevsoi+1) ! d(smp)/d(vol_liq)
    real(r8) :: num                         ! used in calculating qin, qout
    real(r8) :: qin(lbc:ubc,1:nlevsoi+1)    ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(lbc:ubc,1:nlevsoi+1)   ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: s_node                      ! soil wetness
    real(r8) :: s1                          ! "s" at interface of layer
    real(r8) :: s2                          ! k*s**(2b+2)
    real(r8) :: smp(lbc:ubc,1:nlevsoi)      ! soil matrix potential [mm]
    real(r8) :: sdamp                       ! extrapolates soiwat dependence of evaporation
    integer  :: pi                          ! pft index
    real(r8) :: temp(lbc:ubc)               ! accumulator for rootr weighting
    integer  :: jwt(lbc:ubc)                ! index of the soil layer right above the water table (-)
    real(r8) :: smp1,dsmpdw1,wh,wh_zwt,ka
    real(r8) :: dwat2(lbc:ubc,1:nlevsoi+1)
    real(r8) :: dzq                         ! used in calculating qin, qout (difference in equilbirium matric potential)
    real(r8) :: zimm(lbc:ubc,0:nlevsoi)     ! layer interface depth [mm]
    real(r8) :: zq(lbc:ubc,1:nlevsoi+1)     ! equilibrium matric potential for each layer [mm]
    real(r8) :: vol_eq(lbc:ubc,1:nlevsoi+1) ! equilibrium volumetric water content
    real(r8) :: tempi                       ! temp variable for calculating vol_eq
    real(r8) :: temp0                       ! temp variable for calculating vol_eq
    real(r8) :: voleq1                      ! temp variable for calculating vol_eq
    real(r8) :: zwtmm(lbc:ubc)              ! water table depth [mm]
    real(r8) :: icefrac(lbc:ubc,1:nlevsoi)  ! fraction of ice, Added by Shaobo Sun, 2015/5/24
    real(r8) :: vol_ice(lbc:ubc,1:nlevsoi)  ! Add by Shaobo sun

!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    qcharge_EASS           => clm3%g%l%c%cws%qcharge_EASS
    hkdepth                => clm3%g%l%c%cps%hkdepth
    zi_EASS                => clm3%g%l%c%cps%zi_EASS
    zwt_EASS               => clm3%g%l%c%cws%zwt_EASS
    ctype                  => clm3%g%l%c%itype
    npfts                  => clm3%g%l%c%npfts
    z_EASS                 => clm3%g%l%c%cps%z_EASS
    dz_EASS                => clm3%g%l%c%cps%dz_EASS
    smpmin                 => clm3%g%l%c%cps%smpmin
    watsat                 => clm3%g%l%c%cps%watsat
    hksat                  => clm3%g%l%c%cps%hksat
    bsw                    => clm3%g%l%c%cps%bsw
    sucsat                 => clm3%g%l%c%cps%sucsat
    eff_porosity_EASS      => clm3%g%l%c%cps%eff_porosity_EASS
    rootr_col_EASS         => clm3%g%l%c%cps%rootr_column_EASS
    t_soisno_EASS          => clm3%g%l%c%ces%t_soisno_EASS
    h2osoi_liq_EASS        => clm3%g%l%c%cws%h2osoi_liq_EASS
    h2osoi_vol_EASS        => clm3%g%l%c%cws%h2osoi_vol_EASS
    qflx_infl_EASS         => clm3%g%l%c%cwf%qflx_infl_EASS
    fracice_EASS           => clm3%g%l%c%cps%fracice_EASS
    qflx_tran_veg_col_EASS => clm3%g%l%c%cwf%pwf_a%qflx_tran_veg_EASS
    pfti                   => clm3%g%l%c%pfti
    h2osoi_ice_EASS        => clm3%g%l%c%cws%h2osoi_ice_EASS   ! Add by Shaobo Sun

!#    smp_l_EASS             => clm3%g%l%c%cws%smp_l_EASS
!#    hk_l_EASS              => clm3%g%l%c%cws%hk_l_EASS

    ! Assign local pointers to derived type members (pft-level)

    qflx_tran_veg_pft_EASS => clm3%g%l%c%p%pwf%qflx_tran_veg_EASS
    rootr_pft              => clm3%g%l%c%p%pps%rootr
    pwtcol                 => clm3%g%l%c%p%wtcol
    pwtgcell               => clm3%g%l%c%p%wtgcell

    ! Get time step

    dtime = get_step_size()

    ! Because the depths in this routine are in mm, use local
    ! variable arrays instead of pointers
    ! commented by Shaobo Sun, 2015/05/28 
    ! nlevsoi is 10 active soil layers    
    ! num_hydrologyc is the number of the computing grids
    ! filted arrays
    ! 1.e3 = 1000, covert m to mm 

    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          zmm(c,j)  = z_EASS(c,j)*1.e3_r8
          dzmm(c,j) = dz_EASS(c,j)*1.e3_r8
          zimm(c,j) = zi_EASS(c,j)*1.e3_r8
          ! calculate icefrac, added by Shaobo Sun
          vol_ice(c,j) = min(watsat(c,j), h2osoi_ice_EASS(c,j)/(dz_EASS(c,j)*denice))
          icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
       end do
    end do

    do fc = 1, num_hydrologyc 
       c = filter_hydrologyc(fc)
       zimm(c,0) = 0.0_r8
       zwtmm(c)  = zwt_EASS(c)*1.e3_r8
    end do

    ! First step is to calculate the column-level effective rooting
    ! fraction in each soil layer. This is done outside the usual
    ! PFT-to-column averaging routines because it is not a simple
    ! weighted average of the PFT level rootr arrays. Instead, the
    ! weighting depends on both the per-unit-area transpiration
    ! of the PFT and the PFTs area relative to all PFTs.

    temp(:) = 0._r8

    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          rootr_col_EASS(c,j) = 0._r8
       end do
    end do

    do pi = 1,max_pft_per_col
       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (pi <= npfts(c)) then
                p = pfti(c) + pi - 1
                if (pwtgcell(p)>0._r8) then
                   rootr_col_EASS(c,j) = rootr_col_EASS(c,j) + rootr_pft(p,j) &
                                       * qflx_tran_veg_pft_EASS(p) * pwtcol(p)
                end if
             end if
          end do
       end do
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (pi <= npfts(c)) then
             p = pfti(c) + pi - 1
             if (pwtgcell(p)>0._r8) then
                temp(c) = temp(c) + qflx_tran_veg_pft_EASS(p) * pwtcol(p)
             end if
          end if
       end do
    end do

    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (temp(c) /= 0._r8) then
             rootr_col_EASS(c,j) = rootr_col_EASS(c,j)/temp(c)
          end if
       end do
    end do

    !compute jwt index
    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt(c) = nlevsoi
       do j = 2,nlevsoi
          if(zwt_EASS(c) <= zi_EASS(c,j)) then
             jwt(c) = j-1
             exit
          end if
       enddo
    end do

    ! calculate the equilibrium water content based on the water table depth
            
    do j=1,nlevsoi 
       do fc=1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if ((zwtmm(c) .lt. zimm(c,j-1))) then   !fully saturated when wtd is less than the layer top
             vol_eq(c,j) = watsat(c,j)
            
          ! use the weighted average from the saturated part (depth > wtd) and the equilibrium solution for the
          ! rest of the layer

          else if ((zwtmm(c) .lt. zimm(c,j)) .and. (zwtmm(c) .gt. zimm(c,j-1))) then
             tempi = 1.0_r8
             temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
             voleq1 = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j-1))*(tempi-temp0)
             vol_eq(c,j) = (voleq1*(zwtmm(c)-zimm(c,j-1)) + watsat(c,j)*(zimm(c,j)-zwtmm(c)))/(zimm(c,j)-zimm(c,j-1))
             vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
             vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
          else
             tempi = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
             temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
             vol_eq(c,j) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zimm(c,j)-zimm(c,j-1))*(tempi-temp0)
             vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
             vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
          endif
          zq(c,j) = -sucsat(c,j)*(max(vol_eq(c,j)/watsat(c,j),0.01_r8))**(-bsw(c,j))
          zq(c,j) = max(smpmin(c), zq(c,j))
       end do
    end do

    ! If water table is below soil column calculate zq for the 11th layer
    j = nlevsoi
    do fc=1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(jwt(c) == nlevsoi) then 
          tempi = 1._r8
          temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
          vol_eq(c,j+1) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j))*(tempi-temp0)
          vol_eq(c,j+1) = max(vol_eq(c,j+1),0.0_r8)
          vol_eq(c,j+1) = min(watsat(c,j),vol_eq(c,j+1))
          zq(c,j+1) = -sucsat(c,j)*(max(vol_eq(c,j+1)/watsat(c,j),0.01_r8))**(-bsw(c,j))
          zq(c,j+1) = max(smpmin(c), zq(c,j+1))
       end if
    end do

    ! Hydraulic conductivity and soil matric potential and their derivatives

    sdamp = 0._r8
    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          s1 = 0.5_r8*(h2osoi_vol_EASS(c,j) + h2osoi_vol_EASS(c,min(nlevsoi, j+1))) / &
               (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
          s1 = min(1._r8, s1)
          s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)

          hk_EASS(c,j) = (10._r8**(-e_ice*(0.5_r8*(icefrac(c,j)+icefrac(c,min(nlevsoi,j+1))))))*s1*s2  ! add by Shaobo Sun
          
          dhkdw_EASS(c,j) = (10._r8**(-e_ice*(0.5_r8*(icefrac(c,j)+icefrac(c,min(nlevsoi,j+1))))))* &
                         (2._r8*bsw(c,j)+3._r8)*s2*0.5_r8/watsat(c,j)


!         hk_EASS(c,j) = (1._r8-0.5_r8*(fracice_EASS(c,j)+fracice_EASS(c,min(nlevsoi, j+1))))*s1*s2

!         dhkdw_EASS(c,j) = (1._r8-0.5_r8*(fracice_EASS(c,j)+fracice_EASS(c,min(nlevsoi, j+1))))* &
!                       (2._r8*bsw(c,j)+3._r8)*s2*0.5_r8/watsat(c,j)

          s_node = max(h2osoi_vol_EASS(c,j)/watsat(c,j), 0.01_r8)
          s_node = min(1.0_r8, s_node)

          smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
          smp(c,j) = max(smpmin(c), smp(c,j))

          dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))

!#          smp_l_EASS(c,j) = smp(c,j)
!#          hk_l_EASS(c,j) = hk_EASS(c,j)

       end do
    end do

    ! aquifer (11th) layer
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       zmm(c,nlevsoi+1) = 0.5*(1.e3_r8*zwt_EASS(c) + zmm(c,nlevsoi))
       if(jwt(c) < nlevsoi) then
         dzmm(c,nlevsoi+1) = dzmm(c,nlevsoi)
       else
         dzmm(c,nlevsoi+1) = (1.e3_r8*zwt_EASS(c) - zmm(c,nlevsoi))
       end if
    end do

    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1 (top)

    j = 1
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qin(c,j)    = qflx_infl_EASS(c)
       den    = (zmm(c,j+1)-zmm(c,j))
       dzq    = (zq(c,j+1)-zq(c,j))
       num    = (smp(c,j+1)-smp(c,j)) - dzq
       qout(c,j)   = -hk_EASS(c,j)*num/den
       dqodw1(c,j) = -(-hk_EASS(c,j)*dsmpdw(c,j)   + num*dhkdw_EASS(c,j))/den
       dqodw2(c,j) = -( hk_EASS(c,j)*dsmpdw(c,j+1) + num*dhkdw_EASS(c,j))/den
       rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col_EASS(c) * rootr_col_EASS(c,j)
       amx(c,j) =  0._r8
       bmx(c,j) =  dzmm(c,j)*(sdamp+1._r8/dtime) + dqodw1(c,j)
       cmx(c,j) =  dqodw2(c,j)
    end do

    ! Nodes j=2 to j=nlevsoi-1

    do j = 2, nlevsoi - 1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          den    = (zmm(c,j) - zmm(c,j-1))
          dzq    = (zq(c,j)-zq(c,j-1))
          num    = (smp(c,j)-smp(c,j-1)) - dzq
          qin(c,j)    = -hk_EASS(c,j-1)*num/den
          dqidw0(c,j) = -(-hk_EASS(c,j-1)*dsmpdw(c,j-1) + num*dhkdw_EASS(c,j-1))/den
          dqidw1(c,j) = -( hk_EASS(c,j-1)*dsmpdw(c,j)   + num*dhkdw_EASS(c,j-1))/den
          den         = (zmm(c,j+1)-zmm(c,j))
          dzq         = (zq(c,j+1)-zq(c,j))
          num         = (smp(c,j+1)-smp(c,j)) - dzq
          qout(c,j)   = -hk_EASS(c,j)*num/den
          dqodw1(c,j) = -(-hk_EASS(c,j)*dsmpdw(c,j)   + num*dhkdw_EASS(c,j))/den
          dqodw2(c,j) = -( hk_EASS(c,j)*dsmpdw(c,j+1) + num*dhkdw_EASS(c,j))/den
          rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_tran_veg_col_EASS(c)*rootr_col_EASS(c,j)
          amx(c,j)    = -dqidw0(c,j)
          bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
          cmx(c,j)    =  dqodw2(c,j)
       end do
    end do

    ! Node j=nlevsoi (bottom)

    j = nlevsoi
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(j > jwt(c)) then !water table is in soil column
         den    = (zmm(c,j) - zmm(c,j-1))
         dzq    = (zq(c,j)-zq(c,j-1))
         num    = (smp(c,j)-smp(c,j-1)) - dzq
         qin(c,j)    = -hk_EASS(c,j-1)*num/den
         dqidw0(c,j) = -(-hk_EASS(c,j-1)*dsmpdw(c,j-1) + num*dhkdw_EASS(c,j-1))/den
         dqidw1(c,j) = -( hk_EASS(c,j-1)*dsmpdw(c,j)   + num*dhkdw_EASS(c,j-1))/den
         qout(c,j)   =  0._r8
         dqodw1(c,j) =  0._r8
         rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_tran_veg_col_EASS(c)*rootr_col_EASS(c,j)
         amx(c,j)    = -dqidw0(c,j)
         bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
         cmx(c,j)    =  0._r8

         !scs: next set up aquifer layer; hydrologically inactive
         rmx(c,j+1) = 0._r8
         amx(c,j+1) = 0._r8
         bmx(c,j+1) = dzmm(c,j+1)/dtime
         cmx(c,j+1) = 0._r8

       else ! water table is below soil column

         !scs: compute aquifer soil moisture as average of layer 10 and saturation
         s_node = max(0.5*(1.0_r8+h2osoi_vol_EASS(c,j)/watsat(c,j)), 0.01_r8)
         s_node = min(1.0_r8, s_node)

         !scs: compute smp for aquifer layer
         smp1 = -sucsat(c,j)*s_node**(-bsw(c,j))
         smp1 = max(smpmin(c), smp1)

         !scs: compute dsmpdw for aquifer layer
         dsmpdw1 = -bsw(c,j)*smp1/(s_node*watsat(c,j))

         !scs: first set up bottom layer of soil column
         den    = (zmm(c,j) - zmm(c,j-1))
         dzq    = (zq(c,j)-zq(c,j-1))
         num    = (smp(c,j)-smp(c,j-1)) - dzq
         qin(c,j)    = -hk_EASS(c,j-1)*num/den
         dqidw0(c,j) = -(-hk_EASS(c,j-1)*dsmpdw(c,j-1) + num*dhkdw_EASS(c,j-1))/den
         dqidw1(c,j) = -( hk_EASS(c,j-1)*dsmpdw(c,j)   + num*dhkdw_EASS(c,j-1))/den
         den    = (zmm(c,j+1)-zmm(c,j))
         dzq    = (zq(c,j+1)-zq(c,j))
         num    = (smp1-smp(c,j)) - dzq
         qout(c,j)   = -hk_EASS(c,j)*num/den
         dqodw1(c,j) = -(-hk_EASS(c,j)*dsmpdw(c,j)   + num*dhkdw_EASS(c,j))/den
         dqodw2(c,j) = -( hk_EASS(c,j)*dsmpdw1 + num*dhkdw_EASS(c,j))/den

         rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col_EASS(c)*rootr_col_EASS(c,j)
         amx(c,j) = -dqidw0(c,j)
         bmx(c,j) =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
         cmx(c,j) =  dqodw2(c,j)

         !scs: next set up aquifer layer; den/num unchanged, qin=qout
         qin(c,j+1)    = qout(c,j)
         dqidw0(c,j+1) = -(-hk_EASS(c,j)*dsmpdw(c,j) + num*dhkdw_EASS(c,j))/den
         dqidw1(c,j+1) = -( hk_EASS(c,j)*dsmpdw1   + num*dhkdw_EASS(c,j))/den
         qout(c,j+1)   =  0._r8  ! zero-flow bottom boundary condition
         dqodw1(c,j+1) =  0._r8  ! zero-flow bottom boundary condition
         rmx(c,j+1) =  qin(c,j+1) - qout(c,j+1)
         amx(c,j+1) = -dqidw0(c,j+1)
         bmx(c,j+1) =  dzmm(c,j+1)/dtime - dqidw1(c,j+1) + dqodw1(c,j+1)
         cmx(c,j+1) =  0._r8

       endif
    end do

    ! Solve for dwat

    jtop(:) = 1
    call Tridiagonal(lbc, ubc, 1, nlevsoi+1, jtop, num_hydrologyc, filter_hydrologyc, &
                     amx, bmx, cmx, rmx, dwat2 )
    !scs: set dwat
    do fc = 1,num_hydrologyc
       c = filter_hydrologyc(fc)
       do j = 1, nlevsoi
          dwat_EASS(c,j)=dwat2(c,j)
       end do
    end do

    ! Renew the mass of liquid water
    !scs: also compute qcharge from dwat in aquifer layer
    !scs: update in drainage for case jwt < nlevsoi

    do fc = 1,num_hydrologyc
       c = filter_hydrologyc(fc)
       do j = 1, nlevsoi
          h2osoi_liq_EASS(c,j) = h2osoi_liq_EASS(c,j) + dwat2(c,j)*dzmm(c,j)
       end do

       !scs: calculate qcharge for case jwt < nlevsoi
       if(jwt(c) < nlevsoi) then
          wh_zwt = 0._r8   !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat

          s_node = max(h2osoi_vol_EASS(c,jwt(c))/watsat(c,jwt(c)), 0.01_r8)
          s_node = min(1.0_r8, s_node)

          !scs: use average moisture between water table and layer jwt
          s1 = 0.5_r8*(1.0+s_node)
          s1 = min(1._r8, s1)

          !scs: this is the expression for unsaturated hk
          ka = hksat(c,jwt(c))*s1**(2._r8*bsw(c,jwt(c))+3._r8)

          ! Recharge rate qcharge to groundwater (positive to aquifer)
          smp1 = -sucsat(c,jwt(c))*s_node**(-bsw(c,jwt(c)))
          smp1 = max(smpmin(c), smp(c,jwt(c)))
          wh   = smp1 - zq(c,jwt(c))
          qcharge_EASS(c) = -ka * (wh_zwt-wh)  /((zwt_EASS(c)-z_EASS(c,jwt(c)))*1000._r8)

          ! To limit qcharge  (for the first several timesteps)
          qcharge_EASS(c) = max(-10.0_r8/dtime,qcharge_EASS(c))
          qcharge_EASS(c) = min( 10.0_r8/dtime,qcharge_EASS(c))
       else
          !scs: if water table is below soil column, compute qcharge from dwat2(11)
          qcharge_EASS(c) = dwat2(c,nlevsoi+1)*dzmm(c,nlevsoi+1)/dtime
       endif
    end do

  end subroutine SoilWater_EASS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Drainage_EASS
!
! !INTERFACE:
  subroutine Drainage_EASS(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                      num_urbanc, filter_urbanc, vol_liq_EASS, hk_EASS, &
                      icefrac_EASS)
!
! !DESCRIPTION:
! Calculate subsurface drainage
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use clm_time_manager, only : get_step_size
    use clm_varcon  , only : pondmx, tfrz, icol_roof, icol_road_imperv, icol_road_perv, watmin
    use clm_varpar  , only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                     ! column bounds
    integer , intent(in) :: num_hydrologyc               ! number of column soil points in column filter
    integer , intent(in) :: num_urbanc                   ! number of column urban points in column filter
    integer , intent(in) :: filter_urbanc(ubc-lbc+1)     ! column filter for urban points
    integer , intent(in) :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
    real(r8), intent(in) :: vol_liq_EASS(lbc:ubc,1:nlevsoi)   ! partial volume of liquid water in layer
    real(r8), intent(in) :: hk_EASS(lbc:ubc,1:nlevsoi)        ! hydraulic conductivity (mm h2o/s)
    real(r8), intent(in) :: icefrac_EASS(lbc:ubc,1:nlevsoi)   ! fraction of ice in layer
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 4/26/05, Peter Thornton and David Lawrence: Turned off drainage from
! middle soil layers for both wet and dry fractions.
! 04/25/07  Keith Oleson: Completely new routine for CLM3.5 hydrology
! 27 February 2008: Keith Oleson; Saturation excess modification
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: ctype(:)               !column type index
    integer , pointer :: snl_EASS(:)            !number of snow layers
    real(r8), pointer :: qflx_snwcp_liq_EASS(:) !excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd_EASS(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow_EASS(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_EASS(:)  !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: dz_EASS(:,:)           !layer depth (m)
    real(r8), pointer :: bsw(:,:)               !Clapp and Hornberger "b"
    real(r8), pointer :: eff_porosity_EASS(:,:) !effective porosity = porosity - vol_ice
    real(r8), pointer :: t_soisno_EASS(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: hksat(:,:)             !hydraulic conductivity at saturation (mm H2O /s)
    real(r8), pointer :: sucsat(:,:)            !minimum soil suction (mm)
    real(r8), pointer :: z_EASS(:,:)            !layer depth (m)
    real(r8), pointer :: zi_EASS(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: watsat(:,:)            !volumetric soil water at saturation (porosity)
    real(r8), pointer :: hkdepth(:)             !decay factor (m)
    real(r8), pointer :: zwt_EASS(:)            !water table depth (m)
    real(r8), pointer :: wa_EASS(:)             !water in the unconfined aquifer (mm)
    real(r8), pointer :: wt_EASS(:)             !total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: qcharge_EASS(:)        !aquifer recharge rate (mm/s)
    !------------------------
    ! added by Jing Chen
    real(r8), pointer :: qflx_qrgwl(:)          !qflx_surf at glaciers,      wetlands, lakes (mm H2O /s)
    !---------------------
!
! local pointers to original implicit inout arguments
!
    real(r8), pointer :: h2osoi_ice_EASS(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq_EASS(:,:)   !liquid water (kg/m2)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_drain_EASS(:)     !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_irrig(:)          !irrigation flux (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl_EASS(:)     !qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
    real(r8), pointer :: eflx_impsoil_EASS(:)   !implicit evaporation for soil temperature equation
    real(r8), pointer :: qflx_rsub_sat_EASS(:)  !soil saturation excess [mm h2o/s]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,j,fc,i                 !indices
    real(r8) :: dtime                    !land model time step (sec)
    real(r8) :: xs(lbc:ubc)              !water needed to bring soil moisture to watmin (mm)
    real(r8) :: dzmm(lbc:ubc,1:nlevsoi)  !layer thickness (mm)
    integer  :: jwt(lbc:ubc)             !index of the soil layer right above the water table (-)
    real(r8) :: rsub_bot(lbc:ubc)        !subsurface runoff - bottom drainage (mm/s)
    real(r8) :: rsub_top(lbc:ubc)        !subsurface runoff - topographic control (mm/s)
    real(r8) :: fff(lbc:ubc)             !decay factor (m-1)
    real(r8) :: xsi(lbc:ubc)             !excess soil water above saturation at layer i (mm)
    real(r8) :: xsia(lbc:ubc)            !available pore space at layer i (mm)
    real(r8) :: xs1(lbc:ubc)             !excess soil water above saturation at layer 1 (mm)
    real(r8) :: smpfz(1:nlevsoi)         !matric potential of layer right above water table (mm)
    real(r8) :: wtsub                    !summation of hk*dzmm for layers below water table (mm**2/s)
    real(r8) :: rous                     !aquifer yield (-)
    real(r8) :: wh                       !smpfz(jwt)-z(jwt) (mm)
    real(r8) :: wh_zwt                   !water head at the water table depth (mm)
    real(r8) :: ws                       !summation of pore space of layers below water table (mm)
    real(r8) :: s_node                   !soil wetness (-)
    real(r8) :: dzsum                    !summation of dzmm of layers below water table (mm)
    real(r8) :: icefracsum               !summation of icefrac*dzmm of layers below water table (-)
    real(r8) :: fracice_rsub(lbc:ubc)    !fractional impermeability of soil layers (-)
    real(r8) :: ka                       !hydraulic conductivity of the aquifer (mm/s)
    real(r8) :: dza                      !fff*(zwt-z(jwt)) (-)
    real(r8) :: available_h2osoi_liq     !available soil liquid water in a layer
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    ctype               => clm3%g%l%c%itype
!   cgridcell      => clm3%g%l%c%gridcell

    snl_EASS           => clm3%g%l%c%cps%snl_EASS
    dz_EASS            => clm3%g%l%c%cps%dz
    bsw                => clm3%g%l%c%cps%bsw
    t_soisno_EASS      => clm3%g%l%c%ces%t_soisno_EASS
    hksat              => clm3%g%l%c%cps%hksat
    sucsat             => clm3%g%l%c%cps%sucsat
    z_EASS             => clm3%g%l%c%cps%z_EASS
    zi_EASS            => clm3%g%l%c%cps%zi_EASS
    watsat             => clm3%g%l%c%cps%watsat
    hkdepth            => clm3%g%l%c%cps%hkdepth
    zwt_EASS           => clm3%g%l%c%cws%zwt_EASS
    wa_EASS            => clm3%g%l%c%cws%wa_EASS
    wt_EASS            => clm3%g%l%c%cws%wt_EASS
    qcharge_EASS       => clm3%g%l%c%cws%qcharge_EASS
    eff_porosity_EASS  => clm3%g%l%c%cps%eff_porosity_EASS
    qflx_snwcp_liq_EASS=> clm3%g%l%c%cwf%pwf_a%qflx_snwcp_liq_EASS
    qflx_dew_grnd_EASS => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd_EASS
    qflx_dew_snow_EASS => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow_EASS
    qflx_sub_snow_EASS => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow_EASS
    qflx_drain_EASS    => clm3%g%l%c%cwf%qflx_drain_EASS
    qflx_irrig         => clm3%g%l%c%cwf%qflx_irrig
    qflx_qrgwl_EASS    => clm3%g%l%c%cwf%qflx_qrgwl_EASS
    qflx_rsub_sat_EASS => clm3%g%l%c%cwf%qflx_rsub_sat_EASS
    eflx_impsoil_EASS  => clm3%g%l%c%cef%eflx_impsoil_EASS
    h2osoi_liq_EASS    => clm3%g%l%c%cws%h2osoi_liq_EASS
    h2osoi_ice_EASS    => clm3%g%l%c%cws%h2osoi_ice_EASS
    qflx_qrgwl         => clm3%g%l%c%cwf%qflx_qrgwl

    ! Get time step

    dtime = get_step_size()

    ! Convert layer thicknesses from m to mm

    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          dzmm(c,j) = dz_EASS(c,j)*1.e3_r8
       end do
    end do

    ! Initial set

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qflx_drain_EASS(c)    = 0._r8 
       rsub_bot(c)           = 0._r8
       qflx_rsub_sat_EASS(c) = 0._r8
       rsub_top(c)           = 0._r8
       fracice_rsub(c)       = 0._r8
    end do

    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt(c) = nlevsoi
       do j = 2,nlevsoi
          if(zwt_EASS(c) <= zi_EASS(c,j)) then
             jwt(c) = j-1
             exit
          end if
       enddo
    end do

    ! Topographic runoff
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       fff(c)     = 1._r8/ hkdepth(c)
       dzsum      = 0._r8
       icefracsum = 0._r8
       do j = jwt(c), nlevsoi
          dzsum   = dzsum + dzmm(c,j)
          icefracsum = icefracsum + icefrac_EASS(c,j) * dzmm(c,j)
       end do
       fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum)))- exp(-3._r8))/(1.0_r8-exp(-3._r8))
       rsub_top(c)     = (1._r8 - fracice_rsub(c)) * 5.5e-3_r8 * exp(-fff(c)*zwt_EASS(c))
    end do

    rous = 0.2_r8

    ! Water table calculation

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       ! Water storage in aquifer + soil
       wt_EASS(c)  = wt_EASS(c) + (qcharge_EASS(c) - rsub_top(c)) * dtime

       if(jwt(c) == nlevsoi) then          ! water table is below the soil column
          wa_EASS(c)  = wa_EASS(c) + (qcharge_EASS(c) -rsub_top(c)) * dtime
          wt_EASS(c)  = wa_EASS(c)
          zwt_EASS(c)     = (zi_EASS(c,nlevsoi) + 25._r8) - wa_EASS(c)/1000._r8/rous
          h2osoi_liq_EASS(c,nlevsoi) = h2osoi_liq_EASS(c,nlevsoi) + max(0._r8,(wa_EASS(c)-5000._r8))
          wa_EASS(c)  = min(wa_EASS(c), 5000._r8)
       else                                ! water table within soil layers
          if (jwt(c) == nlevsoi-1) then    ! water table within bottom soil layer
             zwt_EASS(c) = zi_EASS(c,nlevsoi)- (wt_EASS(c)-rous*1000._r8*25._r8) /eff_porosity_EASS(c,nlevsoi)/1000._r8
          else                                   ! water table within soil layers 1-9
             ws = 0._r8   ! water used to fill soil air pores regardless of water content
             do j = jwt(c)+2,nlevsoi
               ws = ws + eff_porosity_EASS(c,j) * 1000._r8 * dz_EASS(c,j)
             enddo
             zwt_EASS(c) = zi_EASS(c,jwt(c)+1)-(wt_EASS(c)-rous*1000_r8*25._r8-ws) /eff_porosity_EASS(c,jwt(c)+1)/1000._r8
          endif

          wtsub = 0._r8
          do j = jwt(c)+1, nlevsoi
             wtsub = wtsub + hk_EASS(c,j)*dzmm(c,j)
          end do

          ! Remove subsurface runoff
          do j = jwt(c)+1, nlevsoi 
             h2osoi_liq_EASS(c,j) = h2osoi_liq_EASS(c,j) - rsub_top(c)*dtime*hk_EASS(c,j)*dzmm(c,j)/wtsub
          end do
       end if

       zwt_EASS(c) = max(0.05_r8,zwt_EASS(c))
       zwt_EASS(c) = min(80._r8,zwt_EASS(c))

    end do

    !  excessive water above saturation added to the above unsaturated layer like a bucket
    !  if column fully saturated, excess water goes to runoff

    do j = nlevsoi,2,-1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          xsi(c)    = max(h2osoi_liq_EASS(c,j)-eff_porosity_EASS(c,j)*dzmm(c,j),0._r8)
          h2osoi_liq_EASS(c,j)   = min(eff_porosity_EASS(c,j)*dzmm(c,j), h2osoi_liq_EASS(c,j))
          h2osoi_liq_EASS(c,j-1) = h2osoi_liq_EASS(c,j-1) + xsi(c)
       end do
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       xs1(c)       = max(max(h2osoi_liq_EASS(c,1),0._r8)-max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice_EASS(c,1))),0._r8)
       h2osoi_liq_EASS(c,1) = min(max(0._r8,pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice_EASS(c,1)), h2osoi_liq_EASS(c,1))
       qflx_rsub_sat_EASS(c)= xs1(c) / dtime
    end do

    ! Limit h2osoi_liq to be greater than or equal to watmin.
    ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
    ! If insufficient water in soil layers, get from aquifer water

    do j = 1, nlevsoi-1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (h2osoi_liq_EASS(c,j) < watmin) then
             xs(c) = watmin - h2osoi_liq_EASS(c,j)
          else
             xs(c) = 0._r8
          end if
          h2osoi_liq_EASS(c,j  ) = h2osoi_liq_EASS(c,j  ) + xs(c)
          h2osoi_liq_EASS(c,j+1) = h2osoi_liq_EASS(c,j+1) - xs(c)
       end do
    end do

! Get water for bottom layer from layers above if possible
    j = nlevsoi
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (h2osoi_liq_EASS(c,j) < watmin) then
          xs(c) = watmin-h2osoi_liq_EASS(c,j)
          searchforwater: do i = nlevsoi-1, 1, -1
             available_h2osoi_liq = max(h2osoi_liq_EASS(c,i)-watmin-xs(c),0._r8)
             if (available_h2osoi_liq .ge. xs(c)) then
               h2osoi_liq_EASS(c,j) = h2osoi_liq_EASS(c,j) + xs(c)
               h2osoi_liq_EASS(c,i) = h2osoi_liq_EASS(c,i) - xs(c)
               xs(c) = 0._r8
               exit searchforwater
             else
               h2osoi_liq_EASS(c,j) = h2osoi_liq_EASS(c,j) + available_h2osoi_liq
               h2osoi_liq_EASS(c,i) = h2osoi_liq_EASS(c,i) - available_h2osoi_liq
               xs(c) = xs(c) - available_h2osoi_liq
             end if
          end do searchforwater
       else
          xs(c) = 0._r8
       end if
! Needed in case there is no water to be found
       h2osoi_liq_EASS(c,j) = h2osoi_liq_EASS(c,j) + xs(c)
       wt_EASS(c) = wt_EASS(c) - xs(c)
! Instead of removing water from aquifer where it eventually
! shows up as excess drainage to the ocean, take it back out of 
! drainage
       rsub_top(c) = rsub_top(c) - xs(c)/dtime
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       ! Sub-surface runoff and drainage

       qflx_drain_EASS(c) = qflx_rsub_sat_EASS(c) + rsub_top(c)

       ! Set imbalance for snow capping

       qflx_qrgwl_EASS(c) = qflx_snwcp_liq_EASS(c)

       ! Implicit evaporation term is now zero

       eflx_impsoil_EASS(c) = 0._r8

       ! Renew the ice and liquid mass due to condensation

       if (snl_EASS(c)+1 >= 1) then
          h2osoi_liq_EASS(c,1) = h2osoi_liq_EASS(c,1) + qflx_dew_grnd_EASS(c) * dtime
          h2osoi_ice_EASS(c,1) = h2osoi_ice_EASS(c,1) + (qflx_dew_snow_EASS(c) * dtime)
          if (qflx_sub_snow_EASS(c)*dtime > h2osoi_ice_EASS(c,1)) then
             qflx_sub_snow_EASS(c) = h2osoi_ice_EASS(c,1)/dtime
             h2osoi_ice_EASS(c,1) = 0._r8
          else
             h2osoi_ice_EASS(c,1) = h2osoi_ice_EASS(c,1) - (qflx_sub_snow_EASS(c) * dtime)
          end if
       end if
    end do

    ! No drainage for urban columns (except for pervious road as computed above)

    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       if (ctype(c) /= icol_road_perv) then
         qflx_drain_EASS(c) = 0._r8
         qflx_irrig(c) = 0._r8
         ! This must be done for roofs and impervious road (walls will be zero)
         qflx_qrgwl_EASS(c) = qflx_qrgwl(c)
         eflx_impsoil_EASS(c) = 0._r8
       end if

       ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

       if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then
         if (snl_EASS(c)+1 >= 1) then
            h2osoi_liq_EASS(c,1) = h2osoi_liq_EASS(c,1) + qflx_dew_grnd_EASS(c) * dtime
            h2osoi_ice_EASS(c,1) = h2osoi_ice_EASS(c,1) + (qflx_dew_snow_EASS(c) * dtime)
            if (qflx_sub_snow_EASS(c)*dtime > h2osoi_ice_EASS(c,1)) then
               qflx_sub_snow_EASS(c) = h2osoi_ice_EASS(c,1)/dtime
               h2osoi_ice_EASS(c,1) = 0._r8
            else
               h2osoi_ice_EASS(c,1) = h2osoi_ice_EASS(c,1) - (qflx_sub_snow_EASS(c) * dtime)
            end if
         end if
       end if

    end do

  end subroutine Drainage_EASS

end module SoilHydrology_EASSMod
