module PenmanEqu_EASSMod
!--------------------------------------------------
!BOP
!! MODULE:PenmanEqu_EASSMod

!! DESCRIPTON:
! Calculates the leaf flux and ground flux based on Penman equ
!! USE: 
   use shr_kind_mod, only: r8 => shr_kind_r8
   

!! PUBLIC TYPES:
   implicit none
   save

!! PUBLIC MEMBER FUNCTIONS:
   public :: LeafFlux_EASS   
   public :: CanopyFlux_EASS 
   public :: GrndFlux_EASS   
   public :: LHvaporFunc_EASS

!! RVESION HISTORY:
! 5/21/2012, Jing Chen: created the module
!EOP
!--------------------------------------------------
contains
!-------------------------------------------------
    ! BOP
    !
    !!IROUTINE: LeafFlux_EASS
    !
    !!INTERFACE:
    subroutine LeafFlux_EASS(tak, rhoair, cpair_EASS, Tleaf_EASS, es_EASS, espar_EASS,&
                             G_EASS, frac_cov, qflx_leaf_EASS)

    !! DESCRIPTION:
    ! Calculates the leaf flux based on penman equ.
    !! REVISION HISTROY:
    !! USE:
         use clm_varcon     , only : tfrz

    !! ARGUMENTS:
         implicit none
         real(r8), intent(in) :: tak         ! temperature[K]
         real(r8), intent(in) :: rhoair         ! density of air at 0 degree Celsius[kg m-3]
         real(r8), intent(in) :: cpair_EASS     ! specific heat of dry air in EASS[J kg-1 k-1]
                                                ! depanding on saturated water vapor specific humidity
         real(r8), intent(in) :: Tleaf_EASS     ! leaf temperature
         real(r8), intent(in) :: es_EASS        ! saturated atmopsheric pressure [pa]
         real(r8), intent(in) :: espar_EASS     ! partical atmopsheric pressure [pa]
         real(r8), intent(in) :: G_EASS         ! conductance of sunlit and shaded leaves(overstory, understory)[m s-1], including
                                                ! (1)water conductance for transpiration
                                                ! (2)intercepted water conductance for water(liq.) evaporation and snow of leaves and
                                                ! (3)heat (transfer) conductance for sensible heat
         real(r8), intent(in) :: frac_cov       ! fraction of fall on the leaves
         real(r8), intent(out):: qflx_leaf_EASS ! vapor flux of leaves[w m-2]

         !!OTHER VARIABLES:
         real(r8), parameter :: psy=0.066_r8    ! psychrometer constant[kPa K-1]
                                                ! at a standard pressure of 101.kPa, the value is about 66 Pa/K at 0, and increase to 67 Pa/K at 20C
         real(r8) :: Bowen_const          ! constent in Bowen Ratio
         real(r8) :: de_EASS              ! unsaturated atmopsheric pressure [kpa]
         real(r8) :: slope_EASS           ! slope of vapor pressure curve [kpa K-1]
         !------------------------------------------------------
         ! constent of Bowen ratio
         Bowen_const= rhoair * cpair_EASS / psy      ! unit:1[kg m-3]*1[J kg-1 K-1]*1[kPa-1 oC]=1 J kpa-1 m-2

         ! difference between saturated and partical water vapor pressure at the reference height(de_EASS)
         de_EASS=es_EASS-espar_EASS      ! [pa]
         de_EASS=de_EASS/1000.0_r8       ! unit change: [pa]->[kpa]

         slope_EASS= 4098.0_r8 * (es_EASS/1000.0_r8)/((tak-tfrz+237.3_r8)**2.0_r8)  ! slope of vapor pressure curve[kpa k-1]

         ! leaf flux
         qflx_leaf_EASS=G_EASS*frac_cov*Bowen_const*(de_EASS+slope_EASS*(Tleaf_EASS-tak))
         ! unit: 1[m s-1]*1*1 [J kpa-1 m-3]*(1kpa+1[kpa oC-1]*1oC)=1J s-1 m-2=1w m-2

      end subroutine LeafFlux_EASS

    !-----------------------------------
    ! BOP
    !
    !!IROUTINE: CanopyFlux_EASS
    !
    !!INTERFACE:
     subroutine CanopyFlux_EASS(flux_sunlit,flux_shaded, h, lt_sunlit, lt_shaded, CanopyFlux)
    !
    !! DESCRIPTION
    !  calulate canopy flux using sunlit and shaded leaves
    !!USES:

    !! ARGUMENTS:
        implicit none
        real(r8),intent(in) :: flux_sunlit     ! flux of sunlit leaves, including canopy temperature, heat and water vapor[w m-2]
        real(r8),intent(in) :: flux_shaded     ! flux of shaded leaves, including canopy temperature, heat and water vapor[w m-2]
        real(r8),intent(in) :: h               ! indicate hvap or hsub[J kg-1]
        real(r8),intent(in) :: lt_sunlit       ! elai_sunlit+esai_sunlit
        real(r8),intent(in) :: lt_shaded       ! elai_shaded+esai_shaded
        real(r8),intent(out) :: CanopyFlux     ! elai_shaded+esai_shaded
    !
    !! CALLED FROM:
    ! subroutine main in this module
    !
    !! REVISTION HISTORY:
    !----------------------------------------------------
       CanopyFlux=(flux_sunlit*lt_sunlit+flux_shaded*lt_shaded)/h 
      ! unit: w m-2 * kg J-1 =kg w m-2 J-1=kg m-2 s-1=mmH2O s-1

     end subroutine CanopyFlux_EASS

    !---------------------------------------------------
    ! BOP
    !
    !!IROUTINE: GrndFlux_EASS
    !
    !!INTERFACE:
    subroutine GrndFlux_EASS(cpair_EASS, slope_grnd_EASS, radnet_grnd_EASS, &
                             const, de_grnd_EASS, ra_grnd_EASS, rs_grnd_EASS,&
                             h, frac_cov, qflx_grnd_EASS)

    !! DESCRIPTION:
    ! Calculates the ground flux based on penman equ.
    !! REVISION HISTROY:
    !! USE:

    !!ARGUMENTS:
        implicit none
        real(r8), intent(in) :: cpair_EASS        ! specific heat capacity of air[J kg-1 K-1]
        real(r8), intent(in) :: slope_grnd_EASS   ! [kpa K-1]
        real(r8), intent(in) :: radnet_grnd_EASS   ! [W m-2]
        real(r8), intent(in) :: const
        real(r8), intent(in) :: de_grnd_EASS      ! [kpa]
        real(r8), intent(in) :: ra_grnd_EASS      ! air resitance, ground[s m-1]
        real(r8), intent(in) :: rs_grnd_EASS      ! soil resistant for evaporation[s m-1]
        real(r8), intent(in) :: h                 ! J kg-1
        real(r8), intent(in) :: frac_cov          ! fraction of snow on the ground[-]
        real(r8), intent(out) :: qflx_grnd_EASS   ! vapor flux of leaves[mmH2O s-1]


    !!OTHER VARIABLES:
        real(r8), parameter :: forc_rho=1.292_r8 ! density of air at 0 degree Celsius[kg m-3]
        real(r8), parameter :: psy=0.066_r8      ! psychrometer constant[kPa k-1]
                                              ! at a standard pressure of 101.kPa, the value is about 66 Pa/K at 0, 
                                              ! and increase to 67 Pa/K at 20C
        real(r8) :: heat
        real(r8) :: sl

    !----------------------------------------------------------------
        ! leaf flux
        heat=slope_grnd_EASS*radnet_grnd_EASS*const+cpair_EASS*forc_rho*de_grnd_EASS/ra_grnd_EASS
        ! unit: kpa K-1* W m-2 * * +J kg-1 K-1 *kg m-3 *kpa* s-1 m= w  K-1 * m-2 *kpa*
        sl=slope_grnd_EASS+psy*(1.0_r8+rs_grnd_EASS/ra_grnd_EASS)      ! unit: kpa K-1

        qflx_grnd_EASS= 1.0_r8 / h * ( heat / sl ) * frac_cov    
        ! unit:  [mmH2O s-1](=kg J-1 * w m-2 = J s-1 kg J-1 m-2= kg m-2 s-1)

    end subroutine GrndFlux_EASS

   !---------------------------------------------------
   !BOP

   !!INTERFACE:
     real(r8) function LHvaporFunc_EASS(tak)

   !!DESRIPTION:
   ! calculate Latent heat of vaporiation[J kg-1] based on temperature[K]

   !! CALLED FROM:
   ! subroutine PhotoSynth_EASS in the module PhotoSynth_EASSMod , to calculate vcmax
   ! subroutine Canopy Fluxes , to calculate latent heat

   !!REVISION HISTORY:
   ! 12 January 2012: Jing Chen F90 Revision of ../eass_1031/EASS_region_1028/db_PHOTOSYNTHESIS.c

   !!USES:
       use clm_varcon,     only: tfrz, hvap

   !! ARGUMENTS:
      implicit none
      real (r8), intent(in) :: tak    ![K]

   !EOP
   !-------------------
      LHvaporFunc_EASS = hvap-2370.0_r8*(tak-tfrz)
!      LHvaporFunc_EASS = hvap-2260.0_r8*(tak-tfrz)   !add by Shaobo Sun, 2016/08/14

      ! add heat of fusion for melting ice
      !  if (tak .lt. 273.0) then
      !        LHvaporFunc_EASS = LHvaporFunc_EASS + 333.0
      !  end if
      return

   end function LHvaporFunc_EASS
  !------------------------------------

end module PenmanEqu_EASSMod
