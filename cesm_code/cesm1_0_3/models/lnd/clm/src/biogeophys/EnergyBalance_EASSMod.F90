module EnergyBalance_EASSMod
!--------------------------------------------------
!BOP
!! USE:
   use shr_kind_mod,    only: r8=>shr_kind_r8
   use clm_varctl,      only : iulog

!! PUBLIC TYPES:
   implicit none
   save

!! PUBLIC MEMBER FUNCTIONS:
   public :: EnergyBalance_EASS         ! calculate air resistnace and leaf boundry resistance

!! RVESION HISTORY

!EOP
!--------------------------------------------------
 contains
!-------------------------------------------------
!! BOPf
!! IROUTINE: EnergyBalance_EASS
!! INTERFACE

!#   subroutine EnergyBalance_EASS(fn, filterp, lbp, ubp ) 
   subroutine EnergyBalance_EASS(p, c, g, espar)

!! DESCIRPTION:


!! USES:
      use clmtype
      use clm_varcon,           only : sb   !Stefan-boltzmann constant[W m-2 K-4]= 5.67*1.0e-8
      use SurfaceRadiationMod,  only : SunlitLAI_EASS 
      use PenmanEqu_EASSMod,    only : CanopyFlux_EASS
      use clm_atmlnd,           only : clm_a2l

!! ARGUMENTS:
      implicit none
!#      integer, intent(in) :: fn                  ! size of pft filter
!#      integer, intent(in) :: filterp(fn)         ! pft filter
!#      integer, intent(in) :: lbp, ubp            ! pft bounds
      integer,intent(in) :: p                           ! pft index
      integer,intent(in) :: c                           ! column index
      integer,intent(in) :: g                           ! gridcell index
      real(r8), intent(in) :: espar           ! partiacl water vapro pressure[pa]

!! CALLED FROM:

!!LOCAL VARIABLES:
! local pointer to implicit in arguments
!#      integer , pointer :: pcolumn(:)              ! pft's column index
!#      integer , pointer :: pgridcell(:)            ! pft's gridcell index
      real(r8), pointer :: forc_t(:)               ! atmospheric temperature [Kelvin] at the pft level
      real(r8), pointer :: forc_solar(:)           ! solar (shortwave) radiation [W m-2]		  
      real(r8), pointer :: forc_lwrad(:)           ! downward infrared (longwave) radiation [W m-2]
      real(r8), pointer :: clumping_EASS(:)        ! clumping index in EASS
      real(r8), pointer :: elai_over_EASS(:)       ! elai of overstory
      real(r8), pointer :: elai_under_EASS(:)      ! elai of understroy
      real(r8), pointer :: esai_over_EASS(:)       ! esai of overstory
      real(r8), pointer :: esai_under_EASS(:)      ! esai of understory
      real(r8), pointer :: lt_over_sunlit_EASS(:)  ! elai+slai of sunlit leaves, overstory. To calculate leaf temperature.
      real(r8), pointer :: lt_over_shaded_EASS(:)  ! elai+slai of shaded leaves, overstory. To calculate leaf temperature.
      real(r8), pointer :: lt_under_sunlit_EASS(:) ! elai+slai of sunlit leaves, understory. To calculate leaf temperature.
      real(r8), pointer :: lt_under_shaded_EASS(:) ! elai+slai of shaded leaves, understory. To calculate leaf temperature.
      real(r8), pointer :: Tleaf_over_EASS(:)      ! leaf temperature of overstory[K]
      real(r8), pointer :: Tleaf_under_EASS(:)     ! leaf temperature of understory[K]
      real(r8), pointer :: Tleaf_over_sunlit_EASS(:)   ! leaf temperature of sunlit leaves, overstory [K]
      real(r8), pointer :: Tleaf_under_sunlit_EASS(:)  ! leaf temperature of sunlit leaves, understory [K]
      real(r8), pointer :: Tleaf_over_shaded_EASS(:)   ! leaf temperature of shaded leaves, overstory [K]
      real(r8), pointer :: Tleaf_under_shaded_EASS(:)  ! leaf temperature of shaded leaves, understory [K]
      real(r8), pointer :: forc_t_col(:)            ! atmospheric temperature [Kelvin] at the column level
     !real(r8), pointer :: Tgrnd_EASS(:)            ! ground temperature[K]
      real(r8), pointer :: coszen(:)                ! cosine of solar zenith angle
      real(r8), pointer :: dalb_over_EASS(:)        ! diff. between 1.0 and albedo of overstory
      real(r8), pointer :: dalb_under_EASS(:)       ! diff. between 1.0 and albedo of understory
      real(r8), pointer :: dalb_grnd_EASS(:)        ! diff. between 1.0 and albedo of ground
      real(r8), pointer :: t_grnd(:)                ! ground surface temperature [K]
      integer , pointer :: frac_veg_nosno(:)        ! fraction of vegetation not covered by snow (0 OR 1) [-]
      !------------------------------

      ! local pointer to implicit out arguments
      real(r8), pointer :: radsb_over_sunlit_EASS(:)     ! total solar radiation aborbed by sunlit leave, overstory[W m-2]
      real(r8), pointer :: radsb_over_shaded_EASS(:)     ! total solar radiation aborbed by shaded leave, overstory[W m-2]
      real(r8), pointer :: radsb_under_sunlit_EASS(:)    ! total solar radiation aborbed by sunlit leave, understory[W m-2]
      real(r8), pointer :: radsb_under_shaded_EASS(:)    ! total solar radiation aborbed by shaded leave, understory[W m-2]
      real(r8), pointer :: radnet_over_EASS(:)           ! net radiation of overstory [W m-2]
      real(r8), pointer :: radnet_under_EASS(:)          ! net radiation of understory [W m-2]
      real(r8), pointer :: radnet_grnd_EASS(:)           ! net radiation of ground [W m-2]
      real(r8), pointer :: radnet_EASS(:)                ! net radiation [W m-2](radnet_over_EASS+ radnet_under_EASS+radnet_grnd_EASS)
      real(r8), pointer :: lwnet_EASS(:)                 ! net longwave radiation [W m-2]
      real(r8), pointer :: lwrad_grnd_EASS(:)            ! longwave radiation of grnd
      real(r8), pointer :: swnet_EASS(:)                 ! net solar (shortwave) radiation [W m-2]
      real(r8), pointer :: swnet_over_EASS(:)            ! net solar (shortwave) radiation of overstory[W m-2]
      real(r8), pointer :: swnet_under_EASS(:)           ! net solar (shortwave) radiation of understory[W m-2]
      real(r8), pointer :: swnet_grnd_EASS(:)            ! net solar (shortwave) radiation of grnd[W m-2]
      real(r8), pointer :: par_over_sunlit_EASS(:)       ! average absorbed PAR for sunlit leaves, overstory[W m-2]
      real(r8), pointer :: par_over_shaded_EASS(:)       ! average absorbed PAR for shaded leaves, overstory[W m-2]
      real(r8), pointer :: par_under_sunlit_EASS(:)      ! average absorbed PAR for sunlit leaves, understory[W m-2]
      real(r8), pointer :: par_under_shaded_EASS(:)      ! average absorbed PAR for shaded leaves, understory[W m-2]
      !------------------------

      !! OTHER LOCAL VARIABLES:
!#      integer :: f                           ! filter index
!#      integer :: p                           ! pft index
!#      integer :: c                           ! column index
!#      integer :: g                           ! gridcell index
      real(r8) :: cosz                       ! cosine of solar zenith	  
      real(r8) :: cl                         ! clearness index, also referred as the cloudiness index
      real(r8) :: qita_over
      real(r8) :: qita_under
      real(r8) :: lt_over_EASS               ! elai_over+esai_over
      real(r8) :: lt_under_EASS              ! elai_under+esai_under
      real(r8) :: elai_over_under_EASS       ! elai_over+elai_under
      real(r8) :: flc_elai_over              ! fractional area of sunflecks on horizontal plane below the LAI, overstory
      real(r8) :: flc_elai_under             ! fractional area of sunflecks on horizontal plane below the LAI, understory(divieded by coszen)
      real(r8) :: flc_elai_over_under        ! fractional area of sunflecks on horizontal plane below the(LAI_over+LAI_under)(divieded by coszen)
      real(r8) :: flc_lt_over                ! fractional area of sunflecks on horizontal plane below the(LAI_over+SAI_over)(divieded by coszen)
      real(r8) :: flq_elai_over              ! fractional area of sunflecks on horizontal plane below the LAI, overstory(divieded by qita)
      real(r8) :: flq_elai_under             ! fractional area of sunflecks on horizontal plane below the LAI, understory(divieded by qita)
      real(r8) :: flq_lt_over                ! fractional area of sunflecks on horizontal plane below the(LAI_over+SAI_over)(divieded by qita)
      real(r8) :: flq_lt_under               ! fractional area of sunflecks on horizontal plane below the(LAI_under+SAI_under)(divieded by qita)
      real(r8) :: swdif_EASS                 ! diffused solar radiation [W m-2]
      real(r8) :: swdir_EASS                 ! direct solar radiation [W m-2]
      real(r8) :: epsilon                    ! [-]
      real(r8) :: lwrad_EASS                 ! longwave radiation in EASS[W m-2]
      real(r8) :: swdir_over_EASS            ! direct solar radiation(shortwave) of overstory
      real(r8) :: swdir_under_EASS           ! direct solar radiation(shortwave) of understory
      real(r8) :: swdir_grnd_EASS            ! direct solar radiation(shortwave) of grnd
      real(r8) :: pardir_over_EASS           ! direct solar radiation aborbed by overstory 
      real(r8) :: pardir_under_EASS          ! direct solar radiation aborbed by understory 
      real(r8) :: swdif_over_EASS            ! diffused solar radiation(shortwave) of overstory
      real(r8) :: swdif_under_EASS           ! diffused solar radiation(shortwave) of understory
      real(r8) :: swdif_grnd_EASS            ! diffused solar radiation(shortwave) of grnd
      real(r8) :: pardif_over_EASS           ! diffused solar radiation aborbed by overstory 
      real(r8) :: pardif_under_EASS          ! diffused solar radiation aborbed by understory 
      real(r8) :: temp_over                  ! temporary var. about overstory to calculate longwave radiation
      real(r8) :: temp_under                 ! temporary var. about understory to calculate longwave radiation
      real(r8) :: temp_grnd                  ! temporary var. about ground to calculate longwave radiation
      real(r8) :: lwrad_over_EASS            ! longwave radiation of overstory
      real(r8) :: lwrad_under_EASS           ! longwave radiation of understory
      real(r8) :: h                          ! temporary var.
      real(r8) :: Tgrnd_EASS                 ! ground temperature[K]
      real(r8), parameter :: em_over=0.98_r8    ! overstory emissivity
      real(r8), parameter :: em_under=0.98_r8   ! understory emissivity
      real(r8), parameter :: em_grnd=0.96_r8    ! ground emissivity
      character(len=20):: phase                 ! phase="LAI" or "energy"
      integer :: ib

      !-------------------------------------------------
      ! Assign local pointers to derived subtypes componsents(gridcell-level)
      forc_t                     => clm_a2l%forc_t
      forc_solar                 => clm_a2l%forc_solar
      forc_lwrad                 => clm_a2l%forc_lwrad

      ! Assign local pointers to derived subtypes componsents(column-level)
      forc_t_col                 => clm3%g%l%c%ces%forc_t
     !Tgrnd_EASS                 => clm3%g%l%c%ces%Tgrnd_EASS
      coszen                     => clm3%g%l%c%cps%coszen
      t_grnd                     => clm3%g%l%c%ces%t_grnd

      ! Assign local pointers to derived subtypes componsents(pft-level)
!#      pcolumn                    => clm3%g%l%c%p%column
!#      pgridcell                  => clm3%g%l%c%p%gridcell
      clumping_EASS              => clm3%g%l%c%p%pps%clumping_EASS
      elai_over_EASS             => clm3%g%l%c%p%pps%elai_over_EASS
      elai_under_EASS            => clm3%g%l%c%p%pps%elai_under_EASS
      esai_over_EASS             => clm3%g%l%c%p%pps%esai_over_EASS
      esai_under_EASS            => clm3%g%l%c%p%pps%esai_under_EASS
      lt_over_shaded_EASS        => clm3%g%l%c%p%pps%lt_over_shaded_EASS
      lt_over_sunlit_EASS        => clm3%g%l%c%p%pps%lt_over_sunlit_EASS
      lt_under_shaded_EASS       => clm3%g%l%c%p%pps%lt_under_shaded_EASS
      lt_under_sunlit_EASS       => clm3%g%l%c%p%pps%lt_under_sunlit_EASS
      Tleaf_over_EASS            => clm3%g%l%c%p%pes%Tleaf_over_EASS 
      Tleaf_under_EASS           => clm3%g%l%c%p%pes%Tleaf_under_EASS
      Tleaf_over_sunlit_EASS     => clm3%g%l%c%p%pes%Tleaf_over_sunlit_EASS
      Tleaf_under_sunlit_EASS    => clm3%g%l%c%p%pes%Tleaf_under_sunlit_EASS
      Tleaf_over_shaded_EASS     => clm3%g%l%c%p%pes%Tleaf_over_shaded_EASS
      Tleaf_under_shaded_EASS    => clm3%g%l%c%p%pes%Tleaf_under_shaded_EASS
      radsb_over_sunlit_EASS     => clm3%g%l%c%p%pef%radsb_over_sunlit_EASS
      radsb_over_shaded_EASS     => clm3%g%l%c%p%pef%radsb_over_shaded_EASS
      radsb_under_sunlit_EASS    => clm3%g%l%c%p%pef%radsb_under_sunlit_EASS
      radsb_under_shaded_EASS    => clm3%g%l%c%p%pef%radsb_under_shaded_EASS
      radnet_over_EASS           => clm3%g%l%c%p%pef%radnet_over_EASS
      radnet_under_EASS          => clm3%g%l%c%p%pef%radnet_under_EASS
      radnet_grnd_EASS           => clm3%g%l%c%p%pef%radnet_grnd_EASS
      radnet_EASS                => clm3%g%l%c%p%pef%radnet_EASS
      lwnet_EASS                 => clm3%g%l%c%p%pef%lwnet_EASS
      lwrad_grnd_EASS            => clm3%g%l%c%p%pef%lwrad_grnd_EASS
      swnet_EASS                 => clm3%g%l%c%p%pef%swnet_EASS
      swnet_over_EASS            => clm3%g%l%c%p%pef%swnet_over_EASS
      swnet_under_EASS           => clm3%g%l%c%p%pef%swnet_under_EASS
      swnet_grnd_EASS            => clm3%g%l%c%p%pef%swnet_grnd_EASS
      par_over_sunlit_EASS       => clm3%g%l%c%p%pef%par_over_sunlit_EASS
      par_over_shaded_EASS       => clm3%g%l%c%p%pef%par_over_shaded_EASS
      par_under_sunlit_EASS      => clm3%g%l%c%p%pef%par_under_sunlit_EASS
      par_under_shaded_EASS      => clm3%g%l%c%p%pef%par_under_shaded_EASS
      dalb_over_EASS             => clm3%g%l%c%p%pef%dalb_over_EASS
      dalb_under_EASS            => clm3%g%l%c%p%pef%dalb_under_EASS
      dalb_grnd_EASS             => clm3%g%l%c%p%pef%dalb_grnd_EASS
      frac_veg_nosno            => clm3%g%l%c%p%pps%frac_veg_nosno
!-------------------------------------------------

!#      do f=1, fn
!#         p=filterp(f)
!#         c=pcolumn(p)
!#         g=pgridcell(p)

         !(0.1) clearness index (or cloudiness index, referred from SeparateRadiation.c in 
         ! Distributed Hydrology Soil Vegetation Model)
         ! (Gueymard, 2000. Prediction and performance assessment of mean hourly global radiation,
         ! solar energy, 68(3), 285-303)
         if (coszen(c) .lt. 0.001_r8) then
              cl=0.0_r8
         else
              cl=forc_solar(g)/(1367.0_r8*coszen(c))
         end if

         !(0.2) diffused solar raidation(swdif)and direct solar radiation (swdir)
         if(cl .gt. 0.8_r8) then
              swdif_EASS=0.13_r8*forc_solar(g)
         else
              swdif_EASS=(0.943_r8+0.734_r8*cl-4.9_r8*cl*cl+1.796_r8*cl*cl*cl+2.058_r8*cl*cl*cl*cl)*forc_solar(g)
         end if

         swdif_EASS=max(0.0_r8,swdif_EASS)
         swdif_EASS=min(swdif_EASS, forc_solar(g)) 
         swdir_EASS=forc_solar(g)-swdif_EASS

         !(0.3) longwave radiation                  
         epsilon=1.0_r8-exp(-1.0_r8*((espar/1000.0_r8 *10.0_r8)** (forc_t(g)/1200.0_r8)) )     ! "espar_EASS" unit change [pa]->[Kpa]
         epsilon=min(1.0_r8, epsilon)
         epsilon=max(0.7_r8, epsilon)

         if (forc_lwrad(g) .lt. -100.0_r8) then
             lwrad_EASS=epsilon*sb*forc_t(g)**4.0_r8         ! (W m-2 K-4) * K^4=W m-2
         else 
             lwrad_EASS=forc_lwrad(g)          ![w m-2]
         end if      

         if (frac_veg_nosno(p)==0._r8) then

               !(1) solar radiation
               if ((forc_solar(g) .gt. 1.0e-10_r8) .and. (coszen(c) .gt. 1.0e-10_r8)) then
                      ! direct solar raidation of ground
                      swdir_grnd_EASS=swdir_EASS*dalb_grnd_EASS(p)

                      ! diffused solar raidation of ground
                      swdif_grnd_EASS=swdif_EASS*dalb_grnd_EASS(p)                                     
               else
                      swdir_grnd_EASS=0.0_r8
                      swdif_grnd_EASS=0.0_r8
               end if

               !(2)longwave radiation
               Tgrnd_EASS=t_grnd(c)
               temp_grnd=em_grnd*sb*Tgrnd_EASS**4.0_r8
               lwrad_grnd_EASS(p)=em_grnd*lwrad_EASS-temp_grnd

               !(3) output
               !(3.1)net radiation
               radnet_grnd_EASS(p)  =swdir_grnd_EASS+swdif_grnd_EASS+lwrad_grnd_EASS(p)
               radnet_EASS(p)=radnet_grnd_EASS(p)

               swnet_EASS(p)=swdir_grnd_EASS+swdif_grnd_EASS
               lwnet_EASS(p)=lwrad_grnd_EASS(p)

               !(3.2) Add the following to aviod NaN
               radnet_over_EASS(p)       =0._r8
               radnet_under_EASS(p)      =0._r8
               radsb_over_sunlit_EASS(p) =0._r8
               radsb_over_shaded_EASS(p) =0._r8
               radsb_under_sunlit_EASS(p)=0._r8
               radsb_under_shaded_EASS(p)=0._r8
               par_over_sunlit_EASS(p)   =0._r8
               par_over_shaded_EASS(p)   =0._r8
               par_under_sunlit_EASS(p)  =0._r8
               par_under_shaded_EASS(p)  =0._r8

         else 

               !(0.4) LAI
               lt_over_EASS=elai_over_EASS(p)+esai_over_EASS(p)
               lt_under_EASS=elai_under_EASS(p)+esai_under_EASS(p)
               elai_over_under_EASS=elai_over_EASS(p)+elai_under_EASS(p)

               phase='energy'
               call SunlitLAI_EASS(elai_over_EASS(p), clumping_EASS(p), coszen(c), phase, flc_elai_over)        ! lai_o
               call SunlitLAI_EASS(elai_under_EASS(p), clumping_EASS(p), coszen(c), phase, flc_elai_under)      ! lai_u
               call SunlitLAI_EASS(elai_over_under_EASS, clumping_EASS(p), coszen(c), phase, flc_elai_over_under)      ! lai_o+lai_u
               call SunlitLAI_EASS(lt_over_EASS, clumping_EASS(p), coszen(c), phase, flc_lt_over)               ! lai_o+lai_u

               qita_over=0.537_r8+0.025_r8*elai_over_EASS(p)
               qita_under=0.537_r8+0.025_r8*elai_under_EASS(p)

               call SunlitLAI_EASS(elai_over_EASS(p), clumping_EASS(p), qita_over, phase, flq_elai_over)   ! lai_o
               call SunlitLAI_EASS(elai_under_EASS(p), clumping_EASS(p), qita_under, phase, flq_elai_under)! lai_u
               call SunlitLAI_EASS(lt_over_EASS, clumping_EASS(p), qita_over, phase, flq_lt_over)          ! lai_o+lai_u
               call SunlitLAI_EASS(lt_under_EASS, clumping_EASS(p), qita_under, phase, flq_lt_under)       ! lai_o+lai_u

               !(1) solar radiation
               if ((forc_solar(g) .gt. 1.0e-10_r8) .and. (coszen(c) .gt. 1.0e-10_r8)) then
                      ! direct solar raidation of overstroy
                      swdir_over_EASS=swdir_EASS*(dalb_over_EASS(p)-dalb_under_EASS(p)* flc_elai_over)

                      ! direct solar raidation of understroy
                      swdir_under_EASS=swdir_EASS*flc_elai_over*(dalb_under_EASS(p)-dalb_grnd_EASS(p)*flc_elai_under)

                      ! direct solar raidation of ground
                      swdir_grnd_EASS=swdir_EASS*flc_elai_over_under*dalb_grnd_EASS(p)

                      ! total solar radiation aborbed by overstory 
                      pardir_over_EASS=0.5_r8*swdir_EASS/coszen(c)       ! 0.25 is changed from 0.5
      
                      ! total solar radiation aborbed by understory
                      pardir_under_EASS=0.5_r8*swdir_EASS*flc_lt_over/coszen(c)
                else
                      swdir_over_EASS=0.0_r8
                      swdir_under_EASS=0.0_r8
                      swdir_grnd_EASS=0.0_r8                                                                                                                                              
                      pardir_over_EASS=0.0_r8
                      pardir_under_EASS=0.0_r8
                end if

                !(2) diffused solar radiation
                if ((forc_solar(g) .gt. 1.0e-10_r8) .and. (coszen(c) .gt. 1.0e-3_r8)) then
                      ! diffused solar raidation of overstroy
                      swdif_over_EASS=swdif_EASS*(dalb_over_EASS(p)-dalb_under_EASS(p)*flq_elai_over) &
                               +0.21*swdir_EASS*clumping_EASS(p)*(1.1_r8-0.1_r8*elai_over_EASS(p))*exp(-coszen(c))

                      ! diffused solar raidation of understroy
                      swdif_under_EASS=swdif_EASS*flq_elai_over*(dalb_under_EASS(p)-dalb_grnd_EASS(p)*flq_elai_under) &
                                +0.21_r8*swdir_EASS*clumping_EASS(p)*flc_elai_over*(1.1_r8-0.1_r8*elai_under_EASS(p))*exp(-coszen(c))

                      ! diffused solar raidation of ground
                      swdif_grnd_EASS=swdif_EASS*dalb_grnd_EASS(p)*exp(-0.5_r8*clumping_EASS(p)*elai_over_EASS(p) &
                               /qita_over-0.5_r8*clumping_EASS(p)*elai_under_EASS(p)/qita_under)

                      ! total solar radiation aborbed by overstory 
                      pardif_over_EASS=swdif_EASS*(1.0_r8-flq_lt_over)/lt_over_EASS &
                                +0.07_r8*swdir_EASS*(1.1_r8-0.1_r8*lt_over_EASS)*exp(-coszen(c))

                      ! total solar radiation aborbed by understory
                      pardif_under_EASS=swdif_EASS*flq_elai_over*(1.0_r8-flq_lt_under)/lt_under_EASS &
                                 +0.05_r8*swdir_EASS*flc_elai_over*(1.1_r8-0.1_r8*lt_under_EASS)*exp(-coszen(c))
                else
                      swdif_over_EASS=0.0_r8
                      swdif_under_EASS=0.0_r8
                      swdif_grnd_EASS=0.0_r8
                      pardif_over_EASS=0.0_r8
                      pardif_under_EASS=0.0_r8
                end if

                !(3) longwave radiation
                ! canopy temperature(including overstory and understory) and ground temperature     
                h=lt_over_sunlit_EASS(p)+lt_over_shaded_EASS(p)
                call CanopyFlux_EASS(Tleaf_over_sunlit_EASS(p), Tleaf_over_shaded_EASS(p), h, &
                               lt_over_sunlit_EASS(p), lt_over_shaded_EASS(p), Tleaf_over_EASS(p) )
        
                h=lt_under_sunlit_EASS(p)+lt_under_shaded_EASS(p)
                call CanopyFlux_EASS(Tleaf_under_sunlit_EASS(p), Tleaf_under_shaded_EASS(p), h, &
                               lt_under_sunlit_EASS(p), lt_under_shaded_EASS(p),Tleaf_under_EASS(p) )

                !Tgrnd_EASS(c)=forc_t_col(c)
                Tgrnd_EASS=t_grnd(c)

                temp_over=em_over*sb*Tleaf_over_EASS(p)**4.0_r8
                temp_under=em_under*sb*Tleaf_under_EASS(p)**4.0_r8
                temp_grnd=em_grnd*sb*Tgrnd_EASS**4.0_r8


                ! overstory
                !constent 1.85 changed from 2.0 by chenjinn,Aug 23, 2012 

 !#               lwrad_over_EASS=(em_over*(lwrad_EASS+temp_under*(1.0_r8-flq_elai_under)+temp_grnd*flq_elai_under)&
 !#                        - 2.0_r8*temp_over)*(1.0_r8- flq_elai_over)+(1.0_r8-em_under)*(1.0_r8-flq_elai_under)&
 !#                        * (lwrad_EASS * flq_elai_over + temp_over * (1.0_r8-flq_elai_over))
 
                lwrad_over_EASS=(em_over*(lwrad_EASS+temp_under*(1.0_r8-flq_elai_under)+temp_grnd*flq_elai_under)&
                         - 1.85_r8*temp_over)*(1.0_r8- flq_elai_over)+(1.0_r8-em_under)*(1.0_r8-flq_elai_under)&     !1.85
                         * (lwrad_EASS * flq_elai_over + temp_over * (1.0_r8-flq_elai_over))

   
                ! understory
 !#               lwrad_under_EASS=(em_under*(lwrad_EASS*flq_elai_over+em_over*(1.0_r8-flq_elai_over)+temp_grnd)-2.0_r8*temp_under) &
 !#                         *(1.0_r8-flq_elai_under)+(1.0_r8-em_grnd)*(lwrad_EASS*flq_elai_over+temp_over*(1.0_r8-flq_elai_over) &
 !#                         *flq_elai_under+temp_under*(1.0_r8-flq_elai_under)) &
 !#                         + (1.0_r8-em_over) * (temp_under*(1.0_r8-flq_elai_under)+temp_grnd*flq_elai_under)*(1.0_r8-flq_elai_under)

                lwrad_under_EASS=(em_under*(lwrad_EASS*flq_elai_over+em_over*(1.0_r8-flq_elai_over)+temp_grnd)-1.85_r8*temp_under) &   !1.85
                          *(1.0_r8-flq_elai_under)+(1.0_r8-em_grnd)*(lwrad_EASS*flq_elai_over+temp_over*(1.0_r8-flq_elai_over) &
                          *flq_elai_under+temp_under*(1.0_r8-flq_elai_under)) &
                          + (1.0_r8-em_over) * (temp_under*(1.0_r8-flq_elai_under)+temp_grnd*flq_elai_under)*(1.0_r8-flq_elai_under)

                ! ground
                lwrad_grnd_EASS(p)=em_grnd*(lwrad_EASS*flq_elai_over+temp_over*(1.0_r8-flq_elai_over)*flq_elai_under+temp_under &
                         *(1.0_r8-flq_elai_under))-temp_grnd+(1.0_r8-em_under)*temp_grnd*(1.0_r8-flq_elai_under)

                !------------------------------------
                ! (4) output
                ! (1)net radiation 
                !    (output value only, not used in the model)
                radnet_over_EASS(p)  =swdir_over_EASS+swdif_over_EASS+lwrad_over_EASS
                radnet_under_EASS(p) =swdir_under_EASS+swdif_under_EASS+lwrad_under_EASS
                radnet_grnd_EASS(p)  =swdir_grnd_EASS+swdif_grnd_EASS+lwrad_grnd_EASS(p)
         !#       radnet_EASS(p)=radnet_over_EASS(p)+radnet_under_EASS(p)+radnet_grnd_EASS(p)

         !#       swnet_EASS(p)=swdir_over_EASS+swdif_over_EASS+swdir_under_EASS+swdif_under_EASS+swdir_grnd_EASS+swdif_grnd_EASS
                swnet_over_EASS(p)  =swdir_over_EASS+swdif_over_EASS
                swnet_under_EASS(p) =swdir_under_EASS+swdif_under_EASS
                swnet_grnd_EASS(p)  =swdir_grnd_EASS+swdif_grnd_EASS
                swnet_EASS(p)=swnet_over_EASS(p)+swnet_under_EASS(p)+swnet_grnd_EASS(p)
                lwnet_EASS(p)=lwrad_over_EASS+lwrad_under_EASS+lwrad_grnd_EASS(p)

                ! (2)total solar radiation aborbed by sunlit and shaded leaves, including overstory and understory 
                ! (2.1) sunlit leaves, overstory
                if (lt_over_sunlit_EASS(p) .gt. 0.0_r8 ) then
                       radsb_over_sunlit_EASS(p)= swdir_over_EASS/lt_over_sunlit_EASS(p)+lwrad_over_EASS/lt_over_EASS
                else
                       radsb_over_sunlit_EASS(p)= swdir_over_EASS+lwrad_over_EASS
                end if

                ! (2.2) shaded leaves, overstory
                if (lt_over_shaded_EASS(p) .gt. 0.0_r8 ) then
                       radsb_over_shaded_EASS(p)= swdif_over_EASS/lt_over_shaded_EASS(p)+lwrad_over_EASS/lt_over_EASS
                else
                       radsb_over_shaded_EASS(p)= swdif_over_EASS+lwrad_over_EASS
                end if

                ! (2.3) sunlit leaves, understory
                if (lt_under_sunlit_EASS(p) .gt. 0.0_r8 ) then
                       radsb_under_sunlit_EASS(p)= swdir_under_EASS/lt_under_sunlit_EASS(p)+lwrad_under_EASS/lt_under_EASS
                else
                       radsb_under_sunlit_EASS(p)= swdir_under_EASS+lwrad_under_EASS
                end if

                ! (2.4) shaded leaves, understory
                if (lt_under_shaded_EASS(p) .gt. 0.0_r8 ) then
                       radsb_under_shaded_EASS(p)= swdif_under_EASS/lt_under_shaded_EASS(p)+lwrad_under_EASS/lt_under_EASS
                else
                       radsb_under_shaded_EASS(p)= swdif_under_EASS+lwrad_under_EASS
                end if

                !(3) average absorbed PAR for sunlit and shaded leaves, including overstory and undstory
                par_over_sunlit_EASS(p)=(pardif_over_EASS+pardir_over_EASS)*dalb_over_EASS(p)
                par_over_shaded_EASS(p)= pardif_over_EASS*dalb_over_EASS(p)
                par_under_sunlit_EASS(p)=(pardif_under_EASS+pardir_under_EASS)*dalb_under_EASS(p)
                par_under_shaded_EASS(p)= pardif_under_EASS*dalb_under_EASS(p)

         end if

 !#   end do !(f)

  end subroutine EnergyBalance_EASS

end module EnergyBalance_EASSMod
