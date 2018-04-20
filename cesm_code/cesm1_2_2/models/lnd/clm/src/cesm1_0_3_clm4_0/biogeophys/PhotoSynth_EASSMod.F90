module PhotoSynth_EASSMod 
!--------------------------------------------------
!BOP
!! USE:
    use shr_kind_mod       , only : r8 => shr_kind_r8   
    use clm_varctl         , only : iulog

!! PUBLIC TYPES:
   implicit none
   save

!! PUBLIC MEMBER FUNCTIONS:
   public :: PhotoSynth_EASS         ! calculate air resistnace and leaf boundry resistance

!! PRIVATE MEMBER FUNCTIONS:
   private :: TBOLZ                 !Boltzmann tempertaure distribution for photosynthesis
   private :: SFC_VPD               !computes the relative humidity at the leaf surface
   private :: TEMP_FUNC             !

!! RVESION HISTORY

!EOP
!--------------------------------------------------
   contains
!-------------------------------------------------
!! BOP
!! IROUTINE: PhotoSynth_EASS
!! INTERFACE
   subroutine PhotoSynth_EASS(p, g, c, Tleaf_EASS, &                        ! input var.
                              par_EASS,rb_EASS,dayl_factor, phase,  &        ! input var.
                              ci_EASS,psn_EASS,cs_EASS, rs_EASS, vcmax, jmax)    ! output var.

!! DESRIPTION:
! c3 Leaf stomatal resistance and leaf photosynthesis(net co2 assimilation rate, A) of the fortran version of EASS model
! the clumping indix(clumping_EASS) and Extinction coefficient (kext) is adopted from 
! Chen, J.M., Li, J., Leblanc, S.G., Lacaze, R., Roujean, J.L., 2003. Multi-angular optical remote sensing for assessing vegetation structureand carbon abs.or.ption. Remote Sensing of Environment, 84, 516-525.

!! CALLED FROM:
!  subroutine LeafTemper , to calculate rs, ci, psn

!!REVISION HISTORY:
! 12 January 2012: Jing Chen, F90 Revision of ../eass_1031/EASS_region_1028/db_PHOTOSYNTHESIS.c

!!USE:
        use clmtype
        use clm_varcon      , only : rgas,tfrz   !=8.314*1.0e3 [J K-1 kmol-1]
        use shr_const_mod   , only : SHR_CONST_PSTD,SHR_CONST_PI  !standard pressure =101325.0 [pa]
        use QsatMod         , only : Qsat
        use clm_atmlnd      , only : clm_a2l

!! ARGUMENTS:
       implicit none
       integer,intent(in) :: p                    ! pft index [-]
       integer,intent(in) :: g                    ! grid cell index[-]
       integer,intent(in) :: c                    ! column index[-]
       real(r8),intent(in) :: Tleaf_EASS          ! leaf temperature[K]
       real(r8),intent(in) :: par_EASS            ! absorbed PAR for leaves[w m-2]
       real(r8),intent(in) :: rb_EASS             ! leaf boundary resistance [s m-1]
       real(r8),intent(in) :: dayl_factor         ! scalar (0-1) for daylength
       real(r8),intent(inout) :: ci_EASS          ! intracellular leaf CO2 concentration[pa]
       real(r8),intent(out) :: psn_EASS           ! foliage photosynthesis[umol CO2 m-2 s-1]
       real(r8),intent(out) :: cs_EASS            ! CO2 concentration at the leaf surface(cs_EASS)[pa]
       real(r8),intent(out) :: rs_EASS            ! stomatal resistance(rs)[m s-1] 
       real(r8),intent(out) :: vcmax              ! maximun rate of carbonxylation [umol CO2 m-2 s-1]
       real(r8),intent(out) :: jmax               ! light-saturated rate of electron transport 
                                                  ! in the photosynthetic carbon reduction cycle in leaf cells[umol m-2 s-1]
       character(len=*),intent(in) :: phase       ! 'over_sunlit','over_shaded','under_sunlit','under_shaded'

       !!LOCAL VARIABLES:
       !local pointer to implicait in varialbles
       integer, pointer :: ivt(:)                  ! pft vegetation type
       real(r8), pointer :: forc_pbot(:)           ! atmospheric pressure (Pa)           
       real(r8), pointer :: forc_pco2(:)           ! partial pressure co2 [Pa]
       real(r8), pointer :: forc_po2(:)            ! partial pressure o2  [Pa]
       real(r8), pointer :: coszen(:)              ! cosine of solar zenith angle
       real(r8), pointer :: btran(:)               ! soil water transpiration factor (0 to 1)[-]

       real(r8), pointer :: clumping_EASS(:)       ! clumping index
       real(r8), pointer :: elai_EASS(:)           ! one-sided leaf area index with burying by snow, corrected by clumping index, in EASS
       real(r8), pointer :: vcmax25_EASS(:)        ! maximum capacity of rubisco at 25C-Vcmax[umol CO2 m-2 s-1]  parameter1[36] in readparam.c 
       real(r8), pointer :: Nleaf_EASS(:)          ! leaf Nitrogen content    (mean value + 1 SD g/m2)    parameter1[46] in readparam.c 
       real(r8), pointer :: slope_vcmaxN_EASS(:)   ! slope of Vcmax-N curve                              parameter1[47] in readparam.c
       real(r8), pointer :: qe25(:)        ! quantum efficiency at 25C (umol CO2 / umol photon)
       real(r8), pointer :: forc_q(:)      ! atmospheric specific humidity (kg/kg)
       real(r8), pointer :: qg(:)          ! specific humidity at ground surface [kg/kg]
       real(r8), pointer :: c3psn(:)       ! photosynthetic pathway: 0. = c4, 1. = c3
       real(r8), pointer :: mp_EASS(:)          ! slope of conductance-to-photosynthesis relationship

       !!OHTER LOCAL VARIABLES:
       real(r8), parameter :: knext=0.3            ! Extinction coefficient, constent
       real(r8), parameter :: G_theta_EASS=0.5     ! mean rojection of leaf mormals in the direction of zentih angle, 
                                                   !Takeing the usual value of 0.5 for K under the assumption of random leaf angle distribution
       real(r8), parameter :: evc_EASS=55000.0     ! activation energy for carboxylation[J mol-1]
       real(r8), parameter :: ejm_EASS=55000.0     ! activation energy for electron transport[J mol-1]
       real(r8), parameter :: toptvc_EASS=301.0    ! optimum temperature for maximum carboxylation[K]
       real(r8), parameter :: toptjm_EASS=301.0    ! optimum temperature for maximum electron transport[K]
       real(r8), parameter :: t25k_EASS=298.16     ! absolute temperature at 25 degree celsius[Kelvin]
       real(r8), parameter :: tau25_EASS=2904.12   ! tau coefficient[-]( =specificity factor(102.33)Kco2/Kh2o(28.38) )
                                                   ! Similar number from Dreyer et al. 2001, Tree Physiol, tau= 2710
       real(r8), parameter :: ektau_EASS=-29000.0  ! [J mol-1] (Jordan and Ogren, 1984)
       real(r8), parameter :: kc25_EASS=274.6      ! kinetic coef for CO2 at 25 C[microbars, ubr]     1ubar=10^-3mbar=0.1pa
       real(r8), parameter :: ko25_EASS=419.8      ! kinetic coef for O2 at 25C [millibars, mbr]      1mbar=1hpa=100pa  
       real(r8), parameter :: ekc_EASS=80500.0     ! Activation energy for K of CO2[J mol-1]
       real(r8), parameter :: eko_EASS=14500.0     ! Activation energy for K of O2[J mol-1]
       real(r8), parameter :: erd_EASS=38000.0     ! activation energy for dark respiration, eg Q10=2[J mol-1]
       real(r8) :: kext                            ! Extinction coefficient, depending on clumping index and solar zenith angle
       real(r8) :: expr1                           ! intermediate variable in calculate vcmax
       real(r8) :: expr2                           ! intermediate variable in calculate vcmax
       real(r8) :: expr3                           ! intermediate variable in calculate vcmax    
       real(r8) :: vcmax25                         ! vcmax at 25 degree Celsius [umol CO2 m-2 s-1]
       real(r8) :: ppf                             ! absorb photosynthetic photon flux denisty [umol photons m-2 s-1]
       real(r8) :: PPFD                            ! electon transport  [umol co2 m-2 s-1] 
       real(r8) :: dt25                            ! temperature difference =leaf temperature - 298.16 [Kelvin]
       real(r8) :: tau                             ! to convert for value in solution to that based in air[-]
       real(r8) :: cp                              ! co2 compensation point [pa]
       real(r8) :: kc                              ! co2 michaelis-menten constant [pa]
       real(r8) :: ko                              ! o2 michaelis-menten constant [pa]
       real(r8) :: kenzyme                         ! enzyme kinetics[pa]
       real(r8) :: wc                              ! rubisco limited photosynthesis [umol CO2 m-2 s-1]
       real(r8) :: wj                              ! light limited photosynthesis [umol CO2 m-2 s-1]
       real(r8) :: we                              ! export limited photosynthesis [umol co2/m**2/s]
       real(r8) :: jmax25                          ! Jmax at 25 oC(=2.39parameter1[37]-14.2)  [umol m-2 s-1]
       real(r8) :: j_photon                        ! maximum value of electron transport rate[umol m-2 s-1]
       real(r8) :: rb_mole                         ! leaf boundary resistance [m2 s mol-1]
       real(r8) :: kc25                            ! kinetic coef for co2 at 25 C[pa](=kc25_EASS*0.1)
       real(r8) :: ko25                            ! kinetic coef for o2 at 25C [pa](=ko25_EASS*100)
       real(r8) :: eleaf                           ! vapor pressure on leaf surface [pa]
       real(r8) :: deleafdT                        ! derivative of "el" on "t_veg" [pa/K] ( var. in CLM, not used in EASS)
       real(r8) :: qsatleaf                        ! leaf specific humidity [kg/kg]( var. in CLM, not used in EASS)
       real(r8) :: qsatleafdT                      ! derivative of "qsatl" on "t_veg"( var. in CLM, not used in EASS)
       real(r8) :: cf                              ! s m2/umol -> s/m
       real(r8) :: rs_mole                         ! stomatal resistance [umol m-2 s-1]
       real(r8) :: q_can   ! humidity of canopy canopy air [kg kg-1]
       real(r8) :: es_can  ! constrain of vapor pressure of canopy air(es_can)[pa]
       real(r8) :: cea     ! constrain ea or else model blows up
       real(r8) :: iter    ! iteration index
       real(r8) :: atmp    ! intermediate calculations for rs
       real(r8) :: btmp    ! intermediate calculations for rs
       real(r8) :: ctmp    ! intermediate calculations for rs
       real(r8) :: q       ! intermediate calculations for rs
       real(r8) :: r1,r2   ! roots for rs
       real(r8) :: tc      ! leaf temperature (degree celsius)
       integer , parameter :: niter = 3          ! number of iterations
       real(r8), parameter :: bp=2000._r8        ! minimum leaf conductance (umol/m2/s)
       real(r8), parameter :: rsmax0 = 2.e4_r8   ! maximum stomatal resistance [s/m]
       real(r8), parameter :: akc = 2.1_r8       ! q10 for kc25
       real(r8), parameter :: ako = 1.2_r8       ! q10 for ko25
   !#    real(r8), parameter :: kc25 = 30._r8      ! co2 michaelis-menten constant at 25c [pa]
   !#    real(r8), parameter :: ko25 = 30000._r8   ! o2 michaelis-menten constant at 25c [pa]

       !----------------------------------------------------------------------

       ! Assign local pointers to derived type members (pft-level)
       ivt              => clm3%g%l%c%p%itype
       clumping_EASS    => clm3%g%l%c%p%pps%clumping_EASS
       elai_EASS        => clm3%g%l%c%p%pps%elai_EASS
       btran            => clm3%g%l%c%p%pps%btran

       ! Assign local pointers to derived type members (gridcell-level) 
       forc_pco2        => clm_a2l%forc_pco2
       forc_po2         => clm_a2l%forc_po2
       forc_pbot        => clm_a2l%forc_pbot

       ! Assign local pointers to derived type members (column-level) 
       coszen           => clm3%g%l%c%cps%coszen

       ! Assign local pointers to pft constants
       qe25             => pftcon%qe25
       vcmax25_EASS     => pftcon%vcmax25_EASS
       Nleaf_EASS       => pftcon%Nleaf_EASS
       slope_vcmaxN_EASS=> pftcon%slope_vcmaxN_EASS
       qg               => clm3%g%l%c%cws%qg
       c3psn            => pftcon%c3psn
       mp_EASS          => pftcon%mp_EASS
       forc_q           => clm_a2l%forc_q

!-------------------------------------------------------------------------

        !(1) vcmax (maximum carboxylation rate[umol CO2 m-2 s-1])
        !         c version in eass_m31.c

        kext = G_theta_EASS * clumping_EASS(p) / coszen(c)
        expr1=1.0_r8-exp(-kext*elai_EASS(p))
        expr2=1.0_r8-exp(-(knext+kext)*elai_EASS(p))
        expr3=1.0_r8-exp(-knext*elai_EASS(p))

        if (phase == 'over_sunlit' .or. phase == 'under_sunlit' ) then
           if (expr1 .gt. 0.0_r8) then
                vcmax25=vcmax25_EASS(ivt(p))*Nleaf_EASS(ivt(p))*slope_vcmaxN_EASS(ivt(p)) &   !modefied by Chenjinn
                      *kext*expr2/(kext+knext)/expr1
            else
                vcmax25=vcmax25_EASS(ivt(p))
            end if
        else if (phase == 'over_shaded' .or. phase == 'under_shaded' ) then
            if ( (kext .gt. 0.0_r8) .and. (elai_EASS(p) .gt. (expr1/kext)) ) then
                vcmax25=vcmax25_EASS(ivt(p))*Nleaf_EASS(ivt(p))*slope_vcmaxN_EASS(ivt(p)) &
                        *(expr3/knext-expr2/(kext+knext))/(elai_EASS(p)-expr1/kext)
            else
                        vcmax25=vcmax25_EASS(ivt(p))
            end if
        end if

        vcmax= TBOLZ(vcmax25,evc_EASS,toptvc_EASS,Tleaf_EASS)
        vcmax= vcmax * btran(p) * dayl_factor

        ! (2) electron transport rate(j_photon) [umol m-2 s-1]
        jmax25= 2.39_r8 * vcmax25_EASS(ivt(p))-14.2_r8    !# 1.2 for Sen.test
        jmax=TBOLZ(jmax25, ejm_EASS, toptjm_EASS, Tleaf_EASS)   ![umol CO2 m-2 s-1] 

        ppf=par_EASS*4.55_r8       ! [umol photons m-2 s-1]. par_EASS, par absorbed per unit lai [w m-2(=J s-1 m-2)]                                                   
   !#     PPFD = 0.5_r8 * ppf      ! electon transport  [umol co2 m-2 s-1]
                                 ! 0.5, quantum efficiency at 25C [(umol co2)*(umol photons)-1]
        PPFD=ppf*qe25(ivt(p))

        j_photon= jmax * PPFD/ (PPFD + 2.1_r8*jmax)    ![umol co2 s-1 m-2]
      
   !#     ppf = 4.6_r8 * par_EASS
   !#     j_photon = ppf * qe25(ivt(p))

        !(3) co2 compensation point of without dark respiration(cp)[pa]
 !#       tc = Tleaf_EASS - tfrz
 !#       kc = kc25 * akc **((tc-25.0_r8)/10.0_r8)
 !#       ko = ko25 * ako **((tc-25.0_r8)/10.0_r8)
 !#       kenzyme = kc * (1._r8+forc_po2(g)/ko)
 !#       cp = 0.5_r8*kc/ko*forc_po2(g)*0.21_r8

        dt25 = Tleaf_EASS - t25k_EASS            ![K]
        tau    = TEMP_FUNC(tau25_EASS, ektau_EASS, dt25, t25k_EASS, Tleaf_EASS)    ! unit[-]
        ! unit  tau25_EASS=2904.12  [-]
        !       ektau_EASS=-29000.0 [J mol-1]
        !       dt25=Tleaf_EASS-298.16 [K]
        !       t25k_EASS=298.16    [K]
        !       Tleaf_EASS             [K]

        cp = 0.5_r8 * forc_po2(g) / tau            ![pa] 
        ! to make co2 compensation point consist with the value in CLM, the 0.5 is changed from 500

        !(4) Kenzyme[pa] (a function of enzyme kinetics)
        kc25 = kc25_EASS * 0.1_r8                    ! unit change[ubar]->[pa]
        kc    = TEMP_FUNC(kc25, ekc_EASS, dt25, t25k_EASS, Tleaf_EASS) ! [pa]

        ko25 = ko25_EASS * 100.0_r8                                        ! unit change[mbar]->[pa]
        ko = TEMP_FUNC(ko25, eko_EASS, dt25, t25k_EASS, Tleaf_EASS)        ! [pa]

        kenzyme = kc * (1.0_r8 + forc_po2(g) / ko)       ! [pa]

        !(5) boundary reistance[s m-1]->[s m2 umol-1]
        cf=forc_pbot(g)/(rgas*0.001_r8*Tleaf_EASS)*1.0e6_r8  ![umol pa J-1]=[umol m-3] 
        rb_mole=rb_EASS/cf

        !(6) constrain of vapor pressure of canopy air(es_can)[pa]
        q_can=(forc_q(g)+qg(c))/2._r8
        es_can=forc_pbot(g)*q_can/0.622_r8

        call QSat(Tleaf_EASS, forc_pbot(g), eleaf, deleafdT, qsatleaf, qsatleafdT)
        cea = max(0.25_r8*eleaf*c3psn(ivt(p))+0.40_r8*eleaf*(1._r8-c3psn(ivt(p))), min(es_can,eleaf) ) 
 
        !(7) foliage photosynthesis(psn, umol m-2 s-1)
        !      co2 concentration at leaf surface(cs, pa)
        !      leaf stomatal resistance (rs, s/m)
        !      intracellular leaf CO2 (ci, Pa)

        do iter = 1, niter
              ! (1.6) Rubisco-limited gross photosynthesis rate (wc, [umol m-2 s-1])
              !       light-limited gross photosynthesis rate (wj, [umol m-2 s-1])
              wc = max(ci_EASS - cp,0._r8)* vcmax   / (ci_EASS + kenzyme ) * c3psn(ivt(p)) + vcmax * ( 1.0_r8 - c3psn(ivt(p)) )       ![umol m-2 s-1]
              wj = max(ci_EASS - cp,0._r8)*j_photon /(( ci_EASS + 2.0_r8 * cp )) * c3psn(ivt(p)) + j_photon * ( 1.0_r8 - c3psn(ivt(p)) )
              we = 0.5_r8*vcmax*c3psn(ivt(p)) + 4000._r8*vcmax*ci_EASS/forc_pbot(g)*(1._r8-c3psn(ivt(p))) 
              psn_EASS=min(wc,wj,we)
              cs_EASS = max( forc_pco2(g)-1.37_r8*rb_mole*forc_pbot(g)*psn_EASS/4._r8, 1.0e-6_r8 )  ![pa]
              atmp = mp_EASS(ivt(p))*psn_EASS/4._r8*forc_pbot(g)*cea / (cs_EASS*eleaf) + bp      !# for Sen. 
              btmp = ( mp_EASS(ivt(p))*psn_EASS/4._r8*forc_pbot(g)/cs_EASS + bp ) * rb_mole - 1._r8   !# for Sen
              ctmp = -rb_mole
              if (btmp >= 0._r8) then
                 q = -0.5_r8*( btmp + sqrt(btmp*btmp-4._r8*atmp*ctmp) )
              else
                 q = -0.5_r8*( btmp - sqrt(btmp*btmp-4._r8*atmp*ctmp) )
              end if
              r1 = q/atmp
              r2 = ctmp/q
              rs_mole = max(r1,r2)
              ci_EASS = max( cs_EASS-psn_EASS/4._r8*forc_pbot(g)*1.65_r8*rs_mole, 0._r8 )
         end do

         rs_EASS= min(rsmax0, rs_mole*cf)

  end subroutine PhotoSynth_EASS

!-----------------------------------------------------
!BOP

!!INTERFACE:
   real(r8) function TBOLZ(rate, eakin, topt, tl)

!!DESRIPTION:
! Boltzmann tempertaure distribution for photosynthesis
 
!! CALLED FROM:
!subroutine Stomata in this module, to calculate vcmax

! !REVISION HISTORY:
! 12 January 2012: Jing Chen F90 Revision of ../eass_1031/EASS_region_1028/db_PHOTOSYNTHESIS.c

!!USES:
     use clm_varcon,    only: rgas

!! ARGUMENTS:
     implicit none
     real(r8), intent(in) :: rate         ! [umol CO2 m-2 s-1] 
     real(r8), intent(in) :: eakin        ! activation energy for carboxylation or electron transport[J mol-1]
     real(r8), intent(in) :: topt         ! [K]
     real(r8), intent(in) :: tl           ! [K]

!! LOCAL VARIABLES:
!EOP
     real(r8) :: dtlopt
     real(r8) :: prodt
     real(r8) :: numm
     real(r8) :: denom
     real(r8), parameter :: hkin=2.0e5_r8  !enthalpy term [J mol-1]

     dtlopt  = tl   - topt                                ! [K]
     prodt   = rgas/1000.0 * topt * tl                    ! [J K-1 mole-1] * K* K=J K mol-1
     numm    = hkin * exp(eakin * (dtlopt) / (prodt))     ! J mol-1 * exp(J mol-1 * K /(J K mol-1))=J mol-1
     denom   = hkin - eakin * (1.0_r8 - exp(hkin * (dtlopt) / (prodt)))   ! [J mol-1]- [J mol-1]*(1-exp([J mol-1]* K/J K mol-1))=[J mol-1]
     TBOLZ   = rate * numm / denom                        ! [umol CO2 m-2 s-1]

   end function TBOLZ

!------------------------------------------------------------------
!BOP

!!INTERFACE
   real(r8) function SFC_VPD(tlk, es_leaf, espar, qflu, rb_EASS)

!!DESRIPTION:
! this function computes the relative humidity at the leaf surface for application in the Ball Berry Equation
! latent heat flux, LE, is passed through the function, mol m-2 s-1    and it solves for the humidity at leaf surface
! es_leaf : saturation vapor pressure at leaf temperature.

! the algorithm of abs_hum is based on  http://zh.wikepedia.org , in 'absolute humidity' introduction 
! and the result is same as the old one, but the unit is clearer
! note that the unit change in 'abs_hum 'calculate: 1[J pa-1]=1[m3]
!-----------------------------------------------------------
!- prove:
!- condiation:      J = N * m  ; Pa= N /m^2                - 
!- so:              m^3=N*m^3/N=N*m*m^2/N=J/(N/m^2)=J/Pa   -
!-----------------------------------------------------------

!! CALLED FROM:
!subroutine Stomata in this module, to calculate vcmax

! !REVISION HISTORY:
! 12 January 2012: Jing Chen F90 Revision of ../eass_1031/EASS_region_1028/db_PHOTOSYNTHESIS.c

!!USES:
        use PenmanEqu_EASSMod, only: LHvaporFunc_EASS

!! ARGUMENTS:
        implicit none
        real(r8), intent(in) :: tlk            ! leaf temperature [K]
        real(r8), intent(in) :: es_leaf        ! saturated water vapor pressure on the leaf[pa]                                  
        real(r8), intent(in) :: espar          ! actual (unsaturated) water vapor pressure[pa]
        real(r8), intent(in) :: qflu           ! transpiration from leaves [W m-2]
        real(r8), intent(in) :: rb_EASS        ! leaf boundary resistance [s m-1]

!! LOCAL VARIABLES:
!EOP
        real(r8) :: abs_hum                     ! absolute humidity [kg m-3]
        real(r8) :: rhov_sfc                    ! [kg m-3 ]
        real(r8) :: espar_sfc                   ! unsaturated water vapor pressure on the leaf surface [kg m-3]
        real(r8) :: vpd_sfc                     ! difference between saturated and unsatruated water vapro pressure
                                                ! on the surface [pa]
        real(r8), parameter :: Rw=461.62_r8     ! the gas constant for water vapor [J kg-1 K-1]

!-------------------------------------------------------------------------------------------------------

        abs_hum   = espar/(Rw * tlk)    ! [pa] / ( J kg-1 K-1 *K)=pa kg J-1= kg m-3 

        rhov_sfc  = (qflu / LHvaporFunc_EASS(tlk)) * rb_EASS + abs_hum    !(J s-1 m-2 *J-1 kg)* s m-1 + kg m-3 = kg m-3
        espar_sfc = rhov_sfc * Rw * tlk              ![pa]
        vpd_sfc   = es_leaf - espar_sfc              ![pa]
        SFC_VPD   = 1.0_r8 - vpd_sfc / es_leaf       ! 0 to 1.0
        
        return
     
   end function SFC_VPD
!-----------------------------------------------------
   !BOP

   !!INTERFACE:
   real(r8) function TEMP_FUNC(rate, eact, tprime, tref, tlk)

   !!DESRIPTION:
   ! Arrhennius temperature function, formula for the tempertaure dependence of the reaction rate constant, 
   ! and therefor, rate of a chemical reation 
   ! the equation gives the dependance of the rate cosntant (rate) of chemical reactions on the temperature(tlk,[K]) 
   ! and activation energy (eact)
   ! the unit of "TEMP_FUNC" depanding on that of "rate"


   !! CALLED FROM:
   !subroutine Stomata in this module, to calculate vcmax

   !!REVISION HISTORY:
   ! 12 January 2012: Jing Chen F90 Revision of ../eass_1031/EASS_region_1028/db_PHOTOSYNTHESIS.c

   !!USES:
        use clm_varcon,    only: rgas

   !! ARGUMENTS:
        implicit none
        real(r8), intent(in) :: rate          ! per-exponential factor of simple the perfactor [-]
        real(r8), intent(in) :: eact          ! activation energy[J mol-1]
        real(r8), intent(in) :: tprime        ![K]
        real(r8), intent(in) :: tref          ![K]
        real(r8), intent(in) :: tlk           ![K]
       !-------------------------------------

        TEMP_FUNC = rate * exp(tprime * eact / (tref * rgas/1000.0_r8*tlk))  ! exp(K* J mol-1/(K*J K-1 mol-1)*K)=[-]                                                                     
        return

    end function TEMP_FUNC

end module PhotoSynth_EASSMod
