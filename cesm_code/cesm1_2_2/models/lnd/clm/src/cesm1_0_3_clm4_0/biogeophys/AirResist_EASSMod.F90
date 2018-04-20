module AirResist_EASSMod
!--------------------------------------------------
!BOP
!!MODULE:AirResist_EASSMod
!
!!DESCTIPTION:
!Calculates (1)air resitance(ra) of overstory, understoy and ground 
!(2) leaf boundry resitance(rb) of overstory and understoy, 

!!USE:
    use shr_kind_mod,    only : r8=>shr_kind_r8
    use clm_varctl,      only : iulog

!!PUBLIC TYPES:
   implicit none
   save

!!PUBLIC MEMBER FUNCTIONS:
   public :: AirResist_EASS         ! calculate air resistnace and leaf boundry resistance

!! RVESION HISTORY
!5/19/2012, Jing Chen: The module was created based on ...\eass_1031\EASS_region_1028\ra_o1.c
 
!EOP
!--------------------------------------------------
 contains
!-------------------------------------------------
!!BOP
!!IROUTINE: AirResist_EASS
!!INTERFACE
! subroutine AirResist_EASS(fn, filterp, lbp, ubp)  
 subroutine AirResist_EASS(p, c, g, cpair)  
!!DESCRIPTION:
!Calculates (1)air resitance(ra) of overstory, understoy and ground
!(2) leaf boundry resitance(rb) of overstory and understoy,

!!USES:
   use clm_varcon,       only : vkc, grav, tfrz, zlnd, zsno  
   use clmtype
   use clm_atmlnd,       only : clm_a2l       

!!ARGUMENTS:
   implicit none
 !  integer , intent(in) :: fn           ! size of pft filter
 !  integer , intent(in) :: filterp(fn)  ! pft filter
 !  integer , intent(in) :: lbp, ubp     ! pft bounds
    integer , intent(in) :: p            ! pft index
    integer , intent(in) :: c            ! column index
    integer , intent(in) :: g            ! gridcell index
    real(r8), intent(in) :: cpair        ! specific heat of dry air in EASS, depanding on 
                                         ! saturated water vapor specific humidity [J kg-1 k-1]


!!LOCAL VARIABLES:
! local pointer to implicit in variables
 !  integer , pointer :: pcolumn(:)      ! pft's column index
 !  integer , pointer :: pgridcell(:)    ! pft's gridcell index
   integer , pointer :: ivt(:)               ! pft vegetation type
   real(r8), pointer :: elai_over_EASS(:)    ! elai overstory
   real(r8), pointer :: esai_over_EASS(:)    ! esai understory
   real(r8), pointer :: clumping_EASS(:)     ! clumping index
   real(r8), pointer :: eflx_sh_over_EASS(:) ! sensible heat flux of overstory [w m-2]
   real(r8), pointer :: forc_hgt_u(:)        ! observational height of wind [m]
   real(r8), pointer :: forc_wind(:)         ! atmospheric wind speed[m s-1], calculate in lnd_comp_esmf.F90
   real(r8), pointer :: forc_rho(:)          ! atmospheric density [kg m-3]
   real(r8), pointer :: forc_t_col(:)        ! atmospheric temperature [K]
   real(r8), pointer :: htop_EASS(:)         ! canopy top (m)
   real(r8), pointer :: hunder_EASS(:)       ! the hight of understory[m]
   real(r8), pointer :: frac_snow_grnd_EASS(:)  ! fraction of snow on the ground(0~1)
   integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    real(r8), pointer :: eflx_sh_grnd_EASS(:)       ! sensible heat flux, ground [w m-2]

! local pointer to implicit inout variables
   real(r8), pointer :: rs_over_shaded_EASS(:)  ! stomatal resistant of shaded leaves, overstroy[s m-1]
   real(r8), pointer :: rs_over_sunlit_EASS(:)  ! stomatal resistant of sunlit leaves, overstroy[s m-1]
   real(r8), pointer :: rs_under_shaded_EASS(:) ! stomatal resistant of shaded leaves, understroy[s m-1]
   real(r8), pointer :: rs_under_sunlit_EASS(:) ! stomatal resistant of sunlit leaves, understroy[s m-1]

! local pointer to implicit out variables
   real(r8), pointer :: ra_over_EASS(:)         ! air resistance of overstor [s m-1]
   real(r8), pointer :: ra_under_EASS(:)        ! air resistance of understory [s m-1]
   real(r8), pointer :: rb_over_EASS(:)         ! leaf boundary resistance of overstory [s m-1]
   real(r8), pointer :: rb_under_EASS(:)        ! leaf boundary resistance of understory [s m-1]
   real(r8), pointer :: ra_grnd_EASS(:)         ! air resistance of grnd [s m-1]
   real(r8), pointer :: Gheat_over_EASS(:)      ! heat conductance, overstory[m s-1]
   real(r8), pointer :: Gheat_under_EASS(:)     ! heat conductance, overstory[m s-1]
   real(r8), pointer :: Gevap_over_EASS(:)      ! intercepted water conductance, overstory[m s-1]
   real(r8), pointer :: Gevap_under_EASS(:)     ! intercepted water conductance, overstory[m s-1]
   real(r8), pointer :: Gleaf_over_sunlit_EASS(:) ! leaf conductance of sunlit leaves , overstory[m s-1]
   real(r8), pointer :: Gleaf_over_shaded_EASS(:) ! leaf conductance of shaded leaves , overstory[m s-1]
   real(r8), pointer :: Gleaf_under_sunlit_EASS(:)! leaf conductance of sunlit leaves , understory[m s-1]
   real(r8), pointer :: Gleaf_under_shaded_EASS(:)! leaf conductance of shaded leaves , understory[m s-1]    

!!OTHER LOCAL VARIABLES:
 !  integer :: f                          ! filter index
 !  integer :: p                          ! pft index
 !  integer :: c                          ! column index
 !  integer :: g                          ! gridcell index
   real(r8) :: hdisp                      ! displacement height[m]
   real(r8) :: z0                         ! roughness lengh[m]
   real(r8) :: lt_over_EASS               ! (elai+esai) of overstory
   real(r8) :: ustar_EASS                 ! friction velocity [m s-1]
   real(r8) :: leng_EASS                  ! length[m^-1]
   real(r8) :: psi                        ! [-]
   real(r8) :: kh_o                       ! [m2 s-1]
   real(r8) :: kh_u                       ! [m2 s-1]
   real(r8) :: gamm                       ! temperary var.
   real(r8) :: Le                         ! [-]
   real(r8) :: windsp_hcan                ! wind speed at canopy top[m s-1]
   real(r8) :: windsp_hdisp               ! Wind speed at displancement height[m s-1]
   real(r8) :: windsp_hunder              ! wind speed at understory height [m s-1]
   real(r8) :: visco                      ! kinematic viscosity of air[m2 s-1]
   real(r8) :: Re                         ! Reynold's number[s m-1]
   real(r8) :: Nu                         ! Nusselt number[s m-1]
   real(r8) :: alfaw                      ! ??unit                                               

! Assign local pointers to pft constants
   hunder_EASS         => pftcon%hunder_EASS
   htop_EASS           => pftcon%htop_EASS

! Assign local pointers to derived type members(column-level)
   frac_snow_grnd_EASS   => clm3%g%l%c%cps%frac_snow_grnd_EASS

! Assign local pointers to derived type members(pft-level)
!#   pcolumn             => clm3%g%l%c%p%column
!#   pgridcell           => clm3%g%l%c%p%gridcell
   ivt                 => clm3%g%l%c%p%itype
   elai_over_EASS      => clm3%g%l%c%p%pps%elai_over_EASS
   esai_over_EASS      => clm3%g%l%c%p%pps%esai_over_EASS
   clumping_EASS       => clm3%g%l%c%p%pps%clumping_EASS
   eflx_sh_over_EASS   => clm3%g%l%c%p%pef%eflx_sh_over_EASS
   ra_over_EASS        => clm3%g%l%c%p%pps%ra_over_EASS
   ra_under_EASS       => clm3%g%l%c%p%pps%ra_under_EASS
   rb_over_EASS        => clm3%g%l%c%p%pps%rb_over_EASS
   rb_under_EASS       => clm3%g%l%c%p%pps%rb_under_EASS
   ra_grnd_EASS        => clm3%g%l%c%p%pps%ra_grnd_EASS
  !htop                => clm3%g%l%c%p%pps%htop
   Gheat_over_EASS     => clm3%g%l%c%p%pps%Gheat_over_EASS
   Gheat_under_EASS    => clm3%g%l%c%p%pps%Gheat_under_EASS
   Gevap_over_EASS     => clm3%g%l%c%p%pps%Gevap_over_EASS
   Gevap_under_EASS    => clm3%g%l%c%p%pps%Gevap_under_EASS
   Gleaf_over_sunlit_EASS  => clm3%g%l%c%p%pps%Gleaf_over_sunlit_EASS
   Gleaf_over_shaded_EASS  => clm3%g%l%c%p%pps%Gleaf_over_shaded_EASS
   Gleaf_under_sunlit_EASS => clm3%g%l%c%p%pps%Gleaf_under_sunlit_EASS
   Gleaf_under_shaded_EASS => clm3%g%l%c%p%pps%Gleaf_under_shaded_EASS
   rs_over_shaded_EASS     => clm3%g%l%c%p%pps%rs_over_shaded_EASS
   rs_over_sunlit_EASS     => clm3%g%l%c%p%pps%rs_over_sunlit_EASS
   rs_under_shaded_EASS    => clm3%g%l%c%p%pps%rs_under_shaded_EASS
   rs_under_sunlit_EASS    => clm3%g%l%c%p%pps%rs_under_sunlit_EASS
   frac_veg_nosno          => clm3%g%l%c%p%pps%frac_veg_nosno
   eflx_sh_grnd_EASS       => clm3%g%l%c%p%pef%eflx_sh_grnd_EASS

! Assign local pointers to derived subtypes components(gridcell-level)
   forc_wind           => clm_a2l%forc_wind 
   forc_rho            => clm_a2l%forc_rho
   forc_t_col          => clm3%g%l%c%ces%forc_t
   forc_hgt_u          => clm_a2l%forc_hgt_u

!-------------------------------------------------

 !  do f=1, fn
 !     p=filterp(f)
 !     c=pcolumn(p)
 !     g=pgridcell(p)

      if (frac_veg_nosno(p)==0._r8) then
             if (forc_wind(g)==0.0_r8) then
                   ra_grnd_EASS(p)  =300.0_r8
             else

                   !(1) ra_g , air resistnace of ground
                   hdisp=0            ! displacement height [m]
                  ! if (frac_sno(c) > 0._r8) then
                   if (frac_snow_grnd_EASS(c) > 0._r8) then
                       z0 = zsno         ! roughness lengh [m]
                   else
                       z0 = zlnd
                   end if

                   ustar_EASS=forc_wind(g)*vkc/log((forc_hgt_u(g)-hdisp)/z0)    !friction velocity [m s-1]

                   leng_EASS= -1.0_r8 * (vkc*grav*eflx_sh_grnd_EASS(p))/(forc_rho(g)*cpair*forc_t_col(c)*ustar_EASS**3.0_r8)
                   leng_EASS=max(-2.0_r8, leng_EASS)

                   if(leng_EASS .gt. 0.0_r8) then
                         psi=1.0_r8+5.0_r8*(forc_hgt_u(g)-hdisp)*leng_EASS
                   else
                         psi=(1.0_r8-16.0_r8*(forc_hgt_u(g)-hdisp)*leng_EASS)**(-0.5_r8)
                   end if

                   psi = min(10.0_r8, psi)

                   htop_EASS(ivt(p))=0.001_r8
                   hunder_EASS(ivt(p))=0.001_r8

                   kh_o = 0.41_r8*ustar_EASS*(htop_EASS(ivt(p))-htop_EASS(ivt(p))*0.8_r8)/psi   ! m s-1* m=m2 s-1

                   gamm=4.0_r8
                   ra_grnd_EASS(p)=htop_EASS(ivt(p))/(gamm* kh_o)* (exp(gamm)-exp(gamm*(1.0_r8-hunder_EASS(ivt(p))/htop_EASS(ivt(p)))))   ! m/(m2 s-1)=s m-1

             !     ra_grnd_EASS(p)=htop_EASS(ivt(p))/(gamm* kh_o)* (exp(gamm*(1.0_r8-0.0_r8/htop_EASS(ivt(p))))&
             !               -exp(gamm*(1.0_r8-hunder_EASS(ivt(p))/htop_EASS(ivt(p)))))   ! m/(m2 s-1)=s m-1
             end if
             ra_over_EASS(p)  =0.0_r8
             ra_under_EASS(p) =0.0_r8
             rb_over_EASS(p)  =0.0_r8
             rb_under_EASS(p) =0.0_r8
      else

             if (forc_wind(g)==0.0_r8) then
                   ra_over_EASS(p)  =200.0_r8
                   ra_under_EASS(p) =0.0_r8
                   rb_over_EASS(p)  =200.0_r8
                   rb_under_EASS(p) =200.0_r8
                   ra_grnd_EASS(p)  =300.0_r8
             else
                   !(0) Common variable
                   hdisp   = 0.8_r8*htop_EASS(ivt(p))         ! displacement height [m]
                   z0      = 0.08_r8*htop_EASS(ivt(p))        ! roughness lengh [m]

                   !(1) ra_o, air resistnace of overstroy
                   ustar_EASS=forc_wind(g)*vkc/log((forc_hgt_u(g)-hdisp)/z0)    !friction velocity [m s-1]

                   leng_EASS= -1.0_r8 * (vkc*grav*eflx_sh_over_EASS(p))&
                              /(forc_rho(g)*cpair*forc_t_col(c)*ustar_EASS**3.0_r8)

                   !unit:(m s-2 w m-2)/(kg m-3 J kg-1 k-1 K m3 s-3)=(m s-2 J s-1 m-2)/(J s-3)=m-1
                   leng_EASS=max(-2.0_r8, leng_EASS)

                   ra_over_EASS(p)= 1.0_r8/(vkc*ustar_EASS)* (log((forc_hgt_u(g)-hdisp)/z0)&
                                  +(5.0_r8*(forc_hgt_u(g)-hdisp)*leng_EASS)) ! m-1 s + m * m-1=s m-1
                   ra_over_EASS(p)=max(2.0_r8, ra_over_EASS(p))
                   ra_over_EASS(p)=min(100.0_r8, ra_over_EASS(p))                
            
                   !(2) ra_u , air resistnace of understroy 
                   if(leng_EASS .gt. 0.0_r8) then
                         psi=1.0_r8+5.0_r8*(forc_hgt_u(g)-hdisp)*leng_EASS  
                   else
                         psi=(1.0_r8-16.0_r8*(forc_hgt_u(g)-hdisp)*leng_EASS)**(-0.5_r8)
                   end if
                   psi = min(10.0_r8, psi)

                   kh_o = 0.41_r8*ustar_EASS*(htop_EASS(ivt(p))-htop_EASS(ivt(p))*0.8_r8)/psi   ! m s-1* m=m2 s-1
            
                   lt_over_EASS=elai_over_EASS(p)+esai_over_EASS(p)
                   gamm = 0.1_r8 + lt_over_EASS ** 0.75_r8

                   ra_under_EASS(p) = htop_EASS(ivt(p))*(exp(gamm*(1.0_r8-hunder_EASS(ivt(p))/htop_EASS(ivt(p))))-1.0_r8)/(gamm* kh_o)     !m/(m2 s-1)=s m-1
                   ra_under_EASS(p) = ra_over_EASS(p) + ra_under_EASS(p)

                   !(3) rb_o , leaf boundary resistance of overstory 
                   Le=lt_over_EASS*clumping_EASS(p)         ![-]
                   windsp_hcan=1.1_r8*ustar_EASS/vkc        ! wind speed at canopy top[m s-1]

                   gamm=(0.167_r8+0.179_r8*windsp_hcan)* Le**(1.0_r8/3.0_r8) ! ?unit is [m s-1] or [-]
                   !for the managed sand, in eq B.3, ChenBZ et al.,2007, p296 

                   ! Wind speed at displancement height (=0.8*htop(p)) within the canopy, 
                   ! taking as the mean wind speed inside a stand */
                   windsp_hdisp=windsp_hcan*exp(-gamm*( 1.0_r8-hdisp/htop_EASS(ivt(p)) ))

                   ! kinematic viscosity of air[m2 s-1]
                   visco=(13.3_r8+(forc_t_col(c)-tfrz)*0.07_r8)/1.0e6_r8    ! exprimental or statistcal function

                   !Reynold's number
                   Re=(windsp_hdisp*0.1_r8)/visco
                   ! the constant 0.1 is the varialbe u* in eq5.71(p75, CLM4.0 Technical despription)

                   !Nusselt number
                   Nu=1.0_r8*Re**0.5_r8
                   alfaw=(18.9_r8+(forc_t_col(c)-tfrz)*0.07_r8)/1.0e6_r8

                   ! leaf boundary resistance of overstory
                   rb_over_EASS(p) = 0.5_r8*0.1_r8/(alfaw*Nu)

              !#   rb_over_EASS(p) = min(40.0_r8,rb_over_EASS(p))

                   !(4) rb_u , leaf boundary resistance of understory
                   gamm=0.1_r8+lt_over_EASS**0.75_r8

                   ! wind speed at the height of understory= hunder*0.8
                   windsp_hunder=windsp_hcan*exp(-gamm*( 1.0_r8-hunder_EASS(ivt(p))*0.8_r8/htop_EASS(ivt(p)) ))   ![m s-1]

                   ! Reynold's number 
                   Re=(windsp_hunder*0.1_r8)/visco ! [-]=(m s-1)/(s m-1)

                   !Nusselt number
                   Nu=1.0_r8*Re**0.5_r8          ! [-]

                   !leaf boundary resistance of understory
                   rb_under_EASS(p) = 0.5_r8*0.1_r8/(alfaw*Nu)
               !##  rb_under_EASS(p) = min(40.0_r8,rb_under_EASS(p))

                   !(5) ra_g , air resistnace of ground
                   gamm=4.0_r8

                   !kh_u=kh_o*exp(-gamm*(1.0_r8-hunder_EASS(ivt(p))/htop_EASS(ivt(p))))   ![m2 s-1]
                   !????? kh_u indicate?

                   ra_grnd_EASS(p)=htop_EASS(ivt(p))/(gamm* kh_o)* (exp(gamm*(1.0_r8-0.0_r8/htop_EASS(ivt(p))))&
                            -exp(gamm*(1.0_r8-hunder_EASS(ivt(p))/htop_EASS(ivt(p)))))   ! m/(m2 s-1)=s m-1
                   ra_grnd_EASS(p)=ra_grnd_EASS(p)+ra_under_EASS(p)
                   ra_grnd_EASS(p)=max(120.0_r8,ra_grnd_EASS(p))

             end if
      
             ! conductance of leave [m s-1]
             ! heat conductance of sunlit and shaded leaves of overstory and understory    
             Gheat_over_EASS(p)=1.00_r8/(ra_over_EASS(p)+0.50_r8*rb_over_EASS(p))
             Gheat_under_EASS(p)=1.00_r8/(ra_under_EASS(p)+0.50_r8*rb_under_EASS(p))   

             ! intercepted water conductance of sunlit and shaded leaves of overstory and understory
             Gevap_over_EASS(p)=1.00_r8/(ra_over_EASS(p)+rb_over_EASS(p)+100.0_r8)
             Gevap_under_EASS(p)=1.00_r8/(ra_under_EASS(p)+rb_under_EASS(p)+100.0_r8)

             ! leaf conductance of overstory and understory
             Gleaf_over_sunlit_EASS(p)=1.00_r8/(ra_over_EASS(p)+rb_over_EASS(p)+rs_over_sunlit_EASS(p))
             Gleaf_over_shaded_EASS(p)=1.00_r8/(ra_over_EASS(p)+rb_over_EASS(p)+rs_over_shaded_EASS(p))
             Gleaf_under_sunlit_EASS(p)=1.00_r8/(ra_under_EASS(p)+rb_under_EASS(p)+rs_under_sunlit_EASS(p))
             Gleaf_under_shaded_EASS(p)=1.00_r8/(ra_under_EASS(p)+rb_under_EASS(p)+rs_under_shaded_EASS(p))

      end if
!   end do

 end subroutine AirResist_EASS

end module AirResist_EASSMod
