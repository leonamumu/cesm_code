module Tridiagonal_EASSMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Tridiagonal_EASSMod
!
! !DESCRIPTION:
! Tridiagonal matrix solution
!
! !USE:
  use clm_varctl,   only: iulog

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Tridiagonal_EASS
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tridiagonal_EASS
!
! !INTERFACE:
  subroutine Tridiagonal_EASS ( lbj, ubj, jtop, ci, a, &
                                b, c, r, u)
!
! !DESCRIPTION:
! Tridiagonal matrix solution
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
  !#  integer , intent(in)    :: fc                !  column indices
    integer , intent(in)    :: ci                ! filter
    integer , intent(in)    :: lbj, ubj          ! lbinning and ubing level indices
    integer , intent(in)    :: jtop              ! top level for each column
    real(r8), intent(in)    :: a(lbj:ubj)    ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(lbj:ubj)    ! "b" diagonal column for tridiagonal matrix
    real(r8), intent(in)    :: c(lbj:ubj)    ! "c" right off diagonal tridiagonal matrix
    real(r8), intent(in)    :: r(lbj:ubj)    ! "r" forcing term of tridiagonal matrix
    real(r8), intent(inout) :: u(lbj:ubj)    ! solution
!
! !CALLED FROM:
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine SoilTemperature in module SoilTemperatureMod
! subroutine SoilWater in module HydrologyMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!  1 July 2003: Mariana Vertenstein; modified for vectorization
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: j  !#,ci,fc                   !indices
    real(r8) :: gam(lbj:ubj)      !temporary
    real(r8) :: bet              !temporary
!-----------------------------------------------------------------------

    ! Solve the matrix

    bet = b(jtop)

    do j = lbj, ubj
          if (j >= jtop) then
             if (j == jtop) then
                u(j) = r(j) / bet
             else
                gam(j) = c(j-1) / bet
                bet = b(j) - a(j) * gam(j)
                u(j) = (r(j) - a(j)*u(j-1)) / bet
             end if
          end if
    end do

    do j = ubj-1,lbj,-1
          if (j >= jtop) then
             u(j) = u(j) - gam(j+1) * u(j+1)
          end if
    end do

  end subroutine Tridiagonal_EASS

end module Tridiagonal_EASSMod
