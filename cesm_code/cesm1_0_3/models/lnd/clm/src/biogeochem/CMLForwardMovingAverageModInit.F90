! initialize the module
module CMLForwardMovingAverageModInit
!LAST MODIFICATION DATE: 2013-11-27
!BY CHE Mingliang
use CMLForwardMovingAverageMod
implicit none
contains

subroutine init()
    integer::n=15
    integer::i
    print *,"Allocating..."
    allocate(fma%array(n))
    !initialzing array
    do i=1,n
       fma%array(i)=0
    end do
end subroutine init

subroutine release()
    print *,"Deallocating..."
end subroutine release

end module CMLForwardMovingAverageModInit
