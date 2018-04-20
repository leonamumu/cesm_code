!define the modle
module CMLForwardMovingAverageMod
!LAST MODIFICATION DATE: 2013-11-27
!BY CHE Mingliang
implicit none
type,public:: ForwardMovingAverage
   real,dimension(:),pointer::array ! define a dynamic pointer array
end type ForwardMovingAverage

type(ForwardMovingAverage),public,target::fma

end module CMLForwardMovingAverageMod

