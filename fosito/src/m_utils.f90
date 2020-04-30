module m_utils
  use m_types
  
  implicit none
  private

  public :: UT_twoNorm
  public :: UT_binarySearch
  public :: UT_linInterpList
  public :: UT_swapInt
  public :: UT_swapDBLE
  public :: UT_IxFirstTrue

contains

  subroutine UT_swapInt(a, b)
    integer, intent(INOUT) :: a, b
    integer                :: temp
    temp = a
    a = b
    b = temp   
  end subroutine UT_swapInt

  subroutine UT_swapDBLE(a, b)
    real(dp), intent(INOUT) :: a, b
    real(dp)                :: temp
    temp = a
    a = b
    b = temp   
  end subroutine UT_swapDBLE

  integer function UT_ixFirstTrue(bools)
    logical, intent(in) :: bools(:)
    integer             :: n
    
    UT_ixFirstTrue = -1
    do n = 1, size(bools)
       if (bools(n)) then
          UT_ixfirstTrue = n
          exit
       end if
    end do
  end function UT_ixFirstTrue

  ! Searches list for the interval containing value, such that list(i) <= value < list(i+1), and returns i
  integer function UT_binarySearch(list, value)
    real(dp), dimension(:), intent(IN) :: list
    real(dp), intent(IN)               :: value
    integer                                    :: iMin, iMax, iMiddle

    iMin = lbound(list, 1)
    iMax = ubound(list, 1)
    
    do
       iMiddle = iMin + (iMax - iMin) / 2
       if ( value < list(iMiddle) ) then
          iMax = iMiddle         
       else
          iMin = iMiddle
       end if
       if (iMax - iMin <= 1) exit
    end do
    
    UT_binarySearch = iMin
  end function UT_binarySearch

  subroutine UT_linInterpList(xList, yList, xValue, yValue)
    real(dp), intent(IN)    :: xList(:), yList(:)
    real(dp), intent(IN)    :: xValue
    real(dp), intent(INOUT) :: yValue

    integer                         :: i, iMin, iMax
    real(dp)                :: temp
    
    iMin = lbound(xList, 1)
    iMax = ubound(xList, 1)

    if (xValue <= xList(iMin)) then
       yValue = yList(iMin)
    else if (xValue >= xList(iMax)) then
       yValue = yList(iMax)
    else
       i = UT_binarySearch(xList, xValue)
       temp = (xValue - xList(i)) / (xList(i+1) - xList(i))
       yValue = (1.0D0 - temp) * yList(i) + temp * yList(i+1)  
    end if

  end subroutine UT_linInterpList

  ! Computes the two-norm of a vector
  real(dp) function UT_twoNorm(vec)
    real(dp), intent(IN) :: vec(:)
    UT_twoNorm = sqrt( sum(vec**2) )
  end function UT_twoNorm

end module m_utils
