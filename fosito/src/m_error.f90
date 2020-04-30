module m_error
   implicit none
   private

   public :: ERR_warn
   public :: ERR_show

contains


   subroutine ERR_warn(msg)
      character(len=*), intent(in) :: msg
      print *, "The following WARNING has been issued"
      print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      print *, trim(msg)
      print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
   end subroutine ERR_warn

   subroutine ERR_show(msg)
      character(len=*), intent(in) :: msg
      integer :: set_me_to_zero
      print *, "The following ERROR has occurred"
      print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      print *, trim(msg)
      print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      print *, "Will halt now..."
      print *, "Also, let's try to give you a backtrace. Enter '0':"
      read(*,*) set_me_to_zero
      print *, 1/set_me_to_zero
      print *, "To see the backtrace, enable signalling on division by zero"
      ! Check whether this can be replaced by a real call,
      ! http://gcc.gnu.org/ml/fortran/2012-12/msg00091.html
      stop
   end subroutine ERR_show

end module m_error
