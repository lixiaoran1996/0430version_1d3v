!> Module for the random number generation
module m_random
   use m_types

   implicit none
   private


   INTEGER :: xx=123456789, yy=362436069, zz=521288629, ww=916191069

   DOUBLE PRECISION, PARAMETER :: maxIntPlusOne = DBLE(HUGE(1)) + 1.0D0
   DOUBLE PRECISION, PARAMETER :: toUniformFactor = 0.5D0 / maxIntPlusOne


   public  :: RNG_uniform
   public  :: RNG_normal
   public  :: RNG_poisson
   public  :: RNG_set_seed



contains

   !> Initialize the RNG seed
   subroutine RNG_set_seed(arg_seed)
      integer, intent(IN)  :: arg_seed
      integer              :: n, seed_size
      integer, allocatable :: used_seed(:)

      call random_seed(size = seed_size) ! Get required size of seed
      allocate(used_seed(seed_size))

      do n = 1, seed_size
         used_seed(n) = arg_seed * n ! Fill with some arbritrary values depending on arg_seed
      end do

      call random_seed(put = used_seed)
      deallocate(used_seed)
   end subroutine RNG_set_seed


!    !> Sets the seeds for the KISS RNG
!    SUBROUTINE seedKiss(seeds)
!       INTEGER, INTENT(IN) :: seeds(4)
!       xx = seeds(1)
!       yy = seeds(2)
!       zz = seeds(3)
!       ww = seeds(4)
!    END SUBROUTINE seedKiss

   !> The  KISS (Keep It Simple Stupid) random number generator. Combines:
   !! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
   !! (2) A 3-shift shift-register generator, period 2^32-1,
   !! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
   !!  Overall period>2^123;  Default seeds x,y,z,w.
   !!  Set your own seeds with statement i=kisset(ix,iy,iz,iw).
   INTEGER FUNCTION kiss()
      xx = 69069 * xx + 1327217885
      yy = mm (mm (mm (yy, 13), - 17), 5)
      zz = 18000 * iand (zz, 65535) + ishft (zz, - 16)
      ww = 30903 * iand (ww, 65535) + ishft (ww, - 16)
      kiss = xx + yy + ishft (zz, 16) + ww
   CONTAINS
      !> Small helper function for kiss() RNG.
      INTEGER FUNCTION mm(k, n)
         INTEGER :: k, n
         mm = ieor (k, ishft (k, n) )
      END FUNCTION mm
   END FUNCTION kiss


!> Anbang: The following are two ways to generate the random number

   !> Method 1 
   !> Generate a uniform deviate between 0.0 (inclusive) and 1.0 (exclusive)
   !! RINT / (INT_MAX + 1) lies in [-1, 1), so 0.5 * RINT / (INT_MAX + 1) + 0.5
   !! should be in the interval [0,1)
    DOUBLE PRECISION FUNCTION RNG_uniform()
       RNG_uniform = DBLE(kiss()) * toUniformFactor + 0.5D0
    END FUNCTION RNG_uniform

   !> Method 2
   !> Generate a uniform deviate between 0.0 (inclusive) and 1.0 (exclusive)
!    real(dp) function RNG_uniform()
!       call random_number(RNG_uniform)
!    end function RNG_uniform

   !> Return normal random variate with mean 0 and variance 1
   ! Source: http://en.wikipedia.org/wiki/Marsaglia_polar_method
   real(dp) function RNG_normal()
      logical, save :: have_spare = .false.
      real(dp), save :: spare
      real(dp) :: u, v, s, mul

      if (have_spare) then
         have_spare = .false.
         RNG_normal = spare
      else
         do
            u = RNG_uniform() * 2 - 1
            v = RNG_uniform() * 2 - 1
            s = u * u + v * v
            if (s < 1.0_dp .and. s /= 0.0_dp) exit
         end do

         mul        = sqrt(-2 * log(s) / s)
         RNG_normal = u * mul
         spare      = v * mul
         have_spare = .true.
      end if
   end function RNG_normal

   !> Return Poisson random variate with rate labda. Works well for labda < 30 or so.
   !! For labda >> 1 it can produce wrong results due to roundoff error.
   integer function RNG_poisson(labda)
      real(dp), intent(IN) :: labda
      integer :: k
      real(dp) :: expL, p

      expL = exp(-labda)
      k = 0
      p = RNG_uniform()

      do while (p > expL)
         k = k + 1
         p = p * RNG_uniform()
      end do
      RNG_poisson = k

   end function RNG_poisson

end module m_random
