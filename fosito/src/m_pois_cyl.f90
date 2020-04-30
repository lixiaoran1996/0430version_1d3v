module m_pois_cyl
   use m_types

   implicit none
   private

   public :: PSC_solve_dirichlet
   public :: PSC_solve_free
   public :: PSC_ring_sol

contains

   ! Solve laplace phi = rhs (Poisson's eq.) for a given source distribution with Dirichlet boundary conditions.
   ! rhs(i,j) contains the source terms at r = (i-0.5)*deltas(1) and z = (j-0.5)*deltas(2).
   ! rhs(i_max, :), rhs(:, 1), rhs(:, j_max) contains the boundary values.
   subroutine PSC_solve_dirichlet(rhs, deltas)
      use m_error
      real(dp), intent(inout)     :: rhs(:, :)
      real(dp), intent(in)        :: deltas(2)

      integer                     :: N_rz(2), ierror, n_work
      integer, save               :: n_work_allocated = 0
      real(dp)                    :: lengths(2), dummy(1), pertrb,  zero = 0.0_dp
      real(dp), allocatable, save :: workspace(:)

      interface
         subroutine hwscyl (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
              elmbda, f, idimf, pertrb, ierror, w)
            import
            real(dp) :: a, b, bda(*), bdb(*), c, d, bdc(*), bdd(*), elmbda, f(*), pertrb, w(*)
            integer  :: m, mbdcnd, n, nbdcnd, idimf, ierror
         end subroutine hwscyl
      end interface

      N_rz    = shape(rhs)
      lengths = (N_rz-1) * deltas
      n_work  = 4*(N_rz(2)+1) + (13 + INT(LOG(N_rz(2)+1.0_dp)/log(2.0_dp)))*(N_rz(1)+1)

      if (n_work > n_work_allocated) then
         if (allocated(workspace)) deallocate(workspace)
         allocate(workspace(n_work))
         n_work_allocated = n_work
      end if

      call hwscyl(zero, lengths(1), N_rz(1)-1, 5, dummy, dummy, &
           zero, lengths(2), N_rz(2)-1, 1, dummy, dummy, &
           zero, rhs, N_rz(1), pertrb, ierror, workspace)

      if (ierror /= 0) call ERR_show("PSC_solve_dirichlet error in calling hwscyl")
   end subroutine PSC_solve_dirichlet

   ! Solve laplace phi = rhs (Poisson's eq.) for a given source distribution with free space boundary conditions.
   ! rhs(i,j) contains the source terms at r = (i-0.5)*deltas(1) and z = (j-0.5)*deltas(2).
   subroutine PSC_solve_free(rhs, deltas, arg_n_buff)
      real(dp), intent(inout)       :: rhs(:, :)
      real(dp), intent(in)          :: deltas(2)
      integer, intent(in), optional :: arg_n_buff

      integer                       :: i, j, n, n_buff, N_bnd, M_bnd
      integer                       :: N_r, N_z, M_r, M_z, i_zmin, i_zmax, j_zmin, j_zmax
      real(dp)                      :: H_domain, H_buff, temp

      ! If all of these stay the same between calls, we can avoid doing extra work
      real(dp), allocatable, save   :: bnd_effect(:,:), rhs_buf(:,:), bnd_rho(:), rr(:), zz(:), ww(:)
      logical, save                 :: first_time = .true.
      integer, save                 :: prev_shape(2) = 0, prev_n_buff = 0
      real(dp), save                :: prev_deltas(2) = 0.0_dp

      N_r       = size(rhs, 1)
      N_z       = size(rhs, 2)
      if (present(arg_n_buff)) then
         n_buff = arg_n_buff
      else
         n_buff = max(8, nint(sqrt(real(max(N_r, N_z), dp))))
      end if
      M_r       = N_r + n_buff
      M_z       = N_z + 2 * n_buff
      N_bnd     = 2 * N_r + N_z
      M_bnd     = 2 * M_r + M_z
      H_domain  = (N_z-1) * deltas(2)
      H_buff    = n_buff * deltas(2)
      i_zmin    = N_z
      i_zmax    = i_zmin + N_r
      j_zmin    = M_z
      j_zmax    = j_zmin + M_r

      if (first_time .or. any((/N_r, N_z/) /= prev_shape) .or. &
           any(deltas /= prev_deltas) .or. n_buff /= prev_n_buff) then

         ! Set up storage for the charge on the inner boundary, and a matrix to evaluate this on the outer boundary
         if (.not. first_time) then
            deallocate(rhs_buf)
            deallocate(bnd_effect)
            deallocate(bnd_rho)
            deallocate(rr)
            deallocate(zz)
            deallocate(ww)
         end if

         prev_shape  = (/N_r, N_z/)
         prev_deltas = deltas
         prev_n_buff = n_buff
         first_time  = .false.

         allocate(rhs_buf(M_r, M_z))
         allocate(bnd_effect(N_bnd, M_bnd))
         allocate(bnd_rho(N_bnd))
         allocate(rr(M_r))
         allocate(zz(M_z))
         allocate(ww(M_r))

         ! Set up the 'interaction' matrix
         do n = 1, M_r
            rr(n) = (n-1) * deltas(1)
            ww(n) = product(deltas)
         end do

         ! The cells at r = 0 have half the width, and their center is at dr/4
         rr(1) = deltas(1) * 0.25_dp
         ww(1) = ww(1) * 0.5_dp

         do n = 1, M_z
            zz(n) = (n-1) * deltas(2)
         end do


         ! Set values for r_max at the outside grid, we can get half from symmetry
         bnd_effect = 0.0_dp

         do j = 1, (M_z+1)/2                                                                ! outside grid: r_max
            do i = 1, N_z                                                                   ! inside grid: r_max
               bnd_effect(i,j)              = PSC_ring_sol(rr(N_r), rr(M_r), ((j-n_buff-i)*deltas(2))**2, ww(N_r))
               bnd_effect(N_z-i+1, M_z-j+1) = bnd_effect(i,j)
            end do

            do i = 1, N_r                                                                   ! inside grid: z_min and z_max
               bnd_effect(i_zmin+i,j)        = PSC_ring_sol(rr(i), rr(M_r), (zz(j)-H_buff)**2, ww(i))
               bnd_effect(i_zmax+i,j)        = PSC_ring_sol(rr(i),  rr(M_r), (H_domain+H_buff-zz(j))**2, ww(i))
               bnd_effect(i_zmin+i, M_z-j+1) = bnd_effect(i_zmax+i, j)
               bnd_effect(i_zmax+i, M_z-j+1) = bnd_effect(i_zmin+i, j)
            end do
         end do

         ! Set values for the z_min and z_max at the outside grid. We have to compute only half, because of the symmetry.
         do j = 1, M_r                                                                                 ! outside grid: z_min and z_max
            do i = 1, N_z
               bnd_effect(i,j_zmin+j)       = PSC_ring_sol(rr(N_r), rr(j), (zz(i)+H_buff)**2, ww(N_r)) ! ins: r_max, out: z_min
               bnd_effect(N_z-i+1,j_zmax+j) = bnd_effect(i,j_zmin+j)
            end do

            do i = 1, N_r
               bnd_effect(i_zmin+i,j_zmin+j) = PSC_ring_sol(rr(i), rr(j), H_buff**2, ww(i)) ! ins: z_min, out: z_min
               bnd_effect(i_zmax+i,j_zmin+j) = PSC_ring_sol(rr(i), rr(j), (H_domain+H_buff)**2, ww(i))
               bnd_effect(i_zmin+i,j_zmax+j) = bnd_effect(i_zmax+i,j_zmin+j)                ! ins: z_min, out: z_max
               bnd_effect(i_zmax+i,j_zmax+j) = bnd_effect(i_zmin+i,j_zmin+j)                ! ins: z_max, out: z_max
            end do
         end do

      end if

      ! Store boundary charge and set boundary conditions to zero at boundary
      bnd_rho(1:N_z)           = rhs(N_r,:)
      rhs(N_r,:)               = 0.0_dp
      bnd_rho(i_zmin+1:i_zmax) = rhs(:,1)
      rhs(:, 1)                = 0.0_dp
      bnd_rho(i_zmax+1:N_bnd)  = rhs(:,N_z)
      rhs(:, N_z)              = 0.0_dp

      ! Store a copy of rhs and compute the current approximation to the solution
      rhs_buf(:,:)                        = 0.0_dp
      rhs_buf(1:N_r, n_buff+1:n_buff+N_z) = rhs

      ! Solve system with homog. dirichlet b.c.
      call PSC_solve_dirichlet(rhs, deltas)

      ! Compute minus the 'induced' charge at the boundary
      do n = 1, N_z           ! r_max
         temp = rhs(N_r-1, n)
         bnd_rho(n) = bnd_rho(n) - temp * (rr(N_r) - 0.5_dp * deltas(1))/(rr(N_r)*deltas(1)**2)
      end do

      do n = 1, N_r           ! z_min and z_max
         bnd_rho(i_zmin+n) = bnd_rho(i_zmin+n) - rhs(n, 2) / deltas(2)**2
         bnd_rho(i_zmax+n) = bnd_rho(i_zmax+n) - rhs(n, N_z-1) / deltas(2)**2
      end do

      ! Set the correction potential at the boundary of the outside grid
      do n = 1, M_z
         rhs_buf(M_r, n) = sum(bnd_effect(:,n) * bnd_rho)
      end do

      do n = 1, M_r
         rhs_buf(n, 1)   = sum(bnd_effect(:,j_zmin+n) * bnd_rho)
         rhs_buf(n, M_z) = sum(bnd_effect(:,j_zmax+n) * bnd_rho)
      end do

      ! Solve with the corrected boundary condition (rhs_buf contains the original rhs, see above)
      call PSC_solve_dirichlet(rhs_buf, deltas)

      ! Copy out the interior part of the solution
      rhs = rhs_buf(1:N_r, n_buff+1:n_buff+N_z)

   end subroutine PSC_solve_free

   real(dp) function PSC_ring_sol(r_ring, r_pos, z_dist2, lambda)
      real(dp), intent(in) :: r_ring, r_pos, z_dist2, lambda
      real(dp), parameter  :: inv_pi = 1.0_dp / acos(-1.0_dp)
      real(dp)             :: temp, k_arg

      temp         = 1 / ((r_ring + r_pos)**2 + z_dist2)
      k_arg        = 1 - (4 * r_ring * r_pos * temp)
      PSC_ring_sol = -inv_pi * r_ring * lambda * sqrt(temp) * drf(0.0_dp, k_arg, 1.0_dp)
   end function PSC_ring_sol

   ! Compute the incomplete or complete elliptic integral of the 1st kind
   ! See drf.f on Netlib.org. Error handling has been removed.
   real(dp) function DRF (xx, yy, zz)
      real(dp), intent(in) :: xx, yy, zz

      real(dp), parameter :: errtol = (4 * epsilon(1.0_dp))**(1/6.0_dp)
      real(dp), parameter :: c1 = 1/24.0_dp, c2 = 3/44.0_dp, c3 = 1/14.0_dp
      real(dp), parameter :: one_third = 1 / 3.0_dp
      real(dp) :: e2, e3, mu, lamda, epslon
      real(dp) :: xyz(3), xyz_dev(3), xyz_root(3)

      xyz = (/xx, yy, zz/)

      do
         mu       = sum(xyz) * one_third
         xyz_dev  = 2 - (mu + xyz)/mu
         epslon   = maxval(abs(xyz_dev))
         if (epslon < errtol) exit

         xyz_root = sqrt(xyz)
         lamda    = xyz_root(1) * (xyz_root(2)+xyz_root(3)) + xyz_root(2) * xyz_root(3)
         xyz      = (xyz + lamda) * 0.25_dp
      end do

      e2  = xyz_dev(1) * xyz_dev(2) - xyz_dev(3) * xyz_dev(3)
      e3  = xyz_dev(1) * xyz_dev(2) * xyz_dev(3)
      drf = 1 + (c1 * e2 - 0.1_dp - c2 * e3) * e2 + c3 * e3
      drf = drf/sqrt(mu)
   end function DRF


end module m_pois_cyl
