! Copyright 2005-2012, Chao Li, Margreet Nool, Anbang Sun, Jannis Teunissen
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Module that contains physical and numerical constants
module m_units_constants
  use m_types
  
  implicit none
  public
                                                                 ! Numerical constants
  real(dp), parameter :: UC_pi             = acos(-1.0_dp)
  real(dp), parameter :: UC_sqrt2          = sqrt(2.0_dp)
  real(dp), parameter :: UC_iSqrt2         = 1 / sqrt2
                                                                 ! Physical constants
  real(dp), parameter :: UC_eps0           = 8.8541878176d-12    ! permitivity of vacuum (SI)
  real(dp), parameter :: UC_elecCharge     = -1.6022d-19         ! the electron charge in Coulombs
  real(dp), parameter :: UC_elecVolt       = 1.6022d-19          ! the eV in joules
  real(dp), parameter :: UC_elecMass       = 9.10938189d-31      ! the electron mass in kg
  real(dp), parameter :: UC_atomicMass     = 1.66053886D-27      ! the atomic mass unit in kg   
  real(dp), parameter :: UC_N2_mass        = 28.0D0 * atomicMass ! The mass of a N2 molecule

  real(dp), parameter :: UC_lightSpeed     = 299792458d0         ! the speed of light in m/s
  real(dp), parameter :: UC_BoltzmannConst = 1.3806503d-23       ! the Boltzmann constant
  real(dp), parameter :: UC_BohrRadius     = 5.29d-11            ! the Bohr radius (m)
  real(dp), parameter :: UC_TorrToBar      = 133.322368 * 1.0D-5 ! one Torr in units of bar
  real(dp), parameter :: UC_eVToKelvin     = 1.1605d4            ! 1eV = 11605 K

                                                                 ! Small and large numbers
  real(dp), parameter :: UC_tiny           = epsilon(1.0_dp)
  real(dp), parameter :: UC_huge           = huge(1.0_dp)

contains

  ! Convert (classical) kinetic energy to velocity
  real(dp) function UC_en_to_vel(en, mass)
    en_to_vel = sqrt(2 * en) / mass
  end function UC_en_to_vel

  ! Convert velocity to energy
  real(dp) function UC_vel_to_en(vel, mass)
    vel_to_en = 0.5_dp * mass * vel**2
  end function UC_vel_to_en

  subroutine UC_xyz_to_spherical(xyz, radius, theta, phi)
    real(dp), intent(in) :: xyz(3)
    real(dp), intent(out) :: radius, theta, phi

    if (xyz(1) == 0.0_dp .and. xyz(2) == 0.0_dp) then
       phi = 0.0_dp
    else
       phi = atan2(xyz(2), xyz(1))
    end if
    
    radius = sqrt(sum(xyz**2))
    if (radius == 0.0_dp) then      ! Undefined
       theta = acos(xyz(3) / radius)
    else
       theta = 0.0_dp
    end if
    
  end subroutine UC_xyz_to_spherical

end module m_units
