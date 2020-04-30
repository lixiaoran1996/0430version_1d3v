module m_init_cond_1d
   use m_types
   use m_error

   implicit none
   private

   integer, parameter :: INIT_gaussian_type = 1, &
        INIT_block_type = 2, INIT_background_type = 3
   integer            :: INIT_type          = -1
   real(dp)           :: INIT_location
   real(dp)           :: INIT_width
   real(dp)           :: INIT_elec_dens
   real(dp)           :: INIT_elec_energy
   real(dp)           :: INIT_ion_dens
   real(dp)           :: INIT_ion_energy
   real(dp)           :: INIT_elec_temp
   real(dp)           :: INIT_ion_temp
   logical            :: INIT_use_temp

   public :: INIT_init_cfg
   public :: INIT_get_elec_dens
   public :: INIT_get_ion_dens
   public :: INIT_get_en_dens
   public :: INIT_get_ion_en_dens

contains

   subroutine INIT_init_cfg()
      use m_config
      use m_units_constants
      use m_random
      use m_particle_core

      character(len=name_len) :: init_name
      real(dp)                :: domain_length

      integer                 :: INIT_grid_size , number
      real(dp)                :: INIT_grid_delta_x
      integer                 :: num_part, ion_num_part, n, pp
      real(dp) :: xx
      real(dp) :: elec_dens, ion_dens, en_dens, mean_en, ion_en_dens, ion_mean_en
      real(dp) :: pos(3), vel(3), accel(3)
      real(dp) :: ion_pos(3), ion_vel(3), ion_accel(3)
      real(dp) :: weightZeroTime
      real(dp) :: INIT_transverse_area, INIT_grid_volume, INIT_inv_grid_volume, ionMassRatio

      real(dp)    :: INIT_dbd_len(2), temp
      integer     :: nn

      ! Get the type of initial condition
      call CFG_get("init_cond_name", init_name)

      select case (init_name)
      case ("gaussian")
         INIT_type = INIT_gaussian_type
      case ("block")
         INIT_type = INIT_block_type
      case ("background")
         INIT_type = INIT_background_type
      case default
         call ERR_show("INIT_initialize_with_config: unknown initial condition specified")
      end select


      ionMassRatio      = CFG_get_real("part_ionMassRatio")
      weightZeroTime    = CFG_get_real("initial_weight")
      INIT_grid_delta_x = CFG_get_real("grid_delta_x") 
      INIT_grid_size    = CFG_get_int("grid_num_points")

      domain_length = INIT_grid_delta_x * (INIT_grid_size - 1)
      INIT_transverse_area = CFG_get_real("part_transverse_area")
      INIT_grid_volume     = INIT_grid_delta_x * INIT_transverse_area
      INIT_inv_grid_volume = 1 / INIT_grid_volume

      ! parameters for DBD
      do number = 1, 2
          INIT_dbd_len(number) = CFG_get_real("sim_len_del",number)
      end do    

      INIT_use_temp = CFG_get_logic("init_use_temInK")
      INIT_width = CFG_get_real("init_width")
      INIT_location = CFG_get_real("init_rel_pos") * domain_length
      INIT_elec_dens = CFG_get_real("init_elec_dens")

      INIT_elec_temp = CFG_get_real("init_elec_tem")
      INIT_ion_temp  = CFG_get_real("init_ion_tem")

      if (INIT_use_temp) then  ! use temperature as initial conditions or use energy in eV
          INIT_elec_energy = CFG_get_real("init_elec_tem") * UC_boltzmann_const / UC_elec_volt * 3.d0 /2.d0 ! kB * Te  (in eV)
          INIT_ion_energy = CFG_get_real("init_ion_tem") * UC_boltzmann_const / UC_elec_volt  * 3.d0 /2.d0
      else
          INIT_elec_energy = CFG_get_real("init_elec_energy")
          INIT_ion_energy = CFG_get_real("init_ion_energy")
      end if
      INIT_ion_dens = CFG_get_real("init_elec_dens") ! electrons and ions have the same density


! Initialization of electron and ion density
      do n = 1, INIT_grid_size
         xx = (n-1) * INIT_grid_delta_x + INIT_grid_delta_x * (RNG_uniform())
         !Anbang: do not initialize particles out of the domain
         if (xx >= ((INIT_grid_size -1 ) * INIT_grid_delta_x - INIT_dbd_len(2)) .or. (xx <= INIT_dbd_len(1)))  cycle
           

        call INIT_get_elec_dens(xx, elec_dens)
        call INIT_get_ion_dens(xx, ion_dens)
        call INIT_get_en_dens(xx, en_dens)
        call INIT_get_ion_en_dens(xx, ion_en_dens)
        num_part = nint(elec_dens * INIT_grid_volume / weightZeroTime)
        ion_num_part = nint(ion_dens * INIT_grid_volume / weightZeroTime)


         !> For electrons
         if (elec_dens > 0.0_dp) then
            mean_en = en_dens * UC_elec_volt / elec_dens
         else
            mean_en = 0.0_dp
         end if

        !> For ions
         if (ion_dens > 0.0_dp) then
            ion_mean_en = ion_en_dens * UC_elec_volt / ion_dens
         else
            ion_mean_en = 0.0_dp
         end if
         
         if (num_part >= 1) then   ! We have to make sure there is particle when we initialize it
            !> Creat electrons information: velocity and position 
            pos = (/0.0_dp, 0.0_dp, xx/)
            ! Accel is set after computing electric field
            accel = 0.0_dp
            do pp = 1, num_part
                 ! method 1: from jannis
                ! Set a Maxwellian velocity distribution with the correct mean energy
!                 vel = (/RNG_normal(), RNG_normal(), RNG_normal()/)
!                 vel = vel * sqrt(2 * mean_en / (3 * UC_elec_mass))


                ! method 2
            do nn = 1, 3
                temp =  RNG_uniform()
                if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
                vel(nn) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * INIT_elec_temp / UC_elec_mass ) * &
                            & sin(2.d0 * UC_pi * RNG_uniform())
            end do

             call PC_create_part(pos, vel, accel, weightZeroTime)   
            end do
         end if
            
        if (ion_num_part  >=1) then 
            !> Create the information for ions
            ion_pos = (/0.0_dp, 0.0_dp, xx/)
            ion_accel = 0.0_dp 
            do pp = 1, ion_num_part
                ! Set a Maxwellian velocity distribution with the correct mean energy
!                 ion_vel = (/RNG_normal(), RNG_normal(), RNG_normal()/)
!  
!                 ion_vel = ion_vel * sqrt(2 * ion_mean_en / (3 * UC_atomic_mass * ionMassRatio))  

                ! method 2
            do nn = 1, 3
                temp =  RNG_uniform()
                if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
                ion_vel(nn) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * INIT_ion_temp / (UC_atomic_mass * ionMassRatio) ) * &
                            & sin(2.d0 * UC_pi * RNG_uniform())
            end do

                call PC_create_ion_part(ion_pos, ion_vel, ion_accel, weightZeroTime)   
             end do
         end if

      end do
     
      
   end subroutine INIT_init_cfg

   subroutine INIT_get_elec_dens(xx, elec_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: elec_dens

      select case (INIT_type)
      case (INIT_gaussian_type)
         elec_dens = INIT_elec_dens * exp(-(xx - INIT_location)**2 / (2 * INIT_width**2))
         if (elec_dens < 1.0_dp) elec_dens = 0.0_dp
      case (INIT_block_type)
         if (abs(xx - INIT_location) < INIT_width) then
            elec_dens = INIT_elec_dens
         else
            elec_dens = 0.0_dp
         end if
      case (INIT_background_type)
         elec_dens = INIT_elec_dens
      case default
         elec_dens = 0.0_dp
         call ERR_show("m_initial_cond is being used without proper initialization")
      end select

   end subroutine INIT_get_elec_dens

   subroutine INIT_get_ion_dens(xx, ion_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: ion_dens
      call INIT_get_elec_dens(xx, ion_dens)
   end subroutine INIT_get_ion_dens

   subroutine INIT_get_en_dens(xx, en_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: en_dens
      real(dp) :: elec_dens
      call INIT_get_elec_dens(xx, elec_dens)
      en_dens = INIT_elec_energy * elec_dens
   end subroutine INIT_get_en_dens

   subroutine INIT_get_ion_en_dens(xx, ion_en_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: ion_en_dens
      real(dp) :: ion_dens
      call INIT_get_ion_dens(xx, ion_dens)
      ion_en_dens = INIT_ion_energy * ion_dens
   end subroutine INIT_get_ion_en_dens

end module m_init_cond_1d
