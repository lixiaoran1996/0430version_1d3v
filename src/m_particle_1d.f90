module m_particle_1d
   use m_types
   use m_particle_core

   implicit none
   private

   integer, parameter    :: PM_num_vars = 11
   integer, parameter    :: PM_iv_elec  = 1, PM_iv_ion = 2, PM_iv_en = 3, PM_iv_ion_en = 4, &
                            & PM_iv_vel_elec = 5, PM_iv_vel_ion = 6, &
                            & PM_directed_vel_elec_x = 7, &
                            & PM_directed_vel_elec_y = 8, &
                            & PM_directed_vel_elec_z = 9, &
                            & PM_directed_en_elec    = 10, &
                            & PM_ioni_sour_term     = 11
                       

   integer               :: PM_grid_size
   real(dp), allocatable :: PM_vars(:,:)
   real(dp)              :: PM_delta_x, PM_inv_delta_x
   real(dp)              :: PM_transverse_area
   real(dp)              :: PM_grid_volume, PM_inv_grid_volume
   real(dp)              :: PM_vel_rel_weight, PM_part_per_cell
   
   real(dp)              :: ionMassRatio  ! the mass ratio between ions and standard atoms
   real(dp)              :: atomMassRatio ! the mass ratio between simulated atoms and standard atoms

   real(dp), allocatable :: PM_pos_aver_data(:,:)  ! the average data in a certain time, 
                                                !Anbang:  for comparison with turner's resutls


   real(dp), allocatable :: PM_data_per_peroid(:,:,:)  !here include all the data in one peroid after peroidic steady state
   integer               :: PM_num_per_peroid

   ! Anbang: look for the position of sheath
   real(dp), allocatable :: PM_pos_sheath1(:,:)  ! this is calculated by (ne/ni = 1/2)
   real(dp), allocatable :: PM_pos_sheath2(:,:)  ! this is calculated by (ni - ne)/ ni = 0.01


   ! Anbang: output the average physical data, likes ni, energy of electron in the middle plane and etc
   real(dp), allocatable :: PM_phys_aver_data(:) 
   real(dp)               :: PM_step_num_averger    ! the step nums for avergering the values

!    ! Anbang, here we define the parameters for DBD
   logical               :: PM_useDBD
   real(dp),allocatable  :: PM_delLen(:)

   real(dp)              :: PM_surChargeAtnodes(2), PM_surChargeAtroot(2)

    ! Anbang: here we define the averged number of particles and mean weight
    real(dp), allocatable :: PM_part_num_per_cell(:,:)  ! for electrons and ions
    real(dp), allocatable :: PM_part_weight_per_cell(:,:)
  
  ! parameter for field emission 
   real(dp)              :: PM_field_factor, PM_work_fuction 
  
  ! Public routines
   public :: PM_initialize
   public :: PM_advance
   public :: PM_get_output
   public :: PM_get_coeffs

   public :: PM_get_num_sim_part, PM_get_num_real_part
   public :: PM_get_num_sim_ion_part, PM_get_num_real_ion_part

   public :: PM_ion_advance

   public :: PM_get_aver_para
   public :: PM_output_aver_para
   public :: PM_get_eedf
   public :: PM_output_aver_phys_para
   public :: PM_correct_accel_and_vel
   public :: PM_get_accel

   public :: PM_set_elec_density
   public :: PM_set_ion_density
   public :: PM_update_efield

 ! parameters for MPI
   public :: PM_collectDensitiesAtRoot, PM_shareDens
   public :: PM_collectEEDFAtRoot
   public :: PM_shareEEDF
   public :: PM_adapt_weights, PM_adapt_ion_weights
   public :: PM_get_surface_charge_density, PM_collectsurfChargeAtRoot

   public :: PM_save_surf_charge_info, PM_read_surf_charge_info

   public :: PM_aver_para_one_peroid, PM_output_aver_para_one_peroid, PM_sheath_pos_cal_output

   public   :: PM_get_aver_num_and_weight_per_cell
   public   :: PM_output_aver_num_and_weight_per_cell
   public   :: PM_calculate_curr
   public   :: PM_get_max_dens
      ! parameter for field emission 
   public   :: PM_field_emission

contains
	subroutine PM_field_emission(tt)
        use m_config
        use m_units_constants
        use m_random
        use m_efield_1d
        real(dp)     :: dbd_field(2), emission_current(2)
        real(dp)     :: variable_y(2),t(2),v(2),enhance_field(2),number_elec(2)
        real(dp),intent(in)    :: tt
        real(dp)     :: weight_initial,xx,pos(3),accel(3),vel(3),elec_temp,temp,old_vel
        integer      :: pp
        elec_temp = CFG_get_real("init_elec_tem")
        weight_initial= CFG_get_real("initial_weight")
        PM_field_factor = CFG_get_real("field factor")
        PM_work_fuction = CFG_get_real("work fuction")

        call EF_get_values_dbd(dbd_field)

        if (abs(dbd_field(1))>1.0d9) then 
		!Es=Belta*E
			enhance_field(1) = PM_field_factor * dbd_field(1)
			!y=3.79d-5*sqrt(Es)/PM_work_fuction
			variable_y(1) = 3.79d-5 * sqrt(enhance_field(1)) / PM_work_fuction	
			!t2=1-y^2*(1-Ln(y)/3) v=1+y^2*(1-Ln(y)/9)
			t(1) = 1 - variable_y(1)**2 * (1 - dlog(variable_y(1)) / 3)
			v(1) = 1 + variable_y(1)**2 * (1 - dlog(variable_y(1)) / 9)
			!Jfe calculate
			emission_current(1) = 1.51414d-6 * (enhance_field(1)**2 / PM_work_fuction*t(1))* &
			exp(-6.8309d9 * PM_work_fuction**(1.5_dp) * v(1) / enhance_field(1))
			!handle emission electron  N=J*S*dt/q now N is real paticles number
			number_elec(1) = emission_current(1) * PM_transverse_area * tt / UC_elem_charge
			!due to electron emission,equal to that suface charge positive increase
			PM_surChargeAtroot(1) = PM_surChargeAtroot(1) + number_elec(1)
			!calculate mcroparticles ,i.e.N=N/initial/weight, assume new particle weight is initial weight
			number_elec(1) = number_elec(1) / weight_initial

			if (number_elec(1)>= 1) then   ! We have to make sure there is particle when we create it
				do pp = 1, nint(number_elec(1))!rounding and change to integer in order to use do loop
				!> Creat electrons information: velocity and position, acceleration
					if(mod(pp,2)/=0) then
						xx=PM_dellen(1)+RNG_uniform()*PM_delta_x
						pos = (/0.0_dp, 0.0_dp, xx/)
						accel(1:2) = 0.0_dp
						accel(3) = UC_elec_charge / UC_elec_mass * EF_get_at_pos(xx)
						! maxwell distribution to describe the particles velocity
						temp =  RNG_uniform()
						if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
						vel(1) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& sin(2.d0 * UC_pi * RNG_uniform())
						vel(2) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& cos(2.d0 * UC_pi * RNG_uniform())
						temp =  RNG_uniform()
						if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
						vel(3)=sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& sin(2.d0 * UC_pi * RNG_uniform())
						old_vel=sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& cos(2.d0 * UC_pi * RNG_uniform())
						call PC_create_part(pos, vel, accel, weight_initial)
					endif
					if(mod(pp,2)==0) then
					! initialize second macroparticle
						xx=PM_dellen(1)+RNG_uniform()*PM_delta_x
						pos = (/0.0_dp, 0.0_dp, xx/)
						accel(1:2) = 0.0_dp
						accel(3) = UC_elec_charge / UC_elec_mass * EF_get_at_pos(xx)
						vel(1) = old_vel
						temp =  RNG_uniform()
						if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
						vel(2) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& sin(2.d0 * UC_pi * RNG_uniform())
						vel(3) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& cos(2.d0 * UC_pi * RNG_uniform())
						call PC_create_part(pos, vel, accel, weight_initial)  
					endif 
				end do
			endif	
		end if

        if ( abs(dbd_field(2))>1.0d9) then 
			enhance_field(2)=PM_field_factor*dbd_field(2) 
			variable_y(2)=3.79d-5*sqrt(enhance_field(2))/PM_work_fuction
			t(2)=1-variable_y(2)**2*(1-dlog(variable_y(2))/3)
			v(2)=1+variable_y(2)**2*(1-dlog(variable_y(2))/9)
			emission_current(2)=1.51414d-6*enhance_field(2)**2/(PM_work_fuction*t(2))*&
			exp(-6.8309d9*PM_work_fuction**1.5*v(2)/enhance_field(2))
			number_elec(2)=emission_current(2)*PM_transverse_area*tt/UC_elem_charge
			PM_surChargeAtroot(2) = PM_surChargeAtroot(2) + number_elec(2)
			number_elec(2)=number_elec(2)/weight_initial
			if (number_elec(2)>= 1) then   ! We have to make sure there is particle when we create it
				!> Creat electrons information: velocity and position
				do pp = 1, nint(number_elec(1))!rounding and change to integer in order to use do loop
				!> Creat electrons information: velocity and position, acceleration
					if(mod(pp,2)/=0) then 
							xx=PM_delta_x * (PM_grid_size - 1)-PM_dellen(2)-RNG_uniform()*PM_delta_x
							pos = (/0.0_dp, 0.0_dp, xx/)
						accel(1:2) = 0.0_dp
							accel(3) = UC_elec_charge / UC_elec_mass * EF_get_at_pos(xx)
							! maxwell distribution to describe the particles velocity
							temp =  RNG_uniform()
						if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
						vel(1) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& sin(2.d0 * UC_pi * RNG_uniform())
							vel(2) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& cos(2.d0 * UC_pi * RNG_uniform())
							temp =  RNG_uniform()
						if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
							vel(3)=sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& sin(2.d0 * UC_pi * RNG_uniform())
						old_vel=sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& cos(2.d0 * UC_pi * RNG_uniform())
							call PC_create_part(pos, vel, accel, weight_initial)
					endif
					if(mod(pp,2)==0) then
					! initialize second macroparticle
						xx=PM_delta_x * (PM_grid_size - 1)-PM_dellen(2)-RNG_uniform()*PM_delta_x
						pos = (/0.0_dp, 0.0_dp, xx/)
							accel(1:2) = 0.0_dp
								accel(3) = UC_elec_charge / UC_elec_mass * EF_get_at_pos(xx)
						vel(1) = old_vel
						temp =  RNG_uniform()
						if (temp == 0.d0) temp = temp + UC_tiny  ! it cannot be zero
						vel(2) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& sin(2.d0 * UC_pi * RNG_uniform())
									vel(3) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * elec_temp / UC_elec_mass ) * &
						& cos(2.d0 * UC_pi * RNG_uniform())
							call PC_create_part(pos, vel, accel, weight_initial)  
					endif
				end do 
			endif
		endif
    end subroutine PM_field_emission

   subroutine PM_initialize(cross_secs, cross_secs_for_ions, ntasks)
      use m_config
      use m_random
      use m_init_cond_1d
      use m_units_constants
      use m_cross_sec
      use m_cross_sec_for_ions

      type(CS_type), intent(in) :: cross_secs(:)
      type(CS_ion_type), intent(in) :: cross_secs_for_ions(:)

      integer ::   varSize, number
      
      integer,allocatable  :: bd_elec_type(:), bd_ion_type(:)
      real(dp),allocatable :: coSecEl(:), coElecRefl(:)
      integer, intent(in)  :: ntasks

      PM_grid_size       = CFG_get_int("grid_num_points")
      PM_delta_x         = CFG_get_real("grid_delta_x")
      PM_inv_delta_x     = 1 / PM_delta_x
      PM_transverse_area = CFG_get_real("part_transverse_area")
      PM_grid_volume     = PM_delta_x * PM_transverse_area
      PM_inv_grid_volume = 1 / PM_grid_volume
      PM_part_per_cell = CFG_get_real("apm_part_per_cell")
      PM_vel_rel_weight = CFG_get_real("apm_vel_rel_weight")
      
      ionMassRatio      = CFG_get_real("part_ionMassRatio")
      atomMassRatio     = CFG_get_real("part_atomMassRatio")
      PM_useDBD         = CFG_get_logic("sim_DBD")

      allocate(PM_vars(PM_grid_size, PM_num_vars))
      PM_vars = 0.0_dp

      allocate(PM_pos_aver_data(PM_grid_size, 13))
      allocate(PM_phys_aver_data(6))
      PM_pos_aver_data = 0.0_dp
      PM_phys_aver_data =0.0_dp
      PM_step_num_averger = CFG_get_int("turner_average_stepNum") * CFG_get_real("sim_conv_fac")

      PM_num_per_peroid = CFG_get_int("sim_div_num_one_peroid")
      allocate(PM_data_per_peroid(PM_num_per_peroid, PM_grid_size, 15))
      PM_data_per_peroid = 0.d0

      allocate(PM_pos_sheath1(PM_num_per_peroid,2))
      allocate(PM_pos_sheath2(PM_num_per_peroid,2))
      PM_pos_sheath1 = 0.d0
      PM_pos_sheath2 = 0.d0

      ! the length of the dielectric length on both sides
      varSize = CFG_get_size("sim_len_del")
      allocate(PM_delLen(varSize))
      do number = 1, varSize
          PM_delLen(number) = CFG_get_real("sim_len_del",number)
      end do    

      ! the type of boundary conditions of electrons on both sides
      varSize = CFG_get_size("sim_bd_elec")
      allocate(bd_elec_type(varSize))
      do number = 1, varSize
          bd_elec_type(number) = CFG_get_int("sim_bd_elec",number)
      end do    

      ! the type of boundary conditions of ions on both sides
      varSize = CFG_get_size("sim_bd_ion")
      allocate(bd_ion_type(varSize))
      do number = 1, varSize
          bd_ion_type(number) = CFG_get_int("sim_bd_ion",number)
      end do    

      ! the secondary electron emission coefficients on the both sides
      varSize = CFG_get_size("sim_co_secondary_elec_emmison")
      allocate(coSecEl(varSize))
      do number = 1, varSize
          coSecEl(number) = CFG_get_real("sim_co_secondary_elec_emmison",number)
      end do    

      ! the reflection coefficients of electrons on the both sides
      varSize = CFG_get_size("sim_ref_coeff_elec")
      allocate(coElecRefl(varSize))
      do number = 1, varSize
          coElecRefl(number) = CFG_get_real("sim_ref_coeff_elec",number)
      end do    

      ! the surface charge at dielectrics
      PM_surChargeAtnodes = 0.d0
      PM_surChargeAtroot  = 0.d0

    ! particle number and mean weight per cell
      allocate(PM_part_num_per_cell(2,PM_grid_size -1))
      allocate(PM_part_weight_per_cell(2,PM_grid_size - 1))
      PM_part_num_per_cell = 0.d0
      PM_part_weight_per_cell = 0.d0


      call PC_initialize(UC_elec_mass, UC_atomic_mass * ionMassRatio, CFG_get_int("part_max_number_of") / ntasks, &
           CFG_get_int("ion_part_max_number_of") / ntasks, &
           cross_secs, cross_secs_for_ions, CFG_get_int("part_lkptbl_size"), CFG_get_real("part_max_energy_eV"), &
           PM_delta_x, PM_grid_size, CFG_get_real("gas_temperature"), atomMassRatio * UC_atomic_mass, &
           CFG_get_logic("sim_useBgVel"), CFG_get_int("sim_frame_type"), PM_delLen, bd_elec_type, bd_ion_type,coSecEl, &
            CFG_get_real("sim_secElecEn"), CFG_get_int("sim_rot_scheme"), CFG_get_int("sim_table_create_scheme"), &
            coElecRefl, CFG_get_int("sim_bw_scheme"), CFG_get_int("sim_merge_scheme"))
 
      
      print *, "The PM_initialize is done!"

   end subroutine PM_initialize

   subroutine PM_set_elec_density()
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
		PM_vars(:, PM_iv_elec) = 0.0_dp
    !print *, "Set_elec_density!"
		call PC_loop_part(add_elec_to_dens)
    ! Anbang: attention!!!: Here we set the boundary as 2 times of its density,
    ! because only half cell contributes to the grid
		PM_vars(bd_lr(1), PM_iv_elec) = 2.d0 * PM_vars(bd_lr(1), PM_iv_elec) 
		PM_vars(bd_lr(2), PM_iv_elec) = 2.d0 * PM_vars(bd_lr(2), PM_iv_elec)
    end subroutine PM_set_elec_density

    subroutine PM_set_ion_density()
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
		PM_vars(:, PM_iv_ion) = 0.0_dp
    !print *, "Set_ion_density!"
		call PC_loop_ion_part(add_ion_to_dens)
		PM_vars(bd_lr(1), PM_iv_ion) = 2.d0 * PM_vars(bd_lr(1), PM_iv_ion)
		PM_vars(bd_lr(2), PM_iv_ion) = 2.d0 * PM_vars(bd_lr(2), PM_iv_ion)
    end subroutine PM_set_ion_density

    subroutine set_elec_en_density()
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
    !print *, "Set_elec_en_density!"
		PM_vars(:, PM_iv_en) = 0.0_dp
		call PC_loop_part(add_elec_en_to_dens)
		PM_vars(bd_lr(1), PM_iv_en) = 2.d0 * PM_vars(bd_lr(1), PM_iv_en)
		PM_vars(bd_lr(2), PM_iv_en) = 2.d0 * PM_vars(bd_lr(2), PM_iv_en)
    end subroutine set_elec_en_density

    subroutine set_ion_en_density()
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
    !print *, "Set_ion_en_density!"
		PM_vars(:, PM_iv_ion_en) = 0.0_dp
		call PC_loop_ion_part(add_ion_en_to_dens)
		PM_vars(bd_lr(1), PM_iv_ion_en) = 2.d0 * PM_vars(bd_lr(1), PM_iv_ion_en)
		PM_vars(bd_lr(2), PM_iv_ion_en) = 2.d0 * PM_vars(bd_lr(2), PM_iv_ion_en)
    end subroutine set_ion_en_density

    ! set the total velocity of electrons and ions
    subroutine set_elec_vel_density()
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
		PM_vars(:, PM_iv_vel_elec) = 0.0_dp
		call PC_loop_part(add_elec_vel_to_dens)
		PM_vars(bd_lr(1), PM_iv_vel_elec) = 2.d0 * PM_vars(bd_lr(1), PM_iv_vel_elec)
		PM_vars(bd_lr(2), PM_iv_vel_elec) = 2.d0 * PM_vars(bd_lr(2), PM_iv_vel_elec)
    end subroutine set_elec_vel_density

    subroutine set_ion_vel_density()
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
		PM_vars(:, PM_iv_vel_ion) = 0.0_dp
		call PC_loop_ion_part(add_ion_vel_to_dens)
		PM_vars(bd_lr(1), PM_iv_vel_ion) = 2.d0 * PM_vars(bd_lr(1), PM_iv_vel_ion)
		PM_vars(bd_lr(2), PM_iv_vel_ion) = 2.d0 * PM_vars(bd_lr(2), PM_iv_vel_ion)
    end subroutine set_ion_vel_density

    ! set the directed energy of electrons: m/2 *<u>^2
    subroutine set_elec_directed_energy_dens()
		use m_units_constants
		use m_efield_1d
		integer  :: dbd_type,bd_lr(2)
		if (PM_useDBD) then
			call EF_get_DBD_index_and_type(bd_lr, dbd_type)
		else
			bd_lr(1)=1
			bd_lr(2)=PM_grid_size
		endif
		PM_vars(:, PM_directed_vel_elec_x) = 0.0_dp
		PM_vars(:, PM_directed_vel_elec_y) = 0.0_dp
		PM_vars(:, PM_directed_vel_elec_z) = 0.0_dp
		PM_vars(:, PM_directed_en_elec)    = 0.0_dp
		call PC_loop_part(add_elec_directedVel_to_dens)
		where (PM_vars(:, PM_iv_elec) > 0)
			PM_vars(:, PM_directed_en_elec)  = 0.5d0 * UC_elec_mass * ( PM_vars(:, PM_directed_vel_elec_x)**2 + &
				PM_vars(:, PM_directed_vel_elec_y)**2 + &
				PM_vars(:, PM_directed_vel_elec_z)**2 )
		elsewhere
			PM_vars(:, PM_directed_en_elec) = 0.0_dp
		end where
		PM_vars(bd_lr(1), PM_directed_vel_elec_x) = 2.d0 * PM_vars(bd_lr(1), PM_directed_vel_elec_x)
		PM_vars(bd_lr(2), PM_directed_vel_elec_x) = 2.d0 * PM_vars(bd_lr(2), PM_directed_vel_elec_x)
		PM_vars(bd_lr(1), PM_directed_vel_elec_y) = 2.d0 * PM_vars(bd_lr(1), PM_directed_vel_elec_y)
		PM_vars(bd_lr(2), PM_directed_vel_elec_y) = 2.d0 * PM_vars(bd_lr(2), PM_directed_vel_elec_y)
		PM_vars(bd_lr(1), PM_directed_vel_elec_z) = 2.d0 * PM_vars(bd_lr(1), PM_directed_vel_elec_z)
		PM_vars(bd_lr(2), PM_directed_vel_elec_z) = 2.d0 * PM_vars(bd_lr(2), PM_directed_vel_elec_z)
		PM_vars(bd_lr(1), PM_directed_en_elec) = 2.d0 * PM_vars(bd_lr(1), PM_directed_en_elec)
		PM_vars(bd_lr(2), PM_directed_en_elec) = 2.d0 * PM_vars(bd_lr(2), PM_directed_en_elec)
    !Anbang: the boundary is not important, so we do not set it up
    end subroutine set_elec_directed_energy_dens

   !Anbang: routine to collect the ionization numbers in one time step
   subroutine set_ionization_source_term()
        
        real(dp)   :: ioni_num(PM_grid_size)
        call PC_output_ionization_number(ioni_num)
        PM_vars(:, PM_ioni_sour_term) = ioni_num
        
   end subroutine set_ionization_source_term


   subroutine add_to_dens_cic(amount, zz, dens)
      use m_error
      real(dp), intent(in)    :: amount, zz
      real(dp), intent(inout) :: dens(:)

      integer                 :: low_ix
      real(dp)                :: temp, weight, dens_dif

      ! Index i is at (i-1) * dx
      temp   = zz * PM_inv_delta_x
      low_ix = floor(temp) + 1
      weight = low_ix - temp

      dens_dif = amount * PM_inv_grid_volume

      if (low_ix < 1 .or. low_ix > PM_grid_size-1) then
         print *, "position", zz
         call ERR_warn("Particle module: a particle is outside the computational domain")
      end if
      dens(low_ix)   = dens(low_ix) + weight * dens_dif
      dens(low_ix+1) = dens(low_ix+1) + (1-weight) * dens_dif
   end subroutine add_to_dens_cic

   subroutine add_elec_to_dens(my_part)
      type(PC_part_t), intent(in) :: my_part
      call add_to_dens_cic(real(my_part%weight, dp), my_part%x(3), PM_vars(:, PM_iv_elec))
   end subroutine add_elec_to_dens

   subroutine add_ion_to_dens(my_ion_part)
      type(PC_ion_part_t), intent(in) :: my_ion_part
      call add_to_dens_cic(real(my_ion_part%weight, dp), my_ion_part%x(3), PM_vars(:, PM_iv_ion))
   end subroutine add_ion_to_dens

   subroutine add_elec_en_to_dens(my_part)
      use m_units_constants
      type(PC_part_t), intent(in) :: my_part
      real(dp) :: energy
      energy = 0.5_dp * my_part%weight * UC_elec_mass * sum(my_part%v**2)
      call add_to_dens_cic(energy, my_part%x(3), PM_vars(:, PM_iv_en))
   end subroutine add_elec_en_to_dens

   subroutine add_ion_en_to_dens(my_ion_part)
      use m_units_constants
      type(PC_ion_part_t), intent(in) :: my_ion_part
      real(dp) :: energy
      energy = 0.5_dp * my_ion_part%weight * UC_atomic_mass * ionMassRatio * sum(my_ion_part%v**2)
      call add_to_dens_cic(energy, my_ion_part%x(3), PM_vars(:, PM_iv_ion_en))
   end subroutine add_ion_en_to_dens

  ! Anbang: calculate the velocity density in 1D
   subroutine add_elec_vel_to_dens(my_part)
      use m_units_constants
      type(PC_part_t), intent(in) :: my_part
      real(dp) :: vel
      vel = my_part%weight * my_part%v(3)
      call add_to_dens_cic(vel, my_part%x(3), PM_vars(:, PM_iv_vel_elec))
   end subroutine add_elec_vel_to_dens

   subroutine add_ion_vel_to_dens(my_part)
      use m_units_constants
      type(PC_ion_part_t), intent(in) :: my_part
      real(dp) :: vel
      vel = my_part%weight * my_part%v(3)
      call add_to_dens_cic(vel, my_part%x(3), PM_vars(:, PM_iv_vel_ion))
   end subroutine add_ion_vel_to_dens

  ! Anbang: calculate the directed energy densitys
  ! Here is from Detlef V = u + w
  ! V: total particle velocity; u: directed velocity; w :random velocity
  ! So: <m V^2/2> = m/2 *<u>^2 + m/2 *<w>^2
  ! Here we calculate the second term of the above equation
  ! we note that m/2 *<w>^2 = 3/2 * kB *Te
  ! <u>^2 = <ux>^2 +  <uy>^2  + <uz>^2                     
   subroutine add_elec_directedVel_to_dens(my_part)
      use m_units_constants
      type(PC_part_t), intent(in) :: my_part
      real(dp) :: vel(3)
      vel = my_part%weight * my_part%v
      call add_to_dens_cic(vel(1), my_part%x(3), PM_vars(:, PM_directed_vel_elec_x))
      call add_to_dens_cic(vel(2), my_part%x(3), PM_vars(:, PM_directed_vel_elec_y))
      call add_to_dens_cic(vel(3), my_part%x(3), PM_vars(:, PM_directed_vel_elec_z))
   end subroutine add_elec_directedVel_to_dens

   !> Anbang: Here need to check, if we have attachment
!    subroutine coll_callback(my_part, col_type, col_ix)
!       use m_cross_sec
!       type(PC_part_t), intent(in) :: my_part
!       integer, intent(in) :: col_type, col_ix
! 
!       select case(col_type)
!       case (CS_ionize_t)
!          call add_to_dens_cic(real(my_part%weight, dp), my_part%x(3), PM_vars(:, PM_iv_ion))
!       case (CS_attach_t)
!          call add_to_dens_cic(-real(my_part%weight, dp), my_part%x(3), PM_vars(:, PM_iv_ion))
!       end select
!    end subroutine coll_callback

   subroutine PM_update_efield(time, myrank, root)
      use m_efield_1d
      use m_units_constants

      integer, intent(in) ::  myrank, root

      real(dp),intent(in) :: time

      if ( .not. PM_useDBD) then
          call EF_compute((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge,time, [0.d0,0.d0], myrank, root)
      else
!           print *, "Time and the value of the surcharge:", time, PM_surChargeAtroot
          call EF_compute((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge,time, PM_surChargeAtroot, myrank, root)
      end if

      ! print *, myrank, "PM_update_efield is done!"
   end subroutine PM_update_efield

   subroutine PM_set_accel(my_part, accel)
      use m_efield_1d
      use m_units_constants
      type(PC_part_t), intent(in) :: my_part
      real(dp), intent(out)       :: accel(3)
      real(dp), parameter         :: accel_fac = UC_elec_charge / UC_elec_mass

      ! Only acceleration in the z-direction
      accel(1:2) = 0.0_dp
      accel(3) = accel_fac * EF_get_at_pos(my_part%x(3))
   end subroutine PM_set_accel

   subroutine PM_set_ion_accel(my_ion_part, ion_accel)
      use m_efield_1d
      use m_units_constants
      type(PC_ion_part_t), intent(in) :: my_ion_part
      real(dp), intent(out)       :: ion_accel(3)
      real(dp)                    :: accel_fac 

      accel_fac = UC_elem_charge / (UC_atomic_mass * ionMassRatio)   !> Mass to be checked

      ! Only acceleration in the z-direction
      ion_accel(1:2) = 0.0_dp
      ion_accel(3) = accel_fac * EF_get_at_pos(my_ion_part%x(3))
   end subroutine PM_set_ion_accel

   subroutine PM_advance(dt)
      real(dp), intent(in) :: dt
      call PC_advance(dt)
      
   end subroutine PM_advance

   subroutine PM_correct_accel_and_vel( dt)
        real(dp), intent(in) :: dt
        call PC_correct_new_accel(dt, PM_set_accel)
        call PC_correct_new_ion_accel(dt, PM_set_ion_accel) 
   end subroutine PM_correct_accel_and_vel

  ! get the accel of electrons and ions
   subroutine PM_get_accel()
        call PC_set_accel(PM_set_accel)
        call PC_set_ion_accel(PM_set_ion_accel)
   end subroutine PM_get_accel

    subroutine PM_adapt_weights()
        real(dp) :: v_fac

        v_fac = PM_vel_rel_weight * PM_delta_x / PM_part_per_cell
        call PC_merge_and_split((/.false., .false., .true./), v_fac, .true., &
            get_desired_weight, PC_merge_part_rxv, PC_split_part)
    end subroutine PM_adapt_weights

    real(dp) function get_desired_weight(my_part)
        type(PC_part_t), intent(in) :: my_part
        real(dp)                    :: elec_dens
        elec_dens = get_from_dens_lin_interp(my_part%x(3), PM_vars(:, PM_iv_elec))
        get_desired_weight = elec_dens * PM_grid_volume / PM_part_per_cell
        get_desired_weight = max(1.0_dp, get_desired_weight)
    end function get_desired_weight

    ! for ions
    real(dp) function get_desired_weight_for_ion(my_part)
        type(PC_ion_part_t), intent(in) :: my_part
        real(dp)                    :: ion_dens
        ion_dens = get_from_dens_lin_interp(my_part%x(3), PM_vars(:, PM_iv_ion))
        get_desired_weight_for_ion = ion_dens * PM_grid_volume / PM_part_per_cell
        get_desired_weight_for_ion = max(1.0_dp, get_desired_weight_for_ion)
    end function get_desired_weight_for_ion

    real(dp) function get_from_dens_lin_interp(zz, dens)
        real(dp), intent(in)    :: zz
        real(dp), intent(inout) :: dens(:)

        integer                 :: low_ix
        real(dp)                :: temp, weight

        ! Index i is at (i-1) * dx
        temp   = zz * PM_inv_delta_x
        low_ix = floor(temp) + 1
        weight = low_ix - temp
        get_from_dens_lin_interp = dens(low_ix) * weight + dens(low_ix+1) * (1-weight)
    end function get_from_dens_lin_interp

  !!> Anbang: here we add ion's movements
   subroutine PM_ion_advance(dt)
      real(dp), intent(in) :: dt

      call PC_advance_ions(dt)
      
   end subroutine PM_ion_advance

    subroutine PM_adapt_ion_weights()
        real(dp) :: v_fac
        
        v_fac = PM_vel_rel_weight * PM_delta_x / PM_part_per_cell
        call PC_ion_merge_and_split((/.false., .false., .true./), v_fac, .true., &
            get_desired_weight_for_ion, PC_ion_merge_part_rxv, PC_ion_split_part)
    end subroutine PM_adapt_ion_weights


   subroutine PM_get_output(pos_data, sca_data, data_names, n_pos, n_sca,time, myrank, root)
      use m_efield_1d
      use m_units_constants
      real(dp), intent(out), allocatable                :: pos_data(:,:), sca_data(:)
      character(len=name_len), intent(out), allocatable :: data_names(:)
      integer, intent(out)                              :: n_pos, n_sca
      integer                                           :: n
      real(dp)                                          :: temp_data(PM_grid_size)
      real(dp)                                          :: potential(PM_grid_size)
      real(dp), intent(in)                              :: time
      integer, intent(in)                               :: myrank, root

      n_pos = 13
      n_sca = 2
      allocate(pos_data(PM_grid_size, n_pos))
      allocate(sca_data(n_sca))
      allocate(data_names(n_pos+n_sca))

      do n = 1, PM_grid_size
         temp_data(n) = (n-1) * PM_delta_x
      end do

      data_names(1) = "position (m)"
      pos_data(:,1) = temp_data

      if ( .not. PM_useDBD) then
          call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, &
                               & potential, temp_data,time, [0.d0,0.d0], myrank, root)
      else
          call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, &
            & potential, temp_data,time, PM_surChargeAtroot, myrank, root)
      end if
      
      data_names(2) = "electric potential (V)"
      pos_data(:,2) = potential

      data_names(3) = "electric field (V/m)"
      pos_data(:,3) = temp_data

      call PM_set_elec_density()
      call PM_set_ion_density()
      call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root      
      call PM_shareDens(myrank, root)  ! broadcast to all the nodes
      data_names(4) = "electron density (1/m3)"
      pos_data(:,4) = PM_vars(:, PM_iv_elec)

      
      data_names(5) = "ion density (1/m3)"
      pos_data(:,5) = PM_vars(:, PM_iv_ion)

      call set_elec_en_density()
      call set_ion_en_density()
      call PM_collectenergyDensAtRoot(myrank, root)
      call PM_shareEnergyDens(myrank, root)
      
      data_names(6) = "electron energy density (eV/m3)"
      pos_data(:,6) = PM_vars(:, PM_iv_en) / UC_elec_volt

      data_names(7) = "electron mean energy (eV)"
      where (PM_vars(:, PM_iv_elec) > 0)
         pos_data(:,7) = PM_vars(:, PM_iv_en) / (UC_elec_volt * PM_vars(:, PM_iv_elec))
      elsewhere
         pos_data(:,7) = 0.0_dp
      end where
      
      data_names(8) = "ion energy density (eV/m3)"
      pos_data(:,8) = PM_vars(:, PM_iv_ion_en) / UC_elec_volt

      data_names(9) = "ion mean energy(eV)"
      where (PM_vars(:, PM_iv_ion) > 0)
         pos_data(:,9) = PM_vars(:, PM_iv_ion_en) / (UC_elec_volt * PM_vars(:, PM_iv_ion))
      elsewhere
         pos_data(:,9) = 0.0_dp
      end where

      call set_elec_vel_density()
      call set_ion_vel_density() 
      call PM_collectvelocityDensAtRoot(myrank, root)
      call PM_sharevelocityDens(myrank, root)

      data_names(10) = "electron power density (W/m3)"  ! Pe =(Je*E) = (ne*ve*e*E)
      pos_data(:,10) = PM_vars(:, PM_iv_vel_elec) * UC_elec_charge * pos_data(:,3)
        
      data_names(11) = "ion power density (W/m3)"  ! Pi =(Ji*E) = (ni*vi*e*E)
      pos_data(:,11) = PM_vars(:, PM_iv_vel_ion) * UC_elem_charge * pos_data(:,3)

      call set_elec_directed_energy_dens()
      call set_ionization_source_term()
      call PM_collectDiretEnDensAndSourceTermAtRoot(myrank, root)
      call PM_shareDiretEnDensAndSourceTerm(myrank, root)
      data_names(12) = "electron directed energy (eV)" 
      where (PM_vars(:, PM_iv_elec) > 0)
         pos_data(:,12) = PM_vars(:,PM_directed_en_elec) / UC_elec_volt/ PM_vars(:,PM_iv_elec)**2 ! attention: here
      elsewhere
         pos_data(:,12) = 0.0_dp
      end where
      
      data_names(13) = "the ionization source term (m^-3 * s^-1)"
      pos_data(:,13) = PM_vars(:, PM_ioni_sour_term) * PM_inv_grid_volume
   end subroutine PM_get_output


   !Anbang: Here we set up a routine to average the values during a certain time
   !        this was done by turner in his capactive discharge, we do the same to compare the results
   !        
   subroutine PM_get_aver_para(time, myrank, root)
      use m_efield_1d
      use m_units_constants

      integer                                           :: n
      real(dp)                                          :: temp_data(PM_grid_size)
      real(dp)                                          :: potential(PM_grid_size)
      real(dp),intent(in)                               :: time
      integer                                           :: idnx

      integer, intent(in)                               :: myrank, root


      do n = 1, PM_grid_size
         temp_data(n) = (n-1) * PM_delta_x
      end do
      
      !The middle grid point
      idnx = int(PM_grid_size/2.d0) + 1

      PM_pos_aver_data(:,1) =  temp_data 

      if ( .not. PM_useDBD) then
          call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, &
                               & potential, temp_data,time, [0.d0,0.d0], myrank, root)
      else
          call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, &
            & potential, temp_data,time, PM_surChargeAtroot, myrank, root)
      end if


      PM_pos_aver_data(:,2) = PM_pos_aver_data(:,2) + potential /PM_step_num_averger
      PM_pos_aver_data(:,3) = PM_pos_aver_data(:,3) + temp_data /PM_step_num_averger


      ! density
      call PM_set_elec_density()
      call PM_set_ion_density()
      call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root
      call PM_shareDens(myrank, root)  ! broadcast to all the nodes
      PM_pos_aver_data(:,4) = PM_pos_aver_data(:,4) + PM_vars(:, PM_iv_elec) / PM_step_num_averger
      PM_pos_aver_data(:,5) = PM_pos_aver_data(:,5) + PM_vars(:, PM_iv_ion) /PM_step_num_averger

      ! energy 
      call set_elec_en_density()
      call set_ion_en_density()
      call PM_collectenergyDensAtRoot(myrank, root)
      call PM_shareEnergyDens(myrank, root)
      PM_pos_aver_data(:,6) = PM_pos_aver_data(:,6) + (PM_vars(:, PM_iv_en) / UC_elec_volt) /PM_step_num_averger
      where (PM_vars(:, PM_iv_elec) > 0)
         PM_pos_aver_data(:,7) = PM_pos_aver_data(:,7) + &
                                (PM_vars(:, PM_iv_en) / (UC_elec_volt * PM_vars(:, PM_iv_elec))) / PM_step_num_averger
      elsewhere
         PM_pos_aver_data(:,7) = PM_pos_aver_data(:,7) + 0.0_dp
      end where

      PM_pos_aver_data(:,8) = PM_pos_aver_data(:,8) + &
                            (PM_vars(:, PM_iv_ion_en) / UC_elec_volt) / PM_step_num_averger

      where (PM_vars(:, PM_iv_ion) > 0)
         PM_pos_aver_data(:,9) = PM_pos_aver_data(:,9) + &
                    (PM_vars(:, PM_iv_ion_en) / (UC_elec_volt * PM_vars(:, PM_iv_ion))) /PM_step_num_averger
      elsewhere
         PM_pos_aver_data(:,9) = PM_pos_aver_data(:,9) + 0.0_dp
      end where

      ! velocity
      call set_elec_vel_density()
      call set_ion_vel_density() 
      call PM_collectvelocityDensAtRoot(myrank, root)
      call PM_sharevelocityDens(myrank, root)
      PM_pos_aver_data(:,10) = PM_pos_aver_data(:,10) + &
                              (PM_vars(:,PM_iv_vel_elec) * UC_elec_charge * temp_data(:)) / PM_step_num_averger
      PM_pos_aver_data(:,11) = PM_pos_aver_data(:,11) + &
                              (PM_vars(:,PM_iv_vel_ion) * UC_elem_charge * temp_data(:)) / PM_step_num_averger


      ! directed energy and source term
      call set_elec_directed_energy_dens()
      call set_ionization_source_term()
      call PM_collectDiretEnDensAndSourceTermAtRoot(myrank, root)
      call PM_shareDiretEnDensAndSourceTerm(myrank, root)

      where (PM_vars(:, PM_iv_elec) > 0)
            PM_pos_aver_data(:,12) = PM_pos_aver_data(:,12) + &
                (PM_vars(:,PM_directed_en_elec) / UC_elec_volt / PM_vars(:, PM_iv_elec) **2) / PM_step_num_averger
      elsewhere
            PM_pos_aver_data(:,12) = PM_pos_aver_data(:,12) + 0.0_dp
      end where
     
      PM_pos_aver_data(:,13) = PM_pos_aver_data(:,13) + &
                              PM_vars(:,PM_ioni_sour_term) * PM_inv_grid_volume / PM_step_num_averger

    ! calculate the physical parameters
    ! the ion density in the middle point
    PM_phys_aver_data(1) = PM_phys_aver_data(1) + PM_vars(idnx, PM_iv_ion) /PM_step_num_averger

    ! the electron energy in the middle point
    if (PM_vars(idnx, PM_iv_elec) > 0) then
        PM_phys_aver_data(2) = PM_phys_aver_data(2) + &
                            (PM_vars(idnx, PM_iv_en) / (UC_elec_volt * PM_vars(idnx, PM_iv_elec))) / PM_step_num_averger
    else
        PM_phys_aver_data(2) = PM_phys_aver_data(2)+ 0.0_dp
    end if

   ! line integrated electrical power coulded electrons and ions
   PM_phys_aver_data(3) =  PM_phys_aver_data(3) + sum((PM_vars(:,PM_iv_vel_elec) * UC_elec_charge * temp_data(:)) *&
                           PM_delta_x) / PM_step_num_averger
   PM_phys_aver_data(4) =  PM_phys_aver_data(4) + sum((PM_vars(:,PM_iv_vel_ion) * UC_elem_charge * temp_data(:)) * &
                           PM_delta_x) / PM_step_num_averger
   
   ! the current density collected at two electrode
   PM_phys_aver_data(5) =  PM_phys_aver_data(5) + PM_vars(1,PM_iv_vel_ion) * UC_elem_charge / PM_step_num_averger
   PM_phys_aver_data(6) =  PM_phys_aver_data(6) + PM_vars(PM_grid_size,PM_iv_vel_ion) * UC_elem_charge / PM_step_num_averger

   end subroutine PM_get_aver_para


   !Anbang: Here we output the time evaluation of parameter during one time step after peroidic steady state is accieved.
   subroutine PM_aver_para_one_peroid(num, num_ave, time, myrank, root)
      use m_efield_1d
      use m_units_constants

      integer                                           :: n
      real(dp)                                          :: temp_data(PM_grid_size)
      real(dp)                                          :: potential(PM_grid_size)
      real(dp),intent(in)                               :: time
      integer                                           :: idnx

      integer, intent(in)                               :: myrank, root
      integer, intent(in)                               :: num, num_ave  ! the latter is how many times are averaged


      do n = 1, PM_grid_size
         temp_data(n) = (n-1) * PM_delta_x
      end do
      
      !The middle grid point
      idnx = int(PM_grid_size/2.d0) + 1

      PM_data_per_peroid(num,:,1) =  temp_data 

      ! density
      call PM_set_elec_density()
      call PM_set_ion_density()
      call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root
      call PM_shareDens(myrank, root)  ! broadcast to all the nodes

      if ( .not. PM_useDBD) then
          call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, &
                               & potential, temp_data,time, [0.d0,0.d0], myrank, root)
      else
          call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, &
            & potential, temp_data,time, PM_surChargeAtroot, myrank, root)
      end if

      PM_data_per_peroid(num,:,2) = PM_data_per_peroid(num,:,2) + potential /num_ave
      PM_data_per_peroid(num,:,3) = PM_data_per_peroid(num,:,3) + temp_data /num_ave


      PM_data_per_peroid(num,:,4) = PM_data_per_peroid(num,:,4) + PM_vars(:, PM_iv_elec) / num_ave
      PM_data_per_peroid(num,:,5) = PM_data_per_peroid(num,:,5) + PM_vars(:, PM_iv_ion) /num_ave

      ! energy 
      call set_elec_en_density()
      call set_ion_en_density()
      call PM_collectenergyDensAtRoot(myrank, root)
      call PM_shareEnergyDens(myrank, root)
      PM_data_per_peroid(num,:,6) = PM_data_per_peroid(num,:,6) + (PM_vars(:, PM_iv_en) / UC_elec_volt) /num_ave
      where (PM_vars(:, PM_iv_elec) > 0)
         PM_data_per_peroid(num,:,7) = PM_data_per_peroid(num,:,7) + &
                                (PM_vars(:, PM_iv_en) / (UC_elec_volt * PM_vars(:, PM_iv_elec))) / num_ave
      elsewhere
         PM_data_per_peroid(num,:,7) = PM_data_per_peroid(num,:,7) + 0.0_dp
      end where

      PM_data_per_peroid(num,:,8) = PM_data_per_peroid(num,:,8) + &
                            (PM_vars(:, PM_iv_ion_en) / UC_elec_volt) / num_ave

      where (PM_vars(:, PM_iv_ion) > 0)
         PM_data_per_peroid(num,:,9) = PM_data_per_peroid(num,:,9) + &
                    (PM_vars(:, PM_iv_ion_en) / (UC_elec_volt * PM_vars(:, PM_iv_ion))) /num_ave
      elsewhere
         PM_data_per_peroid(num,:,9) = PM_data_per_peroid(num,:,9) + 0.0_dp
      end where

      ! velocity
      call set_elec_vel_density()
      call set_ion_vel_density() 
      call PM_collectvelocityDensAtRoot(myrank, root)
      call PM_sharevelocityDens(myrank, root)
      PM_data_per_peroid(num,:,10) = PM_data_per_peroid(num,:,10) + &
                              (PM_vars(:,PM_iv_vel_elec) * UC_elec_charge * temp_data(:)) / num_ave
      PM_data_per_peroid(num,:,11) = PM_data_per_peroid(num,:,11) + &
                              (PM_vars(:,PM_iv_vel_ion) * UC_elem_charge * temp_data(:)) / num_ave


      ! directed energy and source term
      call set_elec_directed_energy_dens()
      call set_ionization_source_term()
      call PM_collectDiretEnDensAndSourceTermAtRoot(myrank, root)
      call PM_shareDiretEnDensAndSourceTerm(myrank, root)

      where (PM_vars(:, PM_iv_elec) > 0)
            PM_data_per_peroid(num,:,12) = PM_data_per_peroid(num,:,12) + &
                                    (PM_vars(:,PM_directed_en_elec) / UC_elec_volt / PM_vars(:, PM_iv_elec) **2) / num_ave
      elsewhere
            PM_data_per_peroid(num,:,12) = PM_data_per_peroid(num,:,12) + 0.0_dp
      end where
     
      PM_data_per_peroid(num,:,13) = PM_data_per_peroid(num,:,13) + &
                              PM_vars(:,PM_ioni_sour_term) * PM_inv_grid_volume / num_ave

      ! the current density
      PM_data_per_peroid(num,:,14) = PM_data_per_peroid(num,:,14) + &
                              (PM_vars(:,PM_iv_vel_elec) * UC_elec_charge ) / num_ave
      PM_data_per_peroid(num,:,15) = PM_data_per_peroid(num,:,15) + &
                              (PM_vars(:,PM_iv_vel_ion) * UC_elem_charge ) / num_ave

   end subroutine PM_aver_para_one_peroid

  !Anbang: output the average parameters during one peroid after peroidic steady state
    subroutine PM_output_aver_para_one_peroid(num, sim_name, time, peroid)

        character(len=*), intent(in) :: sim_name
        real(dp), intent(in) :: time, peroid
        character(len=name_len) :: filename
        integer              :: myunit, nn, mm
        integer, intent(in)  :: num
        real(dp)             :: timeOut

        do mm = 1, num
            myunit = 15795 + mm
            write(filename, fmt="(I0,A)"), mm, ".txt"
            filename = "output/" // trim(sim_name) // "_paOnePd_"// trim(filename)
            open(unit = myunit, file = filename , status = 'unknown')
            timeOut = time + (mm - 1) * peroid / num
            do nn = 1, PM_grid_size
                write(myunit, *) timeOut, PM_data_per_peroid(mm,nn,:)
            end do
            close(myunit)
        end do

    end subroutine PM_output_aver_para_one_peroid

   !!Anbang: here i calculate the mean particle number and mean weight per cell
    subroutine PM_get_aver_num_and_weight_per_cell()

        integer :: i
        integer :: tot_elec_sim_num_per_cell(PM_grid_size - 1)
        real(dp) :: tot_elec_mean_weight_per_cell(PM_grid_size - 1)
        integer :: tot_ion_sim_num_per_cell(PM_grid_size - 1)
        real(dp) :: tot_ion_mean_weight_per_cell(PM_grid_size - 1)
        
        tot_elec_mean_weight_per_cell = 0.d0
        tot_elec_sim_num_per_cell = 0

        call PC_cal_elec_num_and_weight_per_cell(tot_elec_sim_num_per_cell, tot_elec_mean_weight_per_cell)
        call PC_cal_ion_num_and_weight_per_cell(tot_ion_sim_num_per_cell, tot_ion_mean_weight_per_cell)

        do i =  1, PM_grid_size -1    
            PM_part_num_per_cell(1,i) = PM_part_num_per_cell(1,i) + (tot_elec_sim_num_per_cell(i) / PM_step_num_averger)
            PM_part_weight_per_cell(1,i) = PM_part_weight_per_cell(1,i) + &
                        & tot_elec_mean_weight_per_cell(i) / PM_step_num_averger
            PM_part_num_per_cell(2,i) = PM_part_num_per_cell(2,i) + (tot_ion_sim_num_per_cell(i) / PM_step_num_averger)
            PM_part_weight_per_cell(2,i) = PM_part_weight_per_cell(2,i) + &
                        & tot_ion_mean_weight_per_cell(i) / PM_step_num_averger                 
        end do
        
    end subroutine PM_get_aver_num_and_weight_per_cell

  !Anbang: output the average values caculated above  for particle number and weight per cell
    subroutine PM_output_aver_num_and_weight_per_cell(sim_name)

        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer              :: myunit, nn
        
        myunit = 521115
        filename = "output/" // trim(sim_name) // "_part_" // "num_weight.txt"
        open(unit = myunit, file = filename , status = 'unknown')
        write(myunit, *) "pos/ elec_num/ ion_num/ elec_mean_weight/ ion_mean_weight"
        do nn = 1, PM_grid_size - 1
            write(myunit, *) (nn-0.5d0) * PM_delta_x, PM_part_num_per_cell(:,nn), PM_part_weight_per_cell(:,nn) 
        end do
        close(myunit)

        print *, "output the averged particle number and mean weight per cell: PM_output_aver_num_and_weight_per_cell!"
    end subroutine PM_output_aver_num_and_weight_per_cell

 !Anbang: calculate the sheath position during one peroid
    subroutine PM_sheath_pos_cal_output(num, sim_name, time, peroid)

        use m_units_constants

        character(len=*), intent(in) :: sim_name
        real(dp), intent(in) :: time, peroid
        character(len=name_len) :: filename
        integer              :: myunit, nn, mm
        integer, intent(in)  :: num
        real(dp),allocatable :: timeOut(:)

       allocate(timeOut(num))
        
        ! find the sheath pos: two methods: 1). ne/ni = 1/2;  2).(ni-ne)/ni =0.01
        ! Anbang: this method is not very efficient, but we only calculated once that is ok
        do mm = 1, num
            timeOut(mm) = time + (mm - 1) * peroid / num

            !Anbang: left sheath, method1
            do nn = 1, int(PM_grid_size/2.d0)
                if ( abs(PM_data_per_peroid(mm,nn,4)/(PM_data_per_peroid(mm,nn,5)+ UC_tiny)) >=0.50) then
                    PM_pos_sheath1(mm,1) = nn * PM_delta_x
                    exit
                end if
            end do

            !Anbang: left sheath, method 2
            do nn = 1, int(PM_grid_size/2.d0)
                if ( abs ((PM_data_per_peroid(mm,nn,5) - PM_data_per_peroid(mm,nn,4))/ &
                        &(PM_data_per_peroid(mm,nn,5)+ UC_tiny)) <=1.d-2 ) then
                    PM_pos_sheath2(mm,1) = nn * PM_delta_x
                    exit
                end if
            end do

            !Anbang: right sheath, method 1
            do nn = PM_grid_size, int(PM_grid_size/2.d0),-1
                if ( abs(PM_data_per_peroid(mm,nn,4)/(PM_data_per_peroid(mm,nn,5)+ UC_tiny)) >=0.50) then
                    PM_pos_sheath1(mm,2) = nn * PM_delta_x
                    exit
                end if
            end do

            !Anbang: right sheath, method 2
            do nn = PM_grid_size, int(PM_grid_size/2.d0), -1
                if ( abs ((PM_data_per_peroid(mm,nn,5) - PM_data_per_peroid(mm,nn,4))/ &
                     &(PM_data_per_peroid(mm,nn,5)+ UC_tiny)) <=1.d-2 ) then
                   PM_pos_sheath2(mm,2) = nn * PM_delta_x
                   exit
                end if
            end do

        end do

        myunit = 79621
        filename = "output/" // trim(sim_name) // "_sheath_pos"// ".txt"
        open(unit = myunit, file = filename , status = 'unknown')
        do mm = 1, num
            write(myunit, *) timeOut(mm), PM_pos_sheath1(mm,:), PM_pos_sheath2(mm,:)
        end do
        close(myunit)

    end subroutine PM_sheath_pos_cal_output

  !Anbang: output the average values caculated above 
    subroutine PM_output_aver_para(sim_name, time)

        character(len=*), intent(in) :: sim_name
        real(dp), intent(in) :: time
        character(len=name_len) :: filename
        integer              :: myunit, nn
        
        myunit = 16789
        filename = "output/" // trim(sim_name) // "_part_" // "aver_para.txt"
        open(unit = myunit, file = filename , status = 'unknown')
        write(myunit, *) "Here we output average value during", PM_step_num_averger, "steps"
!         write(myunit, *) "position (m) / electric potential (V)/ electric field (V/m)/ electron density (1/m3)" &
!                        & "/ ion density (1/m3)/ electron energy density (eV) / electron mean energy densiity (eV/m3) " &
!                        & "/ ion energy density (eV/m3)/ ion mean energy density (eV)"
        do nn = 1, PM_grid_size
            write(myunit, *) PM_pos_aver_data(nn,:)
        end do
        close(myunit)

        print *, "output the average values " // trim(filename) // " at t = ", time
    end subroutine PM_output_aver_para

  !Anbang: output the average physical values caculated above 
    subroutine PM_output_aver_phys_para(sim_name, time)

        character(len=*), intent(in) :: sim_name
        real(dp), intent(in) :: time
        character(len=name_len) :: filename
        integer              :: myunit
        
        myunit = 1639273
        filename = "output/" // trim(sim_name) // "_part_" // "averphys_para.txt"
        open(unit = myunit, file = filename , status = 'unknown')
        write(myunit, *) "Here we output average physical values during", PM_step_num_averger, "steps"
        write(myunit, *) "And the time is:", time
        write(myunit, *) "ion densiy & electron energy in the middle plasne/Se(W m-2)/ Si(W m-2)/ Ji(A m-2) at two electrodes"

        write(myunit, *) PM_phys_aver_data(:)
        write(myunit, *) "write the averaged particle number and weight of all particles"
        write(myunit, *) sum(PM_part_num_per_cell(1,:)) / (PM_grid_size - 1)
        write(myunit, *) sum(PM_part_num_per_cell(2,:)) / (PM_grid_size - 1)
        write(myunit, *) sum(PM_part_weight_per_cell(1,:)) / (PM_grid_size - 1)
        write(myunit, *) sum(PM_part_weight_per_cell(2,:)) / (PM_grid_size - 1)

        close(myunit)

    end subroutine PM_output_aver_phys_para

   subroutine PM_get_coeffs(coeff_data, coeff_names, n_coeffs)
      use m_cross_sec
      use m_lookup_table
      real(dp), intent(out), allocatable :: coeff_data(:,:)
      character(len=name_len), intent(out), allocatable :: coeff_names(:)
      integer, intent(out) :: n_coeffs
      type(PC_coll_t) :: coll_data
      integer :: nn, n_rows

      call PC_get_colls(coll_data)
      n_coeffs = coll_data%num + 2
      n_rows = LT_get_num_rows(coll_data%rate_lt)
      allocate(coeff_data(n_coeffs, n_rows))
      allocate(coeff_names(n_coeffs))

      call LT_get_all_data(coll_data%sum_rate_lt, coeff_data(1, :), coeff_data(2:2,:))
      call LT_get_all_data(coll_data%rate_lt, coeff_data(1, :), coeff_data(3:,:))
      coeff_names(1) = "velocity (m/s)"
      coeff_names(2) = "sum coll_rate (1/s)"
      do nn = 1, coll_data%num
         select case (coll_data%types(nn))
         case (CS_ionize_t)
            coeff_names(2+nn) = "ionization"
         case (CS_attach_t)
            coeff_names(2+nn) = "attachment"
         case (CS_elastic_t)
            coeff_names(2+nn) = "elastic"
         case (CS_excite_t)
            coeff_names(2+nn) = "excitation"
         case default
            coeff_names(2+nn) = "unknown"
         end select
      end do
   end subroutine PM_get_coeffs


  !> New added
    !> for electrons, transfer all the particles from all the nodes
    integer function PM_get_num_sim_part()
        use m_error
        PM_get_num_sim_part = PC_getNumSimelecPartMPI()
        if (PM_get_num_sim_part==0) then
            call ERR_warn("Error: There is no electrons left in the domain.")
            stop
        end if
    end function PM_get_num_sim_part

    real(dp) function PM_get_num_real_part()
        PM_get_num_real_part = PC_getNumRealelecPartMPI()
    end function PM_get_num_real_part
    
     !> For ions
    integer function PM_get_num_sim_ion_part()
        PM_get_num_sim_ion_part = PC_getNumSimionPartMPI()
    end function PM_get_num_sim_ion_part

    real(dp) function PM_get_num_real_ion_part()
        PM_get_num_real_ion_part = PC_getNumRealionPartMPI()
    end function PM_get_num_real_ion_part

    ! For the EEDF/EEPF of electrons column: 1- eV, 2- EEDF, 3-EEPF
    subroutine PM_get_eedf(eedf, energy_range, field_range, n_bins)
        use m_units_constants
        real(dp), intent(inout), allocatable :: eedf(:,:)
        real(dp), intent(in) :: energy_range(2), field_range(2)
        integer, intent(in) :: n_bins

        integer :: ix
        real(dp) :: accel_range(2)
        
        allocate(eedf(3, n_bins))

        do ix = 1, n_bins
        eedf(1, ix) = energy_range(1) + (ix-0.5_dp) * &
                (energy_range(2) - energy_range(1)) / n_bins
        end do

        ! Convert to range of velocity**2
        eedf(1,:) = eedf(1,:) * UC_elec_volt / (0.5_dp * UC_elec_mass)
        accel_range = field_range * UC_elem_charge / UC_elec_mass

        call PC_histogram(get_vel_squared, select_part_by_accel, &
            accel_range, eedf(1,:), eedf(2,:),eedf(3,:))

        ! Convert to energy range
        eedf(1,:) = eedf(1,:) * 0.5_dp * UC_elec_mass / UC_elec_volt

    end subroutine PM_get_eedf

    real(dp) function get_vel_squared(my_part)
        type(PC_part_t), intent(in) :: my_part
        get_vel_squared = sum(my_part%v**2)
    end function get_vel_squared


    logical function select_part_by_accel(my_part, accel_range)
        type(PC_part_t), intent(in) :: my_part
        real(dp), intent(in)        :: accel_range(:)
        real(dp)                    :: accel
        accel                = norm2(my_part%a)
        select_part_by_accel = accel >= accel_range(1) .and. accel <= accel_range(2)
    end function select_part_by_accel


    !> get the charge density on the surface of dielectrics
    !> Anbang: here i need to be checked with markus
    subroutine PM_get_surface_charge_density(dt)
        use m_units_constants
        use m_efield_1d
        
        real(dp),intent(in)     :: dt
        real(dp)                :: surfElecFlux(2)
        real(dp)                :: surfIonFlux(2)
        integer                 :: dbd_index(2),dbd_type
	PM_surChargeAtnodes=0.d0

	call PC_curr_cal(PM_surChargeAtnodes)
	PM_surChargeAtnodes = UC_elem_charge*PM_surChargeAtnodes/PM_transverse_area
        call EF_get_DBD_index_and_type(dbd_index, dbd_type)

    select case (dbd_type)
    case(1)
        PM_surChargeAtnodes(2) = 0.d0
    case(2)
        PM_surChargeAtnodes(1) = 0.d0
    end select
!         print *, "PM_surChargeAtnodes =", PM_surChargeAtnodes
        
    end subroutine PM_get_surface_charge_density


   ! Anbang: here we save PM_surChargeAtroot during intermediate steps
    subroutine PM_save_surf_charge_info(myrank, root, sim_name)
        integer, intent(IN) :: myrank, root
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer                 :: myunit, i

        if (myrank == root) then 

            myunit = 9856
            filename = "midstore/" // trim(sim_name) // "_info_srf"// ".txt"
            open(unit = myunit, file = filename , status = 'unknown')
            write(myunit, *) PM_surChargeAtroot
            close(myunit)
         end if

    end subroutine PM_save_surf_charge_info

   ! Anbang: here we read PM_surChargeAtroot during intermediate steps
    subroutine PM_read_surf_charge_info(myrank, root, sim_name)
        integer, intent(IN) :: myrank, root
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer                 :: myunit, i

        if (myrank == root) then 

            myunit = 9856
            filename = "midstore/" // trim(sim_name) // "_info_srf"// ".txt"
            open(unit = myunit, file = filename , status = 'old')
            write(myunit, *) PM_surChargeAtroot
            close(myunit)
         end if

    end subroutine PM_read_surf_charge_info


!###########################################
!!!!Routines for MPI########################
!!##########################################
   !> Collect the surface charge at root
   !! Important:
   subroutine PM_collectsurfChargeAtRoot(myrank, root)

      include "mpif.h"

      integer, intent(in)           :: myrank, root

      integer                       :: ierr, dataSize
    
      dataSize = 2
      
      if (myrank == root) then   ! collect at root
            call MPI_REDUCE(MPI_IN_PLACE, PM_surChargeAtnodes, dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            PM_surChargeAtroot = PM_surChargeAtroot + PM_surChargeAtnodes   ! the surf charge is accumulated all the time
!             print *, "the surface charge at the root is:", PM_surChargeAtroot

      else 
             call MPI_REDUCE(PM_surChargeAtnodes, PM_surChargeAtnodes, dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

             PM_surChargeAtnodes = 0.d0

      end if

   end subroutine PM_collectsurfChargeAtRoot

   !> Collect electron and ion densities.
   !! Important: for the particle code, we collect all the electron/ion density to root
   subroutine PM_collectDensitiesAtRoot(myrank, root)

      include "mpif.h"

      integer, intent(in)           :: myrank, root

      integer                       :: ierr, dataSize


      dataSize = PM_grid_size

      if (myrank == root) then   ! collect at root
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(:,PM_iv_elec), dataSize, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(:, PM_iv_ion), dataSize, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      else
            call MPI_REDUCE(PM_vars(:,PM_iv_elec), PM_vars(:,PM_iv_elec), dataSize, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(PM_vars(:, PM_iv_ion), PM_vars(:, PM_iv_ion), dataSize, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr) 

            PM_vars(:, PM_iv_ion) = 0.d0
            PM_vars(:, PM_iv_elec) = 0.d0
      end if
   end subroutine PM_collectDensitiesAtRoot

  ! broadcast to all the nodes, be careful here all the nodes have the right densities
   subroutine PM_shareDens(myrank, root)

      include "mpif.h"

      integer, intent(in)  :: myrank, root
      integer              :: dataSize, ierr
      real(dp), allocatable :: temp1(:), temp2(:)

      dataSize = PM_grid_size
      allocate(temp1(dataSize))
      allocate(temp2(dataSize))

      if (myrank == root) then
          temp1 = PM_vars(:,PM_iv_elec)
          temp2 = PM_vars(:,PM_iv_ion)
      end if

      call MPI_BCAST(temp1, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(temp2, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)

      if (myrank /= root) then
            PM_vars(:,PM_iv_elec) = temp1
            PM_vars(:,PM_iv_ion)  = temp2
      end if

      deallocate(temp1,temp2)
 
   end subroutine PM_shareDens
   

   !> Collect the energy densities of electron and ion 
   !! Important: for the particle code, we collect all the electron/ion density to root
   subroutine PM_collectenergyDensAtRoot(myrank, root)

      include "mpif.h"
      integer, intent(in)           :: myrank, root
      integer                       :: ierr, dataSize
    
      dataSize = PM_grid_size
      
      if (myrank == root) then   ! collect at root
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(: ,PM_iv_en), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(: ,PM_iv_ion_en), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
       else
            call MPI_REDUCE(PM_vars(: ,PM_iv_en), PM_vars(: ,PM_iv_en), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(PM_vars(: ,PM_iv_ion_en), PM_vars(: ,PM_iv_ion_en), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr) 

            PM_vars(: ,PM_iv_en) = 0.d0
            PM_vars(: ,PM_iv_ion_en) = 0.d0

      end if

   end subroutine PM_collectenergyDensAtRoot

  ! broadcast to all the nodes, be careful here all the nodes have the right densities
   subroutine PM_shareEnergyDens(myrank, root)

      include "mpif.h"
      integer, intent(in)  :: myrank, root
      integer              :: dataSize, ierr
      real(dp), allocatable :: temp1(:), temp2(:)
    
      dataSize = PM_grid_size

      allocate(temp1(dataSize))
      allocate(temp2(dataSize))

      if (myrank == root) then
          temp1 = PM_vars(:,PM_iv_en)
          temp2 = PM_vars(:,PM_iv_ion_en)
      end if
      
      call MPI_BCAST(temp1, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(temp2, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)

      if (myrank /= root) then
          PM_vars(:,PM_iv_en) = temp1
          PM_vars(:,PM_iv_ion_en) = temp2
      end if
    
       deallocate(temp1,temp2)
   end subroutine PM_shareEnergyDens    


   !> Collect the velocity densities of electron and ion 
   subroutine PM_collectvelocityDensAtRoot(myrank, root)

      include "mpif.h"
      integer, intent(in)           :: myrank, root
      integer                       :: ierr, dataSize
    
      dataSize = PM_grid_size
      
      if (myrank == root) then   ! collect at root
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(: ,PM_iv_vel_elec), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(: ,PM_iv_vel_ion), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
       else
            call MPI_REDUCE(PM_vars(: ,PM_iv_vel_elec), PM_vars(: ,PM_iv_vel_elec), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(PM_vars(: ,PM_iv_vel_ion), PM_vars(: ,PM_iv_vel_ion), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

            PM_vars(: ,PM_iv_vel_elec) = 0.d0
            PM_vars(: ,PM_iv_vel_ion) =0.d0
      end if

   end subroutine PM_collectvelocityDensAtRoot

  ! broadcast to all the nodes, be careful here all the nodes have the right densities
   subroutine PM_sharevelocityDens(myrank, root)

      include "mpif.h"
      integer, intent(in)  :: myrank, root
      integer              :: dataSize, ierr
      real(dp), allocatable:: temp1(:), temp2(:)

      dataSize = PM_grid_size
      allocate(temp1(dataSize))
      allocate(temp2(dataSize))

      if (myrank == root) then
          temp1 = PM_vars(:,PM_iv_vel_elec)
          temp2 = PM_vars(:,PM_iv_vel_ion)
      end if

      call MPI_BCAST(temp1, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(temp2, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)

      if (myrank /= root) then
         PM_vars(:,PM_iv_vel_elec) = temp1
         PM_vars(:,PM_iv_vel_ion) = temp2
      end if
      deallocate(temp1, temp2)

   end subroutine PM_sharevelocityDens    

   !> Collect the directed energy density and the source term
   subroutine PM_collectDiretEnDensAndSourceTermAtRoot(myrank, root)

      include "mpif.h"
      integer, intent(in)           :: myrank, root
      integer                       :: ierr, dataSize
    
      dataSize = PM_grid_size
      
      if (myrank == root) then   ! collect at root
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(: ,PM_directed_en_elec), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, PM_vars(: ,PM_ioni_sour_term), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      else
            call MPI_REDUCE(PM_vars(: ,PM_directed_en_elec), PM_vars(: ,PM_directed_en_elec), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(PM_vars(: ,PM_ioni_sour_term), PM_vars(: ,PM_ioni_sour_term), dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

            PM_vars(: ,PM_directed_en_elec) = 0.d0
            PM_vars(: ,PM_ioni_sour_term) = 0.d0

      end if

   end subroutine PM_collectDiretEnDensAndSourceTermAtRoot

  ! broadcast to all the nodes, be careful here all the nodes have the right densities
   subroutine PM_shareDiretEnDensAndSourceTerm(myrank, root)

      include "mpif.h"
      integer, intent(in)  :: myrank, root
      integer              :: dataSize, ierr
      real(dp),allocatable :: temp1(:), temp2(:)

      dataSize = PM_grid_size
      allocate(temp1(dataSize))
      allocate(temp2(dataSize))

      if (myrank == root) then
          temp1 = PM_vars(:,PM_directed_en_elec)
          temp2 = PM_vars(:,PM_ioni_sour_term)
      end if

      call MPI_BCAST(temp1, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(temp2, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)

      if (myrank /= root) then
          PM_vars(:,PM_directed_en_elec) = temp1
          PM_vars(:,PM_ioni_sour_term) = temp2
      end if

      deallocate(temp1, temp2)
   end subroutine PM_shareDiretEnDensAndSourceTerm 


   !> Collect the eedf at the root
    subroutine PM_collectEEDFAtRoot(eedf, myrank, root, ntasks, n_bins)
        include "mpif.h"
        integer, intent(in) :: myrank, root, ntasks, n_bins
        integer             :: dataSize,ierr
      !  real(dp), allocatable     :: temp1(:,:) 

        real(dp), intent(inout)  :: eedf(3,n_bins)

       ! allocate(temp1(3,n_bins))
        dataSize = 3 * n_bins

        if (myrank == root) then  ! the root collect the eedf
            call MPI_REDUCE(MPI_IN_PLACE, eedf, dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
 !           eedf = temp1
            eedf = eedf / ntasks

        else
            call MPI_REDUCE(eedf, eedf, dataSize, &
                MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

            eedf = 0.d0
    
        end if
    end subroutine PM_collectEEDFAtRoot

    subroutine PM_shareEEDF(eedf, myrank, root, n_bins)
        include "mpif.h"
        integer, intent(in) :: myrank, root, n_bins
        integer             :: dataSize, ierr

        real(dp), intent(inout)  :: eedf(3,n_bins)
        real(dp), allocatable :: temp(:,:)

        dataSize = 3 * n_bins
        allocate(temp(3,n_bins))
        if (myrank == root) temp = eedf

        call MPI_BCAST(temp, dataSize, MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ierr)
        if (myrank /= root) eedf = temp
        deallocate(temp)

    end subroutine PM_shareEEDF

!!! anbang: new added staff, for current
  ! anbang: calculate the conducting,displacement and total current
   subroutine PM_calculate_curr(dis_curr, conduct_curr, tot_curr,tot_curr_density, dt)

        use m_units_constants
        use m_efield_1d

        real(dp), intent(out),allocatable :: dis_curr(:)
        real(dp), intent(out), allocatable :: conduct_curr(:,:), tot_curr(:),tot_curr_density(:)
        real(dp), intent(in)  :: dt
        integer               :: indx(4)
        integer               :: j
        
        call EF_index_for_current(indx)

        allocate(dis_curr(size(indx)))
        allocate(conduct_curr(2, size(indx)))
        allocate(tot_curr(size(indx)))
		allocate(tot_curr_density(size(indx)))
        
        call EF_displace_currtent_at(dt, dis_curr, PM_transverse_area)

        conduct_curr = 0.d0
        tot_curr = 0.d0
		tot_curr_density=0.d0

        do j = 1, size(indx)
            ! for electron
            conduct_curr(1,j) = PM_vars(indx(j), PM_iv_vel_elec) * UC_elec_charge * PM_transverse_area
            ! for ions
            conduct_curr(2,j) = PM_vars(indx(j), PM_iv_vel_ion)  * &
                    UC_elem_charge * PM_transverse_area

            ! total current
            tot_curr(j) = conduct_curr(1,j) + conduct_curr(2,j) + dis_curr(j)
			tot_curr_density(j) = tot_curr(j) / PM_transverse_area
        end do           
   end subroutine PM_calculate_curr

    ! get the maximum density of electrons and ions
    subroutine PM_get_max_dens(max_dens,sur_all)
        real(dp), intent(out) :: max_dens(2)
    real(dp), intent(out) :: sur_all(2)
        
        max_dens(1) = maxval(PM_vars(:,PM_iv_elec))
        max_dens(2) = maxval(PM_vars(:,PM_iv_ion))
    sur_all=PM_surChargeAtroot
    end subroutine PM_get_max_dens

end module m_particle_1d
