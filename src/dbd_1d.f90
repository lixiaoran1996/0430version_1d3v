!> \mainpage
!! This webpage contains documentation for a 1d streamer simulation simulation code.
!! It can simulate these 'ionization waves' or 1d streamer using either a particle or various fluid models.
!!
!! \section secLabel Usage overview
!! Usage with the default configuration is easy: just type "make" to build the program and then execute
!! "dbd_1d" with no arguments. If you want to execute the program with a configuration file,
!! then use "dbd_1d [file]".
!! Input files are located in the "input" directory. Output files are located in the "output" directory.
!! The input file for the particle model contains cross sections.
!! The input file for the fluid models contains transport data (mobility, diffusion, ionization rate etc.).
!! The output files have the name and type of the simulation prepended (e.g., simName_type_000001.txt for the first one).
!! They contain the following columns: position, electric field, electron density, ion density, mean energy and possibly more.


program dbd_1d
   use m_types
   use m_model_choice
   use m_config
   use m_efield_1d
   use m_output_1d
   use m_fluid_dd_1d
   use m_particle_1d
   use m_particle_core
   use m_gas
   use m_error
   use m_init_cond_1d
   use m_transport_data
   use m_cross_sec
   use m_cross_sec_for_ions

   include "mpif.h"

  ! implicit none

   character(len=80) :: sim_name, fluid_input_file, cfg_name, gas_mix_name, outTimeFile

   integer,parameter :: my_unit = 1256, my_unit_flux = 2367 
   integer           :: sim_type
   integer           :: steps, n_steps_apm
   integer           :: output_cntr
   integer           :: info_cntr
   integer           :: prev_apm_elec_part, prev_apm_ion_part
   integer           :: n_elec_sim_sum, n_ion_sim_sum
   real(dp)          :: n_elec_real_sum, n_ion_real_sum,sur_all(6)
   real(dp), allocatable          :: flux_all(:)
   integer           :: timeStepIonOverElec

   integer           :: stepIons

   real(dp)          :: err_goal
   real(dp)          :: sim_time, sim_end_time
   real(dp)          :: dt, max_dt, new_dt
   real(dp)          :: output_dt
   real(dp)          :: time_now, time_start
   real(dp)          :: apm_increase
   
   real(dp)          ::potential_left_bound,potential_gap, potential_dielectric

   logical            :: do_apm, do_apm_ion, output_numerical_para
   logical            :: use_DBD

   ! parameters in turner's code, to average the values in a certain time
   logical           :: doAverage
   integer           :: stepToExcute
   integer           :: stepToAverage

   ! parameters for calculating time evaluation of paramerters during one peroid 
   real(dp)          :: potPreStep
   real(dp)          :: timeStartEva
   logical           :: flagfindTime, outputTimeEvaPra
   integer           :: ave_num_timeEva
   integer           :: sim_voltage_type
   integer           :: div_num_timeEva
   real(dp)          :: peroid 
   real(dp)          :: outputMoment
   integer           :: cnrt

    ! parameters for current output
    real(dp), allocatable  :: dis_curr(:), conduct_curr(:,:), tot_curr(:),tot_curr_density(:)
    integer     :: time_para_output_interval 
   ! for output electric field
    real(dp), allocatable :: efield_ix(:)
    
    !output max electron and ion densities
    real(dp)         :: max_dens(2)

    integer :: nn, n_gas_comp
    character(len=name_len) :: sim_type_name, cs_file, cs_file_for_ions
    character(len=name_len), allocatable :: gas_comp_names(:)
    real(dp), allocatable :: gas_comp_fracs(:)
    type(CS_type), allocatable :: cross_secs(:)
    type(CS_ion_type), allocatable :: cross_secs_for_ions(:)


   ! parameters for MPI
   integer  :: ierr, myrank, ntasks, root = 0

   integer  :: tempElecNum, tempIonNum
   logical  :: didAPMforElec, didAPMforIons
   logical  :: do_verlet_bf_coll

    ! save and read particle info
    logical   :: save_mid_info, read_mid_info
    integer   :: save_step_interval

    ! for convegence test
    real(dp)  :: factor_convegence
    character(len=80) :: read_file_name

    ! Set up MPI environment
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    

  ! call initialize_all() ! Initialize all necessary modules (routine contained in program)

      call create_sim_config()      ! Create default parameters for the simulation (routine contained below)

      call CFG_sort()               ! Sort parameters in config module

      ! Use either no configuration file, or the one specified as argument
      select case (command_argument_count())
      case (0)
         continue
      case (1)
         call get_command_argument(1, cfg_name)

         call CFG_read_file(trim(cfg_name))
         
      case default
         call ERR_show("dbd_1d expects at most 1 argument: the configuration filename")
      end select

      call CFG_get("sim_name", sim_name)
      call CFG_get("sim_type", sim_type_name)

      if (myrank == root) call CFG_write("output/" // trim(sim_name) // "_config.txt")

      call MODEL_initialize(sim_type_name)
      sim_type = MODEL_get_type()    ! Get the type of simulation model

      n_gas_comp = CFG_get_size("gas_component_names")
      if (n_gas_comp /= CFG_get_size("gas_component_fractions")) &
           call ERR_show("dbd_1d: variables gas_component_names and gas_component_fracs have unequal size")
      allocate(gas_comp_names(n_gas_comp))
      allocate(gas_comp_fracs(n_gas_comp))
      if (CFG_get_size("gas_component_names") /= CFG_get_size("gas_component_fractions")) &
           call ERR_show("dbd_1d: variables gas_component_names and gas_component_fracs have unequal size")

      call CFG_get("gas_component_names", gas_comp_names)
      call CFG_get("gas_component_fractions", gas_comp_fracs)

      ! Initialize gas and electric field module
      if (myrank == root) print *, "~~~~GAS initialize"
      call GAS_initialize(gas_comp_names, gas_comp_fracs, &
           CFG_get_real("gas_pressure"), CFG_get_real("gas_temperature"))

      if (myrank == root) print *, "~~~~electric field initialize"
      call EF_initialize()


       if (myrank ==root) print *, "~~~initialize the parameters for output"
       call OUT_init()

        ! Initialize variables
        sim_time     = 0.0_dp
        output_cntr  = 0
        info_cntr    = 0
        steps        = 0
        stepIons     = 0 
        dt           = CFG_get_real("sim_initial_dt")
        max_dt       = CFG_get_real("sim_max_dt")
        output_dt    = CFG_get_real("output_interval")
        sim_end_time = CFG_get_real("sim_end_time")
        err_goal     = CFG_get_real("sim_rel_error_goal")
        n_steps_apm  = CFG_get_int("apm_steps_between")
        apm_increase = CFG_get_real("apm_increase_factor")
        timeStepIonOverElec= CFG_get_int("part_timeStepRatio")

        doAverage    = CFG_get_logic("turner_average_active")
        stepToExcute = CFG_get_int("turner_execute_stepNum")
        stepToAverage= CFG_get_int("turner_average_stepNum")
            
        output_numerical_para = CFG_get_logic("sim_output_numerical_para")
        use_DBD         = CFG_get_logic("sim_DBD")

        !parameters for time evaluations of parameters during one peroid
        outputTimeEvaPra = CFG_get_logic("sim_output_time_eva_para")
        potPreStep       = 0.d0
        timeStartEva     = 0.d0
        flagfindTime     = .true.
        ave_num_timeEva  = CFG_get_int("sim_ave_num_timeEva")
        div_num_timeEva  = CFG_get_int("sim_div_num_one_peroid")
		sim_voltage_type = CFG_get_int("sim_voltage_type")
		select case(sim_voltage_type)
		case(1)!use AC voltage
			peroid   = 1.d0 / CFG_get_real("sim_voltage_AC_fre")
		case(2)!use pulse voltage
			peroid   = 1.d0 / CFG_get_real("sim_voltage_pulse_fre")
		end select
        outputMoment     = 0.d0
        cnrt             = 0

        do_verlet_bf_coll = CFG_get_logic("sim_do_verlet_before_coll")

        ! save and read particle info
        save_mid_info = CFG_get_logic("sim_save_mid_info")
        read_mid_info = CFG_get_logic("sim_read_mid_info")
        save_step_interval = CFG_get_int("sim_save_steps_interval")

! save time depedent parameters
        time_para_output_interval = CFG_get_int("sim_tim_para_out_interval")


        ! for convegence test
        factor_convegence = CFG_get_real("sim_conv_fac")
        call CFG_get("sim_read_file_name", read_file_name)
        dt =  dt  / factor_convegence
        stepToExcute = int(stepToExcute * factor_convegence)
        stepToAverage = int(stepToAverage * factor_convegence)

      select case (sim_type)
      case (MODEL_part)
         if (myrank == root) print *, "Reading crossection data: electrons and ions!"

         if (CFG_get_size("gas_crosssec_files") /= n_gas_comp) &
              call ERR_show("dbd_1d(Electrons): variables gas_crosssec_files and gas_component_fracs have unequal size")

         if (CFG_get_size("gas_crosssec_files_for_ions") /= n_gas_comp) &
              call ERR_show("dbd_1d(Ions): variables gas_crosssec_files and gas_component_fracs have unequal size")

         do nn = 1, n_gas_comp
            cs_file = CFG_get_string("gas_crosssec_files", nn)
            call CS_read_file("input/" // trim(cs_file), trim(gas_comp_names(nn)), 1.0_dp, &
                 gas_comp_fracs(nn) * GAS_get_number_dens(), CFG_get_real("part_max_energy_eV"))

            ! for ions
            cs_file_for_ions = CFG_get_string("gas_crosssec_files_for_ions", nn)
            call CS_read_file_for_ions("input/" // trim(cs_file_for_ions), trim(gas_comp_names(nn)), 1.0_dp, &
                 gas_comp_fracs(nn) * GAS_get_number_dens(), CFG_get_real("part_max_energy_eV"))
         end do

         if (myrank == root)  print *, "geting and wrighting crossection data: electrons!"
         call CS_get_cross_secs(cross_secs)    !Anbang: Here checked later, whether only needs for root
         if (myrank == root) then
            call CS_write_summary("output/" // trim(sim_name) // "_cs_summary.txt")
            call CS_write_all("output/" // trim(sim_name) // "_cs.txt")
         end if

         ! for ions
         if (myrank == root)  print *, "geting and wrighting crossection data:Ions!"
            call CS_get_cross_secs_for_ions(cross_secs_for_ions)
         if (myrank == root) then
            call CS_write_summary_for_ions("output/" // trim(sim_name) // "_cs_ions_summary.txt")
            call CS_write_all_for_ions("output/" // trim(sim_name) // "_cs_ions.txt")
         end if

         if (myrank == root)  print *, "~~~initializing particle module..."
         call PM_initialize(cross_secs, cross_secs_for_ions, ntasks)   
            !Anbang: here i need to remove updata density and efield routines,
            !only gives the velocity and position of electrons and ions
          !only gives the velocity and position of electrons and ions
         if (myrank == root) then               
                if (read_mid_info) then
                    print *, "~~Initialize the seed from an exsited particle file!"
                    call PC_read_elec_part_info(read_file_name)     
                    call PC_read_ion_part_info(read_file_name)
                else
                    print *, "~~~initialize the type of seeds and energy and the position from zero time"
                    call INIT_init_cfg()
                end if
         end if

        ! read necessary numberical parameters 
        !Anbang: i ignore this for convergence testing
        if (read_mid_info) then
            call OUT_read_mid_sim_info(sim_time, steps, stepIons, output_cntr, read_file_name)  
            ! here need to be careful, we need go to the next step
            sim_time = sim_time + dt
            steps = steps + 1
            stepIons = stepIons + int(1.d0/timeStepIonOverElec)
            ! only read at root, the surface charge
            if (use_DBD) call PM_read_surf_charge_info(myrank, root, read_file_name)
            stepToExcute = stepToExcute + steps
        end if
        
        if (myrank == root) print *, "distributing electrons and ions to", ntasks," nodes"        
        call shareParticlesElectrons(myrank, ntasks)
        call shareParticlesIons(myrank, ntasks)

        ! anbang: checking the initial weight and pos of particles
        ! we do not need it normally
        if (myrank == root) print *, "checking the weight and positions of particles"
        call PC_check_weight_particles()
        call PC_check_pos_particles()

        if (myrank == root) print *, "calculating the initial electron and ion density."
        call PM_set_elec_density()
        call PM_set_ion_density()
        call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root
        if (myrank == root) then  !only root calculating the efield, then it will broadcast to all nodes
            print *, "calculating the initial electric field."
        end if
        call PM_update_efield(sim_time, myrank, root)

        if (myrank == root) print *, "geting the initial accel of electron and ions."
        call PM_get_accel()

      case (MODEL_fluid)
         call CFG_get("fluid_input_file", fluid_input_file)
         call CFG_get("gas_mixture_name", gas_mix_name)
         call TD_read_file("input/" // trim(fluid_input_file), trim(gas_mix_name))

         call FL_init_cfg()
      end select
      print *, "The initialize_all is done!"

   if (myrank == root) then
        call OUT_write_coeffs(sim_name, sim_type)
        call cpu_time(time_start)
            
        ! Anbang: Here we write the parameters as a function of time
        outTimeFile = "output/" // trim(sim_name) // "_timePars.txt"
        open(unit = my_unit, file = outTimeFile , status = 'unknown')
        write (my_unit, *) "Here we output parameters information as a function of time"
        write (my_unit, *) "time", "potential"," n_elec", "n_ion","tot_curr(1)", "tot_curr(4)", "max_elec","max_ion", "sur_all"
        
        outTimeFile = "output/" // trim(sim_name) // "_fluxes.txt"
        open(unit = my_unit_flux, file = outTimeFile , status = 'unknown')
        write (my_unit_flux, *) "Here we output flux information as a function of time"
    end if


    prev_apm_elec_part = PC_getNumSimelecPartMPI()
    prev_apm_ion_part = PC_getNumSimionPartMPI() 

    do_apm = .false.
    do_apm_ion =.false.
    didAPMforElec = .false.
    didAPMforIons = .false.

   ! Here the simulation starts
   do
      if (sim_time > sim_end_time) exit

      n_elec_sim_sum = PC_getNumSimelecPartMPI()
      n_ion_sim_sum  = PC_getNumSimionPartMPI()
      n_elec_real_sum   = PC_getNumRealelecPartMPI()
      n_ion_real_sum    = PC_getNumRealionPartMPI()

      if (myrank == root) then
            call cpu_time(time_now)
            if (time_now - time_start > 10 * info_cntr) then
                if (sim_type == MODEL_part) then
                    write(*, "(I8,A,E8.2,A,E8.2,A,E8.2,A,E8.2,A,E8.2)") steps, " -- t = ", sim_time, ", dt = ", dt, &
                        ", wct = ", time_now - time_start, ", n_elecs = ", real(n_elec_sim_sum, dp), &
                        ",n_ions = ", real(n_ion_sim_sum, dp)
                else
                    write(*, "(I8,A,E8.2,A,E8.2,A,E8.2)") steps, " -- t = ", sim_time, ", dt = ", dt, &
                        ", wct = ", time_now - time_start
                end if
                info_cntr = info_cntr + 1
            end if
        end if

      if (sim_time >= output_cntr * output_dt) then   !finish, whether it is right
         output_cntr = output_cntr + 1
         call OUT_write_vars(sim_name, sim_type, output_cntr, sim_time, myrank, root, ntasks)
      end if

      select case (sim_type)
      case (MODEL_part)

        if (myrank == root) then
            !> Here we check whether we should use adaptive particles for ELECTRONS and Ions
            !> The first time step. no apm, ! initially, step/stepion = zero
            if (steps /= 0 .and. stepIons /= 0) then
                do_apm = (n_steps_apm > 0 .and. mod(steps, n_steps_apm) == 0) &
                    .or. n_elec_sim_sum > prev_apm_elec_part * apm_increase
                do_apm_ion = (n_steps_apm > 0 .and. mod(stepIons, n_steps_apm) == 0) &
                    .or. n_ion_sim_sum > prev_apm_ion_part * apm_increase
            end if
        end if

        !> broadcast the flag whether we need to do grid refinement
        call MPI_BCAST(do_apm, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(do_apm_ion, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)

        if (do_verlet_bf_coll) then
            ! here if we use verlet firstly, we need to update pos-accel-vel-collision seperately
            call PC_advance_elec_pos_verlet(dt)
            if (mod(steps,timeStepIonOverElec) == 0 ) then
                call PC_advance_ion_pos_verlet(dt * timeStepIonOverElec)  
                stepIons = stepIons + 1
            end if
            call PM_set_elec_density()
            call PM_set_ion_density()
            call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root
            call PM_update_efield(sim_time, myrank, root)
            call PM_get_accel()
            call PC_advance_elec_vel_verlet(dt)
            if (mod(steps,timeStepIonOverElec) == 0 ) &
               & call PC_advance_ion_vel_verlet(dt * timeStepIonOverElec) 
            call PC_elec_collide(dt)
             if (mod(steps,timeStepIonOverElec) == 0 ) &
               & call PC_ion_collide(dt * timeStepIonOverElec)    
        else
        print *," PM_advance"
            call PM_advance(dt)
            ! move the ions       
            if (mod(steps,timeStepIonOverElec) == 0 ) then
                call PM_ion_advance(dt * timeStepIonOverElec)  
                stepIons = stepIons + 1
            end if
        end if

        call PM_set_elec_density()
        call PM_set_ion_density()
         call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root  
         !Anbang: here sharing dens is very important for particle remaping, as we need the density to get the desired weight of particles
         call PM_shareDens(myrank, root)  ! broadcast to all the nodes

        !Anbang: Her we set if we need to use dbd, we have to calculate the surface charge
        if (use_DBD) then
            call PM_get_surface_charge_density(dt)

            call PM_collectsurfChargeAtRoot(myrank, root)  ! the charge is collected at root
        end if


       ! adjust the weight of electrons
        if (do_apm) then  
            call PM_adapt_weights()
            didAPMforElec = .true. 
        end if

        ! adjust the weight of ions
        if (do_apm_ion) then
            call PM_adapt_ion_weights()
            didAPMforIons = .true.    
        end if

        ! save the info of elec/ion and numberical parameters
        if (save_mid_info .and. mod(steps, save_step_interval) == 0) then
            if (myrank == root) print *, "Saving mid info at steps:", steps
            call PC_save_elec_part_info(myrank, root, ntasks, sim_name)
            call PC_save_ion_part_info(myrank, root, ntasks, sim_name)
            if (myrank == root) call OUT_save_mid_sim_info(sim_time, steps, stepIons, output_cntr, sim_name)
            if (use_DBD) call PM_save_surf_charge_info(myrank, root, sim_name)
            if (myrank == root) print *, "Saving info is done at steps:", steps
        end if


!         ! if there is no apm at all, we calculate the efield and accel directly
!         !BE CAREFUL that when we calculate the electric field, all the densities collected at the root and other nodes are set as zero
         if ((.not. didAPMforElec) .and. (.not. didAPMforIons) ) then
                call PM_update_efield(sim_time, myrank, root)
                !Anbang: if we do collisions between time step, we need to adjust accel here
                if (.not. do_verlet_bf_coll) call PM_correct_accel_and_vel(dt)
         end if


        ! here we check if we did apm for electrons or ions, if we did, we should adjust the accel for electrons and all ions
        if (didAPMforElec .or. didAPMforIons) then
            if (didAPMforElec) then
                n_elec_sim_sum     = PC_getNumSimelecPartMPI()
                n_elec_real_sum   = PC_getNumRealelecPartMPI()
                prev_apm_elec_part = n_elec_sim_sum
                if (myrank == root) write(*,'(A, E12.4, A, I0)') " After merging:(Electrons) real/sim part: ", &
                    n_elec_real_sum, " / ", int(n_elec_sim_sum)
            end if
!
            if (didAPMforIons) then
                n_ion_sim_sum     = PC_getNumSimionPartMPI()
                n_ion_real_sum    = PC_getNumRealionPartMPI()
                prev_apm_ion_part = n_ion_sim_sum
                if (myrank == root) write(*,'(A, E12.4, A, I0)') " After merging:(Ions) real/sim part: ", &
                n_ion_real_sum, " / ", int(n_ion_sim_sum)
            end if

            call PM_set_elec_density()
            call PM_set_ion_density()
            call PM_collectDensitiesAtRoot(myrank, root) ! collect densities at root
            call PM_update_efield(sim_time, myrank, root)
            !Anbang: we here get new accel
            if (do_verlet_bf_coll) then
                call PM_get_accel()
            else
                call PM_correct_accel_and_vel(dt)
            end if
            didAPMforElec = .false.
            do_apm  = .false.
            didAPMforIons = .false.
            do_apm_ion = .false.
        end if
        
        sim_time = sim_time + dt

         ! we here to average the para in a certain time, according to turner's code
         if (doAverage) then

            !Anbang: Here i output the parameters during 1 or 2 peroids, to show the time evaluation of parameters
            !> We sum the paramerters at the same moment during e.g. 100 time steps, then we average them to smooth the data
            if (steps > stepToExcute .and. outputTimeEvaPra .and. flagfindTime) then   ! we arrive the peroid steady state
                ! find the time when we start to calculate the data
                potential_left_bound = EF_get_potential_at_bound()
                if (potPreStep < 0.d0 .and. potential_left_bound >= 0.d0) then
                    timeStartEva = sim_time
                    flagfindTime = .false.
                else
                    potPreStep = potential_left_bound
                end if
            end if

            if (timeStartEva > 0.d0 .and. sim_time <= (timeStartEva + peroid * ave_num_timeEva)) then  ! here we find time when we can start and finish
                 outputMoment = timeStartEva + cnrt * peroid / div_num_timeEva
                 if (sim_time >= outputMoment) then
                    nn = mod(cnrt, div_num_timeEva) + 1
                    call PM_aver_para_one_peroid(nn, ave_num_timeEva, sim_time, myrank, root)
                    cnrt = cnrt + 1
                 end if
            end if


            if (steps >= stepToExcute .and. steps < (stepToExcute + stepToAverage))  then
                call PM_get_aver_para(sim_time, myrank, root)
                call Cal_aver_EEDF_elec(stepToAverage,  myrank, root, ntasks)

               ! calculate the particle number and mean weight per cell
                call PM_get_aver_num_and_weight_per_cell()
            end if

            if (steps == (stepToExcute + stepToAverage) .and. (myrank == root)) then
                call PM_output_aver_para_one_peroid(div_num_timeEva, sim_name, timeStartEva, peroid) !Output parameters during one peroid
                call PM_sheath_pos_cal_output(div_num_timeEva, sim_name, timeStartEva, peroid)  ! output sheath positions

                call PM_output_aver_para(sim_name, sim_time)
                call PM_output_aver_phys_para(sim_name, sim_time)
                call OUT_write_aver_EEDF_elec(sim_name, sim_time)

                ! output the averaged particle number and mean weight per cell
                call PM_output_aver_num_and_weight_per_cell(sim_name)
            end if       
            
            if (steps > (stepToExcute + stepToAverage)) go to 333   ! finish 
         end if
       
        ! we redistribute the particles
        call shareParticlesElectrons(myrank, ntasks)
        call shareParticlesIons(myrank, ntasks)
         
      case (MODEL_fluid)
         call FL_advance(sim_time, dt, max_dt, new_dt, err_goal)
         dt = new_dt
      end select

     !Anbang: here we output the parameters of calculation
      if (myrank == root) then
            if (output_numerical_para) then
                ! anbang: pass the old field, for calculating displacement current
                if (mod(steps + 1, time_para_output_interval)==0 .and. time_para_output_interval /= 1)  call EF_get_old_efield()
                if (mod(steps, time_para_output_interval)==0) then
                    potential_left_bound = EF_get_potential_at_bound()
					call EF_get_Gap_and_dielectric_potential(potential_gap , potential_dielectric)
                    ! calculate current
                    call PM_calculate_curr(dis_curr, conduct_curr, tot_curr,tot_curr_density,dt)
                    if (time_para_output_interval == 1) call EF_get_old_efield()  !if we sameple every time step, we pass the old field here
                    if (steps == 0) then ! if at time zero, we set the dis_curr = 0.0
                        dis_curr = 0.d0
                        tot_curr(:) = conduct_curr(1,:) + conduct_curr(2,:) 
                    end if

                    ! get electric field at index points
                    call EF_get_efield_index(efield_ix)

!                     write(my_unit, *) sim_time, potential_left_bound, efield_ix, &
!                         & dis_curr,tot_curr, conduct_curr(1,1:size(dis_curr)), conduct_curr(2, 1:size(dis_curr))

                    ! output time dependent parameters during a certain time_now
                    if (mod(steps, 1) == 0) then
                        call PM_get_max_dens(max_dens,sur_all)
                        write(my_unit, *) sim_time,potential_left_bound,potential_gap,potential_dielectric,n_elec_real_sum,&
                            & n_ion_real_sum,tot_curr_density(1), tot_curr_density(4), max_dens, sur_all
                            
                        call PM_out_flux(flux_all)
                        write(my_unit_flux,*) flux_all
                    end if

                end if
            end if
       end if

      steps = steps + 1     
      
        if (mod(steps,1000)==0) then
            tempElecNum = PC_getNumSimelecPartMPI()
            tempIonNum  = PC_getNumSimionPartMPI()
            if (myrank == root) then
                print *, "**********************"
                print *, "steps =", steps
                print *, "the nubmer of electrons/ions:", tempElecNum, tempIonNum   
                print *, "**********************"
            end if
        end if
   end do
   333 continue

   if (myrank == root) then
   close(my_unit)
   close(my_unit_flux)
   end if

   call MPI_FINALIZE(ierr)

   print *, "The simulation has finished!"
   stop

contains

   !> Create the parameters and their default values for the simulation
   subroutine create_sim_config()
	! Field emission parameters
      call CFG_add("field factor", 5.0D1, "The local field enhancement factor due to roughness of surface")
      call CFG_add("work fuction", 4.0D0, "The energy(eV) for electron can escape trap") 
      ! General simulation parameters
      call CFG_add("sim_type", "fluid", "The type of simulation to run, options are particle, fluid_min, fluid_ee")
      call CFG_add("sim_end_time", 1.0D-9, "The desired endtime in seconds of the simulation")
      call CFG_add("sim_initial_dt", 1.0D-13, "The initial/fixed timestep in seconds")
      call CFG_add("sim_max_dt", 1.0D-11, "The maximal timestep in seconds")
      call CFG_add("sim_rel_error_goal", 1.0D-4, "Desired relative error in the solution between consecutive steps")
      call CFG_add("sim_name", "my_sim", "The name of the simulation")

      ! Grid parameters
      call CFG_add("grid_num_points", 2000, "The number of grid cells")
      call CFG_add("grid_delta_x", 2.0d-6, "The length of a grid cell")

      ! Gas parameters
      call CFG_add("gas_pressure", 1.0D0, "The gas pressure (bar)")
      call CFG_add("gas_temperature", 293.0D0, "The gas temperature (Kelvin)")
      call CFG_add("gas_mixture_name", "Helium", "The name of the gas mixture used")
      call CFG_add("gas_component_names", (/"He"/), "The names of the gases used in the simulation", .true.)
      call CFG_add("gas_crosssec_files", (/"elecCS_He.txt"/), &
           & "Electrons: the files in which to find cross section data for each gas", .true.)
      call CFG_add("gas_component_fractions", (/1.0_dp/), &
           & "The partial pressure of the gases (as if they were ideal gases)", .true.)
      ! call CFG_add("gas_crosssec_scaling", 1.0D0, "Scale factor for the cross sections in the input data", .true.)

      call CFG_add("gas_crosssec_files_for_ions", (/"ionCS_He.txt"/), &
           & "Ions: the files in which to find cross section data for each gas", .true.)

      ! Electric field parameters
      CALL CFG_add("sim_elec_times", (/0.0d-9/), "The times at which the  voltage/electric field is specified.", .true.)
      CALL CFG_add("sim_applied_efield", (/1.0D7/), "The electric field at the specified times.", .true.)
      call CFG_add("sim_constant_efield", .false., "Whether the electric field is kept constant")

      !Voltage parameters: Anbang: if we want to use voltage as boundary conditions
      call CFG_add("sim_use_voltage", .false., "Whether to use voltage")
      CALL CFG_add("sim_applied_voltage", (/10.0D3/), "The potential used at left boundary at the specified times.", .true.)
      CALL CFG_add("sim_voltage_type", 0, "The type of applied voltage: 0- DC, 1- AC,2- pulse")

      ! amplitude and frequency for AC
      CALL CFG_add("sim_voltage_AC_max ", 4.5d2, "The amplitude of AC, unit: V")
      CALL CFG_add("sim_voltage_AC_fre ", 13.56d6, "The frequency of AC, unit: Hz")
	   
	  ! amplitude and frequency for pulse
	  CALL CFG_add("sim_voltage_pulse_max ", 4.5d2, "The amplitude of pulse, unit: V")
      CALL CFG_add("sim_voltage_pulse_fre ", 2d6, "The frequency of pulse, unit: Hz")
	  CALL CFG_add("sim_voltage_rising_edge ", 1.5d-8, "The amplitude of AC, unit: s")
      CALL CFG_add("sim_voltage_falling_edge ", 1.5d-8, "The frequency of AC, unit: s")
	  CALL CFG_add("sim_voltage_platform ", 1.5d-8, "The amplitude of AC, unit: s")
	  
      ! Initial conditions
      call CFG_add("init_cond_name", "gaussian", "The type of initial condition")
      call CFG_add("init_elec_dens", 1.0d17 , "The number of initial ion pairs")
      call CFG_add("init_elec_energy", 1.0D0 , "The initial energy of the electrons in eV")
      call CFG_add("init_ion_energy", 0.0D0 , "The initial energy of the ions in eV")

      call CFG_add("init_elec_tem", 3.d4 , "The initial temperature of the electrons in K")
      call CFG_add("init_ion_tem", 0.0D0 , "The initial temperature of the ions in K")
      call CFG_add("init_use_temInK", .false., "whether we use temperature as initial parameter")

      call CFG_add("init_rel_pos", 0.5D0, "The relative position of the initial seed")
      call CFG_add("init_width", 25.0d-6, "The standard deviation used for Gaussian initial profiles")
      call CFG_add("init_background_density", 0.0D0, "The background ion and electron density in 1/m^3")
      call CFG_add("initial_weight", 1.0D0, "The initial weight of particles (electrons/ions)")

      ! Output parameters
      call CFG_add("output_interval", 1.0D-10, "The timestep for writing output")

      call CFG_add("apm_part_per_cell", 1.0d2, &
         "The desired number of particles per cell")
      call CFG_add("apm_vel_rel_weight", 1.0d-6, &
         "Relative weight of vel. in the k-d tree compared to position")
      call CFG_add( "apm_steps_between", 100, &
         "Adapt weight every apm_steps_between steps")
      call CFG_add("apm_increase_factor", 1.2_dp, &
         "Adapt weight if #part increases by this factor")

      ! Particle model related parameters
      call CFG_add("part_lkptbl_size", 20*1000, "The size of the lookup table for the collision rates")
      call CFG_add("part_transverse_area", 5.0D-12, "The transverse area of the simulation, to convert to m^3 etc.")
      call CFG_add("part_max_number_of", 25*1000*1000, "The maximum number of particles allowed per task(electrons)")
      call CFG_add("ion_part_max_number_of", 25*1000*1000, "The maximum number of particles allowed per task(ions)")
      call CFG_add("part_max_energy_eV", 1000.0D0, "The maximum energy in eV for particles in the simulation")
      call CFG_add("part_tableSize", 50000, "The table size to use in the particle model for lookup tables")
      call CFG_add("part_timeStepRatio", 1, "The time step ratio between ions and electrons")
      call CFG_add("part_ionMassRatio", 1.d0, "The mass ratio between ions and standard atoms")
      call CFG_add("part_atomMassRatio", 1.d0, "The mass ratio between simulated atoms and standard atoms")
      call CFG_add("sim_useBgVel", .false., "Whether the velocity of bg gas is considered!")
      call CFG_add("sim_frame_type", 1, "The frame we used: 1 - lab frame; 2 - Center of mass frame")
      call CFG_add("sim_merge_scheme", 1, "randomly choose one particle position; 2 - weighted average position")


      ! parameters for comparing with turner's results
      call CFG_add("turner_average_stepNum", 12800, "The step nums for averaging the valuse in turner's code")
      call CFG_add("turner_execute_stepNum", 512000, "The step nums when to excute in turner's code")
      call CFG_add("turner_average_active", .false., "Whether we active average routines")

      ! parameters for outputing EEDF of electrons
      call CFG_add("output_num_bins", 200, "The number of bins to use for histograms")
      call CFG_add("output_eedf_min_max_fields", (/0.0d0, 1.0d10/), "The field range(s) over which to get the EEDF")
      call CFG_add("output_eedf_eV_range", (/0.0d0, 1.0d2/), &
                    "The energy range over which to get the EEDF")
      

      ! General fluid model parameters
      call CFG_add("fluid_use_energy", .false., "Whether to use an energy equation")
      call CFG_add("fluid_use_en_mob", .false., "Whether to use energy dependent mobility")
      call CFG_add("fluid_use_en_dif", .false., "Whether to use energy dependent diffusion coefficient")
      call CFG_add("fluid_use_en_src", .false., "Whether to use energy dependent source term")
      call CFG_add("fluid_input_file", "transport_data_nitrogen.txt" , "The input file for the fluid models")
      call CFG_add("fluid_lkptbl_size", 1000, "The transport data table size in the fluid model")
      call CFG_add("fluid_lkptbl_max_efield", 3.0d7, "The maximum electric field in the fluid model coefficients")
      call CFG_add("fluid_lkptbl_max_energy", 1.0d2, "The maximum mean energy in eV in the fluid model coefficients")
      call CFG_add("fluid_small_density", 1.0d0, "Regularization density to compute a mean energy")

      call CFG_add("fluid_en_fld", "energy[eV]_vs_efield[V/m]", "The name of the energy vs efield list")
      call CFG_add("fluid_en_mob", "energy[eV]_vs_mu[m2/Vs]", "The name of the mobility coefficient")
      call CFG_add("fluid_en_dif", "energy[eV]_vs_dif[m2/s]", "The name of the diffusion coefficient")
      call CFG_add("fluid_en_alpha", "energy[eV]_vs_alpha[1/m]", "The name of the eff. ionization coeff.")
      call CFG_add("fluid_en_eta", "energy[eV]_vs_eta[1/m]", "The name of the eff. attachment coeff.")
      call CFG_add("fluid_en_loss", "energy[eV]_vs_loss[eV/s]", "The name of the energy loss coeff.")

      call CFG_add("fluid_fld_mob", "efield[V/m]_vs_mu[m2/Vs]", "The name of the mobility coefficient")
      call CFG_add("fluid_fld_en", "efield[V/m]_vs_energy[eV]", "The name of the energy(fld) coefficient")
      call CFG_add("fluid_fld_dif", "efield[V/m]_vs_dif[m2/s]", "The name of the diffusion coefficient")
      call CFG_add("fluid_fld_alpha", "efield[V/m]_vs_alpha[1/m]", "The name of the eff. ionization coeff.")
      call CFG_add("fluid_fld_eta", "efield[V/m]_vs_eta[1/m]", "The name of the eff. attachment coeff.")
      call CFG_add("fluid_fld_loss", "efield[V/m]_vs_loss[eV/s]", "The name of the energy loss coeff.")

      ! Add parameters for DBD
      call CFG_add("sim_DBD", .false., "Whether we have dielectrics")
      call CFG_add("sim_len_del", (/0.d0, 0.d0/), "the length of the dielectric in both sides")
      call CFG_add("sim_relPer_del", (/1.d0, 1.d0/), "the relative permittivity of the dielectric in both sides")
      call CFG_add("sim_bd_elec", (/1, 1/), &
            "the boudary conditions for electrons on both sides: 1- absorbtion; 2 - reflection")
      call CFG_add("sim_bd_ion", (/1, 1/), &
            "the db for ions on both sides: 1- absorbtion; 2 - reflection, 3- secondary electton emission")
      call CFG_add("sim_co_secondary_elec_emmison", (/0.d0, 0.d0/), &
            "the secondary electron emission coeffients on the both sides")
      call CFG_add("sim_secElecEn", 2.d0, &
            "the energy of secondary electrons in eV")
      call CFG_add("sim_ref_coeff_elec", (/0.d0, 0.d0/), &
            "the reflection coeffients of electrons on the both sides")

      call CFG_add("sim_output_numerical_para", .false., &
            "whether we output the potential as a funtion of time in the main program.")

      call CFG_add("sim_rot_scheme", 1, &
            "which rotation scheme we choose for scattering: 1- Yousfi, 2- Detlef, 3- Donko")

      call CFG_add("sim_table_create_scheme", 1, &
            "1: we use velocity as the column of the lookup table; 2: we use com energy")

      call CFG_add("sim_bw_scheme", 1, &
            "the scheme for backward scattering, 1 - turner(charge exchange); 2 - transfer between com and lab frame")

      call CFG_add("sim_do_verlet_before_coll", .true., &
            "Whether we do verlet updating completly before do collisions of particles")

       !Anbang: The following are parameters for outputing time evaluation of paramters during one peroid after peroidic steady state
      call CFG_add("sim_div_num_one_peroid", 100, &
            "How many parts are divided of one peroid, which means we output this times of paramters during one peroid")
      call CFG_add("sim_output_time_eva_para", .true., &
            "whether we output time evaluated paramerters after steady state")
        call CFG_add("sim_ave_num_timeEva", 100, &
            "We average this number of peroid to get a good fit")

        ! save and read particle info during simulations
        call CFG_add("sim_save_mid_info", .false.,  &
                & "Here we decide whether we should save particle info to an existed file")
        call CFG_add("sim_read_mid_info", .false.,  &
                & "Here we decide whether we should read particle info from an existed file")
        call CFG_add("sim_save_steps_interval", 100,  &
                & "Here we decide how many simulation time steps we should save the mid particle info")

        ! some parameters for convergence testing
        call CFG_add("sim_read_file_name", "my_sim", "The name of files we want to read for initializing particles")
        call CFG_add("sim_conv_fac", 1.d0, "The factor of time step that want to test convegences")

 !time parameters output interval
         call CFG_add("sim_tim_para_out_interval", 1, &
            "time parameters output interval")

   end subroutine create_sim_config

end program dbd_1d
