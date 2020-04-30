! Authors: Jannis Teunissen, concepts based on work of Chao Li, Margreet Nool
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

!> Core particle module. The user has to supply his own routines for customization.
module m_particle_core
   use m_types
   use m_lookup_table

   implicit none
   private

   !> The number of active particles in the simulation (electrons, ions)
   integer     :: PC_num_part, PC_num_ion_part

   !> The electron type
   type PC_part_t
      real(dp) :: T_left   ! The time until the next timestep
      real(dp) :: x(3)     ! The position
      real(dp) :: v(3)     ! The velocity
      real(dp) :: a(3)     ! The electric field's acceleration
      real(dp) :: weight   ! The weight factor of the (super)particle
      logical  :: live     ! If .false. the particle will not be used (anymore)
   end type PC_part_t

   !****************************************************************
      !> The ion type
   type PC_ion_part_t
      real(dp) :: T_left   ! The time until the next timestep
      real(dp) :: x(3)     ! The position
      real(dp) :: v(3)     ! The velocity
      real(dp) :: a(3)     ! The electric field's acceleration
      real(dp) :: weight   ! The weight factor of the (super)particle
      logical  :: live     ! If .false. the particle will not be used (anymore)
   end type PC_ion_part_t
   !****************************************************************

   !> The mass of the particles in the simulation (electrons, ions)
    real(dp) :: PC_particle_elec_mass, PC_particle_ion_mass, PC_reduced_mass,PC_Q_bd(2)

   !> The maximum number of particles (electrons, ions)
   integer  :: PC_max_num_part, PC_max_num_ion_part



   !> The temperature of gas
   real(dp) :: PC_gas_tem
   
   !> The mass of simulated atoms
   real(dp) :: PC_particle_atom_mass
  
   !> Whether we use background gas velocity
   logical  :: PC_useBgVel
   
   !Anbang: we decide use com/lab frame for scattering
   integer, parameter   ::  PC_lab_frame = 1, PC_com_frame = 2
   integer              ::  PC_frame_type

   !> The length of the domain
   real(dp) :: PC_domainLength

   !> The length of the dielectrics in both sides
   real(dp) :: PC_len_del(2)
   
   !> The positions of left and right boundaries
   real(dp) :: PC_posBD(2)

   !> The type of boundary conditions for electrons and ions
   integer :: PC_bdType_Ion(2), PC_bdType_elec(2)

   !> the coefficient of sendary electrons emmision
   real(dp) :: PC_coSecElec(2)

   !> the electron reflection coefficients on both sides
   real(dp) :: PC_coElecRefl(2)

   !> the energy of electrons
   real(dp) :: PC_secElecEn

   real(dp), allocatable :: PC_ionization_number(:)

   real(dp)         :: PC_delta_x
   integer          :: PC_grid_size

   integer          :: PC_rot_scheme

   integer          :: PC_flagTable
   integer, parameter ::  flag_vel =1, flag_en = 2 

   integer          :: PC_BWflag

  ! Anbang: the meerge shcme of particle positon
   integer              :: PC_merge_scheme

   ! collision for electrons
   type PC_coll_t
      integer                :: num                    ! The number of different collisions
      real(dp)               :: max_rate, inv_max_rate ! Maximum collision rate and inverse
      integer, allocatable   :: types(:)               ! The types of the collisions that can occur
      real(dp), allocatable  :: buff(:)                ! Temporary storage for the collision rates
      real(dp), allocatable  :: spec_values(:)         ! Special value used by this collision
      real(dp), allocatable  :: counts(:)              ! To collect statistics during runtime
      type(LT_type), pointer :: rate_lt, sum_rate_lt   ! Lookup table with collision rates
   end type PC_coll_t

   type(PC_coll_t) :: PC_coll


   ! collision for Ions
   type PC_ion_coll_t
      integer                :: num                    ! The number of different collisions
      real(dp)               :: max_rate, inv_max_rate ! Maximum collision rate and inverse
      integer, allocatable   :: types(:)               ! The types of the collisions that can occur
      real(dp), allocatable  :: buff(:)                ! Temporary storage for the collision rates
      real(dp), allocatable  :: spec_values(:)         ! Special value used by this collision
      real(dp), allocatable  :: counts(:)              ! To collect statistics during runtime
      type(LT_type), pointer :: rate_lt, sum_rate_lt   ! Lookup table with collision rates
   end type PC_ion_coll_t

   type(PC_ion_coll_t) :: PC_ion_coll

   !> The list that contains all the particles in the simulation (electrons, ions)
   type(PC_part_t), allocatable     :: PC_particles(:)
   type(PC_part_t), allocatable :: PC_particles_collect(:)

   type(PC_ion_part_t), allocatable :: PC_ion_particles(:)   
   type(PC_ion_part_t), allocatable :: PC_ion_particles_collect(:)


   !> MPI type for sending the 'electrons/ ions' between tasks.
   integer :: PC_elecpartTypeMPI, PC_ionpartTypeMPI

  !> old accel of particles, useful for verlet scheme
   real(dp), allocatable :: PC_elec_old_accel(:,:), PC_ion_old_accel(:,:)


   interface
      subroutine if_ipart(my_part)
         import
         type(PC_part_t), intent(in) :: my_part
      end subroutine if_ipart

      subroutine if_iint2_ipart(my_part, my_int_1, my_int_2)
         import
         type(PC_part_t), intent(in) :: my_part
         integer, intent(in)         :: my_int_1, my_int_2
      end subroutine if_iint2_ipart

      subroutine if_ovec3(my_vec)
         import
         real(dp), intent(out) :: my_vec(3)
      end subroutine if_ovec3

      subroutine if_ipart_ovec3(my_part, my_vec)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_vec(3)
      end subroutine if_ipart_ovec3
   end interface

!> copy from jannis's new codes: 20140916
  abstract interface

     subroutine p_to_r3_p(my_part, my_vec)
       import
       type(PC_part_t), intent(in) :: my_part
       real(dp), intent(out)       :: my_vec(3)
     end subroutine p_to_r3_p

     function p_to_r3_f(my_part) result(my_vec)
       import
       type(PC_part_t), intent(in) :: my_part
       real(dp)                    :: my_vec(3)
     end function p_to_r3_f

     real(dp) function p_to_r_f(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_r_f

     ! for ions
     real(dp) function p_to_r_f_ion(my_part)
       import
       type(PC_ion_part_t), intent(in) :: my_part
     end function p_to_r_f_ion

     logical function p_to_logic_f(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_logic_f

  end interface

   !!!> For ions
   interface
      subroutine if_ion_ipart(my_ion_part)
         import
         type(PC_ion_part_t), intent(in) :: my_ion_part
      end subroutine if_ion_ipart

      subroutine if_ion_ipart_ovec3(my_ion_part, my_ion_vec)
         import
         type(PC_ion_part_t), intent(in) :: my_ion_part
         real(dp), intent(out)       :: my_ion_vec(3)
      end subroutine if_ion_ipart_ovec3
   end interface

   procedure(if_iint2_ipart), pointer :: PC_pptr_coll_callback => null()
   procedure(if_ovec3), pointer :: PC_pptr_bg_vel_sampler => null()

   ! Types
   public :: PC_part_t, PC_ion_part_t
   public :: PC_coll_t

   ! Procedures
   public :: PC_initialize
   public :: PC_advance, PC_advance_ions
   public :: PC_create_part, PC_create_ion_part
   public :: PC_periodify, PC_ion_periodify
   public :: PC_set_part !, PC_set_ion_part
   public :: PC_get_part !, PC_get_ion_part
   public :: PC_get_num_sim_part, PC_get_num_sim_ion_part
   public :: PC_get_num_real_part, PC_get_num_real_ion_part
   public :: PC_set_accel, PC_set_ion_accel
   public :: PC_correct_new_accel, PC_correct_new_ion_accel
   public :: PC_get_num_colls
   public :: PC_get_colls
   public :: PC_get_coll_count
   public :: PC_reset_coll_count
   public :: PC_loop_part, PC_loop_ion_part
   public :: PC_merge_and_split, PC_ion_merge_and_split
   public :: PC_merge_part_rxv, PC_ion_merge_part_rxv
   public :: PC_split_part, PC_ion_split_part
   public :: PC_check_pos_particles
   public :: PC_check_weight_particles
    public :: PC_histogram,PC_curr_cal
   public :: PC_output_ionization_number

   !MPI
   !public :: shareParticlesIons,  shareParticlesElectrons
   public :: shareParticlesElectrons, shareParticlesIons
   public :: PC_getNumRealelecPartMPI, PC_getNumSimelecPartMPI
   public :: PC_getNumRealionPartMPI, PC_getNumSimionPartMPI

   ! Verlet scheme 2
   public :: PC_advance_elec_pos_verlet, PC_advance_elec_vel_verlet
   public :: PC_advance_ion_pos_verlet, PC_advance_ion_vel_verlet
   public :: PC_elec_collide, PC_ion_collide

   public :: PC_save_elec_part_info, PC_save_ion_part_info
   public :: PC_read_elec_part_info, PC_read_ion_part_info

   public :: PC_cal_elec_num_and_weight_per_cell, PC_cal_ion_num_and_weight_per_cell

contains

   !> Initialization routine for the particle module
   subroutine PC_initialize(mass, massIon, nPartMax, nPartMaxIon, cross_secs, cross_secs_for_ions, lookup_tableSize, &
              max_en_eV, delta_x, grid_size, gas_tem, massAtom, useBgVel, frame_type, len_del, bd_elec, bd_ion, coSecElec, &
                elecEn, rot_scheme, flag_table, coElecRefl, schemeBWSC, scheme_Merge)
      use m_cross_sec
      use m_cross_sec_for_ions
      use m_error
      use m_units_constants

      type(CS_type), intent(in) :: cross_secs(:)
      type(CS_ion_type), intent(in) :: cross_secs_for_ions(:)
      integer, intent(in)       :: nPartMax, lookup_tableSize, nPartMaxIon
      real(dp), intent(in)      :: mass, max_en_eV,massIon
      real(dp), intent(in)      :: delta_x
      integer,  intent(in)      :: grid_size 
      real(dp), intent(in)      :: massAtom, gas_tem
      logical                   :: useBgVel
      integer                   :: frame_type
      real(dp), intent(in)      :: len_del(2)   ! the length of delectrics on bothside
      integer, intent(in)       :: bd_elec(2), bd_ion(2) 
      real(dp), intent(in)      :: coSecElec(2)  ! the secondary electron emisson coefficients
      real(dp), intent(in)      :: elecEn
      integer                   :: rot_scheme, flag_table
      real(dp), intent(in)      :: coElecRefl(2)
      integer                   :: schemeBWSC
      integer, intent(in)       :: scheme_Merge
      
      PC_particle_elec_mass = mass
      PC_max_num_part  = nPartMax
      allocate( PC_particles(PC_max_num_part) )
      PC_num_part      = 0
    PC_Q_bd=0.D0
      PC_gas_tem = gas_tem
      print *, "the temperature of gas:", PC_gas_tem
      PC_particle_atom_mass= massAtom
      PC_useBgVel = useBgVel
      PC_frame_type = frame_type
      PC_secElecEn =  elecEn
      PC_rot_scheme = rot_scheme
      PC_flagTable  = flag_table
      PC_coElecRefl = coElecRefl
      PC_BWflag     = schemeBWSC   ! the flag for backwardscattering: 1 turner; 2 com/lab frame transfer
      PC_merge_scheme = scheme_Merge

     ! print *, "PC_coElecRefl = ", PC_coElecRefl
      
      PC_len_del = len_del
      PC_domainLength = delta_x * (grid_size - 1)
      PC_delta_x  = delta_x
      PC_grid_size = grid_size
      print *, "domainLength =", PC_domainLength      

      allocate(PC_ionization_number(grid_size))
      PC_ionization_number = 0.d0
       
      PC_posBD(1) = PC_len_del(1)
      PC_posBD(2) = PC_domainLength - PC_len_del(2)

     ! print *, "PC_posBD positions are:", PC_posBD

      PC_bdType_elec = bd_elec  
      PC_bdType_Ion  = bd_ion

      PC_coSecElec  = coSecElec
      
      ! for ions
      PC_particle_ion_mass = massIon
      PC_max_num_ion_part  = nPartMaxIon
     ! print *, "The mass of ion and atoms are:", PC_particle_ion_mass, PC_particle_atom_mass

      ! the reduced mass in the com frame
      PC_reduced_mass =  PC_particle_ion_mass * PC_particle_atom_mass / &
                    (PC_particle_ion_mass + PC_particle_atom_mass)
     ! print *, "The reduced mass is:", PC_reduced_mass
     ! print *, "The flag is:", PC_flagTable
    
      allocate( PC_ion_particles(PC_max_num_ion_part))
      PC_num_ion_part = 0

      allocate(PC_elec_old_accel(3, PC_max_num_part))
      allocate(PC_ion_old_accel(3, PC_max_num_ion_part))
      PC_elec_old_accel = 0.d0
      PC_ion_old_accel  = 0.d0

      if (size(cross_secs) > 0) then
         ! Create the lookup table for the collision rates, minimum velocity is 0.0
         call create_coll_rate_table(cross_secs, 0.d0, max_en_eV, lookup_tableSize)
         allocate(PC_coll%counts(PC_coll%num))
         PC_coll%counts(:) = 0

         print *, "Electron: 1/collision rate", PC_coll%inv_max_rate
      else
         PC_coll%num = 0
      end if

      ! for ions as well
      if (size(cross_secs_for_ions) > 0) then
         ! Create the lookup table for the collision rates for ions, minimum energy is 0.0
         !Anbang: Here we creat a table that the collision rate is a function of velocity in the COM frame
         ! BE careful here: !!!!! To compare with turner's benchmark
         select case(PC_flagTable)
         case(flag_vel)
            call create_coll_rate_table_for_ions_input_vel(cross_secs_for_ions, 0.d0, max_en_eV, lookup_tableSize)
         case(flag_en)
            call create_coll_rate_table_for_ions_input_EneV(cross_secs_for_ions, 0.d0, max_en_eV, lookup_tableSize)
         end select

         allocate(PC_ion_coll%counts(PC_ion_coll%num))
         PC_ion_coll%counts(:) = 0

         print *, "Ion: 1/collision rate", PC_ion_coll%inv_max_rate
      else
         PC_ion_coll%num = 0
      end if

!       call createElecPartTypeMPI(PC_elecpartTypeMPI)
!       call createIonPartTypeMPI(PC_ionpartTypeMPI)

   end subroutine PC_initialize

    subroutine PC_curr_cal(grid_curr_out)
    real(dp), intent(out) :: grid_curr_out(2)
    grid_curr_out=PC_Q_bd
    PC_Q_bd=0.d0
    end subroutine PC_curr_cal
   !> Loop over all the particles, and for each set the time left 'dt' for this step.
   !! Then they are fed to the move_and_collide() routine, which advances them in time
   subroutine PC_advance(dt)
      real(dp), intent(IN)      :: dt
      integer                   :: ll
      real(dp)                  :: tempDt
        
      PC_ionization_number = 0.d0

      PC_particles(1:PC_num_part)%t_left = dt

      if (PC_coll%num > 0) then ! There are collisions
         ll = 1
         do while (ll <= PC_num_part)
            call move_and_collide(ll)
            ll = ll + 1
         end do
         call remove_dead_particles()
      else                      ! There are no collisions
         do ll = 1, PC_num_part
            tempDt = PC_particles(ll)%t_left
            call advance_particle(ll, PC_particles(ll)%t_left)
            !> boundary conditions: absorption
            if (isOutOfGas(PC_particles(ll)%x)) then
                call PC_boundary_for_electron(ll, PC_particles(ll)%x(3), tempDt)
            !    print *, "One electron is out:", PC_particles(ll)%x
            end if
         end do
         call remove_dead_particles()
      end if
     ! print *, "Min and Max position =", minval(PC_particles(:)%x(3)), maxval(PC_particles(:)%x(3))

     PC_ionization_number = PC_ionization_number / dt

   end subroutine PC_advance

   !> Anbang: For ions
    !>Loop over all the particles, and for each set the time left 'dt' for this step.
   !! Then they are fed to the move_and_collide() routine, which advances them in time
   subroutine PC_advance_ions(dt)
      real(dp), intent(IN)      :: dt
      integer                   :: ll
      real(dp)                  :: tempDt

      PC_ion_particles(1:PC_num_ion_part)%t_left = dt

      if (PC_ion_coll%num > 0) then ! There are collisions
         ll = 1
         do while (ll <= PC_num_ion_part)
            call move_and_collide_for_ions(ll)
            ll = ll + 1
         end do
         call remove_dead_ion_particles()
      else                      ! There are no collisions
         do ll = 1, PC_num_ion_part
            tempDt = PC_ion_particles(ll)%t_left
            call advance_ion_particle(ll, PC_ion_particles(ll)%t_left)
            
            !> boundary conditions
            if (isOutOfGas(PC_ion_particles(ll)%x)) then
                call PC_boundary_for_ion(ll, PC_ion_particles(ll)%x(3), tempDt)
            end if
         end do
         call remove_dead_ion_particles()
      end if

   end subroutine PC_advance_ions

   !> Perform a collision for an electron, either elastic, excitation, ionizationCollision,
   !! attachment or null.
   subroutine move_and_collide(ll)
      use m_cross_sec
      use m_utils

      integer, intent(in)  :: ll
      integer              :: cIx, cType
      double precision     :: coll_time, new_vel, tempDt


      coll_time = sample_coll_time()

      do while (coll_time < PC_particles(ll)%t_left .and. PC_particles(ll)%live)      
         call advance_particle(ll, coll_time)   ! Set x,v at the collision time

        !Anbang: Here i added another check, because it is possible that the electron is already out of domain here
        !If it is the case, we let it out without following collision
        if (isOutOfGas(PC_particles(ll)%x)) then
            call PC_boundary_for_electron(ll, PC_particles(ll)%x(3), coll_time)
        else
            new_vel = sqrt(sum(PC_particles(ll)%v**2))
            cIx     = get_coll_index(new_vel)

            if (cIx > 0) then
                cType               = PC_coll%types(cIx)
                PC_coll%counts(cIx) = PC_coll%counts(cIx) + PC_particles(ll)%weight

                ! Perform the corresponding collision
                select case (cType)
                case (CS_attach_t)
                    call attach_collision(ll)
                case (CS_elastic_t)
                    call elastic_collision_COM(ll, cIx)
                    !call elastic_collision_COM_rotation(ll, cIx)
                case (CS_excite_t)
                    call excite_collision(ll, cIx, new_vel)
                case (CS_ionize_t)
                    call ionization_collision(ll, cIx, new_vel)
                end select
            end if

            coll_time = sample_coll_time() ! Set the next collision time
        end if

         goto 444   !Anbang: i test if there is only one collision per time step
      end do
      
      444 continue
 
      !> Anbang: here is very important for the bourndaty conditions of electrons: likes reflection
      tempDt = PC_particles(ll)%t_left
 
      ! Update the particle position and velocity to the next timestep
      if (PC_particles(ll)%live) call advance_particle(ll, PC_particles(ll)%t_left)
      
      !> boundary conditions: absorption
      if (isOutOfGas(PC_particles(ll)%x)) then
          call PC_boundary_for_electron(ll, PC_particles(ll)%x(3), tempDt)
          !print *, "One electron is out:", PC_particles(ll)%x
      end if

   end subroutine move_and_collide


   !> Perform a collision for an ion, elastic, backwardscattering and null
   !> Anbang: Here we should pay attention, becasue the mass of ions and melecular is similar
   !>         We should use the relative velocity between them.
   !> Be careful: the cross section in turner's paper is a function of center of mass energy!!!! 
   subroutine move_and_collide_for_ions(ll)
      use m_cross_sec_for_ions
      use m_units_constants

      integer, intent(in)  :: ll
      integer              :: cIx, cType
      double precision     :: coll_time, tempDt
      real(dp)             :: bg_vel(3),com_vel(3), rel_vel
      coll_time = sample_coll_time_for_ions()

      do while (coll_time < PC_ion_particles(ll)%t_left .and. PC_ion_particles(ll)%live)

        call advance_ion_particle(ll, coll_time)   ! Set x,v at the collision time

         ! check if the particle is out or not, if it is out of domain, we delete it without collision
        if (isOutOfGas(PC_ion_particles(ll)%x)) then
            call PC_boundary_for_ion(ll, PC_ion_particles(ll)%x(3), coll_time)
        else
            !Anbang: Here is not clear how turner set the background gas velocities
            if (PC_useBgVel) then
                call creatANeutralVel(bg_vel)
            else
                bg_vel = 0.d0
            end if

            ! Anbang: Here we set the COM velocity
            com_vel = (PC_particle_ion_mass * PC_ion_particles(ll)%v + PC_particle_atom_mass * bg_vel) / &
                    (PC_particle_ion_mass + PC_particle_atom_mass)
            

            !relative velocity
            rel_vel = sqrt(sum((PC_ion_particles(ll)%v - bg_vel)**2))
        
            select case(PC_flagTable)
            case(flag_vel)
                cIx     = get_coll_index_for_ions(rel_vel)
            case(flag_en)
                ! velocity to com energy, to the look up table
                rel_vel = vel_to_en( rel_vel, PC_reduced_mass) / UC_elec_volt
                cIx    = get_coll_index_for_ions_formular(rel_vel)
            end select    
                  
            if (cIx > 0) then
                cType               = PC_ion_coll%types(cIx)
                PC_ion_coll%counts(cIx) = PC_ion_coll%counts(cIx) + PC_ion_particles(ll)%weight

                ! Perform the corresponding collision
                select case (cType)
                case (CS_ion_elastic_t)
                    select case (PC_frame_type)
                    case (PC_lab_frame)
                        call elastic_collision_for_ions_hard_sphere(ll,cIx, bg_vel, com_vel)
                    case (PC_com_frame)
                        call scatter_velocity_for_ions(ll, bg_vel, cType, com_vel)

                    end select
                case (CS_ion_backwardscattering_t)
                !  print *, "before backwardscattering:", PC_ion_particles(ll)%v, &
                !          & vel_to_en(sqrt(sum(PC_ion_particles(ll)%v**2)), PC_particle_ion_mass)/ UC_elec_volt
                    select case (PC_frame_type)
                    case (PC_lab_frame)
                        call scatter_velocity_for_ions(ll, bg_vel,cType, com_vel)
                    case (PC_com_frame)
                        if (PC_BWflag == 1) then
                            call backwardscattering_turner(ll, bg_vel)
                        else if (PC_BWflag == 2) then
                            call scatter_velocity_for_ions(ll, bg_vel,cType, com_vel)
                        else
                            print *, "the scheme of backward scattering is unknow!"
                            stop
                        end if
                    end select
                ! print *, "after backwardscattering:", PC_ion_particles(ll)%v, &
                    !          & vel_to_en(sqrt(sum(PC_ion_particles(ll)%v**2)), PC_particle_ion_mass)/ UC_elec_volt
                end select
            end if
            
            coll_time = sample_coll_time_for_ions() ! Set the next collision time 
         end if
            
         goto 555    ! get out of the routine, to make sure ther is only one collision during one time step
      end do
      
      555 continue

      tempDt = PC_ion_particles(ll)%t_left
      ! Update the particle position and velocity to the next timestep
      if (PC_ion_particles(ll)%live) call advance_ion_particle(ll, PC_ion_particles(ll)%t_left)
      
      !> boundary conditions
      if (isOutOfGas(PC_ion_particles(ll)%x)) then
          call PC_boundary_for_ion(ll, PC_ion_particles(ll)%x(3), tempDt )
      end if
   end subroutine move_and_collide_for_ions

   !> Checks whether pos is outside the computational domain 
   !! and if so returns .TRUE., otherwise returns .FALSE.
   logical function isOutOfGas(pos)
      real(dp), dimension(3), intent(IN)     :: pos
     
      ! Check whether the component of the position is outside of the domain (1D)
      isOutOfGas = (pos(3) <= PC_posBD(1) .or. pos(3) >= PC_posBD(2))

   end function isOutOfGas

   !> Returns a sample from the exponential distribution of the collision times
   ! RNG_uniform() is uniform on [0,1), but log(0) = nan, so we take 1 - RNG_uniform()
   real(dp) function sample_coll_time()
      use m_random
      sample_coll_time = -log(1.0D0 - RNG_uniform()) * PC_coll%inv_max_rate
   end function sample_coll_time

   real(dp) function sample_coll_time_for_ions()
      use m_random
      sample_coll_time_for_ions = -log(1.0D0 - RNG_uniform()) * PC_ion_coll%inv_max_rate
   end function sample_coll_time_for_ions

   !> From the list crosssec(:) select the index of the process that will occur,
   !! or set colIndex = 0 if there is a null collision
   integer function get_coll_index(velocity)
      use m_random
      real(dp), intent(IN) :: velocity
      integer              :: j
      real(dp)             :: rate, rand_rate, sum_rate
      type(LT_loc_type)    :: lt_pos

      call LT_get_location(PC_coll%sum_rate_lt, velocity**2, lt_pos)
      get_coll_index = 0
      call LT_get(PC_coll%sum_rate_lt, lt_pos, 1, sum_rate) ! Get value at first (and only) column
      rand_rate      = RNG_uniform() * PC_coll%max_rate ! Random collision frequency

      if (rand_rate > sum_rate) return ! If we have a null collision we are done

      ! Otherwise, fill the coll_rateList array first with the appropriate values
      call LT_get(PC_coll%rate_lt, lt_pos, PC_coll%buff)

      !Anbang: here i test if we use two sets of random numbers
      rand_rate = RNG_uniform() * sum_rate

      ! Determine the type of collision
      ! This uses that process j has a probability coll_rate(j) / PC_coll%max_rate
      rate = 0.0D0
      do j = 1, PC_coll%num
         rate = rate + PC_coll%buff(j)
         if (rand_rate <= rate) then
            get_coll_index = j
            exit
         end if
      end do

   end function get_coll_index

   integer function get_coll_index_for_ions(velocity)
      use m_random
      real(dp), intent(IN) :: velocity
      integer              :: j
      real(dp)             :: rate, rand_rate, sum_rate
      type(LT_loc_type)    :: lt_pos

      
      call LT_get_location(PC_ion_coll%sum_rate_lt, velocity**2, lt_pos)
      get_coll_index_for_ions = 0
      call LT_get(PC_ion_coll%sum_rate_lt, lt_pos, 1, sum_rate) ! Get value at first (and only) column
      rand_rate      = RNG_uniform() * PC_ion_coll%max_rate ! Random collision frequency

      if (rand_rate > sum_rate) return ! If we have a null collision we are done

      ! Otherwise, fill the coll_rateList array first with the appropriate values
      call LT_get(PC_ion_coll%rate_lt, lt_pos, PC_ion_coll%buff)

      !Anbang: here i test if we use two sets of random numbers
      rand_rate = RNG_uniform() * sum_rate

      ! Determine the type of collision
      ! This uses that process j has a probability coll_rate(j) / PC_ion_coll%max_rate
      rate = 0.0D0
      do j = 1, PC_ion_coll%num
         rate = rate + PC_ion_coll%buff(j)
         if (rand_rate <= rate) then
            get_coll_index_for_ions = j
            exit
         end if
      end do

   end function get_coll_index_for_ions

   !Anbang: Here i test if we use formular of phelps, the results maybe better 
   integer function get_coll_index_for_ions_formular(en)
      use m_random
      use m_gas
      use m_units_constants

      real(dp), intent(IN) :: en
      integer              :: j
      real(dp)             :: rate, rand_rate, sum_rate
      type(LT_loc_type)    :: lt_pos
      real(dp)             :: cs_iso, cs_bw, rate_iso, rate_bw

      ! we should get sum rate here
      get_coll_index_for_ions_formular = 0
      if (en > 0.d0) then
        ! for argon
        ! from lxcat
!         cs_iso = 2.d-19 / (2.d0 * en)**(0.5d0) / (1.d0 + 2.d0 * en) + 3.d-19 * 2.d0 * en / (1.d0 + 2.d0 * en / 3.d0)**(2.3d0)
!         cs_bw  = 0.5d0 * (1.15d-18 * (1.d0 + 0.015/2.d0/en)**(0.6d0) / (2.d0 * en)**(0.1d0) - 2.d-19/(2.d0 * en)**(0.5d0) &
!                 & / (1.d0 + 2.d0 * en) - 3.d-19 * 2.d0 * en/ (1.d0 + 2.d0 * en/3.d0)**(2.3d0))

       ! from http://jilawww.colorado.edu/~avp/collision_data/ionneutral/IONATOM.TXT
!            cs_iso = 21.6d0/(en)**0.5d0 * 1.d-20
!            cs_bw =  52.d0/(en)**0.08d0 /(1.d0 + 0.08d0/en)/(1.d0 + en/1.d3)**0.3d0 * 1.d-20

         ! from the phelps 1994, J. Appl. Phys.   
         ! Remark: this is the same as used by Donko
         cs_iso = 2.d-19 / (en**0.5d0 * (1.d0 + en)) + 3.d-19 * en / (1.d0 + en/3.d0)**2.3d0
         cs_bw = (1.15d-18 * en**(-0.1d0) * (1.d0 + 0.015d0 /en)**0.6d0 - cs_iso) /2.d0
      else
        cs_iso = 0.d0
        cs_bw  = 0.d0
      end if

      rate_iso = cs_iso * GAS_get_number_dens() * en_to_vel(en*UC_elec_volt, PC_reduced_mass)
      rate_bw  = cs_bw  * GAS_get_number_dens() * en_to_vel(en*UC_elec_volt, PC_reduced_mass)
        
      sum_rate = rate_iso + rate_bw

      rand_rate      = RNG_uniform() * PC_ion_coll%max_rate ! Random collision frequency
     ! print *, "some parameters:",en, sum_rate, PC_ion_coll%max_rate, GAS_get_number_dens()

      if (rand_rate > sum_rate) return ! If we have a null collision we are done

!       ! Otherwise, fill the coll_rateList array first with the appropriate values
!       call LT_get(PC_ion_coll%rate_lt, lt_pos, PC_ion_coll%buff)
      !> 1: isotropic; 2: backwarscattering
      PC_ion_coll%buff(1) = rate_iso
      PC_ion_coll%buff(2) = rate_bw

      !Anbang: here i test if we use two sets of random numbers
      rand_rate = RNG_uniform() * sum_rate

      ! Determine the type of collision
      ! This uses that process j has a probability coll_rate(j) / PC_ion_coll%max_rate
      rate = 0.0D0
      do j = 1, PC_ion_coll%num
         rate = rate + PC_ion_coll%buff(j)
         if (rand_rate <= rate) then
            get_coll_index_for_ions_formular = j
            exit
         end if
      end do

   end function get_coll_index_for_ions_formular

   !> Perform an elastic collision for particle 'll'
    subroutine elastic_collision_COM(ll, coll_ix)
        integer, intent(IN)  :: ll, coll_ix
        real(dp)             :: bg_vel(3), com_vel(3)
        

!         if (PC_useBgVel) then
!             call creatANeutralVel(bg_vel)
!         else
!             bg_vel = 0.0_dp
!         end if

!        Anbang: we use cold gas here
        bg_vel = 0.0_dp

        ! Compute center of mass velocity
        com_vel = (PC_coll%spec_values(coll_ix) * PC_particles(ll)%v + bg_vel) / &
            (1 + PC_coll%spec_values(coll_ix))

        ! Scatter in center of mass coordinates
        PC_particles(ll)%v = PC_particles(ll)%v - com_vel
        call scatter_isotropic(ll, sqrt(sum(PC_particles(ll)%v**2)))
        PC_particles(ll)%v = PC_particles(ll)%v + com_vel
    end subroutine elastic_collision_COM

       !> Like Donko did, psst, 2011
       !> the netural gas is cold, only scatter electrons in the lab frame
    subroutine test_electron_elestic(ll)
        integer, intent(IN)  :: ll    
        call scatter_isotropic(ll, sqrt(sum(PC_particles(ll)%v**2)))
    end subroutine test_electron_elestic

 !> Anbang: we test if the rotations are important or not
    subroutine elastic_collision_COM_rotation(ll, coll_ix)
        use m_units_constants
        use m_random

        integer, intent(IN)  :: ll, coll_ix
        real(dp)             :: bg_vel(3), com_vel(3)
        real(dp)                :: chi, psi

        real(dp) ::  theta, phi, vel(3), temp, vel_norm
        real(dp) :: costheta, cosphi, sintheta, sinphi, coschi
        real(dp) :: cospsi, sinchi, sinpsi
        
!         if (PC_useBgVel) then
!             call creatANeutralVel(bg_vel)
!         else
!             bg_vel = 0.0_dp
!         end if

        !anbang: we set the bg as cold gas
        bg_vel = 0.0_dp

        ! Compute center of mass velocity
        com_vel = (PC_coll%spec_values(coll_ix) * PC_particles(ll)%v + bg_vel) / &
            (1 + PC_coll%spec_values(coll_ix))

      
      ! the velocity in the com frame
      vel_norm = sqrt(sum((PC_particles(ll)%v - bg_vel)**2))


      ! we assume isotropic scattering in the center of mass frame
      !Then, we have

       theta  = acos(1.0_dp - 2.0_dp * RNG_uniform())  ! isotropic scattering

       phi    = 2.d0 * UC_pi * RNG_uniform()

       ! energy loss, due to the elastic collision
       vel_norm  = sqrt(vel_norm **2 * ( 1.d0 - 2.d0 * PC_coll%spec_values(coll_ix) * (1.d0 - cos(theta))))

       ! Here i test another method developed by jannis
       vel = PC_particles(ll)%v - bg_vel

       select case(PC_rot_scheme)
        case(1,2)
            call UC_xyz_to_spherical(vel, temp, chi, psi)
        case(3)
            call UC_Donko_xyz_to_spherical(vel,temp,chi,psi)
       end select

      costheta = cos(theta)
      sintheta = sin(theta)
      cosphi   = cos(phi)
      sinphi   = sin(phi)
      coschi   = cos(chi)
      sinchi   = sin(chi)
      cospsi   = cos(psi)
      sinpsi   = sin(psi)

      select case(PC_rot_scheme)
      case(1)
      ! scattering of the relative velocity, See the paper of Yousfi, Physical review E, 49(4), 1994
        vel(1) = - sintheta * sinphi * sinpsi + sintheta * cosphi * coschi * cospsi + costheta * sinchi * cospsi
        vel(2) = sintheta * sinphi * cospsi + sintheta * cosphi * coschi * sinpsi + costheta * sinchi * sinpsi
        vel(3) = - sintheta * cosphi * sinchi + costheta * coschi
      case(2)
     ! from Deltlef
        vel(1) = - sintheta * sinpsi * cosphi - coschi * cospsi * sintheta * sinphi + sinchi * cospsi * costheta
        vel(2) = sintheta * cospsi * cosphi - coschi * sinpsi * sintheta * sinphi + sinchi * sinpsi * costheta
        vel(3) = sinchi * sintheta * sinphi + coschi * costheta
      case(3) 
        ! from Donko, theta is the angle between Vel and x axis
        vel(1) = coschi * costheta - sinchi * sintheta * cosphi
        vel(2) = sinchi * cospsi * costheta + coschi * cospsi * sintheta * cosphi - sinpsi * sintheta * sinphi
        vel(3) = sinchi * sinpsi * costheta + coschi * sinpsi * sintheta * cosphi + cospsi * sintheta * sinphi
      end select


      vel    = vel * vel_norm

      ! now we need to put the velocity back to the laboratory frame
      PC_particles(ll)%v = PC_particle_atom_mass / (PC_particle_elec_mass + PC_particle_atom_mass) * vel + &
                            com_vel

    end subroutine elastic_collision_COM_rotation

 
   !> Perform an excitation-collision for particle 'll'
   subroutine excite_collision(ll, coll_ix, velocity)
      use m_units_constants
      integer, intent(IN)  :: ll, coll_ix
      real(dp), intent(IN) :: velocity
      real(dp)             :: new_vel, energy, old_en

      old_en  = vel_to_en(velocity, PC_particle_elec_mass)
      energy  = max(0.0_dp, old_en - PC_coll%spec_values(coll_ix))
      new_vel = en_to_vel(energy, PC_particle_elec_mass)

      call scatter_isotropic(ll, new_vel)
   end subroutine excite_collision

   !> Perform an ionizing collision for particle 'll'
   subroutine ionization_collision(ll, coll_ix, velocity)
      use m_units_constants
      integer, intent(IN)  :: ll, coll_ix
      real(dp), intent(IN) :: velocity
      real(dp)             :: energy, old_en, old_vel, new_vel
      real(dp)             :: en_s1, en_s2
      integer              :: idx

      PC_num_part = PC_num_part + 1
      call check_num_particles(PC_num_part)

      old_en = vel_to_en(velocity, PC_particle_elec_mass)
      energy = max(0.0_dp, old_en - PC_coll%spec_values(coll_ix))
      en_s1  = 0.5D0 * energy
      en_s2  = 0.5D0 * energy

      PC_particles(PC_num_part) = PC_particles(ll)
      old_vel = en_to_vel(en_s1, PC_particle_elec_mass)
      new_vel = en_to_vel(en_s2, PC_particle_elec_mass)

      call scatter_isotropic(ll, old_vel)
      call scatter_isotropic(PC_num_part, new_vel)
      
      ! a new ion is also added to the list
      PC_num_ion_part = PC_num_ion_part + 1
      call check_num_ion_particles(PC_num_ion_part)
      
      ! how to give a new energy to ions? now we simplily give zero velocity 
      PC_ion_particles(PC_num_ion_part)%x = PC_particles(ll)%x
      PC_ion_particles(PC_num_ion_part)%weight = PC_particles(ll)%weight
      PC_ion_particles(PC_num_ion_part)%T_left = PC_particles(ll)%T_left
      PC_ion_particles(PC_num_ion_part)%a = 0.d0 ! approximately
      PC_ion_particles(PC_num_ion_part)%live = PC_particles(ll)%live

      !PC_ion_particles(PC_num_ion_part)%v = 0.d0         ! we make it's velocity as zero
      ! make it a have a energy of 0.04eV
      call scatter_isotropic_for_ions(PC_num_ion_part, en_to_vel(0.04d0 * UC_elec_Volt, PC_particle_ion_mass))
    ! we can set the new ion velocity as the velocity of the neutral molecule
    !  call creatANeutralVel(PC_ion_particles(PC_num_ion_part)%v)

      ! calculate the ionization source term 
      idx = floor(PC_ion_particles(PC_num_ion_part)%x(3)/ PC_delta_x) + 1
      if (idx < 1) idx = 1
      if (idx > PC_grid_size) idx = PC_grid_size
      PC_ionization_number(idx) = PC_ionization_number(idx) + PC_ion_particles(PC_num_ion_part)%weight
      
   end subroutine ionization_collision

   subroutine PC_output_ionization_number(number)

        real(dp), intent(out) :: number(:)
        number = PC_ionization_number
   end subroutine PC_output_ionization_number

   !> Perform attachment of electron 'll'
   subroutine attach_collision(ll)
      integer, intent(IN) :: ll
      call PC_killElec(ll)
   end subroutine attach_collision

   ! Anbang: Here we creat a neutral molecule, for modeing the collision with ions
   subroutine creatANeutralVel(bg_vel)
      use m_units_constants
      use m_random
      real(dp)             :: vel_norm, temp
      real(dp),intent(out) :: bg_vel(3)
      integer              :: nn

       ! Set a Maxwellian velocity distribution with the correct mean energy
       
      !Method 1: like what jannis did
      !Anbang: devide the total velocity into three dimensions
!         vel_norm = 3.d0 / 2.d0 * UC_boltzmann_const * PC_gas_tem
!         vel_norm = sqrt(2 * vel_norm / (3 * PC_particle_atom_mass))    
!         bg_vel = (/RNG_normal(), RNG_normal(), RNG_normal()/)
!         bg_vel = bg_vel * vel_norm

       
       ! method 2: from a master thesis
       !> Attention: Here if random number equals zero, the log(0.d0) is not right, becareful when use it
        do nn = 1, 3
            temp = RNG_uniform()
            if (temp == 0.d0) temp = temp +UC_tiny
            bg_vel(nn) = sqrt(-log(temp) * 2.d0 * UC_boltzmann_const * PC_gas_tem / PC_particle_atom_mass ) * &
                            & sin(2.d0 * UC_pi * RNG_uniform())
        end do

   end subroutine creatANeutralVel

   subroutine scatter_isotropic(ll, vel_norm)
      use m_units_constants
      use m_random
      integer, intent(in)  :: ll
      real(dp), intent(in) :: vel_norm
      real(dp)             :: theta, phi

      theta             = acos(1.0_dp - 2.0_dp * RNG_uniform())
      phi               = 2 * UC_pi * RNG_uniform()
      PC_particles(ll)%v(1) = vel_norm * sin(theta) * cos(phi)
      PC_particles(ll)%v(2) = vel_norm * sin(theta) * sin(phi)
      PC_particles(ll)%v(3) = vel_norm * cos(theta)
   end subroutine scatter_isotropic


!!!!!!!******************!!!!!!!!!!!!!!
!!!!!! For ions !!!!!!!!!!!!!!!!!
!!!!!!!******************!!!!!!!!!!!!!!
    !> Ions: Perform an elastic collision for particle 'll'
    subroutine elastic_collision_for_ions_COM(ll, com_vel)
      integer, intent(IN)  :: ll
      real(dp), intent(in) :: com_vel(3)

     ! Scatter in center of mass coordinates
        PC_ion_particles(ll)%v = PC_ion_particles(ll)%v - com_vel
        call scatter_isotropic_for_ions(ll, sqrt(sum(PC_ion_particles(ll)%v**2)))
        PC_ion_particles(ll)%v = PC_ion_particles(ll)%v + com_vel

    end subroutine elastic_collision_for_ions_COM

    !>Anbang: The hard sphere collisions See. Vahedi, Surendra. 87 (1995) 179-198, Comput. phys. Commu.
    subroutine elastic_collision_for_ions_hard_sphere(ll, coll_ix, bg_vel, com_vel)
      use m_random
      use m_units_constants

      integer, intent(IN)  :: ll, coll_ix
      real(dp)             :: theta, factorLoss
      real(dp)             :: en_new, en_old, vel_norm, phi
      real(dp),intent(in)  :: bg_vel(3), com_vel(3) 
      
      vel_norm = sqrt(sum((PC_ion_particles(ll)%v - bg_vel)**2))
      en_old = vel_to_en(vel_norm, PC_ion_coll%spec_values(coll_ix) * PC_particle_atom_mass)
      
      theta             = acos(1.0_dp - 2.0_dp * RNG_uniform())

      factorLoss = 2.d0 * PC_ion_coll%spec_values(coll_ix) * PC_particle_atom_mass **2 / &
                &(PC_ion_coll%spec_values(coll_ix) *PC_particle_atom_mass + PC_particle_atom_mass) **2 &
                & * (1.d0 - cos(theta)) 

      en_new = (1.d0 - factorLoss) * en_old
      
      vel_norm = en_to_vel(en_new, PC_ion_coll%spec_values(coll_ix) * PC_particle_atom_mass)

      phi               = 2 * UC_pi * RNG_uniform()
      PC_ion_particles(ll)%v(1) = vel_norm * sin(theta) * cos(phi)
      PC_ion_particles(ll)%v(2) = vel_norm * sin(theta) * sin(phi)
      PC_ion_particles(ll)%v(3) = vel_norm * cos(theta)


    end subroutine elastic_collision_for_ions_hard_sphere


   subroutine scatter_isotropic_for_ions(ll, vel_norm)
      use m_units_constants
      use m_random
      integer, intent(in)  :: ll
      real(dp), intent(in) :: vel_norm
      real(dp)             :: theta, phi

      theta             = acos(1.0_dp - 2.0_dp * RNG_uniform())
      phi               = 2 * UC_pi * RNG_uniform()
      PC_ion_particles(ll)%v(1) = vel_norm * sin(theta) * cos(phi)
      PC_ion_particles(ll)%v(2) = vel_norm * sin(theta) * sin(phi)
      PC_ion_particles(ll)%v(3) = vel_norm * cos(theta)
   end subroutine scatter_isotropic_for_ions


  !Anbang: test tuner's advice, just change the velocity between ion and neutral
  !> this is like charge exchange collision
   subroutine backwardscattering_turner(ll, bg_vel)
      integer, intent(in)  :: ll
      real(dp), intent(in) :: bg_vel(3)
    
      PC_ion_particles(ll)%v = bg_vel
   end subroutine backwardscattering_turner

 
   ! Anbang: To get the scattering angles in the spherical coordinate
   ! chi: polar angle- the angle between the z axis and the velocity
   ! psi: azimuthal angle- the angle tetween the velocity component on the xy plane and x axis
   ! those angles are in the lab frame
   subroutine scatter_angle_for_ions(ll,chi, psi, bg_vel)
        use m_units_constants

        real(dp), intent(out) :: chi, psi
        real(dp), intent(in)  :: bg_vel(3)
        real(dp)              :: vel_norm
        integer,  intent(in)  :: ll
        
        ! the relative velocity, for the elastic collision.
        ! if it is an inelastic collision, we should minus some enery loss
        ! See the paper of Yousfi, Physical review E, 49(4), 1994
        vel_norm = sqrt(sum((PC_ion_particles(ll)%v - bg_vel)**2))

        if (vel_norm == 0.d0) then  ! here we set if the velocity is very small
            chi = 0.d0
            psi = 0.d0
        else
            chi = acos((PC_ion_particles(ll)%v(3) - bg_vel(3)) / vel_norm)

            if ((PC_ion_particles(ll)%v(1) - bg_vel(1)) == 0.d0) then
                if ((PC_ion_particles(ll)%v(2) - bg_vel(2)) == 0.d0) then
                    psi = 0.d0
                else
                    if (abs(vel_norm * sin(chi)) < abs(PC_ion_particles(ll)%v(2) - bg_vel(2))) then
                        psi = sign(asin(1.d0), (PC_ion_particles(ll)%v(2) - bg_vel(2)))   ! 90 degree
                    else
                        psi = asin((PC_ion_particles(ll)%v(2) - bg_vel(2)) / (vel_norm * sin(chi) + UC_tiny))
                    end if
                end if
            else
                psi = atan2((PC_ion_particles(ll)%v(2) - bg_vel(2)), (PC_ion_particles(ll)%v(1) - bg_vel(1)))
            end if
        end if

   end subroutine scatter_angle_for_ions

   subroutine scatter_velocity_for_ions(ll,bg_vel,cType, com_vel)
      use m_units_constants
      use m_random
      use m_cross_sec_for_ions

      real(dp)                :: chi, psi
      real(dp), intent(in)    :: bg_vel(3), com_vel(3)
      integer, intent(in)     :: cType, ll

      real(dp) :: vel_norm, theta, phi, vel(3), temp
      real(dp) :: costheta, cosphi, sintheta, sinphi, coschi
      real(dp) :: cospsi, sinchi, sinpsi
      
      ! the velocity in the center of mass frame
      vel_norm = sqrt(sum((PC_ion_particles(ll)%v - bg_vel)**2))

      ! we assume isotropic scattering / backward scattering in the center of mass frame
      !Then
       select case(cType)
       case(CS_ion_elastic_t)
            theta  = acos(1.0_dp - 2.0_dp * RNG_uniform())  ! isotropic scattering
       case(CS_ion_backwardscattering_t) 
            theta  = acos(-1.d0)        ! backward scattering
       end select
       phi    = 2.d0 * UC_pi * RNG_uniform()


      vel = PC_ion_particles(ll)%v - bg_vel

      select case(PC_rot_scheme)
      case(1,2)
            call scatter_angle_for_ions(ll,chi, psi, bg_vel)
      case(3)
            call UC_Donko_xyz_to_spherical(vel, temp, chi, psi)
      end select

      costheta = cos(theta)
      sintheta = sin(theta)
      cosphi   = cos(phi)
      sinphi   = sin(phi)
      coschi   = cos(chi)
      sinchi   = sin(chi)
      cospsi   = cos(psi)
      sinpsi   = sin(psi)


      select case(PC_rot_scheme)
      case(1)
            ! scattering of the relative velocity, See the paper of Yousfi, Physical review E, 49(4), 1994
            vel(1) = - sintheta * sinphi * sinpsi + sintheta * cosphi * coschi * cospsi + costheta * sinchi * cospsi
            vel(2) = sintheta * sinphi * cospsi + sintheta * cosphi * coschi * sinpsi + costheta * sinchi * sinpsi
            vel(3) = - sintheta * cosphi * sinchi + costheta * coschi
      case(2)
            ! from Deltlef
            vel(1) = - sintheta * sinpsi * cosphi - coschi * cospsi * sintheta * sinphi + sinchi * cospsi * costheta
            vel(2) = sintheta * cospsi * cosphi - coschi * sinpsi * sintheta * sinphi + sinchi * sinpsi * costheta
            vel(3) = sinchi * sintheta * sinphi + coschi * costheta
       case(3)
            ! from Donko, theta is the angle between Vel and x axis
            vel(1) = coschi * costheta - sinchi * sintheta * cosphi
            vel(2) = sinchi * cospsi * costheta + coschi * cospsi * sintheta * cosphi - sinpsi * sintheta * sinphi
            vel(3) = sinchi * sinpsi * costheta + coschi * sinpsi * sintheta * cosphi + cospsi * sintheta * sinphi
      end select

      vel    = vel * vel_norm

      ! now we need to put the velocity back to the laboratory frame
      PC_ion_particles(ll)%v = PC_particle_atom_mass / (PC_particle_ion_mass + PC_particle_atom_mass) * vel + &
                            com_vel
        
   end subroutine scatter_velocity_for_ions

   !> Advance the particle position and velocity over time tt, electrons
   subroutine advance_particle(ll, tt)
      integer, intent(IN)  :: ll
      real(dp), intent(IN) :: tt

      PC_particles(ll)%x      = PC_particles(ll)%x + PC_particles(ll)%v * tt + &
           0.5_dp * PC_particles(ll)%a * tt**2
      PC_particles(ll)%v      = PC_particles(ll)%v + PC_particles(ll)%a * tt
      PC_particles(ll)%t_left = PC_particles(ll)%t_left - tt
   end subroutine advance_particle

   !> Advance the particle position and velocity over time tt, ions
   subroutine advance_ion_particle(ll, tt)
      integer, intent(IN)  :: ll
      real(dp), intent(IN) :: tt

      PC_ion_particles(ll)%x      = PC_ion_particles(ll)%x + PC_ion_particles(ll)%v * tt + &
           0.5_dp * PC_ion_particles(ll)%a * tt**2
      PC_ion_particles(ll)%v      = PC_ion_particles(ll)%v + PC_ion_particles(ll)%a * tt
      PC_ion_particles(ll)%t_left = PC_ion_particles(ll)%t_left - tt
   end subroutine advance_ion_particle

   !> advance_elec_and_ion's position and velocity with verlet scheme
   subroutine PC_advance_elec_pos_verlet(tt)
      integer              :: ll
      real(dp), intent(IN) :: tt

      do ll = 1, PC_num_part
        PC_particles(ll)%x      = PC_particles(ll)%x + PC_particles(ll)%v * tt + &
            0.5_dp * PC_particles(ll)%a * tt**2 
        if (isOutOfGas(PC_particles(ll)%x)) then
            call PC_boundary_for_electron(ll, PC_particles(ll)%x(3), tt)
        end if      
      end do
      call remove_dead_particles()

      ! put the accel to old accel
      do ll = 1, PC_num_part
            PC_elec_old_accel(:,ll) = PC_particles(ll)%a(:)
      end do  
   end subroutine PC_advance_elec_pos_verlet

   subroutine PC_advance_ion_pos_verlet(tt)
      integer              :: ll
      real(dp), intent(IN) :: tt

      do ll = 1, PC_num_ion_part
        PC_ion_particles(ll)%x      = PC_ion_particles(ll)%x + PC_ion_particles(ll)%v * tt + &
            0.5_dp * PC_ion_particles(ll)%a * tt**2 
        if (isOutOfGas(PC_ion_particles(ll)%x)) then
            call PC_boundary_for_ion(ll, PC_ion_particles(ll)%x(3), tt)
        end if       
      end do
      call remove_dead_ion_particles()

      ! put the accel to old accel
      do ll = 1, PC_num_ion_part
            PC_ion_old_accel(:,ll) = PC_ion_particles(ll)%a(:)
      end do  
   end subroutine PC_advance_ion_pos_verlet

  ! calculate the velocity of particles, with Verlet scheme
  subroutine  PC_advance_elec_vel_verlet(tt)
      integer              :: ll
      real(dp), intent(IN) :: tt
      do ll = 1, PC_num_part
        PC_particles(ll)%v(3)      = PC_particles(ll)%v(3) + 0.5d0 * (PC_elec_old_accel(3,ll) + PC_particles(ll)%a(3)) * tt
      end do
  end subroutine PC_advance_elec_vel_verlet

  subroutine  PC_advance_ion_vel_verlet(tt)
      integer              :: ll
      real(dp), intent(IN) :: tt
      do ll = 1, PC_num_ion_part
        PC_ion_particles(ll)%v(3)      = PC_ion_particles(ll)%v(3) + &
                & 0.5d0 * (PC_ion_old_accel(3,ll) + PC_ion_particles(ll)%a(3)) * tt
      end do
  end subroutine PC_advance_ion_vel_verlet

! do collisions of particles seprately
   subroutine PC_elec_collide(tt)
      use m_cross_sec
      use m_utils

      integer              :: ll
      real(dp), intent(in) :: tt
      integer              :: cIx, cType
      real(dp)             :: coll_time, new_vel
      integer              :: elec_num

      PC_ionization_number = 0.d0

     ! we only do collision for orginal particles, new born particles are not done
      elec_num  = PC_num_part
      do ll = 1, elec_num
            coll_time = sample_coll_time()

            if (coll_time < tt)  then   
 
                new_vel = sqrt(sum(PC_particles(ll)%v**2))
                cIx     = get_coll_index(new_vel)

                if (cIx > 0) then
                    cType               = PC_coll%types(cIx)
                    PC_coll%counts(cIx) = PC_coll%counts(cIx) + PC_particles(ll)%weight

                    ! Perform the corresponding collision
                    select case (cType)
                    case (CS_attach_t)
                        call attach_collision(ll)
                    case (CS_elastic_t)
                        call elastic_collision_COM(ll, cIx)
                    case (CS_excite_t)
                        call excite_collision(ll, cIx, new_vel)
                    case (CS_ionize_t)
                        call ionization_collision(ll, cIx, new_vel)
                    end select
                end if
            end if
       end do     

       PC_ionization_number = PC_ionization_number /tt
    end subroutine PC_elec_collide

   subroutine PC_ion_collide(tt)
      use m_cross_sec_for_ions
      use m_utils
      use m_units_constants

      integer              :: ll
      real(dp), intent(in) :: tt
      integer              :: cIx, cType
      double precision     :: coll_time, tempDt
      real(dp)             :: bg_vel(3),com_vel(3), rel_vel
      integer              :: ion_num

      ion_num  = PC_num_ion_part

      do ll = 1, ion_num
            coll_time = sample_coll_time_for_ions()

            if (coll_time < tt) then
                    !Anbang: Here is not clear how turner set the background gas velocities
                    if (PC_useBgVel) then
                        call creatANeutralVel(bg_vel)
                    else
                        bg_vel = 0.d0
                    end if

                    ! Anbang: Here we set the COM velocity
                    com_vel = (PC_particle_ion_mass * PC_ion_particles(ll)%v + PC_particle_atom_mass * bg_vel) / &
                            (PC_particle_ion_mass + PC_particle_atom_mass)
                    

                    !relative velocity
                    rel_vel = sqrt(sum((PC_ion_particles(ll)%v - bg_vel)**2))
                
                    select case(PC_flagTable)
                    case(flag_vel)
                        cIx     = get_coll_index_for_ions(rel_vel)
                    case(flag_en)
                        ! velocity to com energy, to the look up table
                        rel_vel = vel_to_en( rel_vel, PC_reduced_mass) / UC_elec_volt
                        cIx    = get_coll_index_for_ions_formular(rel_vel)
                    end select    
                        
                    if (cIx > 0) then
                        cType               = PC_ion_coll%types(cIx)
                        PC_ion_coll%counts(cIx) = PC_ion_coll%counts(cIx) + PC_ion_particles(ll)%weight

                        ! Perform the corresponding collision
                        select case (cType)
                        case (CS_ion_elastic_t)
                            select case (PC_frame_type)
                            case (PC_lab_frame)
                                call elastic_collision_for_ions_hard_sphere(ll,cIx, bg_vel, com_vel)
                            case (PC_com_frame)
                                call scatter_velocity_for_ions(ll, bg_vel, cType, com_vel)

                            end select
                        case (CS_ion_backwardscattering_t)
                        !  print *, "before backwardscattering:", PC_ion_particles(ll)%v, &
                        !          & vel_to_en(sqrt(sum(PC_ion_particles(ll)%v**2)), PC_particle_ion_mass)/ UC_elec_volt
                            select case (PC_frame_type)
                            case (PC_lab_frame)
                                call scatter_velocity_for_ions(ll, bg_vel,cType, com_vel)
                            case (PC_com_frame)
                                if (PC_BWflag == 1) then
                                    call backwardscattering_turner(ll, bg_vel)
                                else if (PC_BWflag == 2) then
                                    call scatter_velocity_for_ions(ll, bg_vel,cType, com_vel)
                                else
                                    print *, "the scheme of backward scattering is unknow!"
                                    stop
                                end if
                            end select
                        ! print *, "after backwardscattering:", PC_ion_particles(ll)%v, &
                            !          & vel_to_en(sqrt(sum(PC_ion_particles(ll)%v**2)), PC_particle_ion_mass)/ UC_elec_volt
                        end select
                    end if 
            end if

        end do
    end subroutine PC_ion_collide
   
  

   subroutine PC_set_accel(accel_func)
      procedure(if_ipart_ovec3) :: accel_func
      integer                   :: ll
      real(dp)                  :: new_accel(3)

      do ll = 1, PC_num_part
         call accel_func(PC_particles(ll), new_accel)
         PC_particles(ll)%a = new_accel
      end do
   end subroutine PC_set_accel

   subroutine PC_set_ion_accel(accel_func)
      procedure(if_ion_ipart_ovec3) :: accel_func
      integer                   :: ll
      real(dp)                  :: new_accel(3)

      do ll = 1, PC_num_ion_part
         call accel_func(PC_ion_particles(ll), new_accel)
         PC_ion_particles(ll)%a = new_accel
      end do
   end subroutine PC_set_ion_accel


   !> Correct particle velocities for the previous timestep of 'dt'
   !!
   !! During the timestep x,v have been advanced to:
   !! x(t+1) = x(t) + v(t)*dt + 0.5*a(t)*dt^2,
   !! v(t+1) = v(t) + a(t)*dt
   !! But the velocity at t+1 should be v(t+1) = v(t) + 0.5*(a(t) + a(t+1))*dt,
   !! to have a second order leapfrog scheme, so here we set it to that value.
   subroutine PC_correct_new_accel(dt, accel_func)
      use m_units_constants
      real(dp), intent(IN)      :: dt
      procedure(if_ipart_ovec3) :: accel_func
      integer                   :: ll
      real(dp)                  :: new_accel(3)

      do ll = 1, PC_num_part
         call accel_func(PC_particles(ll), new_accel)
         PC_particles(ll)%v = PC_particles(ll)%v + 0.5_dp * (new_accel - PC_particles(ll)%a) * dt
         PC_particles(ll)%a = new_accel
      end do
      !print *, "the maximum/ minmum accel of electrons is:", maxval(PC_particles(:)%a(3)), minval(PC_particles(:)%a(3))
   end subroutine PC_correct_new_accel

!   !> Anbang: same for ions
   subroutine PC_correct_new_ion_accel(dt, accel_func)
      use m_units_constants
      real(dp), intent(IN)          :: dt
      procedure(if_ion_ipart_ovec3) :: accel_func
      integer                       :: ll
      real(dp)                      :: new_accel(3)

      do ll = 1, PC_num_ion_part
         call accel_func(PC_ion_particles(ll), new_accel)
         PC_ion_particles(ll)%v = PC_ion_particles(ll)%v + 0.5_dp * (new_accel - PC_ion_particles(ll)%a) * dt
         PC_ion_particles(ll)%a = new_accel
      end do
   end subroutine PC_correct_new_ion_accel

   !> Remove all particles for which live == .FALSE. from the list, which can happen
   !! due moving outside the computational domain or attachment.
   subroutine remove_dead_particles()
      integer :: ll, ix_end
      ix_end = PC_num_part

      do ll = 1, PC_num_part
         if (.not. PC_particles(ll)%live) then
            ! Find the first alive particle from the end of the list, and place it at index ll
            do while (.not. PC_particles(ix_end)%live .and. ix_end > ll)
               ix_end = ix_end - 1
            end do

            ! If there is no alive particle available in the range [ll+1 : PC_num_part] we have
            ! ll >= ix_end and we exit the routine. Else we put the alive particle at index ll
            if (ll >= ix_end) then
               PC_num_part = ll - 1
               exit
            else
               PC_particles(ll)          = PC_particles(ix_end)
               PC_particles(ix_end)%live = .false.
               ix_end                = ix_end - 1
            end if
         end if
      end do

   end subroutine remove_dead_particles


   subroutine PC_create_part(pos, vel, accel, weight)
      real(dp), intent(IN) :: pos(3), vel(3), accel(3)
      real(dp), intent(IN)  :: weight

      PC_num_part = PC_num_part + 1
      call check_num_particles(PC_num_part)

      PC_particles(PC_num_part)%weight   = weight
      PC_particles(PC_num_part)%live     = .true.
      PC_particles(PC_num_part)%t_left   = 0.0_dp
      PC_particles(PC_num_part)%x        = pos
      PC_particles(PC_num_part)%v        = vel
      PC_particles(PC_num_part)%a    = accel
   end subroutine PC_create_part

   !> Mark particle 'll' as inactive, so that it will be removed from the particle list
   subroutine PC_killElec(ll)
      integer, intent(IN) :: ll
      PC_particles(ll)%live = .false.
   end subroutine PC_killElec

   subroutine PC_periodify(is_periodic, lengths)
      logical, intent(in) :: is_periodic(3)
      real(dp), intent(in) :: lengths(3)
      integer :: ix

      do ix = 1, 3
         if (is_periodic(ix)) then
            PC_particles(1:PC_num_part)%x(ix) = modulo(PC_particles(1:PC_num_part)%x(ix), lengths(ix))
         end if
      end do
   end subroutine PC_periodify

   subroutine PC_get_part(ll, my_part)
      integer, intent(in)          :: ll
      type(PC_part_t), intent(out) :: my_part
      my_part = PC_particles(ll)
   end subroutine PC_get_part

   subroutine PC_set_part(ll, my_part)
      type(PC_part_t), intent(in) :: my_part
      integer, intent(in)         :: ll
      PC_particles(ll) = my_part
   end subroutine PC_set_part

   !> Return the number of real particles
   real(dp) function PC_get_num_real_part()
      PC_get_num_real_part = sum(PC_particles(1:PC_num_part)%weight)
   end function PC_get_num_real_part

   !> Return the number of simulation particles
   integer function PC_get_num_sim_part()
      PC_get_num_sim_part = PC_num_part
   end function PC_get_num_sim_part

   !> Loop over all the particles and call pptr_loop for each of them
   subroutine PC_loop_part(pptr_loop)
      procedure(if_ipart) :: pptr_loop
      integer             :: n
      do n = 1, PC_num_part
         call pptr_loop(PC_particles(n))
      end do
   end subroutine PC_loop_part

   subroutine check_num_particles(nPart)
      use m_error
      integer, intent(in) :: nPart
      if (nPart > PC_max_num_part) call ERR_show("Particle module error (electrons), too many particles")
   end subroutine check_num_particles

   subroutine check_num_ion_particles(nPart)
      use m_error
      integer, intent(in) :: nPart
      if (nPart > PC_max_num_ion_part) then
        print *, "npart, PC_max_num_ion_part", nPart, PC_max_num_ion_part
        call ERR_show("Particle module error (ions), too many particles")
      end if
   end subroutine check_num_ion_particles

   real(dp) function vel_to_en(vel, mass)
      real(dp), intent(in) :: vel, mass
      vel_to_en = 0.5_dp * mass * vel**2
   end function vel_to_en

   real(dp) function en_to_vel(en, mass)
      real(dp), intent(in) :: en, mass
      en_to_vel = sqrt(2 * en / mass)
   end function en_to_vel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!1For ions !!!!!!!!!!!!!!!!
   !> Remove all particles for which live == .FALSE. from the list, which can happen
   !! due moving outside the computational domain or attachment.
   subroutine remove_dead_ion_particles()
      integer :: ll, ix_end
      ix_end = PC_num_ion_part

      do ll = 1, PC_num_ion_part
         if (.not. PC_ion_particles(ll)%live) then
            ! Find the first alive particle from the end of the list, and place it at index ll
            do while (.not. PC_ion_particles(ix_end)%live .and. ix_end > ll)
               ix_end = ix_end - 1
            end do

            ! If there is no alive particle available in the range [ll+1 : PC_num_ion_part] we have
            ! ll >= ix_end and we exit the routine. Else we put the alive particle at index ll
            if (ll >= ix_end) then
               PC_num_ion_part = ll - 1
               exit
            else
               PC_ion_particles(ll)          = PC_ion_particles(ix_end)
               PC_ion_particles(ix_end)%live = .false.
               ix_end                = ix_end - 1
            end if
         end if
      end do

   end subroutine remove_dead_ion_particles

   subroutine PC_create_ion_part(pos, vel, accel, weight)
      real(dp), intent(IN) :: pos(3), vel(3), accel(3)
      real(dp), intent(IN)  :: weight

      PC_num_ion_part = PC_num_ion_part + 1
      call check_num_ion_particles(PC_num_ion_part)

      PC_ion_particles(PC_num_ion_part)%weight   = weight
      PC_ion_particles(PC_num_ion_part)%live     = .true.
      PC_ion_particles(PC_num_ion_part)%t_left   = 0.0_dp
      PC_ion_particles(PC_num_ion_part)%x        = pos
      PC_ion_particles(PC_num_ion_part)%v        = vel
      PC_ion_particles(PC_num_ion_part)%a    = accel
   end subroutine PC_create_ion_part

   !> Mark particle 'll' as inactive, so that it will be removed from the particle list
   subroutine PC_killIon(ll)
      integer, intent(IN) :: ll
      PC_ion_particles(ll)%live = .false.
   end subroutine PC_killIon

   subroutine PC_ion_periodify(is_periodic, lengths)
      logical, intent(in) :: is_periodic(3)
      real(dp), intent(in) :: lengths(3)
      integer :: ix

      do ix = 1, 3
         if (is_periodic(ix)) then
            PC_ion_particles(1:PC_num_ion_part)%x(ix) = modulo(PC_ion_particles(1:PC_num_ion_part)%x(ix), lengths(ix))
         end if
      end do
   end subroutine PC_ion_periodify


   !> Return the number of real particles
   real(dp) function PC_get_num_real_ion_part()
      PC_get_num_real_ion_part = sum(PC_ion_particles(1:PC_num_ion_part)%weight)
   end function PC_get_num_real_ion_part

   !> Return the number of simulation particles
   integer function PC_get_num_sim_ion_part()
      PC_get_num_sim_ion_part = PC_num_ion_part
   end function PC_get_num_sim_ion_part



  !> Loop over all the particles and call pptr_loop for each of them
   subroutine PC_loop_ion_part(pptr_loop)
      procedure(if_ion_ipart) :: pptr_loop
      integer                 :: n
      do n = 1, PC_num_ion_part
         call pptr_loop(PC_ion_particles(n))
      end do
   end subroutine PC_loop_ion_part

   !!!!!!!!!!!!End for ions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   !> Create a lookup table with cross sections for a number of energies
   subroutine create_coll_rate_table(cross_secs, min_eV, max_eV, table_size)
      use m_units_constants
      use m_cross_sec
      use m_utils
      use m_gas
      type(CS_type), intent(in) :: cross_secs(:)
      integer, intent(IN)       :: table_size
      real(dp), intent(in)      :: min_eV, max_eV

      real(dp)                  :: vel_list(table_size), rate_list(table_size), sum_rate_list(table_size)
      real(dp)                  :: vel2_list(table_size)
      integer                   :: i, n, m, iXc
      integer, allocatable      :: collIxs(:)
      real(dp)                  :: min_vel, max_vel, en_eV, temp, min_vel2, max_vel2

      character(len=name_len) :: filename
      integer              :: myunit
      real(dp), allocatable :: tempRate(:,:)

      min_vel = en_to_vel(min_eV * UC_elec_volt, PC_particle_elec_mass)
      max_vel = en_to_vel(max_eV * UC_elec_volt, PC_particle_elec_mass)
        
      max_vel2 = max_vel**2.d0
      min_vel2 = min_vel**2.d0

      PC_coll%num = count(cross_secs(:)%min_energy < max_eV .and. cross_secs(:)%max_energy > min_eV)
      allocate( collIxs(PC_coll%num) )
      allocate( PC_coll%types(PC_coll%num) )
      allocate( PC_coll%spec_values(PC_coll%num) )
      allocate( PC_coll%buff(PC_coll%num) )
      allocate( tempRate(PC_coll%num, table_size))
      tempRate = 0.d0
        
      ! Get the indexes of the collisions in this energy range
      m = 1
      do n = 1, size(cross_secs)
         if (cross_secs(n)%min_energy < max_eV .and. cross_secs(n)%max_energy > min_eV) then
            collIxs(m)  = n
            m           = m + 1
         end if
      end do

      ! Set the velocities for the lookup table
      do i = 1, table_size
         temp = real(i-1, dp) / (table_size-1)
         vel_list(i) = min_vel + temp * (max_vel - min_vel)
         vel2_list(i) = vel_list(i)**2.d0
      end do

      ! Create collision rate table
      PC_coll%rate_lt => LT_create(min_vel2, max_vel2, table_size, n_storage = PC_coll%num)

      sum_rate_list = 0.0_dp
        
!       print *, "PC_coll%num =", PC_coll%num

      do n = 1, PC_coll%num
         iXc = collIxs(n)
         PC_coll%types(n) = cross_secs(iXc)%col_type

         ! Ionizations and excitations use a threshold energy which is lost
         select case (PC_coll%types(n))
         case (CS_ionize_t, CS_excite_t)
            PC_coll%spec_values(n) = cross_secs(iXc)%spec_value * UC_elec_volt
!             print *, "the spec_value of ionization/ excitation of electrons:", PC_coll%spec_values(n)/UC_elec_volt
         case (CS_elastic_t)
            PC_coll%spec_values(n) = cross_secs(iXc)%spec_value
!             print *, "the spec_value of eleastic coll of electrons:", PC_coll%spec_values(n)
         case DEFAULT
            PC_coll%spec_values(n) = 0.0_dp
         end select

         do i = 1, table_size
            en_eV = vel_to_en(vel_list(i), PC_particle_elec_mass) / UC_elec_volt
            call UT_linInterpList(cross_secs(iXc)%en_cs(1, :), &
                 cross_secs(iXc)%en_cs(2, :), en_eV, rate_list(i))
            rate_list(i) = rate_list(i) * vel_list(i)
            if (vel_list(i) /= 0.d0) tempRate(n,i) = rate_list(i) / vel_list(i) /GAS_get_number_dens()
         end do
         call LT_add_column(PC_coll%rate_lt, vel2_list, rate_list)
         sum_rate_list = sum_rate_list + rate_list
      end do
      ! Create a new table with the sums of the rates. This can speed up the null collision method.
      PC_coll%sum_rate_lt => LT_create(minval(vel2_list), maxval(vel2_list), table_size)
      call LT_add_column(PC_coll%sum_rate_lt, vel2_list, sum_rate_list)

      PC_coll%max_rate = maxval(sum_rate_list)
      PC_coll%inv_max_rate = 1 / PC_coll%max_rate


      ! output the cross section in the lookup-table, to check whether it is correct
      myunit = 112233
      filename = "output/" // "output_elec_cs_test.txt"
      open(unit = myunit, file = filename , status = 'unknown')
      write(myunit, *) "Here we output the cross section of electrons. For Testing!!"
      do i = 1, table_size
          write(myunit, *) 0.5d0 * vel2_list(i)* PC_particle_elec_mass / UC_elec_volt, tempRate(:,i)
      end do
      close(myunit)
   end subroutine create_coll_rate_table

   !> Anbang: Create a lookup table with cross sections for IONS
   !> And the mass here should be the COM mass, which is m1 *m2/ (m1+m2). Or we use real mass of ions in the lab frame
   !>  attention, here we use energy instead of velocity
   subroutine create_coll_rate_table_for_ions_input_EneV(cross_secs_for_ions, min_eV, max_eV, table_size)
      use m_units_constants
      use m_cross_sec_for_ions
      use m_utils
      use m_gas
      type(CS_ion_type), intent(in) :: cross_secs_for_ions(:)
      integer, intent(IN)       :: table_size
      real(dp), intent(in)      :: min_eV, max_eV

      real(dp)                  :: en_list(table_size), rate_list(table_size), sum_rate_list(table_size)
      integer                   :: i, n, m, iXc
      integer, allocatable      :: collIxs(:)
      real(dp)                  :: en_eV, temp
      character(len=name_len) :: filename
      integer              :: myunit
      real(dp), allocatable     :: tempRate(:,:)


         
  !     min_vel = en_to_vel(min_eV * UC_elec_volt, PC_particle_ion_mass)
  !     max_vel = en_to_vel(max_eV * UC_elec_volt, PC_particle_ion_mass)

    !Anbang: Here for turner's model, we use reduced mass   
!       min_vel = en_to_vel(min_eV * UC_elec_volt, PC_reduced_mass)
!       max_vel = en_to_vel(max_eV * UC_elec_volt, PC_reduced_mass)

      PC_ion_coll%num = count(cross_secs_for_ions(:)%min_energy < max_eV .and. cross_secs_for_ions(:)%max_energy > min_eV)
!       print *, "PC_ion_coll%num =", PC_ion_coll%num 
      allocate( collIxs(PC_ion_coll%num) )
      allocate( PC_ion_coll%types(PC_ion_coll%num) )
      allocate( PC_ion_coll%spec_values(PC_ion_coll%num) )
      allocate( PC_ion_coll%buff(PC_ion_coll%num) )
      allocate( tempRate(PC_ion_coll%num, table_size))
      tempRate = 0.d0

      ! Get the indexes of the collisions in this energy range
      m = 1
      do n = 1, size(cross_secs_for_ions)
         if (cross_secs_for_ions(n)%min_energy < max_eV .and. cross_secs_for_ions(n)%max_energy > min_eV) then
            collIxs(m)  = n
            m           = m + 1
         end if
      end do

      ! Set the energy for the lookup table
      do i = 1, table_size
         temp = real(i-1, dp) / (table_size-1)
         en_list(i) = min_eV + temp * (max_eV - min_eV)
      end do

      ! Create collision rate table
      PC_ion_coll%rate_lt => LT_create(min_eV, max_eV, table_size, n_storage = PC_ion_coll%num)

      sum_rate_list = 0.0_dp

      do n = 1, PC_ion_coll%num
         iXc = collIxs(n)
         PC_ion_coll%types(n) = cross_secs_for_ions(iXc)%col_type

         ! elastic scattering and backwardscattering
         select case (PC_ion_coll%types(n))
         case (CS_ion_elastic_t, CS_ion_backwardscattering_t)
            PC_ion_coll%spec_values(n) = cross_secs_for_ions(iXc)%spec_value
       !     print *,"spec_values of ions (mass ration between ion and neutral molecular):", PC_ion_coll%spec_values(n)
         case DEFAULT
            PC_ion_coll%spec_values(n) = 0.0_dp
         end select

         do i = 1, table_size
            
  !          en_eV = vel_to_en(vel_list(i), PC_particle_ion_mass) / UC_elec_volt
!            en_eV = vel_to_en(vel_list(i), PC_reduced_mass) / UC_elec_volt
             en_eV = en_list(i)

            call UT_linInterpList(cross_secs_for_ions(iXc)%en_cs(1, :), &
                 cross_secs_for_ions(iXc)%en_cs(2, :), en_eV, rate_list(i))
          !  rate_list(i) = rate_list(i) * vel_list(i)
            rate_list(i) = rate_list(i) * en_to_vel(en_list(i)*UC_elec_volt, PC_reduced_mass)  !!Anbang: here which mass we should use?

            if (en_list(i) /= 0.d0) tempRate(n,i) = rate_list(i) / &
                    & en_to_vel(en_list(i)*UC_elec_volt, PC_reduced_mass) /GAS_get_number_dens()
         end do
         call LT_add_column(PC_ion_coll%rate_lt, en_list, rate_list)
         sum_rate_list = sum_rate_list + rate_list
      end do
      ! Create a new table with the sums of the rates. This can speed up the null collision method.
      PC_ion_coll%sum_rate_lt => LT_create(minval(en_list), maxval(en_list), table_size)
      call LT_add_column(PC_ion_coll%sum_rate_lt, en_list, sum_rate_list)

      PC_ion_coll%max_rate = maxval(sum_rate_list)
      PC_ion_coll%inv_max_rate = 1.d0 / PC_ion_coll%max_rate

      ! output the cross section in the lookup-table, to check whether it is correct
      myunit = 332211
      filename = "output/" // "output_ion_cs_test.txt"
      open(unit = myunit, file = filename , status = 'unknown')
      write(myunit, *) "Here we output the cross section of ions. For Testing!!"
      do i = 1, table_size
          write(myunit, *) en_list(i), tempRate(:,i)
      end do
      close(myunit)    

   end subroutine create_coll_rate_table_for_ions_input_EneV


  ! use velocity as input data instead of energy in eV.  This is the same as electrons
   subroutine create_coll_rate_table_for_ions_input_vel(cross_secs_for_ions, min_eV, max_eV, table_size)
      use m_units_constants
      use m_cross_sec_for_ions
      use m_utils
      use m_gas

      type(CS_ion_type), intent(in) :: cross_secs_for_ions(:)
      integer, intent(IN)       :: table_size
      real(dp), intent(in)      :: min_eV, max_eV

      real(dp)                  :: vel_list(table_size), rate_list(table_size), sum_rate_list(table_size)
      real(dp)                  :: vel2_list(table_size)
      integer                   :: i, n, m, iXc
      integer, allocatable      :: collIxs(:)
      real(dp)                  :: en_eV, temp, min_vel, max_vel, min_vel2, max_vel2
      character(len=name_len)   :: filename
      integer                   :: myunit
      real(dp), allocatable     :: tempRate(:,:)


    !Anbang: Here for turner's model, we use reduced mass   
       min_vel = en_to_vel(min_eV * UC_elec_volt, PC_reduced_mass)
       max_vel = en_to_vel(max_eV * UC_elec_volt, PC_reduced_mass)

       min_vel2 = min_vel * min_vel
       max_vel2 = max_vel * max_vel

      PC_ion_coll%num = count(cross_secs_for_ions(:)%min_energy < max_eV .and. cross_secs_for_ions(:)%max_energy > min_eV)
     ! print *, "PC_ion_coll%num =", PC_ion_coll%num 
      allocate( collIxs(PC_ion_coll%num) )
      allocate( PC_ion_coll%types(PC_ion_coll%num) )
      allocate( PC_ion_coll%spec_values(PC_ion_coll%num) )
      allocate( PC_ion_coll%buff(PC_ion_coll%num) )
      allocate( tempRate(PC_ion_coll%num, table_size))
      tempRate = 0.d0

      ! Get the indexes of the collisions in this energy range
      m = 1
      do n = 1, size(cross_secs_for_ions)
         if (cross_secs_for_ions(n)%min_energy < max_eV .and. cross_secs_for_ions(n)%max_energy > min_eV) then
            collIxs(m)  = n
            m           = m + 1
         end if
      end do

      ! Set the velocities for the lookup table
      do i = 1, table_size
         temp = real(i-1, dp) / (table_size-1)
         vel_list(i) = min_vel + temp * (max_vel - min_vel)
         vel2_list(i) = vel_list(i) **2.d0
      end do

      ! Create collision rate table
      PC_ion_coll%rate_lt => LT_create(min_vel2, max_vel2, table_size, n_storage = PC_ion_coll%num)

      sum_rate_list = 0.0_dp

      do n = 1, PC_ion_coll%num
         iXc = collIxs(n)
         PC_ion_coll%types(n) = cross_secs_for_ions(iXc)%col_type

         ! elastic scattering and backwardscattering
         select case (PC_ion_coll%types(n))
         case (CS_ion_elastic_t, CS_ion_backwardscattering_t)
            PC_ion_coll%spec_values(n) = cross_secs_for_ions(iXc)%spec_value
         !   print *,"spec_values of ions (mass ration between ion and neutral molecular):", PC_ion_coll%spec_values(n)
         case DEFAULT
            PC_ion_coll%spec_values(n) = 0.0_dp
         end select

         do i = 1, table_size
            
            en_eV = vel_to_en(vel_list(i), PC_reduced_mass) / UC_elec_volt
            call UT_linInterpList(cross_secs_for_ions(iXc)%en_cs(1, :), &
                 cross_secs_for_ions(iXc)%en_cs(2, :), en_eV, rate_list(i))
            rate_list(i) = rate_list(i) * vel_list(i)
            if (vel_list(i) /= 0.d0) tempRate(n,i) = rate_list(i)/ vel_list(i)/ GAS_get_number_dens()
         end do
         call LT_add_column(PC_ion_coll%rate_lt, vel2_list, rate_list)
         sum_rate_list = sum_rate_list + rate_list
      end do
      ! Create a new table with the sums of the rates. This can speed up the null collision method.
      PC_ion_coll%sum_rate_lt => LT_create(minval(vel2_list), maxval(vel2_list), table_size)
      call LT_add_column(PC_ion_coll%sum_rate_lt, vel2_list, sum_rate_list)

      PC_ion_coll%max_rate = maxval(sum_rate_list)
      PC_ion_coll%inv_max_rate = 1.d0 / PC_ion_coll%max_rate

      ! output the cross section in the lookup-table, to check whether it is correct
      myunit = 331122
      filename = "output/" // "output_ion_cs_test.txt"
      open(unit = myunit, file = filename , status = 'unknown')
      write(myunit, *) "Here we output the cross section of ion. For Testing!!"
      do i = 1, table_size
          write(myunit, *) 0.5d0 * vel2_list(i) * PC_reduced_mass / UC_elec_volt, tempRate(:,i)
      end do
      close(myunit)

   end subroutine create_coll_rate_table_for_ions_input_vel

   ! Routine to merge and split particles. Input arguments are the coordinate weights, used
   ! to determine the 'distance' between particles. The first three elements of the array are
   ! the weights of the xyz position coordinates, the next three the weights of the xyz
   ! velocity coordinates. Max_distance is the maxium Euclidean distance between particles
   ! to be merged. The weight_func returns the desired weight for a particle, whereas the
   ! pptr_merge and pptr_split procedures merge and split particles.
   subroutine PC_merge_and_split(x_mask, v_fac, use_v_norm, weight_func, &
       pptr_merge, pptr_split)
      use m_mrgrnk
      use kdtree2_module

      real(dp), intent(in)       :: v_fac
      logical, intent(in)        :: x_mask(3), use_v_norm

      interface

         subroutine pptr_merge(part_a, part_b,part_out)
            ! Takes in two particles,remove them and generate one new particle
            import
            type(PC_part_t), intent(in) :: part_a, part_b
            type(PC_part_t), intent(out):: part_out
         end subroutine pptr_merge

         subroutine pptr_split(part_a,w_ratio, part_out,n_part_out)
            ! Split part_a into two, with the other 'half' put into part_out
            import
            type(PC_part_t), intent(in) :: part_a
            real(dp), intent(in)           :: w_ratio
            type(PC_part_t), intent(inout) :: part_out(:)
            integer, intent(inout)         :: n_part_out           
         end subroutine pptr_split
      end interface

      procedure(p_to_r_f)          :: weight_func

      integer, parameter           :: num_neighbors = 1
      integer, parameter           :: n_part_out_max = 16
      real(dp), parameter          :: large_ratio   = 1.5_dp, small_ratio = 1 / large_ratio
      type(kdtree2), pointer       :: kd_tree
      type(kdtree2_result)         :: kd_results(num_neighbors)

      integer                      :: n_x_coord, n_coords
      integer                      :: num_part, num_merge, num_split
      integer                      :: p_min, p_max, n_too_far
      integer               :: o_ix, o_nn_ix
      integer               :: i, ix, neighbor_ix
      integer               :: n_part_out
      logical, allocatable  :: already_merged(:)
      integer, allocatable  :: sorted_ixs(:), coord_ixs(:)
      real(dp), allocatable :: coord_data(:, :), weight_ratios(:)
      type(PC_part_t)       :: part_out(n_part_out_max)

      p_min                                    = 1
      p_max                                    = PC_num_part
      num_part                                 = p_max - p_min + 1

      allocate(weight_ratios(num_part))
      allocate(sorted_ixs(num_part))

      do ix = 1, num_part
         weight_ratios(ix) = PC_particles(p_min+ix-1)%weight / weight_func(PC_particles(p_min+ix-1))
      end do
      num_merge = count(weight_ratios <= small_ratio)
      n_x_coord = count(x_mask)     ! position(1d-3d)
      n_coords  = n_x_coord + 3     ! position + velocity? 
      if (use_v_norm) n_coords = n_coords - 2

      allocate(coord_data(n_coords, num_merge))
      allocate(coord_ixs(n_coords))
      allocate(already_merged(num_merge))
      already_merged = .false.
      n_too_far = 0

     !   print *, "Electron: num_merge = ", num_merge

        ! Sort particles by their relative weight
        call mrgrnk(weight_ratios, sorted_ixs)

        ! Only create a k-d tree if there are enough particles to be merged
   if (num_merge > n_coords) then
        do ix = 1, num_merge
            o_ix = sorted_ixs(ix)
            coord_data(1:n_x_coord, ix) = pack(PC_particles(o_ix)%x, x_mask)
            if (use_v_norm) then
                coord_data(n_x_coord+1, ix) = v_fac * norm2(PC_particles(o_ix)%v)
            else
                coord_data(n_x_coord+1:, ix) = v_fac * PC_particles(o_ix)%v
            end if
        end do

         ! Create k-d tree
         kd_tree => kdtree2_create(coord_data)

       ! Merge particles
       do ix = 1, num_merge
          if (already_merged(ix)) cycle

          call kdtree2_n_nearest_around_point(kd_tree, idxin=ix, &
               nn=num_neighbors, correltime=1, results=kd_results)
          neighbor_ix = kd_results(1)%idx

          if (already_merged(neighbor_ix)) cycle

          ! Get indices in the original particle list
          o_ix = sorted_ixs(ix)
          o_nn_ix = sorted_ixs(neighbor_ix)

          ! Merge, then remove neighborCS_ion_elastic_t
          call pptr_merge(PC_particles(o_ix), PC_particles(o_nn_ix), &
               part_out(1))
          PC_particles(o_ix) = part_out(1)
          call PC_killElec(o_nn_ix)     ! make it 'inactive' that can be removed later
          !PC_particles(o_nn_ix)%weight  = PC_dead_weight   ! give it a large weight that will be remover later
          already_merged((/ix, neighbor_ix/)) = .true.
       end do

       call kdtree2_destroy(kd_tree)
    end if

    ! Split particles.
    num_split = count(weight_ratios >= large_ratio)
  !  print *, "Electron: num_split = ", num_split

    do ix = num_part - num_split + 1, num_part
       o_ix = sorted_ixs(ix)
       call pptr_split(PC_particles(o_ix), weight_ratios(o_ix), part_out, &
            n_part_out)
       PC_particles(o_ix) = part_out(1)
        
    !> After split, the original particle has the information of one new particle;
    !> while the second particle is added to the list
       do i = 2, n_part_out
          PC_num_part = PC_num_part + 1
          PC_particles(PC_num_part) = part_out(i)
       end do

    end do

    
   
   ! Remove dead particles that were marked with 'inactive'
   call remove_dead_particles()

  ! print *, "Electron: the min/ max weight =", minval(PC_particles(1:PC_num_part)%weight),&
   !                     & maxval(PC_particles(1:PC_num_part)%weight)

   end subroutine PC_merge_and_split

! Merge two particles into part_a, should remove part_b afterwards
  subroutine PC_merge_part_rxv(part_a, part_b, part_out)
    use m_random
    type(PC_part_t), intent(in)    :: part_a, part_b
    type(PC_part_t), intent(out) :: part_out
    !type(RNG_t), intent(inout) :: rng

   if (PC_merge_scheme == 1)  then ! select position from one of the particle
        if (RNG_uniform() > part_a%weight / (part_a%weight + part_b%weight)) then
            part_out%x      = part_b%x
            part_out%v      = part_b%v
            part_out%a      = part_b%a
            part_out%T_left = part_b%T_left
            part_out%live   = part_b%live
        else
            part_out%x      = part_a%x
            part_out%v      = part_a%v
            part_out%a      = part_a%a
            part_out%T_left = part_a%T_left
            part_out%live   = part_a%live
        end if
     else if (PC_merge_scheme == 2) then
        if (RNG_uniform() > part_a%weight / (part_a%weight + part_b%weight)) then
            part_out%v      = part_b%v
            part_out%a      = part_b%a
            part_out%T_left = part_b%T_left
            part_out%live   = part_b%live
        else
            part_out%v      = part_a%v
            part_out%a      = part_a%a
            part_out%T_left = part_a%T_left
            part_out%live   = part_a%live
        end if
            !weighted average 
            part_out%x      = (part_a%x * part_a%weight + part_b%x * part_b%weight)/ &
                            &  (part_a%weight + part_b%weight )

     else
        print *, "PC_merge_part_rxv Error: the merge scheme is not defined!"
        stop
     end if
    part_out%weight = part_a%weight + part_b%weight
  end subroutine PC_merge_part_rxv

  subroutine PC_split_part(part_a, w_ratio, part_out, n_part_out)
    type(PC_part_t), intent(in)    :: part_a
    real(dp), intent(in)           :: w_ratio
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(inout)         :: n_part_out
    !type(RNG_t), intent(inout)     :: rng

    n_part_out    = 2
    part_out(1)   = part_a
    part_out(1)%weight = 0.5_dp * part_out(1)%weight
    part_out(2)   = part_out(1)
  end subroutine PC_split_part

   integer function PC_get_num_colls()
      PC_get_num_colls = PC_coll%num
   end function PC_get_num_colls

   subroutine PC_get_colls(out_colls)
      type(PC_coll_t), intent(out) :: out_colls
      out_colls = PC_coll
   end subroutine PC_get_colls

   !> Returns the number of collisions as real, to avoid overflow with standard integers
   subroutine PC_get_coll_count(collCount)
      real(dp), intent(out) :: collCount(:)
      collCount = PC_coll%counts
   end subroutine PC_get_coll_count

   !> Reset the number of collisions
   subroutine PC_reset_coll_count()
      PC_coll%counts(:) = 0
   end subroutine PC_reset_coll_count

    !> Here to check the minimum and maximum position,acceleration of electrons and ions
   subroutine PC_check_pos_particles()
      real(dp)  :: posElecMin, posElecMax, PosIonMin, PosIonMax, accelElecMax, accelIonMax
      posElecMin = minval(PC_particles(:)%x(3))
      posElecMax = maxval(PC_particles(:)%x(3))
      posIonMin  = minval(PC_ion_particles(:)%x(3))
      posIonMax  = maxval(PC_ion_particles(:)%x(3))
      accelElecMax = maxval(abs(PC_particles(:)%a(3)))
      accelIonMax  = maxval(abs(PC_ion_particles(:)%a(3)))
!       print *, "The minimum and maximum positions of electrons are: ", posElecMin, posElecMax
!       print *, "The minimum and maximum positions of ions are: ", posIonMin, posIonMax
! 
!       print *, "accelElecMax and accelIonMax are", accelElecMax, accelIonMax     
   end subroutine PC_check_pos_particles


    !> Here to check the minimum and maximum weight of particles
   subroutine PC_check_weight_particles()

      integer :: ll

!       print *, "The min/max weight of electrons are: ", minval(PC_particles(1:PC_num_part)%weight),&
!                             & maxval(PC_particles(1:PC_num_part)%weight)
!       print *, "The min/max weight of ions are: ", minval(PC_ion_particles(1:PC_num_ion_part)%weight), &
!                             & maxval(PC_ion_particles(1:PC_num_ion_part)%weight)
        
!       print *, "the size compare:", size(PC_particles(:)%weight), PC_num_part
! 
!       ll = 1
!       do while (PC_particles(ll)%weight == 0.d0 .and. (ll <= PC_num_part))
!         print *, "ll = ", ll
!         print *, " PC_particles(ll)%x(3) = ", PC_particles(ll)%x(3) 
!         ll = ll + 1
!       end do
!       

   end subroutine PC_check_weight_particles


    !!!We calculate EEDF for electrons
    subroutine PC_histogram(hist_func, filter_f, &
        filter_args, x_values, y_values, z_values)
        use m_mrgrnk
        use m_units_constants
        procedure(p_to_r_f)   :: hist_func
        real(dp), intent(in)        :: x_values(:)
        real(dp), intent(out)       :: y_values(:)
        real(dp), intent(out)       :: z_values(:) !EEPF
        interface
        logical function filter_f(my_part, real_args)
            import
            type(PC_part_t), intent(in) :: my_part
            real(dp), intent(in) :: real_args(:)
        end function filter_f
        end interface
        real(dp), intent(in)        :: filter_args(:)

        integer                     :: ix, p_ix, o_ix, n_used, num_bins, n_part
        integer, allocatable        :: sorted_ixs(:)
        real(dp)                    :: boundary_value
        real(dp), allocatable       :: part_values(:)
        logical, allocatable        :: part_mask(:)
        real(dp)                    :: total_part_counted

        n_part = PC_num_part
        allocate(part_mask(n_part))
        do ix = 1, n_part
        part_mask(ix) = filter_f(PC_particles(ix), filter_args)
        end do

        n_used = count(part_mask)

        allocate(part_values(n_used))
        allocate(sorted_ixs(n_used))

        p_ix = 0
        do ix = 1, PC_num_part
        if (part_mask(ix)) then
            p_ix = p_ix + 1
            part_values(p_ix) = hist_func(PC_particles(p_ix))
        end if
        end do

        !Anbang: the total number of counted particles
        total_part_counted = 0.d0
        do ix = 1, PC_num_part
            if (part_mask(ix)) total_part_counted = total_part_counted + PC_particles(ix)%weight
        end do

        call mrgrnk(part_values, sorted_ixs)

        num_bins = size(x_values)
        p_ix     = 1
        y_values = 0.d0
        z_values = 0.d0

        outer: do ix = 1, num_bins-1
        boundary_value = 0.5_dp * (x_values(ix) + x_values(ix+1))
        do
            if (p_ix == n_used + 1) exit outer
            o_ix = sorted_ixs(p_ix) ! Index in the 'old' list
            if (part_values(o_ix) > boundary_value) exit

            y_values(ix) = y_values(ix) + PC_particles(o_ix)%weight
            !! EEPF 
            z_values(ix) = z_values(ix) + PC_particles(o_ix)%weight / sqrt(x_values(ix) * 0.5d0 * UC_elec_mass / UC_elec_Volt)
            p_ix         = p_ix + 1
        end do
        end do outer

        ! Fill last bin
        y_values(num_bins) = sum(PC_particles(sorted_ixs(p_ix:n_used))%weight)
        z_values(num_bins) = sum(PC_particles(sorted_ixs(p_ix:n_used))%weight) / &
                            sqrt(x_values(num_bins) * 0.5d0 * UC_elec_mass / UC_elec_Volt)


        if (total_part_counted /= 0.d0) then
            y_values = y_values / total_part_counted
            z_values = z_values / total_part_counted
        else
            y_values = 0.d0
            z_values = 0.d0
        end if
    end subroutine PC_histogram



!!!!!!!!!!!*********************************************************************************
!!!!!!!!!!For ions merge and split *********************************************************
!!!!!!!!!************************************************************************************
   ! Routine to merge and split ion particles. Input arguments are the coordinate weights, used
   ! to determine the 'distance' between particles. The first three elements of the array are
   ! the weights of the xyz position coordinates, the next three the weights of the xyz
   ! velocity coordinates.
   subroutine PC_ion_merge_and_split(x_mask, v_fac, use_v_norm, weight_func, &
       pptr_ion_merge, pptr_ion_split)
      use m_mrgrnk
      use kdtree2_module

      real(dp), intent(in)       :: v_fac
      logical, intent(in)        :: x_mask(3), use_v_norm
        
      interface

         subroutine pptr_ion_merge(part_a, part_b,part_out)
            ! Takes in two particles,remove them and generate one new particle
            import
            type(PC_ion_part_t), intent(in) :: part_a, part_b
            type(PC_ion_part_t), intent(out):: part_out
         end subroutine pptr_ion_merge

         subroutine pptr_ion_split(part_a,w_ratio, part_out,n_part_out)
            ! Split part_a into two, with the other 'half' put into part_out
            import
            type(PC_ion_part_t), intent(in) :: part_a
            real(dp), intent(in)           :: w_ratio
            type(PC_ion_part_t), intent(inout) :: part_out(:)
            integer, intent(inout)         :: n_part_out           
         end subroutine pptr_ion_split
      end interface
      procedure(p_to_r_f_ion)          :: weight_func

      integer, parameter           :: num_neighbors = 1
      integer, parameter           :: n_part_out_max = 16
      real(dp), parameter          :: large_ratio   = 1.5_dp, small_ratio = 1 / large_ratio
      type(kdtree2), pointer       :: kd_tree
      type(kdtree2_result)         :: kd_results(num_neighbors)

      integer                      :: n_x_coord, n_coords
      integer                      :: num_part, num_merge, num_split
      integer                      :: p_min, p_max, n_too_far
      integer               :: o_ix, o_nn_ix
      integer               :: i, ix, neighbor_ix
      integer               :: n_part_out
      logical, allocatable  :: already_merged(:)
      integer, allocatable  :: sorted_ixs(:), coord_ixs(:)
      real(dp), allocatable :: coord_data(:, :), weight_ratios(:)
      type(PC_ion_part_t)       :: part_out(n_part_out_max)

      p_min                                    = 1
      p_max                                    = PC_num_ion_part
      num_part                                 = p_max - p_min + 1

      allocate(weight_ratios(num_part))
      allocate(sorted_ixs(num_part))

      do ix = 1, num_part
         weight_ratios(ix) = PC_ion_particles(p_min+ix-1)%weight / weight_func(PC_ion_particles(p_min+ix-1))
      end do
      num_merge = count(weight_ratios <= small_ratio)
      n_x_coord = count(x_mask)     ! position(1d-3d)
      n_coords  = n_x_coord + 3     ! position + velocity? 
      if (use_v_norm) n_coords = n_coords - 2

      allocate(coord_data(n_coords, num_merge))
      allocate(coord_ixs(n_coords))
      allocate(already_merged(num_merge))
      already_merged = .false.
      n_too_far = 0

       ! print *, "Ions: num_merge = ", num_merge

        ! Sort particles by their relative weight
        call mrgrnk(weight_ratios, sorted_ixs)

        ! Only create a k-d tree if there are enough particles to be merged
   if (num_merge > n_coords) then
        do ix = 1, num_merge
            o_ix = sorted_ixs(ix)
            coord_data(1:n_x_coord, ix) = pack(PC_ion_particles(o_ix)%x, x_mask)
            if (use_v_norm) then
                coord_data(n_x_coord+1, ix) = v_fac * norm2(PC_ion_particles(o_ix)%v)
            else
                coord_data(n_x_coord+1:, ix) = v_fac * PC_ion_particles(o_ix)%v
            end if
        end do
     

         ! Create k-d tree
         kd_tree => kdtree2_create(coord_data)

       ! Merge particles
       do ix = 1, num_merge
          if (already_merged(ix)) cycle

          call kdtree2_n_nearest_around_point(kd_tree, idxin=ix, &
               nn=num_neighbors, correltime=1, results=kd_results)
          neighbor_ix = kd_results(1)%idx

          if (already_merged(neighbor_ix)) cycle

          ! Get indices in the original particle list
          o_ix = sorted_ixs(ix)
          o_nn_ix = sorted_ixs(neighbor_ix)

          ! Merge, then remove neighbor
          call pptr_ion_merge(PC_ion_particles(o_ix), PC_ion_particles(o_nn_ix), &
               part_out(1))
          PC_ion_particles(o_ix) = part_out(1)
          call PC_killIon(o_nn_ix)     ! make it 'inactive' that can be removed later
          already_merged((/ix, neighbor_ix/)) = .true.
       end do

       call kdtree2_destroy(kd_tree)
    end if

    ! Split particles.
    num_split = count(weight_ratios >= large_ratio)
  !  print *, "Ions: num_split = ", num_split

    do ix = num_part - num_split + 1, num_part
       o_ix = sorted_ixs(ix)
       call pptr_ion_split(PC_ion_particles(o_ix), weight_ratios(o_ix), part_out, &
            n_part_out)
       PC_ion_particles(o_ix) = part_out(1)
        
    !> After split, the original particle has the information of one new particle;
    !> while the second particle is added to the list
       do i = 2, n_part_out
          PC_num_ion_part = PC_num_ion_part + 1
          PC_ion_particles(PC_num_ion_part) = part_out(i)
       end do

    end do
   
   ! Remove dead particles that were marked with 'inactive'
   call remove_dead_ion_particles()

  ! print *, "Ion: the min/ max weight =", minval(PC_ion_particles(1:PC_num_ion_part)%weight), &
   !         & maxval(PC_ion_particles(1:PC_num_ion_part)%weight)

   end subroutine PC_ion_merge_and_split

! Merge two particles into part_a, should remove part_b afterwards
  subroutine PC_ion_merge_part_rxv(part_a, part_b, part_out)
    use m_random
    type(PC_ion_part_t), intent(in)    :: part_a, part_b
    type(PC_ion_part_t), intent(out) :: part_out
    !type(RNG_t), intent(inout) :: rng

   if (PC_merge_scheme == 1)  then ! select position from one of the particle
        if (RNG_uniform() > part_a%weight / (part_a%weight + part_b%weight)) then
            part_out%x      = part_b%x
            part_out%v      = part_b%v
            part_out%a      = part_b%a
            part_out%T_left = part_b%T_left
            part_out%live   = part_b%live
        else
            part_out%x      = part_a%x
            part_out%v      = part_a%v
            part_out%a      = part_a%a
            part_out%T_left = part_a%T_left
            part_out%live   = part_a%live
        end if
     else if (PC_merge_scheme == 2) then
        if (RNG_uniform() > part_a%weight / (part_a%weight + part_b%weight)) then
            part_out%v      = part_b%v
            part_out%a      = part_b%a
            part_out%T_left = part_b%T_left
            part_out%live   = part_b%live
        else
            part_out%v      = part_a%v
            part_out%a      = part_a%a
            part_out%T_left = part_a%T_left
            part_out%live   = part_a%live
        end if
            !weighted average 
            part_out%x      = (part_a%x * part_a%weight + part_b%x * part_b%weight)/ &
                            &  (part_a%weight + part_b%weight )

     else
        print *, "PC_merge_part_rxv Error: the merge scheme is not defined!"
        stop
     end if
    part_out%weight = part_a%weight + part_b%weight
  end subroutine PC_ion_merge_part_rxv

  subroutine PC_ion_split_part(part_a, w_ratio, part_out, n_part_out)
    type(PC_ion_part_t), intent(in)    :: part_a
    real(dp), intent(in)           :: w_ratio
    type(PC_ion_part_t), intent(inout) :: part_out(:)
    integer, intent(inout)         :: n_part_out
    !type(RNG_t), intent(inout)     :: rng

    n_part_out    = 2
    part_out(1)   = part_a
    part_out(1)%weight = 0.5_dp * part_out(1)%weight
    part_out(2)   = part_out(1)
  end subroutine PC_ion_split_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!Anbang:Here we develop routines for boundary conditions for electrons and ions!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine PC_boundary_for_electron(ll, pos, dt)
        use m_error
        use m_random

        integer, intent(in)  :: ll
        real(dp), intent(in) :: pos, dt
    
        if (pos <= PC_posBD(1)) then

            select case(PC_bdType_elec(1))
            case(1)   ! absorb
                call PC_killElec(ll)
            PC_Q_bd(1)=PC_Q_bd(1)-PC_particles(ll)%weight
            case(2)   ! partial reflection
                if (RNG_uniform() <= PC_coElecRefl(1)) then
                    !call specularReflectionForElec(ll)
                    call reflectionTskhakayaForElec(ll,dt)
                else
                     call PC_killElec(ll)
                PC_Q_bd(1)=PC_Q_bd(1)-PC_particles(ll)%weight
                end if
            end select
        else if (pos >= PC_posBD(2)) then
            select case(PC_bdType_elec(2))
            case(1)   ! absorb
                call PC_killElec(ll)
            PC_Q_bd(2)=PC_Q_bd(2)-PC_particles(ll)%weight
            case(2)   ! partial reflection
                if (RNG_uniform() <= PC_coElecRefl(2)) then
                    !call specularReflectionForElec(ll)
                    call reflectionTskhakayaForElec(ll,dt)
                else
                    call PC_killElec(ll)
                PC_Q_bd(2)=PC_Q_bd(2)-PC_particles(ll)%weight
                end if
            end select
        else
        call ERR_show("PC_boundary_for_electron: the particle is still in the domain!")
        end if
    end subroutine PC_boundary_for_electron

  !!Anbang: specular reflection bd of electrons, we need to reconsider it,here is not accurate
   subroutine specularReflectionForElec(ll)
        integer,intent(in)   :: ll

        if (abs(PC_particles(ll)%a(3)) > 1.d18) print *, "the accel of aparticle ll is too big", ll, PC_particles(ll)%a(3) 

        if (PC_particles(ll)%x(3) <= PC_posBD(1)) then
            PC_particles(ll)%x(3) = PC_posBD(1) + (PC_posBD(1) - PC_particles(ll)%x(3))
            PC_particles(ll)%v(3) = -PC_particles(ll)%v(3)
        else if (PC_particles(ll)%x(3) >= PC_posBD(2)) then
            PC_particles(ll)%x(3) = PC_posBD(2) - (PC_particles(ll)%x(3) - PC_posBD(2))
            PC_particles(ll)%v(3) = -PC_particles(ll)%v(3)
        else
            print *, "Error(specularReflectionForElec): the electron is still in the domain!"
            stop
        end if

        if (isOutOfGas(PC_particles(ll)%x)) then
            print *, "Error(specularReflectionForElec): the reflected electron is still out of the domain?"
            print *, "The position of the electron:", PC_particles(ll)%x(3), PC_particles(ll)%v(3), PC_particles(ll)%a(3)
            stop
        end if
   end subroutine specularReflectionForElec

  !Anbang: a method provided by Tskhakaya, Matyash and Schneider et al, Contrib. plasma phys. 47, 563-594 (2007) 
   subroutine reflectionTskhakayaForElec(ll,dt)
        use m_random
        
        integer,intent(in)   :: ll
        real(dp),intent(in)  :: dt
        real(dp)             :: velAtWall, timeAtWall, velOld, temp(3)

!          print *, "before reflection, the position/vel is:", PC_particles(ll)%x(3), PC_particles(ll)%v(3), &
!              PC_particles(ll)%a(3) 
        temp(1) = PC_particles(ll)%x(3)
        temp(2) = PC_particles(ll)%v(3)
        temp(3) = PC_particles(ll)%a(3)

        if (abs(PC_particles(ll)%a(3)) > 1.d18) print *, "the accel of aparticle ll is too big", ll, PC_particles(ll)%a(3) 

        if (PC_particles(ll)%x(3) <= PC_posBD(1)) then
            velOld =  PC_particles(ll)%v(3) - PC_particles(ll)%a(3) * dt
            !timeAtWall = (PC_particles(ll)%x(3) - PC_posBD(1)) / velOld    ! Anbang: here is not very accurate? !!!
            timeAtWall  = dt * RNG_uniform()    ! i simply sample a time
!              print *, "timeAtWall/dt =",  timeAtWall, dt
            velAtWall  = velOld + PC_particles(ll)%a(3) * timeAtWall
!              print *, "velOLd/velAtWall =",  velOld, velAtWall

            PC_particles(ll)%v(3) = - velAtWall + PC_particles(ll)%a(3) * (dt - timeAtWall)
            PC_particles(ll)%x(3) = PC_posBD(1) + abs (( - velAtWall *(dt - timeAtWall) + &
                    0.5d0 * PC_particles(ll)%a(3) * (dt - timeAtWall)**2 ))               



        else if (PC_particles(ll)%x(3) >= PC_posBD(2)) then
            velOld =  PC_particles(ll)%v(3) - PC_particles(ll)%a(3) * dt
            timeAtWall = dt * RNG_uniform()
!              print *, "timeAtWall/dt =",  timeAtWall, dt
            velAtWall  = velOld + PC_particles(ll)%a(3) * timeAtWall

            PC_particles(ll)%v(3) = - velAtWall + PC_particles(ll)%a(3) * (dt - timeAtWall)
            PC_particles(ll)%x(3) = PC_posBD(2) - abs( (- velAtWall *(dt - timeAtWall) + &
                    0.5d0 * PC_particles(ll)%a(3) * (dt - timeAtWall)**2 ))

        else
            print *, "Error(reflectionTskhakayaForElec): reflection of electrons, electron is still in the domain? "
            stop
        end if

!          print *, "after reflection: position/vel of the electron:", PC_particles(ll)%x(3), PC_particles(ll)%v(3), &
!                  PC_particles(ll)%a(3) 

        ! further check
        if (isOutOfGas(PC_particles(ll)%x)) then
            print *, "Error(reflectionTskhakayaForElec): the reflected electron is still out of the domain?"
            print *, "The position of the electron:", PC_particles(ll)%x(3), PC_particles(ll)%v(3), PC_particles(ll)%a(3), temp
            stop
        end if
   end subroutine reflectionTskhakayaForElec

    subroutine PC_boundary_for_ion(ll, pos, dt)
        use m_error

        integer, intent(in)  :: ll
        real(dp), intent(in) :: pos, dt
    
        if (pos <= PC_posBD(1)) then

            select case(PC_bdType_Ion(1))
            case(1)   ! absorb
                call PC_killIon(ll)
            PC_Q_bd(1)=PC_Q_bd(1)+PC_ion_particles(ll)%weight
            case(2)   ! reflection
                call specularReflectionForion(ll)
            case(3)  
                call secondaryElectronEmissonForIon(ll)
            end select
        else if (pos >= PC_posBD(2)) then
            select case(PC_bdType_Ion(2))
            case(1)   ! absorb
                call PC_killIon(ll)
            PC_Q_bd(2)=PC_Q_bd(2)+PC_ion_particles(ll)%weight
            case(2)   ! reflection
                call specularReflectionForion(ll)
            case(3)  
                call secondaryElectronEmissonForIon(ll)
            end select
        else
            call ERR_show("PC_boundary_for_ion: the particle is still in the domain!")
        end if
    end subroutine PC_boundary_for_ion

  !!Anbang: specular reflection bd of ions
   subroutine specularReflectionForion(ll)
        integer,intent(in)   :: ll

        if (PC_ion_particles(ll)%x(3) <= PC_posBD(1)) then
            PC_ion_particles(ll)%x(3) = PC_posBD(1) + (PC_posBD(1) - PC_ion_particles(ll)%x(3))
            PC_ion_particles(ll)%v(3) = -PC_ion_particles(ll)%v(3)
        else if (PC_ion_particles(ll)%x(3) >= PC_posBD(2)) then
            PC_ion_particles(ll)%x(3) = PC_posBD(2) - (PC_ion_particles(ll)%x(3) - PC_posBD(2))
            PC_ion_particles(ll)%v(3) = -PC_ion_particles(ll)%v(3)
        else
            print *, "Error(specularReflectionForion): the electron is still in the domain!"
            stop
        end if
   end subroutine specularReflectionForion

  !Anbang: a method provided by Tskhakaya, Matyash and Schneider et al, Contrib. plasma phys. 47, 563-594 (2007) 
   subroutine reflectionTskhakayaForIon(ll,dt)
        use m_random

        integer,intent(in)   :: ll
        real(dp),intent(in)  :: dt
        real(dp)             :: velAtWall, timeAtWall, velOld

        if (PC_ion_particles(ll)%x(3) <= PC_posBD(1)) then
            velOld =  PC_ion_particles(ll)%v(3) - PC_ion_particles(ll)%a(3) * dt
            timeAtWall = dt * RNG_uniform()
            velAtWall  = velOld + PC_ion_particles(ll)%a(3) * timeAtWall
            PC_ion_particles(ll)%v(3) = - velAtWall + PC_ion_particles(ll)%a(3) * (dt - timeAtWall)
            PC_ion_particles(ll)%x(3) = PC_posBD(1) + abs(( - velAtWall *(dt - timeAtWall) + &
                        0.5d0 * PC_ion_particles(ll)%a(3) * (dt - timeAtWall)**2 ))

        else if (PC_ion_particles(ll)%x(3) >= PC_posBD(2)) then
            velOld =  PC_ion_particles(ll)%v(3) - PC_ion_particles(ll)%a(3) * dt
            timeAtWall = dt * RNG_uniform()
            velAtWall  = velOld + PC_ion_particles(ll)%a(3) * timeAtWall
            PC_ion_particles(ll)%v(3) = - velAtWall + PC_ion_particles(ll)%a(3) * (dt - timeAtWall)
            PC_ion_particles(ll)%x(3) = PC_posBD(2) - abs((- velAtWall *(dt - timeAtWall) + &
                    0.5d0 * PC_ion_particles(ll)%a(3) * (dt - timeAtWall)**2 ))
        else
            print *, "Error(reflectionTskhakayaForIon):the ion is still in the domain? "
            stop
        end if

        ! further check
        if (isOutOfGas(PC_ion_particles(ll)%x)) then
            print *, "Error(reflectionTskhakayaForIon): the reflected ion is still outof the domain?"
            print *, "The position of the ion:", PC_ion_particles(ll)%x(3)
            stop
        end if
   end subroutine reflectionTskhakayaForIon

   !Anbang: secondary electron emisson
    subroutine secondaryElectronEmissonForIon(ll)
        use m_random

        integer, intent(in)  :: ll

        if (PC_ion_particles(ll)%x(3) <= PC_posBD(1)) then
            if (RNG_uniform() <= PC_coSecElec(1)) then
                call PC_killIon(ll)
                ! a electron is emmited from the left wall
                call PC_create_second_elec(ll)
            PC_Q_bd(1)=PC_Q_bd(1)+2.d0*PC_ion_particles(ll)%weight
                !print *, "the left boundary has a sendary electron:", PC_particles(PC_num_part)%v
            else
                call PC_killIon(ll)
            PC_Q_bd(1)=PC_Q_bd(1)+PC_ion_particles(ll)%weight
            end if
        else if (PC_ion_particles(ll)%x(3) >= PC_posBD(2)) then
            if (RNG_uniform() <= PC_coSecElec(2)) then
                call PC_killIon(ll)
                ! a electron is emmited from the right wall
                call PC_create_second_elec(ll)
                !print *, "the right boundary has a sendary electron:", PC_particles(PC_num_part)%v
            PC_Q_bd(2)=PC_Q_bd(2)+2.d0*PC_ion_particles(ll)%weight
            else
                call PC_killIon(ll)
            PC_Q_bd(2)=PC_Q_bd(2)+PC_ion_particles(ll)%weight
            end if 
         else
            print *, "Error(secondaryElectronEmissonForIon):the ion is still in the domain? "
            stop
         end if       
        
    end subroutine secondaryElectronEmissonForIon

    !> Creat a secondary electron from the wall, 
    !> scattering angle is [0, pi/2] or [pi/2, pi]
    subroutine PC_create_second_elec(ll)
        use m_units_constants
        use m_random

        integer, intent(IN)   :: ll
        real(dp)              :: en_eV, vel_norm
        
        real(dp)              :: theta_sec, phi_sec

        PC_num_part = PC_num_part + 1
        call check_num_particles(PC_num_part)

        ! give an energy to the secondary electron
        ! in this case, we give it as an energy between 0 and 5 ev, with a maxwellian distribution
        ! the peak energy is PC_secElecEn = 2eV
        do 
            en_eV = RNG_normal() + PC_secElecEn   ! shift to right as PC_secElecEn
            if (en_eV <= 5.d0 .and.en_eV >0.d0) exit
        end do
!         print *, "**********************************************************"
!         print *, "we get a secondary eletron with an energy of (eV):", en_eV
        vel_norm    =  en_to_vel(en_eV * UC_elec_volt, PC_particle_elec_mass)

        PC_particles(PC_num_part)%weight   = PC_ion_particles(ll)%weight
        PC_particles(PC_num_part)%live     = .true.
        PC_particles(PC_num_part)%t_left   = 0.0_dp
        PC_particles(PC_num_part)%a    = - PC_ion_particles(ll)%a * PC_particle_ion_mass / PC_particle_elec_mass 
        ! accel can be decided from ion's accel

        if (PC_ion_particles(ll)%x(3) <= PC_posBD(1)) then
            PC_particles(PC_num_part)%x        = PC_posBD(1) + UC_tiny

            theta_sec = acos(RNG_uniform())  ! left boundary, the scattering angle is between 0 and pi/2
            phi_sec   = 2.d0 * UC_pi * RNG_uniform()
            PC_particles(PC_num_part)%v(1) = vel_norm * sin(theta_sec) * cos(phi_sec)
            PC_particles(PC_num_part)%v(2) = vel_norm * sin(theta_sec) * sin(phi_sec)
            PC_particles(PC_num_part)%v(3) = vel_norm * cos(theta_sec)
        else if (PC_ion_particles(ll)%x(3) >= PC_posBD(2)) then
            PC_particles(PC_num_part)%x        = PC_posBD(2) - UC_tiny

            theta_sec = acos(- RNG_uniform())  ! right boundary, the scattering angle is between pi/2 and pi
            phi_sec   = 2.d0 * UC_pi * RNG_uniform()
            PC_particles(PC_num_part)%v(1) = vel_norm * sin(theta_sec) * cos(phi_sec)
            PC_particles(PC_num_part)%v(2) = vel_norm * sin(theta_sec) * sin(phi_sec)
            PC_particles(PC_num_part)%v(3) = vel_norm * cos(theta_sec)
        end if

!         print *, "the velocity of incident ion is:", PC_ion_particles(ll)%v
!         print *, "the velocity of secondary electron is", PC_particles(PC_num_part)%v
        
    end subroutine PC_create_second_elec

!     !> Creat a secondary electron from the wall, 
!     !> Attention: this routine get a velocity depends on the incident ion
!     subroutine PC_create_second_elec(ll)
!         use m_units_constants
!         use m_random
! 
!         integer, intent(IN)   :: ll
!         real(dp)              :: vel(3), en_eV, vel_norm
!         real(dp)              :: theta_inc, phi_inc, radius
!         
!         real(dp)              :: theta_sec, phi_sec
! 
!         PC_num_part = PC_num_part + 1
!         call check_num_particles(PC_num_part)
! 
!         ! give an energy to the secondary electron
!         ! in this case, we give it as an energy between 0 and 5 ev, with a maxwellian distribution
!         ! the peak energy is PC_secElecEn = 2eV
! 
!         do 
!             en_eV = RNG_normal() + PC_secElecEn   ! shift to right as PC_secElecEn
!             if (en_eV <= 5.d0 .and.en_eV >=0.d0) exit
!         end do
! !         print *, "**********************************************************"
! !         print *, "we get a secondary eletron with an energy of (eV):", en_eV
!         vel_norm    =  en_to_vel(en_eV * UC_elec_volt, PC_particle_elec_mass)
! 
!         ! get the angles of incident ions
!         call UC_xyz_to_spherical(PC_ion_particles(ll)%v, radius, theta_inc, phi_inc)
!         !we get from a paper: A Monte Carlo investigation of secondary electron emission from solid
!             !targets.
!         ! two angles of electrons
!         theta_sec = asin(cos(theta_inc))  ! sin(theta_sec) = cos(theta_inc)
!         phi_sec   = UC_pi + phi_inc    ! phi_sec = phi_inc + pi
! 
!         PC_particles(PC_num_part)%weight   = PC_ion_particles(ll)%weight
!         PC_particles(PC_num_part)%live     = .true.
!         PC_particles(PC_num_part)%t_left   = 0.0_dp
!         PC_particles(PC_num_part)%a    = - PC_ion_particles(ll)%a * PC_particle_ion_mass / PC_particle_elec_mass 
!         ! accel can be decided from ion's accel
! 
!         if (PC_ion_particles(ll)%x(3) <= PC_posBD(1)) then
!             PC_particles(PC_num_part)%x        = PC_posBD(1) + UC_tiny
!             PC_particles(PC_num_part)%v(1) = vel_norm * sin(theta_sec) * cos(phi_sec)
!             PC_particles(PC_num_part)%v(2) = vel_norm * sin(theta_sec) * sin(phi_sec)
!             PC_particles(PC_num_part)%v(3) = vel_norm * cos(theta_sec)
!         else if (PC_ion_particles(ll)%x(3) >= PC_posBD(2)) then
!             PC_particles(PC_num_part)%x        = PC_posBD(2) - UC_tiny
!             PC_particles(PC_num_part)%v(1) = vel_norm * sin(theta_sec) * cos(phi_sec)
!             PC_particles(PC_num_part)%v(2) = vel_norm * sin(theta_sec) * sin(phi_sec)
!             PC_particles(PC_num_part)%v(3) = vel_norm * - cos(theta_sec)
!         end if
! 
! !         print *, "the velocity of incident ion is:", PC_ion_particles(ll)%v
! !         print *, "the velocity of secondary electron is", PC_particles(PC_num_part)%v
!         
!     end subroutine PC_create_second_elec
    
!     !> Creat a secondary electron from the wall, 
!     !> Attention: velocity of emitted electrons from the wall have to be checked!
!     subroutine PC_create_second_elec(ll)
!         use m_units_constants
!         use m_random
! 
!         integer, intent(IN)   :: ll
!         real(dp)              :: vel(3), en_ev
! 
!         PC_num_part = PC_num_part + 1
!         call check_num_particles(PC_num_part)
! 
!         vel     = (/RNG_normal(), RNG_normal(), RNG_normal()/)
!         vel     = vel * sqrt(2 * PC_secElecEn * UC_elec_volt / (3 * UC_elec_mass))
! 
!         PC_particles(PC_num_part)%weight   = PC_ion_particles(ll)%weight
!         PC_particles(PC_num_part)%live     = .true.
!         PC_particles(PC_num_part)%t_left   = 0.0_dp
!         if (PC_ion_particles(ll)%x(3) <= PC_posBD(1)) then
!             PC_particles(PC_num_part)%x        = PC_posBD(1) + UC_tiny
!             PC_particles(PC_num_part)%v        = vel
!             PC_particles(PC_num_part)%v(3) = abs(vel(3))  ! to keep it as positve velocity
!         else if (PC_ion_particles(ll)%x(3) >= PC_posBD(2)) then
!             PC_particles(PC_num_part)%x        = PC_posBD(2) - UC_tiny
!             PC_particles(PC_num_part)%v        = vel
!             PC_particles(PC_num_part)%v(3) = - abs(vel(3))  ! to keep it as negative velocity   
!         end if
! 
!         PC_particles(PC_num_part)%a    = - PC_ion_particles(ll)%a * PC_particle_ion_mass / PC_particle_elec_mass 
!         ! accel can be decided from ion's accel
!         
!     end subroutine PC_create_second_elec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************************
!****************************************************************************************************
   !> Routine for MPI

   !> MPI: Share the electrons evenly over the tasks, with no domain decomposition
    subroutine shareParticlesElectrons(myrank, ntasks)
        include "mpif.h"
        integer, intent(IN) :: myrank, ntasks

        integer :: i_more, i_less, n_send
        integer :: ierr, mean_n_part, i_min, i_max, tag
        integer, dimension(0:ntasks-1) :: n_part_task

!          print *, myrank, "Electron: Before sharing: PC_num_part", PC_num_part
        ! Get the number of particles each task has
        call MPI_ALLGATHER(PC_num_part, 1, MPI_integer, n_part_task, 1, &
            MPI_integer, MPI_COMM_WORLD, ierr)
        mean_n_part = (sum(n_part_task) + ntasks - 1) / ntasks

        i_less = 0
        i_more = 0

        do while (i_more < ntasks)
        if (n_part_task(i_more) > mean_n_part) then

            ! Find first task with less than mean particles
            do while (n_part_task(i_less) >= mean_n_part)
                i_less = i_less + 1
            end do

            ! Determine amount of particles to send
            n_send = min(n_part_task(i_more) - mean_n_part, &
                mean_n_part - n_part_task(i_less))

            if (myrank == i_more) then ! Send
                i_min = n_part_task(i_more) - n_send + 1
                i_max = i_min + n_send - 1
                tag   = i_more
                call send_elec_part_mpi(PC_particles(i_min:i_max), i_less, tag)
            else if (myrank == i_less) then ! Receive
                i_min = n_part_task(i_less) + 1
                i_max = i_min + n_send - 1
                call check_num_particles(i_max)
                tag   = i_more
                call recv_elec_part_mpi(PC_particles(i_min:i_max), i_more, tag)
            end if

            ! Update the number of particles for each task
            n_part_task(i_less) = n_part_task(i_less) + n_send
            n_part_task(i_more) = n_part_task(i_more) - n_send
        end if

        if (n_part_task(i_more) <= mean_n_part) i_more = i_more + 1
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        PC_num_part = n_part_task(myrank)
!          print *, myrank, " Electron: After sharing: PC_num_part", PC_num_part
    end subroutine shareParticlesElectrons


    subroutine send_elec_part_mpi(parts, recver, tag)
        include "mpif.h"
        type(PC_part_t), intent(in) :: parts(:)
        integer, intent(in) :: recver, tag

        integer :: n_send, req_size, rel_size, ierr
        real(dp), allocatable :: r_buf(:)

        ! We use a buffer because we cannot send the part type *neatly* with mpi
        rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
        n_send   = size(parts)
        req_size = n_send * rel_size
        allocate(r_buf(req_size))

        r_buf = transfer(parts, r_buf, req_size)
        call MPI_SEND(r_buf, req_size, MPI_DOUBLE_PRECISION, recver, tag, &
            & MPI_COMM_WORLD, ierr)
    end subroutine send_elec_part_mpi

    subroutine recv_elec_part_mpi(parts, sender, tag)
        include "mpif.h"
        type(PC_part_t), intent(inout) :: parts(:)
        integer, intent(in) :: sender, tag

        integer :: n_send, req_size, rel_size, ierr
        integer :: mpi_status(MPI_STATUS_SIZE)
        real(dp), allocatable :: r_buf(:)

        ! We use a buffer because we cannot send the part type *neatly* with mpi
        n_send   = size(parts)
        rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
        req_size = n_send * rel_size
        allocate(r_buf(req_size))
        call MPI_RECV(r_buf, req_size, MPI_DOUBLE_PRECISION, sender, tag, &
            & MPI_COMM_WORLD, mpi_status, ierr)
        parts = transfer(r_buf, parts, n_send)
    end subroutine recv_elec_part_mpi


   !!!!!!!!!!!!!!!!!!!For ions!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> MPI: Share the IONS evenly over the tasks, with no domain decomposition
    subroutine shareParticlesIons(myrank, ntasks)
        include "mpif.h"
        integer, intent(IN) :: myrank, ntasks

        integer :: i_more, i_less, n_send
        integer :: ierr, mean_n_part, i_min, i_max, tag
        integer, dimension(0:ntasks-1) :: n_part_task

       !  print *, myrank, "Ion: Before sharing: PC_num_ion_part", PC_num_ion_part
        ! Get the number of particles each task has
        call MPI_ALLGATHER(PC_num_ion_part, 1, MPI_integer, n_part_task, 1, &
            MPI_integer, MPI_COMM_WORLD, ierr)
        mean_n_part = (sum(n_part_task) + ntasks - 1) / ntasks

        i_less = 0
        i_more = 0

        do while (i_more < ntasks)
        if (n_part_task(i_more) > mean_n_part) then

            ! Find first task with less than mean particles
            do while (n_part_task(i_less) >= mean_n_part)
                i_less = i_less + 1
            end do

            ! Determine amount of particles to send
            n_send = min(n_part_task(i_more) - mean_n_part, &
                mean_n_part - n_part_task(i_less))

            if (myrank == i_more) then ! Send
                i_min = n_part_task(i_more) - n_send + 1
                i_max = i_min + n_send - 1
                tag   = i_more
                call send_ion_part_mpi(PC_ion_particles(i_min:i_max), i_less, tag)
            else if (myrank == i_less) then ! Receive
                i_min = n_part_task(i_less) + 1
                i_max = i_min + n_send - 1
                call check_num_particles(i_max)
                tag   = i_more
                call recv_ion_part_mpi(PC_ion_particles(i_min:i_max), i_more, tag)
            end if

            ! Update the number of particles for each task
            n_part_task(i_less) = n_part_task(i_less) + n_send
            n_part_task(i_more) = n_part_task(i_more) - n_send
        end if

        if (n_part_task(i_more) <= mean_n_part) i_more = i_more + 1
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        PC_num_ion_part = n_part_task(myrank)
       !  print *, myrank, "Ions: After sharing: PC_num_ion_part", PC_num_ion_part
    end subroutine shareParticlesIons


    subroutine send_ion_part_mpi(parts, recver, tag)
        include "mpif.h"
        type(PC_ion_part_t), intent(in) :: parts(:)
        integer, intent(in) :: recver, tag

        integer :: n_send, req_size, rel_size, ierr
        real(dp), allocatable :: r_buf(:)

        ! We use a buffer because we cannot send the part type *neatly* with mpi
        rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
        n_send   = size(parts)
        req_size = n_send * rel_size
        allocate(r_buf(req_size))

        r_buf = transfer(parts, r_buf, req_size)
        call MPI_SEND(r_buf, req_size, MPI_DOUBLE_PRECISION, recver, tag, &
            & MPI_COMM_WORLD, ierr)
    end subroutine send_ion_part_mpi

    subroutine recv_ion_part_mpi(parts, sender, tag)
        include "mpif.h"
        type(PC_ion_part_t), intent(inout) :: parts(:)
        integer, intent(in) :: sender, tag

        integer :: n_send, req_size, rel_size, ierr
        integer :: mpi_status(MPI_STATUS_SIZE)
        real(dp), allocatable :: r_buf(:)

        ! We use a buffer because we cannot send the part type *neatly* with mpi
        n_send   = size(parts)
        rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
        req_size = n_send * rel_size
        allocate(r_buf(req_size))
        call MPI_RECV(r_buf, req_size, MPI_DOUBLE_PRECISION, sender, tag, &
            & MPI_COMM_WORLD, mpi_status, ierr)
        parts = transfer(r_buf, parts, n_send)
    end subroutine recv_ion_part_mpi

!> Return the number of real particles
    real(dp) function PC_getNumRealelecPartMPI()
        include "mpif.h"
        integer  :: ierr
        real(dp) :: part_sum
        part_sum = sum(PC_particles(1:PC_num_part)%weight)
        call MPI_ALLREDUCE(part_sum, PC_getNumRealelecPartMPI, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
    end function PC_getNumRealelecPartMPI

    !> Return the number of simulation particles
    integer function PC_getNumSimelecPartMPI()
        include "mpif.h"
        integer :: ierr
        call MPI_ALLREDUCE(PC_num_part, PC_getNumSimelecPartMPI, 1, MPI_integer, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
    end function PC_getNumSimelecPartMPI

    !> Return the number of real particles
    real(dp) function PC_getNumRealionPartMPI()
        include "mpif.h"
        integer  :: ierr
        real(dp) :: part_sum
        part_sum = sum(PC_ion_particles(1:PC_num_ion_part)%weight)
        call MPI_ALLREDUCE(part_sum, PC_getNumRealionPartMPI, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
    end function PC_getNumRealionPartMPI

    !> Return the number of simulation particles
    integer function PC_getNumSimionPartMPI()
        include "mpif.h"
        integer :: ierr
        call MPI_ALLREDUCE(PC_num_ion_part, PC_getNumSimionPartMPI, 1, MPI_integer, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
    end function PC_getNumSimionPartMPI

!Anbang: here i save the informaton of particle during a certain time steps, for the convenience of later computing
    subroutine PC_save_elec_part_info(myrank, root, ntasks, sim_name)
        include "mpif.h"
        integer, intent(IN) :: myrank, ntasks, root

        integer :: ierr, total_part, i_min, i_max, tag
        integer, dimension(0:ntasks-1) :: n_part_task
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer                 :: myunit, i

        ! we do this for electrons
        ! Get the number of electrons particles each task has
        call MPI_ALLGATHER(PC_num_part, 1, MPI_integer, n_part_task, 1, &
            MPI_integer, MPI_COMM_WORLD, ierr)

        total_part = sum(n_part_task)
        allocate(PC_particles_collect(total_part))

        tag = 1

        do while (tag < ntasks)   ! transfer particle information to the root
            if (myrank == tag) then
                call send_elec_part_mpi(PC_particles(1:n_part_task(myrank)), root, tag)  ! send to root
            else if (myrank == root) then
                i_min = sum(n_part_task(0:tag -1)) + 1
                i_max = i_min + n_part_task(tag) - 1
                call recv_elec_part_mpi(PC_particles_collect(i_min:i_max), tag, tag)
            end if

            tag = tag + 1
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (myrank == root) then 
            ! direction copy at the root
            PC_particles_collect(1:n_part_task(myrank)) = PC_particles(1:n_part_task(myrank))

            myunit = 299
            filename = "midstore/" // trim(sim_name) // "_info_elec"// ".txt"
            open(unit = myunit, file = filename , status = 'unknown')
            write(myunit, *) total_part
            do i = 1, total_part
                write(myunit, *) PC_particles_collect(i)%weight, PC_particles_collect(i)%x(3), PC_particles_collect(i)%v
            end do
            close(myunit)
         end if

        deallocate(PC_particles_collect)

    end subroutine PC_save_elec_part_info

! save ion particle info
    subroutine PC_save_ion_part_info(myrank, root, ntasks, sim_name)
        include "mpif.h"
        integer, intent(IN) :: myrank, ntasks, root

        integer :: ierr, total_part, i_min, i_max, tag
        integer, dimension(0:ntasks-1) :: n_part_task
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer                 :: myunit, i

        ! we do this for ions
        ! Get the number of ions particles each task has
        call MPI_ALLGATHER(PC_num_ion_part, 1, MPI_integer, n_part_task, 1, &
            MPI_integer, MPI_COMM_WORLD, ierr)

        total_part = sum(n_part_task)
        allocate(PC_ion_particles_collect(total_part))

        tag = 1

        do while (tag < ntasks)   ! transfer particle information to the root
            if (myrank == tag) then
                call send_ion_part_mpi(PC_ion_particles(1:n_part_task(myrank)), root, tag)  ! send to root
            else if (myrank == root) then
                i_min = sum(n_part_task(0:tag -1)) + 1
                i_max = i_min + n_part_task(tag) - 1
                call recv_ion_part_mpi(PC_ion_particles_collect(i_min:i_max), tag, tag)
            end if

            tag = tag + 1
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if (myrank == root) then 
            ! direction copy at the root
            PC_ion_particles_collect(1:n_part_task(myrank)) = PC_ion_particles(1:n_part_task(myrank))

            myunit = 300
            write(filename, fmt="(A)") &
                "midstore/" // trim(sim_name) // "_info_ion"//".txt"
            open(unit = myunit, file = filename , status = 'unknown')
            write(myunit, *) total_part
            do i = 1, total_part
                write(myunit, *) PC_ion_particles_collect(i)%weight, PC_ion_particles_collect(i)%x(3), &
                    PC_ion_particles_collect(i)%v
            end do
            close(myunit)
         end if

        deallocate(PC_ion_particles_collect)

    end subroutine PC_save_ion_part_info

    
    ! read elec info from saved files
    subroutine PC_read_elec_part_info(sim_name)
        integer :: total_part
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer                 :: myunit, i

        
        myunit = 299
        filename = "midstore/" // trim(sim_name) // "_info_elec"// ".txt"
        open(unit = myunit, file = filename , status = 'old')
            read(myunit, *) total_part
            if (total_part > PC_max_num_part) then
                print *, "there are too many electrons when read the particle infor: PC_read_elec_part_info, re-allocate!"
                PC_max_num_part = int(1.5d0 * total_part) 
                deallocate(PC_particles)
                allocate(PC_particles(PC_max_num_part))
            end if

            PC_num_part = total_part
            PC_particles(1:PC_num_part)%live = .true.
            do i = 1,3
                PC_particles(1:PC_num_part)%x(i)    = 0.d0
                PC_particles(1:PC_num_part)%v(i)    = 0.d0
                PC_particles(1:PC_num_part)%a(i)    = 0.d0
            end do
            PC_particles(1:PC_num_part)%T_left = 0.d0

            !anbang: at each node, we read certain amount of particles
            do i = 1, total_part
                read(myunit, *) PC_particles(i)%weight, PC_particles(i)%x(3), PC_particles(i)%v
            end do
        close(myunit)
    end subroutine PC_read_elec_part_info

    ! ions read particle info
    subroutine PC_read_ion_part_info(sim_name)
        integer :: total_part
        character(len=*), intent(in) :: sim_name
        character(len=name_len) :: filename
        integer                 :: myunit, i


        
        myunit = 300
        write(filename, fmt="(A)") &
                "midstore/" // trim(sim_name) // "_info_ion"//".txt"
        open(unit = myunit, file = filename , status = 'old')
            read(myunit, *) total_part    
            if (total_part > PC_max_num_ion_part) then
                print *, "there are too many ions when read the particle infor: PC_read_ion_part_info, re-Allocate!"
                PC_max_num_ion_part = int(1.5d0 * total_part) 
                deallocate(PC_ion_particles)
                allocate(PC_ion_particles(PC_max_num_ion_part))
            end if

            PC_num_ion_part = total_part
            PC_ion_particles(1:PC_num_ion_part)%live = .true.
            do i = 1,3
                PC_ion_particles(1:PC_num_ion_part)%x(i)    = 0.d0
                PC_ion_particles(1:PC_num_ion_part)%v(i)    = 0.d0
                PC_ion_particles(1:PC_num_ion_part)%a(i)    = 0.d0
            end do
            PC_ion_particles(1:PC_num_ion_part)%T_left = 0.d0

            do i = 1, total_part
                read(myunit, *) PC_ion_particles(i)%weight, PC_ion_particles(i)%x(3), PC_ion_particles(i)%v
            end do
        close(myunit)
    end subroutine PC_read_ion_part_info


!!!!Anbang: from here, i calculate the info of  particle number per cell/ mean particle weight per cell
    subroutine PC_cal_elec_num_and_weight_per_cell(tot_elec_sim_num_per_cell, tot_elec_mean_weight_per_cell)
        include "mpif.h"
        integer :: ierr
        integer :: elec_sim_num_per_cell(PC_grid_size - 1)  ! electron number per cell
        integer :: i, low_ix
        real(dp) :: temp, weight, inv_delta_x
        real(dp) :: elec_real_num_per_cell(PC_grid_size - 1), tot_elec_real_num_per_cell(PC_grid_size - 1)
        integer, intent(out) :: tot_elec_sim_num_per_cell(PC_grid_size - 1)
        real(dp), intent(out) :: tot_elec_mean_weight_per_cell(PC_grid_size - 1)
        

        elec_sim_num_per_cell = 0
        elec_real_num_per_cell = 0.d0
        tot_elec_mean_weight_per_cell = 0.d0
        tot_elec_real_num_per_cell = 0.d0
        tot_elec_sim_num_per_cell = 0

        inv_delta_x = 1.d0 / PC_delta_x
        do i = 1,PC_num_part
            ! Index i is at (i-1) * dx
            temp   = PC_particles(i)%x(3) * inv_delta_x
            low_ix = floor(temp) + 1

            if (low_ix < 1 .or. low_ix > PC_grid_size-1) then
                print *, "position", PC_particles(i)%x(3)
                print *, "PC_cal_elec_num_and_weight_per_cell: a particle is outside the computational domain"
            end if
            
            elec_sim_num_per_cell(low_ix) = elec_sim_num_per_cell(low_ix) + 1
            elec_real_num_per_cell(low_ix) = elec_real_num_per_cell(low_ix) + PC_particles(i)%weight          
        end do

        ! collect and dristribute to all nodes
        call MPI_ALLREDUCE(elec_sim_num_per_cell, tot_elec_sim_num_per_cell, PC_grid_size - 1, MPI_integer, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(elec_real_num_per_cell, tot_elec_real_num_per_cell, PC_grid_size - 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, ierr)

        do i = 1, PC_grid_size - 1
            if (tot_elec_sim_num_per_cell(i) /= 0) then
                tot_elec_mean_weight_per_cell(i) = tot_elec_real_num_per_cell(i) / tot_elec_sim_num_per_cell(i)
                !if (tot_elec_mean_weight_per_cell(i) < 1.d0) tot_elec_mean_weight_per_cell(i) = 1.d0
            end if
         end do
    end subroutine PC_cal_elec_num_and_weight_per_cell

    subroutine PC_cal_ion_num_and_weight_per_cell(tot_ion_sim_num_per_cell, tot_ion_mean_weight_per_cell)
        include "mpif.h"
        integer :: ierr
        integer :: ion_sim_num_per_cell(PC_grid_size - 1)  ! ion number per cell
        integer :: i, low_ix
        real(dp) :: temp, weight, inv_delta_x
        real(dp) :: ion_real_num_per_cell(PC_grid_size - 1), tot_ion_real_num_per_cell(PC_grid_size - 1)
        integer, intent(out) :: tot_ion_sim_num_per_cell(PC_grid_size - 1)
        real(dp), intent(out) :: tot_ion_mean_weight_per_cell(PC_grid_size - 1)
        

        ion_sim_num_per_cell = 0
        ion_real_num_per_cell = 0.d0
        tot_ion_mean_weight_per_cell = 0.d0
        tot_ion_real_num_per_cell = 0.d0
        tot_ion_sim_num_per_cell = 0

        inv_delta_x = 1.d0 / PC_delta_x
        do i = 1,PC_num_ion_part
            ! Index i is at (i-1) * dx
            temp   = PC_ion_particles(i)%x(3) * inv_delta_x
            low_ix = floor(temp) + 1

            if (low_ix < 1 .or. low_ix > PC_grid_size-1) then
                print *, "position", PC_ion_particles(i)%x(3)
                print *, "PC_cal_ion_num_and_weight_per_cell: a particle is outside the computational domain"
            end if
            
            ion_sim_num_per_cell(low_ix) = ion_sim_num_per_cell(low_ix) + 1
            ion_real_num_per_cell(low_ix) = ion_real_num_per_cell(low_ix) + PC_ion_particles(i)%weight          
        end do

        ! collect and dristribute to all nodes
        call MPI_ALLREDUCE(ion_sim_num_per_cell, tot_ion_sim_num_per_cell, PC_grid_size - 1, MPI_integer, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ion_real_num_per_cell, tot_ion_real_num_per_cell, PC_grid_size - 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, ierr)

        do i = 1, PC_grid_size - 1
            if (tot_ion_sim_num_per_cell(i) /= 0) then
                tot_ion_mean_weight_per_cell(i) = tot_ion_real_num_per_cell(i) / tot_ion_sim_num_per_cell(i)
                !if (tot_ion_mean_weight_per_cell(i) < 1.d0) tot_ion_mean_weight_per_cell(i) = 1.d0
            end if
         end do
    end subroutine PC_cal_ion_num_and_weight_per_cell


    
end module m_particle_core
